################################################################################
# Cassiopea microbiome project
# Strobilation-focused bacterial community analysis comparing aposymbiotic, mutant, and strobilating polyp outcomes using PCoA, Wilcoxon tests, Cohen’s d, and log2 fold changes..
################################################################################

# ============================================================================
# CASSIOPEA STROBILATION ANALYSIS 
# ============================================================================

required_packages <- c(
  "phyloseq",
  "vegan",
  "ggplot2",
  "dplyr",
  "tidyr",
  "ggpubr",
  "rstatix",
  "patchwork",
  "scales",
  "ggtext",
  "multcompView"
)


# ============================================================================
# SETUP
# ============================================================================

# ---------------------------------------------------------------------------
# Path configuration
# ---------------------------------------------------------------------------
project_root <- "."
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "03_beta_diversity/strobilation_pcoa_effect_sizes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
phyloseq_bact_file <- file.path(data_dir, "bacteria_phyloseq.rds")

cat("Loading phyloseq object...
")
ps_filt <- readRDS(phyloseq_bact_file)

cat("Loaded:", nsamples(ps_filt), "samples,", ntaxa(ps_filt), "ASVs

")

# ============================================================================
# DATA PREPARATION
# ============================================================================

cat("=== DATA PREPARATION ===\n")

create_strobilation_groups <- function(ps_obj) {
  meta <- data.frame(sample_data(ps_obj))
  
  meta$Symbiosis_Outcome <- case_when(
    meta$Treatment == "Polyp_Aposymbiotic" ~ "Aposymbiotic",
    meta$Treatment == "Polyp_Mutant" ~ "Mutant",
    meta$Treatment %in% c("Polyp_Native", 
                          "Polyp_Control",
                          "Polyp_Antibiotic") ~ "Strobilation",
    TRUE ~ "Exclude"
  )
  
  sample_data(ps_obj) <- meta
  
  cat("\nGrouping results:\n")
  print(table(meta$Symbiosis_Outcome))
  
  return(ps_obj)
}

ps_grouped <- create_strobilation_groups(ps_filt)

# ============================================================================
# FILTER DATA
# ============================================================================

cat("\n=== FILTERING ===\n")

ps_polyp <- subset_samples(ps_grouped, Symbiosis_Outcome != "Exclude")

meta_pre <- data.frame(sample_data(ps_polyp))
cat("\nSamples before read filtering:\n")
print(table(meta_pre$Symbiosis_Outcome))

cat("\nSamples with <5000 reads that will be removed:\n")
low_reads <- meta_pre[sample_sums(ps_polyp) < 5000, c("Symbiosis_Outcome")]
if(length(low_reads) > 0) {
  print(table(low_reads))
} else {
  cat("  None\n")
}

ps_polyp <- prune_samples(sample_sums(ps_polyp) >= 5000, ps_polyp)
ps_polyp <- prune_taxa(taxa_sums(ps_polyp) > 0, ps_polyp)

meta <- data.frame(sample_data(ps_polyp))
sample_counts <- table(meta$Symbiosis_Outcome)

cat("\n FINAL sample distribution:\n")
print(sample_counts)
cat("Total samples:", nsamples(ps_polyp), "\n")
cat("Total ASVs:", ntaxa(ps_polyp), "\n\n")

ps_rel <- transform_sample_counts(ps_polyp, function(x) x/sum(x))

# ============================================================================
# PCoA
# ============================================================================

cat("=== PCoA ===\n")

dist_bray <- phyloseq::distance(ps_rel, method = "bray")

pcoa <- cmdscale(dist_bray, eig = TRUE, k = 3)
eigenvals <- pcoa$eig[pcoa$eig > 0]
var_exp <- eigenvals / sum(eigenvals) * 100

pcoa_df <- data.frame(
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2],
  Sample = rownames(pcoa$points)
)
pcoa_df <- merge(pcoa_df, meta, by.x = "Sample", by.y = "row.names")

cat("Variance explained - PC1:", round(var_exp[1], 1), "%, PC2:", round(var_exp[2], 1), "%\n")

set.seed(123)
permanova <- adonis2(dist_bray ~ Symbiosis_Outcome, data = meta, permutations = 999)
cat("\nPERMANOVA Results:\n")
print(permanova)

dispersion <- betadisper(dist_bray, meta$Symbiosis_Outcome)
dispersion_test <- permutest(dispersion, permutations = 999)
cat("\nDispersion test P-value:", round(dispersion_test$tab$`Pr(>F)`[1], 4), "\n")

panel_A <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Symbiosis_Outcome, shape = Symbiosis_Outcome)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = Symbiosis_Outcome), linetype = "dashed", alpha = 0.7, level = 0.95) +
  scale_color_manual(values = c("Aposymbiotic" = "#A23B72", 
                                "Mutant" = "#DE8F05",
                                "Strobilation" = "#029E73")) +
  scale_shape_manual(values = c("Aposymbiotic" = 16, 
                                "Mutant" = 17,
                                "Strobilation" = 15)) +
  labs(x = paste0("Axis1 (", round(var_exp[1], 1), "%)"),
       y = paste0("Axis2 (", round(var_exp[2], 1), "%)"),
       color = "", shape = "") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14))

ggsave(file.path(output_dir, "Panel_A_PCoA.png"), panel_A, width = 8, height = 7, dpi = 300)

# ============================================================================
# DIFFERENTIAL ABUNDANCE ANALYSIS - WITH BOTH EFFECT SIZES
# ============================================================================

cat("=== DIFFERENTIAL ABUNDANCE ANALYSIS (Cohen's d + Log2FC) ===\n")

perform_wilcoxon_analysis <- function(ps_rel, group1, group2, comparison_name, 
                                      min_prevalence = 0.1) {
  
  cat("Comparing:", group1, "vs", group2, "\n")
  
  meta <- data.frame(sample_data(ps_rel))
  samples_g1 <- rownames(meta)[meta$Symbiosis_Outcome == group1]
  samples_g2 <- rownames(meta)[meta$Symbiosis_Outcome == group2]
  
  cat("  Group 1 (", group1, "):", length(samples_g1), "samples\n")
  cat("  Group 2 (", group2, "):", length(samples_g2), "samples\n")
  
  otu_mat <- as(otu_table(ps_rel), "matrix")
  if (taxa_are_rows(ps_rel)) otu_mat <- t(otu_mat)
  
  total_samples <- length(c(samples_g1, samples_g2))
  min_occurrences <- ceiling(min_prevalence * total_samples)
  
  # Calculate pseudocount for log2FC 
  non_zero_values <- otu_mat[otu_mat > 0]
  pseudocount <- min(non_zero_values) / 2
  
  results <- data.frame()
  
  for (asv in colnames(otu_mat)) {
    g1_abund <- otu_mat[samples_g1, asv]
    g2_abund <- otu_mat[samples_g2, asv]
    combined <- c(g1_abund, g2_abund)
    
    if (sum(combined > 0) >= min_occurrences) {
      
      test <- tryCatch({
        wilcox.test(g1_abund, g2_abund, exact = FALSE)
      }, error = function(e) list(p.value = 1.0))
      
      mean1 <- mean(g1_abund)
      mean2 <- mean(g2_abund)
      median1 <- median(g1_abund)
      median2 <- median(g2_abund)
      
      # Cohen's d with pooled SD (accounts for unequal sample sizes)
      pooled_sd <- sqrt(((length(g1_abund)-1)*var(g1_abund) + 
                           (length(g2_abund)-1)*var(g2_abund)) / 
                          (length(g1_abund) + length(g2_abund) - 2))
      cohens_d <- ifelse(pooled_sd > 0, (mean1 - mean2) / pooled_sd, 0)
      
      # Log2 fold-change with pseudocount
      log2fc <- log2((mean1 + pseudocount) / (mean2 + pseudocount))
      
      results <- rbind(results, data.frame(
        ASV = asv,
        Comparison = comparison_name,
        P_value = test$p.value,
        Median_Group1 = median1,
        Median_Group2 = median2,
        Median_Difference = median1 - median2,
        Mean_Group1 = mean1,
        Mean_Group2 = mean2,
        Cohens_D = cohens_d,
        Log2FC = log2fc,
        Prevalence = sum(combined > 0) / total_samples
      ))
    }
  }
  
  results$P_adjusted <- p.adjust(results$P_value, method = "BH")
  
  tax_df <- data.frame(tax_table(ps_rel))
  results <- merge(results, tax_df, by.x = "ASV", by.y = "row.names", all.x = TRUE)
  
  results <- results[order(results$P_adjusted, -abs(results$Cohens_D)), ]
  
  
  return(results)
}

wilcox_strob_vs_apo <- perform_wilcoxon_analysis(ps_rel, "Strobilation", "Aposymbiotic", 
                                                 "Strobilation_vs_Aposymbiotic")
wilcox_strob_vs_mut <- perform_wilcoxon_analysis(ps_rel, "Strobilation", "Mutant", 
                                                 "Strobilation_vs_Mutant")
wilcox_mut_vs_apo <- perform_wilcoxon_analysis(ps_rel, "Mutant", "Aposymbiotic", 
                                               "Mutant_vs_Aposymbiotic")

# ============================================================================
# FILTERING BIOMARKERS - USING COHEN'S D THRESHOLDS
# ============================================================================

filter_biomarkers <- function(results, p_thresh = 0.05, effect_thresh = 0.001, 
                              cohens_thresh = 0.2, log2fc_thresh = NULL) {
  
  sig <- results[
    results$P_adjusted < p_thresh & 
      abs(results$Median_Difference) > effect_thresh &
      abs(results$Cohens_D) > cohens_thresh,
  ]
  
  # also filter by log2FC
  if (!is.null(log2fc_thresh)) {
    sig <- sig[abs(sig$Log2FC) > log2fc_thresh, ]
  }
  
  return(sig)
}

sig_strob_vs_apo <- filter_biomarkers(wilcox_strob_vs_apo)
sig_strob_vs_mut <- filter_biomarkers(wilcox_strob_vs_mut)
sig_mut_vs_apo <- filter_biomarkers(wilcox_mut_vs_apo)

cat("\nSignificant biomarkers:\n")
cat("  Strobilation vs Aposymbiotic:", nrow(sig_strob_vs_apo), "\n")
cat("  Strobilation vs Mutant:", nrow(sig_strob_vs_mut), "\n")
cat("  Mutant vs Aposymbiotic:", nrow(sig_mut_vs_apo), "\n\n")

# Save ALL results 
write.csv(wilcox_strob_vs_apo, file.path(output_dir, "Strobilation_vs_Aposymbiotic_All.csv"), row.names = FALSE)
write.csv(wilcox_strob_vs_mut, file.path(output_dir, "Strobilation_vs_Mutant_All.csv"), row.names = FALSE)
write.csv(wilcox_mut_vs_apo, file.path(output_dir, "Mutant_vs_Aposymbiotic_All.csv"), row.names = FALSE)

# Save SIGNIFICANT results 
write.csv(sig_strob_vs_apo, file.path(output_dir, "Strobilation_vs_Aposymbiotic_Significant.csv"), row.names = FALSE)
write.csv(sig_strob_vs_mut, file.path(output_dir, "Strobilation_vs_Mutant_Significant.csv"), row.names = FALSE)
write.csv(sig_mut_vs_apo, file.path(output_dir, "Mutant_vs_Aposymbiotic_Significant.csv"), row.names = FALSE)

# ============================================================================
# BIOMARKER BOXPLOTS 
# ============================================================================

cat("=== BIOMARKER BOXPLOTS ===\n")

create_enhanced_boxplot <- function(biomarker_results, ps_rel_obj, top_n = 4, 
                                    letter_offset = 0.05) {
  
  if (nrow(biomarker_results) == 0) {
    cat("No significant biomarkers found\n")
    return(NULL)
  }
  
  biomarker_results$Score <- -log10(biomarker_results$P_adjusted) * abs(biomarker_results$Cohens_D)
  top_taxa <- biomarker_results[order(biomarker_results$Score, decreasing = TRUE), ][1:min(top_n, nrow(biomarker_results)), ]
  
  cat("Selected top", nrow(top_taxa), "biomarkers\n")
  
  otu_mat <- as(otu_table(ps_rel_obj), "matrix")
  if (taxa_are_rows(ps_rel_obj)) otu_mat <- t(otu_mat)
  
  meta <- data.frame(sample_data(ps_rel_obj))
  
  plot_data <- data.frame()
  letter_data <- data.frame()
  
  for (i in 1:nrow(top_taxa)) {
    asv <- top_taxa$ASV[i]
    if (asv %in% colnames(otu_mat)) {
      
      family <- ifelse(is.na(top_taxa$Family[i]) | top_taxa$Family[i] == "", 
                       "Unclassified", top_taxa$Family[i])
      genus <- ifelse(is.na(top_taxa$Genus[i]) | top_taxa$Genus[i] == "", 
                      "Unclassified", top_taxa$Genus[i])
      
      asv_data <- data.frame(
        Sample = rownames(otu_mat),
        Abundance = otu_mat[, asv] * 100,
        ASV = asv,
        Family = family,
        Genus = genus,
        Cohens_D = top_taxa$Cohens_D[i],
        Log2FC = top_taxa$Log2FC[i]
      )
      
      asv_data <- merge(asv_data, meta, by.x = "Sample", by.y = "row.names")
      
      asv_data$Symbiosis_Outcome <- factor(asv_data$Symbiosis_Outcome,
                                           levels = c("Aposymbiotic", "Mutant", "Strobilation"))
      
      # Pairwise Wilcoxon tests
      pw_result <- pairwise.wilcox.test(
        asv_data$Abundance, 
        asv_data$Symbiosis_Outcome,
        p.adjust.method = "BH",
        exact = FALSE
      )
      
      # Create full p-value matrix
      groups <- levels(asv_data$Symbiosis_Outcome)
      full_matrix <- matrix(1, nrow = length(groups), ncol = length(groups))
      rownames(full_matrix) <- colnames(full_matrix) <- groups
      
      pw_matrix <- pw_result$p.value
      for (j in 1:(length(groups)-1)) {
        for (k in (j+1):length(groups)) {
          if (!is.na(pw_matrix[k-1, j])) {
            full_matrix[j, k] <- full_matrix[k, j] <- pw_matrix[k-1, j]
          }
        }
      }
      
      # Generate compact letters
      letters_result <- multcompLetters(full_matrix, threshold = 0.05, 
                                        compare = "<", Letters = letters)
      
      # Fixed distance approach
      data_range <- max(asv_data$Abundance) - min(asv_data$Abundance)
      fixed_offset <- data_range * letter_offset
      
      max_vals <- asv_data %>%
        group_by(Symbiosis_Outcome) %>%
        summarize(MaxVal = max(Abundance) + fixed_offset, .groups = "drop")
      
      letters_df <- data.frame(
        Symbiosis_Outcome = names(letters_result$Letters),
        Letter = letters_result$Letters,
        ASV = asv,
        Family = family,
        Genus = genus,
        Cohens_D = top_taxa$Cohens_D[i],
        Log2FC = top_taxa$Log2FC[i]
      )
      letters_df <- merge(letters_df, max_vals, by = "Symbiosis_Outcome")
      
      plot_data <- rbind(plot_data, asv_data)
      letter_data <- rbind(letter_data, letters_df)
    }
  }
  
  # Create facet labels (keeping Cohen's d 
  plot_data$Facet_Label <- paste0(
    plot_data$Family, "\n*", plot_data$Genus, "*\n(d = ", 
    sprintf("%.2f", plot_data$Cohens_D), ")"
  )
  
  letter_data$Facet_Label <- paste0(
    letter_data$Family, "\n*", letter_data$Genus, "*\n(d = ", 
    sprintf("%.2f", letter_data$Cohens_D), ")"
  )
  
  sample_sizes <- table(meta$Symbiosis_Outcome)
  
  # Create plot
  panel_B <- ggplot(plot_data, aes(x = Symbiosis_Outcome, y = Abundance, fill = Symbiosis_Outcome)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
    geom_text(data = letter_data,
              aes(x = Symbiosis_Outcome, y = MaxVal, label = Letter),
              inherit.aes = FALSE, size = 6) +
    facet_wrap(~ Facet_Label, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Aposymbiotic" = "#A23B72", 
                                 "Mutant" = "#DE8F05",
                                 "Strobilation" = "#029E73")) +
    scale_x_discrete(labels = c("Aposymbiotic" = paste0("Aposymbiotic\n(n=", sample_sizes["Aposymbiotic"], ")"),
                                "Mutant" = paste0("Mutant\n(n=", sample_sizes["Mutant"], ")"),
                                "Strobilation" = paste0("Symbiotic\n(n=", sample_sizes["Strobilation"], ")"))) +
    labs(x = "", y = "Abund. (%)") +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          legend.position = "none",
          strip.text = element_markdown(size = 13),
          strip.background = element_rect(fill = "gray95", color = "gray70"))
  
  
  
  return(panel_B)
}

panel_B <- create_enhanced_boxplot(sig_strob_vs_mut, ps_rel, top_n = 4, 
                                   letter_offset = 0.08)

if (!is.null(panel_B)) {
  ggsave(file.path(output_dir, "Panel_B_Biomarkers.png"), 
         panel_B, width = 10, height = 8, dpi = 300)
  
}

# ============================================================================
# COHEN'S D EFFECT SIZES 
# ============================================================================

cat("=== COHEN'S D EFFECT SIZES ===\n")

create_cohens_d_plot <- function(strob_vs_apo, strob_vs_mut) {
  
  sig_apo <- filter_biomarkers(strob_vs_apo)
  sig_mut <- filter_biomarkers(strob_vs_mut)
  
  all_biomarkers <- unique(c(sig_apo$ASV, sig_mut$ASV))
  
  if (length(all_biomarkers) == 0) {
    cat("No significant biomarkers for Cohen's d plot\n")
    return(NULL)
  }
  
  cat("Creating effect size plot for", length(all_biomarkers), "biomarkers\n")
  
  plot_data <- data.frame()
  
  for (asv in all_biomarkers) {
    genus_name <- NA
    if (asv %in% strob_vs_apo$ASV) {
      genus_name <- strob_vs_apo$Genus[strob_vs_apo$ASV == asv][1]
    } else if (asv %in% strob_vs_mut$ASV) {
      genus_name <- strob_vs_mut$Genus[strob_vs_mut$ASV == asv][1]
    }
    
    genus_name <- ifelse(is.na(genus_name) | genus_name == "", "Unclassified", genus_name)
    
    cohens_apo <- ifelse(asv %in% strob_vs_apo$ASV, 
                         strob_vs_apo$Cohens_D[strob_vs_apo$ASV == asv][1], NA)
    cohens_mut <- ifelse(asv %in% strob_vs_mut$ASV, 
                         strob_vs_mut$Cohens_D[strob_vs_mut$ASV == asv][1], NA)
    
    if (!is.na(cohens_apo)) {
      plot_data <- rbind(plot_data, data.frame(
        Genus = genus_name,
        Cohens_D = cohens_apo,
        Comparison = "Strobilation vs Aposymbiotic",
        ASV = asv
      ))
    }
    if (!is.na(cohens_mut)) {
      plot_data <- rbind(plot_data, data.frame(
        Genus = genus_name,
        Cohens_D = cohens_mut,
        Comparison = "Strobilation vs Mutant",
        ASV = asv
      ))
    }
  }
  
  genus_order <- plot_data %>%
    group_by(Genus) %>%
    summarize(Max_D = max(abs(Cohens_D))) %>%
    arrange(desc(Max_D)) %>%
    pull(Genus)
  
  plot_data$Genus <- factor(plot_data$Genus, levels = genus_order)
  
  panel_D <- ggplot(plot_data, aes(x = Cohens_D, y = Genus, fill = Comparison)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("Strobilation vs Aposymbiotic" = "#029E73",
                                 "Strobilation vs Mutant" = "#DE8F05"),
                      labels = c("Strobilation vs\nAposymbiotic", 
                                 "Strobilation vs\nMutant")) +
    labs(x = "Cohen's d Effect Size", y = "", fill = "") +
    theme_classic(base_size = 16) +
    theme(axis.text.y = element_text(size = 13, face = "italic"),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          legend.position = "bottom",
          legend.text = element_text(size = 12))
  
  return(panel_D)
}

panel_D <- create_cohens_d_plot(wilcox_strob_vs_apo, wilcox_strob_vs_mut)

if (!is.null(panel_D)) {
  ggsave(file.path(output_dir, "Panel_D_Cohens_D.png"), panel_D, width = 8, height = 7, dpi = 300)
  
}

# ============================================================================
# LOG2FC PLOT 
# ============================================================================

cat("=== LOG2FC EFFECT SIZES ===\n")

create_log2fc_plot <- function(strob_vs_apo, strob_vs_mut) {
  
  sig_apo <- filter_biomarkers(strob_vs_apo)
  sig_mut <- filter_biomarkers(strob_vs_mut)
  
  all_biomarkers <- unique(c(sig_apo$ASV, sig_mut$ASV))
  
  if (length(all_biomarkers) == 0) {
    cat("No significant biomarkers for Log2FC plot\n")
    return(NULL)
  }
  
  cat("Creating Log2FC plot for", length(all_biomarkers), "biomarkers\n")
  
  plot_data <- data.frame()
  
  for (asv in all_biomarkers) {
    genus_name <- NA
    if (asv %in% strob_vs_apo$ASV) {
      genus_name <- strob_vs_apo$Genus[strob_vs_apo$ASV == asv][1]
    } else if (asv %in% strob_vs_mut$ASV) {
      genus_name <- strob_vs_mut$Genus[strob_vs_mut$ASV == asv][1]
    }
    
    genus_name <- ifelse(is.na(genus_name) | genus_name == "", "Unclassified", genus_name)
    
    log2fc_apo <- ifelse(asv %in% strob_vs_apo$ASV, 
                         strob_vs_apo$Log2FC[strob_vs_apo$ASV == asv][1], NA)
    log2fc_mut <- ifelse(asv %in% strob_vs_mut$ASV, 
                         strob_vs_mut$Log2FC[strob_vs_mut$ASV == asv][1], NA)
    
    if (!is.na(log2fc_apo)) {
      plot_data <- rbind(plot_data, data.frame(
        Genus = genus_name,
        Log2FC = log2fc_apo,
        Comparison = "Strobilation vs Aposymbiotic",
        ASV = asv
      ))
    }
    if (!is.na(log2fc_mut)) {
      plot_data <- rbind(plot_data, data.frame(
        Genus = genus_name,
        Log2FC = log2fc_mut,
        Comparison = "Strobilation vs Mutant",
        ASV = asv
      ))
    }
  }
  
  # Order by absolute Log2FC
  genus_order <- plot_data %>%
    group_by(Genus) %>%
    summarize(Max_Log2FC = max(abs(Log2FC))) %>%
    arrange(desc(Max_Log2FC)) %>%
    pull(Genus)
  
  plot_data$Genus <- factor(plot_data$Genus, levels = genus_order)
  
  panel_log2fc <- ggplot(plot_data, aes(x = Log2FC, y = Genus, fill = Comparison)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("Strobilation vs Aposymbiotic" = "#029E73",
                                 "Strobilation vs Mutant" = "#DE8F05"),
                      labels = c("Strobilation vs\nAposymbiotic", 
                                 "Strobilation vs\nMutant")) +
    labs(x = "Log2 Fold-Change", y = "", fill = "") +
    theme_classic(base_size = 16) +
    theme(axis.text.y = element_text(size = 13, face = "italic"),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          legend.position = "bottom",
          legend.text = element_text(size = 12))
  
  return(panel_log2fc)
}

panel_log2fc <- create_log2fc_plot(wilcox_strob_vs_apo, wilcox_strob_vs_mut)

if (!is.null(panel_log2fc)) {
  ggsave(file.path(output_dir, "Panel_Log2FC_Effect_Sizes.png"), 
         panel_log2fc, width = 8, height = 7, dpi = 300)
 
}

# ============================================================================
# COMBINED FIGURE
# ============================================================================

cat("=== CREATING COMBINED FIGURE ===\n")

if (!is.null(panel_B) && !is.null(panel_D)) {
  
  top_row <- ggarrange(panel_A, panel_D, 
                       labels = c("A", "D"), 
                       font.label = list(size = 18),
                       ncol = 2, nrow = 1,
                       widths = c(1, 1))
  
  bottom_row <- ggarrange(panel_B, 
                          labels = c("B"),
                          font.label = list(size = 18),
                          ncol = 1)
  
  combined_figure <- ggarrange(top_row, bottom_row,
                               ncol = 1, nrow = 2,
                               heights = c(1, 1.2))
  
  ggsave(file.path(output_dir, "Figure4_Combined.png"), 
         combined_figure, width = 16, height = 14, dpi = 300)
  
  }

# ============================================================================
# SUMMARY STATISTICS TABLE
# ============================================================================

cat("=== CREATING SUMMARY TABLE ===\n")

create_summary_table <- function(sig_results, comparison_name) {
  if (nrow(sig_results) == 0) {
    return(data.frame())
  }
  
  summary <- sig_results %>%
    arrange(desc(abs(Cohens_D))) %>%
    select(Genus, Family, Cohens_D, Log2FC, P_adjusted, 
           Mean_Group1, Mean_Group2) %>%
    mutate(
      Comparison = comparison_name,
      Direction = ifelse(Cohens_D > 0, "Enriched_Group1", "Depleted_Group1")
    ) %>%
    head(10)
  
  return(summary)
}

summary_strob_vs_mut <- create_summary_table(sig_strob_vs_mut, "Strobilation_vs_Mutant")
summary_strob_vs_apo <- create_summary_table(sig_strob_vs_apo, "Strobilation_vs_Aposymbiotic")
summary_mut_vs_apo <- create_summary_table(sig_mut_vs_apo, "Mutant_vs_Aposymbiotic")

all_summaries <- rbind(summary_strob_vs_mut, summary_strob_vs_apo, summary_mut_vs_apo)

write.csv(all_summaries, file.path(output_dir, "Top_Biomarkers_Summary.csv"), row.names = FALSE)

