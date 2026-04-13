################################################################################
# Cassiopea microbiome project
# Bray-Curtis PCoA workflow with overall and pairwise PERMANOVA testing
################################################################################

# ==============================================================================
# PCoA Analysis: Bacteria 
# ==============================================================================

required_packages <- c(
  "phyloseq", "ggplot2", "vegan", "dplyr", "tibble",
  "ggpubr", "rstatix", "stringr", "tidyr", "scales", "grid",
  "pairwiseAdonis"
)

# ---------------------------------------------------------------------------
# Path configuration
# ---------------------------------------------------------------------------
project_root <- "."
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "03_beta_diversity/pcoa_pairwise_permanova")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
phyloseq_bact_file <- file.path(data_dir, "bacteria_phyloseq.rds")
grouping_column <- "Treatment"
cat("Output folder:", output_dir, "\n")

# Rarefaction depth & read threshold
rarefaction_depth_bacteria   <- 6000
min_reads_threshold_bacteria <- 5000

# --- Treatment Name Mapping  ---
treatment_name_mapping <- c(
  "Algae_Smic_KB8"       = "Algae_Native",
  "Algae_Bmin_SSB01"     = "Algae_Control",      
  "Algae_antibiotic"     = "Algae_Antibiotic",     
  "Algae_mutant"         = "Algae_Mutant",       
  # POLYPS
  "Polyp_Aposymbiotic"   = "Apo",
  "Polyp_Smic_KB8"       = "Polyp_Native",
  "Polyp_Bmin_SSB01"     = "Polyp_Control",
  "Polyp_antibiotic"     = "Polyp_Antibiotic",
  "Polyp_Mutant"         = "Polyp_Mutant"
  )

apply_treatment_mapping <- function(x) {
  sapply(x, function(name) if (name %in% names(treatment_name_mapping)) treatment_name_mapping[name] else name,
         USE.NAMES = FALSE)
}

# --- Color & Shape Mapping ---
clean_treatment_order <- c(
  "Apo", "Polyp_Native", "Polyp_Control", "Polyp_Antibiotic", "Polyp_Mutant",
  "Algae_Native", "Algae_Control", "Algae_Antibiotic", "Algae_Mutant"
)

publication_colors <- c(
  "Apo"              = "#762A83",
  "Polyp_Native"     = "#D73027",
  "Polyp_Control"    = "#F46D43",
  "Polyp_Antibiotic" = "#FDAE61",
  "Polyp_Mutant"     = "#E6F598",
  "Algae_Native"     = "#1B7837",
  "Algae_Control"    = "#5AAE61",
  "Algae_Antibiotic" = "#7FBC41",
  "Algae_Mutant"     = "#C2E681"
)

publication_shapes <- c(
  "Apo"         = 8,
  "Polyp_Native"  = 18,
  "Polyp_Control"  = 17,
  "Polyp_Antibiotic" = 16,
  "Polyp_Mutant"   = 15,
  "Algae_Native"  = 18,
  "Algae_Control"  = 17,
  "Algae_Antibiotic" = 16,
  "Algae_Mutant"   = 15
)

# Helper function for pairwise PERMANOVA with clear letter output ---
perform_pairwise_permanova <- function(dist_mat, group_vector, p_adjust_method = "fdr") {
  
  # Ensure we have a factor
  if (!is.factor(group_vector)) {
    group_vector <- as.factor(group_vector)
  }
  
  # Create a simple data frame with matching names
  sample_names <- labels(dist_mat)
  pairwise_data <- data.frame(
    Group = group_vector[sample_names],
    row.names = sample_names,
    stringsAsFactors = TRUE
  )
  
  # Double-check alignment
  if (nrow(pairwise_data) != attr(dist_mat, "Size")) {
    cat("Sample mismatch detected, realigning...\n")
    common_samples <- intersect(sample_names, names(group_vector))
    dist_mat <- as.dist(as.matrix(dist_mat)[common_samples, common_samples])
    pairwise_data <- data.frame(
      Group = group_vector[common_samples],
      row.names = common_samples,
      stringsAsFactors = TRUE
    )
  }
  
  cat("    Testing", nlevels(pairwise_data$Group), "groups with", nrow(pairwise_data), "samples\n")
  
  # Pairwise PERMANOVA
  pairwise_result <- tryCatch({
    pairwiseAdonis::pairwise.adonis(dist_mat, pairwise_data$Group, 
                                    p.adjust.m = p_adjust_method, 
                                    perm = 999)
  }, error = function(e) {
    cat("Pairwise PERMANOVA failed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(pairwise_result)) return(NULL)
  
  # Extract p-values and create results
  pairs <- as.character(pairwise_result$pairs)
  p_values <- pairwise_result$p.adjusted
  
  # Create results data frame
  results_df <- data.frame(
    Comparison = pairs,
    R2 = round(pairwise_result$R2, 3),
    F_value = round(pairwise_result$F.Model, 2),
    p_value = round(pairwise_result$p.value, 4),
    p_adjusted = round(p_values, 4),
    Significant = ifelse(p_values < 0.001, "***", 
                         ifelse(p_values < 0.01, "**",
                                ifelse(p_values < 0.05, "*", "ns"))),
    stringsAsFactors = FALSE
  )
  
  # Generate SIMPLE CLD using rcompanion
  group_levels <- levels(pairwise_data$Group)
  n_groups <- length(group_levels)
  
  # Create a full comparison table for cldList
  comparison_table <- data.frame(
    Group1 = character(),
    Group2 = character(),
    p.adj = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(pairwise_result)) {
    pair_names <- strsplit(as.character(pairwise_result$pairs[i]), "_vs_")[[1]]
    if (length(pair_names) == 2 && all(pair_names %in% group_levels)) {
      comparison_table <- rbind(comparison_table, data.frame(
        Group1 = pair_names[1],
        Group2 = pair_names[2],
        p.adj = pairwise_result$p.adjusted[i],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  cat("  Generating compact letter display...\n")
  
  # Check for NA values in p.adj (for debugging)
  if (any(is.na(comparison_table$p.adj))) {
    n_na <- sum(is.na(comparison_table$p.adj))
    }
  
  # CLD GENERATION
  # Build a similarity matrix (1 = not significantly different, 0 = significantly different)
  sim_matrix <- matrix(1, nrow = n_groups, ncol = n_groups)
  rownames(sim_matrix) <- colnames(sim_matrix) <- group_levels
  
  # Fill in the similarity matrix based on p-values
  for (i in 1:nrow(comparison_table)) {
    g1 <- comparison_table$Group1[i]
    g2 <- comparison_table$Group2[i]
    
    # Handle NA p-values explicitly (treat as not significant)
    p_val <- comparison_table$p.adj[i]
    
    # Use isTRUE for safer NA handling
    if (isTRUE(p_val < 0.05)) {
      # Significantly different - mark as 0
      sim_matrix[g1, g2] <- 0
      sim_matrix[g2, g1] <- 0
    }
    # If p_val is NA or >= 0.05, leave as 1 (not significantly different)
  }
  

  letter_list <- vector("list", n_groups)
  names(letter_list) <- group_levels
  
  # Assign letters: groups that are NOT significantly different share letters
  current_letter <- 1  # Start with 'a'
  assigned <- rep(FALSE, n_groups)
  
  for (i in 1:n_groups) {
    if (!assigned[i]) {
      # Find all groups similar to this one (including itself)
      similar_to_i <- which(sim_matrix[i, ] == 1)
      
      # Assign current letter to all similar groups
      for (j in similar_to_i) {
        if (!assigned[j]) {
          letter_list[[j]] <- c(letter_list[[j]], letters[current_letter])
          assigned[j] <- TRUE
        }
      }
      current_letter <- current_letter + 1
    }
  }
  
  # Convert to compact format (combine letters)
  cld_letters <- sapply(letter_list, function(x) paste(sort(x), collapse = ""))
  
    for (i in 1:n_groups) {
    if (cld_letters[i] == "") {
      cld_letters[i] <- letters[i]
    }
  }
  
  # Create final CLD dataframe
  cld <- data.frame(
    Group = group_levels,
    Letter = cld_letters,
    stringsAsFactors = FALSE
  )
  
  cat("CLD generated successfully\n")
  
  return(list(
    pairwise_results = results_df,
    cld = cld,
    comparison_table = comparison_table
  ))
}

# --- Function to print pairwise results with letter display ---
print_pairwise_results <- function(pairwise_obj, method_name) {
  if (is.null(pairwise_obj)) return()
  
  cat(sprintf("\n PAIRWISE COMPARISONS - %s\n", method_name))
  cat(paste(rep("=", 90), collapse = ""), "\n")
  
  results <- pairwise_obj$pairwise_results
  print(results, row.names = FALSE)
  
  if (!is.null(pairwise_obj$cld)) {
    cat("\n📋 COMPACT LETTER DISPLAY (for manual plot addition):\n")
    cat("   Groups sharing letters are NOT significantly different (p.adj > 0.05)\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    cld_df <- pairwise_obj$cld
    # Simple print - just Group and Letter columns
    print(cld_df, row.names = FALSE)
    
    cat(paste(rep("-", 50), collapse = ""), "\n")
  }
  
  cat(paste(rep("=", 90), collapse = ""), "\n\n")
}

# --- Rarefaction helper ---
perform_rarefaction_pcoa <- function(ps, depth, dataset_name) {
  cat("  Performing rarefaction to", depth, "reads per sample...\n")
  sample_sums_vec <- sample_sums(ps)
  samples_above_depth <- sum(sample_sums_vec >= depth)
  cat("  Samples with ≥", depth, "reads:", samples_above_depth, "out of", nsamples(ps), "\n")
  if (samples_above_depth < 5) {
    cat("Too few samples (", samples_above_depth, ") at depth; skipping rarefaction for", dataset_name, "\n")
    return(ps)
  }
  set.seed(12345)
  ps_rare <- tryCatch(
    rarefy_even_depth(ps, sample.size = depth, rngseed = 12345, replace = FALSE, trimOTUs = TRUE, verbose = TRUE),
    error = function(e) { cat("Rarefaction failed:", e$message, "\n"); ps }
  )
  cat("Rarefaction completed. Samples retained:", nsamples(ps_rare), "\n")
  ps_rare
}

# --- Distance helper ---
calculate_distance_safe <- function(ps_obj, method) {
  cat("    Computing", method, "distance...")
  if (method == "bray") {
    otu_mat <- as(otu_table(ps_obj), "matrix"); if (taxa_are_rows(ps_obj)) otu_mat <- t(otu_mat)
    res <- vegan::vegdist(otu_mat, method = "bray")
  } else if (method %in% c("unifrac", "wunifrac")) {
    weighted_flag <- (method == "wunifrac")
    res <- phyloseq::UniFrac(ps_obj, weighted = weighted_flag)
  } else stop("Unknown distance method: ", method)
  cat(" Done.\n"); res
}

# --- Core PCoA function ---
run_enhanced_pcoa_bacteria <- function(phylo_file, min_reads, rare_depth, out_dir) {
  cat("\n=== PROCESSING BACTERIA ===\n")
  if (!file.exists(phylo_file)) { cat("File not found:", phylo_file, "\n"); return(FALSE) }
  
  ps <- readRDS(phylo_file)
  cat("Loaded", nsamples(ps), "samples,", ntaxa(ps), "ASVs\n")
  
    # Read threshold
  ps <- prune_samples(sample_sums(ps) >= min_reads, ps)
  cat("After filtering <", min_reads, "reads:", nsamples(ps), "samples remain\n")
  if (nsamples(ps) < 5) { cat("Too few samples remaining\n"); return(FALSE) }
  
  # Rarefaction
  ps_rare <- perform_rarefaction_pcoa(ps, rare_depth, "Bacteria")
  
  # Prepare labels & palettes (clean names)
  sample_data(ps_rare)[[grouping_column]] <- factor(sample_data(ps_rare)[[grouping_column]])
  orig_levels  <- levels(sample_data(ps_rare)[[grouping_column]])
  clean_levels <- apply_treatment_mapping(orig_levels)
  
  plot_order <- clean_treatment_order[clean_treatment_order %in% clean_levels]
  cols <- publication_colors[plot_order]
  shps <- publication_shapes[plot_order]
  
  # Distance methods to run
  methods_to_try <- c("bray", "unifrac", "wunifrac")
  
  for (method in methods_to_try) {
    cat(sprintf("\n--- %s Distance ---\n", toupper(method)))
    # Need tree for UniFrac
    if (method != "bray" && is.null(phy_tree(ps_rare, errorIfNULL = FALSE))) {
      cat("No phylogeny available: skipping", method, "\n")
      next
    }
    
    # Distance
    dist_mat <- tryCatch(calculate_distance_safe(ps_rare, method),
                         error = function(e) { cat("Dstance error:", e$message, "\n"); NULL })
    if (is.null(dist_mat)) next
    
    # PCoA
    pcoa <- tryCatch(cmdscale(dist_mat, k = min(nsamples(ps_rare) - 1, 10), eig = TRUE),
                     error = function(e) { cat("PcoA error:", e$message, "\n"); NULL })
    if (is.null(pcoa)) next
    
    eig <- pcoa$eig
    var_expl <- eig[eig > 0] / sum(eig[eig > 0]) * 100
    v1 <- round(var_expl[1], 1); v2 <- round(var_expl[2], 1)
    
    
    # Build plotting data
    coords <- as.data.frame(pcoa$points[, 1:2]); names(coords) <- c("Axis1", "Axis2")
    meta <- data.frame(sample_data(ps_rare))
    plot_df <- cbind(coords, meta[rownames(coords), , drop = FALSE])
    
    # Keep original grouping for PERMANOVA
    plot_df$Group_orig <- as.factor(plot_df[[grouping_column]])
    
    # Map to clean names for plotting
    plot_df$Group_clean <- factor(
      apply_treatment_mapping(as.character(plot_df[[grouping_column]])),
      levels = plot_order
    )
    
    # PERMANOVA 
    pm <- tryCatch(
      vegan::adonis2(dist_mat ~ Group_orig, data = plot_df, permutations = 999),
      error = function(e) { cat("PERMANOVA error:", e$message, "\n"); NULL }
    )
    if (!is.null(pm)) {
      r2 <- round(pm$R2[1], 3); pval <- pm[["Pr(>F)"]][1]
      cat(sprintf("R²=%.3f, p=%s\n", r2, format.pval(pval, digits = 3)))
      permanova_text <- paste0("PERMANOVA R² = ", r2, ", p ",
                               ifelse(pval < 0.001, "< 0.001", paste("=", format.pval(pval, digits = 3))))
    } else {
      permanova_text <- "PERMANOVA failed"
    }
    
    # *** PAIRWISE PERMANOVA ***
    group_vec <- setNames(plot_df$Group_orig, rownames(plot_df))
    
    pairwise_obj <- perform_pairwise_permanova(
      dist_mat = dist_mat, 
      group_vector = group_vec,
      p_adjust_method = "fdr"
    )
    print_pairwise_results(pairwise_obj, paste("Bacteria", toupper(method)))
    
    # Label for method
    lab_method <- switch(method,
                         bray = "Bray-Curtis",
                         unifrac = "Unweighted UniFrac",
                         wunifrac = "Weighted UniFrac")
    
    # Plot WITHOUT CLD letters (for manual addition)
    p <- ggplot(plot_df, aes(Axis1, Axis2, color = Group_clean, shape = Group_clean)) +
      geom_point(size = 6, alpha = 0.8) +
      stat_ellipse(aes(group = Group_clean),
                   type = "norm", level = 0.95,
                   linetype = "dashed", linewidth = 0.8,
                   show.legend = FALSE, alpha = 0.7) +
      scale_color_manual(values = cols, name = "Treatment", drop = TRUE) +
      scale_shape_manual(values = shps, name = "Treatment", drop = TRUE) +
      labs(
        x = paste0("Axis1 (", v1, "%)"),
        y = paste0("Axis2 (", v2, "%)"),
        caption = paste0(
          "Data rarefied to ", format(rare_depth, big.mark = ","), " reads/sample. ",
          permanova_text, ". ",
          "Pairwise comparisons with FDR correction. ",
          "Samples with <", format(min_reads, big.mark = ","), " reads excluded."
        )
      ) +
      theme_classic(base_size = 18) +
      theme(
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.text.x  = element_text(size = 28, color = "black"),
        axis.text.y  = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 28, color = "black"),
        axis.title.y = element_text(size = 28, color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 28, color = "black"),
        legend.text  = element_text(size = 28, color = "black"),
        legend.key.size = unit(1.2, "cm"),
        plot.caption = element_text(size = 24, hjust = 0, color = "black"),
        text = element_text(face = "plain", color = "black")
      ) +
      guides(color = guide_legend(override.aes = list(size = 5)))
    
    # Save
    file_safe <- gsub("[^A-Za-z0-9]", "", lab_method)
    png_out <- file.path(out_dir, paste0("Bacteria_", file_safe, "_PCoA_NoLetters.png"))
    pdf_out <- file.path(out_dir, paste0("Bacteria_", file_safe, "_PCoA_NoLetters.pdf"))
    ggsave(png_out, p, width = 12, height = 9, dpi = 300, bg = "white")
    ggsave(pdf_out, p, width = 12, height = 9, units = "in", bg = "white")
    
    
    # Save pairwise results and CLD to CSV
    if (!is.null(pairwise_obj)) {
      csv_out <- file.path(out_dir, paste0("Bacteria_", file_safe, "_PairwiseTests.csv"))
      write.csv(pairwise_obj$pairwise_results, csv_out, row.names = FALSE)
      
      
      if (!is.null(pairwise_obj$cld)) {
        cld_out <- file.path(out_dir, paste0("Bacteria_", file_safe, "_CompactLetters.csv"))
        write.csv(pairwise_obj$cld, cld_out, row.names = FALSE)
        
      }
    }
  }
  
  
# ==============================================================================
# POLYP-ONLY Beta Diversity Analysis
# ==============================================================================

polyp_treatment_order <- c("Apo", "Polyp_Smic", "Polyp_Bmin", "Polyp_Antib", "Polyp_Mut")

polyp_colors <- c(
  "Apo"              = "#762A83",
  "Polyp_Native"     = "#D73027", 
  "Polyp_Control"    = "#F46D43",
  "Polyp_Antibiotic" = "#FDAE61",
  "Polyp_Mutant"     = "#E6F598"
)

polyp_shapes <- c(
  "Apo"         = 8,
  "Polyp_Native"  = 18,
  "Polyp_Control"  = 17, 
  "Polyp_Antibiotic" = 16,
  "Polyp_Mutant"   = 15
)

run_polyp_only_pcoa <- function(phylo_file, min_reads, rare_depth, out_dir) {
  cat("\n=== PROCESSING POLYPS ONLY ===\n")
  if (!file.exists(phylo_file)) { cat("File not found:", phylo_file, "\n"); return(FALSE) }
  
  ps <- readRDS(phylo_file)
  cat("Loaded", nsamples(ps), "samples,", ntaxa(ps), "ASVs\n")
  
  # Rarefaction
  ps_rare <- perform_rarefaction_pcoa(ps, rare_depth, "Polyps_Only")
  
  # Prepare labels & palettes
  sample_data(ps_rare)[[grouping_column]] <- factor(sample_data(ps_rare)[[grouping_column]])
  orig_levels  <- levels(sample_data(ps_rare)[[grouping_column]])
  clean_levels <- apply_treatment_mapping(orig_levels)
  
  plot_order <- polyp_treatment_order[polyp_treatment_order %in% clean_levels]
  cols <- polyp_colors[plot_order]
  shps <- polyp_shapes[plot_order]
  cat("Polyp plot order:", paste(plot_order, collapse = " → "), "\n")
  
  methods_to_try <- c("bray", "unifrac", "wunifrac")
  
  for (method in methods_to_try) {
    cat(sprintf("\n--- %s Distance (Polyps Only) ---\n", toupper(method)))
    
    if (method != "bray" && is.null(phy_tree(ps_rare, errorIfNULL = FALSE))) {
      cat("No phylogeny available: skipping", method, "\n")
      next
    }
    
    dist_mat <- tryCatch(calculate_distance_safe(ps_rare, method),
                         error = function(e) { cat("Distance error:", e$message, "\n"); NULL })
    if (is.null(dist_mat)) next
    
    pcoa <- tryCatch(cmdscale(dist_mat, k = min(nsamples(ps_rare) - 1, 10), eig = TRUE),
                     error = function(e) { cat("PCoA error:", e$message, "\n"); NULL })
    if (is.null(pcoa)) next
    
    eig <- pcoa$eig
    var_expl <- eig[eig > 0] / sum(eig[eig > 0]) * 100
    v1 <- round(var_expl[1], 1); v2 <- round(var_expl[2], 1)
    cat(sprintf("Variance explained: Axis1=%.1f%%, Axis2=%.1f%%\n", v1, v2))
    
    coords <- as.data.frame(pcoa$points[, 1:2]); names(coords) <- c("Axis1", "Axis2")
    meta <- data.frame(sample_data(ps_rare))
    plot_df <- cbind(coords, meta[rownames(coords), , drop = FALSE])
    
    plot_df$Group_orig <- as.factor(plot_df[[grouping_column]])
    plot_df$Group_clean <- factor(
      apply_treatment_mapping(as.character(plot_df[[grouping_column]])),
      levels = plot_order
    )
    
    cat("Running overall PERMANOVA...")
    pm <- tryCatch(
      vegan::adonis2(dist_mat ~ Group_orig, data = plot_df, permutations = 999),
      error = function(e) { cat("PERMANOVA error:", e$message, "\n"); NULL }
    )
    if (!is.null(pm)) {
      r2 <- round(pm$R2[1], 3); pval <- pm[["Pr(>F)"]][1]
      cat(sprintf("R²=%.3f, p=%s\n", r2, format.pval(pval, digits = 3)))
      permanova_text <- paste0("PERMANOVA R² = ", r2, ", p ",
                               ifelse(pval < 0.001, "< 0.001", paste("=", format.pval(pval, digits = 3))))
    } else {
      permanova_text <- "PERMANOVA failed"
    }
    
    group_vec <- setNames(plot_df$Group_orig, rownames(plot_df))
    
    pairwise_obj <- perform_pairwise_permanova(
      dist_mat = dist_mat, 
      group_vector = group_vec,
      p_adjust_method = "fdr"
    )
    print_pairwise_results(pairwise_obj, paste("Polyps", toupper(method)))
    
    lab_method <- switch(method,
                         bray = "Bray-Curtis",
                         unifrac = "Unweighted UniFrac", 
                         wunifrac = "Weighted UniFrac")
    
    p <- ggplot(plot_df, aes(Axis1, Axis2, color = Group_clean, shape = Group_clean)) +
      geom_point(size = 6, alpha = 0.8) +
      stat_ellipse(aes(group = Group_clean),
                   type = "norm", level = 0.95,
                   linetype = "dashed", linewidth = 0.8,
                   show.legend = FALSE, alpha = 0.7) +
      scale_color_manual(values = cols, name = "Treatment", drop = TRUE) +
      scale_shape_manual(values = shps, name = "Treatment", drop = TRUE) +
      labs(
        x = paste0("Axis1 (", v1, "%)"),
        y = paste0("Axis2 (", v2, "%)"),
        title = paste0("Polyp Bacterial Communities - ", lab_method),
        caption = paste0(
          "Polyp samples only. Rarefied to ", format(rare_depth, big.mark = ","), " reads/sample. ",
          permanova_text, ". ",
          "Pairwise comparisons with FDR correction."
        )
      ) +
      theme_classic(base_size = 18) +
      theme(
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.text.x  = element_text(size = 24, color = "black"),
        axis.text.y  = element_text(size = 24, color = "black"),
        axis.title.x = element_text(size = 26, color = "black"),
        axis.title.y = element_text(size = 26, color = "black"),
        plot.title = element_text(size = 28, color = "black", hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 24, color = "black"),
        legend.text  = element_text(size = 22, color = "black"),
        legend.key.size = unit(1.2, "cm"),
        plot.caption = element_text(size = 20, hjust = 0, color = "black"),
        text = element_text(face = "plain", color = "black")
      ) +
      guides(color = guide_legend(override.aes = list(size = 5)))
    
    file_safe <- gsub("[^A-Za-z0-9]", "", lab_method)
    png_out <- file.path(out_dir, paste0("PolypOnly_", file_safe, "_PCoA_NoLetters.png"))
    pdf_out <- file.path(out_dir, paste0("PolypOnly_", file_safe, "_PCoA_NoLetters.pdf"))
    ggsave(png_out, p, width = 12, height = 9, dpi = 300, bg = "white")
    ggsave(pdf_out, p, width = 12, height = 9, units = "in", bg = "white")
    cat("Polyp plots saved (NO LETTERS):\n   PNG:", png_out, "\n   PDF:", pdf_out, "\n")
    
    if (!is.null(pairwise_obj)) {
      csv_out <- file.path(out_dir, paste0("PolypOnly_", file_safe, "_PairwiseTests.csv"))
      write.csv(pairwise_obj$pairwise_results, csv_out, row.names = FALSE)
      cat("Pairwise results saved:", csv_out, "\n")
      
      if (!is.null(pairwise_obj$cld)) {
        cld_out <- file.path(out_dir, paste0("PolypOnly_", file_safe, "_CompactLetters.csv"))
        write.csv(pairwise_obj$cld, cld_out, row.names = FALSE)
        
      }
    }
  }
   }

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

cat("\n STARTING ENHANCED PCoA (Bacteria - All Treatments)\n")
ok <- run_enhanced_pcoa_bacteria(
  phylo_file = phyloseq_bact_file,
  min_reads  = min_reads_threshold_bacteria,
  rare_depth = rarefaction_depth_bacteria,
  out_dir    = output_dir
)

