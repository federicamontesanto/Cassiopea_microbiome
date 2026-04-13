################################################################################
# Cassiopea microbiome project
# Alpha-diversity summary statistics for bacterial phyloseq data, including descriptive statistics and pairwise tests across Cassiopea treatment groups.
################################################################################

required_packages <- c(
  "phyloseq",
  "dplyr",
  "tidyr",
  "rstatix",
  "multcompView",
  "tibble"
)


# ---------------------------------------------------------------------------
# Path configuration
# ---------------------------------------------------------------------------
project_root <- "."
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "02_alpha_diversity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
phyloseq_bact_file <- file.path(data_dir, "bacteria_phyloseq.rds")

# --- 1. Load data ---
ps_relaxed <- readRDS(phyloseq_bact_file)

# Fix trailing spaces in Treatment labels
sample_data(ps_relaxed)$Treatment <- trimws(
  as.character(sample_data(ps_relaxed)$Treatment)
)

# --- 2. Filter and rarefy ---
ps_filtered <- prune_samples(sample_sums(ps_relaxed) >= 5000, ps_relaxed)
set.seed(12345)
ps_rare <- rarefy_even_depth(ps_filtered,
                             sample.size = 6000,
                             rngseed     = 12345,
                             replace     = FALSE,
                             trimOTUs    = TRUE,
                             verbose     = FALSE)

# --- 3. Calculate alpha diversity ---
alpha_div <- estimate_richness(ps_rare,
                               measures = c("Observed", "Fisher",
                                            "Shannon", "InvSimpson")) %>%
  tibble::rownames_to_column("SampleID")

# --- 4. Attach metadata ---
meta <- data.frame(
  SampleID  = sample_names(ps_rare),
  Treatment = as.character(sample_data(ps_rare)$Treatment),
  stringsAsFactors = FALSE
) %>%
  mutate(SampleID = gsub("-", ".", SampleID))

alpha_div <- alpha_div %>%
  left_join(meta, by = "SampleID") %>%
  filter(!Treatment %in% c("Blank_1", "Blank_2", "Blank_3",
                           "Polyp_antifungal", "Water"))

# --- 5. Define groups ---
algae_groups <- c("Algae_Smic_KB8",
                  "Algae_Bmin_SSB01",
                  "Algae_antibiotic",
                  "Algae_mutant")

polyp_groups <- c("Polyp_Aposymbiotic",
                  "Polyp_Smic_KB8",
                  "Polyp_Bmin_SSB01",
                  "Polyp_antibiotic",
                  "Polyp_Mutant")

# Friendly display names for output
name_map <- c(
  "Algae_Smic_KB8"     = "Algae_Native",
  "Algae_Bmin_SSB01"   = "Algae_Control",
  "Algae_antibiotic"   = "Algae_Antibiotic",
  "Algae_mutant"       = "Algae_Mutant",
  "Polyp_Aposymbiotic" = "Aposymbiotic",
  "Polyp_Smic_KB8"     = "Polyp_Native",
  "Polyp_Bmin_SSB01"   = "Polyp_Control",
  "Polyp_antibiotic"   = "Polyp_Antibiotic",
  "Polyp_Mutant"       = "Polyp_Mutant"
)

alpha_div <- alpha_div %>%
  mutate(Group = name_map[Treatment])

# --- 6. Summary statistics function ---
summarise_alpha <- function(df, groups, label) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SUMMARY STATISTICS —", label, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  df_sub <- df %>%
    filter(Treatment %in% groups) %>%
    mutate(Group = factor(name_map[Treatment],
                          levels = name_map[groups]))
  
  metrics <- c("Observed", "Fisher", "Shannon", "InvSimpson")
  
  for (metric in metrics) {
    cat("\n---", metric, "---\n")
    
    # Descriptive stats
    stats <- df_sub %>%
      group_by(Group) %>%
      summarise(
        n      = n(),
        mean   = round(mean(.data[[metric]], na.rm = TRUE), 2),
        se     = round(sd(.data[[metric]], na.rm = TRUE) / sqrt(n()), 2),
        median = round(median(.data[[metric]], na.rm = TRUE), 2),
        min    = round(min(.data[[metric]], na.rm = TRUE), 2),
        max    = round(max(.data[[metric]], na.rm = TRUE), 2),
        .groups = "drop"
      ) %>%
      mutate(mean_se = paste0(mean, " ± ", se),
             range   = paste0(min, "–", max))
    
    cat("Descriptive statistics (mean ± SE, median, range):\n")
    print(stats %>% select(Group, n, mean_se, median, range))
    
    # Normality check per group
    sw <- df_sub %>%
      group_by(Group) %>%
      summarise(p_shapiro = tryCatch(
        shapiro.test(.data[[metric]])$p.value,
        error = function(e) NA),
        .groups = "drop")
    normal <- all(sw$p_shapiro > 0.05, na.rm = TRUE)
    
    if (normal) {
      # ANOVA + Tukey
      fit   <- aov(as.formula(paste(metric, "~ Group")), data = df_sub)
      aov_s <- summary(fit)
      f_val <- round(aov_s[[1]][["F value"]][1], 3)
      p_val <- aov_s[[1]][["Pr(>F)"]][1]
      eta2  <- round(aov_s[[1]][["Sum Sq"]][1] /
                       sum(aov_s[[1]][["Sum Sq"]]), 3)
      
      cat(paste0("ANOVA: F = ", f_val,
                 ", p = ", format.pval(p_val, digits = 3),
                 ", η² = ", eta2, "\n"))
      
      if (p_val < 0.05) {
        tukey <- TukeyHSD(fit)
        tukey_df <- as.data.frame(tukey$Group) %>%
          tibble::rownames_to_column("Comparison") %>%
          filter(`p adj` < 0.05) %>%
          mutate(across(where(is.numeric), ~ round(., 4)))
        cat("Significant Tukey pairwise comparisons (p.adj < 0.05):\n")
        print(tukey_df %>% select(Comparison, diff, `p adj`))
      }
      
    } else {
      # Kruskal-Wallis + Dunn
      kw <- kruskal.test(as.formula(paste(metric, "~ Group")), data = df_sub)
      cat(paste0("Kruskal-Wallis: H = ", round(kw$statistic, 3),
                 ", df = ", kw$parameter,
                 ", p = ", format.pval(kw$p.value, digits = 3), "\n"))
      
      if (kw$p.value < 0.05) {
        dunn <- df_sub %>%
          dunn_test(as.formula(paste(metric, "~ Group")),
                    p.adjust.method = "BH") %>%
          filter(p.adj < 0.05) %>%
          mutate(across(where(is.numeric), ~ round(., 4)))
        cat("Significant Dunn pairwise comparisons (p.adj < 0.05):\n")
        print(dunn %>% select(group1, group2, statistic, p.adj))
      }
    }
  }
}

# --- 7. Run for algae and polyp groups ---
summarise_alpha(alpha_div, algae_groups, "ALGAL CULTURE TREATMENTS")
summarise_alpha(alpha_div, polyp_groups, "POLYP TREATMENTS")

# --- 8. Save full output to text file ---
output_file <- file.path(output_dir, "alpha_diversity_summary_for_results.txt")
sink(output_file)
summarise_alpha(alpha_div, algae_groups, "ALGAL CULTURE TREATMENTS")
summarise_alpha(alpha_div, polyp_groups, "POLYP TREATMENTS")
sink()

cat("\nSummary statistics saved to:\n", output_file, "\n")