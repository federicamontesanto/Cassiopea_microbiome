################################################################################
# Cassiopea microbiome project
# Relative-abundance summary tables for major bacterial families and genera across algae and polyp treatment groups.
################################################################################

required_packages <- c(
  "phyloseq",
  "dplyr",
  "tidyr",
  "tibble"
)


# ---------------------------------------------------------------------------
# Path configuration
# ---------------------------------------------------------------------------
project_root <- "."
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "04_taxonomic_composition")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
phyloseq_bact_file <- file.path(data_dir, "bacteria_phyloseq.rds")

# --- 1. Load and clean ---
ps <- readRDS(phyloseq_bact_file)
sample_data(ps)$Treatment <- trimws(as.character(sample_data(ps)$Treatment))

ps <- subset_samples(ps, !Treatment %in% c("Blank_1", "Blank_2", "Blank_3",
                                           "Water", "Polyp_antifungal", "Algae_antifungal"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

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
sample_data(ps)$Group <- name_map[sample_data(ps)$Treatment]

# --- 2. Confirm Group ---
cat("Group table:\n")
print(table(sample_data(ps)$Group, useNA = "ifany"))

extract_rel_abund <- function(ps_obj, rank, groups, top_n = 12) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("RELATIVE ABUNDANCE —", rank, "level\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  group_vec <- as.character(sample_data(ps_obj)$Group)
  keep      <- group_vec %in% groups
  ps_sub    <- prune_samples(keep, ps_obj)
  ps_sub    <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  cat("Samples in subset:", nsamples(ps_sub), "\n")
  
  ps_glom <- tax_glom(ps_sub, taxrank = rank, NArm = TRUE)
  ps_rel  <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
  
  # Check orientation of otu_table
  cat("taxa_are_rows:", taxa_are_rows(ps_rel), "\n")
  
  # Extract otu table correctly 
  otu_mat <- as.matrix(otu_table(ps_rel))
  if (taxa_are_rows(ps_rel)) {
    otu_mat <- t(otu_mat)  # make samples = rows
  }
  # Now rows = samples, cols = ASVs
  cat("OTU matrix dimensions (samples x ASVs):", dim(otu_mat), "\n")
  cat("Sample rowname examples:", paste(head(rownames(otu_mat), 3), collapse = ", "), "\n")
  
  otu_df          <- as.data.frame(otu_mat)
  otu_df$SampleID <- rownames(otu_df)
  
  tax_df  <- as.data.frame(tax_table(ps_rel), stringsAsFactors = FALSE)
  meta_df <- data.frame(sample_data(ps_rel),  stringsAsFactors = FALSE)
  
  # Standardize SampleIDs — replace hyphens with dots in both
  otu_df$SampleID  <- gsub("-", ".", otu_df$SampleID)
  meta_df$SampleID <- gsub("-", ".", rownames(meta_df))
  
  cat("OTU SampleID examples:", paste(head(otu_df$SampleID, 3), collapse = ", "), "\n")
  cat("Meta SampleID examples:", paste(head(meta_df$SampleID, 3), collapse = ", "), "\n")
  cat("Matching SampleIDs:", sum(otu_df$SampleID %in% meta_df$SampleID), "\n")
  
  df_long <- otu_df %>%
    tidyr::pivot_longer(-SampleID,
                        names_to  = "ASV",
                        values_to = "RelAbund") %>%
    left_join(tax_df %>%
                tibble::rownames_to_column("ASV") %>%
                dplyr::select(ASV, dplyr::all_of(rank)),
              by = "ASV") %>%
    left_join(meta_df %>%
                dplyr::select(SampleID, Group),
              by = "SampleID") %>%
    mutate(Group = factor(Group, levels = groups))
  
  cat("NA Groups after join:", sum(is.na(df_long$Group)), "\n")
  cat("Non-NA Groups:", sum(!is.na(df_long$Group)), "\n")
  
  top_taxa <- df_long %>%
    filter(!is.na(.data[[rank]]), !is.na(Group)) %>%
    group_by(.data[[rank]]) %>%
    summarise(overall_mean = mean(RelAbund, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(overall_mean)) %>%
    slice_head(n = top_n) %>%
    pull(.data[[rank]])
  
  summary_df <- df_long %>%
    filter(.data[[rank]] %in% top_taxa, !is.na(Group)) %>%
    group_by(Group, .data[[rank]]) %>%
    summarise(
      n      = n(),
      mean   = round(mean(RelAbund, na.rm = TRUE), 2),
      se     = round(sd(RelAbund,   na.rm = TRUE) / sqrt(n()), 2),
      median = round(median(RelAbund, na.rm = TRUE), 2),
      min    = round(min(RelAbund,  na.rm = TRUE), 2),
      max    = round(max(RelAbund,  na.rm = TRUE), 2),
      .groups = "drop"
    ) %>%
    mutate(mean_se = paste0(mean, " ± ", se),
           range   = paste0(min,  "–",   max)) %>%
    arrange(Group, desc(mean))
  
  for (grp in groups) {
    cat("\n---", grp, "---\n")
    sub <- summary_df %>%
      filter(Group == grp) %>%
      dplyr::select(.data[[rank]], n, mean_se, median, range)
    print(sub, n = top_n)
  }
  
  return(invisible(summary_df))
}

# Run
cat("\n\n########## ALGAL CULTURE TREATMENTS ##########\n")
algae_family <- extract_rel_abund(ps, "Family", algae_groups)
algae_genus  <- extract_rel_abund(ps, "Genus",  algae_groups)

cat("\n\n########## POLYP TREATMENTS ##########\n")
polyp_family <- extract_rel_abund(ps, "Family", polyp_groups)
polyp_genus  <- extract_rel_abund(ps, "Genus",  polyp_groups)

write.csv(algae_family, file.path(output_dir, "relabund_algae_family.csv"), row.names = FALSE)
write.csv(algae_genus,  file.path(output_dir, "relabund_algae_genus.csv"),  row.names = FALSE)
write.csv(polyp_family, file.path(output_dir, "relabund_polyp_family.csv"), row.names = FALSE)
write.csv(polyp_genus,  file.path(output_dir, "relabund_polyp_genus.csv"),  row.names = FALSE)

cat("\nDone.\n")
# Proportion of reads from shared ASVs across algal cultures
library(phyloseq)

# Get algae subset
keep_algae <- as.character(sample_data(ps)$Group) %in% algae_groups
ps_algae   <- prune_samples(keep_algae, ps)
ps_algae   <- prune_taxa(taxa_sums(ps_algae) > 0, ps_algae)

# Presence/absence per sample
otu_mat <- as.matrix(otu_table(ps_algae))
if (taxa_are_rows(ps_algae)) otu_mat <- t(otu_mat)

# ASVs present in ALL treatment groups
group_vec   <- as.character(sample_data(ps_algae)$Group)
pa_mat      <- (otu_mat > 0) * 1

# For each ASV check if present in at least one sample per group
asv_in_group <- sapply(algae_groups, function(g) {
  idx <- group_vec == g
  colSums(pa_mat[idx, , drop = FALSE]) > 0
})

shared_all  <- rowSums(asv_in_group) == length(algae_groups)
cat("ASVs present in ALL algal groups:", sum(shared_all), "\n")

# Proportion of total reads from shared ASVs
total_reads  <- sum(otu_mat)
shared_reads <- sum(otu_mat[, shared_all])
cat("Proportion of reads from shared ASVs:",
    round(shared_reads / total_reads * 100, 1), "%\n")

# Same for polyp groups
keep_polyp <- as.character(sample_data(ps)$Group) %in% polyp_groups
ps_polyp   <- prune_samples(keep_polyp, ps)
ps_polyp   <- prune_taxa(taxa_sums(ps_polyp) > 0, ps_polyp)

otu_polyp   <- as.matrix(otu_table(ps_polyp))
if (taxa_are_rows(ps_polyp)) otu_polyp <- t(otu_polyp)

group_polyp <- as.character(sample_data(ps_polyp)$Group)
pa_polyp    <- (otu_polyp > 0) * 1

asv_in_polyp <- sapply(polyp_groups, function(g) {
  idx <- group_polyp == g
  colSums(pa_polyp[idx, , drop = FALSE]) > 0
})
# Supplementary table: number of ASVs per family and genus per treatment
library(phyloseq)
library(dplyr)

# --- Algae ---
keep_algae <- as.character(sample_data(ps)$Group) %in% algae_groups
ps_algae   <- prune_samples(keep_algae, ps)
ps_algae   <- prune_taxa(taxa_sums(ps_algae) > 0, ps_algae)

tax_algae  <- as.data.frame(tax_table(ps_algae), stringsAsFactors = FALSE)

asv_per_family_algae <- tax_algae %>%
  group_by(Family) %>%
  summarise(n_ASVs = n(), .groups = "drop") %>%
  arrange(desc(n_ASVs)) %>%
  mutate(Substrate = "Algal cultures")

asv_per_genus_algae <- tax_algae %>%
  group_by(Family, Genus) %>%
  summarise(n_ASVs = n(), .groups = "drop") %>%
  arrange(Family, desc(n_ASVs)) %>%
  mutate(Substrate = "Algal cultures")

# --- Polyps ---
keep_polyp <- as.character(sample_data(ps)$Group) %in% polyp_groups
ps_polyp   <- prune_samples(keep_polyp, ps)
ps_polyp   <- prune_taxa(taxa_sums(ps_polyp) > 0, ps_polyp)

tax_polyp  <- as.data.frame(tax_table(ps_polyp), stringsAsFactors = FALSE)

asv_per_family_polyp <- tax_polyp %>%
  group_by(Family) %>%
  summarise(n_ASVs = n(), .groups = "drop") %>%
  arrange(desc(n_ASVs)) %>%
  mutate(Substrate = "Polyp hosts")

asv_per_genus_polyp <- tax_polyp %>%
  group_by(Family, Genus) %>%
  summarise(n_ASVs = n(), .groups = "drop") %>%
  arrange(Family, desc(n_ASVs)) %>%
  mutate(Substrate = "Polyp hosts")

# --- Combine and save ---
asv_family_combined <- bind_rows(asv_per_family_algae, asv_per_family_polyp)
asv_genus_combined  <- bind_rows(asv_per_genus_algae,  asv_per_genus_polyp)

write.csv(asv_family_combined,
          file.path(base_dir, "SupTable_ASVs_per_family.csv"),
          row.names = FALSE)
write.csv(asv_genus_combined,
          file.path(base_dir, "SupTable_ASVs_per_genus.csv"),
          row.names = FALSE)

# Print summary for key taxa discussed in Results
key_families <- c("Thalassospiraceae", "Cyclobacteriaceae", "Balneolaceae",
                  "Comamonadaceae", "Marinobacteraceae", "Rhodobacteraceae",
                  "Vibrionaceae", "Saprospiraceae", "Flavobacteriaceae",
                  "Oceanospirillaceae", "Phycisphaeraceae")

cat("ASV counts for key families discussed in Results:\n")
cat("\n--- Algal cultures ---\n")
print(asv_per_family_algae %>%
        filter(Family %in% key_families) %>%
        select(Family, n_ASVs))

cat("\n--- Polyp hosts ---\n")
print(asv_per_family_polyp %>%
        filter(Family %in% key_families) %>%
        select(Family, n_ASVs))

cat("\nSupplementary tables saved to:", base_dir, "\n")
shared_polyp  <- rowSums(asv_in_polyp) == length(polyp_groups)
cat("\nASVs present in ALL polyp groups:", sum(shared_polyp), "\n")

total_polyp  <- sum(otu_polyp)
shared_reads_polyp <- sum(otu_polyp[, shared_polyp])
cat("Proportion of reads from shared ASVs:",
    round(shared_reads_polyp / total_polyp * 100, 1), "%\n")