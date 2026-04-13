################################################################################
# Cassiopea microbiome project
# Bacterial 16S amplicon workflow for Cassiopea xamachana, from raw FASTQ files through ASV inference, taxonomy assignment, decontamination, diversity analyses, differential abundance, and selected visualizations.
################################################################################

################################################################################
# Cassiopea xamachana Microbiome Analysis Pipeline
# Authors: Federica Montesanto et al.
# Description: 16S rRNA gene amplicon sequencing analysis of bacterial
#              communities associated with Cassiopea xamachana polyps and
#              Symbiodiniaceae algal cultures during strobilation.
# Data: Available on NCBI BioProject PRJNA1390206
# Last updated: 2026
################################################################################

# NOTES FOR REPRODUCIBILITY:
# R version: 4.3.x or later recommended
#   - Raw paired-end FASTQ files (from NCBI BioProject PRJNA1390206)
#   - SILVA database v138.1 (download from https://www.arb-silva.de/)

################################################################################
# 0. SETUP - Load required packages
################################################################################

required_packages <- c(
  "dada2",
  "phyloseq",
  "DESeq2",
  "decontam",
  "DECIPHER",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "vegan",
  "ape",
  "phangorn",
  "stringr",
  "RColorBrewer",
  "reshape2",
  "multcomp",
  "ggpubr",
  "writexl",
  "GUniFrac",
  "Biostrings",
  "scales",
  "ComplexUpset",
  "patchwork"
)


################################################################################
# 1. SET PATHS
################################################################################

# ---------------------------------------------------------------------------
# Path configuration
# ---------------------------------------------------------------------------
project_root <- "."
fastq_dir <- file.path(project_root, "data", "fastq", "bacteria")
metadata_path <- file.path(project_root, "metadata", "metadata_Cassiopea_bacteria.csv")
silva_train <- file.path(project_root, "reference", "silva_nr99_v138.1_train_set.fa.gz")
silva_species <- file.path(project_root, "reference", "silva_species_assignment_v138.1.fa.gz")
output_dir <- file.path(project_root, "results", "01_bacterial_16S_pipeline")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


################################################################################
# 2. DADA2 PIPELINE - Quality filtering, denoising, and ASV inference
################################################################################

# --- 2.1 File setup ---
path <- fastq_dir
filt_path <- file.path(path, "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)

# Primers used (515F Parada / 806R Apprill)
primers <- list(
  forward = "GTGYCAGCMGCCGCGGTAA",
  reverse = "GGACTACNVGGGTWTCTAAT"
)

# Get forward and reverse reads
files_F <- sort(list.files(path, pattern = "_R1.*fastq.gz$", full.names = TRUE))
files_R <- sort(list.files(path, pattern = "_R2.*fastq.gz$", full.names = TRUE))
cat("Found", length(files_F), "forward and", length(files_R), "reverse read files\n")

# Define filtered file paths
filtFs <- file.path(filt_path, basename(files_F))
filtRs <- file.path(filt_path, basename(files_R))

# --- 2.2 Quality inspection ---
pdf(file.path(output_dir, "quality_profiles.pdf"))
plotQualityProfile(files_F[1:min(2, length(files_F))])
plotQualityProfile(files_R[1:min(2, length(files_R))])
dev.off()

# --- 2.3 Filter and trim ---
# trimLeft removes primer bases: 19 bp forward, 20 bp reverse
# truncLen set based on quality profile inspection
out <- filterAndTrim(
  files_F, filtFs,
  files_R, filtRs,
  truncLen  = c(240, 200),
  maxN      = 0,
  maxEE     = c(2, 2),
  truncQ    = 2,
  trimLeft  = c(19, 20),
  rm.phix   = TRUE,
  compress  = TRUE,
  verbose   = TRUE,
  multithread = FALSE  # Set to TRUE on Mac/Linux
)
cat("Reads passing filter:\n")
print(out)

# --- 2.4 Learn error rates ---
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
saveRDS(errF, file.path(output_dir, "errF.rds"))
saveRDS(errR, file.path(output_dir, "errR.rds"))

# Plot error rates for inspection
pdf(file.path(output_dir, "error_rates.pdf"))
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

# --- 2.5 Sample inference (denoising) ---
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)
saveRDS(dadaFs, file.path(output_dir, "dadaFs.rds"))
saveRDS(dadaRs, file.path(output_dir, "dadaRs.rds"))

# --- 2.6 Merge paired reads ---
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
saveRDS(mergers, file.path(output_dir, "mergers.rds"))

# --- 2.7 Construct sequence table ---
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(output_dir, "seqtab.rds"))

# Inspect amplicon length distribution
cat("Amplicon length distribution:\n")
print(table(nchar(getSequences(seqtab))))

# --- 2.8 Remove chimeras ---
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = FALSE, verbose = TRUE)
saveRDS(seqtab.nochim, file.path(output_dir, "seqtab_nochim.rds"))

# --- 2.9 Track reads through pipeline ---
getN <- function(x) sum(getUniques(x))
track <- data.frame(
  Input      = out[, 1],
  Filtered   = out[, 2],
  Denoised_F = sapply(dadaFs, getN),
  Denoised_R = sapply(dadaRs, getN),
  Merged     = sapply(mergers, getN),
  Nonchim    = rowSums(seqtab.nochim)
)
write.csv(track, file.path(output_dir, "read_tracking.csv"))
cat("\nRead tracking summary:\n")
print(summary(track))

################################################################################
# 3. TAXONOMY ASSIGNMENT
################################################################################

# Assign taxonomy using SILVA v138.1
taxa <- assignTaxonomy(seqtab.nochim, silva_train, multithread = FALSE)

# Add species-level assignment in batches to avoid memory issues
batch_size  <- 500
total_asvs  <- nrow(taxa)
batches     <- split(seq_len(total_asvs),
                     ceiling(seq_along(seq_len(total_asvs)) / batch_size))
species_results <- vector("list", length(batches))

cat("Assigning species in", length(batches), "batches...\n")
for (i in seq_along(batches)) {
  idx <- batches[[i]]
  species_results[[i]] <- addSpecies(taxa[idx, , drop = FALSE], silva_species)
  cat("Batch", i, "of", length(batches), "complete\n")
}
taxa_species <- do.call(rbind, species_results)
saveRDS(taxa_species, file.path(output_dir, "taxa_species.rds"))

################################################################################
# 4. TAXONOMY FILTERING
################################################################################

# Remove non-bacterial, chloroplast, mitochondrial, and unclassified ASVs
taxa_filtered <- taxa_species[
  taxa_species[, "Kingdom"] == "Bacteria" &
    !grepl("Chloroplast", taxa_species[, "Order"],  ignore.case = TRUE) &
    !grepl("Chloroplast", taxa_species[, "Family"], ignore.case = TRUE) &
    !grepl("Mitochondria", taxa_species[, "Family"], ignore.case = TRUE) &
    !is.na(taxa_species[, "Phylum"]) & taxa_species[, "Phylum"] != "" &
    !is.na(taxa_species[, "Genus"])  & taxa_species[, "Genus"]  != "",
]

seqtab.filtered <- seqtab.nochim[, colnames(seqtab.nochim) %in% rownames(taxa_filtered)]
saveRDS(seqtab.filtered, file.path(output_dir, "seqtab_filtered.rds"))
saveRDS(taxa_filtered,   file.path(output_dir, "taxa_filtered.rds"))

cat("ASVs after taxonomy filtering:", ncol(seqtab.filtered), "\n")
cat("Reads after taxonomy filtering:", sum(seqtab.filtered), "\n")

################################################################################
# 5. DECONTAMINATION USING DECONTAM
################################################################################

# Load metadata
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, sep = ";",
                     check.names = FALSE)

# Clean sample IDs to match ASV table rownames
rownames(seqtab.filtered) <- gsub("_.*", "", rownames(seqtab.filtered))

# Match metadata to ASV table order
metadata$`16S_ID` <- trimws(metadata$`16S_ID`)
shared_ids <- intersect(rownames(seqtab.filtered), metadata$`16S_ID`)
asv_matched  <- seqtab.filtered[shared_ids, ]
meta_matched <- metadata[match(shared_ids, metadata$`16S_ID`), ]

# Identify extraction blanks
# Blanks must be labeled consistently in the Treatment column of the metadata
is.blank <- grepl("Blank", meta_matched$Treatment, ignore.case = TRUE)
cat("Number of blank samples:", sum(is.blank), "\n")
cat("Number of biological samples:", sum(!is.blank), "\n")

# Run decontam using prevalence method
# threshold = 0.5: ASVs more prevalent in blanks than samples are flagged
contamdf <- isContaminant(asv_matched,
                          method    = "prevalence",
                          neg       = is.blank,
                          threshold = 0.5)

cat("Total ASVs before decontam:", nrow(contamdf), "\n")
cat("Contaminant ASVs identified:", sum(contamdf$contaminant), "\n")

# Remove contaminant ASVs
seqtab.decontam <- asv_matched[, !contamdf$contaminant]
saveRDS(seqtab.decontam, file.path(output_dir, "seqtab_decontam.rds"))

cat("ASVs remaining after decontam:", ncol(seqtab.decontam), "\n")

# Update taxonomy to match decontaminated ASV table
taxa_decontam <- taxa_filtered[rownames(taxa_filtered) %in% colnames(seqtab.decontam), ]
saveRDS(taxa_decontam, file.path(output_dir, "taxa_decontam.rds"))
# Remove contaminant ASVs
seqtab.decontam <- asv_matched[, !contamdf$contaminant]
saveRDS(seqtab.decontam, file.path(output_dir, "seqtab_decontam.rds"))

cat("ASVs remaining after decontam:", ncol(seqtab.decontam), "\n")

# Update taxonomy to match decontaminated ASV table
taxa_decontam <- taxa_filtered[rownames(taxa_filtered) %in% colnames(seqtab.decontam), ]
saveRDS(taxa_decontam, file.path(output_dir, "taxa_decontam.rds"))

################################################################################
# 5B. CREATE PHYLOSEQ OBJECT
################################################################################

metadata_phyloseq <- metadata[metadata$`16S_ID` %in% rownames(seqtab.decontam), , drop = FALSE]
metadata_phyloseq <- metadata_phyloseq[match(rownames(seqtab.decontam), metadata_phyloseq$`16S_ID`), , drop = FALSE]
rownames(metadata_phyloseq) <- metadata_phyloseq$`16S_ID`

bacteria_phyloseq <- phyloseq(
  otu_table(as.matrix(seqtab.decontam), taxa_are_rows = FALSE),
  tax_table(as.matrix(taxa_decontam)),
  sample_data(metadata_phyloseq)
)

saveRDS(bacteria_phyloseq, file.path(output_dir, "bacteria_phyloseq.rds"))

################################################################################
# 6. METADATA PREPARATION AND SAMPLE GROUPING
################################################################################

# Standardize treatment names
metadata <- metadata %>%
  mutate(Treatment = str_trim(Treatment)) %>%
  mutate(Treatment = case_when(
    Treatment == "Polyp_Aposymbiotic" ~ "Aposymbiotic",
    Treatment == "Polyp_Smic"         ~ "Polyp_Native",
    Treatment == "Polyp_SSB01"        ~ "Polyp_Control",
    Treatment == "Polyp_antibiotic"   ~ "Polyp_Antibiotic",
    Treatment == "Polyp_Mutant"       ~ "Polyp_Mutant",
    Treatment == "Water"              ~ "Water",
    Treatment == "Algae_KB8"          ~ "Algae_Native",
    Treatment == "Algae_SSB01"        ~ "Algae_Control",
    Treatment == "Algae_antibiotic"   ~ "Algae_Antibiotic",
    Treatment == "Algae_mutant"       ~ "Algae_Mutant",
    TRUE ~ Treatment
  )) %>%
  filter(!grepl("Blank", Treatment, ignore.case = TRUE))

# Define group orders for plotting
algae_order <- c("Algae_Native", "Algae_Control", "Algae_Antibiotic", "Algae_Mutant")
polyp_order <- c("Aposymbiotic", "Polyp_Native", "Polyp_Control",
                 "Polyp_Antibiotic", "Polyp_Mutant")
all_groups_order <- c("Water", algae_order, polyp_order)

# Define color schemes
color_list <- c(
  "Water"             = "#0077BB",
  "Algae_Native"      = "#229944",
  "Algae_Control"     = "#44BB99",
  "Algae_Antibiotic"  = "#66CC88",
  "Algae_Mutant"      = "#55AA55",
  "Aposymbiotic"      = "#850000",
  "Polyp_Native"      = "#EE7733",
  "Polyp_Control"     = "#AA4499",
  "Polyp_Antibiotic"  = "#CC3311",
  "Polyp_Mutant"      = "#DDCC77"
)

# Define shapes
shape_list <- c(
  "Water"             = 16,
  "Algae_Native"      = 15,
  "Algae_Control"     = 16,
  "Algae_Antibiotic"  = 17,
  "Algae_Mutant"      = 19,
  "Aposymbiotic"      = 22,
  "Polyp_Native"      = 15,
  "Polyp_Control"     = 16,
  "Polyp_Antibiotic"  = 17,
  "Polyp_Mutant"      = 19
)

################################################################################
# 7. ALPHA DIVERSITY
################################################################################

# Match decontaminated ASV table to cleaned metadata
shared_ids  <- intersect(rownames(seqtab.decontam), metadata$`16S_ID`)
asv_matched <- seqtab.decontam[shared_ids, ]
meta_matched <- metadata[match(shared_ids, metadata$`16S_ID`), ]

# Calculate alpha diversity metrics
alpha_df <- data.frame(
  SampleID   = rownames(asv_matched),
  Group      = meta_matched$Treatment,
  Shannon    = vegan::diversity(asv_matched, index = "shannon"),
  InvSimpson = vegan::diversity(asv_matched, index = "invsimpson"),
  Observed   = vegan::specnumber(asv_matched),
  Fisher     = vegan::fisher.alpha(asv_matched)
)
alpha_df <- alpha_df[complete.cases(alpha_df), ]

# Convert to long format for plotting
alpha_long <- pivot_longer(alpha_df,
                           cols      = c(Shannon, InvSimpson, Observed, Fisher),
                           names_to  = "Metric",
                           values_to = "Value")

# Function to plot alpha diversity with Tukey post-hoc letters
plot_alpha <- function(df, group_order, colors, title, filename) {
  df <- df %>%
    filter(Group %in% group_order) %>%
    mutate(Group = factor(Group, levels = group_order))

  # Compute Tukey post-hoc letters per metric
  letters_df <- df %>%
    group_by(Metric) %>%
    do({
      aov_model  <- aov(Value ~ Group, data = .)
      tukey      <- multcomp::glht(aov_model, linfct = mcp(Group = "Tukey"))
      cld_result <- multcomp::cld(tukey, level = 0.05)
      data.frame(
        Group   = names(cld_result$mcletters$Letters),
        Letters = cld_result$mcletters$Letters,
        Metric  = unique(.$Metric)
      )
    }) %>%
    ungroup()

  df <- df %>% left_join(letters_df, by = c("Group", "Metric"))

  label_pos <- df %>%
    group_by(Group, Metric) %>%
    summarise(
      y       = max(Value, na.rm = TRUE) * 1.05,
      Letters = unique(Letters),
      .groups = "drop"
    )

  p <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    geom_text(data = label_pos,
              aes(x = Group, y = y, label = Letters),
              inherit.aes = FALSE, size = 4) +
    facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = colors[group_order], limits = group_order) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
      strip.text   = element_text(size = 13),
      legend.position = "none"
    ) +
    labs(title = title, x = "Treatment", y = "Diversity metric")

  ggsave(file.path(output_dir, filename),
         plot = p, width = 8, height = 14, dpi = 300, bg = "white")
  return(p)
}

# Generate alpha diversity plots
plot_alpha(alpha_long, algae_order, color_list,
           "Alpha Diversity — Algae treatments", "alpha_algae.png")
plot_alpha(alpha_long, polyp_order, color_list,
           "Alpha Diversity — Polyp treatments", "alpha_polyps.png")
plot_alpha(alpha_long, all_groups_order, color_list,
           "Alpha Diversity — All treatments", "alpha_all.png")

# Statistical tests
cat("\n=== Alpha Diversity Statistics ===\n")
for (metric in c("Shannon", "InvSimpson", "Observed", "Fisher")) {
  cat("\n---", metric, "---\n")
  sub_df <- alpha_df[alpha_df$Group %in% c(algae_order, polyp_order), ]
  sub_df$Group <- factor(sub_df$Group)

  # Shapiro-Wilk normality test per group
  sw <- tapply(sub_df[[metric]], sub_df$Group, function(x) shapiro.test(x)$p.value)
  cat("Shapiro-Wilk p-values:\n"); print(round(sw, 3))

  # ANOVA or Kruskal-Wallis depending on normality
  if (all(sw > 0.05)) {
    cat("Using ANOVA\n")
    print(summary(aov(as.formula(paste(metric, "~ Group")), data = sub_df)))
  } else {
    cat("Using Kruskal-Wallis\n")
    print(kruskal.test(as.formula(paste(metric, "~ Group")), data = sub_df))
  }
}

################################################################################
# 8. BETA DIVERSITY — BRAY-CURTIS
################################################################################

# Helper function to prepare data subset
prepare_subset <- function(asv_table, metadata, groups) {
  meta_sub <- metadata %>%
    filter(Treatment %in% groups) %>%
    mutate(Treatment = factor(Treatment, levels = groups))
  asv_sub <- asv_table[rownames(asv_table) %in% meta_sub$`16S_ID`, ]
  asv_sub <- asv_sub[, colSums(asv_sub) > 0]
  meta_sub <- meta_sub[match(rownames(asv_sub), meta_sub$`16S_ID`), ]
  dist_bc  <- vegdist(asv_sub, method = "bray")
  list(dist = dist_bc, meta = meta_sub, asv = asv_sub)
}

# PCoA plotting function
plot_pcoa <- function(dist_matrix, metadata, colors, shapes, title, filename) {
  pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  var_exp  <- round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)
  pcoa_df  <- data.frame(
    PC1   = pcoa_res$points[, 1],
    PC2   = pcoa_res$points[, 2],
    Group = metadata$Treatment
  )

  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
    geom_point(size = 3) +
    stat_ellipse(linetype = "dotted", linewidth = 0.5) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "right") +
    labs(
      title = title,
      x     = paste0("PC1 (", var_exp[1], "%)"),
      y     = paste0("PC2 (", var_exp[2], "%)")
    )

  ggsave(file.path(output_dir, filename),
         plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  return(p)
}

# Prepare subsets and run PERMANOVA
datasets_bc <- list(
  algae = prepare_subset(seqtab.decontam, metadata, algae_order),
  polyp = prepare_subset(seqtab.decontam, metadata, polyp_order),
  all   = prepare_subset(seqtab.decontam, metadata, all_groups_order)
)

cat("\n=== PERMANOVA — Bray-Curtis ===\n")
for (nm in names(datasets_bc)) {
  cat("\n---", nm, "---\n")
  print(adonis2(datasets_bc[[nm]]$dist ~ Treatment,
                data = datasets_bc[[nm]]$meta, permutations = 999))
  # Homogeneity of dispersion
  bd <- betadisper(datasets_bc[[nm]]$dist, datasets_bc[[nm]]$meta$Treatment)
  cat("Betadisper permutation test:\n")
  print(permutest(bd))
}

# Generate PCoA plots
plot_pcoa(datasets_bc$algae$dist, datasets_bc$algae$meta,
          color_list[algae_order], shape_list[algae_order],
          "Bray-Curtis PCoA — Algae", "pcoa_bc_algae.png")
plot_pcoa(datasets_bc$polyp$dist, datasets_bc$polyp$meta,
          color_list[polyp_order], shape_list[polyp_order],
          "Bray-Curtis PCoA — Polyps", "pcoa_bc_polyps.png")
plot_pcoa(datasets_bc$all$dist, datasets_bc$all$meta,
          color_list[all_groups_order], shape_list[all_groups_order],
          "Bray-Curtis PCoA — All treatments", "pcoa_bc_all.png")

################################################################################
# 9. BETA DIVERSITY — WEIGHTED UNIFRAC
################################################################################

# Build phylogenetic tree from ASV sequences
build_tree <- function(taxa_table, output_dir) {
  tree_path <- file.path(output_dir, "phylo_tree.rds")
  if (file.exists(tree_path)) {
    cat("Loading existing phylogenetic tree...\n")
    return(readRDS(tree_path))
  }
  cat("Building phylogenetic tree (this may take a while)...\n")
  seqs      <- DNAStringSet(rownames(taxa_table))
  names(seqs) <- rownames(taxa_table)
  alignment <- DECIPHER::AlignSeqs(seqs, anchor = NA)
  phang_align <- phangorn::phyDat(as.matrix(alignment), type = "DNA")
  dm      <- phangorn::dist.ml(phang_align)
  treeNJ  <- phangorn::NJ(dm)
  fit     <- phangorn::pml(treeNJ, data = phang_align)
  fitGTR  <- phangorn::optim.pml(fit, model = "GTR", optInv = TRUE,
                                 optGamma = TRUE, rearrangement = "stochastic",
                                 control = phangorn::pml.control(trace = 0))
  saveRDS(fitGTR$tree, tree_path)
  return(fitGTR$tree)
}

tree <- build_tree(taxa_decontam, output_dir)

# Prepare UniFrac distance for a subset
prepare_unifrac <- function(asv_table, metadata, groups, tree) {
  meta_sub <- metadata %>%
    filter(Treatment %in% groups) %>%
    mutate(Treatment = factor(Treatment, levels = groups))
  asv_sub <- asv_table[rownames(asv_table) %in% meta_sub$`16S_ID`, ]
  asv_sub <- asv_sub[, colSums(asv_sub) > 0]
  common  <- intersect(colnames(asv_sub), tree$tip.label)
  tree_sub <- ape::keep.tip(tree, common)
  asv_sub  <- asv_sub[, common]
  meta_sub <- meta_sub[match(rownames(asv_sub), meta_sub$`16S_ID`), ]
  asv_mat  <- as.matrix(asv_sub)
  storage.mode(asv_mat) <- "numeric"
  physeq <- phyloseq(
    otu_table(asv_mat, taxa_are_rows = FALSE),
    phy_tree(tree_sub)
  )
  dist_uf <- UniFrac(physeq, weighted = TRUE, normalized = TRUE)
  list(dist = dist_uf, meta = meta_sub, asv = asv_sub)
}

datasets_uf <- list(
  algae = prepare_unifrac(seqtab.decontam, metadata, algae_order, tree),
  polyp = prepare_unifrac(seqtab.decontam, metadata, polyp_order, tree),
  all   = prepare_unifrac(seqtab.decontam, metadata, all_groups_order, tree)
)

cat("\n=== PERMANOVA — Weighted UniFrac ===\n")
for (nm in names(datasets_uf)) {
  cat("\n---", nm, "---\n")
  print(adonis2(datasets_uf[[nm]]$dist ~ Treatment,
                data = datasets_uf[[nm]]$meta, permutations = 999))
  bd <- betadisper(datasets_uf[[nm]]$dist, datasets_uf[[nm]]$meta$Treatment)
  cat("Betadisper permutation test:\n")
  print(permutest(bd))
}

# Generate UniFrac PCoA plots
plot_pcoa(datasets_uf$algae$dist, datasets_uf$algae$meta,
          color_list[algae_order], shape_list[algae_order],
          "Weighted UniFrac PCoA — Algae", "pcoa_uf_algae.png")
plot_pcoa(datasets_uf$polyp$dist, datasets_uf$polyp$meta,
          color_list[polyp_order], shape_list[polyp_order],
          "Weighted UniFrac PCoA — Polyps", "pcoa_uf_polyps.png")
plot_pcoa(datasets_uf$all$dist, datasets_uf$all$meta,
          color_list[all_groups_order], shape_list[all_groups_order],
          "Weighted UniFrac PCoA — All treatments", "pcoa_uf_all.png")

################################################################################
# 10. TAXONOMIC COMPOSITION — STACKED BAR PLOTS
################################################################################

plot_taxonomy <- function(asv_table, metadata, taxa, groups,
                          tax_level, colors, filename,
                          min_abundance = 0.01) {
  meta_sub <- metadata %>%
    filter(Treatment %in% groups) %>%
    mutate(Treatment = factor(Treatment, levels = groups))
  asv_sub <- asv_table[rownames(asv_table) %in% meta_sub$`16S_ID`, ]
  meta_sub <- meta_sub[match(rownames(asv_sub), meta_sub$`16S_ID`), ]

  rel_abund <- sweep(asv_sub, 1, rowSums(asv_sub), "/")
  tax_abund <- data.frame(Group = meta_sub$Treatment, rel_abund, check.names = FALSE)
  melted    <- melt(tax_abund, id.vars = "Group")
  melted$Taxonomy <- taxa[melted$variable, tax_level]
  melted$Taxonomy[is.na(melted$Taxonomy)] <- "Unclassified"

  tax_summary <- melted %>%
    group_by(Group, Taxonomy) %>%
    summarise(RelAbund = mean(value), .groups = "drop") %>%
    group_by(Taxonomy) %>%
    filter(max(RelAbund) > min_abundance) %>%
    ungroup()

  n_taxa <- length(unique(tax_summary$Taxonomy))
  pal    <- if (n_taxa <= 12) brewer.pal(max(3, n_taxa), "Paired") else
    colorRampPalette(brewer.pal(12, "Paired"))(n_taxa)

  p <- ggplot(tax_summary, aes(x = Group, y = RelAbund, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
      legend.text  = element_text(size = 9),
      legend.title = element_text(size = 10)
    ) +
    labs(x = "Treatment", y = "Relative abundance",
         title = paste("Taxonomic composition —", tax_level, "level"))

  ggsave(file.path(output_dir, filename),
         plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  return(p)
}

# Family-level composition
plot_taxonomy(seqtab.decontam, metadata, taxa_decontam, algae_order,
              "Family", color_list, "taxonomy_family_algae.png")
plot_taxonomy(seqtab.decontam, metadata, taxa_decontam, polyp_order,
              "Family", color_list, "taxonomy_family_polyps.png")

# Genus-level composition
plot_taxonomy(seqtab.decontam, metadata, taxa_decontam, algae_order,
              "Genus", color_list, "taxonomy_genus_algae.png")
plot_taxonomy(seqtab.decontam, metadata, taxa_decontam, polyp_order,
              "Genus", color_list, "taxonomy_genus_polyps.png")

################################################################################
# 11. DESEQ2 — DIFFERENTIAL ABUNDANCE ANALYSIS
################################################################################

# DESeq2 is run on unrarefied counts using its internal size-factor normalization.
# Pairwise contrasts are extracted between all treatment groups.
# ASVs are considered differentially abundant if:
#   - |log2 fold change| >= 1
#   - Benjamini-Hochberg adjusted p-value < 0.05

run_deseq2 <- function(asv_table, metadata, groups, contrast,
                        taxa_table, label, lfc_threshold = 1) {
  # Subset to groups of interest
  meta_sub <- metadata %>%
    filter(Treatment %in% groups) %>%
    mutate(Treatment = factor(Treatment, levels = groups))
  asv_sub <- asv_table[rownames(asv_table) %in% meta_sub$`16S_ID`, ]
  asv_sub <- asv_sub[, colSums(asv_sub) > 0]
  meta_sub <- meta_sub[match(rownames(asv_sub), meta_sub$`16S_ID`), ]

  # Build DESeq2 object — counts must be integers (round to be safe)
  dds <- DESeqDataSetFromMatrix(
    countData = t(round(asv_sub)),
    colData   = meta_sub,
    design    = ~ Treatment
  )

  # Remove ASVs with very low counts across all samples
  dds <- dds[rowSums(counts(dds)) > 10, ]

  # Run DESeq2
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")

  # Extract results for the specified contrast
  res <- results(dds,
                 contrast  = c("Treatment", contrast[1], contrast[2]),
                 alpha     = 0.05,
                 pAdjustMethod = "BH")
  res_df <- as.data.frame(res)
  res_df$ASV <- rownames(res_df)

  # Add taxonomy
  res_df <- res_df %>%
    left_join(
      data.frame(ASV    = rownames(taxa_table),
                 Family = taxa_table[, "Family"],
                 Genus  = taxa_table[, "Genus"],
                 stringsAsFactors = FALSE),
      by = "ASV"
    )

  # Filter by significance and effect size
  res_sig <- res_df %>%
    filter(!is.na(padj),
           padj < 0.05,
           abs(log2FoldChange) >= lfc_threshold) %>%
    arrange(log2FoldChange)

  cat("\nContrast:", contrast[1], "vs", contrast[2],
      "— significant ASVs:", nrow(res_sig), "\n")

  # Save full results
  write.csv(res_df,
            file.path(output_dir, paste0("deseq2_", label, "_full.csv")),
            row.names = FALSE)
  write.csv(res_sig,
            file.path(output_dir, paste0("deseq2_", label, "_significant.csv")),
            row.names = FALSE)

  return(list(full = res_df, significant = res_sig, dds = dds))
}

# --- 11.1 Pairwise contrasts among algal culture treatments ---
deseq_native_vs_control <- run_deseq2(
  seqtab.decontam, metadata,
  groups   = c("Algae_Native", "Algae_Control"),
  contrast = c("Algae_Native", "Algae_Control"),
  taxa_table = taxa_decontam,
  label    = "algae_native_vs_control"
)

deseq_control_vs_antibiotic <- run_deseq2(
  seqtab.decontam, metadata,
  groups   = c("Algae_Control", "Algae_Antibiotic"),
  contrast = c("Algae_Control", "Algae_Antibiotic"),
  taxa_table = taxa_decontam,
  label    = "algae_control_vs_antibiotic"
)

deseq_control_vs_mutant <- run_deseq2(
  seqtab.decontam, metadata,
  groups   = c("Algae_Control", "Algae_Mutant"),
  contrast = c("Algae_Control", "Algae_Mutant"),
  taxa_table = taxa_decontam,
  label    = "algae_control_vs_mutant"
)

# --- 11.2 Pairwise contrasts among polyp treatments by developmental outcome ---
# Group polyps by outcome: Strobilation (Native, Control, Antibiotic),
# Aposymbiotic, Mutant

# Add outcome grouping to metadata
metadata <- metadata %>%
  mutate(Outcome = case_when(
    Treatment %in% c("Polyp_Native", "Polyp_Control", "Polyp_Antibiotic") ~ "Strobilation",
    Treatment == "Aposymbiotic" ~ "Aposymbiotic",
    Treatment == "Polyp_Mutant" ~ "Mutant",
    TRUE ~ NA_character_
  ))

deseq_strob_vs_apo <- run_deseq2(
  seqtab.decontam,
  metadata %>% mutate(Treatment = Outcome),
  groups   = c("Strobilation", "Aposymbiotic"),
  contrast = c("Strobilation", "Aposymbiotic"),
  taxa_table = taxa_decontam,
  label    = "outcome_strobilation_vs_aposymbiotic"
)

deseq_strob_vs_mutant <- run_deseq2(
  seqtab.decontam,
  metadata %>% mutate(Treatment = Outcome),
  groups   = c("Strobilation", "Mutant"),
  contrast = c("Strobilation", "Mutant"),
  taxa_table = taxa_decontam,
  label    = "outcome_strobilation_vs_mutant"
)

# --- 11.3 Heatmap of differentially abundant genera by outcome ---
# Combine significant ASVs from both outcome contrasts
sig_combined <- bind_rows(
  deseq_strob_vs_apo$significant   %>% mutate(Contrast = "Strob vs Apo"),
  deseq_strob_vs_mutant$significant %>% mutate(Contrast = "Strob vs Mutant")
) %>%
  filter(!is.na(Genus)) %>%
  distinct(ASV, .keep_all = TRUE)

if (nrow(sig_combined) > 0) {
  # Calculate mean relative abundance per genus per outcome group
  outcome_groups <- c("Strobilation", "Aposymbiotic", "Mutant")
  meta_outcome <- metadata %>%
    filter(!is.na(Outcome)) %>%
    mutate(Treatment = Outcome)
  asv_outcome <- seqtab.decontam[rownames(seqtab.decontam) %in% meta_outcome$`16S_ID`, ]
  meta_outcome <- meta_outcome[match(rownames(asv_outcome), meta_outcome$`16S_ID`), ]

  rel_abund <- sweep(asv_outcome, 1, rowSums(asv_outcome), "/") * 100
  sig_asvs  <- intersect(sig_combined$ASV, colnames(rel_abund))
  rel_sig   <- rel_abund[, sig_asvs, drop = FALSE]

  heatmap_df <- data.frame(
    Outcome = meta_outcome$Outcome,
    rel_sig,
    check.names = FALSE
  ) %>%
    group_by(Outcome) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    pivot_longer(-Outcome, names_to = "ASV", values_to = "RelAbund") %>%
    left_join(sig_combined %>% select(ASV, Genus, log2FoldChange), by = "ASV") %>%
    filter(!is.na(Genus)) %>%
    mutate(Outcome = factor(Outcome, levels = outcome_groups))

  p_heatmap <- ggplot(heatmap_df,
                      aes(x = Outcome, y = reorder(Genus, log2FoldChange),
                          fill = RelAbund)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low  = "#2166AC",
                         mid  = "white",
                         high = "#D6604D",
                         midpoint = 0,
                         name = "Mean relative\nabundance (%)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text.y  = element_text(face = "italic"),
          panel.grid   = element_blank()) +
    labs(title = "Differentially abundant genera by developmental outcome",
         x = NULL, y = NULL)

  ggsave(file.path(output_dir, "deseq2_outcome_heatmap.png"),
         plot = p_heatmap, width = 8, height = 10, dpi = 300, bg = "white")
}

################################################################################
# 12. UPSET PLOTS — ASV SHARING ACROSS TREATMENTS
################################################################################


make_upset <- function(asv_table, metadata, groups, label,
                       min_prev = 0.3) {
  # min_prev: minimum proportion of samples within a group in which
  # an ASV must be present to be counted as present in that group

  meta_sub <- metadata %>%
    filter(Treatment %in% groups) %>%
    mutate(Treatment = factor(Treatment, levels = groups))
  asv_sub <- asv_table[rownames(asv_table) %in% meta_sub$`16S_ID`, ]
  asv_sub <- asv_sub[, colSums(asv_sub) > 0]
  meta_sub <- meta_sub[match(rownames(asv_sub), meta_sub$`16S_ID`), ]

  # Presence/absence matrix
  pa <- (asv_sub > 0) * 1

  # Per-group prevalence: is ASV present in >= min_prev of samples in that group?
  group_presence <- sapply(groups, function(g) {
    idx <- meta_sub$Treatment == g
    if (sum(idx) == 0) return(rep(FALSE, ncol(pa)))
    colMeans(pa[idx, , drop = FALSE]) >= min_prev
  })
  group_presence <- as.data.frame(group_presence)
  rownames(group_presence) <- colnames(pa)

  # Keep only ASVs present in at least one group
  group_presence <- group_presence[rowSums(group_presence) > 0, ]

  # Convert to logical for ComplexUpset
  group_presence[] <- lapply(group_presence, as.logical)

  # Print intersection summary
  cat("\nUpSet summary for", label, ":\n")
  cat("Total ASVs included:", nrow(group_presence), "\n")
  cat("ASVs per group:\n")
  print(colSums(group_presence))

  # Save ASV membership table
  write.csv(group_presence,
            file.path(output_dir, paste0("upset_membership_", label, ".csv")))

  # Plot
  p <- upset(
    group_presence,
    intersect   = groups,
    base_annotations = list(
      "Intersection size" = intersection_size(
        counts = TRUE,
        text   = list(size = 3)
      )
    ),
    set_sizes = upset_set_size(geom = geom_bar(fill = "#4477AA")),
    width_ratio = 0.2,
    name = label
  )

  ggsave(file.path(output_dir, paste0("upset_", label, ".png")),
         plot = p, width = 12, height = 6, dpi = 300, bg = "white")

  return(group_presence)
}

# UpSet for algal culture treatments
upset_algae <- make_upset(
  seqtab.decontam, metadata,
  groups = algae_order,
  label  = "algae_cultures"
)

# UpSet for polyp treatments
upset_polyp <- make_upset(
  seqtab.decontam, metadata,
  groups = polyp_order,
  label  = "polyp_treatments"
)

# UpSet for developmental outcome groups
metadata_outcome <- metadata %>%
  filter(!is.na(Outcome)) %>%
  mutate(Treatment = Outcome)

upset_outcome <- make_upset(
  seqtab.decontam, metadata_outcome,
  groups = c("Strobilation", "Aposymbiotic", "Mutant"),
  label  = "developmental_outcome"
)




