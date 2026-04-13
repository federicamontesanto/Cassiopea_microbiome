################################################################################
# Cassiopea microbiome project
#   DESeq2-normalized microbiome heatmaps across Cassiopea treatment subsets.
#
################################################################################
 

required_packages <- c(
  "phyloseq",
  "DESeq2",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "grid",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer",
  "viridis"
)

# ---------------------------------------------------------------------------
# Path configuration
# ---------------------------------------------------------------------------
project_root <- "."
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "05_differential_abundance/publication_heatmaps_deseq2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
phyloseq_bact_file <- file.path(data_dir, "bacteria_phyloseq.rds")


# ----------------------------------------------------------
# Treatment name mapping for cleaner display - CORRECTED
# ----------------------------------------------------------
treatment_name_mapping <- c(
  "Water"              = "Water",
  # ALGAE 
  "Algae_Smic_KB8"     = "Algae_Native",
  "Algae_Bmin_SSB01"   = "Algae_Control",      
  "Algae_antibiotic"   = "Algae_Antibiotic",
  "Algae_mutant"       = "Algae_Muantt",
  # POLYPS
  "Polyp_Aposymbiotic" = "Aposymbiotic",
  "Polyp_Smic_KB8"     = "Polyp_Native",
  "Polyp_Bmin_SSB01"   = "Polyp_Control",
  "Polyp_antibiotic"   = "Polyp_Antibiotic",
  "Polyp_Mutant"       = "Polyp_Mutant",
)

# Extract treatment groups from mapping keys
algae_treatments <- names(treatment_name_mapping)[grepl("^Algae_", names(treatment_name_mapping))]
polyp_treatments <- names(treatment_name_mapping)[grepl("^Polyp_", names(treatment_name_mapping))]
blank_treatments <- names(treatment_name_mapping)[grepl("^Blank_", names(treatment_name_mapping))]



# Function to apply name mapping
apply_treatment_mapping <- function(names) {
  sapply(names, function(name) {
    ifelse(name %in% names(treatment_name_mapping), 
           treatment_name_mapping[name], 
           name)
  }, USE.NAMES = FALSE)
}

# ----------------------------------------------------------
# Publication-Quality Settings
# ----------------------------------------------------------
heatmap_colors <- colorRampPalette(c("#313695","#4575B4","#74ADD1",
                                     "#ABD9E9","#E0F3F8","#FFFFBF",
                                     "#FEE090","#FDAE61","#F46D43",
                                     "#D73027","#A50026"))(100)
heatmap_colors_viridis <- viridis(100, option="B", direction=-1)

# Color palette
group_colors_annotation <- c(
  "Water" = "#4292c6",      # Blue for Water
  "Polyp" = "#fd8d3c",      # Orange for Polyp  
  "Algae" = "#41ab5d",      # Green for Algae
)

# Treatment order
treatment_order <- c(
  "Water",
  algae_treatments,
  polyp_treatments
)

# ----------------------------------------------------------
# Core Heatmap Function with Group Filtering 
# ----------------------------------------------------------
create_publication_heatmap <- function(
    normalized_matrix,
    taxonomy_df,
    metadata_df,
    taxonomic_level     = "Family",
    top_n               = 25,
    output_prefix       = "Microbiome",
    exclude_water       = FALSE,
    group_filter        = NULL,  
    min_prevalence      = 0.1,
    scale_method        = "row",
    clustering_distance = "euclidean",
    clustering_method   = "complete",
    show_names          = TRUE,
    fig_width           = 12,
    fig_height          = 10,
    color_scheme        = "RdBu"
) {
  
  # Determine suffix based on filtering
  if (!is.null(group_filter)) {
    filter_suffix <- paste0(toupper(substring(group_filter, 1, 1)), 
                            substring(group_filter, 2), "Only")
  } else if (exclude_water) {
    filter_suffix <- "NoWater"
  } else {
    filter_suffix <- "AllGroups"
  }
  
  cat(sprintf("\n--- %s: %s-level heatmap (top %d) [%s] ---\n", 
              output_prefix, taxonomic_level, top_n, filter_suffix))
  
  # 1) Data preparation
  norm_long <- as.data.frame(normalized_matrix) %>%
    rownames_to_column("ASV") %>%
    pivot_longer(-ASV, names_to="SampleID", values_to="Abundance") %>%
    left_join(
      if("SampleID" %in% colnames(metadata_df)) metadata_df
      else metadata_df %>% rownames_to_column("SampleID"),
      by="SampleID"
    )
  
  # Standardize SampleIDs - phyloseq converts hyphens to periods
  norm_long <- norm_long %>%
    mutate(SampleID = gsub("-", ".", SampleID))
  
  # Filter by group if specified
  if (!is.null(group_filter)) {
    if (tolower(group_filter) == "algae") {
      norm_long <- norm_long %>% filter(Treatment %in% algae_treatments)
      cat("Filtered to ALGAE treatments only\n")
    } else if (tolower(group_filter) == "polyp") {
      norm_long <- norm_long %>% filter(Treatment %in% polyp_treatments)
      cat("Filtered to POLYP treatments only\n")
    }
  }
  
  # Taxonomy joining
  tax_for_join <- taxonomy_df %>%
    dplyr::select(ASV_Seq, !!taxonomic_level) %>%
    dplyr::rename(ASV = ASV_Seq, TaxLevel = !!taxonomic_level)
  
  norm_long <- norm_long %>%
    left_join(tax_for_join, by = "ASV") %>%
    filter(!is.na(TaxLevel),
           !TaxLevel %in% c("Unclassified","Unknown",""))
  
  # 2) Prevalence filter
  taxon_prev <- norm_long %>%
    group_by(TaxLevel) %>%
    summarise(prev = sum(Abundance>0)/n_distinct(SampleID), .groups="drop") %>%
    filter(prev >= min_prevalence)
  cat(sprintf("Taxa ≥ %.0f%% prevalence: %d\n", min_prevalence*100, nrow(taxon_prev)))
  norm_long <- norm_long %>%
    filter(TaxLevel %in% taxon_prev$TaxLevel)
  
  # 3) Mean abundance per treatment
  abund_sum <- norm_long %>%
    group_by(Treatment, TaxLevel) %>%
    summarise(MeanAbundance = mean(Abundance, na.rm=TRUE), .groups="drop")
  mat <- abund_sum %>%
    pivot_wider(names_from=Treatment, values_from=MeanAbundance, values_fill=0) %>%
    column_to_rownames("TaxLevel")
  cat(sprintf("Matrix dims: %dx%d\n", nrow(mat), ncol(mat)))
  
  # 4) Select top variable taxa
  if(ncol(mat) > 1) {
    stats <- data.frame(
      taxon = rownames(mat),
      mean  = rowMeans(mat),
      sd    = apply(mat,1,sd),
      var   = apply(mat,1,var)
    )
    stats$cv <- with(stats, ifelse(mean>0, sd/mean, 0))
    stats <- subset(stats, var>1e-10)
    stats <- stats[order(stats$cv, decreasing=TRUE), ]
    sel <- head(stats$taxon, max(2, min(top_n, nrow(stats))))
    mat <- mat[sel, , drop=FALSE]
    cat(sprintf("Selected %d taxa\n", nrow(mat)))
  }
  if(nrow(mat) < 2) {
    warning("Need ≥2 taxa after filtering. Skipping this heatmap.")
    return(invisible(NULL))
  }
  
  # 5) Order & exclude water
  cols <- intersect(treatment_order, colnames(mat))
  mat <- mat[, cols, drop=FALSE]
  if(exclude_water && "Water" %in% colnames(mat)) {
    mat <- mat[, setdiff(colnames(mat),"Water"), drop=FALSE]
    cat("Excluded Water\n")
  }
  
  
  # 6) APPLY TREATMENT NAME MAPPING to column names
  original_colnames <- colnames(mat)
  display_colnames <- apply_treatment_mapping(original_colnames)
  
  cat("Treatment name mapping applied:\n")
  for(i in seq_along(original_colnames)) {
    if(original_colnames[i] != display_colnames[i]) {
      cat(sprintf("  %s -> %s\n", original_colnames[i], display_colnames[i]))
    }
  }
  
  # Rename the matrix columns for display
  colnames(mat) <- display_colnames
  
  # 7) Annotation setup 
  ann <- data.frame(
    Group = factor(
      case_when(
        display_colnames=="Water"                  ~ "Water",
        grepl("^Algae", display_colnames)          ~ "Algae",
        grepl("^Polyp", display_colnames)          ~ "Polyp",
        display_colnames=="Apo"                    ~ "Polyp",
        TRUE                                       ~ NA_character_
      ), levels = if (!is.null(group_filter)) {
        if (tolower(group_filter) == "algae") "Algae" else "Polyp"
      } else if (exclude_water) {
        c("Polyp","Algae")
      } else {
        c("Water","Polyp","Algae")
      }
    ), row.names = display_colnames
  )
  
  # Remove columns without annotation
  valid_cols <- !is.na(ann$Group)
  mat <- mat[, valid_cols, drop = FALSE]
  ann <- ann[valid_cols, , drop = FALSE]
  
  cat("Group order applied:\n")
  cat("  Levels:", paste(levels(ann$Group), collapse = " → "), "\n")
  cat("  Groups in data:", paste(unique(ann$Group), collapse = ", "), "\n")
  
  # 8) Color palette selection
  cols_pal <- switch(color_scheme,
                     RdBu    = heatmap_colors,
                     viridis = heatmap_colors_viridis,
                     colorRampPalette(brewer.pal(9, color_scheme))(100)
  )
  
  # 9) Scaling
  plot_mat <- if(scale_method=="row") t(scale(t(mat))) else mat
  
  # 10) Handle NA/Inf
  plot_mat[!is.finite(plot_mat)] <- NA
  plot_mat <- plot_mat[complete.cases(plot_mat), , drop=FALSE]
  if(nrow(plot_mat) < 2) {
    warning("Insufficient finite taxa after NA removal. Skipping this heatmap.")
    return(invisible(NULL))
  }
  
  # 11) Legend breaks & mapping
  if(scale_method=="row") {
    brk <- c(-2,-1,0,1,2); lbl <- as.character(brk)
  } else {
    brk <- pretty(range(plot_mat), 5); lbl <- formatC(brk, format="f", digits=1)
  }
  cols_at <- cols_pal[round(seq(1, length(cols_pal), length.out=length(brk)))]
  col_fun <- colorRamp2(brk, cols_at)
  
  # 12) Build heatmap with renamed column labels
  ht <- Heatmap(
    matrix         = plot_mat,
    name           = if(scale_method=="row") "Z-score" else "Abundance",
    col            = col_fun,
    cluster_rows   = TRUE,
    clustering_distance_rows = clustering_distance,
    clustering_method_rows   = clustering_method,
    cluster_columns=FALSE,
    show_row_names=show_names,
    show_column_names=show_names,
    row_names_gp   = gpar(fontsize=14, fontface = if(taxonomic_level=="Genus") "italic" else "plain"),
    column_names_gp = gpar(fontsize=14),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize=12, fontface="bold"),
    width          = unit(ncol(plot_mat)*8, "mm"),
    height         = unit(nrow(plot_mat)*8, "mm"),
    column_split   = if(length(unique(ann$Group)) > 1) ann$Group else NULL,
    top_annotation = HeatmapAnnotation(
      Group = ann$Group,
      col   = list(Group=group_colors_annotation),
      annotation_name_gp = gpar(fontsize=12, fontface="bold")
    ),
    heatmap_legend_param = list(
      at     = brk,
      labels = lbl,
      title  = if(scale_method=="row") "Z-score" else "Abundance",
      legend_height = unit(4, "cm"),
      title_gp      = gpar(fontsize=12, fontface="bold"),
      labels_gp     = gpar(fontsize=10)
    ),
    column_gap = unit(2, "mm")
  )
  
  # 13) Save 
  base_filename <- sprintf("%s_%s_Top%d_%s",
                           output_prefix, taxonomic_level, top_n, filter_suffix)
  
  pdf(paste0(base_filename, ".pdf"),
      width = fig_width, height = fig_height)
  draw(ht, heatmap_legend_side="right")
  dev.off()
  
  png(paste0(base_filename, ".png"),
      width = fig_width, height = fig_height,
      units = "in", res = 300, bg = "white")
  draw(ht, heatmap_legend_side="right")
  dev.off()
  
  cat(sprintf("Saved: %s.pdf/png\n", base_filename))
  
  invisible(plot_mat)
}

# ----------------------------------------------------------
# Generate Heatmaps
# ----------------------------------------------------------
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("LOADING AND PREPARING DATA\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Load bacteria phyloseq
ps_bacteria <- readRDS(phyloseq_bact_file)

# Clean treatment names
sample_data(ps_bacteria)$Treatment <- trimws(as.character(sample_data(ps_bacteria)$Treatment))


cat("\n=== Filtering samples ===\n")
cat("Before filtering:", nsamples(ps_bacteria), "samples\n")


keep_idx <- !(sample_data(ps_bacteria)$Treatment %in% treatments_to_exclude)
ps_bacteria <- prune_samples(keep_idx, ps_bacteria)

cat("After filtering:", nsamples(ps_bacteria), "samples\n")

# Extract components
otu_b <- as(otu_table(ps_bacteria), "matrix")
if(taxa_are_rows(ps_bacteria)) otu_b <- t(otu_b)

tax_b <- as(tax_table(ps_bacteria), "matrix") %>% 
  as.data.frame() %>% 
  rownames_to_column("ASV_Seq")

meta_b <- as(sample_data(ps_bacteria), "data.frame")

# Standardize SampleIDs in metadata
meta_b <- meta_b %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = gsub("-", ".", SampleID)) %>%
  column_to_rownames("SampleID")

# Standardize rownames in otu table
rownames(otu_b) <- gsub("-", ".", rownames(otu_b))

# Remove singleton treatments for DESeq2
sing <- names(table(meta_b$Treatment)[table(meta_b$Treatment)<=1])
if(length(sing)) { 
  keep <- !meta_b$Treatment %in% sing
  otu_b <- otu_b[keep, ]
  meta_b <- meta_b[keep, ]
  cat("Removed singleton treatments:", paste(sing, collapse = ", "), "\n")
}

# DESeq2 normalization
cat("\n=== Running DESeq2 normalization ===\n")
dds_b <- DESeqDataSetFromMatrix(countData=t(otu_b), colData=meta_b, design=~Treatment)
dds_b <- tryCatch(
  estimateSizeFactors(dds_b, type="poscounts"), 
  error=function(e) estimateSizeFactors(dds_b)
)
vsd_b <- tryCatch(
  vst(dds_b, blind=FALSE), 
  error=function(e) tryCatch(
    varianceStabilizingTransformation(dds_b, blind=FALSE), 
    error=function(e2) rlog(dds_b, blind=FALSE)
  )
)
norm_b <- assay(vsd_b)

cat("Normalization complete. Matrix dimensions:", dim(norm_b), "\n")

# ----------------------------------------------------------
# Generate all heatmap combinations
# ----------------------------------------------------------
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("GENERATING HEATMAPS\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# 1. ALGAE-ONLY HEATMAPS 
cat("\n### ALGAE-ONLY HEATMAPS ###\n")

create_publication_heatmap(
  norm_b, tax_b, meta_b, 
  taxonomic_level = "Family", 
  top_n = 20, 
  output_prefix = "Bacteria", 
  group_filter = "algae",
  min_prevalence = 0.1
)

create_publication_heatmap(
  norm_b, tax_b, meta_b, 
  taxonomic_level = "Genus", 
  top_n = 20, 
  output_prefix = "Bacteria", 
  group_filter = "algae",
  min_prevalence = 0.05
)

# 2. POLYP-ONLY HEATMAPS
cat("\n### POLYP-ONLY HEATMAPS ###\n")

create_publication_heatmap(
  norm_b, tax_b, meta_b, 
  taxonomic_level = "Family", 
  top_n = 20, 
  output_prefix = "Bacteria", 
  group_filter = "polyp",
  min_prevalence = 0.1
)

create_publication_heatmap(
  norm_b, tax_b, meta_b, 
  taxonomic_level = "Genus", 
  top_n = 20, 
  output_prefix = "Bacteria", 
  group_filter = "polyp",
  min_prevalence = 0.05
)

# 3. ALL GROUPS (with Water)
cat("\n### ALL GROUPS (WITH WATER) ###\n")

create_publication_heatmap(
  norm_b, tax_b, meta_b, 
  taxonomic_level = "Family", 
  top_n = 20, 
  output_prefix = "Bacteria", 
  min_prevalence = 0.1
)

create_publication_heatmap(
  norm_b, tax_b, meta_b, 
  taxonomic_level = "Genus", 
  top_n = 20, 
  output_prefix = "Bacteria", 
  min_prevalence = 0.05
)

# 4. NO WATER (Algae + Polyp combined)
cat("\n### NO WATER (ALGAE + POLYP COMBINED) ###\n")

create_publication_heatmap(
  norm_b, tax_b, meta_b, 
  taxonomic_level = "Family", 
  top_n = 20, 
  output_prefix = "Bacteria", 
  exclude_water = TRUE, 
  min_prevalence = 0.1
)

create_publication_heatmap(
  norm_b, tax_b, meta_b, 
  taxonomic_level = "Genus", 
  top_n = 20, 
  output_prefix = "Bacteria", 
  exclude_water = TRUE, 
  min_prevalence = 0.05
)

