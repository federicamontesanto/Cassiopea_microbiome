################################################################################
# Cassiopea microbiome project
#   UpSet plot workflow summarizing ASV overlap patterns across bacterial treatment groups.
#
################################################################################

# ===============================================================================
# Bacterial ASV UpSet Plot Analysis 
# ===============================================================================

# ---------------------------------------------------------------------------
# Path configuration
# ---------------------------------------------------------------------------
project_root <- "."
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "05_differential_abundance/upset_asv_overlap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
phyloseq_bact_file <- file.path(data_dir, "bacteria_phyloseq.rds")

grouping_column <- "Treatment"
output_prefix <- "Bacterial_ASV_UpSet"

### Rarefaction params
min_reads_threshold <- 5000
rarefaction_depth   <- 6000
rarefaction_seed    <- 12345

required_packages <- c("phyloseq")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install required packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}
invisible(lapply(required_packages, library, character.only = TRUE))

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  message("Package `VennDiagram` not available; Venn diagram output will be skipped.")
}

# ===============================================================================
# Treatment name mapping and color scheme
# ===============================================================================

treatment_name_mapping <- c(
  "Water"              = "Water",
  # ALGAE
  "Algae_Smic_KB8"     = "Algae_Native",
  "Algae_Bmin_SSB01"   = "Algae_Control",
  "Algae_antibiotic"   = "Algae_Antibiotic",
  "Algae_mutant"       = "Algae_Mutant",
  # POLYPS
  "Polyp_Aposymbiotic" = "Aposymbiotic",
  "Polyp_Smic_KB8"     = "Polyp_Native",
  "Polyp_Bmin_SSB01"   = "Polyp_Control",
  "Polyp_antibiotic"   = "Polyp_Antibiotic",
  "Polyp_Mutant"       = "Polyp_Mutant",
 )

custom_colors <- c(
  "Water"             = "#2166AC",  # Deep blue
  "Algae_Native"      = "#1B7837",  # Forest green
  "Algae_Control"     = "#5AAE61",  # Medium green
  "Algae_Antibiotic"  = "#7FBC41",  # Light green
  "Algae_Mutant"      = "#C2E681",  # Very light green
  "Aposymbiotic"      = "#762A83",  # Deep purple
  "Polyp_Native"      = "#D73027",  # Red
  "Polyp_Control"     = "#F46D43",  # Orange-red
  "Polyp_Antibiotic"  = "#FDAE61",  # Orange
  "Polyp_Mutat"       = "#E6F598"   # Light yellow-green
)

# AUTOMATIC: Extract treatment groups from mapping keys
algae_treatments <- names(treatment_name_mapping)[grepl("^Algae_", names(treatment_name_mapping))]
polyp_treatments <- names(treatment_name_mapping)[grepl("^Polyp_", names(treatment_name_mapping))]
blank_treatments <- names(treatment_name_mapping)[grepl("^Blank_", names(treatment_name_mapping))]

# Define specific exclusions
exclude_treatments <- c("Polyp_antifungal")

# Remove excluded treatments from polyp group
polyp_treatments <- setdiff(polyp_treatments, exclude_treatments)

# Print what was extracted for verification
cat("\n=== Automatically extracted treatment groups ===\n")
cat("Algae treatments:", paste(algae_treatments, collapse = ", "), "\n")
cat("Polyp treatments:", paste(polyp_treatments, collapse = ", "), "\n")
cat("Blank treatments:", paste(blank_treatments, collapse = ", "), "\n")
cat("Excluded treatments:", paste(exclude_treatments, collapse = ", "), "\n\n")

# Define order 
polyp_order <- c("Apo", "Polyp_Native", "Polyp_Control", "Polyp_Antibiotic", "Polyp_Mutant")
algae_order <- c("Algae_Native", "Algae_Control", "Algae_Antibiotic", "Algae_Mutant")

# Function to get colors for display names
get_colors <- function(group_names) {
  colors <- sapply(group_names, function(group) {
    ifelse(group %in% names(custom_colors), custom_colors[group], "#999999")
  })
  return(unname(colors))
}

# Function to convert original names to display names
convert_names <- function(original_names) {
  display_names <- sapply(original_names, function(name) {
    ifelse(name %in% names(treatment_name_mapping), treatment_name_mapping[name], name)
  })
  return(unname(display_names))
}

# Function to order groups according to specified order
order_groups <- function(group_names, polyps_first = TRUE) {
  # Separate algae and polyp groups
  algae_groups <- intersect(algae_order, group_names)
  polyp_groups <- intersect(polyp_order, group_names)
  other_groups <- setdiff(group_names, c(algae_groups, polyp_groups))
  
  if (polyps_first) {
    return(c(polyp_groups, algae_groups, other_groups))
  } else {
    return(c(algae_groups, polyp_groups, other_groups))
  }
}

# ===============================================================================
# Load, filter, rarefy, and process phyloseq data
# ===============================================================================
cat("Loading phyloseq object from:", phyloseq_bact_file, "\n")
ps_bact <- readRDS(phyloseq_bact_file)
cat("Phyloseq object loaded successfully!\n")
cat("Number of samples:", nsamples(ps_bact), "\n")
cat("Number of ASVs:", ntaxa(ps_bact), "\n")

# Clean treatment names
sample_data(ps_bact)$Treatment <- trimws(as.character(sample_data(ps_bact)$Treatment))

cat("\n=== Filtering samples ===\n")
treatment_col <- sample_data(ps_bact)[[grouping_column]]
keep_idx <- !(treatment_col %in% c(blank_treatments, exclude_treatments))
if (any(!keep_idx)) {
  cat("Excluding", sum(!keep_idx), "samples (blanks + Polyp_antifungal)\n")
  cat("Keeping Algae_antifungal samples\n")
}
ps_bact <- prune_samples(keep_idx, ps_bact)

# Pre-filter by read threshold before rarefaction
reads <- sample_sums(ps_bact)
keep_reads <- reads >= min_reads_threshold
if (any(!keep_reads)) {
  cat("Excluding", sum(!keep_reads), "samples with reads <", min_reads_threshold, "\n")
}
ps_bact <- prune_samples(keep_reads, ps_bact)

# Rarefaction (reproducible)
set.seed(rarefaction_seed)
cat("Rarefying to", rarefaction_depth, "reads per sample...\n")
ps_bact <- rarefy_even_depth(
  ps_bact,
  sample.size = rarefaction_depth,
  rngseed = rarefaction_seed,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = TRUE
)
cat("Post-rarefaction samples:", nsamples(ps_bact), " ASVs:", ntaxa(ps_bact), "\n")

# Extract metadata from ps_bact 
meta <- data.frame(
  SampleID = sample_names(ps_bact),
  Treatment = as.character(sample_data(ps_bact)$Treatment),
  stringsAsFactors = FALSE
)

# Clean treatment names
meta$Treatment <- trimws(meta$Treatment)

# Standardize SampleIDs - phyloseq converts hyphens to periods
meta$SampleID <- gsub("-", ".", meta$SampleID)
rownames(meta) <- meta$SampleID

cat("\n=== SampleID Standardization ===\n")
cat("First 5 SampleIDs after standardization:\n")
print(head(meta$SampleID, 5))

# Now extract the sample_df for analysis
sample_df <- meta

cat("Sample data dimensions:", dim(sample_df), "\n")

# Debug: Print all unique groups found
cat("\n=== GROUP DETECTION ===\n")
all_groups_found <- unique(sample_df$Treatment)
cat("All groups found in data:\n")
for (i in seq_along(all_groups_found)) {
  cat(sprintf("  %d. %s\n", i, all_groups_found[i]))
}

# ===============================================================================
# Extract OTU table 
# ===============================================================================

otu_table_obj <- otu_table(ps_bact)
if (taxa_are_rows(ps_bact)) {
  cat("Taxa are rows, transposing for analysis\n")
  otu_df <- as.data.frame(t(otu_table_obj))
} else {
  cat("Samples are rows, using as is\n")
  otu_df <- as.data.frame(otu_table_obj)
}

# Ensure OTU table rownames match standardized SampleIDs
rownames(otu_df) <- gsub("-", ".", rownames(otu_df))

cat("\n=== OTU Table Info ===\n")
cat("OTU table dimensions:", dim(otu_df), "\n")
cat("First 5 OTU table rownames:\n")
print(head(rownames(otu_df), 5))

# Verify matching between sample_df and otu_df
cat("\n=== Verifying SampleID Matching ===\n")
matching_samples <- sum(rownames(sample_df) %in% rownames(otu_df))
cat("Matching samples:", matching_samples, "out of", nrow(sample_df), "\n")

if (matching_samples == 0) {
  cat("ERROR: No matching samples between metadata and OTU table!\n")
  stop("Cannot proceed without matching samples")
}

# Get groups and exclude Water
groups <- unique(sample_df$Treatment)
groups_no_water <- groups[!grepl("^Water$", groups, ignore.case = TRUE)]

cat("\nGroups excluding Water:", paste(groups_no_water, collapse = ", "), "\n")

# AUTOMATIC: Filter to get algae and polyp sets that exist in data
algae_sets_original <- intersect(algae_treatments, groups_no_water)
polyp_sets_original <- intersect(polyp_treatments, groups_no_water)

combined_algae_polyp <- c(polyp_sets_original, algae_sets_original)  


# ===============================================================================
# Create ASV lists for each group 
# ===============================================================================
cat("\n=== CREATING ASV LISTS ===\n")
asv_list <- list()

for (group in groups_no_water) {
  cat("\nProcessing group:", group, "\n")
  
  # Get samples in this group
  samples_in_group <- rownames(sample_df)[sample_df$Treatment == group]
  
  if (length(samples_in_group) > 0) {
    cat("  Samples in group:", length(samples_in_group), "\n")
    
    # Check for matching samples in OTU table
    matching_samples <- intersect(samples_in_group, rownames(otu_df))
    cat("  Matching samples in OTU table:", length(matching_samples), "\n")
    
    if (length(matching_samples) > 0) {
      # Extract OTU data for this group
      group_otus <- otu_df[matching_samples, , drop = FALSE]
      
      # Get ASVs present in this group (sum > 0 across samples)
      present_asvs <- colnames(group_otus)[colSums(group_otus) > 0]
      asv_list[[group]] <- present_asvs
      cat("  ASVs found:", length(present_asvs), "\n")
    } else {
      cat("  ERROR: No matching samples found in OTU table\n")
      asv_list[[group]] <- character(0)
    }
  } else {
    cat("  ERROR: No samples found for this group\n")
    asv_list[[group]] <- character(0)
  }
}

# Verify ASV lists
cat("\n=== ASV LIST VERIFICATION ===\n")
cat("\nALGAE GROUPS:\n")
for (group in algae_sets_original) {
  display_name <- convert_names(group)
  asv_count <- length(asv_list[[group]])
  cat(sprintf("  %-25s -> %-15s: %d ASVs\n", group, display_name, asv_count))
}
cat("\nPOLYP GROUPS:\n")
for (group in polyp_sets_original) {
  display_name <- convert_names(group)
  asv_count <- length(asv_list[[group]])
  cat(sprintf("  %-25s -> %-15s: %d ASVs\n", group, display_name, asv_count))
}

# ===============================================================================
# Custom plots 
# ===============================================================================

create_custom_upset_plot <- function(asv_list, set_names, filename, title_text) {
  cat("\n=== UPSET PLOT CREATION START ===\n")
  cat("Input set_names:", paste(set_names, collapse = ", "), "\n")
  
  # Convert to display names for ordering
  display_names_for_ordering <- convert_names(set_names)
  ordered_display <- order_groups(display_names_for_ordering, polyps_first = FALSE)  # Algae first for algae-only
  
  # Map back to original names
  ordered_set_names <- sapply(ordered_display, function(dn) {
    names(treatment_name_mapping)[treatment_name_mapping == dn][1]
  })
  ordered_set_names <- ordered_set_names[!is.na(ordered_set_names)]
  
  cat("After ordering:", paste(ordered_set_names, collapse = ", "), "\n")
  
  # Remove groups with 0 ASVs
  ordered_set_names <- ordered_set_names[sapply(ordered_set_names, function(x) length(asv_list[[x]]) > 0)]
  
  if (length(ordered_set_names) < 2) {
    cat("ERROR: Insufficient groups with data\n")
    return(NULL)
  }
  
  n_sets <- length(ordered_set_names)
  display_names <- convert_names(ordered_set_names)
  all_asvs <- unique(unlist(asv_list[ordered_set_names]))
  
  cat("Final sets:", n_sets, "| Total ASVs:", length(all_asvs), "\n")
  
  # Create presence/absence matrix
  presence_matrix <- matrix(0, nrow = length(all_asvs), ncol = n_sets)
  colnames(presence_matrix) <- display_names
  
  for (i in 1:n_sets) {
    presence_matrix[all_asvs %in% asv_list[[ordered_set_names[i]]], i] <- 1
  }
  
  # Find unique intersection patterns
  patterns <- unique(presence_matrix)
  pattern_counts <- apply(patterns, 1, function(p) {
    sum(apply(presence_matrix, 1, function(x) all(x == p)))
  })
  
  # Sort and keep top 40
  ord <- order(pattern_counts, decreasing = TRUE)
  patterns <- patterns[ord, , drop = FALSE]
  pattern_counts <- pattern_counts[ord]
  
  n_show <- min(40, length(pattern_counts))
  patterns <- patterns[1:n_show, , drop = FALSE]
  pattern_counts <- pattern_counts[1:n_show]
  
  # Create plot
  png(filename, width = 35, height = 20, units = "in", res = 300)
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.7, 0.3))
  
  # Top panel: Bar chart
  par(mar = c(1, 25, 4, 2), family = "Arial")
  bp <- barplot(pattern_counts, 
                col = "#2C3E50",
                border = NA,
                main = title_text,
                cex.main = 4.0,
                ylim = c(0, max(pattern_counts) * 1.1),
                axes = FALSE,
                axisnames = FALSE)
  
  mtext("Intersection Size", side = 2, line = 4.5, cex = 3.5, family = "Arial")
  
  y_ticks <- pretty(c(0, max(pattern_counts)))
  axis(2, at = y_ticks, labels = y_ticks, 
       cex.axis = 2.5, family = "Arial", las = 1)
  
  text(bp, pattern_counts + max(pattern_counts) * 0.02, 
       labels = pattern_counts, pos = 3, cex = 3.0, family = "Arial")
  
  if (max(pattern_counts) > 0) {
    abline(h = seq(0, max(pattern_counts), length.out = 5), 
           col = "grey90", lty = 2, lwd = 2)
  }
  
  # Bottom panel: Dot matrix
  par(mar = c(4, 25, 0, 2), family = "Arial")
  plot(0, 0, type = "n",
       xlim = c(0.5, n_show + 0.5),
       ylim = c(0.5, n_sets + 0.5),
       xlab = "", ylab = "", axes = FALSE)
  
  for (i in 1:n_show) {
    active_sets <- which(patterns[i,] == 1)
    if (length(active_sets) > 0) {
      points(rep(i, length(active_sets)), active_sets, 
             pch = 19, cex = 5.5, col = "#34495E")
      if (length(active_sets) > 1) {
        lines(rep(i, 2), range(active_sets), 
              col = "#34495E", lwd = 6)
      }
    }
    inactive_sets <- which(patterns[i,] == 0)
    if (length(inactive_sets) > 0) {
      points(rep(i, length(inactive_sets)), inactive_sets, 
             pch = 1, cex = 5.5, col = "grey70", lwd = 3)
    }
  }
  
  # Add y-axis labels
  tryCatch({
    axis(2, at = 1:n_sets, labels = FALSE, las = 2, cex.axis = 2, family = "Arial", xpd = TRUE)
    text(x = par("usr")[1] - 0.02 * diff(par("usr")[1:2]),
         y = 1:n_sets, 
         labels = display_names, 
         adj = 1, 
         xpd = TRUE, 
         cex = 3.5,
         family = "Arial",
         srt = 0)
  }, error = function(e) {
    axis(2, at = 1:n_sets, labels = display_names, las = 2, cex.axis = 5.0, xpd = TRUE, family = "Arial")
  })
  
  dev.off()
  cat("Created UpSet plot:", filename, "\n")
}

create_set_size_plot <- function(asv_list, set_names, filename, title_text) {
  # Order using display names
  display_names_for_ordering <- convert_names(set_names)
  ordered_display <- order_groups(display_names_for_ordering, polyps_first = FALSE)  # Algae first
  
  ordered_set_names <- sapply(ordered_display, function(dn) {
    names(treatment_name_mapping)[treatment_name_mapping == dn][1]
  })
  ordered_set_names <- ordered_set_names[!is.na(ordered_set_names)]
  ordered_set_names <- ordered_set_names[sapply(ordered_set_names, function(x) length(asv_list[[x]]) > 0)]
  
  set_sizes <- sapply(asv_list[ordered_set_names], length)
  display_names <- convert_names(ordered_set_names)
  set_colors <- get_colors(display_names)
  
  png(filename, width = 38, height = 22, units = "in", res = 300)
  par(mar = c(5, 25, 4, 2), family = "Arial")
  
  bp <- barplot(set_sizes, names.arg = display_names, horiz = TRUE, 
                col = set_colors, border = NA,
                main = title_text,
                cex.main = 4.0,
                cex.names = 4.0,
                las = 1,
                xlim = c(0, max(set_sizes) * 1.15),
                axes = FALSE)
  
  mtext("Tot ASVs", side = 1, line = 3.5, cex = 3.5, family = "Arial")
  
  x_ticks <- pretty(c(0, max(set_sizes)))
  axis(1, at = x_ticks, labels = x_ticks, 
       cex.axis = 2, family = "Arial")
  
  text(set_sizes + max(set_sizes) * 0.02, bp, 
       labels = set_sizes, pos = 4, cex = 3.0, font = 2, family = "Arial")
  
  abline(v = seq(0, max(set_sizes), length.out = 5), 
         col = "grey90", lty = 2, lwd = 2)
  
  dev.off()
  cat("Created set size plot:", filename, "\n")
}

create_summary_barplots <- function(asv_list, set_names, filename_prefix) {
  # Order using display names
  display_names_for_ordering <- convert_names(set_names)
  ordered_display <- order_groups(display_names_for_ordering, polyps_first = FALSE)  # Algae first
  
  ordered_set_names <- sapply(ordered_display, function(dn) {
    names(treatment_name_mapping)[treatment_name_mapping == dn][1]
  })
  ordered_set_names <- ordered_set_names[!is.na(ordered_set_names)]
  ordered_set_names <- ordered_set_names[sapply(ordered_set_names, function(x) length(asv_list[[x]]) > 0)]
  
  total_asvs <- sapply(asv_list[ordered_set_names], length)
  
  if (all(total_asvs == 0)) {
    cat("Warning: No ASVs found\n")
    return(NULL)
  }
  
  display_names <- convert_names(ordered_set_names)
  names(total_asvs) <- display_names
  
  treatment_specific <- sapply(ordered_set_names, function(treat) {
    if (length(asv_list[[treat]]) == 0) return(0)
    other_sets <- setdiff(ordered_set_names, treat)
    current_asvs <- asv_list[[treat]]
    other_asvs <- unique(unlist(asv_list[other_sets]))
    length(setdiff(current_asvs, other_asvs))
  })
  names(treatment_specific) <- display_names
  
  set_colors <- get_colors(display_names)
  
  png(paste0(filename_prefix, "_Summary_Barplots.png"), 
      width = 24, height = 16, units = "in", res = 300)
  
  par(mfrow = c(1, 2), mar = c(18, 8, 4, 2), family = "Arial")
  
  # Total ASVs plot
  max_total <- max(total_asvs)
  if (max_total == 0) max_total <- 1
  
  bp1 <- barplot(total_asvs, 
                 col = set_colors,
                 border = NA,
                 main = "Total ASVs per Treatment",
                 cex.main = 3.2,
                 ylim = c(0, max_total * 1.15),
                 names.arg = rep("", length(display_names)),
                 axes = FALSE)
  
  mtext("Set size", side = 2, line = 4.5, cex = 3.0, family = "Arial")
  
  y_ticks <- pretty(c(0, max_total))
  axis(2, at = y_ticks, labels = y_ticks, 
       cex.axis = 2.5, family = "Arial", las = 1)
  
  text(bp1, -max_total * 0.02, labels = display_names,
       srt = 45, adj = 1, xpd = TRUE, cex = 3.2, family = "Arial")
  
  text(bp1, total_asvs + max_total * 0.02, 
       labels = total_asvs, pos = 3, cex = 2.5, font = 2, family = "Arial")
  
  # Treatment-specific ASVs plot
  max_specific <- max(treatment_specific)
  if (max_specific == 0) max_specific <- 1
  
  bp2 <- barplot(treatment_specific, 
                 col = set_colors,
                 border = NA,
                 main = "Treatment-Specific ASVs",
                 cex.main = 3.2,
                 ylim = c(0, max_specific * 1.15),
                 names.arg = rep("", length(display_names)),
                 axes = FALSE)
  
  mtext("Unique ASVs", side = 2, line = 4.5, cex = 3.0, family = "Arial")
  
  y_ticks2 <- pretty(c(0, max_specific))
  axis(2, at = y_ticks2, labels = y_ticks2, 
       cex.axis = 2.5, family = "Arial", las = 1)
  
  text(bp2, -max_specific * 0.02, labels = display_names,
       srt = 45, adj = 1, xpd = TRUE, cex = 3.2, family = "Arial")
  
  text(bp2, treatment_specific + max_specific * 0.02, 
       labels = treatment_specific, pos = 3, cex = 2.5, font = 2, family = "Arial")
  
  dev.off()
  cat("Created summary bar plots:", paste0(filename_prefix, "_Summary_Barplots.png"), "\n")
}

# ===============================================================================
# Create all plots
# ===============================================================================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("CREATING PLOTS\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# 1. ALGAE-ONLY ANALYSIS
if (length(algae_sets_original) >= 2) {
  cat("\n--- Creating plots: ALGAE-ONLY ---\n")
  
  create_custom_upset_plot(
    asv_list,
    set_names = algae_sets_original,
    filename = file.path(output_dir, paste0(output_prefix, "_AlgaeOnly_NoWater_UpSet.png")),
    title_text = "Bacterial ASV Overlap: Algae-Associated Communities Only"
  )
  
  create_set_size_plot(
    asv_list,
    set_names = algae_sets_original,
    filename = paste0(output_prefix, "_AlgaeOnly_NoWater_SetSizes.png"),
    title_text = "Set Sizes: Algae-Associated Communities Only"
  )
  
  create_summary_barplots(
    asv_list,
    set_names = algae_sets_original,
    filename_prefix = paste0(output_prefix, "_AlgaeOnly_NoWater")
  )
}

# 2. Polyp-only analysis
if (length(polyp_sets_original) >= 2) {
  cat("\n--- Creating plots: Polyp-Only ---\n")
  
  create_custom_upset_plot(
    asv_list,
    set_names = polyp_sets_original,
    filename = file.path(output_dir, paste0(output_prefix, "_PolypOnly_NoWater_UpSet.png")),
    title_text = "Bacterial ASV Overlap: Polyp-Associated Communities Only"
  )
  
  create_set_size_plot(
    asv_list,
    set_names = polyp_sets_original,
    filename = paste0(output_prefix, "_PolypOnly_NoWater_SetSizes.png"),
    title_text = "Set Sizes: Polyp-Associated Communities Only"
  )
  
  create_summary_barplots(
    asv_list,
    set_names = polyp_sets_original,
    filename_prefix = paste0(output_prefix, "_PolypOnly_NoWater")
  )
}

# 3. Combined Algae + Polyp
if (length(combined_algae_polyp) >= 2) {
  cat("\n--- Creating plots: Algae + Polyp Combined ---\n")
  
  create_custom_upset_plot(
    asv_list,
    set_names = combined_algae_polyp,
    filename = file.path(output_dir, paste0(output_prefix, "_AlgaePolyp_NoWater_UpSet.png")),
    title_text = "Bacterial ASV Overlap: Algae + Polyp Communities"
  )
  
  create_set_size_plot(
    asv_list,
    set_names = combined_algae_polyp,
    filename = paste0(output_prefix, "_AlgaePolyp_NoWater_SetSizes.png"),
    title_text = "Set Sizes: Algae + Polyp Communities"
  )
  
  create_summary_barplots(
    asv_list,
    set_names = combined_algae_polyp,
    filename_prefix = paste0(output_prefix, "_AlgaePolyp_NoWater")
  )
}

# 4. All treatments 
if (length(groups_no_water) >= 2) {
  cat("\n--- Creating plots: All Treatments ---\n")
  
  create_custom_upset_plot(
    asv_list,
    set_names = groups_no_water,
    filename = file.path(output_dir, paste0(output_prefix, "_All_NoWater_UpSet.png")),
    title_text = "Bacterial ASV Overlap Across All Treatments"
  )
  
  create_set_size_plot(
    asv_list,
    set_names = groups_no_water,
    filename = paste0(output_prefix, "_All_NoWater_SetSizes.png"),
    title_text = "Set Sizes: All Treatments"
  )
  
  create_summary_barplots(
    asv_list,
    set_names = groups_no_water,
    filename_prefix = paste0(output_prefix, "_NoWater")
  )
}

