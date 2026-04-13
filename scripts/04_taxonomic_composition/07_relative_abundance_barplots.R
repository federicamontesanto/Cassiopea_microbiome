################################################################################
# Cassiopea microbiome project
# Barplots and summary tables for bacterial and fungal taxonomic composition across treatment subsets.
################################################################################

required_packages <- c(
  "phyloseq",
  "ggplot2",
  "dplyr",
  "tidyr",
  "RColorBrewer",
  "viridis",
  "scales"
)

# ---------------------------------------------------------------------------
# Path configuration
# ---------------------------------------------------------------------------
project_root <- "."
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "04_taxonomic_composition/relative_abundance_barplots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
bacteria_phyloseq_file <- file.path(data_dir, "bacteria_phyloseq_relaxed.rds")


# --- Treatment name mapping ---
treatment_name_mapping <- c(
  "Water" = "Water",
  "Algae_Smic_KB8" = "Algae_native",
  "Algae_Antibiotic" = "Algae_Antibiotic",
  "Algae_SSB01_Bmin" = "Algae_Control",
  "Algae_Mutant" = "Algae_Mutant",
  "Polyp_Aposymbiotic" = "Aposymbiotic",
  "Polyp_Smic_KB8" = "Polyp_Native",
  "Polyp_SSB01_Bmin" = "Polyp_Control",
  "Polyp_Antibiotic" = "Polyp_Antibiotic",
  "Polyp_Mutant" = "Polyp_Mutant"
)

# Function to convert original names to display names
convert_treatment_names <- function(original_names) {
  display_names <- sapply(original_names, function(name) {
    ifelse(name %in% names(treatment_name_mapping), treatment_name_mapping[name], name)
  })
  return(unname(display_names))
}

# Define color palettes
get_scientific_color_palette <- function(n_colors) {
  if (n_colors <= 8) {
    
    colors <- brewer.pal(max(3, min(8, n_colors)), "Dark2")
    return(colors[1:n_colors])
  } else if (n_colors <= 11) {
    
    colors <- brewer.pal(max(3, min(11, n_colors)), "Set3")
    return(colors[1:n_colors])
  } else if (n_colors <= 20) {
    
    return(viridis(n_colors - 1, option = "D", end = 0.9))
  } else {
   
    return(plasma(n_colors - 1, end = 0.9))
  }
}

plot_top_taxa_clean <- function(ps_obj, 
                                rank = "Family", 
                                top_n = 10, 
                                group_var = "Treatment",
                                color_palette = NULL,
                                add_separators = FALSE,
                                separator_positions = NULL) {
  
  # Agglomerate at taxonomic rank
  ps_glom <- tax_glom(ps_obj, taxrank = rank, NArm = TRUE)
  
  # Transform to relative abundance
  ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
  
  # Melt to long format
  df <- psmelt(ps_rel)
  
  df <- df %>% filter(!grepl("Antifungal", Treatment))
  
  # Clean rank column
  df[[rank]] <- as.character(df[[rank]])
  df[[rank]][is.na(df[[rank]]) | df[[rank]] == ""] <- "Unclassified"
  
  # Shorten long taxa names for better readability
  df[[rank]] <- gsub("^(.{0,20}).*", "\\1", df[[rank]])  
  df[[rank]] <- gsub("_$", "", df[[rank]])
  
  # Calculate top taxa overall
  top_taxa <- df %>%
    group_by(.data[[rank]]) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    pull(.data[[rank]])
  
  # Collapse less abundant taxa into 'Other'
  df$Taxon <- ifelse(df[[rank]] %in% top_taxa, df[[rank]], "Other")
  
  # Order taxa by abundance (Other always last)
  taxa_order <- df %>%
    filter(Taxon != "Other") %>%
    group_by(Taxon) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>%
    pull(Taxon)
  
  df$Taxon <- factor(df$Taxon, levels = c(taxa_order, "Other"))
  
  # Convert treatment names to display names
  df$Treatment_Display <- convert_treatment_names(df[[group_var]])
  
  # Calculate mean and SE by treatment 
  df_summary <- df %>%
    group_by(Treatment_Display, Taxon) %>%
    summarise(
      MeanAbundance = mean(Abundance, na.rm = TRUE),
      SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Order treatments for logical grouping using display names
  treatment_levels_display <- unique(df_summary$Treatment_Display)
  
  # Create custom treatment order with specific user-requested order (excluding antifungal)
  water_treatments <- treatment_levels_display[grepl("Water", treatment_levels_display)]
  
  # Specify exact order for polyp and algae treatments 
  desired_order <- c(
    "Aposymbiotic",           
    "Polyp_Native",    
    "Polyp_Control",
    "Polyp_Antibiotic", 
    "Polyp_Mutant",     
    "Algae_Native",    
    "Algae_Control",
    "Algae_Antibiotic",
    "Algae_Mutant"      
  )
  
  # Filter to only include treatments present in the data
  present_treatments <- desired_order[desired_order %in% treatment_levels_display]
  
  # Add any other treatments 
  other_treatments <- treatment_levels_display[!treatment_levels_display %in% c(water_treatments, desired_order)]
  
    ordered_treatments <- c(water_treatments, present_treatments, other_treatments)
  df_summary$Treatment_Display <- factor(df_summary$Treatment_Display, levels = ordered_treatments)
  
  # Set color palette for taxa
  n_colors <- length(levels(df$Taxon))
  if (is.null(color_palette)) {
    color_palette <- get_scientific_color_palette(n_colors)
    # Ensure "Other" gets a gray color
    if ("Other" %in% levels(df$Taxon)) {
      color_palette[length(color_palette)] <- "#808080"  
    }
  }
  
  # Create manual legend labels to handle italics
  taxon_levels <- levels(df$Taxon)
  
  # Create custom labels for legend 
  taxon_levels <- levels(df$Taxon)
  custom_labels <- sapply(taxon_levels, function(x) {
    if (x == "Other") {
      return(x)  
    } else {
      return(bquote(italic(.(x))))  #
    }
  }, simplify = FALSE)
  
  # Create plot 
  p <- ggplot(df_summary, aes(x = Treatment_Display, y = MeanAbundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +  
    scale_fill_manual(values = color_palette, name = rank, labels = custom_labels) +
    labs(x = "Treatment", 
         y = "Relative abundance (%)") 
    theme_classic(base_size = 20, base_family = "Arial") +  
    theme(
      # Text styling with Arial and large fonts
      text = element_text(color = "black", family = "Arial"),
      
      # Axis styling with large fonts (REMOVED bold formatting)
      axis.title.x = element_text(size = 24, margin = margin(t = 15), family = "Arial"),  
      axis.title.y = element_text(size = 24, margin = margin(r = 15), family = "Arial"),  
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18, color = "black", family = "Arial"),   
      axis.text.y = element_text(size = 18, color = "black", family = "Arial"),
      legend.title = element_text(size = 20, family = "Arial", face = "plain"),
      legend.text = element_text(size = 16, family = "Arial"),  
      legend.position = "right",
      legend.key.size = unit(1.0, "cm"),  
      legend.spacing.y = unit(0.2, "cm"),
      legend.box.spacing = unit(0.5, "cm"),
      
      # Panel styling
      panel.border = element_rect(fill = NA, color = "black", size = 1.2),  
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", size = 1.0),  
      
      # Plot margins
      plot.margin = unit(c(0.5, 1.5, 0.5, 0.5), "cm")  
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    guides(fill = guide_legend(
      ncol = 1, 
      override.aes = list(color = "black", size = 0.5)
    ))
  
  # Add vertical separators 
  if (add_separators && !is.null(separator_positions)) {
    p <- p + 
      geom_vline(xintercept = separator_positions, 
                 linetype = "dashed", 
                 color = "gray30", 
                 size = 1.2,  
                 alpha = 0.8)
  }
  
  return(p)
}

# Load your data 
ps <- readRDS(bacteria_phyloseq_file)


# ===== 1. PLOTS WITHOUT WATER =====
ps_no_water <- subset_samples(ps, Treatment != "Water")
ps_no_water <- prune_taxa(taxa_sums(ps_no_water) > 0, ps_no_water)



separator_positions_no_water <- c(5.5)  

# Create plots without water with separators
cat("Creating plots without water...\n")
p_bac_fam_no_water <- plot_top_taxa_clean(ps_no_water, 
                                          rank = "Family", 
                                          top_n = 10,
                                          add_separators = TRUE,
                                          separator_positions = separator_positions_no_water)

p_bac_gen_no_water <- plot_top_taxa_clean(ps_no_water, 
                                          rank = "Genus", 
                                          top_n = 10,
                                          add_separators = TRUE,
                                          separator_positions = separator_positions_no_water)


# Save no water plots with larger dimensions for big fonts
ggsave(file.path(output_dir, "Bacteria_Family_NoWater.png"), plot = p_bac_fam_no_water, 
       width = 14, height = 10, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Bacteria_Genus_NoWater.png"), plot = p_bac_gen_no_water, 
       width = 14, height = 10, dpi = 600, bg = "white")


# ===== 2. ALGAE ONLY PLOTS =====
cat("Creating algae-only plots...\n")
# Get all treatments that start with "Algae_" 
algae_treatments <- unique(sample_data(ps)$Treatment)
algae_treatments <- algae_treatments[grepl("^Algae_", algae_treatments)]

ps_algae <- subset_samples(ps, Treatment %in% algae_treatments)
ps_algae <- prune_taxa(taxa_sums(ps_algae) > 0, ps_algae)


# Create algae-only plots (no separators)
p_bac_fam_algae <- plot_top_taxa_clean(ps_algae, rank = "Family", top_n = 10)
p_bac_gen_algae <- plot_top_taxa_clean(ps_algae, rank = "Genus", top_n = 10)


# Save algae plots
ggsave(file.path(output_dir, "Bacteria_Family_AlgaeOnly.png"), plot = p_bac_fam_algae, 
       width = 12, height = 10, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Bacteria_Genus_AlgaeOnly.png"), plot = p_bac_gen_algae, 
       width = 12, height = 10, dpi = 600, bg = "white")

# ===== 3. POLYPS ONLY PLOTS =====
cat("Creating polyps-only plots...\n")
# Get all treatments that start with "Polyp_" PLUS Aposymbiotic 
polyp_treatments <- unique(sample_data(ps)$Treatment)
polyp_treatments <- c(polyp_treatments[grepl("^Polyp_", polyp_treatments)], "Aposymbiotic")

ps_polyps <- subset_samples(ps, Treatment %in% polyp_treatments)
ps_polyps <- prune_taxa(taxa_sums(ps_polyps) > 0, ps_polyps)


# Create polyps-only plots (no separators)
p_bac_fam_polyps <- plot_top_taxa_clean(ps_polyps, rank = "Family", top_n = 10)
p_bac_gen_polyps <- plot_top_taxa_clean(ps_polyps, rank = "Genus", top_n = 10)

# Save polyps plots
ggsave(file.path(output_dir, "Bacteria_Family_PolypsOnly.png"), plot = p_bac_fam_polyps, 
       width = 12, height = 10, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Bacteria_Genus_PolypsOnly.png"), plot = p_bac_gen_polyps, 
       width = 12, height = 10, dpi = 600, bg = "white")

# ===== 4. ALL TREATMENTS =====
cat("Creating all treatments plots...\n")

separator_positions_all <- c(1.5, 5.5)

# Create plots with all treatments and separators
p_bac_fam_all <- plot_top_taxa_clean(ps, 
                                     rank = "Family", 
                                     top_n = 10,
                                     add_separators = TRUE,
                                     separator_positions = separator_positions_all)

p_bac_gen_all <- plot_top_taxa_clean(ps, 
                                     rank = "Genus", 
                                     top_n = 10,
                                     add_separators = TRUE,
                                     separator_positions = separator_positions_all)


# Save all treatment plots with largest dimensions
ggsave(file.path(output_dir, "Bacteria_Family_All.png"), plot = p_bac_fam_all, 
       width = 16, height = 10, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Bacteria_Genus_All.png"), plot = p_bac_gen_all, 
       width = 16, height = 10, dpi = 600, bg = "white")


# Create summary table of relative abundances with display names
cat("\n=== CREATING SUMMARY TABLES WITH DISPLAY NAMES ===\n")

# Function to create summary statistics with display names
create_summary_stats <- function(ps_obj, rank = "Family", top_n = 10) {
  # Agglomerate at taxonomic rank
  ps_glom <- tax_glom(ps_obj, taxrank = rank, NArm = TRUE)
  
  # Transform to relative abundance
  ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
  
  # Melt to long format
  df <- psmelt(ps_rel)
  
  # Clean rank column
  df[[rank]] <- as.character(df[[rank]])
  df[[rank]][is.na(df[[rank]]) | df[[rank]] == ""] <- "Unclassified"
  
  # Calculate top taxa
  top_taxa <- df %>%
    group_by(.data[[rank]]) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    pull(.data[[rank]])
  
  # Convert treatment names to display names
  df$Treatment_Display <- convert_treatment_names(df$Treatment)
  
  # Create summary by treatment 
  summary_stats <- df %>%
    filter(.data[[rank]] %in% top_taxa) %>%
    group_by(Treatment_Display, .data[[rank]]) %>%
    summarise(
      Mean = round(mean(Abundance, na.rm = TRUE), 2),
      SD = round(sd(Abundance, na.rm = TRUE), 2),
      SE = round(sd(Abundance, na.rm = TRUE) / sqrt(n()), 2),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(Treatment_Display, desc(Mean))
  
  return(summary_stats)
}

# Create summary tables
bac_fam_summary <- create_summary_stats(ps, rank = "Family", top_n = 10)
bac_gen_summary <- create_summary_stats(ps, rank = "Genus", top_n = 10)


# Save summary tables
write.csv(bac_fam_summary, file.path(output_dir, "Bacteria_Family_Summary_Stats_DisplayNames.csv"), row.names = FALSE)
write.csv(bac_gen_summary, file.path(output_dir, "Bacteria_Genus_Summary_Stats_DisplayNames.csv"), row.names = FALSE)





