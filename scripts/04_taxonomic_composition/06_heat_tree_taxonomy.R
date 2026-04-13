################################################################################
# Cassiopea microbiome project
#   Heat-tree visualizations of abundant bacterial taxa with updated taxonomy mappings for algae and polyp treatment subsets.
#
################################################################################

# ---------------------------------------------------------------------------
# Path configuration
# ---------------------------------------------------------------------------
project_root <- "."
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "04_taxonomic_composition/heat_trees")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
phyloseq_bact_file <- file.path(data_dir, "bacteria_phyloseq_complete.rds")

required_packages <- c(
  "phyloseq",
  "metacoder",
  "ggplot2",
  "dplyr",
  "tidyr",
  "RColorBrewer"
)


# Treatment name mapping
treatment_name_mapping <- c(
  "Water" = "Water",
  "Algae_Smic_KB8" = "Algae_Native",
  "Algae_Antibiotic" = "Algae_Antibiotic",
  "Algae_SSB01_Bmin" = "Algae_Control",
  "Algae_Mutant" = "Algae_Mutant",
  "Polyp_Aposymbiotic" = "Aposymbiotic",
  "Polyp_Smic_KB8" = "Polyp_Native",
  "Polyp_SSB01_Bmin" = "Polyp_COntrol",
  "Polyp_Antibiotic" = "Polyp_Antibiotic",
  "Polyp_Mutant" = "Polyp_Mutant"
)

# ===== TAXONOMIC NAME MAPPING =====
# Polyp-specific updates
polyp_taxonomy_mapping <- c(
  "Exiguobacteriaceae" = "Bacillaceae",
  "Exiguobacterales" = "Exiguobacteriales",
  "Firmicutes" = "Bacillota",
  "Marinicaulis" = "Hyphococcus",
  "Muricauda" = "Flagellimonas",
  "Myxococcia" = "Deltaproteobacteria",
  "Myxococcota" = "Pseudomonadota",
  "Phormidesmiales" = "Phormidesmidales",
  "Proteobacteria" = "Pseudomonadota",
  "Rhizobiales" = "Hyphomicrobiales",
  "Rhodobacteraceae" = "Paracoccaceae",
  "Thalassobaculales" = "Rhodospirillales"
)

# Algae-specific updates
algae_taxonomy_mapping <- c(
  "Actinobacteria" = "Actinomycetota",
  "Actinobacteriota" = "Actinomycetota",
  "Labrenzia" = "Roseibium",
  "Marinobacteraceae" = "Alteromonadaceae",
  "Micrococcales" = "Micrococcales",
  "Rhizobiales" = "Hyphomicrobiales",
  "Thalassobaculales" = "Rhodospirillales"
)

# Function to update taxonomic names based on group type
update_taxonomy_name <- function(name, group_type) {
  if (tolower(group_type) == "polyp") {
    mapping <- polyp_taxonomy_mapping
  } else {
    mapping <- algae_taxonomy_mapping
  }
  
  # Return updated name if exists in mapping, otherwise return original
  if (name %in% names(mapping)) {
    return(unname(mapping[name]))
  } else {
    return(name)
  }
}

convert_treatment_names <- function(x) {
  unname(ifelse(x %in% names(treatment_name_mapping), treatment_name_mapping[x], x))
}

# WORKING version - uses n_obs method 
create_working_heat_tree <- function(ps_obj, group_type, min_abundance = 5.0, 
                                     output_dir, max_taxa = 12, save_pdf = TRUE) {
  

  
  tryCatch({
    # Get sample information
    sdata <- as.data.frame(sample_data(ps_obj))
    sdata$Treatment_Display <- convert_treatment_names(sdata$Treatment)
    
    # Select group samples 
    if (tolower(group_type) == "algae") {
      group_samples <- rownames(sdata)[grepl("^Algae_", sdata$Treatment_Display)]
      color_scheme <- c("lightgreen", "mediumseagreen", "darkgreen", "forestgreen", "darkslategray")
    } else {
      group_samples <- rownames(sdata)[grepl("^Polyp_|^Apo$", sdata$Treatment_Display)]
      color_scheme <- c("lightpink", "hotpink", "darkorange", "darkorchid", "darkviolet")
    }
    
    if(length(group_samples) == 0) {
      return(NULL)
    }
    
    
    # Transform to relative abundance
    ps_rel <- transform_sample_counts(ps_obj, function(x) x / sum(x) * 100)
    
    # Filter to group samples
    ps_group <- prune_samples(sample_names(ps_rel) %in% group_samples, ps_rel)
    ps_group <- prune_taxa(taxa_sums(ps_group) > 0, ps_group)
    
    # Get taxonomy table
    tax_df <- as.data.frame(as.matrix(tax_table(ps_group)))
    
    # Calculate abundance for each taxon
    abundance_data <- data.frame(
      taxon_id = taxa_names(ps_group),
      abundance = taxa_sums(ps_group),
      stringsAsFactors = FALSE
    )
    
    # Remove NAs and zero abundance
    abundance_data <- abundance_data[!is.na(abundance_data$abundance), ]
    abundance_data <- abundance_data[abundance_data$abundance > 0, ]
    
    # Filter by minimum abundance and select top taxa
    abundant_taxa_df <- abundance_data[abundance_data$abundance >= min_abundance, ]
    abundant_taxa_df <- abundant_taxa_df[order(abundant_taxa_df$abundance, decreasing = TRUE), ]
    
    if(nrow(abundant_taxa_df) == 0) {
      return(NULL)
    }
    
    if(nrow(abundant_taxa_df) > max_taxa * 2) {  
      abundant_taxa_df <- abundant_taxa_df[1:(max_taxa * 2), ]
    }
    
   
    # Get taxonomy for abundant taxa
    tax_clean <- tax_df[abundant_taxa_df$taxon_id, , drop = FALSE]
    
    # Create full taxonomic hierarchy 
   
    lineage_data <- data.frame(
      lineage = character(nrow(tax_clean)),
      abundance = abundant_taxa_df$abundance,
      otu_count = rep(1, nrow(tax_clean)),
      clean_name = character(nrow(tax_clean)),
      stringsAsFactors = FALSE
    )
    
    for(i in 1:nrow(tax_clean)) {
      tax_row <- tax_clean[i, ]
      
      # Build complete taxonomic hierarchy
      hierarchy <- c()
      clean_names <- c()
      
      ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      
      for(rank in ranks) {
        if(rank %in% colnames(tax_row)) {
          value <- as.character(tax_row[[rank]])
          if(!is.na(value) && value != "" && value != "NA" && 
             !grepl("^__$|^NA$|^unknown$|^unidentified$", value, ignore.case = TRUE)) {
            
            # Clean the taxonomic name
            clean_value <- gsub("^[a-z]__", "", value)
            clean_value <- gsub("_", " ", clean_value)
            clean_value <- trimws(clean_value)
            
            # ===== APPLY TAXONOMIC NAME UPDATES HERE =====
            clean_value <- update_taxonomy_name(clean_value, group_type)
            
            if(nchar(clean_value) > 0) {
              hierarchy <- c(hierarchy, clean_value)
              clean_names <- c(clean_names, clean_value)
            }
          } else {
            placeholder <- paste0("Unknown_", rank)
            hierarchy <- c(hierarchy, placeholder)
            clean_names <- c(clean_names, placeholder)
          }
        } else {
          placeholder <- paste0("Unknown_", rank)
          hierarchy <- c(hierarchy, placeholder)
          clean_names <- c(clean_names, placeholder)
        }
      }
      
      if(length(hierarchy) < 3) {
        hierarchy <- c("Bacteria", "Unknown_Phylum", paste0("Taxon_", i))
        clean_names <- hierarchy
      }
      
      lineage_data$lineage[i] <- paste(hierarchy, collapse = ";")
      lineage_data$clean_name[i] <- tail(clean_names, 1)
    }
    
    # Show the full taxonomic hierarchy we extracted
    
    for(i in 1:min(3, nrow(lineage_data))) {
      cat("   ", i, ":", lineage_data$clean_name[i], "\n")
      cat("      Full path:", lineage_data$lineage[i], "\n")
      cat("      Abundance:", round(lineage_data$abundance[i], 2), "%\n\n")
    }
    
   
    # Create taxmap object
    
    
    obj <- parse_tax_data(lineage_data,
                          class_cols = "lineage", 
                          class_sep = ";")
    
    if(is.null(obj)) {
      return(NULL)
    }
    
    
    # Add the abundance and OTU data to the taxmap
    obj$data$tax_data$bacterial_name <- lineage_data$clean_name
    obj$data$tax_data$rel_abundance <- lineage_data$abundance
    obj$data$tax_data$otu_count <- lineage_data$otu_count
    
    # Calculate abundance for intermediate nodes
    tryCatch({
      obj <- obj %>% 
        mutate_obs("rel_abundance_sum", 
                   sum(rel_abundance, na.rm = TRUE),
                   groups = supertaxa, na.rm = TRUE)
      
      obj <- obj %>% 
        mutate_obs("otu_count_sum", 
                   sum(otu_count, na.rm = TRUE), 
                   groups = supertaxa, na.rm = TRUE)
      
        }, error = function(e) {
      obj$data$tax_data$rel_abundance_sum <- obj$data$tax_data$rel_abundance
      obj$data$tax_data$otu_count_sum <- obj$data$tax_data$otu_count
    })
    
    # Add label for root node
    root_taxon_id <- names(obj$taxon_names())[obj$taxon_names() == "Bacteria" | 
                                                grepl("^Bacteria$|^Unknown_Kingdom", obj$taxon_names())]
    if(length(root_taxon_id) > 0) {
      obj$data$tax_data$bacterial_name[obj$data$tax_data$taxon_id == root_taxon_id[1]] <- "Bacteria"
    }
    
    # Create heat trees
    set.seed(42)
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    results <- list()
    
    
    
    #  HIERARCHICAL heat tree
    if(!is.null(results$labeled)) {
        hierarchical_result <- tryCatch({
        heat_tree(obj,
                  node_label = taxon_names,
                  node_size = n_obs,
                  node_color = n_obs,
                  node_color_range = color_scheme,
                  node_color_axis_label = "Number of Descendants",
                  node_size_axis_label = "Number of Descendants",
                  layout = "davidson-harel",
                  initial_layout = "reingold-tilford",
                  title = paste(group_type, "Microbiome - Full Taxonomic Hierarchy (Updated)"),
                  node_size_range = c(0.015, 0.05),
                  node_label_size_range = c(0.02, 0.035),
                  edge_color = "gray40",
                  edge_size_range = c(0.001, 0.003),
                  make_node_legend = TRUE,
                  make_edge_legend = FALSE,
                  node_label_color = "black",
                  overlap_avoidance = 0.9,
                  margin_size = c(0.1, 0.1, 0.1, 0.15),
                  output_file = NULL)
      }, error = function(e) {
         return(NULL)
      })
      
      if(!is.null(hierarchical_result)) {
        hierarchical_file <- file.path(output_dir, paste0("IMPROVED_HIERARCHICAL_Heat_Tree_", group_type, "_", timestamp, ".png"))
        ggsave(hierarchical_file, hierarchical_result, width = 24, height = 20, dpi = 300, bg = "white")
        results$hierarchical <- hierarchical_result
        
        if(save_pdf) {
          hierarchical_pdf <- file.path(output_dir, paste0("IMPROVED_HIERARCHICAL_Heat_Tree_", group_type, "_", timestamp, ".pdf"))
          ggsave(hierarchical_pdf, hierarchical_result, width = 24, height = 20, dpi = 300, bg = "white")
          }
      }
    }
    
    # Return best available result
    if(!is.null(results$hierarchical)) {
      return(list(plot = results$hierarchical, obj = obj, data = lineage_data, method = "hierarchical"))
    } else if(!is.null(results$labeled)) {
      return(list(plot = results$labeled, obj = obj, data = lineage_data, method = "labeled"))
    } else if(!is.null(results$basic)) {
      return(list(plot = results$basic, obj = obj, data = lineage_data, method = "basic"))
    } else {
      return(NULL)
    }
    
  }, error = function(main_error) {
    return(NULL)
  })
}

# Load data
if(!exists("ps")) {
  if (file.exists(phyloseq_bact_file)) {
    ps <- readRDS(phyloseq_bact_file)
    
  } else {
    stop("Cannot find `data/bacteria_phyloseq_complete.rds`.")
  }
}



# Heat trees with updated taxonomy

algae_working <- create_working_heat_tree(ps, "Algae", min_abundance = 5.0, output_dir, max_taxa = 12)


polyp_working <- create_working_heat_tree(ps, "Polyp", min_abundance = 5.0, output_dir, max_taxa = 12)

# Try lower thresholds if needed
if(is.null(algae_working)) {
  
  algae_working <- create_working_heat_tree(ps, "Algae", min_abundance = 5.0, output_dir, max_taxa = 8)
}

if(is.null(polyp_working)) {
  
  polyp_working <- create_working_heat_tree(ps, "Polyp", min_abundance = 5.0, output_dir, max_taxa = 8)
}

# Final results

working_files <- list.files(output_dir, pattern = "(WORKING|LABELED_WORKING|HIERARCHICAL|IMPROVED_HIERARCHICAL)_Heat_Tree.*\\.(png|pdf)$")


if(length(working_files) > 0) {
    for(file in working_files) {
    file_size <- file.size(file.path(output_dir, file))
    cat("   •", file, "(", round(file_size/1024), "KB )\n")
  }
}

