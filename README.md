# Cassiopea Microbiome

This repository contains the analysis code for the Cassiopea microbiome project, focused on characterizing bacterial community variation across host and algal treatment groups associated with different symbiotic and developmental outcomes.

The workflow includes preprocessing of raw 16S rRNA amplicon data, taxonomic assignment, contaminant removal, alpha and beta diversity analyses, taxonomic composition analyses, differential abundance testing, and figure generation.

## Workflow

### 01_pre-processing

Preprocessing and core microbiome pipeline scripts, including sequence filtering, ASV inference, taxonomy assignment, decontamination, and creation of phyloseq objects.

### 02_alpha_diversity

Scripts for calculation and statistical analysis of alpha diversity metrics such as Observed richness, Fisher diversity, Shannon diversity, and Inverse Simpson diversity.

### 03_beta_diversity

Scripts for multivariate community analyses, including Bray–Curtis and phylogeny-based distance analyses, ordination, PERMANOVA, and related visualizations.

### 04_taxonomic_composition

Scripts for relative abundance summaries, taxonomic composition plots, and heat-tree visualization of abundant bacterial taxa.

### 05_differential_abundance

Scripts for DESeq2-based differential abundance analyses, heatmaps, and ASV overlap visualizations such as UpSet plots.

## Data availability

Raw sequencing data are deposited under NCBI BioProject PRJNA1390206.

Code archive DOI: [10.5281/zenodo.19553144](https://doi.org/10.5281/zenodo.19553144)

