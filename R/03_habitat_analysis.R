# Load necessary packages
library(tidyverse)
library(here)
library(ggplot2)
library(patchwork)

# Function to load abundance data
load_abundance_data <- function(file_path) {
  read.csv(here(file_path))
}

# Function to load metadata
load_metadata <- function(file_path) {
  read.csv(here(file_path))
}

# Function to plot environment analysis (customize as needed)
plot_environment <- function(data) {
  # Normalize data if necessary
  normalized_data <- scale(data)
  ggplot(normalized_data, aes(x = factor(Habitat), y = Abundance)) +
    geom_violin() +
    theme_minimal() +
    labs(title = "Environment Plot Analysis")
}

# Function for habitat abundance visualizations
plot_habitat_abundance <- function(data) {
  ggplot(data, aes(x = factor(Habitat), y = Abundance)) +
    geom_violin() +
    theme_minimal() +
    labs(title = "Habitat Abundance Visualization")
}

# Statistical tests
perform_statistical_tests <- function(data) {
  kruskal.test(Abundance ~ Habitat, data = data)
}

# Summary statistics by habitat
summary_statistics_by_habitat <- function(data) {
  data %>%
    group_by(Habitat) %>%
    summarise(median_abundance = median(Abundance), sd_abundance = sd(Abundance))
}

# Directory structure using here::here()
data_dir <- here("data")
raw_data_dir <- here("data/raw")
metadata_dir <- here("data/metadata")
plots_dir <- here("plots")
results_dir <- here("results")

# Ensure directories exist
dir.create(data_dir, showWarnings = FALSE)
dir.create(raw_data_dir, showWarnings = FALSE)
dir.create(metadata_dir, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)

# Main analysis function
main_analysis <- function() {
  abundance_data <- load_abundance_data("data/raw/abundance.csv")
  metadata <- load_metadata("data/metadata/metadata.csv")
  
  # Perform analysis
  plot_env <- plot_environment(abundance_data)
  plot_habitat <- plot_habitat_abundance(abundance_data)
  stat_results <- perform_statistical_tests(abundance_data)
  summary_stats <- summary_statistics_by_habitat(abundance_data)
  
  # Save plots
  ggsave(here(plots_dir, "environment_plot.png"), plot_env)
  ggsave(here(plots_dir, "habitat_abundance_plot.png"), plot_habitat)

  # Return results
  list(statistical_results = stat_results, summary_statistics = summary_stats)
}

# Run the analysis
results <- main_analysis()
