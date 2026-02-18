# Habitat Analysis Script

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(here)

# Function to create environment plots
create_env_plot <- function(data) {
  ggplot(data, aes(x = environment_variable, y = dependent_variable)) + 
    geom_point() + 
    theme_minimal() + 
    labs(title = "Environment Plot")
}

# Function to create habitat violin plots
create_violin_plots <- function(data) {
  ggplot(data, aes(x = habitat, y = nreads)) + 
    geom_violin() + 
    theme_minimal() + 
    labs(title = "Violin Plot of nreads")

  ggplot(data, aes(x = habitat, y = relative_abundance)) + 
    geom_violin() + 
    theme_minimal() + 
    labs(title = "Violin Plot of Relative Abundance")
}

# Function for Kruskal-Wallis test
run_kruskal_test <- function(data){
  kruskal.test(dependent_variable ~ habitat, data = data)
}

# Function for summary statistics
get_summary_stats <- function(data){
  data %>% 
    group_by(habitat) %>% 
    summarise(
      mean_nreads = mean(nreads),
      mean_relative_abundance = mean(relative_abundance),
      .groups = 'drop'
    )
}

# Load your data
# Make sure to adjust the path according to your directory structure
habitat_data <- read.csv(here("data", "habitat_data.csv"))

# Generate plots
env_plot <- create_env_plot(habitat_data)
violin_plots <- create_violin_plots(habitat_data)

# Run statistical tests
kruskal_result <- run_kruskal_test(habitat_data)

# Get summary statistics
summary_stats <- get_summary_stats(habitat_data)

# Print results
print(env_plot)
print(violin_plots)
print(kruskal_result)
print(summary_stats)