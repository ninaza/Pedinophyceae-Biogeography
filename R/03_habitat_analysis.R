# Load required libraries
library(ggplot2)
library(dplyr)
library(here)

# Load data
load_data <- function() {
  tryCatch({
    abundance_table <- read.csv(here("pedinos_filtered_abundance.csv"))
    metadata <- read.csv(here("pedinos_edit_metadata.csv"))
    metadata_all <- read.csv(here("eukbank_18S_V4_samples.tsv"), sep = "\t")
    
    return(list(abundance = abundance_table, metadata = metadata, metadata_all = metadata_all))
  }, error = function(e) {
    stop("Error loading data: ", e$message)
  })
}

# Plot environment analysis
plot_environment <- function(data) {
  tryCatch({
    ggplot(data, aes(x = factor(variable), y = value)) +
      geom_bar(stat = "identity") +
      labs(title = "Environment Plot Analysis") +
      theme_minimal()
  }, error = function(e) {
    stop("Error in environment plot: ", e$message)
  })
}

# Habitat abundance analysis
plot_habitat_abundance <- function(data) {
  tryCatch({
    ggplot(data, aes(x = habitat, y = nreads)) +
      geom_violin() +
      labs(title = "Habitat Abundance Analysis (Nreads)") +
      theme_minimal()
    ggplot(data, aes(x = habitat, y = relative_abundance)) +
      geom_violin() +
      labs(title = "Habitat Abundance Analysis (Relative Abundance)") +
      theme_minimal()
  }, error = function(e) {
    stop("Error in habitat abundance plot: ", e$message)
  })
}

# Kruskal-Wallis test
perform_kruskal_test <- function(data) {
  tryCatch({
    kruskal.test(nreads ~ habitat, data = data)
    kruskal.test(relative_abundance ~ habitat, data = data)
  }, error = function(e) {
    stop("Error performing Kruskal-Wallis test: ", e$message)
  })
}

# Main analysis function
run_analysis <- function() {
  data <- load_data()
  plot_environment(data$abundance)
  plot_habitat_abundance(data$abundance)
  perform_kruskal_test(data$abundance)
}

# Execute the analysis
run_analysis()