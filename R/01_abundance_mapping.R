# ============================================================================
# Pedinophyceae Biogeography Analysis - Abundance Map Visualization
# ============================================================================
# Script for plotting relative abundance of Pedinophyceae across global samples
# Uses metabarcoding data from EukBank

# Session Info and Dependencies =============================================
library(tidyverse)
library(ggplot2)
library(maps)
library(here)

# Print session info for reproducibility
cat("Script executed on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sessionInfo()

# ============================================================================
# Configuration
# ============================================================================

# Set data directories (using here() for reproducibility)
data_dir <- here("data")
output_dir <- here("plots")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Data Loading and Preparation
# ============================================================================

prepare_abundance_data <- function(abundance_file, metadata_file) {
  """
  Load and prepare abundance and metadata files
  
  Args:
    abundance_file: Path to CSV file with abundance data
    metadata_file: Path to TSV file with metadata
  
  Returns:
    Tibble with combined and filtered data
  """
  
  # Load data files
  cat("Loading abundance data from:", abundance_file, "\n")
abundance_table <- read.table(abundance_file, 
                                sep = ",", 
                                header = TRUE, 
                                dec = ".",
                                stringsAsFactors = FALSE)
  
  cat("Loading metadata from:", metadata_file, "\n")
  metadata_all <- read.delim(metadata_file, 
                             header = TRUE, 
                             dec = ".",
                             stringsAsFactors = FALSE)
  
  # Combine datasets
  cat("Merging abundance and metadata...\n")
pedinoAbunMeta <- full_join(abundance_table, metadata_all, by = "sample")
  
  # Prepare and filter data
  pedinoAbunMetaFilter <- pedinoAbunMeta %>%
    select(sample, longitude, latitude, abundance_percent, amplicon) %>%
    # Create marker for samples with/without abundance data
    mutate(marker = ifelse(is.na(abundance_percent), "absent", "present"),
           # Ensure numeric types
           longitude = as.numeric(longitude),
           latitude = as.numeric(latitude),
           abundance_percent = as.numeric(abundance_percent)) %>%
    # Remove rows with missing coordinates
    filter(!is.na(longitude) & !is.na(latitude))
  
  cat("Data preparation complete.\n")
  cat("  Total samples with coordinates:", nrow(pedinoAbunMetaFilter), "\n")
  cat("  Samples with Pedinophyceae:", 
      sum(pedinoAbunMetaFilter$marker == "present", na.rm = TRUE), "\n")
  
  return(pedinoAbunMetaFilter)
}

# ============================================================================
# Visualization Function
# ============================================================================

create_abundance_map <- function(data, title = "Pedinophyceae Global Distribution") {
  """
  Create a world map showing Pedinophyceae abundance across samples
  
  Args:
    data: Tibble with longitude, latitude, and abundance_percent columns
    title: Title for the plot
  
  Returns:
    ggplot2 object
  """
  
  # Separate data for plotting (samples with and without abundance)
  data_present <- data %>% filter(marker == "present")
data_absent <- data %>% filter(marker == "absent")
  
  # Create the map
  map <- ggplot() +
    # Base world map
    borders("world", 
            colour = "gray70", 
            fill = "gray90",
            size = 0.3) +
    
    # Points for samples without detected Pedinophyceae (X marks)
    geom_point(data = data_absent, 
               aes(x = longitude, y = latitude), 
               shape = 4, 
               size = 1.5, 
               color = "gray50",
               alpha = 0.6,
               na.rm = TRUE) +
    
    # Points for samples with detected Pedinophyceae (sized by abundance)
    geom_point(data = data_present, 
               aes(x = longitude, y = latitude, 
                   size = abundance_percent, 
                   colour = abundance_percent), 
               alpha = 0.6,
               na.rm = TRUE) +
    
    # Size scale for abundance
    scale_size_continuous(name = "Abundance (%)",
                          range = c(1, 8),
                          breaks = c(0.1, 0.5, 1, 5, 10),
                          labels = c("0.1%", "0.5%", "1%", "5%", "10%")) +
    
    # Color scale for abundance (viridis for colorblind-friendly visualization)
    scale_color_viridis_c(name = "Abundance (%)",
                          option = "viridis",
                          breaks = c(0.1, 0.5, 1, 5, 10),
                          labels = c("0.1%", "0.5%", "1%", "5%", "10%")) +
    
    # Coordinate limits and breaks
    scale_y_continuous(limits = c(-90, 90), 
                       breaks = seq(-90, 90, by = 30)) +
    scale_x_continuous(breaks = seq(-180, 180, by = 60)) +
    
    # Use quickmap projection
    coord_quickmap() +
    
    # Minimal theme
    theme_minimal() +
    
    # Labels and title
    labs(x = NULL, 
         y = NULL, 
         title = title,
         subtitle = "Based on 18S rRNA V4 metabarcoding data from EukBank",
         caption = "X marks: samples where Pedinophyceae were not detected") +
    
    # Theme customization
    theme(
      axis.text = element_text(size = 7, color = "gray40"),
      axis.ticks = element_line(color = "gray70", size = 0.3),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8, face = "bold"),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray50"),
      plot.caption = element_text(size = 6, color = "gray60", hjust = 1),
      panel.grid.major = element_line(color = "gray95", size = 0.2),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(map)
}

# ============================================================================
# Main Execution
# ============================================================================

# Define file paths
abundance_file <- here("data", "raw", "pedinos_filtered_abundance.csv")
metadata_file <- here("data", "metadata", "eukbank_18S_V4_samples.tsv")

# Check if files exist
if (!file.exists(abundance_file)) {
  stop("Abundance file not found at:", abundance_file)
}
if (!file.exists(metadata_file)) {
  stop("Metadata file not found at:", metadata_file)
}

# Prepare data
pedinoAbunMetaFilter <- prepare_abundance_data(abundance_file, metadata_file)

# Create and display the map
pedinoMap <- create_abundance_map(pedinoAbunMetaFilter)
print(pedinoMap)

# Save the plot
output_path <- file.path(output_dir, "pedino_worldmap.pdf")
cat("Saving plot to:", output_path, "\n")

ggsave(plot = pedinoMap, 
       filename = output_path, 
       dpi = 300, 
       height = 6, 
       width = 11, 
       bg = "white",
       device = "pdf")

cat("✓ Plot saved successfully!\n")

# ============================================================================
# Optional: Summary Statistics
# ============================================================================

summary_stats <- pedinoAbunMetaFilter %>%
  filter(marker == "present") %>%
  summarise(
    n_samples = n(),
    mean_abundance = mean(abundance_percent, na.rm = TRUE),
    median_abundance = median(abundance_percent, na.rm = TRUE),
    max_abundance = max(abundance_percent, na.rm = TRUE),
    min_abundance = min(abundance_percent, na.rm = TRUE),
    sd_abundance = sd(abundance_percent, na.rm = TRUE)
  )

cat("\n=== Summary Statistics for Pedinophyceae ===\n")
print(summary_stats)