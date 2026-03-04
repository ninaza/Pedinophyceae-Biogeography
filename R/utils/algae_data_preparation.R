library(tidyverse)
library(dplyr)

## prepare metadata
# assign categories
# Oceans
metadata <- read.delim("data/raw/eukbank_18S_V4_samples.tsv", header = T, dec = ".")
samples_longhurst <- read.table("data/edited/eukbank_18S_V4_samples_with_longhurst.csv", sep = ",", header = T, dec = ".")

oceans <- samples_longhurst %>%
  select(sample, name) %>%
  rename(ocean = name)

metadata <- inner_join(metadata, oceans)

# Categorize habitats
metadata <- metadata %>%
  mutate(habitat = case_when(
    envplot %in% c("marine_sediment", "marine_water", "marine_organism") ~ "Marine",
    envplot %in% c("land_freshwater") ~ "Freshwater",
    envplot %in% c("land_soil", "land_sediment", "land_organism", "land_water") ~ "Terrestrial",
    envplot %in% c("none") ~ "None",
  ))

## size fraction
metadata <- metadata %>%
  mutate(
    # Ensure limits are numeric
    size_fraction_lower_threshold = as.numeric(size_fraction_lower_threshold),
    size_fraction_upper_threshold = as.numeric(ifelse(is.na(size_fraction_upper_threshold), Inf, size_fraction_upper_threshold)),
    
    # Refined Categorization
    size_fraction = case_when(
      # Bulk samples (Capture everything from a small pore size up)
      size_fraction_lower_threshold <= 0.22 & size_fraction_upper_threshold == Inf ~ "Bulk_Total",
      
      # Specific Fractions
      size_fraction_lower_threshold >= 0.2  & size_fraction_upper_threshold <= 3   ~ "Pico",
      size_fraction_lower_threshold >= 3    & size_fraction_upper_threshold <= 20  ~ "Nano",
      size_fraction_lower_threshold >= 20   & size_fraction_upper_threshold < 200  ~ "Micro",
      size_fraction_lower_threshold >= 200   & size_fraction_upper_threshold < Inf  ~ "Meso",
      
      # Catch-all for overlapping or non-standard ranges
      TRUE ~ "Non_Standard"
    )
  )

## salinity
metadata <- metadata %>%
  mutate(salinity_cat = case_when(
    salinity <= 0.5 ~ "freshwater",
    salinity > 0.5 & salinity <= 30 ~ "brackish",
    salinity > 30 & salinity <= 50 ~ "saline",
    salinity > 50 ~ "briny",
    is.na(salinity) ~ "unknown"
  ))

## depth
metadata <- metadata %>%
  mutate(ocean_layer_fine = case_when(
    depth <= 200 & depth ~ "photic",
    depth > 200 & depth <= 1000 ~ "twilight",
    depth > 1000 ~ "aphotic",
    TRUE ~ NA_character_  # Assign NA if it doesn't fit into any category
  ))

metadata <- metadata %>%
  mutate(ocean_layer = case_when(
    depth <= 200 & depth ~ "sunlit",
    depth > 200 ~ "dark",
    TRUE ~ NA_character_  # Assign NA if it doesn't fit into any category
  ))

write.csv(metadata, file = "data/edited/eukbank_18S_V4_samples_edited.csv", row.names = FALSE)

## subsetting abundance table
abundance_table <- read.table("data/edited/pedinos_abundance.csv", sep = ",", header = T, dec = ".")

# filter out all amplicons appearing in two or less samples
abundance_table <- abundance_table %>% 
  group_by(amplicon) %>%
  filter(sum(sample > 0) > 2) %>%
  ungroup()

write.csv(abundance_table, file = "data/edited/pedinos_filtered_abundance.csv", row.names = FALSE)

# Ensure metadata and ASV table samples match
metadata <- read.table("data/edited/eukbank_18S_V4_samples_edited.csv", sep = ",", header = T, dec = ".")
common_samples <- intersect(abundance_table$sample, metadata$sample)
abundance_table <- abundance_table[abundance_table$sample %in% common_samples, ]
metadata <- metadata[metadata$sample %in% common_samples, ]
write.csv(metadata, file = "data/edited/pedinos_samples_edited.csv", row.names = FALSE)
