library(tidyverse)
library(ggplot2)
library(dplyr)
library(vegan)
library(UpSetR)
library(MetBrewer)

palette <- met.brewer("Renoir")

## PREP DATA
abundance_table <- read.table("data/edited/pedinos_abundance.csv", 
                              sep = ",", header = T, dec = ".")
metadata <- read.table("data/edited/pedinos_edit_filtered_metadata.csv", 
                       sep = ",", header = T, dec = ".")
## SIZE FRACTION
join_marine <- inner_join(metadata, abundance_table) %>%
  filter(habitat == "Marine")

# plot nreads per size fraction
size_violin <- ggplot(join_marine, aes(x = size_fraction, y = nreadsPedino)) +
  geom_violin(show.legend=FALSE, aes(fill = size_fraction), alpha = 0.5) +
  geom_jitter(show.legend=FALSE, width = 0.1, alpha = 0.4) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  #scale_fill_manual(values = palette[c(6,7,8)]) +
  theme_minimal() +
  labs(x=NULL,
       y = "nreads Pedinophyceae")
size_violin

ggsave(plot = size_violin, 
       filename = "plots/pedino_nreads_size_violin.png", 
       dpi = 300, height = 4, width = 6, bg = "white")

## plot size range per sample
size_range <- ggplot(join_marine, aes(y = sample)) +
  geom_segment(aes(x = size_fraction_lower_threshold, 
                   xend = ifelse(size_fraction_upper_threshold == Inf, 200, size_fraction_upper_threshold), 
                   yend = sample, color = size_fraction), size = 2) +
  scale_x_log10() +
  labs(x = "Size Pore (µm, log scale)", y = NULL) +
  theme_minimal()
size_range

ggsave(plot = size_range, 
       filename = "plots/pedino_sample_sizerange.png", 
       dpi = 300, height = 12, width = 16, bg = "white")

## filter only strict size fractions and check ASV overlap and distribution over size fractions
target_fractions <- c("Pico", "Nano", "Micro")

size_filtered <- join_marine %>%
  filter(size_fraction %in% target_fractions)

n_asvs_total <- n_distinct(join_marine$amplicon)
n_asvs_filtered <- n_distinct(size_filtered$amplicon)

## Upset plot of ASV overlap across size fractions
size_filtered <- size_filtered %>%
  filter(nreadsPedino > 0) # Keep only detections

# Determine presence per Size Fraction
# An ASV is 'present' in a fraction if it appears in at least one sample of that fraction
upset_input <- size_filtered %>%
  group_by(amplicon, size_fraction) %>%
  summarise(Present = 1, .groups = "drop") %>%
  pivot_wider(names_from = size_fraction, 
              values_from = Present, 
              values_fill = 0) %>%
  as.data.frame()

# Remove the ID column for the upset function
rownames(upset_input) <- upset_input$amplicon
upset_input_final <- upset_input[,-1]

# Create UpSet plot
png("plots/pedinos_size_upset.png", width = 2400, height = 1800, res = 300)

upset(upset_input_final, 
      sets = colnames(upset_input_final), 
      main.bar.color = "steelblue", 
      sets.bar.color = "darkgray",
      order.by = "freq", 
      empty.intersections = "on")

dev.off()

## Ternary plot of ASV distribution across size fractions
# Prepare a matrix of Relative Abundance per Size Fraction
matrix_data <- size_filtered %>%
  group_by(amplicon, size_fraction) %>%
  summarise(Sum_Reads = sum(nreadsPedino), .groups = "drop") %>%
  group_by(amplicon) %>%
  mutate(Rel_Abund = Sum_Reads / sum(Sum_Reads)) %>% # Normalize per ASV to see "preference"
  select(-Sum_Reads) %>%
  pivot_wider(names_from = size_fraction, values_from = Rel_Abund, values_fill = 0) %>%
  column_to_rownames("amplicon") %>%
  as.matrix()

ternary_df <- as.data.frame(matrix_data)

# Join total reads back to your ternary data
total_reads <- size_filtered %>% 
  group_by(amplicon) %>% 
  summarise(Total_Abundance = log10(sum(nreadsPedino)))

ternary_df_with_size <- ternary_df %>%
  rownames_to_column("amplicon") %>%
  left_join(total_reads, by = "amplicon")

ternary_plot <- ggtern(data = ternary_df_with_size, aes(x = Pico, y = Nano, z = Micro)) +
  geom_point(aes(size = Total_Abundance), alpha = 0.4, color = "darkblue") +
  theme_rgbw()
ternary_plot

ggsave(plot = ternary_plot, 
       filename = "plots/pedino_size_ternary.png", 
       dpi = 300, height = 4, width = 6, bg = "white")

## investigate the phylogenetic placement of miro-leaning ASVs
# Define a threshold for 'Micro-associated'
# For example, ASVs where more than 30% of their reads are in the Micro fraction
potential_symbionts <- ternary_df %>%
  rownames_to_column("amplicon") %>%
  filter(Micro > 0.30) %>% 
  arrange(desc(Micro))
potential_symbionts

## PERMANOVA to test if ASV distribution across size fractions is significant
# Pivot to wide format (Community Matrix)
community_matrix <- size_filtered %>%
  group_by(sample, amplicon) %>%
  summarise(Reads = sum(nreadsPedino), .groups = "drop") %>%
  pivot_wider(names_from = amplicon, values_from = Reads, values_fill = 0) %>%
  column_to_rownames("sample")

# Align Metadata
# Ensure metadata contains only the samples present in the matrix and in the same order
metadata_final <- size_filtered %>%
  select(-amplicon, -nreadsPedino, -relative_abundance, -abundance_percent) %>%
  distinct() %>%
  filter(sample %in% rownames(community_matrix)) %>%
  arrange(match(sample, rownames(community_matrix)))

# Run PERMANOVA: Does 'fraction_type' explain the variation in ASV composition?
permanova_result <- adonis2(community_matrix ~ size_fraction, 
                            data = metadata_final, 
                            permutations = 999, 
                            method = "bray")

print(permanova_result)

# Calculate multivariate dispersion
dist_matrix <- vegdist(community_matrix, method = "bray")
dispersion <- betadisper(dist_matrix, metadata_final$size_fraction)

# Test if dispersion is different between groups
permutest(dispersion)
