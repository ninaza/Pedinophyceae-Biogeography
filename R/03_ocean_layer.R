library(tidyverse)
library(ggplot2)
library(dplyr)
library(vegan)
library(MetBrewer)

palette <- met.brewer("Renoir")

abundance_table <- read.table("data/edited/pedinos_filtered_abundance.csv", 
                              sep = ",", header = T, dec = ".")
metadata <- read.table("data/edited/pedinos_edit_filtered_metadata.csv", 
                       sep = ",", header = T, dec = ".")

join_marine <- inner_join(metadata, abundance_table) %>%
  filter(habitat == "Marine") %>%
  filter(!is.na(ocean_layer))

## filter only strict size fractions and check ASV overlap and distribution over size fractions
target_fractions <- c("Pico", "Nano", "Micro")

size_filtered <- join_marine %>%
  filter(size_fraction %in% target_fractions)

# Summarize the total abundance of your taxon per sample
depth_data <- size_filtered %>%
  group_by(sample) %>%
  summarise(Total_Taxon_Abundance = sum(nreadsPedino)) %>%
  left_join(join_marine, by = "sample")

ggplot(depth_data, aes(x = Total_Taxon_Abundance, y = depth)) +
  geom_point(show.legend=FALSE, alpha = 0.6) +
  scale_y_reverse() + # Put 0 at the top
  geom_hline(yintercept = 200, linetype = "dashed", color = "#5480B5FF") +
  labs(x = "Relative Abundance", y = "Depth (m)") +
  theme_minimal()

## abundance by ocean layer
ocena_layer <- ggplot(depth_data, aes(x = ocean_layer, y = Total_Taxon_Abundance, fill = ocean_layer)) +
  geom_violin(alpha = 0.5, scale = "width", trim = FALSE, color = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, 
               color = "black", linewidth = 0.5) +
  facet_wrap(~size_fraction) +
  scale_y_log10() +
  scale_fill_manual(values = c("sunlit" = "#FFD700", "dark" = "#2F4F4F")) +
  theme_minimal()
ocena_layer

ggsave(plot = ocena_layer, 
       filename = "plots/pedino_nreads_oceanlayer_violin.png", 
       dpi = 300, height = 6, width = 9, bg = "white")

# Filter for only your green algae taxon
green_algae_upset <- size_filtered %>%
  filter(nreadsPedino > 0) %>%
  group_by(amplicon, ocean_layer) %>%
  summarise(Present = 1, .groups = "drop") %>%
  pivot_wider(names_from = ocean_layer, values_from = Present, values_fill = 0) %>%
  as.data.frame()

png("plots/pedinos_oceanlayer_upset.png", width = 2400, height = 1800, res = 300)
upset(green_algae_upset, 
      sets = c("sunlit", "dark"), 
      main.bar.color = "steelblue", 
      sets.bar.color = "darkgray",
      order.by = "freq", 
      empty.intersections = "on")
dev.off()

# Statistical test for Sunlit vs Dark abundance (Pico fraction)
pico_data <- depth_data %>% filter(size_fraction == "Pico")
wilcox.test(Total_Taxon_Abundance ~ ocean_layer, data = pico_data)

nano_data <- depth_data %>% filter(size_fraction == "Nano")
wilcox.test(Total_Taxon_Abundance ~ ocean_layer, data = nano_data)
