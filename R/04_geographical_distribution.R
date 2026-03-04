library(vegan)
library(tidyverse)
library(ggplot2)
library(patchwork)

## PREP DATA
abundance_table <- read.table("data/edited/pedinos_abundance.csv", 
                              sep = ",", header = T, dec = ".")
metadata_pedino <- read.table("data/edited/pedinos_samples_edited.csv", 
                           sep = ",", header = T, dec = ".")
metadata_all <- read.table("data/edited/eukbank_18S_V4_samples_edited.csv", 
                       sep = ",", header = T, dec = ".")

## filter only strict size fractions and check ASV overlap and distribution over size fractions
## only use sunlit ocean samples
target_fractions <- c("Pico", "Nano", "Micro")

metadata_filtered <- metadata_all %>%
  filter(habitat == "Marine") %>%
  filter(size_fraction %in% target_fractions) %>%
  filter(ocean != "") %>%
  filter(ocean_layer == "sunlit") %>%
  select(sample, latitude, longitude, ocean_layer, size_fraction, temperature) %>%
  mutate(abs_lat = abs(latitude))

size_filtered <- inner_join(metadata_pedino, abundance_table) %>%
  filter(habitat == "Marine") %>%
  filter(ocean_layer == "sunlit") %>%
  filter(size_fraction %in% target_fractions) %>%
  filter(ocean != "")

pedinos <- size_filtered %>%
  group_by(sample) %>%
  summarise(
    Total_Reads = sum(nreadsPedino),
    Richness = n_distinct(amplicon),
    .groups = "drop"
  )

full_analysis_df <- metadata_filtered %>%
  left_join(pedinos, by = "sample") %>%
  mutate(
    # If it didn't match, it means 0 reads and 0 richness
    Total_Reads = replace_na(Total_Reads, 0),
    Richness = replace_na(Richness, 0),
    # Create the binary 1/0 column for Presence/Absence
    Present = ifelse(Total_Reads > 0, 1, 0)
  )

# Create 15-degree bins for a clear visual
occupancy_plot_df <- full_analysis_df %>%
  mutate(lat_bin = cut(abs_lat, breaks = seq(0, 90, by = 10), include.lowest = TRUE)) %>%
  group_by(lat_bin, size_fraction) %>%
  summarise(
    Occupancy_Rate = mean(Present), 
    Sample_Count = n(),
    .groups = "drop"
  )

occupancy <- ggplot(occupancy_plot_df, aes(x = lat_bin, y = Occupancy_Rate, fill = size_fraction)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Degrees from Equator",
       y = "Proportion of Samples where Present")
occupancy

ggsave(plot = occupancy, 
       filename = "plots/pedino_lat_occupancy.png", 
       dpi = 300, height = 6, width = 9, bg = "white")

#--------------------------------------------------------------------------------------------------
#### PERMUTATION TEST FOR OCCUPANCY DIFFERENCES ACROSS LATITUDE, TEMPERATURE AND SIZE FRACTION #### 
#--------------------------------------------------------------------------------------------------

# Set reference level (largest/most balanced fraction as reference)
df_sub <- full_analysis_df %>%
  mutate(size_fraction = relevel(factor(size_fraction), ref = "Pico"))

# CHECK COLLINEARITY BETWEEN LATITUDE AND TEMPERATURE
cor_test <- cor.test(df_sub$temperature, df_sub$abs_lat, method = "spearman")
cat(sprintf("Spearman correlation (temp ~ abs_lat): rho = %.3f, p = %.4f\n",
            cor_test$estimate, cor_test$p.value))

# Visualise the relationship
ggplot(df_sub, aes(x = abs_lat, y = temperature)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "loess", color = "tomato") +
  labs(x = "Absolute Latitude (°)", y = "Temperature (°C)",
       title = "Collinearity check: Temperature vs. Latitude") +
  theme_classic(base_size = 13)

# PERMUTATION TEST FUNCTION
permutation_test_predictor <- function(data, formula, n_perm = 999, seed = 42) {
  set.seed(seed)
  
  fit_obs      <- glm(formula, data = data, family = binomial)
  obs_deviance <- fit_obs$null.deviance - fit_obs$deviance
  obs_coefs    <- coef(fit_obs)
  
  perm_deviances <- numeric(n_perm)
  perm_coefs     <- matrix(NA, nrow = n_perm, ncol = length(obs_coefs),
                           dimnames = list(NULL, names(obs_coefs)))
  
  for (i in seq_len(n_perm)) {
    # Permute presence within size fraction to preserve marginal prevalence
    data_perm <- data %>%
      group_by(size_fraction) %>%
      mutate(Present = sample(Present)) %>%
      ungroup()
    
    fit_perm <- tryCatch(
      suppressWarnings(glm(formula, data = data_perm, family = binomial)),
      error = function(e) NULL
    )
    
    if (!is.null(fit_perm)) {
      perm_deviances[i] <- fit_perm$null.deviance - fit_perm$deviance
      # Guard: in case a permuted model drops a coefficient
      shared <- intersect(names(obs_coefs), names(coef(fit_perm)))
      perm_coefs[i, shared] <- coef(fit_perm)[shared]
    }
  }
  
  p_overall <- mean(perm_deviances >= obs_deviance, na.rm = TRUE)
  p_terms   <- sapply(names(obs_coefs), function(term) {
    mean(abs(perm_coefs[, term]) >= abs(obs_coefs[[term]]), na.rm = TRUE)
  })
  
  list(
    formula           = deparse(formula),
    observed_deviance = obs_deviance,
    observed_coefs    = obs_coefs,
    perm_deviances    = perm_deviances,
    perm_coefs        = perm_coefs,
    p_overall         = p_overall,
    p_terms           = p_terms,
    fit_observed      = fit_obs,
    McFadden_R2       = 1 - (fit_obs$deviance / fit_obs$null.deviance)
  )
}

# RUN PERMUTATION TEST AND COMPARE MODELS
# Model 1: latitude only (your existing result)
perm_lat  <- permutation_test_predictor(
  df_sub,
  Present ~ abs_lat * size_fraction
)

# Model 2: temperature only
perm_temp <- permutation_test_predictor(
  df_sub,
  Present ~ temperature * size_fraction
)

# McFadden's pseudo-R² tells you how much variance each model explains
# A higher R² = better explanatory power
# If temp R² >> lat R², temperature is the more proximal driver

model_comparison_reduced <- tibble(
  model       = c("Latitude only", "Temperature only"),
  McFadden_R2 = c(perm_lat$McFadden_R2,
                  perm_temp$McFadden_R2),
  deviance_explained = c(perm_lat$observed_deviance,
                         perm_temp$observed_deviance),
  p_overall   = c(perm_lat$p_overall,
                  perm_temp$p_overall)
) %>%
  mutate(
    # How much of the latitude model's deviance does temperature capture?
    pct_of_lat = deviance_explained / max(deviance_explained) * 100
  )

print(model_comparison_reduced)

# VARIANCE PARTITIONING
# ============================================================
# Even though you can't put both in one model, you can ask:
# "How much of the latitude effect disappears when you use temperature?"
# by comparing the latitude coefficient before and after accounting for temp
# via a partial correlation approach

# Partial correlation: effect of latitude on presence after removing temp signal
# Step 1: regress presence on temperature -> get residuals
# Step 2: regress latitude on temperature -> get residuals  
# Step 3: correlate both sets of residuals
# A near-zero partial correlation means temp fully mediates the lat effect

df_pico <- df_sub %>% filter(size_fraction == "Pico")

res_presence <- residuals(glm(Present ~ temperature,
                              data = df_pico, family = binomial))
res_lat      <- residuals(lm(abs_lat ~ temperature, data = df_pico))

partial_cor  <- cor.test(res_lat, res_presence, method = "spearman")
cat(sprintf(
  "Partial correlation (lat ~ presence | temperature): rho = %.3f, p = %.4f\n",
  partial_cor$estimate, partial_cor$p.value
))

# 3. VISUALISE MODEL COMPARISON
# ============================================================
# Side-by-side: predicted probability vs latitude AND vs temperature
# per size fraction — lets the reader visually judge which fits better

fracs  <- levels(df_sub$size_fraction)
colors <- c("Nano" = "#2196F3", "Micro" = "#FF5722", "Pico" = "#4CAF50")

plot_comparison <- lapply(fracs, function(frac) {
  df_frac <- df_sub %>% filter(size_fraction == frac)
  
  fit_lat  <- glm(Present ~ abs_lat,     data = df_frac, family = binomial)
  fit_temp <- glm(Present ~ temperature, data = df_frac, family = binomial)
  
  # Latitude plot
  lat_seq  <- data.frame(abs_lat = seq(0, 90, length.out = 200))
  pred_lat <- predict(fit_lat, newdata = lat_seq,
                      type = "response", se.fit = TRUE)
  lat_seq  <- lat_seq %>%
    mutate(prob  = pred_lat$fit,
           lower = pred_lat$fit - 1.96 * pred_lat$se.fit,
           upper = pred_lat$fit + 1.96 * pred_lat$se.fit)
  
  p1 <- ggplot() +
    geom_jitter(data = df_frac,
                aes(x = abs_lat, y = Present),
                height = 0.02, alpha = 0.3,
                color  = colors[[frac]]) +
    geom_ribbon(data = lat_seq,
                aes(x = abs_lat, ymin = lower, ymax = upper),
                fill = colors[[frac]], alpha = 0.2) +
    geom_line(data = lat_seq,
              aes(x = abs_lat, y = prob),
              color = colors[[frac]], linewidth = 1) +
    annotate("text", x = 5, y = 0.95, hjust = 0, size = 3.5,
             label = sprintf("R² = %.3f", 1 - fit_lat$deviance /
                               fit_lat$null.deviance)) +
    labs(title = if (frac == fracs[1]) "Latitude" else NULL,
         x = "Absolute Latitude (°)", y = paste(frac, "— P(Present)")) +
    scale_y_continuous(limits = c(-0.05, 1.05)) +
    theme_classic(base_size = 12)
  
  # Temperature plot
  temp_seq  <- data.frame(temperature = seq(min(df_sub$temperature, na.rm = TRUE),
                                            max(df_sub$temperature, na.rm = TRUE),
                                            length.out = 200))
  pred_temp <- predict(fit_temp, newdata = temp_seq,
                       type = "response", se.fit = TRUE)
  temp_seq  <- temp_seq %>%
    mutate(prob  = pred_temp$fit,
           lower = pred_temp$fit - 1.96 * pred_temp$se.fit,
           upper = pred_temp$fit + 1.96 * pred_temp$se.fit)
  
  p2 <- ggplot() +
    geom_jitter(data = df_frac,
                aes(x = temperature, y = Present),
                height = 0.02, alpha = 0.3,
                color  = colors[[frac]]) +
    geom_ribbon(data = temp_seq,
                aes(x = temperature, ymin = lower, ymax = upper),
                fill = colors[[frac]], alpha = 0.2) +
    geom_line(data = temp_seq,
              aes(x = temperature, y = prob),
              color = colors[[frac]], linewidth = 1) +
    annotate("text", x = -Inf, y = 0.95, hjust = -0.1, size = 3.5,
             label = sprintf("R² = %.3f", 1 - fit_temp$deviance /
                               fit_temp$null.deviance)) +
    labs(title = if (frac == fracs[1]) "Temperature" else NULL,
         x = "Temperature (°C)", y = NULL) +
    scale_y_continuous(limits = c(-0.05, 1.05)) +
    theme_classic(base_size = 12)
  
  list(p1, p2)
})

# Arrange in grid: rows = size fractions, columns = lat | temp
all_plots <- unlist(plot_comparison, recursive = FALSE)
model_comp_plot <- wrap_plots(all_plots, ncol = 2, byrow = TRUE)
print(model_comp_plot)
ggsave("plots/lat_vs_temp_comparison.png"
       , model_comp_plot, width = 10, height = 12)

# Summarise and visualise R² comparison per size fraction
# Build a tidy comparison table from your existing model objects
r2_comparison <- tibble(
  size_fraction = rep(levels(df_sub$size_fraction), each = 2),
  predictor     = rep(c("Latitude", "Temperature"),
                      times = length(levels(df_sub$size_fraction))),
  McFadden_R2   = map2_dbl(
    rep(levels(df_sub$size_fraction), each = 2),
    rep(c("abs_lat", "temperature"), times = length(levels(df_sub$size_fraction))),
    function(frac, pred) {
      df_frac <- df_sub %>% filter(size_fraction == frac)
      fit     <- glm(reformulate(pred, "Present"),
                     data = df_frac, family = binomial)
      1 - fit$deviance / fit$null.deviance
    }
  )
)

print(r2_comparison)

# Visualise as a grouped bar chart
r2_plot <- ggplot(r2_comparison,
       aes(x = size_fraction, y = McFadden_R2,
           fill = predictor)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", McFadden_R2)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Latitude"    = "#2196F3",
                               "Temperature" = "#FF5722")) +
  labs(x    = "Size fraction",
       y    = "McFadden's pseudo-R²",
       fill = "Predictor",
       title    = "Latitude vs. Temperature as predictors of presence/absence") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")
r2_plot
ggsave("plots/lat_vs_temp_R2.png"
       , r2_plot, width = 12, height = 10)

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

## Niche breadth: Calculate the range of temperatures where Pedinos is present for each size fraction
# 1. Create Temperature Bins
niche_data <- full_analysis_df %>%
  filter(Present == 1) %>%
  mutate(temp_bin = cut(temperature, breaks = seq(0, 35, by = 5)))

# 2. Calculate proportions (p_i) per fraction
# Abundance-Weighted Levins Breadth
weighted_breadth <- full_analysis_df %>%
  filter(Present == 1) %>%
  mutate(temp_bin = cut(temperature, breaks = seq(0, 35, by = 5))) %>%
  group_by(size_fraction, temp_bin) %>%
  # Use Total_Reads instead of n()
  summarise(sum_reads = sum(Total_Reads), .groups = "drop") %>%
  group_by(size_fraction) %>%
  mutate(p_i = sum_reads / sum(sum_reads)) %>%
  summarise(
    B_weighted = 1 / sum(p_i^2),
    n_bins = n_distinct(temp_bin),
    BA_weighted = (B_weighted - 1) / (n_bins - 1)
  )
print(weighted_breadth)

temp_niche <- ggplot(subset(full_analysis_df, Present == 1), aes(x = temperature, fill = size_fraction)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Thermal Niche Occupancy Profiles",
       subtitle = "Narrower peaks = lower niche breadth (Specialists)",
       x = "Temperature (°C)", y = "Density of Detections")
temp_niche

ggsave("plots/temp_niche_breadth.png"
       , temp_niche, width = 12, height = 10)
