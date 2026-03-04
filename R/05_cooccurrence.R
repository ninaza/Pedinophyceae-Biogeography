# ============================================================
# Script 2: Co-occurrence analysis
# Input:  - algal abundance table (long format: sample, asv, count)
#         - algal metadata table  (one row per sample)
#         - dino_presence.csv     (output of Script 1)
# Output: - asv_cooccurrence_results.csv
#         - cooccurrence_partial_correlations.csv
#         - cooccurrence_rarefied_vs_unrarefied.pdf
#         - cooccurrence_confound_correction.png
#         - asv_cooccurrence_volcano.pdf
#         - asv_cooccurrence_rho_distribution.pdf
#         - asv_prevalence_vs_rho.pdf
#         - sig_asvs_dino_cooccurrence.pdf
#         - sig_asvs_size_fractions.pdf
# ============================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(vegan)

# ============================================================
# LOAD DATA
# ============================================================
abund_algae   <- read.csv("data/edited/pedinos_filtered_abundance.csv")  # sample | asv | count
metadata      <- read.csv("data/edited/eukbank_18S_V4_samples_edited.csv")   # sample | ocean_layer | habitat | size_fraction | ...
dino_presence <- read.csv("data/edited/dino_presence.csv")          # sample | dino_present | n_dino_asvs | ...

cat("Algal abundance table:", nrow(abund_algae), "rows\n")
cat("Metadata table:", nrow(metadata), "rows\n")
cat("Unique samples in abundance:", n_distinct(abund_algae$sample), "\n")
cat("Unique ASVs in abundance:", n_distinct(abund_algae$amplicon), "\n")

# ============================================================
# FILTER METADATA TO TARGET SAMPLES
# ============================================================
metadata_filtered <- metadata %>%
  filter(
    ocean_layer   == "sunlit",
    habitat       == "Marine",
    size_fraction %in% c("Pico", "Nano", "Micro")
  )

cat(sprintf("\nSamples after metadata filtering: %d\n", nrow(metadata_filtered)))
cat(sprintf("Size fraction breakdown:\n"))
print(table(metadata_filtered$size_fraction))

# ============================================================
# FILTER ABUNDANCE TABLE TO TARGET SAMPLES AND ASVs
# ============================================================
# Step 1: keep only samples that pass metadata filters
abund_filtered <- abund_algae %>%
  filter(sample %in% metadata_filtered$sample,
         nreadsPedino> 0)                      # presence = at least 1 read

# Step 2: keep only ASVs detected in at least one of the filtered samples
#         (this gives you your 567 ASVs)
target_asvs <- abund_filtered %>%
  distinct(amplicon) %>%
  pull(amplicon)

cat(sprintf("\nTarget algal ASVs after filtering: %d\n", length(target_asvs)))
# Expecting 597 — if different, check your metadata filter values

# ============================================================
# COLLAPSE TO SAMPLE-LEVEL PRESENCE/ABSENCE
# ============================================================
# Collapse across size fractions:
# a sample (= unique location × time) is positive if the algal group
# was detected in ANY size fraction at that station

# First: get unique station-level sample IDs
# If sample IDs already represent unique location×time (not per size fraction),
# use them directly. If sample IDs include size fraction, you need to
# extract the station ID — adjust the group_by accordingly.

algae_by_sample <- abund_filtered %>%
  filter(amplicon %in% target_asvs) %>%
  group_by(sample) %>%
  summarise(
    algae_present = 1L,
    n_algal_asvs  = n_distinct(amplicon),
    .groups = "drop"
  ) %>%
  # Add samples from metadata with no algae detected
  right_join(
    metadata_filtered %>% distinct(sample),
    by = "sample"
  ) %>%
  mutate(
    algae_present = replace_na(algae_present, 0L),
    n_algal_asvs  = replace_na(n_algal_asvs,  0L)
  )

cat(sprintf("\nSamples with algae present: %d\n", sum(algae_by_sample$algae_present)))
cat(sprintf("Samples with algae absent:  %d\n", sum(algae_by_sample$algae_present == 0)))

# ============================================================
# JOIN ALGAE AND DINOFLAGELLATE PRESENCE
# ============================================================
# Summary of dino ASV richness by algae presence

cooccur_df <- algae_by_sample %>%
  inner_join(dino_presence, by = "sample")

cooccur_df %>%
  group_by(algae_present) %>%
  summarise(
    n_samples       = n(),
    median_dino_asvs = median(n_dino_asvs),
    mean_dino_asvs   = mean(n_dino_asvs),
    sd_dino_asvs     = sd(n_dino_asvs),
    .groups = "drop"
  ) %>%
  mutate(algae_present = factor(algae_present,
                                labels = c("Algae absent",
                                           "Algae present"))) %>%
  print()

# ============================================================
# 1. WILCOXON TEST: DINOFLAGELLATE ASV RICHNESS
# ============================================================
# Non-parametric — appropriate for count data that is likely
# right-skewed and not normally distributed

wilcox_asvs <- wilcox.test(n_dino_asvs ~ algae_present,
                           data    = cooccur_df,
                           exact   = FALSE,
                           conf.int = TRUE)

cat(sprintf(
  "Wilcoxon test (n dino ASVs ~ algae presence):\n  W = %.1f\n  p = %.4f\n  95%% CI of difference: %.2f to %.2f\n",
  wilcox_asvs$statistic,
  wilcox_asvs$p.value,
  wilcox_asvs$conf.int[1],
  wilcox_asvs$conf.int[2]
))

# Effect size: rank-biserial correlation
# r = 0.1 small, 0.3 medium, 0.5 large
n_absent  <- sum(cooccur_df$algae_present == 0)
n_present <- sum(cooccur_df$algae_present == 1)
r_effect  <- 1 - (2 * wilcox_asvs$statistic) / (n_absent * n_present)
cat(sprintf("  Rank-biserial r = %.3f (0.1=small, 0.3=medium, 0.5=large)\n",
            r_effect))

# ============================================================
# CONTROL FOR TEMPERATURE (AND SEQUENCING DEPTH)
# Partial Spearman correlation:
# Is algae presence associated with dino richness BEYOND
# what is explained by shared temperature preference?
# ============================================================
# Replace n_dino_asvs with rarefied version in cooccur_df
dino_richness_rarefied <- read.csv("data/edited/dino_richness_rarefied.csv")

cooccur_df_rarefied <- cooccur_df %>%
  left_join(dino_richness_rarefied, by = "sample") %>%
  filter(!is.na(n_dino_asvs_rarefied))  # drop samples lost to rarefaction

# Re-run Wilcoxon with rarefied richness
wilcox_rarefied <- wilcox.test(n_dino_asvs_rarefied ~ algae_present,
                               data     = cooccur_df_rarefied,
                               exact    = FALSE,
                               conf.int = TRUE)

cat(sprintf(
  "\nWilcoxon test after rarefaction:\n  W = %.1f\n  p = %.4f\n",
  wilcox_rarefied$statistic,
  wilcox_rarefied$p.value
))

n_absent_r  <- sum(cooccur_df_rarefied$algae_present == 0)
n_present_r <- sum(cooccur_df_rarefied$algae_present == 1)
r_effect_r  <- 1 - (2 * wilcox_rarefied$statistic) / (n_absent_r * n_present_r)
cat(sprintf("  Rank-biserial r = %.3f (0.1=small, 0.3=medium, 0.5=large)\n",
            r_effect_r))

# ============================================================
# VISUALISE: RAREFIED VS UNRAREFIED SIDE BY SIDE
# ============================================================

# Unrarefied plot
p_unrarefied <- cooccur_df %>%
  mutate(algae_label = factor(algae_present,
                              labels = c("Algae absent", "Algae present"))) %>%
  ggplot(aes(x = algae_label, y = n_dino_asvs, fill = algae_label)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 1.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.2) +
  annotate("segment",
           x = 1, xend = 2,
           y    = max(cooccur_df$n_dino_asvs) * 1.05,
           yend = max(cooccur_df$n_dino_asvs) * 1.05) +
  annotate("text",
           x     = 1.5,
           y     = max(cooccur_df$n_dino_asvs) * 1.09,
           label = ifelse(wilcox_asvs$p.value < 0.001, "p < 0.001",
                          sprintf("p = %.3f", wilcox_asvs$p.value)),
           size = 3.5, fontface = "italic") +
  scale_fill_manual(values = c("Algae absent"  = "grey80",
                               "Algae present" = "#2196F3")) +
  labs(x        = NULL,
       y        = "Dinoflagellate ASV richness",
       title    = "Unrarefied",
       subtitle = sprintf("n = %d samples", nrow(cooccur_df))) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# Rarefied plot
p_rarefied <- cooccur_df_rarefied %>%
  mutate(algae_label = factor(algae_present,
                              labels = c("Algae absent", "Algae present"))) %>%
  ggplot(aes(x = algae_label, y = n_dino_asvs_rarefied, fill = algae_label)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 1.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.2) +
  annotate("segment",
           x = 1, xend = 2,
           y    = max(cooccur_df_rarefied$n_dino_asvs_rarefied) * 1.05,
           yend = max(cooccur_df_rarefied$n_dino_asvs_rarefied) * 1.05) +
  annotate("text",
           x     = 1.5,
           y     = max(cooccur_df_rarefied$n_dino_asvs_rarefied) * 1.09,
           label = ifelse(wilcox_rarefied$p.value < 0.001, "p < 0.001",
                          sprintf("p = %.3f", wilcox_rarefied$p.value)),
           size = 3.5, fontface = "italic") +
  scale_fill_manual(values = c("Algae absent"  = "grey80",
                               "Algae present" = "#2196F3")) +
  labs(x        = NULL,
       y        = NULL,
       title    = "Rarefied",
       subtitle = sprintf("n = %d samples retained", nrow(cooccur_df_rarefied))) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# Combine
final_plot <- p_unrarefied + p_rarefied +
  plot_annotation(
    title    = "Dinoflagellate ASV richness in algae-present vs. absent samples",
    subtitle = "Left: raw richness | Right: after rarefaction to equal sequencing depth",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 11, color = "grey40"))
  )

print(final_plot)
ggsave("plots/cooccurrence_rarefied_vs_unrarefied.pdf",
       final_plot, width = 10, height = 6)

# ============================================================
# ADD TEMPERATURE TO cooccur_df
# ============================================================
# Temperature should already be in your metadata table
cooccur_df_rarefied <- cooccur_df_rarefied %>%
  left_join(metadata %>% distinct(sample, temperature),
            by = "sample")

# Sanity check
cat(sprintf("Samples with temperature data: %d / %d\n",
            sum(!is.na(cooccur_df_rarefied$temperature)),
            nrow(cooccur_df_rarefied)))

# ============================================================
# CHECK THE CONFOUND VISUALLY
# ============================================================
# Is dino richness correlated with temperature?
# Is algae presence correlated with temperature?
# If both are yes, temperature is a confound

p_dino_temp <- ggplot(cooccur_df_rarefied,
                      aes(x = temperature, y = n_dino_asvs_rarefied)) +
  geom_point(alpha = 0.3, size = 1.5, color = "steelblue") +
  geom_smooth(method = "loess", color = "tomato") +
  labs(x = "Temperature (°C)", y = "Dinoflagellate ASV richness (rarefied)",
       title = "Dino richness ~ Temperature") +
  theme_classic(base_size = 12)

p_algae_temp <- ggplot(cooccur_df_rarefied,
                       aes(x = temperature,
                           y = algae_present,
                           color = factor(algae_present))) +
  geom_jitter(height = 0.02, alpha = 0.3, size = 1.5) +
  geom_smooth(aes(group = 1), method = "glm",
              method.args = list(family = binomial),
              color = "tomato", se = TRUE) +
  scale_color_manual(values = c("0" = "grey60", "1" = "#2196F3")) +
  labs(x = "Temperature (°C)", y = "Algae present (0/1)",
       title = "Algae presence ~ Temperature") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

print(p_dino_temp + p_algae_temp)

# Quantify both correlations
cor_dino_temp  <- cor.test(cooccur_df_rarefied$n_dino_asvs_rarefied,
                           cooccur_df_rarefied$temperature,
                           method = "spearman", exact = FALSE)
cor_algae_temp <- cor.test(cooccur_df_rarefied$algae_present,
                           cooccur_df_rarefied$temperature,
                           method = "spearman", exact = FALSE)

cat(sprintf("Dino richness ~ temperature:  rho = %.3f, p = %.4f\n",
            cor_dino_temp$estimate,  cor_dino_temp$p.value))
cat(sprintf("Algae presence ~ temperature: rho = %.3f, p = %.4f\n",
            cor_algae_temp$estimate, cor_algae_temp$p.value))

# ============================================================
# PARTIAL CORRELATION CONTROLLING FOR TEMPERATURE
# ============================================================
# control for temperature only
res_algae_temp <- residuals(lm(algae_present        ~ temperature,
                               data = cooccur_df_rarefied))
res_dino_temp  <- residuals(lm(n_dino_asvs_rarefied ~ temperature,
                               data = cooccur_df_rarefied))

partial_temp <- cor.test(res_algae_temp, res_dino_temp,
                         method = "spearman", exact = FALSE)

# ============================================================
# SUMMARY TABLE — UNCORRECTED VS TEMPERATURE CORRECTED
# ============================================================
uncorrected <- cor.test(cooccur_df_rarefied$algae_present,
                        cooccur_df_rarefied$n_dino_asvs_rarefied,
                        method = "spearman", exact = FALSE)

results_summary <- tibble(
  model   = c("Uncorrected",
              "Temperature corrected"),
  rho     = c(uncorrected$estimate,
              partial_temp$estimate),
  p_value = c(uncorrected$p.value,
              partial_temp$p.value)
) %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ "n.s."
  ))

print(results_summary)
write.csv(results_summary, "results/cooccurrence_partial_correlations.csv",
          row.names = FALSE)

# ============================================================
# VISUALISE THE KEY COMPARISON
# ============================================================
ggplot(results_summary,
       aes(x = model, y = rho, fill = sig)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = sprintf("rho=%.2f\n%s", rho, sig)),
            vjust = ifelse(results_summary$rho >= 0, -0.3, 1.3),
            size  = 3.5) +
  scale_fill_manual(values = c("***"  = "#D32F2F",
                               "**"   = "#FF5722",
                               "*"    = "#FF9800",
                               "n.s." = "grey70")) +
  scale_x_discrete(limits = results_summary$model) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(x     = NULL,
       y     = "Spearman rho",
       fill  = "Significance",
       title = "Co-occurrence signal after confound correction") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none")

ggsave("plots/cooccurrence_confound_correction.png", width = 8, height = 6)

# ============================================================
# 2. ASV-LEVEL: BUILD SAMPLE × ASV MATRIX
# ============================================================
# One row per sample, one column per ASV
# Collapsed across size fractions: present = detected in ANY fraction
algae_mat <- abund_filtered %>%
  filter(amplicon %in% target_asvs) %>%
  group_by(sample, amplicon) %>%
  summarise(presence = 1L, .groups = "drop") %>%
  # Fill in zeros for sample × ASV combinations with no detection
  complete(sample = metadata_filtered %>% distinct(sample) %>% pull(sample),
           amplicon,
           fill = list(presence = 0L)) %>%
  pivot_wider(names_from  = amplicon,
              values_from = presence,
              values_fill = 0L) %>%
  arrange(sample) %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Get samples that survived rarefaction
samples_rarefied <- dino_richness_rarefied %>%
  pull(sample)

cat(sprintf("Samples in algae_mat:             %d\n", nrow(algae_mat)))
cat(sprintf("Samples that survived rarefaction: %d\n", length(samples_rarefied)))
cat(sprintf("Samples dropped by rarefaction:    %d\n",
            nrow(algae_mat) - length(samples_rarefied)))

# Subset algae_mat to only rarefied samples
algae_mat_sub <- algae_mat[rownames(algae_mat) %in% samples_rarefied, ]

# Align dino_richness_vec to the same sample order as algae_mat_sub
dino_richness_vec <- dino_richness_rarefied %>%
  filter(sample %in% rownames(algae_mat_sub)) %>%
  arrange(sample) %>%
  pull(n_dino_asvs_rarefied)

# Ensure algae_mat_sub is sorted in the same order
algae_mat_sub <- algae_mat_sub[order(rownames(algae_mat_sub)), ]

# Verify alignment
cat(sprintf("\nSamples in algae_mat_sub:   %d\n", nrow(algae_mat_sub)))
cat(sprintf("Length of dino_richness_vec: %d\n", length(dino_richness_vec)))
cat(sprintf("Samples aligned correctly:   %s\n",
            all(rownames(algae_mat_sub) == dino_richness_rarefied %>%
                  filter(sample %in% rownames(algae_mat_sub)) %>%
                  arrange(sample) %>%
                  pull(sample))))

# ============================================================
# ASV-LEVEL: SPEARMAN CORRELATION WITH DINO PRESENCE
# ============================================================
cat("Running ASV-level co-occurrence...\n")

asv_results <- map_dfr(colnames(algae_mat_sub), function(amplicon) {
  asv_vec   <- algae_mat_sub[, amplicon]
  n_present <- sum(asv_vec)
  
  if (n_present < 3) {
    return(tibble(amplicon      = amplicon,
                  n_present   = n_present,
                  rho         = NA_real_,
                  p_value     = NA_real_,
                  skip_reason = "fewer than 3 occurrences"))
  }
  
  if (sd(asv_vec) == 0) {
    return(tibble(amplicon      = amplicon,
                  n_present   = n_present,
                  rho         = NA_real_,
                  p_value     = NA_real_,
                  skip_reason = "zero variance in ASV vector"))
  }
  
  ct <- suppressWarnings(
    cor.test(asv_vec, dino_richness_vec, method = "spearman", exact = FALSE)
  )
  tibble(amplicon      = amplicon,
         n_present   = n_present,
         rho         = ct$estimate,
         p_value     = ct$p.value,
         skip_reason = NA_character_)
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "n.s."
    )
  ) %>%
  arrange(desc(rho))

n_tested <- sum(!is.na(asv_results$rho))
n_sig    <- sum(asv_results$p_adj < 0.05 & asv_results$rho > 0, na.rm = TRUE)
n_skip   <- sum(is.na(asv_results$rho))

cat(sprintf("ASVs tested:   %d\n", n_tested))
cat(sprintf("ASVs skipped (< 3 occurrences): %d\n", n_skip))
cat(sprintf("Significantly co-occurring (FDR < 0.05, rho > 0): %d\n", n_sig))

write.csv(asv_results, "results/asv_cooccurrence_results.csv", row.names = FALSE)

# ============================================================
# ASV-LEVEL: VISUALISE
# ============================================================

# Volcano plot
# Add prevalence as proportion of samples
asv_results <- asv_results %>%
  mutate(prevalence = n_present / nrow(algae_mat_sub))

p_volcano <- ggplot(asv_results %>% filter(!is.na(rho)),
                    aes(x     = rho,
                        y     = -log10(p_adj),
                        color = prevalence,
                        size  = prevalence)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40") +
  scale_color_viridis_c(option    = "plasma",
                        name      = "Prevalence\n(prop. samples)",
                        labels    = scales::percent,
                        direction = -1) +
  scale_size_continuous(range  = c(1, 4),
                        name   = "Prevalence\n(prop. samples)",
                        labels = scales::percent) +
  guides(color = guide_legend(), size = guide_legend()) +
  labs(x        = "Spearman rho (co-occurrence with dinoflagellates)",
       y        = "-log10(FDR-adjusted p-value)",
       title    = "ASV-level co-occurrence with dinoflagellates",
       subtitle = sprintf("%d / %d ASVs significantly co-occur (FDR < 0.05, rho > 0)",
                          n_sig, n_tested)) +
  theme_classic(base_size = 13)
print(p_volcano)

ggsave("plots/asv_cooccurrence_volcano.pdf", p_volcano, width = 8, height = 6)

# Rho distribution
p_rho <- ggplot(asv_results %>% filter(!is.na(rho)),
                aes(x = rho, fill = sig)) +
  geom_histogram(bins = 40, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("***"  = "#D32F2F",
                               "**"   = "#FF5722",
                               "*"    = "#FF9800",
                               "n.s." = "grey70")) +
  labs(x     = "Spearman rho",
       y     = "Number of ASVs",
       fill  = "Significance",
       title = "Distribution of co-occurrence strength across ASVs") +
  theme_classic(base_size = 13)
print(p_rho)

ggsave("plots/asv_cooccurrence_rho_distribution.pdf", p_rho, width = 7, height = 5)

# ============================================================
# SUPPLEMENTARY: rho vs prevalence directly
# shows the parabola explicitly
# ============================================================
p_prev <- ggplot(asv_results %>% filter(!is.na(rho)),
                 aes(x = prevalence, y = rho, color = sig)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "loess", color = "black", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = c("***"  = "#D32F2F",
                                "**"   = "#FF5722",
                                "*"    = "#FF9800",
                                "n.s." = "grey70")) +
  scale_x_continuous(labels = scales::percent) +
  labs(x        = "ASV prevalence (proportion of samples)",
       y        = "Spearman rho",
       color    = "Significance",
       title    = "Co-occurrence strength vs. ASV prevalence") +
  theme_classic(base_size = 13)

print(p_prev)
ggsave("plots/asv_prevalence_vs_rho.pdf", p_prev, width = 8, height = 6)

# ============================================================
# EXTRACT SIGNIFICANT ASVs
# ============================================================
sig_asvs <- asv_results %>%
  filter(p_adj < 0.05, rho > 0, !is.na(rho)) %>%
  arrange(desc(rho))

cat(sprintf("Algal ASVs significantly co-occurring with dinos: %d / %d\n",
            nrow(sig_asvs), sum(!is.na(asv_results$rho))))

# ============================================================
# EXPLORE RHO DISTRIBUTION BEFORE CHOOSING A CUTOFF
# ============================================================
cat("Rho distribution for significant ASVs (p_adj < 0.05, rho > 0):\n")
print(summary(sig_asvs$rho))

# Visualise the distribution
hist(sig_asvs$rho,
     breaks = 30,
     col    = "#2196F3",
     border = "white",
     main   = "Rho distribution — significant algal ASVs",
     xlab   = "Spearman rho",
     ylab   = "Number of ASVs")
abline(v = c(0.2, 0.3, 0.4, 0.5),
       col = c("orange", "red", "darkred", "black"),
       lty = 2, lwd = 1.5)
legend("topright",
       legend = c("rho = 0.2", "rho = 0.3", "rho = 0.4", "rho = 0.5"),
       col    = c("orange", "red", "darkred", "black"),
       lty    = 2, lwd = 1.5, cex = 0.8)

# Count ASVs retained at each threshold
thresholds <- c(0.2, 0.3, 0.4, 0.5)
map_dfr(thresholds, function(t) {
  tibble(
    rho_threshold = t,
    n_asvs        = sum(sig_asvs$rho > t),
    pct_retained  = sprintf("%.0f%%",
                            100 * sum(sig_asvs$rho > t) / nrow(sig_asvs))
  )
}) %>% print()

# Apply chosen threshold — adjust rho_min as appropriate
rho_min <- 0.25

sig_asvs_filtered <- sig_asvs %>%
  filter(rho > rho_min)

cat(sprintf("ASVs retained (p_adj < 0.05, rho > %.1f): %d / %d\n",
            rho_min, nrow(sig_asvs_filtered), nrow(sig_asvs)))

sig_asvs_annotated <- sig_asvs_filtered %>%
  # Add sample-level metadata for each ASV occurrence
  left_join(
    # Summarise metadata across all samples where ASV is present
    abund_filtered %>%
      filter(amplicon %in% sig_asvs_filtered$amplicon) %>%
      left_join(metadata_filtered, by = "sample") %>%
      group_by(amplicon) %>%
      summarise(
        n_samples         = n_distinct(sample),
        size_fractions    = paste(sort(unique(size_fraction)),
                                  collapse = ", "),
        dominant_fraction = names(sort(table(size_fraction),
                                       decreasing = TRUE))[1],
        mean_temperature  = mean(temperature, na.rm = TRUE),
        .groups           = "drop"
      ),
    by = c("amplicon" = "amplicon")
  )

# ============================================================
# VISUALISE — ranked by rho
# ============================================================
p_sig <- sig_asvs_annotated %>%
  mutate(amplicon = fct_reorder(amplicon, rho)) %>%
  ggplot(aes(x = rho, y = amplicon,
             color = dominant_fraction,
             size  = n_present)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_size_continuous(name = "N samples\npresent") +
  labs(x        = "Spearman rho",
       y        = "Algal ASV",
       color    = "Dominant\nsize fraction",
       title    = "Algal ASVs significantly co-occurring with dinoflagellates",
       subtitle = sprintf("%d ASVs, FDR < 0.05, rho > %.1f",
                          nrow(sig_asvs_annotated), rho_min)) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 7))

print(p_sig)
ggsave("plots/sig_asvs_dino_cooccurrence.pdf", p_sig,
       width  = 10,
       height = 5)


# ============================================================
# SIZE FRACTION PRESENCE PER SIGNIFICANT ASV
# ============================================================

# Expand to one row per ASV × size fraction combination
asv_fraction <- abund_filtered %>%
  filter(amplicon %in% sig_asvs_annotated$amplicon) %>%
  left_join(metadata_filtered, by = "sample") %>%
  group_by(amplicon, size_fraction) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop")

# Join rho for ordering
asv_fraction <- asv_fraction %>%
  left_join(sig_asvs_annotated %>% select(amplicon, rho),
            by = "amplicon") %>%
  mutate(amplicon = fct_reorder(amplicon, rho))

# Tile plot — presence in each size fraction
p_fractions <- ggplot(asv_fraction,
                      aes(x    = size_fraction,
                          y    = amplicon,
                          fill = n_samples)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma",
                       name   = "N samples",
                       trans  = "log1p") +
  labs(x        = "Size fraction",
       y        = "Algal ASV (ordered by rho)",
       title    = "Size fraction occurrence of co-occurring algal ASVs",
       subtitle = sprintf("%d ASVs significantly co-occurring with dinos (FDR < 0.05, rho > %.1f)",
                          nrow(sig_asvs_annotated), rho_min)) +
  theme_classic(base_size = 12) +
  theme(axis.text.y  = element_text(size = 7),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.grid   = element_blank())
print(p_fractions)

ggsave("plots/sig_asvs_size_fractions.pdf",
       p_fractions, width  = 10, height = 5)
