# ============================================================
# Script 1: Prepare dinoflagellate presence/absence data
# Input:  - full abundance table (long format: sample, ASV, count)
#         - taxonomy table (one row per ASV: amplicon, taxogroup1, ...)
# Output: - dino_presence.csv (one row per sample, presence/absence)
# ============================================================

library(tidyverse)
library(vegan)

# ============================================================
# 1. LOAD DATA
# ============================================================
abund_full <- read.delim("data/raw/eukbank_18S_V4_counts.tsv", sep = "\t", header = T, dec = ".")   # sample | asv | count
taxonomy   <- read.csv("data/raw/eukbank_18S_V4_asvs.tsv", sep = "\t", header = T, dec = ".")          # amplicon | taxogroup1 | ...

# Quick checks
cat("Full abundance table dimensions:", nrow(abund_full), "rows\n")
cat("Taxonomy table dimensions:", nrow(taxonomy), "rows\n")
cat("Unique ASVs in abundance table:", n_distinct(abund_full$amplicon), "\n")
cat("Unique samples in abundance table:", n_distinct(abund_full$sample), "\n")

# ============================================================
# 2. FILTER TAXONOMY TO DINOFLAGELLATES
# ============================================================
dino_asvs <- taxonomy %>%
  filter(taxogroup1 == "Dinoflagellata") %>%
  pull(amplicon)

cat(sprintf("Dinoflagellate ASVs in taxonomy table: %d\n", length(dino_asvs)))

# ============================================================
# 3. SUBSET ABUNDANCE TABLE TO DINOFLAGELLATES
# ============================================================
abund_dino <- abund_full %>%
  filter(amplicon %in% dino_asvs,
         nreads > 0)          # presence = at least 1 read

cat(sprintf("Dinoflagellate abundance rows (count > 0): %d\n", nrow(abund_dino)))

# ============================================================
# 4. COLLAPSE TO SAMPLE-LEVEL PRESENCE/ABSENCE
# ============================================================
# One row per sample: was ANY dinoflagellate ASV detected?
dino_presence <- abund_dino %>%
  group_by(sample) %>%
  summarise(
    dino_present  = 1L,                  # if it appears here, at least one was present
    n_dino_asvs   = n_distinct(amplicon),     # how many unique dino ASVs detected
    total_dino_reads = sum(nreads),
    .groups = "drop"
  ) %>%
  # Add samples where NO dinoflagellate was detected (presence = 0)
  right_join(
    tibble(sample = unique(abund_full$sample)),
    by = "sample"
  ) %>%
  mutate(
    dino_present     = replace_na(dino_present, 0L),
    n_dino_asvs      = replace_na(n_dino_asvs, 0L),
    total_dino_reads = replace_na(total_dino_reads, 0L)
  )

cat(sprintf("\nSamples with dinoflagellates present:  %d\n",
            sum(dino_presence$dino_present)))
cat(sprintf("Samples with dinoflagellates absent:   %d\n",
            sum(dino_presence$dino_present == 0)))

write.csv(dino_presence, "data/edited/dino_presence.csv", row.names = FALSE)
# ============================================================
# OPTION 2: RAREFY TO EQUAL SEQUENCING DEPTH
# Then recount dinoflagellate ASV richness
# ============================================================

# Work in long format — only pivot one chunk at a time
samples_all <- unique(abund_full$sample)
chunk_size  <- 200   # adjust down if still hitting memory limit
chunks      <- split(samples_all, ceiling(seq_along(samples_all) / chunk_size))

cat(sprintf("Processing %d samples in %d chunks of ~%d\n",
            length(samples_all), length(chunks), chunk_size))

# Determine rarefaction depth from total reads per sample (no wide matrix needed)
depth_per_sample <- abund_full %>%
  group_by(sample) %>%
  summarise(total_reads = sum(nreads), .groups = "drop")

rare_depth <- floor(quantile(depth_per_sample$total_reads, 0.1))
cat(sprintf("Rarefaction depth (10th percentile): %d reads\n", rare_depth))

# Rarefy chunk by chunk
dino_richness_rarefied <- map_dfr(seq_along(chunks), function(i) {
  cat(sprintf("  Processing chunk %d / %d...\r", i, length(chunks)))
  
  chunk_samples <- chunks[[i]]
  
  # Build wide matrix for this chunk only
  chunk_wide <- abund_full %>%
    filter(sample %in% chunk_samples,
           amplicon    %in% dino_asvs) %>%
    pivot_wider(names_from  = amplicon,
                values_from = nreads,
                values_fill = 0L) %>%
    column_to_rownames("sample") %>%
    as.matrix()
  
  # Drop samples below rarefaction depth
  keep <- rowSums(chunk_wide) >= rare_depth
  if (sum(keep) == 0) return(tibble())
  
  chunk_wide <- chunk_wide[keep, , drop = FALSE]
  
  # Rarefy
  chunk_rarefied <- suppressWarnings(
    rrarefy(chunk_wide, sample = rare_depth)
  )
  
  tibble(
    sample               = rownames(chunk_rarefied),
    n_dino_asvs_rarefied = rowSums(chunk_rarefied > 0)
  )
})

cat(sprintf("\nSamples retained after rarefaction: %d / %d\n",
            nrow(dino_richness_rarefied), length(samples_all)))

write.csv(dino_richness_rarefied, "data/edited/dino_richness_rarefied.csv", row.names = FALSE)

