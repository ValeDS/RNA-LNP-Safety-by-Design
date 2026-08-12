# 07_single_cell_human_validation_GSE171964.R
# Human retrospective validation using GSE171964

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

set.seed(123)

source(file.path("scripts", "config.R"))

# 1. Load synthetic signatures

synthetic_signatures <- read.csv(
  file.path(input_dir, "03_synthetic_pathway_signature_scores.csv")
)

expected_signatures <- list(
  interferon_response = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "RSAD2"),
  inflammatory_response = c("IL6", "IL1B", "TNF", "CXCL8", "CXCL10", "CCL2", "NFKBIA", "PTGS2", "SOCS3", "STAT1"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "CIITA", "HLA-DRA"),
  adaptive_immunity = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "MS4A1", "CD19", "IGHG1", "IGHM", "PRDM1"),
  stress_response = c("HSPA1A", "HSPA1B", "DDIT3", "ATF4", "XBP1", "DNAJB1", "HMOX1", "FOS", "JUN", "EGR1")
)

predicted_summary <- synthetic_signatures %>%
  select(-sample_id) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "signature",
    values_to = "predicted_activation_score"
  )

# 2. Load GSE171964 files

gse_dir <- normalizePath(
  Sys.getenv("GSE171964_DATA_DIR", unset = file.path(project_dir, "data", "GSE171964")),
  mustWork = FALSE
)

matrix_file <- file.path(gse_dir, "GSE171964_countsmatrix_v2.mtx.gz")
features_file <- file.path(gse_dir, "GSE171964_feats_v2.tsv.gz")
barcodes_file <- file.path(gse_dir, "GSE171964_barcodes_v2.tsv.gz")
pheno_file <- file.path(gse_dir, "GSE171964_geo_pheno_v2.csv.gz")

counts <- readMM(matrix_file)

features <- read.delim(features_file, header = FALSE)
barcodes <- read.delim(barcodes_file, header = FALSE)
pheno <- read.csv(pheno_file)

# Remove header-like first row if present
features <- features[-1, , drop = FALSE]
barcodes <- barcodes[-1, , drop = FALSE]

gene_names <- sub("^[0-9]+\\s+", "", features$V1)
cell_barcodes <- sub("^[0-9]+\\s+", "", barcodes$V1)

rownames(counts) <- gene_names
colnames(counts) <- cell_barcodes

# Align metadata to matrix columns
pheno <- pheno %>%
  filter(barcode %in% colnames(counts))

counts <- counts[, pheno$barcode]

stopifnot(all(colnames(counts) == pheno$barcode))

# 3. Define validation time points
# day 0 = baseline.
# We compare day 1 and day 2 as early innate response time points.
# Day 7 is retained as adaptive/late response.

pheno <- pheno %>%
  mutate(
    day = as.numeric(day),
    time_group = case_when(
      day == 0 ~ "Day0_baseline",
      day %in% c(1, 2) ~ "Early_post_vaccination",
      day == 7 ~ "Day7_post_vaccination",
      TRUE ~ paste0("Day", day)
    )
  )

# 4. Pseudo-bulk aggregation per patient × day

sample_groups <- pheno %>%
  mutate(pseudobulk_id = paste(pt_id, day, sep = "_Day")) %>%
  select(barcode, pseudobulk_id, pt_id, day, time_group, clustnm, sample_id)

unique_groups <- unique(sample_groups$pseudobulk_id)

message("Creating pseudo-bulk matrix...")

pseudobulk_mat <- sapply(unique_groups, function(g) {
  cells <- sample_groups$barcode[sample_groups$pseudobulk_id == g]
  Matrix::rowSums(counts[, cells, drop = FALSE])
})

pseudobulk_mat <- as.matrix(pseudobulk_mat)

pseudobulk_meta <- sample_groups %>%
  distinct(pseudobulk_id, pt_id, day, time_group)

rownames(pseudobulk_meta) <- pseudobulk_meta$pseudobulk_id
pseudobulk_meta <- pseudobulk_meta[colnames(pseudobulk_mat), ]

# CPM-like normalization + log1p
lib_size <- colSums(pseudobulk_mat)
pseudobulk_log <- log1p(t(t(pseudobulk_mat) / lib_size * 1e6))

# 5. Compute real signature scores

compute_signature_scores <- function(expr_mat, signatures) {
  out <- lapply(names(signatures), function(sig) {
    genes <- intersect(signatures[[sig]], rownames(expr_mat))

    if (length(genes) == 0) {
      score <- rep(NA_real_, ncol(expr_mat))
    } else {
      score <- colMeans(expr_mat[genes, , drop = FALSE], na.rm = TRUE)
    }

    data.frame(
      pseudobulk_id = colnames(expr_mat),
      signature = sig,
      real_signature_score = score,
      n_genes_found = length(genes)
    )
  }) %>% bind_rows()

  return(out)
}

real_scores_long <- compute_signature_scores(
  pseudobulk_log,
  expected_signatures
)

real_scores_long <- real_scores_long %>%
  left_join(pseudobulk_meta, by = "pseudobulk_id")

write.csv(
  real_scores_long,
  file.path(core_output_dir, "07_GSE171964_real_signature_scores_long.csv"),
  row.names = FALSE
)

# 6. Observed activation: post-vaccination minus baseline

baseline_scores <- real_scores_long %>%
  filter(time_group == "Day0_baseline") %>%
  group_by(signature) %>%
  summarise(
    baseline_mean = mean(real_signature_score, na.rm = TRUE),
    .groups = "drop"
  )

observed_activation <- real_scores_long %>%
  filter(time_group %in% c("Early_post_vaccination", "Day7_post_vaccination")) %>%
  group_by(signature, time_group) %>%
  summarise(
    post_mean = mean(real_signature_score, na.rm = TRUE),
    n_genes_found = max(n_genes_found, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(baseline_scores, by = "signature") %>%
  mutate(
    observed_activation_score = post_mean - baseline_mean
  )

write.csv(
  observed_activation,
  file.path(core_output_dir, "07_GSE171964_observed_signature_activation.csv"),
  row.names = FALSE
)

# 7. Predicted vs observed concordance

validation_comparison <- observed_activation %>%
  left_join(predicted_summary, by = "signature")

write.csv(
  validation_comparison,
  file.path(core_output_dir, "07_GSE171964_predicted_vs_observed_signature_comparison.csv"),
  row.names = FALSE
)

concordance_summary <- validation_comparison %>%
  group_by(time_group) %>%
  summarise(
    spearman_concordance = cor(
      predicted_activation_score,
      observed_activation_score,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

write.csv(
  concordance_summary,
  file.path(core_output_dir, "07_GSE171964_biological_concordance_score.csv"),
  row.names = FALSE
)

# 8. Cell composition dynamics

cell_composition <- pheno %>%
  count(pt_id, day, time_group, clustnm, name = "n_cells") %>%
  group_by(pt_id, day, time_group) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  cell_composition,
  file.path(core_output_dir, "07_GSE171964_cell_composition_by_patient_day.csv"),
  row.names = FALSE
)

cell_composition_summary <- cell_composition %>%
  group_by(day, time_group, clustnm) %>%
  summarise(
    mean_fraction = mean(fraction, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  cell_composition_summary,
  file.path(core_output_dir, "07_GSE171964_cell_composition_summary.csv"),
  row.names = FALSE
)

# 9. Figures

p_sig <- ggplot(
  validation_comparison,
  aes(
    x = predicted_activation_score,
    y = observed_activation_score,
    label = signature
  )
) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3.3) +
  facet_wrap(~ time_group) +
  labs(
    title = "GSE171964 human scRNA-seq validation",
    x = "Predicted synthetic signature activation",
    y = "Observed signature activation vs baseline"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "07_GSE171964_predicted_vs_observed_signatures.png"),
  p_sig,
  width = 10,
  height = 5,
  dpi = 300
)

p_time <- ggplot(
  real_scores_long,
  aes(x = factor(day), y = real_signature_score)
) +
  geom_boxplot() +
  facet_wrap(~ signature, scales = "free_y") +
  labs(
    title = "Temporal immune signature dynamics in GSE171964",
    x = "Day",
    y = "Pseudo-bulk signature score"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(core_output_dir, "07_GSE171964_signature_dynamics_by_day.png"),
  p_time,
  width = 11,
  height = 7,
  dpi = 300
)

top_cells <- cell_composition_summary %>%
  group_by(clustnm) %>%
  summarise(mean_fraction_all = mean(mean_fraction), .groups = "drop") %>%
  arrange(desc(mean_fraction_all)) %>%
  slice_head(n = 12) %>%
  pull(clustnm)

p_cells <- cell_composition_summary %>%
  filter(clustnm %in% top_cells) %>%
  ggplot(aes(x = factor(day), y = mean_fraction, group = clustnm)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ clustnm, scales = "free_y") +
  labs(
    title = "Major immune-cell composition dynamics in GSE171964",
    x = "Day",
    y = "Mean cell fraction"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(core_output_dir, "07_GSE171964_cell_composition_dynamics.png"),
  p_cells,
  width = 12,
  height = 8,
  dpi = 300
)

message("Human single-cell validation completed using GSE171964.")
