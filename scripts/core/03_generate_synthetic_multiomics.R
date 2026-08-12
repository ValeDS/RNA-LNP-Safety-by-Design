# 03_generate_synthetic_multiomics.R
# Safety-by-design RNA-LNP workflow
# Module 3: Synthetic multi-omics generation

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

set.seed(123)

source(file.path("scripts", "config.R"))
dir.create(core_output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Module 2 output

responses <- read.csv(file.path(core_output_dir, "02_patient_lnp_synthetic_immune_responses.csv"))

# 2. Select representative subset for multi-omics generation
# We avoid generating high-dimensional omics for all 10,000 pairs.
# This keeps the workflow lightweight and reproducible.

n_omics_samples <- 500

omics_samples <- responses %>%
  slice_sample(n = n_omics_samples) %>%
  mutate(sample_id = paste0("SAMPLE_", sprintf("%04d", 1:n())))

# 3. Define transcriptomic gene signatures

gene_sets <- list(
  interferon_response = c(
    "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "MX2",
    "OAS1", "OAS2", "OAS3", "RSAD2"
  ),
  inflammatory_response = c(
    "IL6", "IL1B", "TNF", "CXCL8", "CXCL10",
    "CCL2", "NFKBIA", "PTGS2", "SOCS3", "STAT1"
  ),
  antigen_presentation = c(
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1",
    "TAP2", "PSMB8", "PSMB9", "CIITA", "HLA-DRA"
  ),
  adaptive_immunity = c(
    "CD3D", "CD3E", "CD4", "CD8A", "CD8B",
    "MS4A1", "CD19", "IGHG1", "IGHM", "PRDM1"
  ),
  stress_response = c(
    "HSPA1A", "HSPA1B", "DDIT3", "ATF4", "XBP1",
    "DNAJB1", "HMOX1", "FOS", "JUN", "EGR1"
  )
)

all_genes <- unique(unlist(gene_sets))

# 4. Generate synthetic transcriptomic matrix

generate_transcriptomics <- function(samples, genes, seed = 123) {
  set.seed(seed)

  expr <- matrix(
    rnorm(nrow(samples) * length(genes), mean = 6, sd = 0.4),
    nrow = length(genes),
    ncol = nrow(samples)
  )

  rownames(expr) <- genes
  colnames(expr) <- samples$sample_id

  for (i in seq_len(nrow(samples))) {
    s <- samples[i, ]

    # Interferon genes track innate activation and cytokine burden
    expr[gene_sets$interferon_response, i] <-
      expr[gene_sets$interferon_response, i] +
      0.9 * s$innate_activation +
      0.5 * s$cytokine_burden

    # Inflammatory genes track cytokine burden and baseline inflammation
    expr[gene_sets$inflammatory_response, i] <-
      expr[gene_sets$inflammatory_response, i] +
      1.0 * s$cytokine_burden +
      0.7 * s$baseline_inflammation

    # Antigen presentation tracks delivery and adaptive activation
    expr[gene_sets$antigen_presentation, i] <-
      expr[gene_sets$antigen_presentation, i] +
      0.8 * s$delivery_efficiency +
      0.5 * s$adaptive_activation

    # Adaptive immunity tracks adaptive activation
    expr[gene_sets$adaptive_immunity, i] <-
      expr[gene_sets$adaptive_immunity, i] +
      0.9 * s$adaptive_activation

    # Stress response tracks off-target activation
    expr[gene_sets$stress_response, i] <-
      expr[gene_sets$stress_response, i] +
      0.7 * s$off_target_activation
  }

  expr <- pmax(expr, 0)
  return(expr)
}

transcriptomics_matrix <- generate_transcriptomics(
  omics_samples,
  all_genes,
  seed = 123
)

# 5. Generate synthetic cytokine layer

generate_cytokine_layer <- function(samples, seed = 124) {
  set.seed(seed)

  cytokines <- samples %>%
    transmute(
      sample_id,
      IL6   = pmax(0, 5 + 4.0 * cytokine_burden + 2.0 * baseline_inflammation + rnorm(n(), 0, 1)),
      TNFa  = pmax(0, 4 + 3.5 * innate_activation + 1.5 * baseline_inflammation + rnorm(n(), 0, 1)),
      IL1b  = pmax(0, 3 + 3.0 * innate_activation + 2.0 * cytokine_burden + rnorm(n(), 0, 1)),
      IFNg  = pmax(0, 4 + 3.0 * adaptive_activation + 1.0 * innate_activation + rnorm(n(), 0, 1)),
      IFNa  = pmax(0, 3 + 3.5 * innate_activation + 1.5 * cytokine_burden + rnorm(n(), 0, 1)),
      CXCL10 = pmax(0, 4 + 4.0 * innate_activation + 2.0 * cytokine_burden + rnorm(n(), 0, 1))
    )

  return(cytokines)
}

cytokine_layer <- generate_cytokine_layer(omics_samples, seed = 124)

# 6. Generate immune-cell dynamics layer

generate_cell_layer <- function(samples, seed = 125) {
  set.seed(seed)

  cells <- samples %>%
    transmute(
      sample_id,
      CD4_T_cells = pmax(0, 100 + 30 * adaptive_activation + rnorm(n(), 0, 8)),
      CD8_T_cells = pmax(0, 90 + 35 * adaptive_activation + 10 * delivery_efficiency + rnorm(n(), 0, 8)),
      B_cells = pmax(0, 80 + 25 * adaptive_activation + rnorm(n(), 0, 8)),
      Plasma_cells = pmax(0, 20 + 18 * adaptive_activation + rnorm(n(), 0, 5)),
      Monocytes = pmax(0, 120 + 25 * innate_activation + 20 * baseline_inflammation + rnorm(n(), 0, 10)),
      NK_cells = pmax(0, 70 + 20 * innate_activation + rnorm(n(), 0, 7))
    )

  return(cells)
}

cell_layer <- generate_cell_layer(omics_samples, seed = 125)

# 7. Compute pathway-level synthetic signatures

compute_signature_scores <- function(expr, gene_sets) {
  scores <- sapply(gene_sets, function(gs) {
    colMeans(expr[gs, , drop = FALSE])
  })

  scores <- as.data.frame(scores)
  scores$sample_id <- rownames(scores)
  scores <- scores %>% relocate(sample_id)

  return(scores)
}

signature_scores <- compute_signature_scores(
  transcriptomics_matrix,
  gene_sets
)

# 8. Export outputs

write.csv(
  omics_samples,
  file.path(core_output_dir, "03_multiomics_sample_metadata.csv"),
  row.names = FALSE
)

write.csv(
  transcriptomics_matrix,
  file.path(core_output_dir, "03_synthetic_transcriptomics_matrix.csv")
)

write.csv(
  cytokine_layer,
  file.path(core_output_dir, "03_synthetic_cytokine_layer.csv"),
  row.names = FALSE
)

write.csv(
  cell_layer,
  file.path(core_output_dir, "03_synthetic_immune_cell_layer.csv"),
  row.names = FALSE
)

write.csv(
  signature_scores,
  file.path(core_output_dir, "03_synthetic_pathway_signature_scores.csv"),
  row.names = FALSE
)

# 9. Visualizations

# Heatmap of pathway signatures
sig_mat <- signature_scores
rownames(sig_mat) <- sig_mat$sample_id
sig_mat$sample_id <- NULL

pheatmap(
  t(scale(sig_mat)),
  show_colnames = FALSE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  filename = file.path(core_output_dir, "03_pathway_signature_heatmap.png"),
  width = 9,
  height = 5
)

# Cytokine burden vs IFNa
plot_df <- omics_samples %>%
  select(sample_id, cytokine_burden, innate_activation) %>%
  left_join(cytokine_layer, by = "sample_id")

p_ifna <- ggplot(plot_df, aes(x = cytokine_burden, y = IFNa)) +
  geom_point(alpha = 0.5) +
  labs(
    title = "Synthetic cytokine layer: IFNa vs cytokine burden",
    x = "Cytokine burden",
    y = "IFNa"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "03_ifna_vs_cytokine_burden.png"),
  p_ifna,
  width = 8,
  height = 6,
  dpi = 300
)

# Immune-cell dynamics
cell_long <- cell_layer %>%
  tidyr::pivot_longer(
    cols = -sample_id,
    names_to = "cell_type",
    values_to = "abundance"
  )

p_cells <- ggplot(cell_long, aes(x = cell_type, y = abundance)) +
  geom_boxplot() +
  coord_flip() +
  labs(
    title = "Synthetic immune-cell dynamics",
    x = "Immune cell type",
    y = "Synthetic abundance"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "03_synthetic_immune_cell_dynamics.png"),
  p_cells,
  width = 8,
  height = 6,
  dpi = 300
)
