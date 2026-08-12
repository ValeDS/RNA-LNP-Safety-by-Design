# 04_compute_safety_by_design_score_consistent.R
# Consistent design-space-relative Safety-by-Design normalization
#
# Endpoint-specific min/max parameters are estimated once from the complete
# simulated patient–formulation design space and reused for every downstream
# representation, including the synthetic multi-omics subset.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source(file.path("scripts", "config.R"))
output_dir <- Sys.getenv(
  "RNA_LNP_SBD_OUTPUT_DIR",
  unset = file.path(project_dir, "results", "core")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

find_input <- function(filename) {
  candidates <- c(
    file.path(core_output_dir, filename),
    file.path(input_dir, filename),
    file.path(primary_dir, filename)
  )
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop(
      paste0(
        "Input file '", filename, "' not found. Expected one of:\n- ",
        paste(candidates, collapse = "\n- ")
      )
    )
  }

  existing[[1]]
}

response_file <- find_input(
  "02_patient_lnp_synthetic_immune_responses.csv"
)
metadata_file <- find_input("03_multiomics_sample_metadata.csv")
signature_file <- find_input("03_synthetic_pathway_signature_scores.csv")
cytokine_file <- find_input("03_synthetic_cytokine_layer.csv")
cell_file <- find_input("03_synthetic_immune_cell_layer.csv")

responses <- read.csv(response_file, stringsAsFactors = FALSE)
omics_samples <- read.csv(metadata_file, stringsAsFactors = FALSE)
signature_scores <- read.csv(signature_file, stringsAsFactors = FALSE)
cytokine_layer <- read.csv(cytokine_file, stringsAsFactors = FALSE)
cell_layer <- read.csv(cell_file, stringsAsFactors = FALSE)

endpoints <- c(
  "delivery_efficiency",
  "adaptive_activation",
  "innate_activation",
  "cytokine_burden",
  "off_target_activation"
)

missing_endpoints <- setdiff(endpoints, names(responses))

if (length(missing_endpoints) > 0) {
  stop(
    "Missing response endpoints: ",
    paste(missing_endpoints, collapse = ", ")
  )
}

# 1. Estimate normalization parameters once

normalization_parameters <- data.frame(
  endpoint = endpoints,
  minimum = vapply(
    responses[endpoints],
    min,
    numeric(1),
    na.rm = TRUE
  ),
  maximum = vapply(
    responses[endpoints],
    max,
    numeric(1),
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(range = maximum - minimum)

if (any(normalization_parameters$range <= 0)) {
  stop("Every endpoint must have a non-zero normalization range.")
}

scale_endpoint <- function(x, endpoint, parameters) {
  row <- parameters[parameters$endpoint == endpoint, , drop = FALSE]

  if (nrow(row) != 1) {
    stop("No unique normalization parameters found for: ", endpoint)
  }

  scaled <- (x - row$minimum) / row$range

  # Clipping protects against small numerical drift or future subset values
  # marginally outside the calibration design space.
  pmin(1, pmax(0, scaled))
}

add_normalized_endpoints <- function(data, parameters) {
  data %>%
    mutate(
      delivery_norm = scale_endpoint(
        delivery_efficiency,
        "delivery_efficiency",
        parameters
      ),
      adaptive_norm = scale_endpoint(
        adaptive_activation,
        "adaptive_activation",
        parameters
      ),
      innate_norm = scale_endpoint(
        innate_activation,
        "innate_activation",
        parameters
      ),
      cytokine_norm = scale_endpoint(
        cytokine_burden,
        "cytokine_burden",
        parameters
      ),
      off_target_norm = scale_endpoint(
        off_target_activation,
        "off_target_activation",
        parameters
      )
    )
}

add_sbd_score <- function(data) {
  data %>%
    mutate(
      safety_by_design_score =
        0.30 * delivery_norm +
        0.20 * adaptive_norm +
        0.20 * (1 - innate_norm) +
        0.15 * (1 - cytokine_norm) +
        0.15 * (1 - off_target_norm),
      risk_class = case_when(
        safety_by_design_score >= 0.75 ~ "Optimal",
        safety_by_design_score >= 0.55 ~ "Acceptable",
        safety_by_design_score >= 0.35 ~ "Caution",
        TRUE ~ "High-risk"
      )
    )
}

# 2. Score the complete simulated design space

responses_scored <- responses %>%
  add_normalized_endpoints(normalization_parameters) %>%
  add_sbd_score()

# 3. Score the multi-omics subset using the same parameters

multiomics_scored <- omics_samples %>%
  select(
    sample_id,
    pair_id,
    patient_id,
    lnp_id,
    all_of(endpoints),
    baseline_inflammation,
    inflammatory_risk_class
  ) %>%
  left_join(signature_scores, by = "sample_id") %>%
  left_join(cytokine_layer, by = "sample_id") %>%
  left_join(cell_layer, by = "sample_id") %>%
  add_normalized_endpoints(normalization_parameters) %>%
  add_sbd_score()

# 4. Formulation-level summaries

lnp_ranking <- responses_scored %>%
  group_by(
    lnp_id,
    size_nm,
    charge_mV,
    peg_mol_percent,
    ionizable_lipid_fraction,
    cholesterol_fraction,
    helper_lipid_fraction,
    targeting_ligand
  ) %>%
  summarise(
    mean_safety_score = mean(safety_by_design_score, na.rm = TRUE),
    sd_safety_score = sd(safety_by_design_score, na.rm = TRUE),
    mean_delivery = mean(delivery_efficiency, na.rm = TRUE),
    mean_adaptive_activation = mean(
      adaptive_activation,
      na.rm = TRUE
    ),
    mean_innate_activation = mean(innate_activation, na.rm = TRUE),
    mean_cytokine_burden = mean(cytokine_burden, na.rm = TRUE),
    mean_off_target_activation = mean(
      off_target_activation,
      na.rm = TRUE
    ),
    n_simulations = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_safety_score))

top20_lnp <- lnp_ranking %>%
  slice_head(n = 20)

# 5. Reproducibility metadata and exports

score_metadata <- data.frame(
  item = c(
    "normalization_scope",
    "calibration_records",
    "unique_virtual_patients",
    "unique_lnp_formulations",
    "score_weights",
    "random_seed"
  ),
  value = c(
    "empirical min-max over the complete simulated design space",
    nrow(responses_scored),
    n_distinct(responses_scored$patient_id),
    n_distinct(responses_scored$lnp_id),
    paste0(
      "delivery=0.30; adaptive=0.20; innate penalty=0.20; ",
      "cytokine penalty=0.15; off-target penalty=0.15"
    ),
    123
  )
)

write.csv(
  normalization_parameters,
  file.path(output_dir, "04_sbd_normalization_parameters.csv"),
  row.names = FALSE
)

write.csv(
  score_metadata,
  file.path(output_dir, "04_sbd_score_metadata.csv"),
  row.names = FALSE
)

write.csv(
  responses_scored,
  file.path(
    output_dir,
    "04_patient_lnp_safety_by_design_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  multiomics_scored,
  file.path(output_dir, "04_multiomics_safety_by_design_scores.csv"),
  row.names = FALSE
)

write.csv(
  lnp_ranking,
  file.path(output_dir, "04_lnp_ranking_by_safety_score.csv"),
  row.names = FALSE
)

write.csv(
  top20_lnp,
  file.path(
    output_dir,
    "04_top20_lnp_candidates_safety_by_design.csv"
  ),
  row.names = FALSE
)

# 6. Diagnostic figures

p_score <- ggplot(
  responses_scored,
  aes(x = safety_by_design_score)
) +
  geom_histogram(bins = 40, fill = "#3B6EA8", color = "white") +
  labs(
    title = "Distribution of Safety-by-Design scores",
    x = "Safety-by-Design score",
    y = "Simulated patient–formulation records"
  ) +
  theme_minimal(base_size = 13)

p_risk <- ggplot(
  responses_scored,
  aes(x = risk_class)
) +
  geom_bar(fill = "#4C956C") +
  labs(
    title = "Risk-class distribution",
    x = "Risk class",
    y = "Simulated records"
  ) +
  theme_minimal(base_size = 13)

p_delivery <- ggplot(
  responses_scored,
  aes(x = delivery_efficiency, y = safety_by_design_score)
) +
  geom_point(alpha = 0.30, color = "#3B6EA8") +
  labs(
    title = "Delivery efficiency and Safety-by-Design score",
    x = "Delivery efficiency",
    y = "Safety-by-Design score"
  ) +
  theme_minimal(base_size = 13)

p_top <- ggplot(
  top20_lnp,
  aes(
    x = reorder(lnp_id, mean_safety_score),
    y = mean_safety_score
  )
) +
  geom_col(fill = "#3B6EA8") +
  coord_flip() +
  labs(
    title = "Top 20 formulations by mean Safety-by-Design score",
    x = "RNA-LNP formulation",
    y = "Mean Safety-by-Design score"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(output_dir, "04_safety_score_distribution.png"),
  p_score,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(output_dir, "04_risk_class_distribution.png"),
  p_risk,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(output_dir, "04_delivery_vs_safety_score.png"),
  p_delivery,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(output_dir, "04_top20_lnp_candidates.png"),
  p_top,
  width = 8,
  height = 7,
  dpi = 300
)

message("Consistent Safety-by-Design scoring completed.")
message("Response input: ", normalizePath(response_file))
message("Output directory: ", normalizePath(output_dir))
message(
  "Score range: ",
  sprintf("%.4f", min(responses_scored$safety_by_design_score)),
  "-",
  sprintf("%.4f", max(responses_scored$safety_by_design_score))
)
