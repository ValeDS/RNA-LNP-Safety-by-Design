# 08_biological_program_concordance_index_corrected.R
# Correct implementation of the Biological Program Concordance Index
#
# BPCI = 0.4 R + 0.4 D + 0.2 T
#   R = rescaled Spearman rank concordance, (rho + 1) / 2
#   D = activation-direction agreement
#   T = overlap between predicted and observed top-3 programs / 3

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source(file.path("scripts", "config.R"))
output_dir <- Sys.getenv("RNA_LNP_BPCI_OUTPUT_DIR", unset = revision_output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

top_k <- 3
weight_rank <- 0.4
weight_direction <- 0.4
weight_top_k <- 0.2

if (
  abs(weight_rank + weight_direction + weight_top_k - 1) >
    .Machine$double.eps^0.5
) {
  stop("BPCI weights must sum to 1.")
}

ranking_candidates <- c(file.path(input_dir, "07c_GSE171964_ranking_comparison_by_day.csv"))

existing_inputs <- ranking_candidates[file.exists(ranking_candidates)]

if (length(existing_inputs) == 0) {
  stop(
    paste0(
      "Ranking-comparison input not found. Expected one of:\n- ",
      paste(ranking_candidates, collapse = "\n- ")
    )
  )
}

input_file <- existing_inputs[[1]]
ranking_comparison <- read.csv(input_file, stringsAsFactors = FALSE)

required_columns <- c(
  "day",
  "signature",
  "predicted_activation_score",
  "delta_vs_baseline"
)

missing_columns <- setdiff(required_columns, names(ranking_comparison))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

pretty_signature <- function(x) {
  dplyr::recode(
    x,
    interferon_response = "Interferon response",
    inflammatory_response = "Inflammatory response",
    antigen_presentation = "Antigen presentation",
    adaptive_immunity = "Adaptive immunity",
    stress_response = "Stress response"
  )
}

# Synthetic activation scores represent post-administration activation states.
# Their direction is therefore positive unless a score is exactly zero.
# Observed direction is defined from the deviation relative to Day 0.
comparison_ranked <- ranking_comparison %>%
  group_by(day) %>%
  mutate(
    predicted_rank = rank(
      -predicted_activation_score,
      ties.method = "average"
    ),
    observed_rank = rank(
      -delta_vs_baseline,
      ties.method = "average"
    ),
    predicted_top_k = predicted_rank <= top_k,
    observed_top_k = observed_rank <= top_k,
    predicted_direction = sign(predicted_activation_score),
    observed_direction = sign(delta_vs_baseline),
    direction_match = predicted_direction == observed_direction
  ) %>%
  ungroup()

bpci_components <- comparison_ranked %>%
  group_by(day) %>%
  summarise(
    n_programs = n_distinct(signature),
    spearman_rank_concordance = cor(
      predicted_activation_score,
      delta_vs_baseline,
      method = "spearman",
      use = "complete.obs"
    ),
    top_k_overlap_count = length(
      intersect(
        signature[predicted_top_k],
        signature[observed_top_k]
      )
    ),
    top3_agreement = top_k_overlap_count / top_k,
    direction_agreement = mean(direction_match, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    spearman_scaled = (spearman_rank_concordance + 1) / 2,
    BPCI =
      weight_rank * spearman_scaled +
      weight_direction * direction_agreement +
      weight_top_k * top3_agreement
  ) %>%
  select(
    day,
    n_programs,
    spearman_rank_concordance,
    spearman_scaled,
    direction_agreement,
    top_k_overlap_count,
    top3_agreement,
    BPCI
  )

bpci_summary <- bpci_components %>%
  summarise(
    n_time_points = n(),
    mean_spearman = mean(spearman_rank_concordance, na.rm = TRUE),
    mean_rank_component = mean(spearman_scaled, na.rm = TRUE),
    mean_direction_agreement = mean(direction_agreement, na.rm = TRUE),
    mean_top3_agreement = mean(top3_agreement, na.rm = TRUE),
    mean_BPCI = mean(BPCI, na.rm = TRUE),
    sd_BPCI = sd(BPCI, na.rm = TRUE),
    min_BPCI = min(BPCI, na.rm = TRUE),
    max_BPCI = max(BPCI, na.rm = TRUE)
  )

bpci_long <- bpci_components %>%
  select(
    day,
    spearman_scaled,
    direction_agreement,
    top3_agreement,
    BPCI
  ) %>%
  pivot_longer(
    cols = -day,
    names_to = "component",
    values_to = "value"
  ) %>%
  mutate(
    component = recode(
      component,
      spearman_scaled = "Rank concordance",
      direction_agreement = "Direction agreement",
      top3_agreement = "Top-3 agreement",
      BPCI = "BPCI"
    )
  )

mean_rank_df <- comparison_ranked %>%
  mutate(signature_label = pretty_signature(signature)) %>%
  group_by(signature, signature_label) %>%
  summarise(
    predicted_rank_mean = mean(predicted_rank, na.rm = TRUE),
    observed_rank_mean = mean(observed_rank, na.rm = TRUE),
    mean_delta_vs_baseline = mean(delta_vs_baseline, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(observed_rank_mean)

# Rank agreement must be calculated within each day.
agreement_df <- comparison_ranked %>%
  group_by(day) %>%
  mutate(
    rank_distance = abs(predicted_rank - observed_rank),
    rank_agreement_score = 1 - rank_distance / (n_distinct(signature) - 1)
  ) %>%
  ungroup() %>%
  mutate(signature_label = pretty_signature(signature))

metadata <- data.frame(
  item = c(
    "input_file",
    "top_k",
    "weight_rank",
    "weight_direction",
    "weight_top_k",
    "direction_definition"
  ),
  value = c(
    normalizePath(input_file),
    top_k,
    weight_rank,
    weight_direction,
    weight_top_k,
    paste0(
      "sign(predicted activation score) compared with ",
      "sign(observed delta from Day 0)"
    )
  )
)

write.csv(
  bpci_components,
  file.path(output_dir, "08_GSE171964_BPCI_by_day.csv"),
  row.names = FALSE
)

write.csv(
  bpci_summary,
  file.path(output_dir, "08_GSE171964_BPCI_summary.csv"),
  row.names = FALSE
)

write.csv(
  bpci_long,
  file.path(output_dir, "08_GSE171964_BPCI_components_long.csv"),
  row.names = FALSE
)

write.csv(
  mean_rank_df,
  file.path(
    output_dir,
    "08_GSE171964_mean_predicted_observed_ranking.csv"
  ),
  row.names = FALSE
)

write.csv(
  agreement_df,
  file.path(
    output_dir,
    "08_GSE171964_signature_level_rank_agreement.csv"
  ),
  row.names = FALSE
)

write.csv(
  metadata,
  file.path(output_dir, "08_GSE171964_BPCI_metadata.csv"),
  row.names = FALSE
)

message("Corrected BPCI analysis completed.")
message("Input: ", normalizePath(input_file))
message("Output directory: ", normalizePath(output_dir))
message(
  "BPCI range: ",
  sprintf("%.3f", bpci_summary$min_BPCI),
  "-",
  sprintf("%.3f", bpci_summary$max_BPCI),
  "; mean: ",
  sprintf("%.3f", bpci_summary$mean_BPCI)
)
