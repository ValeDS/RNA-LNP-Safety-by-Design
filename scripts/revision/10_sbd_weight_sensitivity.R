suppressPackageStartupMessages(library(dplyr))

source(file.path("scripts", "config.R"))
output_dir <- revision_output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scored <- read.csv(
  file.path(primary_dir, "04_patient_lnp_safety_by_design_scores.csv"),
  stringsAsFactors = FALSE
)

components <- scored %>%
  transmute(
    lnp_id,
    delivery = delivery_norm,
    adaptive = adaptive_norm,
    innate_safety = 1 - innate_norm,
    cytokine_safety = 1 - cytokine_norm,
    off_target_safety = 1 - off_target_norm
  ) %>%
  group_by(lnp_id) %>%
  summarise(
    across(
      c(delivery, adaptive, innate_safety, cytokine_safety, off_target_safety),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_records = n(),
    .groups = "drop"
  )

original_weights <- c(0.30, 0.20, 0.20, 0.15, 0.15)
component_names <- c(
  "delivery", "adaptive", "innate_safety", "cytokine_safety", "off_target_safety"
)

component_matrix <- as.matrix(components[, component_names])
rownames(component_matrix) <- components$lnp_id
original_scores <- as.numeric(component_matrix %*% original_weights)
names(original_scores) <- components$lnp_id
original_rank <- rank(-original_scores, ties.method = "min")

# Complete 0.05-step simplex with every component weight >= 0.05.
units_total <- 20
min_units <- 1
grid_units <- expand.grid(
  delivery = min_units:(units_total - 4 * min_units),
  adaptive = min_units:(units_total - 4 * min_units),
  innate_safety = min_units:(units_total - 4 * min_units),
  cytokine_safety = min_units:(units_total - 4 * min_units),
  stringsAsFactors = FALSE
) %>%
  mutate(off_target_safety = units_total - delivery - adaptive - innate_safety - cytokine_safety) %>%
  filter(off_target_safety >= min_units) %>%
  distinct()

weight_grid <- grid_units / units_total
weight_grid$weight_set <- seq_len(nrow(weight_grid))
weight_matrix <- as.matrix(weight_grid[, component_names])

score_matrix <- component_matrix %*% t(weight_matrix)
rank_matrix <- apply(-score_matrix, 2, rank, ties.method = "min")

original_top10 <- names(sort(original_scores, decreasing = TRUE))[1:10]
original_top20 <- names(sort(original_scores, decreasing = TRUE))[1:20]
original_best <- original_top10[1]

weight_set_summary <- data.frame(
  weight_set = weight_grid$weight_set,
  spearman_rank_vs_original = apply(rank_matrix, 2, function(x) cor(x, original_rank, method = "spearman")),
  top10_overlap = apply(rank_matrix, 2, function(x) sum(names(original_rank)[x <= 10] %in% original_top10)),
  top20_overlap = apply(rank_matrix, 2, function(x) sum(names(original_rank)[x <= 20] %in% original_top20)),
  original_best_rank = rank_matrix[match(original_best, rownames(rank_matrix)), ],
  best_lnp_id = apply(score_matrix, 2, function(x) rownames(score_matrix)[which.max(x)]),
  stringsAsFactors = FALSE
) %>%
  left_join(weight_grid, by = "weight_set")

candidate_summary <- data.frame(
  lnp_id = components$lnp_id,
  original_score = original_scores[components$lnp_id],
  original_rank = original_rank[components$lnp_id],
  median_rank = apply(rank_matrix, 1, median),
  min_rank = apply(rank_matrix, 1, min),
  max_rank = apply(rank_matrix, 1, max),
  top1_frequency = rowMeans(rank_matrix == 1),
  top10_frequency = rowMeans(rank_matrix <= 10),
  top20_frequency = rowMeans(rank_matrix <= 20),
  median_score = apply(score_matrix, 1, median),
  min_score = apply(score_matrix, 1, min),
  max_score = apply(score_matrix, 1, max),
  stringsAsFactors = FALSE
) %>%
  arrange(original_rank)

overall_summary <- data.frame(
  metric = c(
    "number_of_weight_sets",
    "minimum_rank_correlation",
    "median_rank_correlation",
    "maximum_rank_correlation",
    "minimum_top10_overlap",
    "median_top10_overlap",
    "minimum_top20_overlap",
    "median_top20_overlap",
    "number_of_distinct_best_formulations",
    "original_best_formulation",
    "original_best_worst_rank"
  ),
  value = c(
    nrow(weight_grid),
    min(weight_set_summary$spearman_rank_vs_original),
    median(weight_set_summary$spearman_rank_vs_original),
    max(weight_set_summary$spearman_rank_vs_original),
    min(weight_set_summary$top10_overlap),
    median(weight_set_summary$top10_overlap),
    min(weight_set_summary$top20_overlap),
    median(weight_set_summary$top20_overlap),
    length(unique(weight_set_summary$best_lnp_id)),
    original_best,
    max(weight_set_summary$original_best_rank)
  ),
  stringsAsFactors = FALSE
)

best_frequency <- weight_set_summary %>%
  count(best_lnp_id, name = "number_of_weight_sets") %>%
  mutate(frequency = number_of_weight_sets / nrow(weight_set_summary)) %>%
  arrange(desc(frequency))

write.csv(weight_grid, file.path(output_dir, "10_SbD_weight_grid.csv"), row.names = FALSE)
write.csv(weight_set_summary, file.path(output_dir, "10_SbD_weight_set_summary.csv"), row.names = FALSE)
write.csv(candidate_summary, file.path(output_dir, "10_SbD_candidate_robustness.csv"), row.names = FALSE)
write.csv(overall_summary, file.path(output_dir, "10_SbD_weight_sensitivity_overall.csv"), row.names = FALSE)
write.csv(best_frequency, file.path(output_dir, "10_SbD_best_candidate_frequency.csv"), row.names = FALSE)

metadata <- data.frame(
  item = c(
    "analysis_scope",
    "weight_step",
    "minimum_component_weight",
    "number_of_weight_sets",
    "ranking_level",
    "original_weights"
  ),
  value = c(
    "Safety-by-Design component-weight sensitivity",
    0.05,
    0.05,
    nrow(weight_grid),
    "mean normalized endpoint components per formulation",
    "delivery=0.30; adaptive=0.20; innate safety=0.20; cytokine safety=0.15; off-target safety=0.15"
  )
)
write.csv(metadata, file.path(output_dir, "10_SbD_weight_sensitivity_metadata.csv"), row.names = FALSE)

message("SbD sensitivity completed with ", nrow(weight_grid), " weight sets.")
