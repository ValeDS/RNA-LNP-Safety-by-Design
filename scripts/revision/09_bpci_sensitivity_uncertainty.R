suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

seed <- 20260812
set.seed(seed)

source(file.path("scripts", "config.R"))
output_dir <- revision_output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ranking <- read.csv(
  file.path(input_dir, "07c_GSE171964_ranking_comparison_by_day.csv"),
  stringsAsFactors = FALSE
)
scores <- read.csv(
  file.path(input_dir, "07_GSE171964_real_signature_scores_long.csv"),
  stringsAsFactors = FALSE
)
cells <- read.csv(
  file.path(input_dir, "07_GSE171964_cell_composition_by_patient_day.csv"),
  stringsAsFactors = FALSE
)

top_k <- 3
original_weights <- c(rank = 0.4, direction = 0.4, top3 = 0.2)

compute_bpci <- function(data, weights = original_weights) {
  stopifnot(abs(sum(weights) - 1) < 1e-8)
  data %>%
    group_by(day) %>%
    mutate(
      predicted_rank = rank(-predicted_activation_score, ties.method = "average"),
      observed_rank = rank(-delta_vs_baseline, ties.method = "average"),
      predicted_top_k = predicted_rank <= top_k,
      observed_top_k = observed_rank <= top_k,
      predicted_direction = sign(predicted_activation_score),
      observed_direction = sign(delta_vs_baseline)
    ) %>%
    summarise(
      spearman = cor(
        predicted_activation_score,
        delta_vs_baseline,
        method = "spearman",
        use = "complete.obs"
      ),
      rank_component = (spearman + 1) / 2,
      direction_component = mean(
        predicted_direction == observed_direction,
        na.rm = TRUE
      ),
      top3_component = length(
        intersect(signature[predicted_top_k], signature[observed_top_k])
      ) / top_k,
      .groups = "drop"
    ) %>%
    mutate(
      BPCI =
        weights[["rank"]] * rank_component +
        weights[["direction"]] * direction_component +
        weights[["top3"]] * top3_component
    )
}

# 1. Observed BPCI and exact shuffled-label null (5! = 120 permutations)

permutations <- function(x) {
  if (length(x) == 1) return(matrix(x, nrow = 1))
  do.call(
    rbind,
    lapply(seq_along(x), function(i) {
      rest <- permutations(x[-i])
      cbind(x[i], rest)
    })
  )
}

observed <- compute_bpci(ranking)
programs <- sort(unique(ranking$signature))
perm_matrix <- permutations(programs)

null_results <- lapply(unique(ranking$day), function(day_value) {
  day_data <- ranking %>%
    filter(day == day_value) %>%
    arrange(match(signature, programs))

  bind_rows(lapply(seq_len(nrow(perm_matrix)), function(i) {
    permuted <- day_data
    source_order <- match(perm_matrix[i, ], programs)
    permuted$predicted_activation_score <-
      day_data$predicted_activation_score[source_order]
    compute_bpci(permuted) %>%
      mutate(permutation_id = i)
  }))
}) %>% bind_rows()

null_summary <- null_results %>%
  group_by(day) %>%
  summarise(
    null_n = n(),
    null_mean = mean(BPCI),
    null_sd = sd(BPCI),
    null_q025 = quantile(BPCI, 0.025),
    null_median = median(BPCI),
    null_q975 = quantile(BPCI, 0.975),
    .groups = "drop"
  ) %>%
  left_join(observed %>% select(day, observed_BPCI = BPCI), by = "day") %>%
  left_join(
    null_results %>%
      left_join(observed %>% select(day, observed_BPCI = BPCI), by = "day") %>%
      group_by(day) %>%
      summarise(
        exact_upper_tail_p = mean(BPCI >= observed_BPCI - 1e-12),
        observed_percentile = mean(BPCI <= observed_BPCI + 1e-12),
        .groups = "drop"
      ),
    by = "day"
  )

global_null <- null_results %>%
  group_by(permutation_id) %>%
  summarise(mean_BPCI_across_days = mean(BPCI), .groups = "drop")

global_observed_mean <- mean(observed$BPCI)
global_null_summary <- data.frame(
  observed_mean_BPCI = global_observed_mean,
  null_n = nrow(global_null),
  null_mean = mean(global_null$mean_BPCI_across_days),
  null_sd = sd(global_null$mean_BPCI_across_days),
  null_q025 = quantile(global_null$mean_BPCI_across_days, 0.025),
  null_median = median(global_null$mean_BPCI_across_days),
  null_q975 = quantile(global_null$mean_BPCI_across_days, 0.975),
  exact_upper_tail_p = mean(
    global_null$mean_BPCI_across_days >= global_observed_mean - 1e-12
  ),
  observed_percentile = mean(
    global_null$mean_BPCI_across_days <= global_observed_mean + 1e-12
  )
)

# 2. Weight sensitivity over a plausible simplex (weights >= 0.10)

weight_grid <- expand.grid(
  weight_rank = seq(0.10, 0.80, 0.05),
  weight_direction = seq(0.10, 0.80, 0.05),
  stringsAsFactors = FALSE
) %>%
  mutate(weight_top3 = 1 - weight_rank - weight_direction) %>%
  filter(weight_top3 >= 0.10 - 1e-10, weight_top3 <= 0.80 + 1e-10) %>%
  mutate(
    weight_top3 = round(weight_top3, 10),
    weight_set = row_number()
  )

weight_sensitivity <- bind_rows(lapply(seq_len(nrow(weight_grid)), function(i) {
  w <- c(
    rank = weight_grid$weight_rank[i],
    direction = weight_grid$weight_direction[i],
    top3 = weight_grid$weight_top3[i]
  )
  compute_bpci(ranking, w) %>%
    mutate(
      weight_set = weight_grid$weight_set[i],
      weight_rank = w[["rank"]],
      weight_direction = w[["direction"]],
      weight_top3 = w[["top3"]]
    )
}))

original_day_order <- observed %>%
  arrange(desc(BPCI), day) %>%
  pull(day)

weight_summary <- weight_sensitivity %>%
  group_by(day) %>%
  summarise(
    n_weight_sets = n(),
    min_BPCI = min(BPCI),
    median_BPCI = median(BPCI),
    max_BPCI = max(BPCI),
    q025_BPCI = quantile(BPCI, 0.025),
    q975_BPCI = quantile(BPCI, 0.975),
    .groups = "drop"
  ) %>%
  left_join(observed %>% select(day, original_BPCI = BPCI), by = "day")

weight_order_stability <- weight_sensitivity %>%
  group_by(weight_set) %>%
  summarise(
    spearman_day_order_vs_original = cor(
      BPCI,
      observed$BPCI[match(day, observed$day)],
      method = "spearman"
    ),
    .groups = "drop"
  )

# 3. Donor/cell counts and paired-donor bootstrap uncertainty

cell_by_donor_day <- cells %>%
  group_by(day, pt_id) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

sample_sizes <- cell_by_donor_day %>%
  group_by(day) %>%
  summarise(
    n_donors = n_distinct(pt_id),
    total_cells = sum(n_cells),
    min_cells_per_donor = min(n_cells),
    median_cells_per_donor = median(n_cells),
    max_cells_per_donor = max(n_cells),
    .groups = "drop"
  )

paired_deltas <- scores %>%
  select(pt_id, day, signature, real_signature_score) %>%
  left_join(
    scores %>%
      filter(day == 0) %>%
      select(pt_id, signature, baseline_score = real_signature_score),
    by = c("pt_id", "signature")
  ) %>%
  filter(day != 0) %>%
  mutate(delta_vs_paired_baseline = real_signature_score - baseline_score)

predicted <- ranking %>%
  select(signature, predicted_activation_score) %>%
  distinct()

n_boot <- 2000
fast_bpci <- function(predicted_values, observed_values) {
  predicted_rank <- rank(-predicted_values, ties.method = "average")
  observed_rank <- rank(-observed_values, ties.method = "average")
  rank_component <- (cor(predicted_values, observed_values, method = "spearman") + 1) / 2
  direction_component <- mean(sign(predicted_values) == sign(observed_values))
  top3_component <- length(intersect(which(predicted_rank <= top_k), which(observed_rank <= top_k))) / top_k
  original_weights[["rank"]] * rank_component +
    original_weights[["direction"]] * direction_component +
    original_weights[["top3"]] * top3_component
}

bootstrap_results <- bind_rows(lapply(sort(unique(paired_deltas$day)), function(day_value) {
  day_data <- paired_deltas %>% filter(day == day_value)
  donor_ids <- sort(unique(day_data$pt_id))
  signature_ids <- sort(unique(day_data$signature))
  delta_matrix <- matrix(
    day_data$delta_vs_paired_baseline[
      match(
        as.vector(outer(donor_ids, signature_ids, paste, sep = "::")),
        paste(day_data$pt_id, day_data$signature, sep = "::")
      )
    ],
    nrow = length(donor_ids),
    ncol = length(signature_ids),
    byrow = FALSE
  )
  # Rebuild explicitly to avoid dependence on input row order.
  delta_matrix <- sapply(signature_ids, function(sig) {
    sig_data <- day_data[day_data$signature == sig, , drop = FALSE]
    sig_data$delta_vs_paired_baseline[match(donor_ids, sig_data$pt_id)]
  })
  predicted_values <- predicted$predicted_activation_score[
    match(signature_ids, predicted$signature)
  ]
  sampled_indices <- matrix(
    sample(seq_along(donor_ids), n_boot * length(donor_ids), replace = TRUE),
    nrow = n_boot
  )
  bpci_values <- apply(sampled_indices, 1, function(idx) {
    fast_bpci(predicted_values, colMeans(delta_matrix[idx, , drop = FALSE]))
  })
  data.frame(day = day_value, iteration = seq_len(n_boot), BPCI = bpci_values)
}))

bootstrap_summary <- bootstrap_results %>%
  group_by(day) %>%
  summarise(
    bootstrap_n = n(),
    bootstrap_valid_n = sum(!is.na(BPCI)),
    bootstrap_mean_BPCI = mean(BPCI, na.rm = TRUE),
    bootstrap_sd_BPCI = sd(BPCI, na.rm = TRUE),
    bootstrap_median_BPCI = median(BPCI, na.rm = TRUE),
    bootstrap_q025_BPCI = quantile(BPCI, 0.025, na.rm = TRUE),
    bootstrap_q975_BPCI = quantile(BPCI, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(sample_sizes, by = "day") %>%
  left_join(observed %>% select(day, original_unpaired_BPCI = BPCI), by = "day")

# 4. Null and bootstrap summary figure

p_null <- ggplot(null_results, aes(x = BPCI)) +
  geom_histogram(binwidth = 0.04, boundary = 0, fill = "#B8C4CE", color = "white") +
  geom_vline(
    data = observed,
    aes(xintercept = BPCI),
    color = "#B33A3A",
    linewidth = 0.7
  ) +
  facet_wrap(~ day, nrow = 2) +
  labs(
    title = "Observed BPCI relative to the exact shuffled-label null",
    subtitle = "Red lines show observed values; null distributions contain all 5! program-label permutations",
    x = "BPCI",
    y = "Number of permutations"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(output_dir, "09_BPCI_exact_null_by_day.png"),
  p_null,
  width = 10,
  height = 6.5,
  dpi = 350
)

# 5. Exports

write.csv(observed, file.path(output_dir, "09_BPCI_observed_components.csv"), row.names = FALSE)
write.csv(null_results, file.path(output_dir, "09_BPCI_exact_null_all_permutations.csv"), row.names = FALSE)
write.csv(null_summary, file.path(output_dir, "09_BPCI_exact_null_summary.csv"), row.names = FALSE)
write.csv(global_null, file.path(output_dir, "09_BPCI_exact_null_global_distribution.csv"), row.names = FALSE)
write.csv(global_null_summary, file.path(output_dir, "09_BPCI_exact_null_global_summary.csv"), row.names = FALSE)
write.csv(weight_grid, file.path(output_dir, "09_BPCI_weight_grid.csv"), row.names = FALSE)
write.csv(weight_sensitivity, file.path(output_dir, "09_BPCI_weight_sensitivity_all.csv"), row.names = FALSE)
write.csv(weight_summary, file.path(output_dir, "09_BPCI_weight_sensitivity_summary.csv"), row.names = FALSE)
write.csv(weight_order_stability, file.path(output_dir, "09_BPCI_weight_order_stability.csv"), row.names = FALSE)
write.csv(sample_sizes, file.path(output_dir, "09_GSE171964_sample_sizes_by_day.csv"), row.names = FALSE)
write.csv(bootstrap_summary, file.path(output_dir, "09_BPCI_donor_bootstrap_summary.csv"), row.names = FALSE)

metadata <- data.frame(
  item = c(
    "analysis_scope",
    "random_seed",
    "null_method",
    "null_permutations_per_day",
    "weight_sensitivity_domain",
    "bootstrap_method",
    "bootstrap_iterations",
    "temporal_interpretation"
  ),
  value = c(
    "BPCI null distribution, weight sensitivity, sample sizes, and uncertainty",
    seed,
    "exact enumeration of all label permutations for five programs",
    nrow(perm_matrix),
    "simplex grid in 0.05 increments; every component weight >= 0.10",
    "paired-donor resampling within each post-vaccination day; donor-specific Day 0 baseline",
    n_boot,
    "fixed synthetic five-program vector compared with observed day-specific deviations from baseline"
  )
)
write.csv(metadata, file.path(output_dir, "09_BPCI_revision_analysis_metadata.csv"), row.names = FALSE)

message("BPCI sensitivity analyses written to ", normalizePath(output_dir))
