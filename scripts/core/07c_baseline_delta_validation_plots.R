# 07c_baseline_delta_validation_plots.R
# Baseline-delta validation plots for GSE171964

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source(file.path("scripts", "config.R"))

real_scores_long <- read.csv(
  first_existing(file.path(core_output_dir, "07_GSE171964_real_signature_scores_long.csv"), file.path(input_dir, "07_GSE171964_real_signature_scores_long.csv"))
)

validation_comparison <- read.csv(
  first_existing(file.path(core_output_dir, "07_GSE171964_predicted_vs_observed_signature_comparison.csv"), file.path(input_dir, "07_GSE171964_predicted_vs_observed_signature_comparison.csv"))
)

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

minmax_norm <- function(x) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    return(rep(0, length(x)))
  }
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# 1. Compute observed delta from Day 0 baseline

signature_day_means <- real_scores_long %>%
  mutate(signature_label = pretty_signature(signature)) %>%
  group_by(day, signature, signature_label) %>%
  summarise(
    mean_signature_score = mean(real_signature_score, na.rm = TRUE),
    .groups = "drop"
  )

baseline_means <- signature_day_means %>%
  filter(day == 0) %>%
  select(signature, baseline_mean = mean_signature_score)

delta_df <- signature_day_means %>%
  left_join(baseline_means, by = "signature") %>%
  mutate(
    delta_vs_baseline = mean_signature_score - baseline_mean
  )

write.csv(
  delta_df,
  file.path(core_output_dir, "07c_GSE171964_signature_delta_vs_baseline_by_day.csv"),
  row.names = FALSE
)

# 2. Heatmap of observed delta from baseline

delta_heatmap_df <- delta_df %>%
  group_by(signature) %>%
  mutate(
    delta_zscore = as.numeric(scale(delta_vs_baseline))
  ) %>%
  ungroup()

write.csv(
  delta_heatmap_df,
  file.path(core_output_dir, "07c_GSE171964_signature_delta_heatmap_data.csv"),
  row.names = FALSE
)

p_delta_heat <- ggplot(
  delta_heatmap_df,
  aes(
    x = factor(day),
    y = signature_label,
    fill = delta_zscore
  )
) +
  geom_tile(color = "white", linewidth = 0.3) +
  labs(
    title = "Observed immune-signature deviation from baseline in GSE171964",
    x = "Day",
    y = "Immune signature",
    fill = "Delta\nz-score"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 11)
  )

ggsave(
  file.path(core_output_dir, "07c_GSE171964_signature_delta_vs_baseline_heatmap.png"),
  p_delta_heat,
  width = 10,
  height = 5.5,
  dpi = 300
)

# 3. Line plot of observed delta from baseline

p_delta_line <- ggplot(
  delta_df,
  aes(
    x = day,
    y = delta_vs_baseline,
    group = signature_label,
    color = signature_label
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  labs(
    title = "Observed pathway-level deviations from baseline",
    x = "Day",
    y = "Mean signature score − Day 0 baseline",
    color = "Immune signature"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(core_output_dir, "07c_GSE171964_signature_delta_vs_baseline_lineplot.png"),
  p_delta_line,
  width = 10,
  height = 6,
  dpi = 300
)

# 4. Ranking-based validation against synthetic predictions

predicted_summary <- validation_comparison %>%
  select(signature, predicted_activation_score) %>%
  distinct()

ranking_comparison <- delta_df %>%
  filter(day != 0) %>%
  left_join(predicted_summary, by = "signature") %>%
  group_by(day) %>%
  mutate(
    observed_rank = rank(-delta_vs_baseline, ties.method = "average"),
    predicted_rank = rank(-predicted_activation_score, ties.method = "average")
  ) %>%
  ungroup()

write.csv(
  ranking_comparison,
  file.path(core_output_dir, "07c_GSE171964_ranking_comparison_by_day.csv"),
  row.names = FALSE
)

ranking_concordance <- ranking_comparison %>%
  group_by(day) %>%
  summarise(
    spearman_rank_concordance = cor(
      predicted_activation_score,
      delta_vs_baseline,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

write.csv(
  ranking_concordance,
  file.path(core_output_dir, "07c_GSE171964_rank_concordance_by_day.csv"),
  row.names = FALSE
)

p_rank <- ggplot(
  ranking_concordance,
  aes(x = factor(day), y = spearman_rank_concordance)
) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Ranking concordance between synthetic predictions and observed baseline deviations",
    x = "Day",
    y = "Spearman rank concordance"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "07c_GSE171964_rank_concordance_by_day.png"),
  p_rank,
  width = 8,
  height = 5,
  dpi = 300
)

# 5. Predicted vs observed delta scatter by day

scatter_df <- ranking_comparison %>%
  mutate(
    signature_label = pretty_signature(signature)
  ) %>%
  group_by(day) %>%
  mutate(
    predicted_norm = minmax_norm(predicted_activation_score),
    observed_delta_norm = minmax_norm(delta_vs_baseline)
  ) %>%
  ungroup()

write.csv(
  scatter_df,
  file.path(core_output_dir, "07c_GSE171964_predicted_vs_observed_delta_scatter_data.csv"),
  row.names = FALSE
)

p_delta_scatter <- ggplot(
  scatter_df,
  aes(
    x = predicted_norm,
    y = observed_delta_norm,
    label = signature_label
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2.8) +
  geom_text(vjust = -0.7, size = 3, check_overlap = TRUE) +
  facet_wrap(~ day) +
  coord_cartesian(
    xlim = c(-0.08, 1.12),
    ylim = c(-0.08, 1.12),
    clip = "off"
  ) +
  labs(
    title = "Predicted activation versus observed deviation from baseline",
    x = "Normalized predicted synthetic activation",
    y = "Normalized observed delta from baseline"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 45, 10, 45)
  )

ggsave(
  file.path(core_output_dir, "07c_GSE171964_predicted_vs_observed_delta_scatter_by_day.png"),
  p_delta_scatter,
  width = 12,
  height = 7,
  dpi = 300
)

message("Baseline-delta validation outputs generated.")
