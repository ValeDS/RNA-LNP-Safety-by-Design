# 06_multiobjective_lnp_optimization.R
# Module 6: Multi-objective LNP candidate prioritization

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source(file.path("scripts", "config.R"))

set.seed(123)

data <- read.csv(
  first_existing(
    file.path(core_output_dir, "04_patient_lnp_safety_by_design_scores.csv"),
    file.path(primary_dir, "04_patient_lnp_safety_by_design_scores.csv")
  )
)

# Aggregate performance per LNP formulation

lnp_summary <- data %>%
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
    mean_delivery = mean(delivery_efficiency, na.rm = TRUE),
    mean_adaptive_activation = mean(adaptive_activation, na.rm = TRUE),
    mean_innate_activation = mean(innate_activation, na.rm = TRUE),
    mean_cytokine_burden = mean(cytokine_burden, na.rm = TRUE),
    mean_off_target_activation = mean(off_target_activation, na.rm = TRUE),
    .groups = "drop"
  )

# Pareto front helper function
# Maximize delivery and adaptive activation
# Minimize innate activation, cytokine burden, off-target activation

is_dominated <- function(i, df) {
  candidate <- df[i, ]

  others <- df[-i, ]

  dominated <- any(
    others$mean_delivery >= candidate$mean_delivery &
      others$mean_adaptive_activation >= candidate$mean_adaptive_activation &
      others$mean_innate_activation <= candidate$mean_innate_activation &
      others$mean_cytokine_burden <= candidate$mean_cytokine_burden &
      others$mean_off_target_activation <= candidate$mean_off_target_activation &
      (
        others$mean_delivery > candidate$mean_delivery |
          others$mean_adaptive_activation > candidate$mean_adaptive_activation |
          others$mean_innate_activation < candidate$mean_innate_activation |
          others$mean_cytokine_burden < candidate$mean_cytokine_burden |
          others$mean_off_target_activation < candidate$mean_off_target_activation
      )
  )

  return(dominated)
}

pareto_flags <- sapply(seq_len(nrow(lnp_summary)), is_dominated, df = lnp_summary)

lnp_summary <- lnp_summary %>%
  mutate(
    pareto_optimal = !pareto_flags
  ) %>%
  arrange(desc(pareto_optimal), desc(mean_safety_score))

write.csv(lnp_summary, file.path(core_output_dir, "06_lnp_multiobjective_summary.csv"), row.names = FALSE)

pareto_candidates <- lnp_summary %>%
  filter(pareto_optimal) %>%
  arrange(desc(mean_safety_score))

write.csv(pareto_candidates, file.path(core_output_dir, "06_pareto_optimal_lnp_candidates.csv"), row.names = FALSE)

top_pareto <- pareto_candidates %>%
  slice_head(n = 20)

write.csv(top_pareto, file.path(core_output_dir, "06_top20_pareto_lnp_candidates.csv"), row.names = FALSE)

# Visualizations

p_pareto <- ggplot(
  lnp_summary,
  aes(
    x = mean_delivery,
    y = mean_cytokine_burden,
    shape = pareto_optimal
  )
) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Multi-objective LNP prioritization",
    x = "Mean delivery efficiency",
    y = "Mean cytokine burden",
    shape = "Pareto optimal"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(core_output_dir, "06_pareto_delivery_vs_cytokine.png"), p_pareto, width = 8, height = 6, dpi = 300)

p_safety_offtarget <- ggplot(
  lnp_summary,
  aes(
    x = mean_safety_score,
    y = mean_off_target_activation,
    shape = pareto_optimal
  )
) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Safety score versus off-target activation",
    x = "Mean safety-by-design score",
    y = "Mean off-target activation",
    shape = "Pareto optimal"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(core_output_dir, "06_safety_vs_offtarget_pareto.png"), p_safety_offtarget, width = 8, height = 6, dpi = 300)

p_top <- ggplot(
  top_pareto,
  aes(x = reorder(lnp_id, mean_safety_score), y = mean_safety_score)
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top Pareto-optimal LNP candidates",
    x = "LNP candidate",
    y = "Mean safety-by-design score"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(core_output_dir, "06_top_pareto_candidates.png"), p_top, width = 8, height = 7, dpi = 300)
