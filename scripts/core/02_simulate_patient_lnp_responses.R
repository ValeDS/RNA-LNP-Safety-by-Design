# 02_simulate_patient_lnp_responses.R
# Safety-by-design RNA-LNP workflow
# Module 2: Patient × LNP sampling + synthetic immune response

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

set.seed(123)

source(file.path("scripts", "config.R"))
dir.create(core_output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Module 1 outputs

virtual_patients <- read.csv(file.path(core_output_dir, "01_virtual_patient_cohort.csv"))
lnp_space <- read.csv(file.path(core_output_dir, "01_lnp_design_space.csv"))

# 2. Sample patient × LNP combinations

generate_patient_lnp_pairs <- function(
    patients,
    lnp,
    n_pairs = 10000,
    seed = 123
) {
  set.seed(seed)

  pairs <- data.frame(
    pair_id = paste0("PAIR_", sprintf("%05d", 1:n_pairs)),
    patient_id = sample(patients$patient_id, n_pairs, replace = TRUE),
    lnp_id = sample(lnp$lnp_id, n_pairs, replace = TRUE)
  )

  pairs <- pairs %>%
    left_join(patients, by = "patient_id") %>%
    left_join(lnp, by = "lnp_id")

  return(pairs)
}

patient_lnp_pairs <- generate_patient_lnp_pairs(
  virtual_patients,
  lnp_space,
  n_pairs = 10000,
  seed = 123
)

# 3. Synthetic immune response model

simulate_immune_response <- function(df, seed = 123) {
  set.seed(seed)

  response <- df %>%
    mutate(
      # Delivery efficiency: optimal around 90 nm, near-neutral charge,
      # moderate PEGylation, targeting ligand, balanced lipid composition.
      delivery_efficiency =
        1.8 -
        0.0008 * (size_nm - 90)^2 -
        0.015  * (charge_mV)^2 +
        0.9    * peg_mol_percent +
        0.7    * targeting_ligand +
        0.5    * ionizable_lipid_fraction +
        rnorm(n(), mean = 0, sd = 0.10),

      # Innate activation: increases with positive charge,
      # high innate reactivity, baseline inflammation, and suboptimal PEG.
      innate_activation =
        0.5 +
        0.08 * pmax(charge_mV, 0) +
        1.2  * baseline_inflammation +
        0.9  * innate_reactivity -
        0.6  * peg_mol_percent +
        rnorm(n(), mean = 0, sd = 0.12),

      # Adaptive activation: driven by delivery and individual adaptive potential.
      adaptive_activation =
        0.4 +
        0.9 * delivery_efficiency +
        0.8 * adaptive_potential +
        0.3 * targeting_ligand +
        rnorm(n(), mean = 0, sd = 0.12),

      # Cytokine burden: safety-related endpoint.
      cytokine_burden =
        0.3 +
        0.8 * innate_activation +
        0.4 * baseline_inflammation +
        0.2 * pmax(charge_mV, 0) +
        rnorm(n(), mean = 0, sd = 0.15),

      # Off-target activation: penalized by targeting,
      # increased by high charge and inflammatory background.
      off_target_activation =
        0.4 +
        0.5 * abs(charge_mV) +
        0.8 * baseline_inflammation -
        0.6 * targeting_ligand +
        rnorm(n(), mean = 0, sd = 0.15)
    )

  response <- response %>%
    mutate(
      across(
        c(
          delivery_efficiency,
          innate_activation,
          adaptive_activation,
          cytokine_burden,
          off_target_activation
        ),
        ~ pmax(.x, 0)
      )
    )

  return(response)
}

simulated_responses <- simulate_immune_response(patient_lnp_pairs, seed = 123)

# 4. Export simulated responses

write.csv(
  simulated_responses,
  file.path(core_output_dir, "02_patient_lnp_synthetic_immune_responses.csv"),
  row.names = FALSE
)

# 5. Exploratory plots

p_delivery <- ggplot(
  simulated_responses,
  aes(x = delivery_efficiency)
) +
  geom_histogram(bins = 40) +
  labs(
    title = "Distribution of predicted delivery efficiency",
    x = "Delivery efficiency",
    y = "Number of patient-LNP pairs"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "02_delivery_efficiency_distribution.png"),
  p_delivery,
  width = 8,
  height = 5,
  dpi = 300
)

p_safety <- ggplot(
  simulated_responses,
  aes(x = innate_activation, y = cytokine_burden)
) +
  geom_point(alpha = 0.4) +
  labs(
    title = "Innate activation and cytokine burden",
    x = "Innate activation",
    y = "Cytokine burden"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "02_innate_activation_vs_cytokine_burden.png"),
  p_safety,
  width = 8,
  height = 6,
  dpi = 300
)

p_lnp_response <- ggplot(
  simulated_responses,
  aes(x = size_nm, y = delivery_efficiency)
) +
  geom_point(alpha = 0.35) +
  labs(
    title = "Effect of LNP size on predicted delivery efficiency",
    x = "Particle size (nm)",
    y = "Delivery efficiency"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "02_lnp_size_vs_delivery_efficiency.png"),
  p_lnp_response,
  width = 8,
  height = 6,
  dpi = 300
)
