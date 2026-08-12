# 01_generate_virtual_patients_and_lnp_space.R
# Safety-by-design RNA-LNP workflow
# Module 1: Virtual patients + LNP design space

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

set.seed(123)

source(file.path("scripts", "config.R"))
dir.create(core_output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Generate virtual patient cohort

generate_virtual_patients <- function(n_patients = 200, seed = 123) {
  set.seed(seed)

  patients <- data.frame(
    patient_id = paste0("VP_", sprintf("%03d", 1:n_patients)),

    age = round(runif(n_patients, min = 18, max = 85)),

    sex = sample(
      c("Female", "Male"),
      size = n_patients,
      replace = TRUE,
      prob = c(0.52, 0.48)
    ),

    baseline_inflammation = rbeta(n_patients, shape1 = 2, shape2 = 6),

    immune_responsiveness = rnorm(n_patients, mean = 1, sd = 0.20),

    innate_reactivity = rnorm(n_patients, mean = 1, sd = 0.25),

    adaptive_potential = rnorm(n_patients, mean = 1, sd = 0.20)
  )

  # Avoid biologically implausible negative values
  patients <- patients %>%
    mutate(
      immune_responsiveness = pmax(immune_responsiveness, 0.2),
      innate_reactivity = pmax(innate_reactivity, 0.2),
      adaptive_potential = pmax(adaptive_potential, 0.2)
    )

  # Simple risk class useful for later stratified analyses
  patients <- patients %>%
    mutate(
      inflammatory_risk_class = case_when(
        baseline_inflammation < 0.25 ~ "Low",
        baseline_inflammation < 0.50 ~ "Intermediate",
        TRUE ~ "High"
      )
    )

  return(patients)
}


# 2. Generate RNA-LNP design space

generate_lnp_design_space <- function(n_lnp = 500, seed = 123) {
  set.seed(seed)

  lnp <- data.frame(
    lnp_id = paste0("LNP_", sprintf("%04d", 1:n_lnp)),

    size_nm = runif(n_lnp, min = 50, max = 150),

    charge_mV = runif(n_lnp, min = -10, max = 10),

    peg_mol_percent = runif(n_lnp, min = 0.1, max = 0.5),

    ionizable_lipid_fraction = runif(n_lnp, min = 0.35, max = 0.60),

    cholesterol_fraction = runif(n_lnp, min = 0.25, max = 0.45),

    helper_lipid_fraction = runif(n_lnp, min = 0.05, max = 0.20),

    targeting_ligand = sample(
      c(0, 1),
      size = n_lnp,
      replace = TRUE,
      prob = c(0.60, 0.40)
    )
  )

  # Normalize lipid fractions to sum to 1, excluding PEG and targeting
  lipid_sum <- lnp$ionizable_lipid_fraction +
    lnp$cholesterol_fraction +
    lnp$helper_lipid_fraction

  lnp <- lnp %>%
    mutate(
      ionizable_lipid_fraction = ionizable_lipid_fraction / lipid_sum,
      cholesterol_fraction = cholesterol_fraction / lipid_sum,
      helper_lipid_fraction = helper_lipid_fraction / lipid_sum
    )

  return(lnp)
}


# 3. Generate objects

virtual_patients <- generate_virtual_patients(n_patients = 200, seed = 123)
lnp_space <- generate_lnp_design_space(n_lnp = 500, seed = 123)


# 4. Export tables

write.csv(
  virtual_patients,
  file = file.path(core_output_dir, "01_virtual_patient_cohort.csv"),
  row.names = FALSE
)

write.csv(
  lnp_space,
  file = file.path(core_output_dir, "01_lnp_design_space.csv"),
  row.names = FALSE
)


# 5. Quick exploratory plots

p_age <- ggplot(virtual_patients, aes(x = age)) +
  geom_histogram(bins = 25) +
  labs(
    title = "Virtual patient cohort: age distribution",
    x = "Age",
    y = "Number of virtual patients"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "01_virtual_patients_age_distribution.png"),
  p_age,
  width = 8,
  height = 5,
  dpi = 300
)

p_inflammation <- ggplot(
  virtual_patients,
  aes(x = inflammatory_risk_class, y = baseline_inflammation)
) +
  geom_boxplot() +
  labs(
    title = "Baseline inflammatory status by virtual patient risk class",
    x = "Inflammatory risk class",
    y = "Baseline inflammation score"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "01_virtual_patients_inflammation_risk.png"),
  p_inflammation,
  width = 8,
  height = 5,
  dpi = 300
)

p_lnp <- ggplot(lnp_space, aes(x = size_nm, y = charge_mV)) +
  geom_point(alpha = 0.6) +
  labs(
    title = "RNA-LNP design space",
    x = "Particle size (nm)",
    y = "Surface charge (mV)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(core_output_dir, "01_lnp_design_space_size_charge.png"),
  p_lnp,
  width = 8,
  height = 6,
  dpi = 300
)
