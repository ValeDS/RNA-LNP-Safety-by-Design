core_scripts <- c(
  "01_generate_virtual_patients_and_lnp_space.R",
  "02_simulate_patient_lnp_responses.R",
  "03_generate_synthetic_multiomics.R",
  "04_compute_safety_by_design_score_consistent.R",
  "05_train_random_forest_grouped_cv.R",
  "06_multiobjective_lnp_optimization.R"
)

for (script in core_scripts) {
  source(file.path("scripts", "core", script), local = new.env(parent = globalenv()))
}
