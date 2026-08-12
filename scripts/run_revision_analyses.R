revision_scripts <- c(
  "09_bpci_sensitivity_uncertainty.R",
  "10_sbd_weight_sensitivity.R",
  "11_inflammatory_risk_robustness.R",
  "12_pareto_candidate_audit.R",
  "13_noise_covariance_audit.R",
  "14_regenerate_independent_multiomics_layers.R",
  "15_sampling_grouped_cv_audit.R",
  "16_balanced_full_factorial_sensitivity.R",
  "17_generate_supplementary_revision_figures.R"
)

for (script in revision_scripts) {
  source(file.path("scripts", "revision", script), local = new.env(parent = globalenv()))
}
