suppressPackageStartupMessages(library(dplyr))

source(file.path("scripts", "config.R"))
out <- revision_output_dir
dir.create(out, recursive = TRUE, showWarnings = FALSE)

responses <- read.csv(file.path(primary_dir, "04_patient_lnp_safety_by_design_scores.csv"))
meta <- read.csv(file.path(input_dir, "03_multiomics_sample_metadata.csv"))
cyt <- read.csv(file.path(input_dir, "03_synthetic_cytokine_layer.csv"))
cells <- read.csv(file.path(input_dir, "03_synthetic_immune_cell_layer.csv"))

# Recover additive residuals before any pmax truncation. None of the realized
# cytokine/cell values are expected to be at zero, but this is checked below.
r_response <- responses %>% transmute(
  delivery = delivery_efficiency - (1.8 - 0.0008*(size_nm-90)^2 - 0.015*charge_mV^2 + 0.9*peg_mol_percent + 0.7*targeting_ligand + 0.5*ionizable_lipid_fraction),
  innate = innate_activation - (0.5 + 0.08*pmax(charge_mV,0) + 1.2*baseline_inflammation + 0.9*innate_reactivity - 0.6*peg_mol_percent),
  adaptive = adaptive_activation - (0.4 + 0.9*delivery_efficiency + 0.8*adaptive_potential + 0.3*targeting_ligand),
  cytokine = cytokine_burden - (0.3 + 0.8*innate_activation + 0.4*baseline_inflammation + 0.2*pmax(charge_mV,0)),
  off_target = off_target_activation - (0.4 + 0.5*abs(charge_mV) + 0.8*baseline_inflammation - 0.6*targeting_ligand)
)

d <- meta %>% select(sample_id, baseline_inflammation, delivery_efficiency, innate_activation, adaptive_activation, cytokine_burden) %>%
  left_join(cyt, by="sample_id") %>% left_join(cells, by="sample_id")
r_layers <- d %>% transmute(
  IL6=IL6-(5+4*cytokine_burden+2*baseline_inflammation),
  TNFa=TNFa-(4+3.5*innate_activation+1.5*baseline_inflammation),
  IL1b=IL1b-(3+3*innate_activation+2*cytokine_burden),
  IFNg=IFNg-(4+3*adaptive_activation+innate_activation),
  IFNa=IFNa-(3+3.5*innate_activation+1.5*cytokine_burden),
  CXCL10=CXCL10-(4+4*innate_activation+2*cytokine_burden),
  CD4=CD4_T_cells-(100+30*adaptive_activation),
  CD8=CD8_T_cells-(90+35*adaptive_activation+10*delivery_efficiency),
  B=B_cells-(80+25*adaptive_activation),
  Plasma=Plasma_cells-(20+18*adaptive_activation),
  Monocytes=Monocytes-(120+25*innate_activation+20*baseline_inflammation),
  NK=NK_cells-(70+20*innate_activation)
)

summ <- function(x, expected_sd, layer) data.frame(
  layer=layer, variable=names(x), expected_sd=expected_sd,
  empirical_mean=sapply(x, mean), empirical_sd=sapply(x, sd),
  empirical_variance=sapply(x, var), stringsAsFactors=FALSE
)
noise_summary <- bind_rows(
  summ(r_response, c(.10,.12,.12,.15,.15), "response_endpoint"),
  summ(r_layers[,1:6], rep(1,6), "cytokine_layer"),
  summ(r_layers[,7:12], c(8,8,8,5,10,7), "cell_layer")
)

cor_response <- cor(r_response)
cor_layers <- cor(r_layers)
write.csv(noise_summary, file.path(out,"13_noise_distribution_summary.csv"), row.names=FALSE)
write.csv(cbind(variable=rownames(cor_response), as.data.frame(cor_response)), file.path(out,"13_response_residual_correlations.csv"), row.names=FALSE)
write.csv(cbind(variable=rownames(cor_layers), as.data.frame(cor_layers)), file.path(out,"13_multiomics_residual_correlations.csv"), row.names=FALSE)

cross_pairs <- data.frame(
  cytokine=c("IL6","TNFa","IL1b","IFNg","IFNa","CXCL10"),
  cell=c("CD4","CD8","B","Plasma","Monocytes","NK")
) %>% rowwise() %>% mutate(residual_correlation=cor(r_layers[[cytokine]],r_layers[[cell]])) %>% ungroup()
write.csv(cross_pairs,file.path(out,"13_same_seed_cross_layer_pairs.csv"),row.names=FALSE)

checks <- data.frame(
  metric=c("max_abs_response_recovered_residual_correlation_after_clipping","max_abs_within_cytokine_residual_correlation","max_abs_within_cell_residual_correlation","max_abs_cross_layer_residual_correlation","zero_delivery","zero_innate","zero_adaptive","zero_cytokine_burden","zero_off_target","zero_cytokine_layer_values","zero_cell_layer_values"),
  value=c(max(abs(cor_response[upper.tri(cor_response)])), max(abs(cor_layers[1:6,1:6][upper.tri(cor_layers[1:6,1:6])])), max(abs(cor_layers[7:12,7:12][upper.tri(cor_layers[7:12,7:12])])), max(abs(cor_layers[1:6,7:12])), sum(responses$delivery_efficiency==0),sum(responses$innate_activation==0),sum(responses$adaptive_activation==0),sum(responses$cytokine_burden==0),sum(responses$off_target_activation==0), sum(cyt[,-1]==0), sum(cells[,-1]==0))
)
write.csv(checks,file.path(out,"13_noise_audit_checks.csv"),row.names=FALSE)
message("Noise covariance audit completed.")
