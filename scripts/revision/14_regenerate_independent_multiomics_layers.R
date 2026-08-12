suppressPackageStartupMessages(library(dplyr))

source(file.path("scripts", "config.R"))
out <- revision_output_dir
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- read.csv(file.path(input_dir, "03_multiomics_sample_metadata.csv"))

generate_cytokine_layer <- function(samples, seed = 124) {
  set.seed(seed)
  samples %>% transmute(
    sample_id,
    IL6=pmax(0,5+4*cytokine_burden+2*baseline_inflammation+rnorm(n(),0,1)),
    TNFa=pmax(0,4+3.5*innate_activation+1.5*baseline_inflammation+rnorm(n(),0,1)),
    IL1b=pmax(0,3+3*innate_activation+2*cytokine_burden+rnorm(n(),0,1)),
    IFNg=pmax(0,4+3*adaptive_activation+innate_activation+rnorm(n(),0,1)),
    IFNa=pmax(0,3+3.5*innate_activation+1.5*cytokine_burden+rnorm(n(),0,1)),
    CXCL10=pmax(0,4+4*innate_activation+2*cytokine_burden+rnorm(n(),0,1))
  )
}

generate_cell_layer <- function(samples, seed = 125) {
  set.seed(seed)
  samples %>% transmute(
    sample_id,
    CD4_T_cells=pmax(0,100+30*adaptive_activation+rnorm(n(),0,8)),
    CD8_T_cells=pmax(0,90+35*adaptive_activation+10*delivery_efficiency+rnorm(n(),0,8)),
    B_cells=pmax(0,80+25*adaptive_activation+rnorm(n(),0,8)),
    Plasma_cells=pmax(0,20+18*adaptive_activation+rnorm(n(),0,5)),
    Monocytes=pmax(0,120+25*innate_activation+20*baseline_inflammation+rnorm(n(),0,10)),
    NK_cells=pmax(0,70+20*innate_activation+rnorm(n(),0,7))
  )
}

cyt <- generate_cytokine_layer(samples)
cells <- generate_cell_layer(samples)
write.csv(cyt,file.path(out,"14_synthetic_cytokine_layer_independent.csv"),row.names=FALSE)
write.csv(cells,file.path(out,"14_synthetic_immune_cell_layer_independent.csv"),row.names=FALSE)

d <- samples %>% select(sample_id,baseline_inflammation,delivery_efficiency,innate_activation,adaptive_activation,cytokine_burden) %>% left_join(cyt,by="sample_id") %>% left_join(cells,by="sample_id")
r <- d %>% transmute(
  IL6=IL6-(5+4*cytokine_burden+2*baseline_inflammation), TNFa=TNFa-(4+3.5*innate_activation+1.5*baseline_inflammation), IL1b=IL1b-(3+3*innate_activation+2*cytokine_burden), IFNg=IFNg-(4+3*adaptive_activation+innate_activation), IFNa=IFNa-(3+3.5*innate_activation+1.5*cytokine_burden), CXCL10=CXCL10-(4+4*innate_activation+2*cytokine_burden),
  CD4=CD4_T_cells-(100+30*adaptive_activation), CD8=CD8_T_cells-(90+35*adaptive_activation+10*delivery_efficiency), B=B_cells-(80+25*adaptive_activation), Plasma=Plasma_cells-(20+18*adaptive_activation), Monocytes=Monocytes-(120+25*innate_activation+20*baseline_inflammation), NK=NK_cells-(70+20*innate_activation)
)
cormat <- cor(r)
write.csv(cbind(variable=rownames(cormat),as.data.frame(cormat)),file.path(out,"14_independent_layer_residual_correlations.csv"),row.names=FALSE)

summary <- data.frame(
  metric=c("cytokine_seed","cell_seed","max_abs_within_cytokine_correlation","max_abs_within_cell_correlation","max_abs_cross_layer_correlation","mean_abs_cross_layer_correlation","zero_cytokine_values","zero_cell_values"),
  value=c(124,125,max(abs(cormat[1:6,1:6][upper.tri(cormat[1:6,1:6])])),max(abs(cormat[7:12,7:12][upper.tri(cormat[7:12,7:12])])),max(abs(cormat[1:6,7:12])),mean(abs(cormat[1:6,7:12])),sum(cyt[,-1]==0),sum(cells[,-1]==0))
)
write.csv(summary,file.path(out,"14_independent_layer_audit_summary.csv"),row.names=FALSE)

dist <- data.frame(variable=names(r), empirical_mean=sapply(r,mean), empirical_sd=sapply(r,sd), empirical_variance=sapply(r,var))
write.csv(dist,file.path(out,"14_independent_layer_noise_distributions.csv"),row.names=FALSE)
message("Independent multi-omics layers regenerated and audited.")
