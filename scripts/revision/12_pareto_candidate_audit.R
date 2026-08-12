suppressPackageStartupMessages(library(dplyr))
source(file.path("scripts", "config.R"))
output_dir <- revision_output_dir
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

all_lnp <- read.csv(file.path(primary_dir,"06_lnp_multiobjective_summary.csv"),stringsAsFactors=FALSE)
ranking <- read.csv(file.path(primary_dir,"04_lnp_ranking_by_safety_score.csv"),stringsAsFactors=FALSE)

audited <- all_lnp %>% left_join(ranking %>% select(lnp_id,n_simulations,sd_safety_score),by="lnp_id")
pareto <- audited %>% filter(pareto_optimal) %>% arrange(desc(mean_safety_score)) %>% mutate(pareto_sbd_rank=row_number())

stopifnot(nrow(pareto)==21, n_distinct(pareto$lnp_id)==21, all(complete.cases(pareto)))

publication_table <- pareto %>% transmute(
  Rank=pareto_sbd_rank, Formulation=lnp_id, `Size (nm)`=round(size_nm,2), `Charge (mV)`=round(charge_mV,2),
  `PEG (mol%)`=round(peg_mol_percent,3), `Ionizable lipid fraction`=round(ionizable_lipid_fraction,3),
  `Cholesterol fraction`=round(cholesterol_fraction,3), `Helper lipid fraction`=round(helper_lipid_fraction,3),
  `Targeting ligand`=targeting_ligand, Delivery=round(mean_delivery,3), Adaptive=round(mean_adaptive_activation,3),
  Innate=round(mean_innate_activation,3), Cytokine=round(mean_cytokine_burden,3), `Off-target`=round(mean_off_target_activation,3),
  `SbD score`=round(mean_safety_score,3), `SbD SD`=round(sd_safety_score,3), `Records`=n_simulations
)

parameters <- c("size_nm","charge_mV","peg_mol_percent","ionizable_lipid_fraction","cholesterol_fraction","helper_lipid_fraction","targeting_ligand")
endpoints <- c("mean_safety_score","mean_delivery","mean_adaptive_activation","mean_innate_activation","mean_cytokine_burden","mean_off_target_activation")
comparison <- bind_rows(lapply(c(parameters,endpoints),function(v){
  data.frame(variable=v,
    pareto_mean=mean(all_lnp[[v]][all_lnp$pareto_optimal]), pareto_median=median(all_lnp[[v]][all_lnp$pareto_optimal]),
    pareto_min=min(all_lnp[[v]][all_lnp$pareto_optimal]),pareto_max=max(all_lnp[[v]][all_lnp$pareto_optimal]),
    nonpareto_mean=mean(all_lnp[[v]][!all_lnp$pareto_optimal]),nonpareto_median=median(all_lnp[[v]][!all_lnp$pareto_optimal]),
    standardized_mean_difference=(mean(all_lnp[[v]][all_lnp$pareto_optimal])-mean(all_lnp[[v]][!all_lnp$pareto_optimal]))/sd(all_lnp[[v]])
  )
}))

# The composite score ranks the non-dominated candidates but is not a Pareto objective.
vals <- as.matrix(all_lnp[,c("mean_delivery","mean_adaptive_activation","mean_innate_activation","mean_cytokine_burden","mean_off_target_activation")])
maximize <- c(TRUE,TRUE,FALSE,FALSE,FALSE)
dominated <- logical(nrow(vals))
for(i in seq_len(nrow(vals))) for(j in seq_len(nrow(vals))) if(i!=j){
  be <- ifelse(maximize,vals[j,]>=vals[i,],vals[j,]<=vals[i,]); st <- ifelse(maximize,vals[j,]>vals[i,],vals[j,]<vals[i,])
  if(all(be)&&any(st)){dominated[i]<-TRUE;break}
}
verification <- data.frame(lnp_id=all_lnp$lnp_id,stored_pareto=all_lnp$pareto_optimal,recomputed_pareto=!dominated,match=all_lnp$pareto_optimal==!dominated)

summary <- data.frame(metric=c("total_formulations","pareto_formulations","all_unique","all_complete","pareto_flag_matches_recomputation","targeted_fraction_pareto","targeted_fraction_nonpareto","records_min","records_median","records_max"),value=c(nrow(all_lnp),nrow(pareto),n_distinct(pareto$lnp_id)==nrow(pareto),all(complete.cases(pareto)),all(verification$match),mean(pareto$targeting_ligand),mean(all_lnp$targeting_ligand[!all_lnp$pareto_optimal]),min(pareto$n_simulations),median(pareto$n_simulations),max(pareto$n_simulations)))

write.csv(pareto,file.path(output_dir,"12_pareto_21_candidates_full.csv"),row.names=FALSE)
write.csv(publication_table,file.path(output_dir,"12_pareto_21_candidates_publication_table.csv"),row.names=FALSE)
write.csv(comparison,file.path(output_dir,"12_pareto_vs_nonpareto_comparison.csv"),row.names=FALSE)
write.csv(verification,file.path(output_dir,"12_pareto_flag_verification.csv"),row.names=FALSE)
write.csv(summary,file.path(output_dir,"12_pareto_candidate_audit_summary.csv"),row.names=FALSE)
message("Pareto audit completed: 21 candidates verified.")
