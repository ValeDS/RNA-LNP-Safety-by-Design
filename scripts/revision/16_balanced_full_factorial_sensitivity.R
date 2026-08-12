suppressPackageStartupMessages(library(dplyr))

source(file.path("scripts", "config.R"))
out <- revision_output_dir
dir.create(out, recursive = TRUE, showWarnings = FALSE)
seed <- 20260812

patients <- read.csv(file.path(input_dir,"01_virtual_patient_cohort.csv"))
lnp <- read.csv(file.path(input_dir,"01_lnp_design_space.csv"))
primary <- read.csv(file.path(primary_dir,"04_patient_lnp_safety_by_design_scores.csv"))
norm <- read.csv(file.path(input_dir,"04_sbd_normalization_parameters.csv"))

# Complete Cartesian product: one realization for every patient-formulation pair.
d <- merge(patients,lnp,by=NULL,sort=FALSE) %>% arrange(patient_id,lnp_id) %>% mutate(pair_id=paste0("BAL_",sprintf("%06d",row_number())))
set.seed(seed)
d <- d %>% mutate(
  delivery_efficiency=1.8-0.0008*(size_nm-90)^2-0.015*charge_mV^2+0.9*peg_mol_percent+0.7*targeting_ligand+0.5*ionizable_lipid_fraction+rnorm(n(),0,.10),
  delivery_efficiency=pmax(delivery_efficiency,0),
  innate_activation=0.5+0.08*pmax(charge_mV,0)+1.2*baseline_inflammation+0.9*innate_reactivity-0.6*peg_mol_percent+rnorm(n(),0,.12),
  innate_activation=pmax(innate_activation,0),
  adaptive_activation=0.4+0.9*delivery_efficiency+0.8*adaptive_potential+0.3*targeting_ligand+rnorm(n(),0,.12),
  adaptive_activation=pmax(adaptive_activation,0),
  cytokine_burden=0.3+0.8*innate_activation+0.4*baseline_inflammation+0.2*pmax(charge_mV,0)+rnorm(n(),0,.15),
  cytokine_burden=pmax(cytokine_burden,0),
  off_target_activation=0.4+0.5*abs(charge_mV)+0.8*baseline_inflammation-0.6*targeting_ligand+rnorm(n(),0,.15),
  off_target_activation=pmax(off_target_activation,0)
)

getrow <- function(ep) norm[norm$endpoint==ep,,drop=FALSE]
scale_locked <- function(x,ep) { z <- (x-getrow(ep)$minimum)/getrow(ep)$range; pmin(1,pmax(0,z)) }
d <- d %>% mutate(
  delivery_norm=scale_locked(delivery_efficiency,"delivery_efficiency"), adaptive_norm=scale_locked(adaptive_activation,"adaptive_activation"),
  innate_norm=scale_locked(innate_activation,"innate_activation"), cytokine_norm=scale_locked(cytokine_burden,"cytokine_burden"), off_target_norm=scale_locked(off_target_activation,"off_target_activation"),
  safety_by_design_score=.30*delivery_norm+.20*adaptive_norm+.20*(1-innate_norm)+.15*(1-cytokine_norm)+.15*(1-off_target_norm)
)

params <- c("lnp_id","size_nm","charge_mV","peg_mol_percent","ionizable_lipid_fraction","cholesterol_fraction","helper_lipid_fraction","targeting_ligand")
summarize_lnp <- function(x,groups=character()) x %>% group_by(across(all_of(c(groups,params)))) %>% summarise(
  n_records=n(), mean_safety_score=mean(safety_by_design_score),sd_safety_score=sd(safety_by_design_score),mean_delivery=mean(delivery_efficiency),mean_adaptive_activation=mean(adaptive_activation),mean_innate_activation=mean(innate_activation),mean_cytokine_burden=mean(cytokine_burden),mean_off_target_activation=mean(off_target_activation),.groups="drop")
bal_global <- summarize_lnp(d) %>% arrange(desc(mean_safety_score)) %>% mutate(sbd_rank=row_number())
bal_risk <- summarize_lnp(d,"inflammatory_risk_class") %>% group_by(inflammatory_risk_class) %>% arrange(desc(mean_safety_score),.by_group=TRUE) %>% mutate(sbd_rank=row_number()) %>% ungroup()

objectives <- c("mean_delivery","mean_adaptive_activation","mean_innate_activation","mean_cytokine_burden","mean_off_target_activation")
pareto_flag <- function(x) {
  y <- as.matrix(x[,objectives]); n <- nrow(y); ans <- rep(TRUE,n)
  for (i in seq_len(n)) { others <- y[-i,,drop=FALSE]; cand <- y[i,]; weak <- others[,1]>=cand[1] & others[,2]>=cand[2] & others[,3]<=cand[3] & others[,4]<=cand[4] & others[,5]<=cand[5]; strict <- others[,1]>cand[1] | others[,2]>cand[2] | others[,3]<cand[3] | others[,4]<cand[4] | others[,5]<cand[5]; ans[i] <- !any(weak & strict) }
  ans
}
bal_global$pareto_optimal <- pareto_flag(bal_global)
bal_risk <- bal_risk %>% group_by(inflammatory_risk_class) %>% group_modify(~mutate(.x,pareto_optimal=pareto_flag(.x))) %>% ungroup()

primary_global <- primary %>% group_by(lnp_id) %>% summarise(mean_safety_score=mean(safety_by_design_score),.groups="drop") %>% arrange(desc(mean_safety_score)) %>% mutate(primary_rank=row_number())
rank_compare <- bal_global %>% select(lnp_id,balanced_rank=sbd_rank,balanced_score=mean_safety_score,balanced_pareto=pareto_optimal) %>% left_join(primary_global,by="lnp_id") %>% arrange(balanced_rank)
primary_risk <- primary %>% group_by(inflammatory_risk_class,lnp_id) %>% summarise(score=mean(safety_by_design_score),.groups="drop") %>% group_by(inflammatory_risk_class) %>% arrange(desc(score),.by_group=TRUE) %>% mutate(primary_rank=row_number()) %>% ungroup()
risk_rank_compare <- bal_risk %>% select(inflammatory_risk_class,lnp_id,balanced_rank=sbd_rank,balanced_score=mean_safety_score,balanced_pareto=pareto_optimal,n_records) %>% left_join(primary_risk,by=c("inflammatory_risk_class","lnp_id"))

risk_summary <- d %>% group_by(inflammatory_risk_class) %>% summarise(patients=n_distinct(patient_id),records=n(),formulations=n_distinct(lnp_id),records_per_formulation=records/formulations,delivery=mean(delivery_efficiency),adaptive=mean(adaptive_activation),innate=mean(innate_activation),cytokine=mean(cytokine_burden),off_target=mean(off_target_activation),sbd=mean(safety_by_design_score),.groups="drop")
pareto_counts <- bal_risk %>% count(inflammatory_risk_class,wt=pareto_optimal,name="pareto_formulations")

top_overlap <- function(a,b,k) length(intersect(a[1:k],b[1:k]))
global_metrics <- data.frame(metric=c("records","unique_pairs","global_rank_spearman","global_top10_overlap","global_top20_overlap","global_pareto_count","primary_top_candidate_balanced_rank","primary_top_candidate_balanced_score"),value=c(nrow(d),n_distinct(paste(d$patient_id,d$lnp_id)),cor(rank_compare$balanced_rank,rank_compare$primary_rank,method="spearman"),top_overlap(rank_compare$lnp_id,primary_global$lnp_id,10),top_overlap(rank_compare$lnp_id,primary_global$lnp_id,20),sum(bal_global$pareto_optimal),rank_compare$balanced_rank[rank_compare$lnp_id=="LNP_0462"],rank_compare$balanced_score[rank_compare$lnp_id=="LNP_0462"]))
risk_metrics <- risk_rank_compare %>% filter(!is.na(primary_rank)) %>% group_by(inflammatory_risk_class) %>% summarise(common_formulations=n(),rank_spearman=cor(balanced_rank,primary_rank,method="spearman"),top20_overlap=length(intersect(lnp_id[order(balanced_rank)][1:20],lnp_id[order(primary_rank)][1:20])),.groups="drop")

range_checks <- data.frame(endpoint=norm$endpoint,below_primary_min=sapply(norm$endpoint,function(ep) sum(d[[ep]]<getrow(ep)$minimum)),above_primary_max=sapply(norm$endpoint,function(ep) sum(d[[ep]]>getrow(ep)$maximum)))
metadata <- data.frame(item=c("analysis_role","records","patients","formulations","seed","normalization","noise_realizations_per_pair"),value=c("balanced full-factorial sensitivity analysis",nrow(d),n_distinct(d$patient_id),n_distinct(d$lnp_id),seed,"locked to primary 10,000-record min/max with clipping",1))

write.csv(d,file.path(out,"16_balanced_100000_scored_records.csv"),row.names=FALSE)
write.csv(bal_global,file.path(out,"16_balanced_global_formulation_ranking.csv"),row.names=FALSE)
write.csv(bal_risk,file.path(out,"16_balanced_risk_formulation_rankings_pareto.csv"),row.names=FALSE)
write.csv(risk_summary,file.path(out,"16_balanced_risk_endpoint_summary.csv"),row.names=FALSE)
write.csv(pareto_counts,file.path(out,"16_balanced_risk_pareto_counts.csv"),row.names=FALSE)
write.csv(rank_compare,file.path(out,"16_balanced_vs_primary_global_ranking.csv"),row.names=FALSE)
write.csv(risk_rank_compare,file.path(out,"16_balanced_vs_primary_risk_rankings.csv"),row.names=FALSE)
write.csv(global_metrics,file.path(out,"16_balanced_global_robustness_metrics.csv"),row.names=FALSE)
write.csv(risk_metrics,file.path(out,"16_balanced_risk_robustness_metrics.csv"),row.names=FALSE)
write.csv(range_checks,file.path(out,"16_balanced_primary_normalization_range_checks.csv"),row.names=FALSE)
write.csv(metadata,file.path(out,"16_balanced_analysis_metadata.csv"),row.names=FALSE)
message("Balanced full-factorial sensitivity analysis completed.")
