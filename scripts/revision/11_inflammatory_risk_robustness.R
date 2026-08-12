suppressPackageStartupMessages(library(dplyr))

source(file.path("scripts", "config.R"))
output_dir <- revision_output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(primary_dir, "04_patient_lnp_safety_by_design_scores.csv"), stringsAsFactors = FALSE)
d$inflammatory_risk_class <- factor(d$inflammatory_risk_class, levels = c("Low", "Intermediate", "High"))

endpoints <- c("delivery_efficiency", "adaptive_activation", "innate_activation", "cytokine_burden", "off_target_activation", "safety_by_design_score")

cohort_summary <- d %>% group_by(inflammatory_risk_class) %>% summarise(
  n_patients = n_distinct(patient_id), n_records = n(), n_formulations = n_distinct(lnp_id),
  median_records_per_formulation = median(as.numeric(table(lnp_id))),
  min_records_per_formulation = min(as.numeric(table(lnp_id))),
  max_records_per_formulation = max(as.numeric(table(lnp_id))),
  .groups = "drop"
)

endpoint_summary <- d %>% group_by(inflammatory_risk_class) %>% summarise(
  across(all_of(endpoints), list(mean = ~mean(.x, na.rm=TRUE), sd = ~sd(.x, na.rm=TRUE), q025 = ~quantile(.x,.025,na.rm=TRUE), median = ~median(.x,na.rm=TRUE), q975 = ~quantile(.x,.975,na.rm=TRUE))),
  .groups = "drop"
)

formulation_risk <- d %>% group_by(inflammatory_risk_class, lnp_id, size_nm, charge_mV, peg_mol_percent, ionizable_lipid_fraction, cholesterol_fraction, helper_lipid_fraction, targeting_ligand) %>% summarise(
  n_records = n(), n_patients = n_distinct(patient_id),
  across(all_of(endpoints), \(x) mean(x, na.rm = TRUE)), .groups="drop"
)

rank_table <- formulation_risk %>% group_by(inflammatory_risk_class) %>% mutate(sbd_rank = rank(-safety_by_design_score, ties.method="min")) %>% ungroup()
rank_wide <- rank_table %>% select(inflammatory_risk_class, lnp_id, sbd_rank, safety_by_design_score, n_records, n_patients) %>% tidyr::pivot_wider(names_from=inflammatory_risk_class, values_from=c(sbd_rank,safety_by_design_score,n_records,n_patients))

complete_ranks <- rank_wide %>% filter(!is.na(sbd_rank_Low), !is.na(sbd_rank_Intermediate), !is.na(sbd_rank_High))
rank_concordance <- data.frame(
  comparison=c("Low_vs_Intermediate","Low_vs_High","Intermediate_vs_High"),
  n_common=c(sum(complete.cases(rank_wide$sbd_rank_Low,rank_wide$sbd_rank_Intermediate)),sum(complete.cases(rank_wide$sbd_rank_Low,rank_wide$sbd_rank_High)),sum(complete.cases(rank_wide$sbd_rank_Intermediate,rank_wide$sbd_rank_High))),
  spearman=c(cor(rank_wide$sbd_rank_Low,rank_wide$sbd_rank_Intermediate,use="complete.obs",method="spearman"),cor(rank_wide$sbd_rank_Low,rank_wide$sbd_rank_High,use="complete.obs",method="spearman"),cor(rank_wide$sbd_rank_Intermediate,rank_wide$sbd_rank_High,use="complete.obs",method="spearman"))
)

top20 <- rank_table %>% filter(sbd_rank <= 20) %>% select(inflammatory_risk_class, lnp_id, sbd_rank)
top20_overlap <- expand.grid(group1=levels(d$inflammatory_risk_class),group2=levels(d$inflammatory_risk_class),stringsAsFactors=FALSE) %>% filter(group1 < group2) %>% rowwise() %>% mutate(overlap=length(intersect(top20$lnp_id[top20$inflammatory_risk_class==group1],top20$lnp_id[top20$inflammatory_risk_class==group2]))) %>% ungroup()

# Separate Pareto fronts. Results are descriptive because subgroup replication is sparse.
pareto_one <- function(x) {
  vals <- as.matrix(x[,c("delivery_efficiency","adaptive_activation","innate_activation","cytokine_burden","off_target_activation")])
  maximize <- c(TRUE,TRUE,FALSE,FALSE,FALSE)
  dominated <- logical(nrow(vals))
  for (i in seq_len(nrow(vals))) {
    for (j in seq_len(nrow(vals))) if (i != j) {
      better_equal <- ifelse(maximize, vals[j,] >= vals[i,], vals[j,] <= vals[i,])
      strictly <- ifelse(maximize, vals[j,] > vals[i,], vals[j,] < vals[i,])
      if (all(better_equal) && any(strictly)) { dominated[i] <- TRUE; break }
    }
  }
  x$pareto_optimal <- !dominated
  x
}
pareto_by_risk <- bind_rows(lapply(split(formulation_risk, formulation_risk$inflammatory_risk_class), pareto_one))
pareto_summary <- pareto_by_risk %>% group_by(inflammatory_risk_class) %>% summarise(n_evaluated=n(), n_pareto=sum(pareto_optimal), median_records=median(n_records), min_records=min(n_records), .groups="drop")

write.csv(cohort_summary,file.path(output_dir,"11_risk_cohort_summary.csv"),row.names=FALSE)
write.csv(endpoint_summary,file.path(output_dir,"11_risk_endpoint_summary.csv"),row.names=FALSE)
write.csv(rank_table,file.path(output_dir,"11_formulation_performance_by_risk.csv"),row.names=FALSE)
write.csv(rank_wide,file.path(output_dir,"11_formulation_rank_robustness_by_risk.csv"),row.names=FALSE)
write.csv(rank_concordance,file.path(output_dir,"11_risk_rank_concordance.csv"),row.names=FALSE)
write.csv(top20_overlap,file.path(output_dir,"11_risk_top20_overlap.csv"),row.names=FALSE)
write.csv(pareto_by_risk,file.path(output_dir,"11_pareto_candidates_by_risk.csv"),row.names=FALSE)
write.csv(pareto_summary,file.path(output_dir,"11_pareto_by_risk_summary.csv"),row.names=FALSE)
message("Risk subgroup analysis completed.")
