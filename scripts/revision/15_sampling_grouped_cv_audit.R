suppressPackageStartupMessages(library(dplyr))

source(file.path("scripts", "config.R"))
out <- revision_output_dir
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(file.path(primary_dir,"04_patient_lnp_safety_by_design_scores.csv"))
folds <- read.csv(file.path(input_dir,"05_lnp_grouped_fold_assignments.csv"))

pair_counts <- d %>% count(patient_id,lnp_id,name="n_records")
form_counts <- d %>% count(lnp_id,name="n_records")
patient_counts <- d %>% count(patient_id,name="n_records")

expected_unique <- 200*500*(1-(1-1/(200*500))^nrow(d))
summary <- data.frame(
  metric=c("records","possible_patient_formulation_pairs","unique_realized_pairs","design_space_coverage_fraction","duplicate_records_beyond_first","pairs_observed_more_than_once","maximum_repeats_of_one_pair","expected_unique_pairs_under_independent_sampling","unique_patients","unique_formulations","records_per_formulation_min","records_per_formulation_median","records_per_formulation_max","records_per_patient_min","records_per_patient_median","records_per_patient_max"),
  value=c(nrow(d),200*500,nrow(pair_counts),nrow(pair_counts)/(200*500),nrow(d)-nrow(pair_counts),sum(pair_counts$n_records>1),max(pair_counts$n_records),expected_unique,n_distinct(d$patient_id),n_distinct(d$lnp_id),min(form_counts$n_records),median(form_counts$n_records),max(form_counts$n_records),min(patient_counts$n_records),median(patient_counts$n_records),max(patient_counts$n_records))
)

# Verify formulation-grouped folds and quantify patient overlap.
ff <- folds %>% distinct(lnp_id,fold)
form_fold_check <- folds %>% group_by(lnp_id) %>% summarise(n_folds=n_distinct(fold),.groups="drop")
fold_summary <- folds %>% group_by(fold) %>% summarise(records=n(),patients=n_distinct(patient_id),formulations=n_distinct(lnp_id),.groups="drop")
patient_fold_presence <- folds %>% distinct(patient_id,fold) %>% count(patient_id,name="folds_present")
cv_checks <- data.frame(
  metric=c("formulations_assigned_to_exactly_one_fold","minimum_formulations_per_fold","maximum_formulations_per_fold","minimum_records_per_fold","maximum_records_per_fold","minimum_patients_per_fold","maximum_patients_per_fold","patients_present_in_all_five_folds"),
  value=c(sum(form_fold_check$n_folds==1),min(fold_summary$formulations),max(fold_summary$formulations),min(fold_summary$records),max(fold_summary$records),min(fold_summary$patients),max(fold_summary$patients),sum(patient_fold_presence$folds_present==5))
)

write.csv(summary,file.path(out,"15_sampling_audit_summary.csv"),row.names=FALSE)
write.csv(pair_counts,file.path(out,"15_patient_formulation_pair_counts.csv"),row.names=FALSE)
write.csv(form_counts,file.path(out,"15_records_per_formulation.csv"),row.names=FALSE)
write.csv(patient_counts,file.path(out,"15_records_per_patient.csv"),row.names=FALSE)
write.csv(fold_summary,file.path(out,"15_grouped_cv_fold_summary.csv"),row.names=FALSE)
write.csv(cv_checks,file.path(out,"15_grouped_cv_checks.csv"),row.names=FALSE)
message("Sampling and grouped-CV audit completed.")
