# 05_train_random_forest_grouped_cv.R
# Reproducible Random Forest analysis of the Safety-by-Design score
#
# Primary validation:
#   5-fold cross-validation grouped by RNA-LNP formulation (lnp_id).
#   All records from a formulation are assigned to the same fold.
#
# Sensitivity analysis:
#   5-fold cross-validation grouped by virtual patient (patient_id).
#
# The model uses 11 predictors:
#   - four latent immune-state variables;
#   - seven RNA-LNP formulation variables.
#
# Age, sex, UISS-derived descriptors, MIIS, and synthetic multi-omics
# variables are intentionally excluded.

suppressPackageStartupMessages({
  library(dplyr)
  library(ranger)
})

seed <- 123
n_folds <- 5
n_trees <- 300

source(file.path("scripts", "config.R"))
output_dir <- Sys.getenv(
  "RNA_LNP_OUTPUT_DIR",
  unset = file.path(project_dir, "results", "core")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_candidates <- c(
  file.path(core_output_dir, "04_patient_lnp_safety_by_design_scores.csv"),
  file.path(primary_dir, "04_patient_lnp_safety_by_design_scores.csv")
)

existing_inputs <- input_candidates[file.exists(input_candidates)]

if (length(existing_inputs) == 0) {
  stop(
    paste0(
      "Input file not found. Expected one of:\n- ",
      paste(input_candidates, collapse = "\n- ")
    )
  )
}

input_file <- existing_inputs[[1]]
data <- read.csv(input_file, stringsAsFactors = FALSE)

predictors <- c(
  "baseline_inflammation",
  "immune_responsiveness",
  "innate_reactivity",
  "adaptive_potential",
  "size_nm",
  "charge_mV",
  "peg_mol_percent",
  "ionizable_lipid_fraction",
  "cholesterol_fraction",
  "helper_lipid_fraction",
  "targeting_ligand"
)

required_columns <- c(
  "pair_id",
  "patient_id",
  "lnp_id",
  predictors,
  "safety_by_design_score"
)

missing_columns <- setdiff(required_columns, names(data))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

model_data <- data %>%
  select(all_of(required_columns)) %>%
  filter(if_all(all_of(c(predictors, "safety_by_design_score")), ~ !is.na(.x)))

if (nrow(model_data) == 0) {
  stop("No complete records are available for model training.")
}

make_grouped_folds <- function(data, group_column, k = 5, seed = 123) {
  set.seed(seed)

  group_ids <- sample(unique(data[[group_column]]))
  group_to_fold <- setNames(
    rep(seq_len(k), length.out = length(group_ids)),
    group_ids
  )

  fold_id <- unname(group_to_fold[as.character(data[[group_column]])])

  data.frame(
    pair_id = data$pair_id,
    patient_id = data$patient_id,
    lnp_id = data$lnp_id,
    validation_group = group_column,
    fold = fold_id,
    stringsAsFactors = FALSE
  )
}

metric_table <- function(observed, predicted) {
  residuals <- observed - predicted
  sse <- sum(residuals^2)
  sst <- sum((observed - mean(observed))^2)

  data.frame(
    n = length(observed),
    RMSE = sqrt(mean(residuals^2)),
    MAE = mean(abs(residuals)),
    R2 = 1 - sse / sst,
    correlation_R2 = cor(observed, predicted)^2
  )
}

run_grouped_cv <- function(
    data,
    predictors,
    group_column,
    k = 5,
    n_trees = 300,
    seed = 123
) {
  assignments <- make_grouped_folds(
    data = data,
    group_column = group_column,
    k = k,
    seed = seed
  )

  fold_id <- assignments$fold
  prediction_list <- vector("list", k)
  performance_list <- vector("list", k)

  for (fold in seq_len(k)) {
    train_data <- data[fold_id != fold, , drop = FALSE]
    test_data <- data[fold_id == fold, , drop = FALSE]

    formula <- reformulate(
      termlabels = predictors,
      response = "safety_by_design_score"
    )

    model <- ranger(
      formula = formula,
      data = train_data,
      num.trees = n_trees,
      importance = "permutation",
      seed = seed + fold,
      num.threads = 1
    )

    predicted <- predict(model, data = test_data)$predictions

    prediction_list[[fold]] <- data.frame(
      pair_id = test_data$pair_id,
      patient_id = test_data$patient_id,
      lnp_id = test_data$lnp_id,
      fold = fold,
      observed = test_data$safety_by_design_score,
      predicted = predicted
    )

    performance_list[[fold]] <- cbind(
      data.frame(fold = fold),
      metric_table(
        observed = test_data$safety_by_design_score,
        predicted = predicted
      )
    )
  }

  predictions <- bind_rows(prediction_list)
  performance_by_fold <- bind_rows(performance_list)

  pooled_performance <- cbind(
    data.frame(
      validation = paste0(k, "-fold grouped CV by ", group_column),
      aggregation = "pooled out-of-fold predictions"
    ),
    metric_table(
      observed = predictions$observed,
      predicted = predictions$predicted
    )
  )

  mean_performance <- performance_by_fold %>%
    summarise(
      validation = paste0(k, "-fold grouped CV by ", group_column),
      aggregation = "mean across folds",
      n = sum(n),
      RMSE = mean(RMSE),
      MAE = mean(MAE),
      R2 = mean(R2),
      correlation_R2 = mean(correlation_R2)
    )

  sd_performance <- performance_by_fold %>%
    summarise(
      validation = paste0(k, "-fold grouped CV by ", group_column),
      aggregation = "SD across folds",
      n = NA_integer_,
      RMSE = sd(RMSE),
      MAE = sd(MAE),
      R2 = sd(R2),
      correlation_R2 = sd(correlation_R2)
    )

  list(
    assignments = assignments,
    predictions = predictions,
    performance_by_fold = performance_by_fold,
    performance_summary = bind_rows(
      pooled_performance,
      mean_performance,
      sd_performance
    )
  )
}

# Primary formulation-grouped validation

lnp_cv <- run_grouped_cv(
  data = model_data,
  predictors = predictors,
  group_column = "lnp_id",
  k = n_folds,
  n_trees = n_trees,
  seed = seed
)

# Patient-grouped sensitivity analysis

patient_cv <- run_grouped_cv(
  data = model_data,
  predictors = predictors,
  group_column = "patient_id",
  k = n_folds,
  n_trees = n_trees,
  seed = seed
)

# Final model trained on all records for interpretation

final_formula <- reformulate(
  termlabels = predictors,
  response = "safety_by_design_score"
)

final_model <- ranger(
  formula = final_formula,
  data = model_data,
  num.trees = n_trees,
  importance = "permutation",
  seed = seed,
  num.threads = 1,
  keep.inbag = TRUE
)

importance_df <- data.frame(
  feature = names(final_model$variable.importance),
  importance = unname(final_model$variable.importance),
  row.names = NULL
) %>%
  arrange(desc(importance))

analysis_metadata <- data.frame(
  item = c(
    "input_file",
    "records",
    "unique_virtual_patients",
    "unique_lnp_formulations",
    "predictors",
    "trees",
    "primary_validation",
    "sensitivity_validation",
    "random_seed"
  ),
  value = c(
    normalizePath(input_file),
    nrow(model_data),
    n_distinct(model_data$patient_id),
    n_distinct(model_data$lnp_id),
    length(predictors),
    n_trees,
    paste0(n_folds, "-fold grouped CV by lnp_id"),
    paste0(n_folds, "-fold grouped CV by patient_id"),
    seed
  )
)

# Export

write.csv(
  lnp_cv$predictions,
  file.path(output_dir, "05_observed_vs_predicted.csv"),
  row.names = FALSE
)

write.csv(
  lnp_cv$performance_by_fold,
  file.path(output_dir, "05_model_performance_by_fold.csv"),
  row.names = FALSE
)

write.csv(
  lnp_cv$performance_summary,
  file.path(output_dir, "05_model_performance.csv"),
  row.names = FALSE
)

write.csv(
  lnp_cv$assignments,
  file.path(output_dir, "05_lnp_grouped_fold_assignments.csv"),
  row.names = FALSE
)

write.csv(
  patient_cv$predictions,
  file.path(
    output_dir,
    "05_patient_grouped_observed_vs_predicted.csv"
  ),
  row.names = FALSE
)

write.csv(
  patient_cv$performance_by_fold,
  file.path(
    output_dir,
    "05_patient_grouped_performance_by_fold.csv"
  ),
  row.names = FALSE
)

write.csv(
  patient_cv$performance_summary,
  file.path(
    output_dir,
    "05_patient_grouped_performance.csv"
  ),
  row.names = FALSE
)

write.csv(
  patient_cv$assignments,
  file.path(output_dir, "05_patient_grouped_fold_assignments.csv"),
  row.names = FALSE
)

write.csv(
  importance_df,
  file.path(output_dir, "05_feature_importance.csv"),
  row.names = FALSE
)

write.csv(
  analysis_metadata,
  file.path(output_dir, "05_analysis_metadata.csv"),
  row.names = FALSE
)

saveRDS(
  final_model,
  file.path(output_dir, "05_random_forest_safety_model.rds")
)

message("Random Forest grouped-cross-validation analysis completed.")
message("Input: ", normalizePath(input_file))
message("Output directory: ", normalizePath(output_dir))
message(
  "Primary pooled R2: ",
  round(lnp_cv$performance_summary$R2[[1]], 4),
  "; RMSE: ",
  round(lnp_cv$performance_summary$RMSE[[1]], 4)
)
message(
  "Patient-grouped pooled R2: ",
  round(patient_cv$performance_summary$R2[[1]], 4),
  "; RMSE: ",
  round(patient_cv$performance_summary$RMSE[[1]], 4)
)
