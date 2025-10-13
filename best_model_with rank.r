## Extract Best Machine Learning Models from Results
## Analyzes all trained models and identifies the best performer for each protein
## Generates publication-ready summary tables and detailed model parameters

# Configuration
cwd <- getwd()

# Required packages
required_packages <- c("caret", "dplyr", "tidyr")
new.packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required_packages, library, character.only = TRUE)

# Protein list (same as training script)
prot_list <- c("OspC", "OspA", "DbpA", "BBK32", "RevA", "P66", "BB_0406")

# Feature subsets
feature_subsets <- c("fs055", "ImpQ1", "ImpQ2", "ImpQ3", "VarQ1", "VarQ2", "VarQ3")

# Model types
models <- c('PLS', "RF", "SVM", 'GLM', "PCANN")

# Validation methods
val_methods <- c("LOOfit", "bsfit")

# Initialize results storage
all_results <- list()
best_models_summary <- data.frame()

# Function to extract performance metrics from a single model
extract_model_performance <- function(model_obj, model_name, protein, feature_set, val_method) {
  if (is.null(model_obj)) return(NULL)
  
  # Extract test performance metrics
  test_perf <- model_obj$Test.Perf
  train_perf <- model_obj$Train.Perf
  
  # Helper function to safely extract column
  safe_extract <- function(matrix, col_name) {
    if (col_name %in% colnames(matrix)) {
      return(matrix[1, col_name])
    } else {
      return(NA)
    }
  }
  
  # Calculate balanced accuracy manually from Sens and Spec
  sens_val <- safe_extract(test_perf, "Sens")
  spec_val <- safe_extract(test_perf, "Spec")
  balanced_acc <- if (!is.na(sens_val) && !is.na(spec_val)) {
    (sens_val + spec_val) / 2
  } else {
    NA
  }
  
  # Create summary data frame
  perf_df <- data.frame(
    Protein = protein,
    FeatureSet = feature_set,
    Model = model_name,
    ValidationMethod = val_method,
    # Test metrics
    Test_Accuracy = safe_extract(test_perf, "Accuracy"),
    Test_Kappa = safe_extract(test_perf, "Kappa"),
    Test_ROC = safe_extract(test_perf, "ROC"),
    Test_Sens = sens_val,
    Test_Spec = spec_val,
    Test_BalancedAccuracy = balanced_acc,
    Test_Precision = safe_extract(test_perf, "Precision"),
    Test_Recall = safe_extract(test_perf, "Recall"),
    Test_F = safe_extract(test_perf, "F"),
    Test_AUC = safe_extract(test_perf, "AUC"),
    Test_PRAUC = safe_extract(test_perf, "PRAUC"),
    Test_Dist = safe_extract(test_perf, "Dist"),
    # Training metrics
    Train_ROC = if (!is.null(train_perf$TrainROC)) train_perf$TrainROC else NA,
    Train_Sens = if (!is.null(train_perf$TrainSens)) train_perf$TrainSens else NA,
    Train_Spec = if (!is.null(train_perf$TrainSpec)) train_perf$TrainSpec  else NA,
    Train_Accuracy = if (!is.null(train_perf$TrainAccuracy)) train_perf$TrainAccuracy else NA,
    stringsAsFactors = FALSE
  )
  
  return(perf_df)
}

# Function to extract model parameters
extract_model_params <- function(fit_obj) {
  if (is.null(fit_obj)) return(NULL)
  
  best_params <- fit_obj$bestTune
  final_model <- fit_obj$finalModel
  
  # Extract method-specific parameters
  params_list <- list(
    method = fit_obj$method,
    bestTune = best_params,
    tuneLength = length(fit_obj$results[,1]),
    finalModel_class = class(final_model)[1]
  )
  
  return(params_list)
}

# Main extraction loop
cat("=== Extracting Model Results ===\n\n")

for (prot in prot_list) {
  cat("Processing protein:", prot, "\n")
  
  results_path <- file.path(cwd, "results", prot)
  
  if (!dir.exists(results_path)) {
    warning(paste("Results directory not found for protein:", prot))
    next
  }
  
  # Loop through each feature subset
  for (ss in feature_subsets) {
    result_file <- file.path(results_path, paste0(ss, "_results.RData"))
    
    if (!file.exists(result_file)) {
      warning(paste("Result file not found:", result_file))
      next
    }
    
    # Load results
    cat("  Loading:", ss, "\n")
    env <- new.env()
    load(result_file, envir = env)
    
    # Extract all models
    for (model in models) {
      for (val_method in val_methods) {
        obj_name <- paste0(ss, "_", model, "_O_", val_method)
        
        if (exists(obj_name, envir = env)) {
          model_list <- get(obj_name, envir = env)
          
          # Process each fold
          for (fold_idx in seq_along(model_list)) {
            fold_result <- model_list[[fold_idx]]
            
            if (!is.null(fold_result)) {
              perf_df <- extract_model_performance(
                fold_result, model, prot, ss, 
                ifelse(val_method == "LOOfit", "LOOCV", "Bootstrap")
              )
              
              if (!is.null(perf_df)) {
                perf_df$Fold <- fold_idx
                best_models_summary <- rbind(best_models_summary, perf_df)
              }
            }
          }
        }
      }
    }
  }
}

# Calculate average performance across folds
cat("\n=== Calculating Average Performance ===\n")

avg_performance <- best_models_summary %>%
  group_by(Protein, FeatureSet, Model, ValidationMethod) %>%
  summarise(
    N_Folds = n(),
    Avg_Test_Accuracy = mean(Test_Accuracy, na.rm = TRUE),
    SD_Test_Accuracy = sd(Test_Accuracy, na.rm = TRUE),
    Avg_Test_ROC = mean(Test_ROC, na.rm = TRUE),
    SD_Test_ROC = sd(Test_ROC, na.rm = TRUE),
    Avg_Test_Sens = mean(Test_Sens, na.rm = TRUE),
    SD_Test_Sens = sd(Test_Sens, na.rm = TRUE),
    Avg_Test_Spec = mean(Test_Spec, na.rm = TRUE),
    SD_Test_Spec = sd(Test_Spec, na.rm = TRUE),
    Avg_Test_BalancedAccuracy = mean(Test_BalancedAccuracy, na.rm = TRUE),
    SD_Test_BalancedAccuracy = sd(Test_BalancedAccuracy, na.rm = TRUE),
    Avg_Test_F = mean(Test_F, na.rm = TRUE),
    SD_Test_F = sd(Test_F, na.rm = TRUE),
    Avg_Test_AUC = mean(Test_AUC, na.rm = TRUE),
    SD_Test_AUC = sd(Test_AUC, na.rm = TRUE),
    Avg_Test_PRAUC = mean(Test_PRAUC, na.rm = TRUE),
    SD_Test_PRAUC = sd(Test_PRAUC, na.rm = TRUE),
    
    Avg_Train_Accuracy = mean(Train_Accuracy, na.rm = TRUE),
    SD_Train_Accuracy = sd(Train_Accuracy, na.rm = TRUE),
    Avg_Train_ROC = mean(Train_ROC, na.rm = TRUE),
    SD_Train_ROC = sd(Train_ROC, na.rm = TRUE),
    Avg_Train_Sens = mean(Train_Sens, na.rm = TRUE),
    SD_Train_Sens = sd(Train_Sens, na.rm = TRUE),
    Avg_Train_Spec = mean(Train_Spec, na.rm = TRUE),
    SD_Train_Spec = sd(Train_Spec, na.rm = TRUE),
 
    .groups = 'drop'
  )

# Identify best model for each protein and validation method
cat("\n=== Identifying Best Models ===\n")

# Calculate composite score based on ROC, Sensitivity, and Specificity
avg_performance <- avg_performance %>%
  mutate(
    # Composite score: equal weight to ROC, Sens, and Spec
    Composite_Score = (Avg_Test_ROC + Avg_Test_Sens + Avg_Test_Spec) / 3
  )

best_models <- avg_performance %>%
  group_by(Protein, ValidationMethod) %>%
  arrange(desc(Composite_Score), desc(Avg_Test_ROC)) %>%
  slice(1) %>%
  ungroup()


# Extract detailed parameters for best models
best_model_details <- list()

for (i in 1:nrow(best_models)) {
  prot <- best_models$Protein[i]
  ss <- best_models$FeatureSet[i]
  model <- best_models$Model[i]
  val_method <- ifelse(best_models$ValidationMethod[i] == "LOOCV", "LOOfit", "bsfit")
  
  result_file <- file.path(cwd, "results", prot, paste0(ss, "_results.RData"))
  
  cat("Loading:", prot, "-", ss, "-", model, "-", val_method, "\n")
  
  if (file.exists(result_file)) {
    env <- new.env()
    load(result_file, envir = env)
    
    obj_name <- paste0(ss, "_", model, "_O_", val_method)
    cat("  Looking for object:", obj_name, "\n")
    
    if (exists(obj_name, envir = env)) {
      model_list <- get(obj_name, envir = env)
      cat("  Found object with", length(model_list), "folds\n")
      
      # Get first non-null fold's model for parameter details
      fit_obj <- NULL
      for (fold_idx in seq_along(model_list)) {
        if (!is.null(model_list[[fold_idx]]) && !is.null(model_list[[fold_idx]]$fit)) {
          fit_obj <- model_list[[fold_idx]]$fit
          cat("  Using fold", fold_idx, "for parameters\n")
          break
        }
      }
      
      if (!is.null(fit_obj)) {
        params <- extract_model_params(fit_obj)
        
        # Create unique key for the list
        detail_key <- paste(prot, best_models$ValidationMethod[i], sep="_")
        
        best_model_details[[detail_key]] <- list(
          protein = prot,
          feature_set = ss,
          model_type = model,
          validation = best_models$ValidationMethod[i],
          parameters = params,
          performance = best_models[i, ]
        )
        cat("  Successfully extracted parameters\n")
      } else {
        cat("  WARNING: No valid fit object found in any fold\n")
      }
    } else {
      cat("  WARNING: Object", obj_name, "not found in", result_file, "\n")
      cat("  Available objects:", paste(ls(envir = env), collapse=", "), "\n")
    }
  } else {
    cat("  WARNING: File not found:", result_file, "\n")
  }
}


# Create publication-ready table
cat("\n=== Creating Publication Table ===\n")
# First, calculate ranks for each scoring system within each Protein-ValidationMethod group
avg_performance_with_ranks <- avg_performance %>%
  group_by(Protein, ValidationMethod) %>%
  mutate(
    Rank = rank(desc(Composite_Score), ties.method = "min"),
    TestROC_formatted = sprintf("%.3f ± %.3f", Avg_Test_ROC, SD_Test_ROC),
    TestAcc_formatted = sprintf("%.3f ± %.3f", Avg_Test_Accuracy, SD_Test_Accuracy),
    TestSens_formatted = sprintf("%.3f ± %.3f", Avg_Test_Sens, SD_Test_Sens),
    TestSpec_formatted = sprintf("%.3f ± %.3f", Avg_Test_Spec, SD_Test_Spec),
    TrainROC_formatted = sprintf("%.3f ± %.3f", Avg_Train_ROC, SD_Train_ROC),
    TrainAcc_formatted = sprintf("%.3f ± %.3f", Avg_Train_Accuracy, SD_Train_Accuracy),
    TrainSens_formatted = sprintf("%.3f ± %.3f", Avg_Train_Sens, SD_Train_Sens),
    TrainSpec_formatted = sprintf("%.3f ± %.3f", Avg_Train_Spec, SD_Train_Spec),
  ) %>%
  ungroup()

# Create a lookup table with all the extra columns we need
rank_lookup <- avg_performance_with_ranks %>%
  select(Protein, FeatureSet, Model, ValidationMethod, 
         Rank, 
         Composite_Score)

# Join the ranks back to best_models
publication_table <- best_models %>%
  left_join(rank_lookup, by = c("Protein", "FeatureSet", "Model", "ValidationMethod")) %>%
  mutate(
    # Composite_Score = round(Composite_Score, 3),
    ROC_formatted = sprintf("%.3f ± %.3f", Avg_Test_ROC, SD_Test_ROC),
    Sens_formatted = sprintf("%.3f ± %.3f", Avg_Test_Sens, SD_Test_Sens),
    Spec_formatted = sprintf("%.3f ± %.3f", Avg_Test_Spec, SD_Test_Spec),
    Accuracy_formatted = sprintf("%.3f ± %.3f", Avg_Test_Accuracy, SD_Test_Accuracy),
    BalancedAcc_formatted = sprintf("%.3f ± %.3f", Avg_Test_BalancedAccuracy, SD_Test_BalancedAccuracy),
    F_Score_formatted = sprintf("%.3f ± %.3f", Avg_Test_F, SD_Test_F)
  )


# Add model parameters to publication table
publication_table$Model_Parameters <- NA_character_

for (i in 1:nrow(publication_table)) {
  prot <- publication_table$Protein[i]
  val_method <- publication_table$ValidationMethod[i]
  detail_key <- paste(prot, val_method, sep="_")
  
  if (detail_key %in% names(best_model_details)) {
    params <- best_model_details[[detail_key]]$parameters$bestTune
    if (!is.null(params) && nrow(params) > 0) {
      # Format parameters as "param1=value1; param2=value2"
      param_str <- paste(names(params), params[1,], sep="=", collapse="; ")
      publication_table$Model_Parameters[i] <- param_str
    }
  }
}

# Create a clean version for CSV with formatted values
publication_table_csv <- publication_table %>%
  select(
    Protein,
    Model,
    FeatureSet,
    ValidationMethod,
    N_Folds,
    Model_Parameters,
    ROC_formatted,
    Sens_formatted,
    Spec_formatted,
    BalancedAcc_formatted,
    F_Score_formatted
  )

# Save results
output_dir <- file.path(cwd, "results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save all results
save(
  best_models_summary,
  avg_performance,
  best_models,
  best_model_details,
  publication_table,
  file = file.path(output_dir, "best_models_summary.RData")
)

# Save CSV files
write.csv(avg_performance_with_ranks, 
          file.path(output_dir, "all_models_average_performance.csv"), 
          row.names = FALSE)

write.csv(publication_table_csv, 
          file.path(output_dir, "best_models_publication_table.csv"), 
          row.names = FALSE)

# Print summary
cat("\n=== BEST MODELS SUMMARY ===\n\n")
print(publication_table[, c("Protein", "Model", "FeatureSet", "ValidationMethod", 
                            "ROC_formatted", "Sens_formatted", "Spec_formatted", 
                            "Accuracy_formatted")])

# Print detailed parameters for each best model
cat("\n\n=== DETAILED MODEL PARAMETERS ===\n\n")
for (prot in names(best_model_details)) {
  cat("Protein:", prot, "\n")
  cat("Model:", best_model_details[[prot]]$model_type, "\n")
  cat("Feature Set:", best_model_details[[prot]]$feature_set, "\n")
  cat("Validation:", best_model_details[[prot]]$validation, "\n")
  cat("Best Hyperparameters:\n")
  print(best_model_details[[prot]]$parameters$bestTune)
  cat("\n")
}

cat("\nResults saved to:", output_dir, "\n")
cat("Files created:\n")
cat("  - best_models_summary.RData (complete R objects)\n")
cat("  - all_models_average_performance.csv (all model comparisons)\n")
cat("  - best_models_publication_table.csv (publication-ready table)\n")

cat("\n=== Extraction Complete ===\n")