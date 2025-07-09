## Machine Learning Model Training with Caret 
## Trains multiple ML models (PLS, RF, SVM, GLM, PCANN) on one-hot encoding protein sequence data
## Uses parallel processing for efficiency

# Configuration
ncores <- 36  # Number of cores for parallel processing
cwd <- getwd()

### Requited Packages####
# Core ML packages
required_packages <- c(
  "caret",           # Main ML framework
  "pls",             # PLS method
  "randomForest",    # Random Forest method
  "kernlab",         # SVM method
  "glmnet",          # GLM method
  "nnet",            # Neural Network method
  "doParallel",      # Parallel processing
  "foreach"          # Parallel loops
)

# Install missing packages
new.packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
lapply(required_packages, library, character.only = TRUE)

print("Successfully loaded minimal required packages", quote = FALSE)

##########
# Custom Summary Function for Model Evaluation
MySummary <- function(data, lev = c("YES","NO"), model = NULL) {
  # Combine multiple performance metrics
  a1 <- defaultSummary(data, lev, model)
  b1 <- twoClassSummary(data, lev, model)
  c1 <- prSummary(data, lev, model)
  out <- c(a1, b1, c1)
  
  # Calculate balanced accuracy and distance from perfect performance
  # Perfect model has sensitivity = 1 and specificity = 1
  coords <- matrix(c(1, 1, out["Spec"], out["Sens"]), 
                   ncol = 2, 
                   byrow = TRUE)
  colnames(coords) <- c("Spec", "Sens")
  rownames(coords) <- c("Best", "Current")
  
  # Fixed balanced accuracy calculation
  BAcc <- (out["Spec"] + out["Sens"]) / 2
  
  c(out, BalancedAccuracy = BAcc, Dist = dist(coords)[1])
}

########
# Protein list for analysis
prot_list <- c("OspC", "OspA", "DbpA", "BBK32", "RevA", "P66", "BB_0406")

# Main processing loop for each prot
for (p in seq_along(prot_list)) {
  prot <- prot_list[p]
  
  cat("\n=== Processing protein:", prot, "===\n")
  
  # Load data
  data_path <- file.path(cwd, "data", prot)
  results_path <- file.path(cwd, "results", prot)
  
  # Check if directories exist
  if (!dir.exists(data_path)) {
    warning(paste("Data directory not found for protein:", prot))
    next
  }
  if (!dir.exists(results_path)) {
    dir.create(results_path, recursive = TRUE)
  }
  
  setwd(data_path)
  
  # Load input data
  input_files <- list.files(pattern = "input_data.RData")
  if (length(input_files) == 0) {
    warning(paste("No input_data.RData file found for proteim:", prot))
    next
  }
  
  load(input_files[1])
  
  # Check if required objects exist
  required_objects <- c("X_OHE", "Y", "fs_X_OHE", "X_O_ImpQ1", "X_O_ImpQ2", 
                        "X_O_ImpQ3", "O_VarQ1", "O_VarQ2", "O_VarQ3")
  missing_objects <- required_objects[!sapply(required_objects, exists)]
  
  if (length(missing_objects) > 0) {
    warning(paste("Missing required objects for protein", prot, ":", 
                  paste(missing_objects, collapse = ", ")))
    next
  }
  
  # Apply SMOTE for class balancing
  cat("Applying SMOTE for class balancing...\n")
  new_X <- upSample(X_OHE, Y)
  new_Y <- new_X$Class

  
  # Remove the Class column from features
  new_X$Class <- NULL
  
  # Set up different feature subsets
  X <- list(
    fs055 = new_X[, colnames(fs_X_OHE), drop = FALSE], 
    ImpQ1 = new_X[, colnames(X_O_ImpQ1), drop = FALSE], 
    ImpQ2 = new_X[, colnames(X_O_ImpQ2), drop = FALSE], 
    ImpQ3 = new_X[, colnames(X_O_ImpQ3), drop = FALSE],
    VarQ1 = new_X[, colnames(O_VarQ1), drop = FALSE],
    VarQ2 = new_X[, colnames(O_VarQ2), drop = FALSE], 
    VarQ3 = new_X[, colnames(O_VarQ3), drop = FALSE]
  )
  Y <- new_Y
  
  # Create train/test split indexes using repeated cross-validation
  set.seed(100)
  idx <- createMultiFolds(Y, k = 4, times = 25)
  # Train models for each feature subset
  for (ss in names(X)) {
    cat("\n--- Processing feature subset:", ss, "---\n")
    
    start_time <- Sys.time()
    start_ml_time <- format(start_time, "%Y_%m_%d_%H_%M_%S")
    
    # Prepare training and testing data for all folds
    cat("Preparing data splits...\n")
    O_Training.X = foreach (i = names(idx) ) %do%  {X[[ss]][idx[[i]],]}
    O_Testing.X  = foreach (i = names(idx) ) %do% {X[[ss]][-idx[[i]],]}
    Training.Y = foreach (i = names(idx) ) %do% {Y[idx[[i]]]}
    Testing.Y = foreach (i = names(idx) ) %do% {Y[-idx[[i]]]}
    
    # Set up training controls
    ctrloo <- trainControl(
      method = 'LOOCV', 
      classProbs = TRUE, 
      savePredictions = TRUE,
      returnResamp = "all", 
      summaryFunction = MySummary
    )
    
    ctrbs <- trainControl(
      method = 'boot', 
      number = 100, 
      classProbs = TRUE, 
      savePredictions = TRUE,
      returnResamp = "all", 
      summaryFunction = MySummary
    )
    
    # Model specifications
    models <- c('PLS', "RF", "SVM", 'GLM', "PCANN")
    m_methods <- c('pls', 'rf', 'svmRadialSigma', 'glmnet', "pcaNNet")
    
    # Register parallel backend
    registerDoParallel(ncores)
    # Train each model type
    for (m in seq_along(models)) {
      cat("Training", models[m], "model...\n")
      
      # LOOCV training
      O_LOOfit <- foreach(i = 1:length(idx), .packages = required_packages) %dopar% {
        tryCatch({
          # Train model
          fit <- train(
            O_Training.X[[i]], Training.Y[[i]],
            method = m_methods[m],
            tuneLength = 15,
            preProc = NULL,
            metric = "ROC",
            trControl = ctrloo
          )
          
          # Get training performance
          train_perf <- getTrainPerf(fit)
          
          # Get test performance
          test_pred <- predict(fit, O_Testing.X[[i]], type = "prob")
          test_pred_class <- predict(fit, O_Testing.X[[i]])
          
          test_perf <- t(MySummary(
            data.frame(
              test_pred,
              pred = test_pred_class,
              obs = Testing.Y[[i]]
            )
          ))
          
          # Variable importance
          vip_scores <- varImp(fit, scale = TRUE)
          vip_scores <- vip_scores[["importance"]]
          
          list(
            fit = fit,
            Train.Perf = train_perf,
            Test.Perf = test_perf,
            VIP_scores = vip_scores
          )
        }, error = function(e) {
          warning(paste("Error in LOOCV for", models[m], "fold", i, ":", e$message))
          return(NULL)
        })
      }
      
      # Bootstrap training
      O_Bsfit <- foreach(i = 1:length(idx), .packages = required_packages) %dopar% {
        tryCatch({
          # Train model
          fit <- train(
            O_Training.X[[i]], Training.Y[[i]],
            method = m_methods[m],
            preProc = NULL,
            metric = "ROC",
            trControl = ctrbs
          )
          
          # Get training performance
          train_perf <- getTrainPerf(fit)
          
          # Get test performance
          test_pred <- predict(fit, O_Testing.X[[i]], type = "prob")
          test_pred_class <- predict(fit, O_Testing.X[[i]])
          
          test_perf <- t(MySummary(
            data.frame(
              test_pred,
              pred = test_pred_class,
              obs = Testing.Y[[i]]
            )
          ))
          
          # Variable importance
          vip_scores <- varImp(fit, scale = TRUE)
          vip_scores <- vip_scores[["importance"]]
          
          list(
            fit = fit,
            Train.Perf = train_perf,
            Test.Perf = test_perf,
            VIP_scores = vip_scores
          )
        }, error = function(e) {
          warning(paste("Error in Bootstrap for", models[m], "fold", i, ":", e$message))
          return(NULL)
        })
      }
      
      # Save results
      assign(paste0(ss, "_", models[m], "_O_LOOfit"), O_LOOfit)
      assign(paste0(ss, "_", models[m], "_O_bsfit"), O_Bsfit)
      
      # Clean up intermediate objects
      rm(O_LOOfit, O_Bsfit)
    }
    
    # Stop parallel processing
    stopImplicitCluster()
    
    # Record end time
    end_time <- Sys.time()
    end_ml_time <- format(end_time, "%Y_%m_%d_%H_%M_%S")
    
    # Print timing information
    cat("\nFeature subset:", ss, "\n")
    cat("Start time:", start_ml_time, "\n")
    cat("End time:", end_ml_time, "\n")
    cat("Total duration:", format(end_time - start_time), "\n")
    
    # Save results
    save.image(paste0(results_path,"/",ss, "_results.RData"))
    
    # Clean up model objects to free memory
    rm(list = ls(pattern = "_O_bsfit|_O_LOOfit"))
  }
  
  # Clean up for next protein
  rm(list = setdiff(ls(), c("prot_list", "p", "ncores", "cwd", "required_packages", "MySummary")))
}

cat("\n=== All proteins processed successfully ===\n")