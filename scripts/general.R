# Libraries ----
required_packages <- c("tidyverse", "readxl", "moments", "doParallel", "caret", "BGLR", "gbm", "coda", "xgboost")
#pak::pak(required_packages)
pacman::p_load(char = required_packages, install = FALSE)

## Raw data ----
genotype_table <- 
  read_excel("data/papaya_sample_meta_data.xlsx",
             sheet = "meta_data")

calibration_table <- read_excel("data/Results_22_23.xlsx", "densities") |> 
  select(Hexanal:Linalool) |> 
  t() |> as.data.frame() |> rownames_to_column("species") |> 
  set_names(c("species", "densities", "OT", "purity"))

res_sheets <- c("Aug22", "Feb23", "Feb24")

VOCs <- res_sheets |> 
  # Read and combine sheets
  map(~read_excel("data/Results_22_23.xlsx", .x)  |> 
        pivot_longer(cols = (Hexanal:BITC),
                     names_to = "species",
                     values_to = "concentration") |>
        mutate(assay = .x)) |>
  list_rbind() |> 
# join data
left_join(calibration_table) |>
  left_join(genotype_table) |>
  
  # Calculate concentrations
  mutate(
    conc_ppb = (concentration * densities * purity) * (0.00425/0.0042), # ug/L -> ug/KgFW = ug/L * volume(L)/mass(Kg)
    species = str_replace_all(species, "-", "_")) |> 
  select(-c(concentration, OT, purity, densities, label, Genotype)) |> 
  pivot_wider(names_from = species, values_from = conc_ppb) |> 
  rename(Genotype = Genotype2)

# cut VOC variables
# TECHNICAL NOTE: no_aroma variables were determined by calculating the odour activity values using published odour thresholds. Only species with odour thresholds >1 in one or more samples were used for modelling. We calculated these values based on the matrix corrected concentrations. Interestingly, when using the raw instrument concentrations more VOCs fell below the odour activity threshold: "P_cymene", "Beta_Cyclocitral","Terpinene","Methyl_octanoate", "Benzaldehyde", "Phenylacetaldehyde"
no_aroma <- c("Methyl_octanoate", "Benzaldehyde", "Phenylacetaldehyde")
# cut sensory variables
poor_lexicon_use <- c("prickly_AT", "metallic_AT")

# Matrix corrected VOC concentrations
all_vars_corrected_vocs <- readRDS("data/all_vars.RDS")

# Raw instrument VOC conentrations
all_vars_raw_voc <- all_vars_corrected_vocs |> 
  select(-c(Hexanal:BITC, label, Genotype2)) |> 
  left_join(VOCs) |> 
  relocate(Hexanal:BITC, .before = aroma_intensity_AR) |> 
  select(-(all_of(c(no_aroma, poor_lexicon_use)))) |> 
  mutate(across(where(is.numeric), ~ replace_na(.x, 0.000001)))

saveRDS(all_vars_raw_voc, "data/all_raw_data.RDS")

## Modelling data sets ----
sensory_data <- readRDS("app_data_files/all_raw_data.RDS")

# Consumer liking scores with matrix corrected vars removed
consumer_liking <- readRDS("data/model_data.RDS") |> 
  rownames_to_column(var = "ID") |> 
  select(ID, emmean) |> 
  rename("Y" = "emmean")

# check data distribution
# Bayesian regression is parametric
# GBM and XGBoost is non-parametric
# Therefore only the variable exploration requires transformations related to parametric assumptions
# however the feature selection is primarily based on the GBM so potentially this step can be simplified
# so the bayes data should be different to the gbm data


bayes_data <- sensory_data |>
  mutate(across(where(is.numeric), \(x) {
    x <- replace_na(x, 0.000001)
    s <- skewness(x)
    if (s > 0.5 & s <= 1) {
      sqrt(x)
    } else if (s > 1) {
      log2(x)
    } else {
      as.numeric(x)
    }
  }))

gbm_data <- sensory_data
# 
# # Consumer liking data
# liking_model_data <- left_join(consumer_liking, sensory_data)

## Feature selection ----
# sensory features that most influence liking may be included alongside liking
# metabolites associated with important sensory features can be modeled more accurately because of greater sample sizes
# metabolites associated with liking have limited sample size but can be assessed alongside sensory attribute predictions

# setup parallel processing for efficient computation
doParallel::registerDoParallel(6)
#doParallel::stopImplicitCluster()
## might help with reproducibility
#xpectr::set_test_seed(1)

# Liking feature importance based on sensory
# Data matrix set up

meta_cols <- c("ID", "assay", "Genotype", "Flesh")

sensory_names <- sensory_data |> select(aroma_intensity_AR:sweet_AT) |> names()
metabolite_names <- sensory_data |> select(Glucose:BITC) |> names()

bayes_carX <- left_join(consumer_liking, bayes_data) |> 
  select(-(all_of(meta_cols))) |> 
  mutate(across(where(is.numeric), ~ as.numeric(scale(.)))) |> 
  as.matrix()

gbm_carX <- left_join(consumer_liking, gbm_data) |> select(-(all_of(meta_cols))) 

liking_gbm_data <- list(mets = gbm_carX[,!colnames(gbm_carX) %in% sensory_names],
                        sens = gbm_carX[,!colnames(gbm_carX) %in% metabolite_names])
red_index <- which(left_join(consumer_liking, bayes_data)[,"Flesh"] == "Red")

# Model grid search parameters
model_params <- list(
  n_trees = c(250, 500, 750, 1000),
  int_depth = c(5, 7, 10),
  shrinkage = c(0.01, 0.02, 0.05),
  min_obs = c(3, 5),
  cv_folds = 10,
  bglr_iter = 30000,
  bglr_burnin = 10000
)
# GBM
grid <- expand.grid(
  .n.trees = model_params$n_trees,
  .interaction.depth = model_params$int_depth,
  .shrinkage = model_params$shrinkage,
  .n.minobsinnode = model_params$min_obs
)
# set seeds for reproducibility
set.seed(0)
seeds <- vector(mode = "list", length = 51)
for(i in 1:50) seeds[[i]] <- sample.int(1000, 18)
seeds[[51]] <- sample.int(1000,1)

mets_gbm_model <- train(
  Y ~ .,
  data = gbm_carX,
  method = "gbm",
  trControl = trainControl(method = "repeatedcv",
                           number = model_params$cv_folds,
                           repeats = 5,
                           seeds = seeds),
  tuneGrid = grid
)
# Caret outputs for model analysis
getTrainPerf(mets_gbm_model)
mets_gbm_model$results
mets_gbm_model$bestTune
mets_gbm_model$resample
mets_gbm_model$finalModel
varImp(mets_gbm_model)$importance 


# BGLR
bglr_model <- BGLR::BGLR(
  y = bayes_carX[,"Y"],
  ETA = list(list(X = bayes_carX[,!colnames(bayes_carX) == "Y" & !colnames(bayes_carX) %in% metabolite_names], model = 'BayesA')),
  nIter = model_params$bglr_iter,
  burnIn = model_params$bglr_burnin,
  thin = 100,
  saveAt = 'data/',
  df0 = 5, # Allows big β’s if data support them; strong shrinkage of small effects
  R2 = 0.5, # variables expected to explain 50% of y
  verbose = F
)
plot_data <- data.frame(x = bglr_model$yHat,y = bglr_model$y)
plot(plot_data$x, plot_data$y)
plot(bglr_model$ETA[[1]]$b)
plot(bglr_model$yHat)
plot(bglr_model$varE)

# each iteration of regression model applies a new intercept, new regression coefficients for each variable + variances of coefficients and a residual variance term. Each term is sampled from a probable distribution (chi-squared distribution) that minimises the probability of high terms and is initailly guided by the R2 term that suggests how much we expect the y term to be influence by the variables.
varE_samples <- read.table("data/varE.dat", header = FALSE)
mu_samples <- read.table("data/mu.dat", header = FALSE)
# Convert to mcmc objects for coda diagnostics
varE_mcmc <- mcmc(varE_samples)
mu_mcmc <- mcmc(mu_samples)
# look for difference between start and end of chains: >2 indicate non-convergence
geweke.diag(varE_mcmc)
geweke.diag(mu_mcmc)
# can be improved by increasing burn-in (not needed here)
# more precise test to show stability across iterations
heidel.diag(varE_mcmc)
heidel.diag(mu_mcmc)
# stationary test indicates convergence
# half width test indicates precision (degree of monte carlo error)
# can be improved by increasing effective samples (e.g. lower thin)
# half-width of mu can fail to pass because of sample size
# If stationarity passes and trace/ACF plots are good, half‑width failures for μ can be safely ignored.

# The first plot shows the variation across iterations. Does it look like a “fat hairy caterpillar” (good)? Is there drift or trending (bad)?
# The second plot shows the most probable error term for the model (single peak). Two peaks means that the model hasn't reached a single most probable solution
plot(varE_mcmc, main = "Trace: Residual Variance (varE)", 
     xlab = "Iteration", ylab = "varE", col = "steelblue", type = "l")
# this plot shows the model error term across the MCMC sampling. lag 0 compares each values with itself, lag 1 compares each value with the previous value (e.g., sample 2 with sample 1) , lag 2 compare sample 3 with 1 etc. High autocorrelation is bad because there is limited effective sampling meaning that the results are biased and probabilities are limited.
acf(varE_samples[,1], main = "Autocorrelation: varE", 
    col = "darkgreen", lwd = 2)

# explore the same themes for the model intercept

plot(mu_mcmc, main = "Trace: Intercept (mu)", 
     xlab = "Iteration", ylab = "mu", col = "steelblue", type = "l")

acf(mu_samples[,1], main = "Autocorrelation: mu", 
    col = "darkgreen", lwd = 2)

# Run models for multiple combinations ------
meta_cols <- c("ID", "assay", "Genotype", "Flesh")

sensory_names <- sensory_data |> select(aroma_intensity_AR:sweet_AT) |> names()

metabolite_names <- sensory_data |> select(Glucose:BITC) |> names()

required_datasets <- c("mets", "sens", "mets_red", "sens_red")

liking_data <- setNames(nm = required_datasets) |>
  map(function(data_set) {
    
    if (str_detect(data_set, "red")) {
      bayes_df <- left_join(consumer_liking, bayes_data) |> filter(Flesh == "Red")
      gbm_df <- left_join(consumer_liking, gbm_data) |> filter(Flesh == "Red")
    } else {
      bayes_df <- left_join(consumer_liking, bayes_data)
      gbm_df <- left_join(consumer_liking, gbm_data)
    }
    
    if (str_detect(data_set, "mets")) {
      bayes <- bayes_df |>
        select(-c(all_of(meta_cols), all_of(sensory_names))) |>
        mutate(across(where(is.numeric), ~ as.numeric(scale(.)))) |>
        as.matrix()
      
      gbm <- gbm_df |> select(-c(all_of(meta_cols), all_of(sensory_names)))
    } else {
      bayes <- bayes_df |>
        select(-c(all_of(meta_cols), all_of(metabolite_names))) |>
        mutate(across(where(is.numeric), ~ as.numeric(scale(.)))) |>
        as.matrix()
      
      gbm <- gbm_df |> select(-c(all_of(meta_cols), all_of(metabolite_names)))
    }
    
    combined_gbm_bayes <- list(bayes = bayes, gbm = gbm)
    
    return(combined_gbm_bayes)
  })

# Model grid search parameters
model_params <- list(
  n_trees = c(250, 500, 750, 1000),
  int_depth = c(5, 7, 10),
  shrinkage = c(0.01, 0.02, 0.05),
  min_obs = c(3, 5),
  cv_folds = 10,
  bglr_iter = 30000,
  bglr_burnin = 10000
)
# GBM
grid <- expand.grid(
  .n.trees = model_params$n_trees,
  .interaction.depth = model_params$int_depth,
  .shrinkage = model_params$shrinkage,
  .n.minobsinnode = model_params$min_obs
)
# set seeds for reproducibility
set.seed(0)
seeds <- vector(mode = "list", length = 51)
for(i in 1:50) seeds[[i]] <- sample.int(1000, 18)
seeds[[51]] <- sample.int(1000,1)

var_imp_models <- setNames(nm = required_datasets) |>
  map(function(data_set) {
    if (str_detect(data_set, "red")) {
      mets <- liking_data$mets_red
      sens <- liking_data$sens_red
    } else {
      mets <- liking_data$mets
      sens <- liking_data$sens
    }
    
    if (str_detect(data_set, "mets")) {
      gbm_df <- mets$gbm
      bayes_X <- mets$bayes
    } else {
      gbm_df <- sens$gbm
      bayes_X <- sens$bayes
    }
    
    # Run both models
    gbm_model <- train(
      Y ~ .,
      data = gbm_df,
      method = "gbm",
      trControl = trainControl(
        method = "repeatedcv",
        number = model_params$cv_folds,
        repeats = 5,
        seeds = seeds, 
        verboseIter = F
      ),
      tuneGrid = grid,
      verbose = F
    )
    
    # BGLR
    bglr_model <- BGLR::BGLR(
      y = bayes_X[, "Y"],
      ETA = list(list(X = bayes_X[, !colnames(bayes_X) == "Y"], model = 'BayesA')),
      nIter = model_params$bglr_iter,
      burnIn = model_params$bglr_burnin,
      thin = 100,
      saveAt = 'data/',
      df0 = 5,
      R2 = 0.5,
      verbose = F
    )
    
    results <- bind_rows(
      varImp(gbm_model)$importance %>%
        rownames_to_column("var") %>%
        mutate(model = paste0(data_set, "_gbm")),
      tibble(
        var = names(bglr_model$ETA[[1]]$b),
        Overall = bglr_model$ETA[[1]]$b,
        model = paste0(data_set, "_beta")
      )
    ) %>%
      rename(coord = Overall)
    
    
    return(results)
  })

## plot liking vars ------
# plot_liking_imp <- var_imp_models |> 
#   list_rbind() |> 
#   left_join(read_xlsx("data/VOC_meta.xlsx"))

## Liking model functions ------

fit_papaya_liking_model <- function(boost_data,
                                    seed = 11,
                                    test_prop = 0.2,
                                    cv_folds = 5) {
  
  # Data preparation
  model_data <- boost_data 
  
  # Train-test split
  set.seed(seed)
  train_index <- createDataPartition(model_data$Y, p = 1 - test_prop, list = FALSE)[, 1]
  
  X_train_full <- as.matrix(model_data[train_index, -1])
  X_test <- as.matrix(model_data[-train_index, -1])
  y_train_full <- model_data$Y[train_index]
  y_test <- model_data$Y[-train_index]
  # 
  # # Expanded hyperparameter grid
  # param_grid <- expand.grid(
  #   eta = c(0.01, 0.05, 0.1, 0.15),
  #   max_depth = c(1:3), #c(2, 3, 4)
  #   subsample = c(5:7/10), #c(0.7, 0.8, 0.9)
  #   colsample_bytree = c(5:7/10), #c(0.7, 0.8, 0.9)
  #   lambda = c(1,5,10), #c(0, 0.1, 1)
  #   alpha = c(1,5,10), #c(0, 0.1, 0.5)
  #   min_child_weight = c(3,5,7), #c(1, 3, 5)
  #   nrounds = c(20, 50, 100) #c(50, 100, 200)
  # )
  
  # Expanded hyperparameter grid
  param_grid <- expand.grid(
    eta = c(0.01, 0.05, 0.1, 0.15),
    max_depth = c(2, 3, 4), 
    subsample = c(0.7, 0.8, 0.9), 
    colsample_bytree = c(0.7, 0.8, 0.9), 
    nrounds = c(50, 100, 200)
  )
  
  # Storage for results
  cv_results_table <- tibble()
  
  cat("Testing", nrow(param_grid), "parameter combinations with", cv_folds, "-fold CV...\n")
  
  for (i in seq_len(nrow(param_grid))) {
    params <- list(
      objective = "reg:squarederror",
      eta = param_grid$eta[i],
      max_depth = param_grid$max_depth[i],
      subsample = param_grid$subsample[i],
      colsample_bytree = param_grid$colsample_bytree[i]#,
      # lambda = param_grid$lambda[i],
      # alpha = param_grid$alpha[i],
      # min_child_weight = param_grid$min_child_weight[i]
    )
    
    # Perform cross-validation and get averaged metrics
    cv_metrics <- perform_cv_validation(
      X_train_full = X_train_full,
      y_train_full = y_train_full,
      params = params,
      nrounds = param_grid$nrounds[i],
      cv_folds = cv_folds,
      seed = seed
    )
    
    # Calculate percentage differences between train and validation
    rmse_diff_pct <- abs(cv_metrics$avg_val_rmse - cv_metrics$avg_train_rmse) / cv_metrics$avg_train_rmse * 100
    r2_diff_pct <- abs(cv_metrics$avg_val_r2 - cv_metrics$avg_train_r2) / max(cv_metrics$avg_train_r2, 0.001) * 100
    
    # Store comprehensive results
    current_result <- tibble(
      combination = i,
      eta = param_grid$eta[i],
      max_depth = param_grid$max_depth[i],
      subsample = param_grid$subsample[i],
      colsample_bytree = param_grid$colsample_bytree[i],
      # lambda = param_grid$lambda[i],
      # alpha = param_grid$alpha[i],
      # min_child_weight = param_grid$min_child_weight[i],
      nrounds = param_grid$nrounds[i],
      avg_best_nrounds = cv_metrics$avg_best_nrounds,
      
      # CV Training metrics (averaged across folds)
      cv_train_rmse = cv_metrics$avg_train_rmse,
      cv_train_r2 = cv_metrics$avg_train_r2,
      cv_train_mae = cv_metrics$avg_train_mae,
      cv_train_corr = cv_metrics$avg_train_corr,
      
      # CV Validation metrics (averaged across folds)
      cv_val_rmse = cv_metrics$avg_val_rmse,
      cv_val_r2 = cv_metrics$avg_val_r2,
      cv_val_mae = cv_metrics$avg_val_mae,
      cv_val_corr = cv_metrics$avg_val_corr,
      
      # Standard deviations across folds
      cv_train_rmse_sd = cv_metrics$train_rmse_sd,
      cv_val_rmse_sd = cv_metrics$val_rmse_sd,
      cv_train_r2_sd = cv_metrics$train_r2_sd,
      cv_val_r2_sd = cv_metrics$val_r2_sd,
      
      # Differences between train and validation
      rmse_diff = cv_metrics$avg_val_rmse - cv_metrics$avg_train_rmse,
      r2_diff = cv_metrics$avg_val_r2 - cv_metrics$avg_train_r2,
      rmse_diff_pct = rmse_diff_pct,
      r2_diff_pct = r2_diff_pct,
      
      # Overfitting indicators
      overfitting_ratio = cv_metrics$avg_val_rmse / cv_metrics$avg_train_rmse,
      
      # Stability indicators (lower CV = more stable)
      rmse_stability = cv_metrics$val_rmse_sd / cv_metrics$avg_val_rmse,
      r2_stability = cv_metrics$val_r2_sd / max(abs(cv_metrics$avg_val_r2), 0.001)
    )
    
    cv_results_table <- bind_rows(cv_results_table, current_result)
    
    if (i %% 50 == 0) cat("Completed", i, "combinations\n")
  }
  
  # Add ranking columns for easy sorting
  cv_results_table <- cv_results_table |>
    mutate(
      rank_val_r2 = rank(-cv_val_r2),
      rank_val_rmse = rank(cv_val_rmse),
      rank_stability = rank(rmse_stability + r2_stability),
      rank_overfitting = rank(overfitting_ratio)
    ) |>
    arrange(rank_val_r2, rank_val_rmse)
  
  # Return results table and training function
  list(
    cv_results = cv_results_table,
    train_final_model = function(selected_row) {
      train_final_model_with_params(
        selected_row = selected_row,
        X_train_full = X_train_full,
        X_test = X_test,
        y_train_full = y_train_full,
        y_test = y_test,
        seed = seed
      )
    },
    data_splits = list(
      train_idx = train_index,
      test_idx = setdiff(seq_len(nrow(model_data)), train_index)
    )
  )
}

# Enhanced cross-validation function with comprehensive metrics
perform_cv_validation <- function(X_train_full, y_train_full, params, nrounds, cv_folds, seed) {
  
  cv_results <- list()
  
  # Create balanced folds
  set.seed(seed)
  cv_folds_indices <- createFolds(y_train_full, k = cv_folds, list = TRUE)
  
  for (fold in 1:cv_folds) {
    # Split data for this fold
    val_indices <- cv_folds_indices[[fold]]
    train_indices <- setdiff(1:length(y_train_full), val_indices)
    
    X_train_cv <- X_train_full[train_indices, ]
    X_val_cv <- X_train_full[val_indices, ]
    y_train_cv <- y_train_full[train_indices]
    y_val_cv <- y_train_full[val_indices]
    
    # Create DMatrix objects
    dtrain_cv <- xgb.DMatrix(data = X_train_cv, label = y_train_cv)
    dval_cv <- xgb.DMatrix(data = X_val_cv, label = y_val_cv)
    
    # Train model with early stopping
    watchlist <- list(train = dtrain_cv, val = dval_cv)
    
    set.seed(seed + fold)
    model_cv <- xgb.train(
      params = params,
      data = dtrain_cv,
      nrounds = nrounds,
      evals = watchlist,
      early_stopping_rounds = 15,
      verbose = 0
    )
    
    # Generate predictions
    train_pred_cv <- predict(model_cv, dtrain_cv)
    val_pred_cv <- predict(model_cv, dval_cv)
    
    best_iter <- attributes(model_cv)$early_stop$best_iteration
    if (is.null(best_iter)) {
      best_iter <- nrounds
    }
    
    # Calculate comprehensive metrics for this fold
    fold_metrics <- calculate_fold_metrics(
      train_y_true = y_train_cv,
      train_y_pred = train_pred_cv,
      val_y_true = y_val_cv,
      val_y_pred = val_pred_cv,
      best_nrounds = best_iter
    )
    
    cv_results[[fold]] <- fold_metrics
  }
  
  # Calculate averages and standard deviations across folds
  avg_metrics <- list(
    avg_train_rmse = mean(sapply(cv_results, function(x) x$train_rmse)),
    avg_train_r2 = mean(sapply(cv_results, function(x) x$train_r2)),
    avg_train_mae = mean(sapply(cv_results, function(x) x$train_mae)),
    avg_train_corr = mean(sapply(cv_results, function(x) x$train_corr)),
    
    avg_val_rmse = mean(sapply(cv_results, function(x) x$val_rmse)),
    avg_val_r2 = mean(sapply(cv_results, function(x) x$val_r2)),
    avg_val_mae = mean(sapply(cv_results, function(x) x$val_mae)),
    avg_val_corr = mean(sapply(cv_results, function(x) x$val_corr)),
    
    train_rmse_sd = sd(sapply(cv_results, function(x) x$train_rmse)),
    train_r2_sd = sd(sapply(cv_results, function(x) x$train_r2)),
    val_rmse_sd = sd(sapply(cv_results, function(x) x$val_rmse)),
    val_r2_sd = sd(sapply(cv_results, function(x) x$val_r2)),
    
    avg_best_nrounds = round(mean(sapply(cv_results, function(x) x$best_nrounds)))
  )
  
  return(avg_metrics)
}

# Helper function to calculate metrics for a single fold
calculate_fold_metrics <- function(train_y_true, train_y_pred, val_y_true, val_y_pred, best_nrounds) {
  
  # Training metrics
  train_residuals <- train_y_true - train_y_pred
  train_rmse <- sqrt(mean(train_residuals^2))
  train_r2 <- 1 - sum(train_residuals^2) / sum((train_y_true - mean(train_y_true))^2)
  train_mae <- mean(abs(train_residuals))
  train_corr <- cor(train_y_pred, train_y_true)
  
  # Validation metrics
  val_residuals <- val_y_true - val_y_pred
  val_rmse <- sqrt(mean(val_residuals^2))
  val_r2 <- 1 - sum(val_residuals^2) / sum((val_y_true - mean(val_y_true))^2)
  val_mae <- mean(abs(val_residuals))
  val_corr <- cor(val_y_pred, val_y_true)
  
  list(
    train_rmse = train_rmse,
    train_r2 = train_r2,
    train_mae = train_mae,
    train_corr = train_corr,
    val_rmse = val_rmse,
    val_r2 = val_r2,
    val_mae = val_mae,
    val_corr = val_corr,
    best_nrounds = best_nrounds
  )
}

calculate_comprehensive_metrics <- function(train_y_true, train_y_pred,
                                            test_y_true, test_y_pred,
                                            n_predictors) {
  
  # Helper function for metrics calculation
  calc_metrics <- function(actual, predicted, dataset_name) {
    residuals <- actual - predicted
    n <- length(actual)
    
    # Core metrics
    mae <- mean(abs(residuals))
    mse <- mean(residuals^2)
    rmse <- sqrt(mse)
    r2 <- 1 - sum(residuals^2) / sum((actual - mean(actual))^2)
    
    # Adjusted R-squared
    adj_r2 <- if (n > (n_predictors + 1)) {
      1 - ((1 - r2) * (n - 1) / (n - n_predictors - 1))
    } else {
      NA_real_
    }
    
    # Additional metrics
    nrmse <- rmse / diff(range(actual))
    corr <- cor(predicted, actual)
    rpd <- sd(actual) / rmse
    
    tibble(
      n = n,
      mae = mae,
      rmse = rmse,
      r2 = r2,
      adj_r2 = adj_r2,
      nrmse = nrmse,
      corr = corr,
      rpd = rpd,
      mean_residual = mean(residuals),
      residual_std = sd(residuals)
    )
  }
  
  # Calculate for both datasets
  train_metrics <- calc_metrics(train_y_true, train_y_pred, "train")
  test_metrics <- calc_metrics(test_y_true, test_y_pred, "test")
  
  # Combine with clear naming
  bind_cols(
    train_metrics |> rename_with(~ paste0("train_", .x)),
    test_metrics |> rename_with(~ paste0("test_", .x))
  )
}

# Function to train final model with selected parameters
train_final_model_with_params <- function(selected_row, X_train_full, X_test, y_train_full, y_test, seed) {
  
  # Extract parameters from selected row
  final_params <- list(
    objective = "reg:squarederror",
    eta = selected_row$eta,
    max_depth = selected_row$max_depth,
    subsample = selected_row$subsample,
    colsample_bytree = selected_row$colsample_bytree#,
    # lambda = selected_row$lambda,
    # alpha = selected_row$alpha,
    # min_child_weight = selected_row$min_child_weight
  )
  
  # Create DMatrix objects
  dtrain_full <- xgb.DMatrix(data = X_train_full, label = y_train_full)
  dtest <- xgb.DMatrix(data = X_test, label = y_test)
  
  # Train final model
  set.seed(seed)
  final_model <- xgb.train(
    params = final_params,
    data = dtrain_full,
    nrounds = selected_row$avg_best_nrounds,
    verbose = 0
  )
  
  # Final predictions and metrics
  train_pred <- predict(final_model, dtrain_full)
  test_pred <- predict(final_model, dtest)
  
  final_metrics <- calculate_comprehensive_metrics(
    train_y_true = y_train_full,
    train_y_pred = train_pred,
    test_y_true = y_test,
    test_y_pred = test_pred,
    n_predictors = ncol(X_train_full)
  )
  
  # Feature importance
  importance <- xgb.importance(
    feature_names = colnames(X_train_full),
    model = final_model
  )
  
  list(
    model = final_model,
    metrics = final_metrics,
    params = final_params,
    nrounds = selected_row$avg_best_nrounds,
    importance = importance,
    predictions = list(train = train_pred, test = test_pred)
  )
}


## Fit liking models -----


# Create all combinations of datasets and top variable counts
combo_results3 <- expand_grid(
  # dataset = c("mets", "mets_red"),
  dataset = c("mets_red"),
  top_variables = c(2, 3, 4, 5)
) |>
  pmap(function(dataset, top_variables) {
    
    # Get important variables for this dataset
    important_variables <- var_imp_models[[dataset]] |> 
      filter(str_detect(model, "(?=.*mets)(?=.*gbm)")) |> 
      arrange(desc(coord))
    
    # Select top N variables
    explanatory_vars <- important_variables |>
      slice_head(n = top_variables) |>
      pull(var)
    
    # Get the appropriate data
    if (dataset == "mets") {
      boost_data <- liking_gbm_data$mets |> 
        select(Y, all_of(explanatory_vars)) 
    } else {
      boost_data <- liking_gbm_data$mets[red_index,] |> 
        select(Y, all_of(explanatory_vars)) 
    }
    
    # Fit model
    xgb_results <- fit_papaya_liking_model(boost_data, seed = 11, cv_folds = 5)
    
    # Select best model
    best_model <- xgb_results$cv_results %>%
      filter(
        rmse_diff_pct <= quantile(rmse_diff_pct, 0.8),
        overfitting_ratio <= quantile(overfitting_ratio, 0.8)
      ) %>%
      mutate(
        # Add small epsilon to prevent division by zero
        val_rmse_norm = (max(cv_val_rmse) - cv_val_rmse) / 
          (max(cv_val_rmse) - min(cv_val_rmse) + 1e-10),
        overfitting_norm = (max(rmse_diff_pct) - rmse_diff_pct) / 
          (max(rmse_diff_pct) - min(rmse_diff_pct) + 1e-10),
        stability_norm = (max(cv_val_rmse_sd) - cv_val_rmse_sd) / 
          (max(cv_val_rmse_sd) - min(cv_val_rmse_sd) + 1e-10),
        
        final_score = 0.60 * val_rmse_norm + 0.25 * overfitting_norm + 0.15 * stability_norm
      ) %>%
      arrange(desc(final_score)) %>%
      slice_head(n = 1)
    
    # Train final model
    final_model_results <- xgb_results$train_final_model(best_model)
    
    # Format results
    performance_table <- final_model_results$metrics %>%
      pivot_longer(cols = everything(),
                   names_to = "metric",
                   values_to = "value") %>%
      mutate(
        dataset_type = ifelse(str_starts(metric, "train_"), "train", "test"),
        metric_clean = str_remove(metric, "^(train_|test_)")
      ) %>%
      select(-metric) %>%
      pivot_wider(names_from = metric_clean, values_from = value) %>%
      select(dataset_type, n, mae, rmse, r2, adj_r2, nrmse, corr, rpd, 
             mean_residual, residual_std) |>
      mutate(
        across(where(is.numeric), ~ round(.x, 2)),
        dataset = dataset,
        n_variables = top_variables,  
        model_vars = paste(explanatory_vars, collapse = ", "),
        final_model = list(final_model_results)
      )
    
    return(performance_table)
  }) |> 
  list_rbind()

mets_model <- combo_results$final_model[[1]]$model
mets_model_params <- combo_results$final_model[[1]]$params
mets_model_importance <- combo_results$final_model[[1]]$importance

combo_results2 # less conservative hyperparameters
combo_results3 # least conservative hyperparamters

# mets model using conservative parameters. mets_red model from less conservative parameters
final_models <- bind_rows(
  combo_results |> filter(model_vars == "Sulcatone, Linalool_oxide"),
  combo_results3 |> filter(model_vars == "Linalool_oxide, Sulcatone, Brix")
)

saveRDS(final_models, "app_data_files/xgb_combos.RDS")
 
# ## Final Liking model -----
# boost_data <- model_data |>
#   select(Y, Sulcatone, Brix, Beta_ionone) |>
#   scale() |>
#   as.data.frame()
# 
# final_xgb_results <- fit_papaya_liking_model(boost_data, seed = 11, cv_folds = 5)
# 
# best_model <- final_xgb_results$cv_results %>%
#   # Step 1: Remove models with excessive overfitting (worst 20%)
#   filter(
#     rmse_diff_pct <= quantile(rmse_diff_pct, 0.8),
#     overfitting_ratio <= quantile(overfitting_ratio, 0.8)
#   ) %>%
#   # Step 2: Create balanced score
#   mutate(
#     # Normalise metrics to 0-1 scale
#     val_rmse_norm = (max(cv_val_rmse) - cv_val_rmse) / (max(cv_val_rmse) - min(cv_val_rmse)),
#     overfitting_norm = (max(rmse_diff_pct) - rmse_diff_pct) / (max(rmse_diff_pct) - min(rmse_diff_pct)),
#     stability_norm = (max(cv_val_rmse_sd) - cv_val_rmse_sd) / (max(cv_val_rmse_sd) - min(cv_val_rmse_sd)),
#     
#     # Weighted score: 60% validation accuracy, 25% overfitting, 15% stability
#     final_score = 0.60 * val_rmse_norm + 0.25 * overfitting_norm + 0.15 * stability_norm
#   ) %>%
#   # Step 3: Select top 10 models
#   arrange(desc(final_score)) %>%
#   slice_head(n = 1)
# 
# # use best model parameters to train final model
# selected_row <- best_model
# 
# selected_row[1:6] |>
#   mutate(across(where(is.numeric), ~round(.x, digits = 2))) |> 
#   gt() |> 
#   tab_caption("Final XGBoost model parameters")
# 
# # Train final model with selected parameters
# final_model_results <- final_xgb_results$train_final_model(selected_row)
# 
# # Make predictions on test set
# predictions <- final_model_results$predictions$test
# 
# performance_table <- final_model_results$metrics %>%
#   pivot_longer(cols = everything(),
#                names_to = "metric",
#                values_to = "value") %>%
#   mutate(
#     dataset = ifelse(str_starts(metric, "train_"), "train", "test"),
#     metric_clean = str_remove(metric, "^(train_|test_)")
#   ) %>%
#   select(-metric) %>%
#   pivot_wider(names_from = metric_clean, values_from = value) %>%
#   select(dataset,
#          n,
#          mae,
#          rmse,
#          r2,
#          adj_r2,
#          nrmse,
#          corr,
#          rpd,
#          mean_residual,
#          residual_std) |>
#   mutate(across(where(is.numeric), ~ round(.x, 2)))
# 
# # OLD Function to run models for a specific flesh type -------
# run_models <- function(data, flesh_type, met_cols, sens_cols, y_col = 41) {
#   # 
#   # Filter and prepare data
#   pap <- data %>%
#     filter(grepl(flesh_type, Flesh, perl = TRUE)) %>%
#     select(-Flesh) %>% 
#     scale()
#   
#   # Define model parameters based on flesh type
#   model_params <- if(flesh_type == "Yellow") {
#     list(
#       n_trees = c(500),
#       int_depth = c(5),
#       shrinkage = c(0.01),
#       min_obs = c(1),
#       cv_folds = 5,
#       bglr_iter = 30000,
#       bglr_burnin = 10000
#     )
#   } else {
#     list(
#       n_trees = c(150, 250, 500, 750, 1000),
#       int_depth = c(5, 7, 10),
#       shrinkage = c(0.01, 0.02, 0.05),
#       min_obs = c(3, 5),
#       cv_folds = 10,
#       bglr_iter = 30000,
#       bglr_burnin = 10000
#     )
#   }
#   
#   # Function to run a single model set (metabolites or sensory)
#   run_model_set <- function(cols, type) {
#     # Prepare data
#     carX <- as.matrix(cbind(Y = pap[,y_col], pap[,cols]))
#     
#     # GBM
#     grid <- expand.grid(
#       .n.trees = model_params$n_trees,
#       .interaction.depth = model_params$int_depth,
#       .shrinkage = model_params$shrinkage,
#       .n.minobsinnode = model_params$min_obs
#     )
#     
#     gbm_model <- train(
#       Y ~ .,
#       data = carX,
#       method = "gbm",
#       trControl = trainControl(method = "cv", number = model_params$cv_folds),
#       tuneGrid = grid
#     )
#     
#     # BGLR
#     bglr_model <- BGLR::BGLR(
#       y = pap[,y_col],
#       ETA = list(list(X = pap[,cols], model = 'BayesA')),
#       nIter = model_params$bglr_iter,
#       burnIn = model_params$bglr_burnin,
#       thin = 100,
#       saveAt = 'data/',
#       df0 = 5,
#       R2 = 0.5
#     )
#     
#     # Combine results
#     bind_rows(
#       varImp(gbm_model)$importance %>%
#         rownames_to_column("var") %>%
#         mutate(model = paste0(flesh_type, "_", type, "_gbm")),
#       tibble(
#         var = colnames(pap[,cols]),
#         Overall = bglr_model$ETA[[1]]$b,
#         model = paste0(flesh_type, "_", type, "_beta")
#       )
#     ) %>%
#       rename(coord = Overall)
#   }
#   
#   # Run both model sets and combine results
#   bind_rows(
#     run_model_set(met_cols, "mets"),
#     run_model_set(sens_cols, "sens")
#   )
# }
# 
# # run functions
# set.seed(11)
# # flesh_types <- c("Red", "Yellow", "Red_Yellow")
# flesh_types <- c("Red", "Yellow", "Red|Yellow")
# met_cols <- 1:22
# sens_cols <- 23:40
# 
# results <- flesh_types %>%
#   map(~run_models(model_data, ., met_cols, sens_cols)) %>%
#   list_rbind() %>% 
#   select(var, model, coord) %>% 
#   mutate(model = sub("|", "_", model, fixed = TRUE))
# save(results, file = "data/liking_var_model.RData")
# 
# # Top sensory characteristics modeled for metabolite importance
# 
# # Liking feature importance based on metabolites
# 
# 
# # Final model

## edit image -----
library(magick)
img <- image_read("app_data_files/Funding_logo.png")
img_no_bg <- image_transparent(img, "white", fuzz = 15)
image_write(img_no_bg, "app_data_files/Funding_logo.png")

# launch app ----
rsconnect::deployApp(                  
  appDir = ".",                             
  appPrimaryDoc = "scripts/Shiny_app_script.R",
  appFiles = c(
    "scripts/Shiny_app_script.R",
    "app_data_files/xgb_combos.RDS",
    "app_data_files/model_data.RDS",
    "app_data_files/all_raw_data.RDS",
    "app_data_files/Funding_logo.png",
    "app_data_files/Univiersity_logo.png",
    "app_data_files/papaya_flavour_wheel.png"
  ),
  appName = "papaya-flavour-predictor"
)
