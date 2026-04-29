#great ref for leave one (group) out cross-validation: https://users.aalto.fi/~ave/CV-FAQ.html#hierarchical

library(rstanarm)
library(bayesplot)
library(pROC)
library(scoring)
library(lme4)

prepDistancesForAnalysis <- function( distRf, distHdg, capForPrediction = TRUE ) {
  distRfInverseScale <- 7.212275
  minTransformedDistRf <- -1.3929
  distHdgInverseScale <- 11.26151
  minTransformedDistHdg <- -1.148933
  
  if (capForPrediction) {
    anDistRf <- ifelse( is.na(distRf) | distRf > 30, 30, distRf)
    anDistHdg <- ifelse( is.na(distHdg) | distHdg > 50, 50, distHdg)
  } else {
    anDistRf <- ifelse( is.na(distRf), max(distRf, na.rm = TRUE), distRf)
    anDistHdg <- ifelse( is.na(distHdg), max(distHdg, na.rm = TRUE), distHdg)
  }
  
  anDistRf <- anDistRf / distRfInverseScale + minTransformedDistRf
  anDistHdg <- anDistHdg / distHdgInverseScale + minTransformedDistHdg
  return(list( anDistRf = anDistRf, anDistHdg = anDistHdg ) )
}

cellingToAnalysisFormat <- function(input) {
  # Check if input is a character string (file name) or a data frame
  if(is.character(input)) {
    dF <- read.csv(input)
  } else if(is.data.frame(input)) {
    dF <- input
  } else {
    stop("Input must be either a file name or a data frame")
  }
  
  habitatNamesInFile <- c( "building", "lawn", "bushes", "trees", "cultivated.bed", "hedge", "artificial.surface", "water" )
  habitatNamesOutput <- c( "roof", "grass", "bush", "tree", "cult", "hedge", "artif", "water" )
    
  # Drop all cells with no parts inside boundary
  dF <- dF[ dF$boundary > 0,]

  # Make more convenient shorter garden ids
  dFS <- data.frame( id = sub( "-.*", "", dF$garden.guid ) )

  # Add each variable divided by boundary to dFS
  for (i in 1:length(habitatNamesInFile)) {
    dFS[[habitatNamesOutput[i]]] <- dF[[habitatNamesInFile[i]]] / dF$boundary
  }

  # Count cultivated as bush
  dFS$bushAndCult <- dFS$cult + dFS$bush
  
  # Add sparrows and distances
  dFS$sparrow <- ifelse( dF$sparrows.proportion > 0, 1, 0 )
  dFS$distRf <- dF$building.distance
  dFS$distHdg <- dF$hedge.distance
  distances <- prepDistancesForAnalysis( dFS$distRf, dFS$distHdg, capForPrediction = FALSE )
  dFS$anDistRf <- distances$anDistRf
  dFS$anDistHdg <- distances$anDistHdg
  
  return( dFS )
}

readGlasgowData <- function() {
  dFG <- read.csv("H:/My Drive/Projects/Environmental psychology/PRO-COAST/KLAMS - Central/Data/rspb20182911_si_002.csv")
  return( data.frame(
    batch = rep("glasgow",nrow(dFG)),
    id = dFG$Colony,
    sparrow = dFG$Sparrow,
    hedge = dFG$Hedge,
    grass = dFG$Grass,
    tree = dFG$Tree,
    roof = dFG$Roof,
    artif = dFG$Artificial.Surface,
    bushAndCult = dFG$Bush,
    anDistRf = dFG$DistRf,
    anDistHdg = dFG$DistHdg
  ) )
}

formStandardUncorr <- sparrow ~ hedge + grass + tree + roof + artif + bushAndCult + anDistHdg + anDistRf + ( 1 | id ) +
  ( 0 + hedge + grass + tree + roof + artif + bushAndCult + anDistHdg + anDistRf || id )


verbose_fit_check <- function( mod, dat ) {
  # Extract fixed effects only predictions
  # re.form=NA tells the function to ignore random effects
  posterior_probs_fixed_only <- posterior_epred(mod, 
                                               newdata = dat, offset = caseControlOffset,
                                               re.form = NA)  # This is a key part - ignore random effects

  # Average across posterior samples
  mean_predicted_probs_fixed <- colMeans(posterior_probs_fixed_only)

  # Calculate AUC with fixed effects only
  roc_obj_fixed <- roc(dat$sparrow, mean_predicted_probs_fixed)
  fixed_only_auc <- auc(roc_obj_fixed)

  # Compare with full model AUC
  posterior_probs_full <- posterior_epred(mod, newdata = dat, offset = caseControlOffset)
  mean_predicted_probs_full <- colMeans(posterior_probs_full)
  roc_obj_full <- roc(dat$sparrow, mean_predicted_probs_full)
  full_model_auc <- auc(roc_obj_full)

  # Print both AUCs
  print(paste("AUC with random effects:", round(full_model_auc, 3)))
  print(paste("AUC without random effects (fixed effects only):", round(fixed_only_auc, 3)))

  # Plot both ROC curves for comparison
  plot(roc_obj_full, main="ROC Curves: With vs. Without Random Effects", 
       col="blue", lwd=2)
  plot(roc_obj_fixed, add=TRUE, col="red", lwd=2)
  legend("bottomright", legend=c("With Random Effects", "Fixed Effects Only"), 
         col=c("blue", "red"), lwd=2)

  # Confusion matrix for fixed effects only model
  create_confusion_matrix <- function(probs, threshold) {
    predicted_classes <- ifelse(probs > threshold, 1, 0)
    conf_matrix <- table(Predicted = factor(predicted_classes, levels=c(0,1)), 
                        Observed = factor(dat$sparrow, levels=c(0,1)))
    
    # Calculate metrics
    sensitivity <- conf_matrix[2,2] / sum(conf_matrix[,2])
    specificity <- conf_matrix[1,1] / sum(conf_matrix[,1])
    accuracy <- (conf_matrix[1,1] + conf_matrix[2,2]) / sum(conf_matrix)
    
    return(list(
      confusion_matrix = conf_matrix,
      metrics = c(threshold = threshold,
                 sensitivity = sensitivity, 
                 specificity = specificity,
                 accuracy = accuracy)
    ))
  }

  # Standard threshold of 0.5
  fixed_results_0.5 <- create_confusion_matrix(mean_predicted_probs_fixed, 0.5)
  print("Confusion Matrix - Fixed Effects Only (threshold = 0.5):")
  print(fixed_results_0.5$confusion_matrix)
  print(round(fixed_results_0.5$metrics, 3))

  # Find optimal threshold for fixed effects model
  # Get the coordinates for best threshold (maximizing sensitivity + specificity)
  coords_result <- coords(roc_obj_fixed, "best")
  # Extract the first threshold value if multiple are returned
  optimal_threshold_fixed <- coords_result$threshold[1]
  print(paste("Optimal threshold (fixed effects only):", round(optimal_threshold_fixed, 3)))

  # Calculate confusion matrix with optimal threshold
  fixed_results_optimal <- create_confusion_matrix(mean_predicted_probs_fixed, optimal_threshold_fixed)
  print("Confusion Matrix - Fixed Effects Only (optimal threshold):")
  print(fixed_results_optimal$confusion_matrix)
  print(round(fixed_results_optimal$metrics, 3))

  # Optional: Calculate and print additional performance metrics
  # Calculate prediction improvement with optimal threshold vs. 0.5
  improvement <- fixed_results_optimal$metrics["accuracy"] - fixed_results_0.5$metrics["accuracy"]
  print(paste("Accuracy improvement with optimal threshold:", round(improvement, 3)))
}

merge_glasgow_with_REDS_data <- function( glas, reds, caseControlOffset ) {
  reds$batch = "new"
  glas$batch = "old"
  glas$prevalenceOffset <- caseControlOffset
  reds$prevalenceOffset <- 0
  return( data.frame(
    batch = c(glas$batch, reds$batch),
    id = c(glas$id, reds$id),
    prevalenceOffset = c(glas$prevalenceOffset, reds$prevalenceOffset),
    sparrow = c(glas$sparrow, reds$sparrow),
    roof = c(glas$roof, reds$roof),
    grass = c(glas$grass, reds$grass),
    bushAndCult = c(glas$bushAndCult, reds$bushAndCult),
    tree = c(glas$tree, reds$tree),
    hedge = c(glas$hedge, reds$hedge),
    artif = c(glas$artif, reds$artif),
    anDistRf = c(glas$anDistRf, reds$anDistRf),
    anDistHdg = c(glas$anDistHdg, reds$anDistHdg)
  ) )
}


model_prediction_AUC_confusion <- function(model, new_data, outcome_var = "sparrow", threshold = 0.5, offset = 0 ) {
  # Extract outcome vector
  true_outcomes <- new_data[[outcome_var]]
  
  # Get fixed-effects-only predictions (ignoring random effects)
  posterior_probs <- posterior_epred(model, newdata = new_data, offset = offset, re.form = NA)
  mean_predicted_probs <- colMeans(posterior_probs)

  # Calculate AUC
  roc_obj <- roc(true_outcomes, mean_predicted_probs, quiet=TRUE)
  model_auc <- auc(roc_obj)
  
  # Calculate confusion matrix
  predicted_classes <- ifelse(mean_predicted_probs > threshold, 1, 0)
  conf_matrix <- table(
    Predicted = factor(predicted_classes, levels = c(0, 1)),
    Observed = factor(true_outcomes, levels = c(0, 1))
  )
  
  # Print results (just AUC and confusion matrix as requested)
  cat(sprintf("Fixed-Effects-Only AUC: %.3f\n", model_auc))
  cat("Confusion Matrix (threshold =", threshold, "):\n")
  print(conf_matrix)
  
  cat("Don't forget, confusion matrix can be bad with perfect AUC because the threshold used for confusion may be poor.\n")
  
  # Return just the AUC
  return(model_auc)
}

multi_garden_prediction_AUC_confusion <- function(model, data, batch_col = "batch", id_col = "id", 
                           outcome_var = "sparrow", batch_value = "new",
                           threshold = 0.5, offset = 0 ) {
  # Filter to only new data
  new_data <- data[data[[batch_col]] == batch_value, ]
  
  # Get unique IDs
  unique_ids <- unique(new_data[[id_col]])
  
  # Initialize results vector
  results <- numeric(length(unique_ids))
  names(results) <- unique_ids
  
  # Loop through each ID
  for (i in seq_along(unique_ids)) {
    current_id <- unique_ids[i]
    id_data <- new_data[new_data[[id_col]] == current_id, ]
    
    # Count observations
    n_obs <- nrow(id_data)
    
    # Check if both classes are present
    if (length(unique(id_data[[outcome_var]])) < 2) {
      cat(sprintf("Skipping ID %s: Less than two classes present\n", current_id))
      results[i] = NA
      next
    }
    
    # Print ID header
    cat("\n========================================\n")
    cat(sprintf("Evaluating ID: %s (n = %d)\n", current_id, n_obs))
    cat("========================================\n")
    
    # Evaluate this ID and store AUC
    results[i] <- model_prediction_AUC_confusion(model, id_data, outcome_var, threshold, offset)
  }
  
  # Filter out NA values
  valid_results <- results[!is.na(results)]
  
  # Overall summary
  cat("\n========================================\n")
  cat("OVERALL SUMMARY\n")
  cat("========================================\n")
  cat(sprintf("Number of IDs evaluated: %d\n", length(valid_results)))
  cat(sprintf("Mean AUC across IDs: %.3f\n", mean(valid_results)))
  cat(sprintf("Median AUC across IDs: %.3f\n", median(valid_results)))
  cat(sprintf("Range of AUC values: %.3f - %.3f\n", 
              min(valid_results), max(valid_results)))
  
  # Return all results
  return(results)
}

runCrossValidationModels <- function(dat, form) {
  modOld <- stan_glmer(form, data = dat[dat$batch == "old", ], offset = prevalenceOffset, chains = 1, family = binomial)
  newIds <- unique(dat[dat$batch == "new", "id"])
  
  # Create a named list to store results
  results <- vector("list", length(newIds))
  names(results) <- newIds
  
  for (i in seq_along(newIds)) {
    currentId <- newIds[i]
    currentTrainData <- dat[dat$id != currentId, ]
    currentTestData <- dat[dat$id == currentId, ]
    nObs <- nrow(currentTestData)
    
    cat("\n========================================\n")
    cat(sprintf("Evaluating ID: %s (n = %d)\n", currentId, nObs ))
    cat("========================================\n")
    
    modCurrent <- stan_glmer(form, data = currentTrainData, offset = prevalenceOffset, chains = 1, family = binomial)
    
    currentTestData$newProbs <- colMeans(posterior_epred(modCurrent, newdata = currentTestData, re.form = NA, offset = caseControlOffset))
    currentTestData$oldProbs <- colMeans(posterior_epred(modOld, newdata = currentTestData, re.form = NA, offset = caseControlOffset))
    
    bSNew <- brierscore(sparrow ~ newProbs, currentTestData, decomp = TRUE)
    bSOld <- brierscore(sparrow ~ oldProbs, currentTestData, decomp = TRUE)
    
    results[[i]] <- list(mod = modCurrent, bSNew = bSNew, bSOld = bSOld, nObs <- nObs )
  }
  
  return(results)
}

assessCrossValidationPermutation <- function(res, n_perms = 100 ) {
  all_i <- vector()
  all_rawDiffScores <- vector()
  all_discrimOld <- vector()
  all_discrimNew <- vector()
  all_miscalOld <- vector()
  all_miscalNew <- vector()
  for(i in seq_along(res)) {
    rawDiffScores <- res[[i]]$bSNew$rawscores - res[[i]]$bSOld$rawscores
    all_i <- c(all_i, rep(i, length(rawDiffScores)))
    all_rawDiffScores <- c(all_rawDiffScores, rawDiffScores)
    all_discrimOld <- c( all_discrimOld, res[[i]]$bSOld$decomp$components["discrim",1] )
    all_discrimNew <- c( all_discrimNew, res[[i]]$bSNew$decomp$components["discrim",1] )
    all_miscalOld <- c( all_miscalOld, res[[i]]$bSOld$decomp$components["miscal",1] )
    all_miscalNew <- c( all_miscalNew, res[[i]]$bSNew$decomp$components["miscal",1] )
  }
  all_discrimDiffScores <- all_discrimNew - all_discrimOld
  all_miscalDiffScores <- all_miscalNew - all_miscalOld
  df <- data.frame(i = factor(all_i), rawDiffScores = all_rawDiffScores)
  observed_model <- lmer(rawDiffScores ~ 1 + (1|i), data = df)
  observed_intercept <- fixef(observed_model)[1]
  perm_intercepts <- numeric(n_perms)
  for(p in 1:n_perms) {
    perm_df <- df
    # For testing the intercept, flip signs of all observations randomly
    perm_df$rawDiffScores <- sample(c(-1, 1), nrow(df), replace = TRUE) * df$rawDiffScores
    suppressMessages({ #We get warnings because flipping signs may mean random effects go to near zero, but that's no issue for the intercept
      perm_model <- lmer(rawDiffScores ~ 1 + (1|i), data = perm_df)
    })
    perm_intercepts[p] <- fixef(perm_model)[1]
  }
  p_value <- mean(abs(perm_intercepts) >= abs(observed_intercept))
  return(list(
    fat = list( 
      data = df,
      observed_model = observed_model,
      permutation_distribution = perm_intercepts,
      bS = list(
        discrimOld = all_discrimOld,
        discrimNew = all_discrimNew,
        discrimDiffScores = all_discrimDiffScores,
        miscalOld = all_miscalOld,
        miscalNew = all_miscalNew,
        miscalDiffScores = all_miscalDiffScores
      )
    ),
    thin = list(
      observed_intercept = observed_intercept,
      permutation_p_value = p_value
    )
  ))
}
  
assessCrossValidationBayesian <- function(res, seed = 12345, chains = 4, iter = 2000, ...) {
  all_i <- vector()
  all_rawDiffScores <- vector()
  all_discrimOld <- vector()
  all_discrimNew <- vector()
  all_miscalOld <- vector()
  all_miscalNew <- vector()
  for(i in seq_along(res)) {
    rawDiffScores <- res[[i]]$bSNew$rawscores - res[[i]]$bSOld$rawscores
    all_i <- c(all_i, rep(i, length(rawDiffScores)))
    all_rawDiffScores <- c(all_rawDiffScores, rawDiffScores)
    all_discrimOld <- c(all_discrimOld, res[[i]]$bSOld$decomp$components["discrim",1])
    all_discrimNew <- c(all_discrimNew, res[[i]]$bSNew$decomp$components["discrim",1])
    all_miscalOld <- c(all_miscalOld, res[[i]]$bSOld$decomp$components["miscal",1])
    all_miscalNew <- c(all_miscalNew, res[[i]]$bSNew$decomp$components["miscal",1])
  }
  all_discrimDiffScores <- all_discrimNew - all_discrimOld
  all_miscalDiffScores <- all_miscalNew - all_miscalOld
  
  # Main data frame for raw score differences
  df <- data.frame(i = factor(all_i), rawDiffScores = all_rawDiffScores)
  
  # Data frames for component scores
  discrim_df <- data.frame(i = factor(1:length(all_discrimDiffScores)), 
                          discrimDiffScores = all_discrimDiffScores)
  miscal_df <- data.frame(i = factor(1:length(all_miscalDiffScores)), 
                         miscalDiffScores = all_miscalDiffScores)
  
  # Fit Bayesian mixed model for raw score differences
  raw_model <- stan_glmer( rawDiffScores ~ 1 + (1|i), data = df  )
  
  # Fit Bayesian models for component differences (no random effects needed)
  discrim_model <- stan_glm( discrimDiffScores ~ 1, data = discrim_df )
  miscal_model <- stan_glm( miscalDiffScores ~ 1, data = miscal_df )
  
  # Calculate posterior probabilities
  raw_samples <- as.matrix(raw_model)
  raw_intercept <- raw_samples[, "(Intercept)"]
  prob_raw_negative <- mean(raw_intercept < 0)
  
  discrim_samples <- as.matrix(discrim_model)
  discrim_intercept <- discrim_samples[, "(Intercept)"]
  prob_discrim_negative <- mean(discrim_intercept < 0)
  
  miscal_samples <- as.matrix(miscal_model)
  miscal_intercept <- miscal_samples[, "(Intercept)"]
  prob_miscal_negative <- mean(miscal_intercept < 0)
  
  return(list(
    fat = list(
      data = df,
      raw_model = raw_model,
      discrim_model = discrim_model,
      miscal_model = miscal_model,
      bS = list(
        discrimOld = all_discrimOld,
        discrimNew = all_discrimNew,
        discrimDiffScores = all_discrimDiffScores,
        miscalOld = all_miscalOld,
        miscalNew = all_miscalNew,
        miscalDiffScores = all_miscalDiffScores
      )
    ),
    thin = list(
      raw_intercept = fixef(raw_model)[1],
      discrim_intercept = coef(discrim_model)[1],
      miscal_intercept = coef(miscal_model)[1],
      raw_credible_interval = posterior_interval(raw_model, prob = 0.95)[1,],
      discrim_credible_interval = posterior_interval(discrim_model, prob = 0.95)[1,],
      miscal_credible_interval = posterior_interval(miscal_model, prob = 0.95)[1,],
      prob_raw_negative = prob_raw_negative,
      prob_discrim_negative = prob_discrim_negative,
      prob_miscal_negative = prob_miscal_negative
    )
  ))
}