#################################################################################################
# PREDICTION METRIC FUNCTIONS
#################################################################################################

# This is needed despite the pre-registered intention to use scoring::brierscore, because:
# 1. That function crashes on one of the gardens, and
# 2. By default it makes the decomposition groups by rounding to one d.p.,
# which is inappropriate for the probabilities here, which tend to be low as prevalence is low
brier_decomposition_own_rolled <- function(observed, predicted) {
  rawscores <- (predicted - observed)^2
  overall_mean <- mean(observed)

  # Create 10 equal-width bins
  pred_groups <- data.frame(
    predicted = predicted,
    observed = observed,
    bin = cut(predicted, breaks = 10, include.lowest = TRUE, labels = FALSE)
  ) %>%
    group_by(bin) %>%
    summarise(
      n = n(),
      mean_predicted = mean(predicted),
      mean_obs = mean(observed),
      .groups = 'drop'
    ) %>%
    filter(n > 0)  # Remove empty bins

  # Uncertainty
  uncertainty <- overall_mean * (1 - overall_mean)

  # Discrimination (Resolution)
  discrimination <- sum(pred_groups$n * (pred_groups$mean_obs - overall_mean)^2) / length(observed)

  # Miscalibration (Reliability)
  miscalibration <- sum(pred_groups$n * (pred_groups$mean_predicted - pred_groups$mean_obs)^2) / length(observed)

  # Verification: Brier score should equal uncertainty + miscalibration - discrimination
  brier_score <- mean(rawscores)
  decomp_brier <- uncertainty + miscalibration - discrimination

  list(
    rawscores = rawscores,
    brier_score = brier_score,
    decomp = list(
      components = matrix(
        c(discrimination, miscalibration, uncertainty),
        nrow = 3,
        dimnames = list(c("discrim", "miscal", "uncertainty"), "")
      ),
      decomp_sum = decomp_brier,
      matches = abs(brier_score - decomp_brier) < 1e-10
    )
  )
}

# Extract predictions from model (works with stan_glmer and INLA)
compute_predictions <- function(model, data, offset = 0, useRandom = FALSE, model_type = "stan_glmer", original_data = NULL) {
  if(model_type == "stan_glmer") {
    if(useRandom) {
      reform <- NULL
    } else {
      reform <- NA
    }
    predicted_probs_matrix <- posterior_epred(model, newdata = data, re.form = reform, offset = offset)
    predicted_probs <- colMeans(predicted_probs_matrix)
    return(predicted_probs)
  } else if(model_type == "inla") {
    if(is.null(original_data)) {
      stop("original_data must be provided when model_type is 'inla'")
    }
    pred_data <- data[, c("id", "prevalenceOffset", pred_names)]
    pred_data$sparrow <- NA
    if(!useRandom) {
      pred_data$id <- "not_in_training_data"
    }
    original_subset <- original_data[, c("id", "sparrow", "prevalenceOffset", pred_names)]
    combined_data <- rbind(original_subset, pred_data)
    new_model <- create_inla_model(combined_data)
    n_original <- nrow(original_data)
    n_pred <- nrow(pred_data)
    pred_indices <- (n_original + 1):(n_original + n_pred)
    predicted_probs <- new_model$summary.fitted.values$mean[pred_indices]
    return(predicted_probs)
  } else if(model_type == "inla_coef_formula") {
    fixed_effects <- model$summary.fixed$mean
    names(fixed_effects) <- rownames(model$summary.fixed)
    linear_pred <- fixed_effects["(Intercept)"] +
                   data$hedge * fixed_effects["hedge"] +
                   data$grass * fixed_effects["grass"] +
                   data$tree * fixed_effects["tree"] +
                   data$roof * fixed_effects["roof"] +
                   data$artif * fixed_effects["artif"] +
                   data$bushAndCult * fixed_effects["bushAndCult"] +
                   data$anDistHdg * fixed_effects["anDistHdg"] +
                   data$anDistRf * fixed_effects["anDistRf"] +
                   offset
    predicted_probs <- plogis(linear_pred)
    return(predicted_probs)
  } else {
    stop("model_type must be 'stan_glmer', 'inla', or 'inla_coef_formula'")
  }
}

# Calculate confusion matrix for each cluster
calculate_cluster_confusion_matrices <- function(data_with_pred) {
  cluster_ids <- unique(data_with_pred$id)
  all_matrices <- list()
  for(i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[i]
    cluster_data <- data_with_pred %>%
      filter(id == cluster_id)
    conf_matrix <- table(
      Predicted = factor(cluster_data$predicted_binary, levels=c(0,1)),
      Observed = factor(cluster_data$sparrow, levels=c(0,1))
    )
    confusion_matrix <- prop.table(conf_matrix) * 100
    dimnames(confusion_matrix) <- list(
      Predicted = c("No Sparrow", "Sparrow"),
      Observed = c("No Sparrow", "Sparrow")
    )
    all_matrices[[i]] <- confusion_matrix
  }
  return(all_matrices)
}

# Optimize threshold using the metric calculation function
optimize_threshold <- function(model, data, offset = 0, useRandom = FALSE,
                               metric = "youden", lower = 0.05, upper = 0.95, tol = 0.01,
                               model_type = "stan_glmer", original_data = NULL, prevalence = NULL
) {
  # Auto-add prevalenceOffset if missing
  if(!"prevalenceOffset" %in% names(data)) {
    data$prevalenceOffset <- offset
  }
  predicted_probs <- compute_predictions( model, data, offset, useRandom, model_type, original_data )
  data_with_pred <- data %>%
    mutate(predicted_probs = predicted_probs)

  objective_function <- function(threshold) {
    predicted_binary <- ifelse(data_with_pred$predicted_probs > threshold, 1, 0)
    temp_data <- data_with_pred %>%
      mutate(predicted_binary = predicted_binary)

    # Calculate cluster-level confusion matrices for all metrics
    all_matrices <- calculate_cluster_confusion_matrices(temp_data)
    mean_confusion_matrix <- Reduce("+", all_matrices) / length(all_matrices)

    # For prevalence metric, calculate absolute difference from target
    if(metric == "prevalence") {
      if(is.null(prevalence)) {
        stop("prevalence argument must be provided when metric = 'prevalence'")
      }
      predicted_prevalence <- sum(mean_confusion_matrix[2, ]) / 100  # Convert percentage to proportion
      value <- abs(predicted_prevalence - prevalence)
      rm(temp_data, all_matrices, mean_confusion_matrix)
      return(value)  # Minimize absolute difference
    }

    # For other metrics, use confusion matrix approach
    value <- calculate_confusion_metrics(mean_confusion_matrix, metric = metric)
    rm(temp_data, all_matrices, mean_confusion_matrix)
    return(-value)  # Minimize negative
  }

  result <- optimize(f = objective_function, interval = c(lower, upper), tol = tol)
  optimal_threshold <- result$minimum

  # Calculate final predicted prevalence for reporting
  predicted_binary <- ifelse(data_with_pred$predicted_probs > optimal_threshold, 1, 0)
  temp_data <- data_with_pred %>%
    mutate(predicted_binary = predicted_binary)
  all_matrices <- calculate_cluster_confusion_matrices(temp_data)
  mean_confusion_matrix <- Reduce("+", all_matrices) / length(all_matrices)
  predicted_prev <- sum(mean_confusion_matrix[2, ]) / 100  # Convert percentage to proportion

  if(metric == "prevalence") {
    optimal_value <- result$objective  # This is the absolute difference
    cat("\n=== OPTIMIZATION RESULTS ===\n")
    cat("Metric:", metric, "\n")
    cat("Target prevalence:", round(prevalence, 3), "\n")
    cat("Optimal threshold:", round(optimal_threshold, 3), "\n")
    cat("Predicted prevalence at optimal threshold:", round(predicted_prev, 3), "\n")
    cat("Absolute difference:", round(optimal_value, 3), "\n")
  } else {
    optimal_value <- -result$objective
    cat("\n=== OPTIMIZATION RESULTS ===\n")
    cat("Metric:", metric, "\n")
    cat("Optimal threshold:", round(optimal_threshold, 3), "\n")
    cat("Metric value at optimal threshold:", round(optimal_value, 3), "\n")
    cat("Predicted prevalence:", round(predicted_prev, 3), "\n")
  }

  return(optimal_threshold)
}

# Cluster as in garden or colony
generate_cluster_prediction_summaries <- function(model, data, threshold = 0.5, offset = 0,
                                                  useRandom = FALSE, outputEach = FALSE,
                                                  own_rolled_brier = TRUE, model_type = "stan_glmer",
                                                  original_data = NULL ) {
  # Auto-add prevalenceOffset if missing (for backwards compatibility)
  if(!"prevalenceOffset" %in% names(data)) {
    data$prevalenceOffset <- offset
  }

  predicted_probs <- compute_predictions(model, data, offset, useRandom, model_type, original_data)
  predicted_binary <- ifelse(predicted_probs > threshold, 1, 0)

  data_with_pred <- data %>%
    mutate(
      predicted_probs = predicted_probs,
      predicted_binary = predicted_binary
    )

  cluster_ids <- unique(data$id)
  cluster_statistics <- list()
  brier_raw_scores <- data.frame()
  brier_decomp_scores <- data.frame()
  cluster_mean_sparrow_probs <- numeric(length(cluster_ids))

  for(i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[i]
    cluster_data <- data_with_pred %>%
      filter(id == cluster_id)

    cluster_mean_sparrow_probs[i] <- mean(cluster_data$predicted_probs)

    n_per_class <- table(cluster_data$sparrow)
    has_both_classes <- length(n_per_class) == 2

    conf_matrix <- table(Predicted = factor(cluster_data$predicted_binary, levels=c(0,1)),
                        Observed = factor(cluster_data$sparrow, levels=c(0,1)))
    total <- sum(conf_matrix)
    confusion_matrix <- prop.table(conf_matrix) * 100
    dimnames(confusion_matrix) <- list(
      Predicted = c("No Sparrow", "Sparrow"),
      Observed = c("No Sparrow", "Sparrow")
    )

    cluster_accuracy <- calculate_confusion_metrics(confusion_matrix, metric = "accuracy")

    if(has_both_classes) {
      sink(nullfile())
      cluster_brier <- if(!own_rolled_brier) {
        brierscore(sparrow ~ predicted_probs, cluster_data, decomp = TRUE)
      } else {
        brier_decomposition_own_rolled(observed = cluster_data$sparrow, predicted = cluster_data$predicted_probs)
      }
      sink()

      discrimination <- cluster_brier$decomp$components["discrim", 1]
      miscalibration <- cluster_brier$decomp$components["miscal", 1]
      mean_brier <- mean(cluster_brier$rawscores)
    } else {
      basic_brier_scores <- (cluster_data$predicted_probs - cluster_data$sparrow)^2
      mean_brier <- mean(basic_brier_scores)
      discrimination <- NA
      miscalibration <- NA
      cluster_brier <- list(rawscores = basic_brier_scores)
    }

    cluster_raw_scores <- data.frame(cluster_id = cluster_id, brier_score = cluster_brier$rawscores)
    brier_raw_scores <- rbind(brier_raw_scores, cluster_raw_scores)

    cluster_decomp_scores <- data.frame(
      cluster_id = cluster_id,
      discrimination = discrimination,
      miscalibration = miscalibration
    )
    brier_decomp_scores <- rbind(brier_decomp_scores, cluster_decomp_scores)

    cluster_statistics[[paste0("Cluster_", cluster_id)]] <- list(
      cluster_id = cluster_id,
      n_observations = total,
      confusion_matrix = confusion_matrix,
      accuracy = cluster_accuracy,
      mean_brier_score = mean_brier,
      discrimination = discrimination,
      miscalibration = miscalibration,
      mean_sparrow_probability = mean(cluster_data$predicted_probs)
    )

    if( outputEach ) {
      cat("\n=== Cluster", cluster_id, "===\n")
      cat("Number of observations:", total, "\n")
      cat("Mean sparrow probability:", round(mean(cluster_data$predicted_probs), 3), "\n")
      cat("Confusion Matrix (%):\n")
      print(round(confusion_matrix, 1))
      cat("Accuracy:", round(cluster_accuracy, 1), "%\n")
      cat("Mean Brier Score:", round(mean_brier, 3), "\n")
      if(!is.na(discrimination)) {
        cat("Discrimination:", round(discrimination, 3), "\n")
        cat("Miscalibration:", round(miscalibration, 3), "\n")
      } else {
        cat("Discrimination: NA (single outcome class)\n")
        cat("Miscalibration: NA (single outcome class)\n")
      }
      cat("------------------------\n")
    }
  }

  # Calculate mean confusion matrix across all clusters
  all_matrices <- lapply(cluster_statistics, function(x) x$confusion_matrix)
  mean_confusion_matrix <- Reduce("+", all_matrices) / length(all_matrices)

  overall_brier_scores <- sapply(cluster_statistics, function(x) x$mean_brier_score)
  discrimination_scores <- sapply(cluster_statistics, function(x) x$discrimination)
  miscalibration_scores <- sapply(cluster_statistics, function(x) x$miscalibration)

  mean_overall_brier <- mean(overall_brier_scores, na.rm = TRUE)
  mean_discrimination <- mean(discrimination_scores, na.rm = TRUE)
  mean_miscalibration <- mean(miscalibration_scores, na.rm = TRUE)

  overall_mean_sparrow_prob <- mean(cluster_mean_sparrow_probs)
  overall_median_sparrow_prob <- median(cluster_mean_sparrow_probs)
  overall_sparrow_prob_iqr <- IQR(cluster_mean_sparrow_probs)

  plr <- calculate_confusion_metrics(mean_confusion_matrix, metric = "plr")

  # Build structured return object
  result <- list(
    threshold = threshold,
    clusters = cluster_statistics,
    brier_raw_scores = brier_raw_scores,
    brier_decomp_scores = brier_decomp_scores,
    mean_confusion_matrix = mean_confusion_matrix,
    overall_mean_brier_score = mean_overall_brier,
    overall_mean_discrimination = mean_discrimination,
    overall_mean_miscalibration = mean_miscalibration,
    overall_mean_sparrow_probability = overall_mean_sparrow_prob,
    overall_median_sparrow_probability = overall_median_sparrow_prob
  )

  cat("\n=== SUMMARY ===\n")
  cat("Mean Confusion Matrix (%):\n")
  print(round(mean_confusion_matrix, 1))
  cat("\nPositive Likelihood Ratio:", round(plr, 3), "\n")
  cat("\nMean Brier Scores Across All Clusters:\n")
  cat("Overall Brier Score:", round(mean_overall_brier, 3), "\n")
  cat("Discrimination:", round(mean_discrimination, 3), "\n")
  cat("Miscalibration:", round(mean_miscalibration, 3), "\n")

  invisible(result)
}

# Calculate metrics from a confusion matrix
calculate_confusion_metrics <- function(conf_matrix, metric = "youden") {
  # Extract confusion matrix values
  TN <- conf_matrix["No Sparrow", "No Sparrow"]
  FP <- conf_matrix["Sparrow", "No Sparrow"]
  FN <- conf_matrix["No Sparrow", "Sparrow"]
  TP <- conf_matrix["Sparrow", "Sparrow"]
  # Calculate base metrics
  TPR <- TP / (TP + FN)
  TNR <- TN / (TN + FP)
  FPR <- FP / (FP + TN)
  precision <- TP / (TP + FP)
  # Calculate requested metric
  value <- switch(metric,
    "youden" = TPR + TNR - 1,
    "f1" = {
      if (precision + TPR == 0) 0 else 2 * (precision * TPR) / (precision + TPR)
    },
    "mcc" = {
      denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
      if (denom == 0) 0 else (TP * TN - FP * FN) / denom
    },
    "plr" = {
      if (FPR == 0) 1e10 else TPR / FPR
    },
    "dor" = {
      if (FP == 0 || FN == 0) 1e10 else log((TP * TN) / (FP * FN))
    },
    "balanced_acc" = (TPR + TNR) / 2,
    "accuracy" = (TP + TN) / (TP + TN + FP + FN) * 100,
    stop("Unknown metric. Choose: youden, f1, mcc, plr, dor, balanced_acc, accuracy")
  )
  return(value)
}

##################################################################################################
# For comparing two posteriors
##################################################################################################

inla_compare_fixed_effects <- function(model1, model2, param_name,
                                       ci_level = 0.95, n_samples = 100000) {
  if (!(param_name %in% names(model1$marginals.fixed))) {
    stop(paste("Parameter", param_name, "not found in model1"))
  }
  if (!(param_name %in% names(model2$marginals.fixed))) {
    stop(paste("Parameter", param_name, "not found in model2"))
  }
  # Extract marginals for the parameter
  marginal1 <- model1$marginals.fixed[[param_name]]
  marginal2 <- model2$marginals.fixed[[param_name]]
  # Sample from the marginals
  samples1 <- inla.rmarginal(n_samples, marginal1)
  samples2 <- inla.rmarginal(n_samples, marginal2)
  # Calculate differences
  difference <- samples1 - samples2
  # Compute credible interval
  alpha <- 1 - ci_level
  ci_probs <- c(alpha/2, 1 - alpha/2)
  ci_difference <- quantile(difference, probs = ci_probs)
  mean_difference <- mean(difference)
  median_difference <- median(difference)
  # Calculate directional posterior probability
  if (median_difference > 0) {
    prob_direction <- mean(difference > 0)
    direction <- "positive"
  } else if (median_difference < 0) {
    prob_direction <- mean(difference < 0)
    direction <- "negative"
  } else {
    # Exactly zero median (very rare)
    prob_direction <- 0.5
    direction <- "zero"
  }
  result <- list(
    mean_difference = mean_difference,
    median_difference = median_difference,
    ci_lower = ci_difference[1],
    ci_upper = ci_difference[2],
    ci_level = ci_level,
    prob_direction = prob_direction,
    direction = direction,
    samples = difference
  )
  cat("Comparison of fixed effect:", param_name, "\n")
  cat("Mean difference:", mean_difference, "\n")
  cat("Median difference:", median_difference, "\n")
  cat(paste0(ci_level * 100, "% Credible Interval: ["),
      ci_difference[1], ",", ci_difference[2], "]\n")
  cat(paste0("Posterior probability difference is ", direction, ": ",
             round(prob_direction * 100, 2), "%\n"))
  return(invisible(result))
}
