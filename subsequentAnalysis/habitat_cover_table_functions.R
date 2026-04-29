##################################################################
# HABITAT SUMMARY TABLE FUNCTIONS
##################################################################

create_habitat_summary_table <- function(data, data_for_hedge_distance_mask=NULL) {
  habitat_vars <- c("sparrow", "roof", "grass", "bushAndCult", "tree", "hedge", "artif", "distRf", "distHdg")
  if(!is.null(data_for_hedge_distance_mask)) {
    n <- length(data_for_hedge_distance_mask)
    tail_indices <- (nrow(data) - n + 1):nrow(data)
    data$distHdg[tail_indices] <- ifelse(
      data_for_hedge_distance_mask,
      NA,
      data$distHdg[tail_indices]
    )
  }
  # Section 1: Garden Advice only (batch == "new")
  garden_advice <- data %>%
    filter(batch == "new")
  section1 <- calculate_garden_summaries(garden_advice, habitat_vars, "GA")
  # Section 2: Both batches, no sparrows, analyzed separately
  no_sparrows <- data %>%
    filter(sparrow == 0)
  # Garden Advice (new batch, no sparrows)
  ga_nosparrow <- no_sparrows %>%
    filter(batch == "new")
  section2_ga <- calculate_garden_summaries(ga_nosparrow, habitat_vars, "GA_NoSparr")
  # GHSP (old batch, no sparrows)
  ghsp_nosparrow <- no_sparrows %>%
    filter(batch == "old")
  section2_ghsp <- calculate_garden_summaries(ghsp_nosparrow, habitat_vars, "GHSP_NoSparr")
  # Combine into wide format
  result <- section1 %>%
    left_join(section2_ga, by = "Variable", suffix = c("", "_temp")) %>%
    left_join(section2_ghsp, by = "Variable")

  # Add a custom print method
  class(result) <- c("habitat_summary", class(result))

  return(result)
}

calculate_garden_summaries <- function(data, vars, prefix) {
  summaries <- lapply(vars, function(var) {
    # Determine if this is a distance variable
    is_distance <- var %in% c("distRf", "distHdg")
    # Calculate per-garden means
    garden_stats <- data %>%
      group_by(id) %>%
      summarise(
        garden_mean = mean(.data[[var]], na.rm = TRUE),
        has_nonzero = any(.data[[var]] != 0, na.rm = TRUE),
        has_nonNA = any(!is.na(.data[[var]])),
        .groups = "drop"
      )
    # Grand summaries: mean and median of garden means
    grand_mean <- mean(garden_stats$garden_mean, na.rm = TRUE)
    grand_median <- median(garden_stats$garden_mean, na.rm = TRUE)
    # Proportion of gardens: non-zero for most vars, non-NA for distance vars
    if (is_distance) {
      prop_gardens <- mean(garden_stats$has_nonNA, na.rm = TRUE)
    } else {
      prop_gardens <- mean(garden_stats$has_nonzero, na.rm = TRUE)
    }
    data.frame(
      Variable = var,
      Mean = grand_mean,
      Median = grand_median,
      PropGardens = prop_gardens
    )
  })
  result <- do.call(rbind, summaries)
  # Add prefix to column names (except Variable)
  names(result)[-1] <- paste0(prefix, "_", names(result)[-1])
  return(result)
}

# Custom print method for habitat_summary objects
print.habitat_summary <- function(x, ...) {
  # Create a copy for formatting
  formatted <- x

  # Get distance variable row indices (rows 8 and 9)
  distance_rows <- which(formatted$Variable %in% c("distRf", "distHdg"))
  cover_rows <- setdiff(1:nrow(formatted), distance_rows)

  # Format PropGardens columns as percentages (no decimals) for all rows
  propgarden_cols <- grep("PropGardens", names(formatted))
  for (col in propgarden_cols) {
    formatted[[col]] <- paste0(round(formatted[[col]] * 100), "%")
  }

  # Format Mean and Median columns
  mean_median_cols <- grep("Mean|Median", names(formatted))

  for (col in mean_median_cols) {
    # For cover rows (1-7): percentages with no decimals
    formatted[[col]][cover_rows] <- paste0(round(x[[col]][cover_rows] * 100), "%")
    # For distance rows (8-9): numbers with 1 decimal place
    formatted[[col]][distance_rows] <- sprintf("%.1f", x[[col]][distance_rows])
  }

  # Print the formatted data frame
  print(as.data.frame(formatted), row.names = FALSE)

  invisible(x)
}

invertPrepDistancesForAnalysis <- function(dF) {
  distRfInverseScale <- 7.212275
  minTransformedDistRf <- -1.3929
  distHdgInverseScale <- 11.26151
  minTransformedDistHdg <- -1.148933
  dF$distRf <- (dF$anDistRf - minTransformedDistRf) * distRfInverseScale
  dF$distHdg <- (dF$anDistHdg - minTransformedDistHdg) * distHdgInverseScale
  return(dF)
}

create_model_comparison_table <- function(model1, model2, filename,
                                          metrics = c("plr"),
                                          pct_digits = 1,
                                          fpr_tpr_digits = 1,
                                          plr_digits = 2,
                                          brier_digits = 2) {

  # Helper function to process a single model
  process_model <- function(model) {
    # Transpose confusion matrix so rows = observed, cols = predicted
    cm <- t(model$mean_confusion_matrix)

    # Extract values
    tn <- cm["No Sparrow", "No Sparrow"]  # True Negative
    fp <- cm["No Sparrow", "Sparrow"]      # False Positive
    fn <- cm["Sparrow", "No Sparrow"]      # False Negative
    tp <- cm["Sparrow", "Sparrow"]         # True Positive

    # Calculate rates
    fpr <- fp / (fp + tn)  # False Positive Rate
    tpr <- tp / (tp + fn)  # True Positive Rate (Sensitivity/Recall)

    # Calculate all requested metrics
    metric_values <- sapply(metrics, function(m) {
      calculate_confusion_metrics(model$mean_confusion_matrix, metric = m)
    })

    # Return results
    list(
      threshold = model$threshold,
      cm_values = c(tn, fp, fn, tp),
      fpr = fpr,
      tpr = tpr,
      metric_values = metric_values,
      brier = model$overall_mean_brier_score,
      discrimination = model$overall_mean_discrimination,
      miscalibration = model$overall_mean_miscalibration
    )
  }

  # Process both models
  m1 <- process_model(model1)
  m2 <- process_model(model2)

  # Create metric labels with nice names
  metric_labels <- sapply(metrics, function(m) {
    switch(m,
      "youden" = "Youden's J statistic",
      "f1" = "F1 score",
      "mcc" = "Matthews correlation coefficient (MCC)",
      "plr" = "Positive likelihood ratio (PLR=TPR/FPR)",
      "dor" = "Diagnostic odds ratio (log scale)",
      "balanced_acc" = "Balanced accuracy",
      "accuracy" = "Accuracy",
      toupper(m)  # fallback
    )
  })

  # Build the initial rows
  metrics_col <- c("", "", "Discrimination threshold",
                   "False positive rate (FPR)", "True positive rate (TPR)")
  observed_col <- c("No sparrow", "Sparrow", rep("", 3))

  m1_nosparrow <- c(
    sprintf(paste0("%.", pct_digits, "f%%"), m1$cm_values[1]),
    sprintf(paste0("%.", pct_digits, "f%%"), m1$cm_values[3]),
    "",
    "",
    ""
  )

  m1_sparrow <- c(
    sprintf(paste0("%.", pct_digits, "f%%"), m1$cm_values[2]),
    sprintf(paste0("%.", pct_digits, "f%%"), m1$cm_values[4]),
    sprintf(paste0("%.", brier_digits, "f"), m1$threshold),
    sprintf(paste0("%.", fpr_tpr_digits, "f%%"), m1$fpr * 100),
    sprintf(paste0("%.", fpr_tpr_digits, "f%%"), m1$tpr * 100)
  )

  m2_nosparrow <- c(
    sprintf(paste0("%.", pct_digits, "f%%"), m2$cm_values[1]),
    sprintf(paste0("%.", pct_digits, "f%%"), m2$cm_values[3]),
    "",
    "",
    ""
  )

  m2_sparrow <- c(
    sprintf(paste0("%.", pct_digits, "f%%"), m2$cm_values[2]),
    sprintf(paste0("%.", pct_digits, "f%%"), m2$cm_values[4]),
    sprintf(paste0("%.", brier_digits, "f"), m2$threshold),
    sprintf(paste0("%.", fpr_tpr_digits, "f%%"), m2$fpr * 100),
    sprintf(paste0("%.", fpr_tpr_digits, "f%%"), m2$tpr * 100)
  )

  # Add metric rows
  for (i in seq_along(metrics)) {
    metrics_col <- c(metrics_col, metric_labels[i])
    observed_col <- c(observed_col, "")
    m1_nosparrow <- c(m1_nosparrow, "")
    m1_sparrow <- c(m1_sparrow, sprintf(paste0("%.", plr_digits, "f"), m1$metric_values[i]))
    m2_nosparrow <- c(m2_nosparrow, "")
    m2_sparrow <- c(m2_sparrow, sprintf(paste0("%.", plr_digits, "f"), m2$metric_values[i]))
  }

  # Add Brier score rows
  metrics_col <- c(metrics_col,
                   "Overall Brier score (lower is better)",
                   "Discrimination decomposition (higher is better)",
                   "Miscalibration decomposition (lower is better)")
  observed_col <- c(observed_col, "", "", "")

  m1_nosparrow <- c(m1_nosparrow, "", "", "")
  m1_sparrow <- c(m1_sparrow,
                  sprintf(paste0("%.", brier_digits, "f"), m1$brier),
                  sprintf(paste0("%.", brier_digits, "f"), m1$discrimination),
                  sprintf(paste0("%.", brier_digits, "f"), m1$miscalibration))

  m2_nosparrow <- c(m2_nosparrow, "", "", "")
  m2_sparrow <- c(m2_sparrow,
                  sprintf(paste0("%.", brier_digits, "f"), m2$brier),
                  sprintf(paste0("%.", brier_digits, "f"), m2$discrimination),
                  sprintf(paste0("%.", brier_digits, "f"), m2$miscalibration))

  # Create the output data frame
  output <- data.frame(
    Metrics = metrics_col,
    Observed = observed_col,
    Model1_NoSparrow = m1_nosparrow,
    Model1_Sparrow = m1_sparrow,
    Model2_NoSparrow = m2_nosparrow,
    Model2_Sparrow = m2_sparrow,
    stringsAsFactors = FALSE
  )

  # Write to TSV
  write.table(output, file = file.path("threshold_statistics", filename ), sep = "\t", row.names = FALSE,
              quote = FALSE, col.names = c("Metrics", "Observed",
                                           "No sparrow", "Sparrow",
                                           "No sparrow", "Sparrow"))
  return(output)
}
