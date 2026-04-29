########################################################################################################
# MODEL MAKING FUNCTIONS, INCLUDING THE ONES FOR LEAVE-ONE-OUT CROSS VALIDATION
########################################################################################################

# Helper function to add ID columns for INLA
# (for some reason it needs a different random variable for each random slope, even though it's all the same cluster ID
add_id_columns <- function(data) {
  data$id.hedge <- data$id
  data$id.grass <- data$id
  data$id.tree <- data$id
  data$id.roof <- data$id
  data$id.artif <- data$id
  data$id.bushAndCult <- data$id
  data$id.anDistHdg <- data$id
  data$id.anDistRf <- data$id
  return(data)
}

# Helper function to create INLA formula
create_inla_formula <- function(prior_spec) {
  sparrow ~ hedge + grass + tree + roof + artif + bushAndCult +
    anDistHdg + anDistRf +
    offset(prevalenceOffset) +
    f(id, model = "iid", hyper = prior_spec) +
    f(id.hedge, hedge, model = "iid", hyper = prior_spec) +
    f(id.grass, grass, model = "iid", hyper = prior_spec) +
    f(id.tree, tree, model = "iid", hyper = prior_spec) +
    f(id.roof, roof, model = "iid", hyper = prior_spec) +
    f(id.artif, artif, model = "iid", hyper = prior_spec) +
    f(id.bushAndCult, bushAndCult, model = "iid", hyper = prior_spec) +
    f(id.anDistHdg, anDistHdg, model = "iid", hyper = prior_spec) +
    f(id.anDistRf, anDistRf, model = "iid", hyper = prior_spec)
}

# Helper function to fit INLA model
fit_inla_model <- function(formula, data) {
  inla(formula,
       data = data,
       family = "binomial",
       control.fixed = list(
         mean.intercept = 0, prec.intercept = 0.001,
         mean = 0, prec = 0.001
       ),
       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
       control.predictor = list(compute = TRUE, link = 1))
}

# Main function to create INLA models
create_inla_model <- function(data, u=10, alpha=0.5) {
  data <- add_id_columns(data)

  # PC prior: P(SD > U) = alpha
  prior_spec <- list(prec = list(prior = "pc.prec", param = c(u, alpha)))

  formINLA <- create_inla_formula(prior_spec)
  model <- fit_inla_model(formINLA, data)

  return(model)
}

# This function is simply a different path to the same result, kept around for testing purposes
create_inla_model_with_loco_legacy <- function(data, holdout_cluster_id, u=10, alpha=0.5, own_rolled_brier = TRUE) {
  if (!holdout_cluster_id %in% data$id) {
    stop(paste("Cluster ID", holdout_cluster_id, "not found in data"))
  }

  modified_data <- data
  holdout_mask <- data$id == holdout_cluster_id
  modified_data$sparrow[holdout_mask] <- NA

  modified_data <- add_id_columns(modified_data)

  # PC prior: P(SD > U) = alpha
  prior_spec <- list(prec = list(prior = "pc.prec", param = c(u, alpha)))

  formINLA <- create_inla_formula(prior_spec)
  model <- fit_inla_model(formINLA, modified_data)

  holdout_data <- data[holdout_mask, ]
  holdout_predictions <- model$summary.fitted.values[holdout_mask, ]
  temp_data <- data.frame(sparrow = holdout_data$sparrow, predicted_probs = holdout_predictions$mean)

  # Check if cluster has both outcome classes for decomposition
  sparrow_counts <- table(holdout_data$sparrow)
  has_both_classes <- length(sparrow_counts) == 2
  sink(nullfile())
  if(has_both_classes) {
    if( !own_rolled_brier ) {
      brier_result <- scoring::brierscore(sparrow ~ predicted_probs, temp_data, decomp = TRUE)
    } else {
      brier_result <- brier_decomposition_own_rolled( temp_data$sparrow, temp_data$predicted_probs )
    }
    brier_scores <- brier_result$rawscores
    discrimination <- brier_result$decomp$components["discrim", 1]
    miscalibration <- brier_result$decomp$components["miscal", 1]
  } else {
    brier_scores <- scoring::brierscore(sparrow ~ predicted_probs, temp_data)
    discrimination <- NA
    miscalibration <- NA
  }
  sink()
  return(list(
    model = model,
    holdout_data = holdout_data,
    holdout_predictions = holdout_predictions,
    holdout_cluster_id = holdout_cluster_id,
    brier_scores = brier_scores,
    mean_brier_score = mean(brier_scores),
    discrimination = discrimination,
    miscalibration = miscalibration
  ))
}

create_inla_model_with_loco <- function(data, holdout_cluster_id, u=10, alpha=0.5,
                                        own_rolled_brier = TRUE, prediction_method = "coef_formula"
) {
  if (!holdout_cluster_id %in% data$id) {
    stop(paste("Cluster ID", holdout_cluster_id, "not found in data"))
  }
  if (!prediction_method %in% c("holdout", "coef_formula")) {
    stop("prediction_method must be either 'holdout' or 'coef_formula'")
  }
  holdout_mask <- data$id == holdout_cluster_id
  if (prediction_method == "coef_formula") {
    training_data <- data[!holdout_mask, ]
    model <- create_inla_model(training_data, u = u, alpha = alpha)
    holdout_data <- data[holdout_mask, ]
    holdout_predictions_mean <- compute_predictions(
      model = model,
      data = holdout_data,
      offset = holdout_data$prevalenceOffset,
      model_type = "inla_coef_formula"
    )
    holdout_predictions <- data.frame(mean = holdout_predictions_mean)
  } else {  # prediction_method == "holdout" (fitting and predicting at the same time in that weird INLA way)
    modified_data <- data
    modified_data$sparrow[holdout_mask] <- NA
    modified_data <- add_id_columns(modified_data)
    prior_spec <- list(prec = list(prior = "pc.prec", param = c(u, alpha))) # PC prior: P(SD > U) = alpha
    formINLA <- create_inla_formula(prior_spec)
    model <- fit_inla_model(formINLA, modified_data)
    holdout_data <- data[holdout_mask, ]
    holdout_predictions <- model$summary.fitted.values[holdout_mask, ]
  }
  temp_data <- data.frame(sparrow = holdout_data$sparrow, predicted_probs = holdout_predictions$mean)
  sparrow_counts <- table(holdout_data$sparrow)
  has_both_classes <- length(sparrow_counts) == 2
  sink(nullfile())
    if(has_both_classes) {
      if( !own_rolled_brier ) {
        brier_result <- scoring::brierscore(sparrow ~ predicted_probs, temp_data, decomp = TRUE)
      } else {
        brier_result <- brier_decomposition_own_rolled( temp_data$sparrow, temp_data$predicted_probs )
      }
      brier_scores <- brier_result$rawscores
      discrimination <- brier_result$decomp$components["discrim", 1]
      miscalibration <- brier_result$decomp$components["miscal", 1]
    } else {
      brier_scores <- scoring::brierscore(sparrow ~ predicted_probs, temp_data)
      discrimination <- NA
      miscalibration <- NA
    }
  sink()
  return(list(
    model = model,
    holdout_data = holdout_data,
    holdout_predictions = holdout_predictions,
    holdout_cluster_id = holdout_cluster_id,
    brier_scores = brier_scores,
    mean_brier_score = mean(brier_scores),
    discrimination = discrimination,
    miscalibration = miscalibration
  ))
}

record_garden_results <- function(garden, baseline_results, enhanced_results,
                                   baseline_col_name, enhanced_col_name,
                                   method_dir) {
  # Create subdirectory for individual model files
  models_dir <- file.path(method_dir, "models")
  if (!dir.exists(models_dir)) {
    dir.create(models_dir, recursive = TRUE)
  }

  # Save the models separately
  baseline_model_file <- file.path(models_dir, paste0(garden, "_baseline.rds"))
  enhanced_model_file <- file.path(models_dir, paste0(garden, "_enhanced.rds"))
  saveRDS(baseline_results$model, baseline_model_file)
  saveRDS(enhanced_results$model, enhanced_model_file)

  n_obs <- length(baseline_results$brier_scores)

  garden_brier_data <- data.frame(
    garden = rep(garden, n_obs),
    baseline = baseline_results$brier_scores,
    enhanced = enhanced_results$brier_scores
  )
  names(garden_brier_data)[2:3] <- c(baseline_col_name, enhanced_col_name)

  garden_decomp_data <- data.frame(
    garden = garden,
    baseline_discrimination = baseline_results$discrimination,
    enhanced_discrimination = enhanced_results$discrimination,
    baseline_miscalibration = baseline_results$miscalibration,
    enhanced_miscalibration = enhanced_results$miscalibration
  )
  names(garden_decomp_data)[2:5] <- c(
    paste0(baseline_col_name, "_discrimination"),
    paste0(enhanced_col_name, "_discrimination"),
    paste0(baseline_col_name, "_miscalibration"),
    paste0(enhanced_col_name, "_miscalibration")
  )

  # Return everything except the models, but include predictions
  return(list(
    brier_scores = garden_brier_data,
    decomp_scores = garden_decomp_data,
    holdout_predictions = list(
      baseline = baseline_results$holdout_predictions,
      enhanced = enhanced_results$holdout_predictions
    ),
    model_files = list(
      baseline = baseline_model_file,
      enhanced = enhanced_model_file
    )
  ))
}

# for prediction_method:
# "holdout" is using the standard INLA mechanisms for making predictions
# "legacy" is also doing that, but with a different code pathway I was using at one time
# "coef_formula" is using a different method, predicting by constructing the formula from extracted coefficients
# (which is the method used for generating predictions in-app, but that is still based on MCMC results,
# which work well for that).
# I don't use coef_formula because it gives different and less plausible results for reasons I have not yet had time
# to investigate (pretty much the same pattern of results just with overall lower Brier scores)
run_loco_comparison <- function(dFO, dF, caseControlOffset, evaluate_all_data = FALSE, output_dir = "loco_results", own_rolled_brier = TRUE, prediction_method = "coef_formula" ) {
  # Create method-specific subdirectory
  method_dir <- file.path(output_dir, prediction_method)
  if (!dir.exists(method_dir)) {
    dir.create(method_dir, recursive = TRUE)
  }
  brier_file <- file.path(method_dir, "brier_scores.rds")
  decomp_file <- file.path(method_dir, "decomp_scores.rds")
  predictions_file <- file.path(method_dir, "holdout_predictions.rds")

  if (file.exists(brier_file) && file.exists(decomp_file)) {
    all_brier_scores <- readRDS(brier_file)
    decomp_scores <- readRDS(decomp_file)
    all_holdout_predictions <- if(file.exists(predictions_file)) readRDS(predictions_file) else list()
    if (nrow(all_brier_scores) > 0) {
      processed_gardens <- unique(all_brier_scores$garden)
      cat("Found existing results with", length(processed_gardens), "gardens already processed\n")
    } else {
      processed_gardens <- character(0)
      cat("Found empty result files - starting from beginning\n")
    }
  } else {
    all_brier_scores <- data.frame()
    decomp_scores <- data.frame()
    all_holdout_predictions <- list()
    processed_gardens <- character(0)
    cat("No existing results found - starting fresh\n")
  }

  dFBoth <- merge_glasgow_with_REDS_data(dFO, dF, caseControlOffset)
  garden_ids <- if (evaluate_all_data) {
    unique(c(dFO$id, dF$id))
  } else {
    unique(dF$id)
  }
  remaining_gardens <- setdiff(garden_ids, processed_gardens)
  if (length(remaining_gardens) == 0) {
    cat("All gardens already processed!\n")
    return(list(brier_scores = all_brier_scores,
                decomp_scores = decomp_scores,
                holdout_predictions = all_holdout_predictions))
  }
  cat("Processing", length(remaining_gardens), "remaining gardens out of", length(garden_ids), "total\n")

  for(i in seq_along(remaining_gardens)) {
    garden <- remaining_gardens[i]
    cat("Processing garden:", garden, "(", i, "of", length(remaining_gardens), "remaining)\n")
    if (evaluate_all_data && garden %in% dFO$id) {
      baseline_data <- merge_glasgow_with_REDS_data(dFO, data.frame(), caseControlOffset)
    } else {
      baseline_data <- merge_glasgow_with_REDS_data(dFO, dF[dF$id == garden, ], caseControlOffset)
    }
    if( prediction_method == "legacy" ) {
      baseline_results <- create_inla_model_with_loco_legacy(baseline_data, garden, own_rolled_brier = own_rolled_brier )
      enhanced_results <- create_inla_model_with_loco_legacy(dFBoth, garden, own_rolled_brier = own_rolled_brier )
    } else {
      baseline_results <- create_inla_model_with_loco(baseline_data, garden, own_rolled_brier = own_rolled_brier, prediction_method = prediction_method )
      enhanced_results <- create_inla_model_with_loco(dFBoth, garden, own_rolled_brier = own_rolled_brier, prediction_method = prediction_method )
    }
    results <- record_garden_results(garden, baseline_results, enhanced_results,
                                      "old_only", "all_data", method_dir)

    # Accumulate results
    all_brier_scores <- rbind(all_brier_scores, results$brier_scores)
    decomp_scores <- rbind(decomp_scores, results$decomp_scores)
    all_holdout_predictions[[garden]] <- results$holdout_predictions

    # Save all result structures (models are already saved separately)
    saveRDS(all_brier_scores, brier_file)
    saveRDS(decomp_scores, decomp_file)
    saveRDS(all_holdout_predictions, predictions_file)

    gc()
  }
  cat("Processing complete! Total gardens processed:", length(unique(all_brier_scores$garden)), "\n")
  return(list(brier_scores = all_brier_scores,
              decomp_scores = decomp_scores,
              holdout_predictions = all_holdout_predictions))
}

analyze_brier_comparison <- function(brier_data, decomp_data, run_quickly = TRUE ) {
  brierRaws <- brier_data
  brierRaws$improvement <- brierRaws$old_only - brierRaws$all_data
  brierRawMeans <- aggregate(improvement ~ garden, data = brierRaws, FUN = mean)
  hist(brierRawMeans$improvement)
  if( !run_quickly ) {
    modBrierRaws <- stan_glmer( improvement ~ 1 + (1|garden), data = brierRaws )
    posterior_interval(modBrierRaws, pars = "(Intercept)", prob = 0.95)
    mcmc_areas(modBrierRaws, pars = "(Intercept)", prob = 0.95)
  } else {
    # Run quick model on one data point per garden (same as decompositions)
    modBrierRaws <- stan_glm(improvement ~ 1, data = brierRawMeans, family = gaussian() )
  }
  decomp_data$discrimination_improvement <- decomp_data$all_data_discrimination - decomp_data$old_only_discrimination
  decomp_data$miscalibration_improvement <- decomp_data$old_only_miscalibration - decomp_data$all_data_miscalibration
  hist(decomp_data$discrimination_improvement)
  hist(decomp_data$miscalibration_improvement)
  modDiscrimination <- stan_glm(discrimination_improvement ~ 1, data = decomp_data, family = gaussian() )
  modMiscalibration <- stan_glm(miscalibration_improvement ~ 1, data = decomp_data, family = gaussian() )
  return(list(
    brierRaws = brierRaws,
    brierRawMeans = brierRawMeans,
    modBrierRaws = modBrierRaws,
    decompRaws = decomp_data,
    modDiscrimination = modDiscrimination,
    modMiscalibration = modMiscalibration
  ))
}

########################################################################################################
# FIXED EFFECTS PROBABILITY DISTRIBUTION PLOTS. CAN DO A MORE THAN THE SINGLE PLOT IN THE MANUSCRIPT
# e.g., can also plot priors
########################################################################################################

# Function to create multi-panel comparison plots for Bayesian models
plot_model_comparison <- function(models,
                                model_types,
                                distribution_types = NULL,
                                model_names = NULL,
                                line_thickness = NULL,
                                line_colors = NULL,
                                fill_colors = NULL,
                                fill_alpha = NULL,
                                title = "Model Comparison",
                                subtitle = NULL,
                                parNames = c("(Intercept)", "hedge", "grass", "tree", "roof", "artif", "bushAndCult", "anDistHdg", "anDistRf"),
                                save_plot = NULL) {

  # Hard-coded panel order and labels
  panel_order <- c("(Intercept)", "anDistHdg", "anDistRf", "tree", "hedge", "roof", "bushAndCult", "grass", "artif")
  panel_labels <- c("(Intercept)" = "Intercept",
                   "anDistHdg" = "Hedge proximity",
                   "anDistRf" = "Roof proximity",
                   "tree" = "Tree",
                   "hedge" = "Hedge",
                   "bushAndCult" = "Bush",
                   "grass" = "Grass",
                   "roof" = "Roof",
                   "artif" = "Artificial")

  # Validate inputs
  n_models <- length(models)
  if(length(model_types) != n_models) {
    stop("model_types must have same length as models")
  }

  # Set defaults
  if(is.null(distribution_types)) distribution_types <- rep("posterior", n_models)
  if(is.null(model_names)) model_names <- paste("Model", 1:n_models)
  if(is.null(line_thickness)) line_thickness <- rep(1.2, n_models)
  if(is.null(line_colors)) line_colors <- rainbow(n_models)
  if(is.null(fill_colors)) fill_colors <- rainbow(n_models)
  if(is.null(fill_alpha)) fill_alpha <- rep(0.7, n_models)

  # Initialize list to store data for each model
  all_data <- list()

  for(i in 1:n_models) {
    model <- models[[i]]
    model_type <- model_types[i]
    dist_type <- distribution_types[i]
    model_name <- model_names[i]

    if(model_type == "MCMC") {
      if(dist_type == "posterior") {
        # Extract MCMC posterior samples
        post_samples <- as.matrix(model, pars = parNames)

        # Convert to long format
        model_data <- as.data.frame(post_samples) %>%
          mutate(model = model_name) %>%
          pivot_longer(-model, names_to = "parameter", values_to = "value") %>%
          mutate(
            # Reverse sign for distance variables to represent proximity
            value = ifelse(parameter %in% c("anDistHdg", "anDistRf"), -value, value),
            model_type = "MCMC",
            line_thickness = line_thickness[i],
            line_color = line_colors[i],
            fill_color = fill_colors[i],
            fill_alpha = fill_alpha[i],
            x = NA_real_,  # Add consistent columns
            y = NA_real_
          ) %>%
          mutate(display_parameter = panel_labels[parameter],
                 display_parameter = factor(display_parameter, levels = panel_labels[panel_order]))

      } else if(dist_type == "prior") {
        # Extract analytical prior distributions
        prior_info <- prior_summary(model)
        marginals_data <- list()

        for(j in 1:length(parNames)) {
          param <- parNames[j]

          # Get prior parameters
          if(param == "(Intercept)") {
            # Intercept: no adjustment, use original scale
            prior_mean <- prior_info$prior_intercept$location
            prior_scale <- prior_info$prior_intercept$scale
          } else {
            # Coefficients: use adjusted_scale, location is always 0
            fixef_pos <- which(names(fixef(model)) == param)
            coef_pos <- fixef_pos - 1  # Remove intercept position

            prior_mean <- prior_info$prior$location[coef_pos]  # Always 0
            prior_scale <- prior_info$prior$adjusted_scale[coef_pos]
          }

          # Create grid of x values and compute analytical density
          x_range <- c(prior_mean - 4*prior_scale, prior_mean + 4*prior_scale)
          x_vals <- seq(x_range[1], x_range[2], length.out = 200)
          y_vals <- dnorm(x_vals, mean = prior_mean, sd = prior_scale)

          # Reverse sign for distance variables to represent proximity
          if(param %in% c("anDistHdg", "anDistRf")) {
            x_vals <- -x_vals
            # Reverse order to maintain proper density alignment
            x_vals <- rev(x_vals)
            y_vals <- rev(y_vals)
          }

          param_data <- data.frame(
            parameter = param,
            model = model_name,
            model_type = "MCMC",
            line_thickness = line_thickness[i],
            line_color = line_colors[i],
            fill_color = fill_colors[i],
            fill_alpha = fill_alpha[i],
            x = x_vals,
            y = y_vals,
            value = NA_real_  # Add consistent columns
          ) %>%
          mutate(display_parameter = panel_labels[parameter],
                 display_parameter = factor(display_parameter, levels = panel_labels[panel_order]))

          marginals_data[[param]] <- param_data
        }

        model_data <- do.call(rbind, marginals_data)
      }

    } else if(model_type == "INLA") {
      # Extract analytical marginals (INLA only does posteriors)
      marginals_data <- list()

      for(param in parNames) {
        # Handle potential intercept name mismatch
        inla_param_name <- param
        if(param == "(Intercept)" && !param %in% names(model$marginals.fixed)) {
          if("X.Intercept." %in% names(model$marginals.fixed)) {
            inla_param_name <- "X.Intercept."
          }
        }

        marginal <- model$marginals.fixed[[inla_param_name]]

        # Reverse sign for distance variables to represent proximity
        x_vals <- marginal[, 1]
        y_vals <- marginal[, 2]
        if(param %in% c("anDistHdg", "anDistRf")) {
          x_vals <- -x_vals
          # Reverse order to maintain proper density alignment
          order_idx <- order(x_vals)
          x_vals <- x_vals[order_idx]
          y_vals <- y_vals[order_idx]
        }

        param_data <- data.frame(
          parameter = param,
          model = model_name,
          model_type = "INLA",
          line_thickness = line_thickness[i],
          line_color = line_colors[i],
          fill_color = fill_colors[i],
          fill_alpha = fill_alpha[i],
          x = x_vals,
          y = y_vals,
          value = NA_real_  # Add consistent columns
        ) %>%
        mutate(display_parameter = panel_labels[parameter],
               display_parameter = factor(display_parameter, levels = panel_labels[panel_order]))

        marginals_data[[param]] <- param_data
      }

      model_data <- do.call(rbind, marginals_data)
    }

    all_data[[i]] <- model_data
  }

  # Combine all data
  combined_data <- do.call(rbind, all_data)

  # Calculate x-axis limits based only on posterior data (ignore priors)
  posterior_data <- combined_data[combined_data$model_type == "MCMC" &
                                  !is.na(combined_data$value), ]

  xlim_scales <- NULL
  if(nrow(posterior_data) > 0) {
    # Calculate limits per parameter based on posterior data only
    xlims <- posterior_data %>%
      group_by(display_parameter) %>%
      summarise(
        xmin = quantile(value, 0.001, na.rm = TRUE),
        xmax = quantile(value, 0.999, na.rm = TRUE),
        .groups = 'drop'
      )

    # Create a list of scale_x_continuous for each panel
    xlim_scales <- list()
    for(i in 1:nrow(xlims)) {
      param_name <- xlims$display_parameter[i]
      xlim_scales[[i]] <- as.formula(paste0("display_parameter == '", param_name, "' ~ scale_x_continuous(limits = c(", xlims$xmin[i], ", ", xlims$xmax[i], "))"))
    }
  }

  # Determine if we have mixed model types
  has_mcmc <- any(combined_data$model_type == "MCMC")
  has_inla <- any(combined_data$model_type == "INLA")
  mixed_types <- has_mcmc && has_inla

  # Create base plot
  p <- ggplot() +
    geom_vline(xintercept = 0, color = "black", linetype = "solid", alpha = 0.7) +
    geom_hline(yintercept = 0, color = "black", linetype = "solid", alpha = 0.7) +
    facet_wrap(~display_parameter, scales = "free") +
    theme_minimal() +
    labs(title = title,
         subtitle = subtitle,
         x = "Posterior parameter (log odds per variable unit)",
         y = ifelse(mixed_types || has_inla, "Density", "Density")) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),  # Remove legend title
          legend.justification.bottom = "left",
          legend.location = "plot",
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.5, "cm"),
          plot.background = element_rect(fill = "white", color = NA),  # White background
          panel.background = element_rect(fill = "white", color = NA),  # White panel background
          text = element_text(size = 14),  # 1.3x base text size (from ~11 to 14.3)
          axis.title = element_text(size = 11),  # Axis titles
          axis.text = element_text(size = 10),  # Axis tick labels
          plot.title = element_text(size = 16),  # Plot title
          plot.subtitle = element_text(size = 11),  # Plot subtitle
          legend.text = element_text(size = 10),  # Legend text
          strip.text = element_text(size = 11))  # Facet panel labels

  # Apply per-panel x-limits if we have posterior data
  if(!is.null(xlim_scales)) {
    p <- p + ggh4x::facetted_pos_scales(x = xlim_scales)
  }

  # Add layers for each model in order
  for(i in 1:n_models) {
    model_data <- all_data[[i]]
    model_name <- model_names[i]
    dist_type <- distribution_types[i]

    if(model_types[i] == "MCMC" && dist_type == "posterior") {
      # Check if line_color is NA for this model
      if(is.na(line_colors[i])) {
        # Add density with fill only (no outline)
        p <- p +
          geom_density(data = model_data,
                      aes(x = value, fill = model),
                      alpha = fill_alpha[i],
                      color = NA,  # No outline
                      linewidth = 0)
      } else {
        # Add density with both fill and line
        p <- p +
          geom_density(data = model_data,
                      aes(x = value, fill = model),
                      alpha = fill_alpha[i],
                      color = line_colors[i],
                      linewidth = line_thickness[i])
      }
    } else if((model_types[i] == "MCMC" && dist_type == "prior") || model_types[i] == "INLA") {
      # Add line for analytical curves (MCMC priors or INLA posteriors)
      p <- p +
        geom_line(data = model_data,
                 aes(x = x, y = y, color = model),
                 linewidth = line_thickness[i])
    }
  }

  # Set colors
  mcmc_posterior_models <- model_names[model_types == "MCMC" & distribution_types == "posterior"]
  analytical_models <- model_names[(model_types == "MCMC" & distribution_types == "prior") | model_types == "INLA"]

  if(length(mcmc_posterior_models) > 0) {
    p <- p + scale_fill_manual(name = "MCMC",  # Add title
                              values = setNames(fill_colors[model_types == "MCMC" & distribution_types == "posterior"],
                                                mcmc_posterior_models))
  }
  if(length(analytical_models) > 0) {
    analytical_indices <- which((model_types == "MCMC" & distribution_types == "prior") | model_types == "INLA")
    p <- p + scale_color_manual(name = "INLA",  # Add title
                               values = setNames(line_colors[analytical_indices],
                                                 analytical_models))
  }

  # Add to theme section to position titles above legend items
  p <- p + theme(
    legend.title = element_text(size=11, hjust = 0.5),  # Center the titles
    legend.title.position = "top",  # Position title above legend items
    legend.justification.bottom = "center",  # Center the legend box itself
    legend.box.just = "center",  # Center the legend box
  )

  # Save plot if filename provided
  if(!is.null(save_plot)) {
    # A4 width in inches: 8.27 inches (210mm)
    plot_width <- 6.27
    plot_height <- 6

    ggsave(filename = save_plot,
           plot = p,
           width = plot_width,
           height = plot_height,
           units = "in",
           dpi = 300,
           bg = "white")

    message(paste("Plot saved to:", save_plot))
  }

  return(p)
}
