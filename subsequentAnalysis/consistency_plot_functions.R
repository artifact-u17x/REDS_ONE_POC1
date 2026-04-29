# ============================================================================
# CONSISTENCY PLOT SYSTEM - note that internally this is called reliability,
# which was the original name
# ============================================================================

# Define patch type colors at top level
patch_colours <- c(
  "hedge" = "#FFA500",        # Hedge
  "grass" = "#00FF00",        # Grass (lawn)
  "tree" = "#8000FF",         # Trees
  "roof" = "#FF0000",         # Roof (building)
  "artif" = "#FF00FF",        # Artificial surface
  "bushAndCult" = "#FFFF00",  # Bush (bushes/cultivated bed)
  "anDistHdg" = "#FFA500",    # Hedge proximity (same as hedge)
  "anDistRf" = "#FF0000"      # Roof proximity (same as roof/building)
)

# Display name mappings
display_name_map <- c(
  "artif" = "Artificial (A)",
  "bushAndCult" = "Bush (B)",
  "grass" = "Grass (G)",
  "hedge" = "Hedge (H)",
  "anDistHdg" = "Hedge\nproximity (HP)",
  "anDistRf" = "Roof\nproximity (RP)",
  "roof" = "Roof (R)",
  "tree" = "Tree (T)"
)

label_map <- c(
  "artif" = "A",
  "bushAndCult" = "B",
  "grass" = "G",
  "hedge" = "H",
  "anDistHdg" = "HP",
  "anDistRf" = "RP",
  "roof" = "R",
  "tree" = "T"
)

# Calculate effect reliability data from model
calculate_effect_reliability <- function(model,
                                        data = NULL,
                                        use_empirical_sds = FALSE,
                                        id_col = "id",
                                        variance_threshold = 1e-10) {

  # Validate empirical SD requirements
  if (use_empirical_sds && is.null(data)) {
    stop("data must be provided when use_empirical_sds = TRUE")
  }

  # Extract posterior draws
  posterior_draws <- as.matrix(model)

  # Get fixed effects (excluding intercept)
  fixed_effects <- c("hedge", "grass", "tree", "roof", "artif", "bushAndCult", "anDistHdg", "anDistRf")
  fixed_posteriors <- posterior_draws[, fixed_effects]

  # Get random slope SDs - either empirical or from hyperparameters
  if (use_empirical_sds) {
    # Calculate unbiased empirical SDs
    empirical_results <- calculate_unbiased_random_sds(model, data, fixed_effects,
                                                       id_col, variance_threshold)

    # Extract empirical SD posteriors and organize into matrix
    random_slope_posteriors <- sapply(fixed_effects, function(var) {
      empirical_results[[var]]$empirical_sd_full
    })

    colnames(random_slope_posteriors) <- fixed_effects

  } else {
    # Use hyperparameters (original behavior)
    # Get random slope SDs (excluding the intercept)
    sigma_cols <- grep("Sigma\\[id:(?!\\(Intercept\\)).*,.*\\]",
                       colnames(posterior_draws),
                       value = TRUE,
                       perl = TRUE)
    random_slope_posteriors <- posterior_draws[, sigma_cols]

    # Clean up names
    colnames(random_slope_posteriors) <- gsub("Sigma\\[id:", "", colnames(random_slope_posteriors))
    colnames(random_slope_posteriors) <- gsub(",.*\\]", "", colnames(random_slope_posteriors))
    random_slope_posteriors <- random_slope_posteriors[, fixed_effects]
  }

  # Calculate pooled within-cluster SDs for standardization
  within_sds <- model.frame(model) %>%
    group_by(id) %>%
    summarise(across(c(hedge, grass, tree, roof, artif, bushAndCult, anDistHdg, anDistRf),
                     ~sd(., na.rm=TRUE), .names = "sd_{.col}")) %>%
    summarise(across(starts_with("sd_"), ~sqrt(mean(.^2, na.rm=TRUE))))

  # Create a named vector for easier matching
  within_sd_vec <- as.numeric(within_sds[1, ])
  names(within_sd_vec) <- gsub("sd_", "", names(within_sds))

  # Create a data frame with all posterior draws
  plot_data <- lapply(fixed_effects, function(var) {
    fixed_vals <- fixed_posteriors[, var]
    random_vals <- random_slope_posteriors[, var]

    # Standardize fixed effect by within-cluster SD (x-axis)
    fixed_vals_std <- fixed_vals * within_sd_vec[var]

    # Standardize random slope SD by within-cluster SD
    random_vals_std <- random_vals * within_sd_vec[var]

    # Reverse x-axis for distance variables
    if (var %in% c("anDistHdg", "anDistRf")) {
      fixed_vals_std <- -fixed_vals_std
    }

    # Calculate percentage of clusters with POSITIVE effect
    pct_positive <- pnorm(fixed_vals_std / random_vals_std) * 100

    data.frame(
      fixed_effect = fixed_vals_std,
      pct_positive = pct_positive,
      variable = var
    )
  }) %>% bind_rows()

  # Add metadata
  plot_data$linetype <- ifelse(plot_data$variable %in% c("anDistHdg", "anDistRf"), "dotted", "solid")
  plot_data$display_name <- display_name_map[plot_data$variable]
  plot_data$label <- label_map[plot_data$variable]

  return(list(
    plot_data = plot_data,
    fixed_effects = fixed_effects,
    use_empirical_sds = use_empirical_sds
  ))
}

# Calculate credible intervals from effect reliability data
calculate_credible_intervals_for_reliability <- function(plot_data, credible_level = 0.95) {

  lower_prob <- (1 - credible_level) / 2
  upper_prob <- 1 - lower_prob

  interval_data <- plot_data %>%
    group_by(variable) %>%
    summarise(
      fixed_median = median(fixed_effect),
      fixed_lower = quantile(fixed_effect, lower_prob),
      fixed_upper = quantile(fixed_effect, upper_prob),
      pct_median = median(pct_positive),
      pct_lower = quantile(pct_positive, lower_prob),
      pct_upper = quantile(pct_positive, upper_prob),
      .groups = "drop"
    )

  # Add metadata
  interval_data$linetype <- ifelse(interval_data$variable %in% c("anDistHdg", "anDistRf"), "dotted", "solid")
  interval_data$display_name <- display_name_map[interval_data$variable]
  interval_data$label <- label_map[interval_data$variable]

  return(interval_data)
}

# Print credible interval summary statistics
print_credible_intervals_for_reliability <- function(plot_data, interval_data, credible_level, use_empirical_sds) {

  lower_prob <- (1 - credible_level) / 2
  upper_prob <- 1 - lower_prob

  cat(sprintf("\n%.0f%% Credible Intervals (Standardized):\n", credible_level * 100))
  if (use_empirical_sds) {
    cat("Using empirical random slope SDs (informative gardens only)\n")
  } else {
    cat("Using hyperparameter random slope SDs\n")
  }
  cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

  fixed_effects <- unique(plot_data$variable)

  for (var in fixed_effects) {
    var_data <- plot_data %>% filter(variable == var)

    # Calculate credible intervals for fixed effect
    fixed_ci <- quantile(var_data$fixed_effect, probs = c(lower_prob, upper_prob))
    fixed_median <- median(var_data$fixed_effect)

    # Calculate credible intervals for percentage positive
    pct_ci <- quantile(var_data$pct_positive, probs = c(lower_prob, upper_prob))
    pct_median <- median(var_data$pct_positive)

    cat(sprintf("%-15s:\n", var))
    cat(sprintf("  Fixed Effect:      %7.3f [%7.3f, %7.3f]\n",
                fixed_median, fixed_ci[1], fixed_ci[2]))
    cat(sprintf("  %% Positive:        %7.2f%% [%6.2f%%, %6.2f%%]\n\n",
                pct_median, pct_ci[1], pct_ci[2]))
  }
  cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")
}

# Create effect reliability plot
create_effect_reliability_plot <- function(plot_data,
                                          interval_data,
                                          credible_level = 0.95,
                                          display_method = "interval_rectangle",
                                          x_min = NULL, x_max = NULL,
                                          y_min = NULL, y_max = NULL) {

  # Validate display_method parameter
  if (!display_method %in% c("credible_region", "interval_rectangle", "error_bars")) {
    stop("display_method must be 'credible_region', 'interval_rectangle', or 'error_bars'")
  }

  # Calculate contour break from credible level
  contour_break <- 1 - credible_level

  # Format credible level for axis labels
  ci_pct <- sprintf("%.0f%%", credible_level * 100)

  # Define linetype mapping for each variable
  linetype_map <- c(
    "hedge" = "solid",
    "grass" = "solid",
    "tree" = "solid",
    "roof" = "solid",
    "artif" = "solid",
    "bushAndCult" = "solid",
    "anDistHdg" = "31",    # Proximity variables get dashed
    "anDistRf" = "31"      # Proximity variables get dashed
  )

  # Create the base plot
  p <- ggplot(interval_data, aes(x = fixed_median, y = pct_median,
                                  color = variable, linetype = variable))

  # Add regions based on display_method
  if (display_method == "credible_region") {
    # Use plot_data for density contours
    p <- ggplot(plot_data, aes(x = fixed_effect, y = pct_positive,
                               fill = variable, color = variable, linetype = variable)) +
      stat_density_2d(geom = "polygon", aes(alpha = after_stat(level)),
                      contour_var = "ndensity", breaks = contour_break,
                      linewidth = 0.5) +
      scale_alpha_continuous(range = c(0.4, 0.4), guide = "none")

  } else if (display_method == "interval_rectangle") {
    # Rectangular credible intervals
    p <- p +
      geom_rect(data = interval_data,
                aes(xmin = fixed_lower, xmax = fixed_upper,
                    ymin = pct_lower, ymax = pct_upper,
                    fill = variable),
                alpha = 0.4, linewidth = 0.5)

  } else if (display_method == "error_bars") {
    # Error bars (no fill)
    p <- p +
      geom_errorbar(aes(ymin = pct_lower, ymax = pct_upper),
                    width = 0, linewidth = 0.8, alpha = 0.8) +
      geom_errorbarh(aes(xmin = fixed_lower, xmax = fixed_upper),
                     height = 0, linewidth = 0.8, alpha = 0.8) +
      geom_point(size = 2.5, alpha = 0.9)
  }

  # Add reference lines - now solid black and narrower
  p <- p +
    geom_hline(yintercept = 50, linetype = "solid", color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3)

  # Apply axis limits only if explicitly provided
  if (!is.null(x_min) || !is.null(x_max)) {
    p <- p + scale_x_continuous(limits = c(x_min, x_max))
  }
  if (!is.null(y_min) || !is.null(y_max)) {
    p <- p + scale_y_continuous(limits = c(y_min, y_max))
  }

  # Add labels differently for error bars (less cluttered)
  if (display_method == "error_bars") {
    p <- p +
      geom_text(data = interval_data,
                aes(x = fixed_median, y = pct_median, label = label, color = variable),
                size = 3, fontface = "bold",
                nudge_x = -0.045, nudge_y = 1.4,
                show.legend = FALSE)
  } else {
    # For filled methods, use white background lozenges
    p <- p +
      geom_label(data = interval_data,
                 aes(x = fixed_median, y = pct_median, label = label),
                 fill = "white", color = NA, alpha = 0.8,
                 label.size = 0, label.padding = unit(0.10, "lines"),
                 show.legend = FALSE) +
      geom_text(data = interval_data,
                aes(x = fixed_median, y = pct_median, label = label, color = variable),
                size = 3, fontface = "plain",
                show.legend = FALSE)
  }

  # Add labels and theme
  if (display_method == "error_bars") {
    p <- p +
      labs(
        x = sprintf("Effect size (standardised posterior median and %s credible interval)", ci_pct),
        y = sprintf("Effect consistency (%% gardens with positive effect)\n(posterior median and %s credible interval)", ci_pct),
        color = "Predictor",
        linetype = "Predictor"
      )
  } else {
    p <- p +
      labs(
        x = sprintf("Effect size (standardised posterior median and %s credible interval)", ci_pct),
        y = sprintf("Effect consistency (%% gardens with positive effect)\n(posterior median and %s credible interval)", ci_pct),
        fill = "Predictor",
        color = "Predictor",
        linetype = "Predictor"
      )
  }

  p <- p +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "white", linewidth = 0.3),
      panel.background = element_rect(fill = "grey85", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "grey85", color = NA),
      legend.key = element_rect(fill = "grey85", color = NA),
      legend.key.height = unit(0.9, "cm"),  # Increased spacing between legend items
      #legend.spacing.y = unit(0.3, "cm"),   # Additional vertical spacing
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 11, lineheight = 0.9),  # Added lineheight for multi-line text
      legend.title = element_text(size = 11)
    ) +
    scale_fill_manual(values = patch_colours, labels = display_name_map) +
    scale_color_manual(values = patch_colours, labels = display_name_map) +
    scale_linetype_manual(values = linetype_map, labels = display_name_map)

  # Merge the legends
  if (display_method %in% c("credible_region", "interval_rectangle")) {
    p <- p + guides(
      fill = guide_legend(override.aes = list(alpha = 0.4, linewidth = 0.8)),
      color = guide_legend(override.aes = list(alpha = 0.4, linewidth = 0.8)),
      linetype = guide_legend(override.aes = list(alpha = 0.4, linewidth = 0.8))
    )
  } else {
    p <- p + guides(
      color = guide_legend(override.aes = list(linewidth = 0.8)),
      linetype = guide_legend(override.aes = list(linewidth = 0.8))
    )
  }

  return(p)
}

# Plot effect consistency with credible intervals
# Again note that the internal functions refer to this as effect reliability
plot_effect_consistency <- function(model,
                                    data = NULL,
                                    use_empirical_sds = FALSE,
                                    id_col = "id",
                                    variance_threshold = 1e-10,
                                    credible_level = 0.95,
                                    display_method = "interval_rectangle",
                                    plot_file_name = "effect_reliability_plot.png",
                                    x_min = NULL, x_max = NULL,
                                    y_min = NULL, y_max = NULL) {

  # Calculate effect reliability data
  calc_results <- calculate_effect_reliability(
    model = model,
    data = data,
    use_empirical_sds = use_empirical_sds,
    id_col = id_col,
    variance_threshold = variance_threshold
  )

  # Calculate credible intervals
  interval_data <- calculate_credible_intervals_for_reliability(
    plot_data = calc_results$plot_data,
    credible_level = credible_level
  )

  # Print summary statistics
  print_credible_intervals_for_reliability(
    plot_data = calc_results$plot_data,
    interval_data = interval_data,
    credible_level = credible_level,
    use_empirical_sds = calc_results$use_empirical_sds
  )

  # Create plot
  p <- create_effect_reliability_plot(
    plot_data = calc_results$plot_data,
    interval_data = interval_data,
    credible_level = credible_level,
    display_method = display_method,
    x_min = x_min,
    x_max = x_max,
    y_min = y_min,
    y_max = y_max
  )

  # Display and save
  print(p)
  ggsave(plot_file_name, plot = p, bg = "white",
         width = 6.27, height = 5.0, units = "in", dpi = 300)

  return(p)
}

compare_random_slope_variability <- function(model, data, var_names,
                                             id_col = "id",
                                             variance_threshold = 1e-10) {
  # Call the core calculation function
  results <- calculate_unbiased_random_sds(model, data, var_names,
                                           id_col, variance_threshold)

  # Create summary table from results
  summary_table <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      variable = x$variable,
      n_informative = x$n_informative_gardens,
      n_total = x$n_total_gardens,
      pct_informative = x$pct_informative,
      empirical_sd_median = round(x$empirical_sd_median, 3),
      empirical_sd_95CI = sprintf("[%.3f, %.3f]", x$empirical_sd_q025, x$empirical_sd_q975),
      hyperpar_sd_median = round(x$hyperpar_sd_median, 3),
      hyperpar_sd_95CI = sprintf("[%.3f, %.3f]", x$hyperpar_sd_q025, x$hyperpar_sd_q975)
    )
  }))

  return(list(
    summary = summary_table,
    detailed_results = results
  ))
}

calculate_unbiased_random_sds <- function(model, data, var_names,
                                           id_col = "id",
                                           variance_threshold = 1e-10) {
  # Extract all random effects
  ranef_posterior <- as.matrix(model, regex_pars = "^b\\[")

  # Get all unique garden IDs from data
  all_ids <- unique(data[[id_col]])

  # Initialize results list
  results <- list()

  for (var_name in var_names) {
    # Identify informative gardens for this variable
    # Use very small threshold to exclude only floating-point error cases
    informative_gardens <- data %>%
      group_by(!!sym(id_col)) %>%
      summarise(
        var_within = var(!!sym(var_name), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(var_within > variance_threshold) %>%
      pull(!!sym(id_col))

    n_informative <- length(informative_gardens)

    # Extract columns for this variable using exact pattern
    var_cols <- grep(paste0("^b\\[", var_name, " ", id_col, ":"),
                     colnames(ranef_posterior),
                     value = TRUE)

    if (length(var_cols) == 0) {
      warning(paste("No random effects found for variable:", var_name))
      next
    }

    var_slopes <- ranef_posterior[, var_cols]

    # Extract garden IDs from column names like "b[anDistHdg id:01b26bec]"
    col_garden_ids <- gsub(paste0("^b\\[", var_name, " ", id_col, ":([^]]+)\\]$"),
                           "\\1",
                           colnames(var_slopes))

    # Convert underscores back to spaces in posterior IDs to match data IDs
    col_garden_ids <- gsub("_", " ", col_garden_ids)

    # Get columns corresponding to informative gardens
    informative_cols <- which(col_garden_ids %in% informative_gardens)

    if (length(informative_cols) == 0) {
      warning(paste("No informative gardens found for variable:", var_name))
      next
    }

    var_slopes_informative <- var_slopes[, informative_cols, drop = FALSE]

    # Calculate empirical SD across informative gardens for each posterior draw
    empirical_sd <- apply(var_slopes_informative, 1, sd)

    # Get hyperparameter using exact pattern: "Sigma[id:varname,varname]"
    hyperpar_name <- paste0("Sigma[", id_col, ":", var_name, ",", var_name, "]")
    posterior_df <- as.data.frame(model)
    hyperpar_sd <- NULL

    if (hyperpar_name %in% names(posterior_df)) {
      hyperpar_sd <- posterior_df[[hyperpar_name]]
    } else {
      warning(paste("Hyperparameter not found for variable:", var_name,
                   "(expected:", hyperpar_name, ")"))
    }

    # Store results
    results[[var_name]] <- list(
      variable = var_name,
      n_informative_gardens = n_informative,
      n_total_gardens = length(all_ids),
      pct_informative = round(100 * n_informative / length(all_ids), 1),
      empirical_sd_median = median(empirical_sd),
      empirical_sd_q025 = quantile(empirical_sd, 0.025),
      empirical_sd_q975 = quantile(empirical_sd, 0.975),
      hyperpar_sd_median = if (!is.null(hyperpar_sd)) median(hyperpar_sd) else NA,
      hyperpar_sd_q025 = if (!is.null(hyperpar_sd)) quantile(hyperpar_sd, 0.025) else NA,
      hyperpar_sd_q975 = if (!is.null(hyperpar_sd)) quantile(hyperpar_sd, 0.975) else NA,
      empirical_sd_full = empirical_sd,
      hyperpar_sd_full = hyperpar_sd
    )
  }

  return(results)
}