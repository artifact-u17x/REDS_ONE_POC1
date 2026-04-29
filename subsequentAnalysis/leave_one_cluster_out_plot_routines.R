##############################################################################
# LEAVE-ONE-OUT COMPARISON PLOTTING SYSTEM
##############################################################################

# Extract credible intervals and medians from LOCO results
extract_credible_intervals_and_medians_for_loco_plot <- function(loco_results) {

  get_ci <- function(model, prob = 0.95) {
    ci <- posterior_interval(model, pars = "(Intercept)", prob = prob)
    ci[1, ]  # Returns c(lower, upper)
  }

  get_median <- function(model) {
    posterior <- as.matrix(model)
    median(posterior[, "(Intercept)"])
  }

  list(
    ci = list(
      discrimination = list(
        old = get_ci(loco_results$old$modDiscrimination),
        new = get_ci(loco_results$new$modDiscrimination),
        all = get_ci(loco_results$all$modDiscrimination)
      ),
      miscalibration = list(
        old = get_ci(loco_results$old$modMiscalibration),
        new = get_ci(loco_results$new$modMiscalibration),
        all = get_ci(loco_results$all$modMiscalibration)
      ),
      brier = list(
        old = get_ci(loco_results$old$modBrierRaws),
        new = get_ci(loco_results$new$modBrierRaws),
        all = get_ci(loco_results$all$modBrierRaws)
      )
    ),
    medians = list(
      discrimination = list(
        old = get_median(loco_results$old$modDiscrimination),
        new = get_median(loco_results$new$modDiscrimination),
        all = get_median(loco_results$all$modDiscrimination)
      ),
      miscalibration = list(
        old = get_median(loco_results$old$modMiscalibration),
        new = get_median(loco_results$new$modMiscalibration),
        all = get_median(loco_results$all$modMiscalibration)
      ),
      brier = list(
        old = get_median(loco_results$old$modBrierRaws),
        new = get_median(loco_results$new$modBrierRaws),
        all = get_median(loco_results$all$modBrierRaws)
      )
    )
  )
}

# Prepare decomposition data for plotting
prepare_decomp_data_for_loco_plot <- function(decomp_raws, glasgow_gardens) {
  decomp_raws %>%
    mutate(data_source = ifelse(garden %in% glasgow_gardens, "glasgow", "garden_advice")) %>%
    select(garden, data_source,
           old_only_discrimination, all_data_discrimination,
           old_only_miscalibration, all_data_miscalibration) %>%
    pivot_longer(
      cols = c(old_only_discrimination, all_data_discrimination,
               old_only_miscalibration, all_data_miscalibration),
      names_to = "metric_model",
      values_to = "value"
    ) %>%
    separate(metric_model, into = c("model", "metric"),
             sep = "_(?=discrimination|miscalibration)") %>%
    mutate(model = factor(model, levels = c("old_only", "all_data")))
}

# Prepare Brier score data for plotting
prepare_brier_data_for_loco_plot <- function(brier_raws, glasgow_gardens) {
  brier_raws %>%
    mutate(data_source = ifelse(garden %in% glasgow_gardens, "glasgow", "garden_advice")) %>%
    group_by(garden, data_source) %>%
    summarise(
      old_only = mean(old_only, na.rm = TRUE),
      all_data = mean(all_data, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(old_only, all_data),
      names_to = "model",
      values_to = "value"
    ) %>%
    mutate(
      metric = "brier",
      model = factor(model, levels = c("old_only", "all_data"))
    )
}

# Calculate means by data source and overall for plotting
calculate_means_for_loco_plot <- function(all_data) {
  # Means by data source
  means_by_source <- all_data %>%
    group_by(model, metric, data_source) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(data_source_plot = recode(data_source,
                                     "glasgow" = "Glasgow",
                                     "garden_advice" = "Garden Advice"))

  # Overall means across all gardens
  means_overall <- all_data %>%
    group_by(model, metric) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      data_source = "both",
      data_source_plot = "Both"
    )

  # Combine and set factor levels
  bind_rows(means_by_source, means_overall) %>%
    mutate(data_source_plot = factor(data_source_plot,
                                     levels = c("Glasgow", "Garden Advice", "Both")))
}

# Create summary table of means for display
create_means_table_for_loco_plot <- function(all_data) {
  # Recode data sources for table
  all_data_recoded <- all_data %>%
    mutate(data_source = recode(data_source,
                                 "glasgow" = "old",
                                 "garden_advice" = "new"))

  # Calculate means by source
  means_by_source <- all_data_recoded %>%
    group_by(model, metric, data_source) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

  # Calculate overall means
  means_overall <- all_data_recoded %>%
    group_by(model, metric) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(data_source = "old_and_new")

  # Combine and format
  bind_rows(means_by_source, means_overall) %>%
    pivot_wider(names_from = metric, values_from = mean_value) %>%
    mutate(
      model_training_data = recode(as.character(model),
                                    "old_only" = "Original",
                                    "all_data" = "Updated"),
      model_evaluation_data = recode(data_source,
                                      "old" = "Glasgow",
                                      "new" = "Garden Advice",
                                      "old_and_new" = "Both")
    ) %>%
    select(model_training_data, model_evaluation_data,
           brier_overall = brier,
           brier_decomposition_discrimination = discrimination,
           brier_decomposition_miscalibration = miscalibration) %>%
    arrange(model_training_data,
            factor(model_evaluation_data, levels = c("Glasgow", "Garden Advice", "Both")))
}

# Apply consistent model and data source labels for plotting
apply_plot_labels_for_loco_plot <- function(data, GHSP_label, updated_label) {
  data %>%
    mutate(
      data_source_plot = recode(data_source,
                                "glasgow" = "Glasgow",
                                "garden_advice" = "Garden Advice",
                                "both" = "Both"),  # Add this line
      model_plot = recode(as.character(model),
                         "old_only" = GHSP_label,
                         "all_data" = updated_label),
      model_plot = factor(model_plot, levels = c(GHSP_label, updated_label))
    )
}

############ MAIN PLOTTING FUNCTIONS

# Create a single metric plot with credible intervals
create_metric_plot_for_loco_plot <- function(data, means, metric_name, title, invert_y,
                                              offsets,
                                              GHSP_label, updated_label,
                                              GHSP_blue, GA_orange) {

  # Filter data for this metric
  metric_data_df <- data %>% filter(metric == metric_name)

  # Prepare means for this metric
  means_colored <- means %>% filter(metric == metric_name, data_source != "both")
  means_both <- means %>% filter(metric == metric_name, data_source == "both")

  # Create base plot with garden-level data
  p <- ggplot(metric_data_df, aes(x = model_plot, y = value,
                                   group = garden, color = data_source_plot)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_line(alpha = 0.3, linewidth=0.3) +
    geom_point(alpha = 0.5, size=1) +
    # Overlay colored means lines (Glasgow and Garden Advice)
    geom_line(data = means_colored,
              aes(x = model_plot, y = mean_value, group = data_source_plot,
                  color = data_source_plot),
              linewidth = 1.2, alpha = 0.9, inherit.aes = FALSE) +
    geom_point(data = means_colored,
               aes(x = model_plot, y = mean_value, color = data_source_plot),
               size = 3.5, alpha = 0.9, inherit.aes = FALSE) +
    # Overlay black "Both" line
    geom_line(data = means_both,
              aes(x = model_plot, y = mean_value, group = data_source_plot,
                  color = data_source_plot),
              linewidth = 1.2, alpha = 0.9, inherit.aes = FALSE) +
    geom_point(data = means_both,
               aes(x = model_plot, y = mean_value, color = data_source_plot),
               size = 3.5, alpha = 0.9, inherit.aes = FALSE) +
    labs(title = title, x = NULL, y = NULL, color = "Dataset for prediction:  ") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9.5), # ADJUST panel header text size
      axis.text.x = element_text(size = 9), # ADJUST x-axis label text size
      axis.text.y = element_text(size = 8.5), # ADJUST y-axis label text size
      axis.title.x = element_text(size = 10, margin = margin(t = 8)), # ADJUST x-axis title text size
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor.y = element_line(color = "grey95"),
      legend.text = element_text(size = 10), # ADJUST legend item text size
      legend.title = element_text(size = 10) # ADJUST legend title text size
    ) +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.2, 0.2))) +
    scale_color_manual(
      values = c("Glasgow" = GHSP_blue, "Garden Advice" = GA_orange, "Both" = "black"),
      breaks = c("Glasgow", "Garden Advice", "Both"),
      labels = c("GHSP", "GA", "All data"),
      guide = guide_legend(order = 1)
    )

  # Add appropriate y-axis scale
  if (invert_y) {
    p <- p + scale_y_reverse(limits = c(NA, 0))
  } else {
    p <- p + scale_y_continuous(limits = c(0, NA))
  }

  p
}

# Create change interval plot for a single metric
create_change_plot_for_loco_plot <- function(ci_old, ci_new, ci_all,
                                              median_old = NULL, median_new = NULL, median_all = NULL,
                                              title,
                                              flip_sign = FALSE,
                                              x_expand = c(0.2, 0.2),
                                              GHSP_blue, GA_orange) {

  # Flip signs if needed (for brier and miscalibration)
  if (flip_sign) {
    ci_old <- -ci_old
    ci_new <- -ci_new
    ci_all <- -ci_all
    if (!is.null(median_old)) median_old <- -median_old
    if (!is.null(median_new)) median_new <- -median_new
    if (!is.null(median_all)) median_all <- -median_all
  }

  # Prepare data frame for change intervals
  change_data <- data.frame(
    data_source = factor(c("Glasgow", "Garden Advice", "Both"),
                        levels = c("Glasgow", "Garden Advice", "Both")),
    lower = c(ci_old[1], ci_new[1], ci_all[1]),
    upper = c(ci_old[2], ci_new[2], ci_all[2]),
    median = c(
      if(!is.null(median_old)) median_old else mean(ci_old),
      if(!is.null(median_new)) median_new else mean(ci_new),
      if(!is.null(median_all)) median_all else mean(ci_all)
    )
  )

  # Create plot
  p <- ggplot(change_data, aes(x = data_source, y = median, color = data_source)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 1) +
    geom_point(size = 3) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9.5), # ADJUST panel header text size
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8.5), # ADJUST y-axis label text size
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 10, margin = margin(t = 8)), # ADJUST x-axis title text size
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor.y = element_line(color = "grey95"),
      legend.position = "none"
    ) +
    scale_color_manual(
      values = c("Glasgow" = GHSP_blue, "Garden Advice" = GA_orange, "Both" = "black")
    ) +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = x_expand))

  # Reverse y-axis for metrics where we flipped the sign
  # This makes improvement (now negative after flipping) appear upward
  if (flip_sign) {
    p <- p + scale_y_reverse()
  }

  p
}

# Create arrow plot for shared y-axis label
create_arrow_plot_for_loco_plot <- function() {
  ggplot() +
    annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = 1,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             linewidth = 1.2) +
    annotate("text", x = 0.5, y = 0.5, label = "Prediction performance",
             angle = 90, vjust = -1.5, size = 3.5) + # ADJUST "Prediction performance" text size
    theme_void() +
    coord_cartesian(clip = "off")
}

# Combine individual metric plots with change plots, arrow and shared legend
combine_plots_for_loco_plot <- function(p1, p2, p3, c1, c2, c3, arrow_plot) {
  # Apply horizontal legend to metric plots
  legend_theme <- theme(legend.direction = "horizontal")

  p1 <- p1 + legend_theme
  p2 <- p2 + legend_theme + xlab("Before and after model update")
  p3 <- p3 + legend_theme

  # Add shared x-axis label to middle change panel
  c2 <- c2 + xlab("Change after updating")

  # Define layout:
  # Row 1: arrow, three metric plots, three change plots
  # Row 2: guide area spanning all columns
  design <- "
    A123456
    .GGGGGG
  "

  patchwork::wrap_plots(
    A = arrow_plot,
    `1` = p1,
    `2` = p2,
    `3` = p3,
    `4` = c1,
    `5` = c2,
    `6` = c3,
    G = patchwork::guide_area(),
    design = design,
    # ADJUST THESE to change panel widths: c(arrow, abs1, abs2, abs3, change1, change2, change3)
    # Currently change panels (0.6) are 60% the width of absolute panels (1.0)
    # To make change panels wider, increase their values (e.g., 0.8 or 1.0)
    # First value (0.5) controls arrow panel width - ADJUST to reduce space on right side of arrow
    widths = c(0.35, .8, .8, .8, 0.7, 0.7, 0.75),
    heights = c(1, 0.1),
    guides = "collect"
  )
}

# MAIN INTERFACE FUNCTION
# Plot LOCO (Leave-One-Cohort-Out) results
plot_loco <- function(loco_results, glasgow_gardens ) {

  # Define model labels
  GHSP_label <- "Before"
  updated_label <- "After"

  # Extract results for main plotting
  results <- loco_results$all

  # Extract credible intervals and medians
  ci_and_medians <- extract_credible_intervals_and_medians_for_loco_plot(loco_results)
  credible_intervals <- ci_and_medians$ci
  medians <- ci_and_medians$medians

  # Prepare data for plotting
  decomp_data <- prepare_decomp_data_for_loco_plot(results$decompRaws, glasgow_gardens)
  brier_data <- prepare_brier_data_for_loco_plot(results$brierRaws, glasgow_gardens)
  all_data <- bind_rows(decomp_data, brier_data)

  # Calculate means and apply plot labels
  means_for_plot <- calculate_means_for_loco_plot(all_data)
  plot_data <- apply_plot_labels_for_loco_plot(all_data, GHSP_label, updated_label)
  means_for_plot <- apply_plot_labels_for_loco_plot(means_for_plot, GHSP_label, updated_label)

  # Create individual metric plots (without CI annotations on secondary axis)
  # Order: Brier, Discrimination, Miscalibration
  p1 <- create_metric_plot_for_loco_plot(
    plot_data, means_for_plot, "brier", "\nBrier score",
    invert_y = TRUE,
    offsets = list(),
    GHSP_label, updated_label, GHSP_blue, GA_orange
  )

  p2 <- create_metric_plot_for_loco_plot(
    plot_data, means_for_plot, "discrimination", "\nDiscrimination",
    invert_y = FALSE,
    offsets = list(),
    GHSP_label, updated_label, GHSP_blue, GA_orange
  )

  p3 <- create_metric_plot_for_loco_plot(
    plot_data, means_for_plot, "miscalibration", "\nMiscalibration",
    invert_y = TRUE,
    offsets = list(),
    GHSP_label, updated_label, GHSP_blue, GA_orange
  )

  # Create change interval plots
  c1 <- create_change_plot_for_loco_plot(
    ci_old = credible_intervals$brier$old,
    ci_new = credible_intervals$brier$new,
    ci_all = credible_intervals$brier$all,
    median_old = medians$brier$old,
    median_new = medians$brier$new,
    median_all = medians$brier$all,
    title = "Δ\nBrier score",
    flip_sign = TRUE,
    x_expand = c(0.2, 0.2),
    GHSP_blue, GA_orange
  )

  c2 <- create_change_plot_for_loco_plot(
    ci_old = credible_intervals$discrimination$old,
    ci_new = credible_intervals$discrimination$new,
    ci_all = credible_intervals$discrimination$all,
    median_old = medians$discrimination$old,
    median_new = medians$discrimination$new,
    median_all = medians$discrimination$all,
    title = "Δ\nDiscrimination",
    flip_sign = FALSE,
    x_expand = c(0.2, 0.2),
    GHSP_blue, GA_orange
  )

  c3 <- create_change_plot_for_loco_plot(
    ci_old = credible_intervals$miscalibration$old,
    ci_new = credible_intervals$miscalibration$new,
    ci_all = credible_intervals$miscalibration$all,
    median_old = medians$miscalibration$old,
    median_new = medians$miscalibration$new,
    median_all = medians$miscalibration$all,
    title = "Δ\nMiscalibration",
    flip_sign = TRUE,
    x_expand = c(0.2, 0.2),
    GHSP_blue, GA_orange
  )

  # Combine all plots with arrow and shared legend
  arrow_plot <- create_arrow_plot_for_loco_plot()
  combined_plot <- combine_plots_for_loco_plot(p1, p2, p3, c1, c2, c3, arrow_plot)

  # Create summary table
  means_table <- create_means_table_for_loco_plot(all_data)

  # Display output
  print(combined_plot)
  cat("\n=== Means by Model and Evaluation Data ===\n")
  print(means_table, row.names = FALSE)

  ga_brier <- means_table %>% filter(model_evaluation_data == "Garden Advice")
  original_brier_ga <- ga_brier$brier_overall[ga_brier$model_training_data == "Original"]
  updated_brier_ga  <- ga_brier$brier_overall[ga_brier$model_training_data == "Updated"]
  delta_brier_ga    <- updated_brier_ga - original_brier_ga
  delta_rmse_ga     <- sqrt(updated_brier_ga) - sqrt(original_brier_ga)

  cat(sprintf("\nΔ Brier score (Garden Advice): %.4f\n", delta_brier_ga))
  cat(sprintf("Δ RMSE / prob-scale improvement (Garden Advice): %.4f\n", delta_rmse_ga))
  
  # Save plot
  ggsave("brier_plot.png",
         plot = combined_plot,
         width = 6.27,  # Increased width to accommodate change panels
         height = 4.5,
         dpi = 300,
         bg = "white")
  cat("\nPlot saved as 'brier_plot.png'\n")

  # Return invisibly for potential further use
  invisible(list(
    plot = combined_plot,
    means_table = means_table
  ))
}
