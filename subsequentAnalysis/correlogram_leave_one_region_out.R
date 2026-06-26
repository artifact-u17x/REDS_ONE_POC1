########################################################################################################
# LEAVE-ONE-BLOCK-OUT CROSS-VALIDATION  +  LOGO vs LOBO COMPARISON PLOT
#
# Conducts leave-one-region-out (LOBO) cross-validation and compares Brier scores with
# results from the leave-one-garden-out (LOGO) analysis.
#
#   PART 1–5 : leave-one-block-out (LOBO) cross-validation
#   PART 6–8 : plot comparing updated-model performance under LOGO vs LOBO
#
# Output goes to loco_results_blocked/ to avoid overwriting existing results.
#
# Requires in environment (and therefore must be run after analysis.R):
#   dFO, dF, finalCaseControlOffset
#   loco_results        — output of the main LOGO script (for the comparison plot)
#   GHSP_blue, GA_orange — colour constants
#   merge_glasgow_with_REDS_data(), add_id_columns(),
#   create_inla_formula(), fit_inla_model(), brier_decomposition_own_rolled()
#   garden_coordinates_rural_urban_region.csv in working directory
#   dFBothFinalOffset, modUpdatedDay5FinalOffset, compute_predictions()
########################################################################################################

library(INLA)
library(ggplot2)

LOBO_purple <- "#7B4F9E"

########################################################################################################
# CORRELOGRAM: spatial autocorrelation of fixed-effect residuals
########################################################################################################

library(rstanarm)
library(ncf)

# --- Coordinates ---
coords_corr <- read.csv("garden_coordinates.csv", stringsAsFactors = FALSE)
coords_corr$id <- sub("-.*", "", coords_corr$garden.guid)

# --- Fixed-effect residuals ---
p_fixed <- compute_predictions(
  model      = modUpdatedDay5FinalOffset,
  data       = dFBothFinalOffset,
  offset     = dFBothFinalOffset$prevalenceOffset,
  useRandom  = FALSE,
  model_type = "stan_glmer"
)
dFBothFinalOffset$residual <- dFBothFinalOffset$sparrow - as.numeric(p_fixed)

# --- Aggregate to location level ---
loc_data <- merge(
  aggregate(residual ~ id, data = dFBothFinalOffset, FUN = mean),
  coords_corr[, c("id", "epsg3857.x", "epsg3857.y")],
  by = "id"
)

# --- Distance summaries ---
coords_mat     <- as.matrix(loc_data[, c("epsg3857.x", "epsg3857.y")])
dist_mat       <- as.matrix(dist(coords_mat))
diag(dist_mat) <- NA
median_nn      <- median(apply(dist_mat, 1, min, na.rm = TRUE))
dist_75        <- quantile(dist_mat, 0.75, na.rm = TRUE)
cat(sprintf("Bin width (median nearest-neighbour distance): %.2f km\n", median_nn / 1000))
cat(sprintf("Cutoff (75th pct pairwise distance):           %.2f km\n", dist_75   / 1000))

# --- Correlogram ---
corr_fit <- correlog(
  x         = loc_data$epsg3857.x,
  y         = loc_data$epsg3857.y,
  z         = loc_data$residual,
  increment = median_nn,
  latlon    = FALSE,
  resamp    = 999,
  quiet     = TRUE
)
corr_df <- data.frame(
  distance_km = corr_fit$mean.of.class / 1000,
  correlation = corr_fit$correlation,
  n_pairs     = corr_fit$n
)
corr_df <- corr_df[!is.na(corr_df$distance_km) &
                     corr_df$distance_km <= dist_75 / 1000, ]

# --- Plot ---
ggplot(corr_df, aes(x = distance_km, y = correlation)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_smooth(aes(weight = n_pairs), method = "lm", se = FALSE,
              colour = "firebrick", linewidth = 0.8) +
  geom_line(colour = "steelblue") +
  geom_point(aes(size = n_pairs), shape = 21,
             colour = "steelblue", fill = "steelblue", alpha = 0.8) +
  scale_size_continuous(name = "Pairs (n)", range = c(2, 6)) +
  labs(
    x = "Distance (km)",
    y = "Moran's I"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title       = element_text(size = 13)
  )
ggsave("correlogram.png", width = 8, height = 5, dpi = 150)
cat("Saved: correlogram.png\n")


########################################################################################################
# PART 1: Load spatial blocks from region_band classification
########################################################################################################

coords <- read.csv("garden_coordinates_rural_urban_region.csv", stringsAsFactors = FALSE)
coords$id <- sub("-.*", "", coords$garden.guid)
stopifnot(length(unique(coords$id)) == nrow(coords))
stopifnot("region_band" %in% names(coords))

blocks_df <- data.frame(
  id    = coords$id,
  block = coords$region_band,
  stringsAsFactors = FALSE
)

{
  cat("=== Spatial blocks (region_band classification) ===\n")
  sz <- table(blocks_df$block)
  print(sz)
  cat(sprintf("Min: %d  |  Median: %.0f  |  Max: %d  |  Total locations: %d\n",
              min(sz), median(sz), max(sz), sum(sz)))
}

ggplot(merge(blocks_df, coords, by = "id"),
       aes(x = epsg3857.x / 1000, y = epsg3857.y / 1000,
           colour = factor(block), label = block)) +
  geom_point(size = 3) +
  geom_text(nudge_y = 15, size = 2.5) +
  labs(title    = "Spatial blocks for leave-one-block-out CV",
       subtitle = "region_band classification  |  colour = block",
       x = "Easting (km)", y = "Northing (km)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")


########################################################################################################
# PART 2: Block-aware model fitting function
#
# Fits the updated model (Glasgow + REDS) with all gardens in the block masked,
# then returns per-garden Brier scores and decompositions for the held-out block.
########################################################################################################

fit_lobo_block <- function(dFBoth, holdout_garden_ids, u = 10, alpha = 0.5,
                            own_rolled_brier = TRUE) {

  stopifnot(all(holdout_garden_ids %in% dFBoth$id))

  modified_data <- dFBoth
  modified_data$sparrow[dFBoth$id %in% holdout_garden_ids] <- NA
  modified_data <- add_id_columns(modified_data)

  prior_spec <- list(prec = list(prior = "pc.prec", param = c(u, alpha)))
  model      <- fit_inla_model(create_inla_formula(prior_spec), modified_data)

  garden_results <- list()

  for (garden in holdout_garden_ids) {
    garden_mask   <- dFBoth$id == garden
    holdout_data  <- dFBoth[garden_mask, ]
    holdout_preds <- model$summary.fitted.values[garden_mask, ]

    temp_data <- data.frame(
      sparrow         = holdout_data$sparrow,
      predicted_probs = holdout_preds$mean
    )

    has_both_classes <- length(table(holdout_data$sparrow)) == 2

    sink(nullfile())
    if (has_both_classes) {
      if (own_rolled_brier) {
        brier_result  <- brier_decomposition_own_rolled(temp_data$sparrow, temp_data$predicted_probs)
      } else {
        brier_result  <- scoring::brierscore(sparrow ~ predicted_probs, temp_data, decomp = TRUE)
      }
      brier_scores   <- brier_result$rawscores
      discrimination <- brier_result$decomp$components["discrim", 1]
      miscalibration <- brier_result$decomp$components["miscal",  1]
    } else {
      brier_scores   <- scoring::brierscore(sparrow ~ predicted_probs, temp_data)
      discrimination <- NA
      miscalibration <- NA
    }
    sink()

    garden_results[[garden]] <- list(
      brier_scores   = brier_scores,
      discrimination = discrimination,
      miscalibration = miscalibration
    )
  }

  return(garden_results)
}


########################################################################################################
# PART 3: Block-level loop
########################################################################################################

run_lobo_comparison <- function(dFO, dF, caseControlOffset,
                                 blocks_df,
                                 output_dir        = "loco_results_blocked",
                                 own_rolled_brier  = TRUE,
                                 prediction_method = "holdout") {

  method_dir <- file.path(output_dir, prediction_method)
  if (!dir.exists(method_dir)) dir.create(method_dir, recursive = TRUE)

  brier_file  <- file.path(method_dir, "brier_scores.rds")
  decomp_file <- file.path(method_dir, "decomp_scores.rds")

  # Resume support
  if (file.exists(brier_file) && file.exists(decomp_file)) {
    all_brier_scores <- readRDS(brier_file)
    decomp_scores    <- readRDS(decomp_file)
    if (nrow(all_brier_scores) > 0) {
      processed_gardens <- unique(all_brier_scores$garden)
      cat("Found existing results with", length(processed_gardens), "gardens already processed\n")
    } else {
      processed_gardens <- character(0)
      cat("Found empty result files - starting from beginning\n")
    }
  } else {
    all_brier_scores  <- data.frame()
    decomp_scores     <- data.frame()
    processed_gardens <- character(0)
    cat("No existing results found - starting fresh\n")
  }

  dFBoth <- merge_glasgow_with_REDS_data(dFO, dF, caseControlOffset)

  gardens_to_process <- setdiff(unique(dF$id), processed_gardens)
  if (length(gardens_to_process) == 0) {
    cat("All gardens already processed!\n")
    return(list(brier_scores = all_brier_scores, decomp_scores = decomp_scores))
  }

  blocks_to_process <- unique(blocks_df$block[blocks_df$id %in% gardens_to_process])
  cat("Processing", length(blocks_to_process), "blocks covering",
      length(gardens_to_process), "remaining gardens\n")

  for (bi in seq_along(blocks_to_process)) {
    block_id      <- blocks_to_process[bi]
    block_gardens <- blocks_df$id[blocks_df$block == block_id]
    block_gardens_to_eval <- intersect(block_gardens, gardens_to_process)
    if (length(block_gardens_to_eval) == 0) next

    cat(sprintf("Processing block '%s' (%d of %d)  |  gardens: %s\n",
                block_id, bi, length(blocks_to_process),
                paste(block_gardens_to_eval, collapse = ", ")))

    block_results <- fit_lobo_block(
      dFBoth, block_gardens_to_eval, own_rolled_brier = own_rolled_brier
    )

    for (garden in block_gardens_to_eval) {
      res   <- block_results[[garden]]
      n_obs <- length(res$brier_scores)

      all_brier_scores <- rbind(all_brier_scores, data.frame(
        garden   = rep(garden, n_obs),
        all_data = res$brier_scores
      ))

      decomp_scores <- rbind(decomp_scores, data.frame(
        garden                  = garden,
        all_data_discrimination = res$discrimination,
        all_data_miscalibration = res$miscalibration
      ))
    }

    saveRDS(all_brier_scores, brier_file)
    saveRDS(decomp_scores,    decomp_file)
    gc()
  }

  cat("Processing complete! Total gardens processed:",
      length(unique(all_brier_scores$garden)), "\n")

  return(list(brier_scores = all_brier_scores, decomp_scores = decomp_scores))
}


########################################################################################################
# PART 4: Analysis function
#
# Simpler version of analyze_brier_comparison() for LOBO: models the raw all_data
# scores directly rather than an improvement score (no old_only column exists here).
# Produces a lobo_results$new with the same named model slots
# (modBrierRaws, modDiscrimination, modMiscalibration) that plot_logo_vs_lobo() expects.
########################################################################################################

analyze_brier_scores_lobo <- function(brier_data, decomp_data, run_quickly = TRUE) {

  brier_means <- aggregate(all_data ~ garden, data = brier_data, FUN = mean)

  if (!run_quickly) {
    mod_brier <- stan_glmer(all_data ~ 1 + (1 | garden), data = brier_data)
  } else {
    mod_brier <- stan_glm(all_data ~ 1, data = brier_means, family = gaussian())
  }

  mod_discrimination <- stan_glm(all_data_discrimination ~ 1,
                                  data = decomp_data, family = gaussian())
  mod_miscalibration <- stan_glm(all_data_miscalibration ~ 1,
                                  data = decomp_data, family = gaussian())

  return(list(
    brier_means       = brier_means,
    modBrierRaws      = mod_brier,
    modDiscrimination = mod_discrimination,
    modMiscalibration = mod_miscalibration
  ))
}


########################################################################################################
# PART 5: Run the LOBO cross-validation
########################################################################################################

loco_output_dir_blocked <- "loco_results_blocked"
prediction_method       <- "holdout"

{
  if (runSlowModels) {
    cat("Starting leave-one-block-out cross-validation...\n")
    lobo_cv_results <- run_lobo_comparison(
      dFO,
      dF,
      finalCaseControlOffset,
      blocks_df         = blocks_df,
      own_rolled_brier  = TRUE,
      output_dir        = loco_output_dir_blocked,
      prediction_method = prediction_method
    )

    cat("\nAnalyzing results...\n")
    lobo_analysis <- analyze_brier_scores_lobo(
      lobo_cv_results$brier_scores,
      lobo_cv_results$decomp_scores,
      run_quickly = run_quickly
    )

    lobo_results <- list(
      new          = lobo_analysis,
      full_results = lobo_cv_results
    )

    saveRDS(lobo_results,
            file.path(loco_output_dir_blocked,
                      paste0("lobo_analysis_", prediction_method, ".rds")))
  } else {
    lobo_results <- readRDS(
      file.path(loco_output_dir_blocked,
                paste0("lobo_analysis_", prediction_method, ".rds"))
    )
  }
}


########################################################################################################
# PART 6: HELPER — posterior median + 95% CI from a stan_glm / stan_glmer model
########################################################################################################

extract_ci_and_median <- function(model, prob = 0.95) {
  ci  <- posterior_interval(model, pars = "(Intercept)", prob = prob)
  med <- median(as.matrix(model)[, "(Intercept)"])
  list(lower = ci[1], upper = ci[2], median = med)
}


########################################################################################################
# PART 7: HELPER — build a tidy data frame of estimates for one CV scheme
#
# For LOGO, modBrierRaws models the improvement score rather than raw performance,
# so we refit a quick stan_glm on the all_data column of brierRaws directly.
# (Discrimination and Miscalibration are already raw all_data models on both sides.)
########################################################################################################

extract_loXo_estimates <- function(results, scheme_label) {

  if (scheme_label == "LOGO") {
    brier_means <- aggregate(all_data ~ garden,
                             data = results$brierRaws, FUN = mean)
    brier_mod   <- stan_glm(all_data ~ 1, data = brier_means, family = gaussian(),
                             refresh = 0)
    brier <- extract_ci_and_median(brier_mod)
  } else {
    brier <- extract_ci_and_median(results$modBrierRaws)
  }

  discr  <- extract_ci_and_median(results$modDiscrimination)
  miscal <- extract_ci_and_median(results$modMiscalibration)

  data.frame(
    scheme = scheme_label,
    metric = c("brier", "discrimination", "miscalibration"),
    median = c(brier$median,  discr$median,  miscal$median),
    lower  = c(brier$lower,   discr$lower,   miscal$lower),
    upper  = c(brier$upper,   discr$upper,   miscal$upper),
    stringsAsFactors = FALSE
  )
}


########################################################################################################
# PART 8: MAIN PLOT FUNCTION + CALL
#
# Compares updated-model performance (REDS gardens only) under LOGO vs LOBO.
########################################################################################################

plot_logo_vs_lobo <- function(loco_results, lobo_results,
                               GHSP_blue, GA_orange,
                               LOBO_purple = "#7B4F9E",
                               save_plot   = "brier_plot_logo_vs_lobo.png") {

  logo_est <- extract_loXo_estimates(loco_results$new, "LOGO")
  lobo_est <- extract_loXo_estimates(lobo_results$new, "LOBO")

  plot_data <- rbind(logo_est, lobo_est)

  metric_labels <- c(
    brier          = "\nBrier score",
    discrimination = "\nDiscrimination",
    miscalibration = "\nMiscalibration"
  )
  plot_data$metric_label <- factor(
    metric_labels[plot_data$metric],
    levels = metric_labels
  )
  plot_data$scheme <- factor(plot_data$scheme, levels = c("LOGO", "LOBO"))

  scheme_colours <- c(LOGO = GA_orange, LOBO = LOBO_purple)

  p <- ggplot(plot_data,
              aes(x = scheme, y = median,
                  ymin = lower, ymax = upper,
                  colour = scheme)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_pointrange(linewidth = 0.7, fatten = 3) +
    facet_wrap(~ metric_label, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = scheme_colours, guide = "none") +
    labs(
      title    = "Updated model performance: LOGO vs LOBO (REDS gardens)",
      subtitle = "Posterior median \u00b1 95% CI",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text       = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 11, face = "bold"),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      plot.background  = element_rect(fill = "white", colour = NA)
    )

  print(p)

  {
    cat("\n=== Updated model performance: LOGO vs LOBO (REDS gardens) ===\n")
    for (m in c("brier", "discrimination", "miscalibration")) {
      cat(sprintf("\n%s:\n", toupper(m)))
      for (s in c("LOGO", "LOBO")) {
        row <- plot_data[plot_data$metric == m & plot_data$scheme == s, ]
        cat(sprintf("  %s:  %.4f  [%.4f, %.4f]\n",
                    s, row$median, row$lower, row$upper))
      }
    }
    cat("\n")
  }

  if (!is.null(save_plot)) {
    ggsave(save_plot, plot = p,
           width = 6.27, height = 3.8, dpi = 300, bg = "white")
    cat(sprintf("Plot saved as '%s'\n", save_plot))
  }

  invisible(list(plot = p, data = plot_data))
}


plot_logo_vs_lobo(loco_results, lobo_results,
                  GHSP_blue   = GHSP_blue,
                  GA_orange   = GA_orange,
                  LOBO_purple = LOBO_purple)