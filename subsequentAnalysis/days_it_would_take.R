library(lubridate)
library(ggplot2)

setwd("C:/Users/benke/Documents/gardenREDSmodelUpdateCode/subsequentAnalysis")

files <- c(
  "../fullProofOfConcept/consoleOutput/2025-05-20_162105_consoleOutput.txt",
  "../fullProofOfConcept/consoleOutput/2025-05-21_130005_consoleOutput.txt",
  "../fullProofOfConcept/consoleOutput/2025-05-22_130006_consoleOutput.txt",
  "../fullProofOfConcept/consoleOutput/2025-05-23_130002_consoleOutput.txt",
  "../fullProofOfConcept/consoleOutput/2025-05-26_130006_consoleOutput.txt"
)

results <- data.frame(
  start_time = ymd_hms(),
  end_time = ymd_hms(),
  n_gardens = integer()
)

for (i in seq_along(files)) {
  content <- readLines(files[i])
  
  start_time <- ymd_hms(sub("Start time: ", "", grep("^Start time:", content, value = TRUE)))
  end_time <- ymd_hms(sub("End time: ", "", grep("^End time:", content, value = TRUE)))
  n_gardens <- as.integer(sub("Number of unique gardens: ", "", grep("^Number of unique gardens:", content, value = TRUE)))
  
  results[i, ] <- list(start_time, end_time, n_gardens)
}

results$n_gardens <- results$n_gardens + 32

# Calculate elapsed hours with correction for May 22nd, when the workstation slept for 18 minutes 33 seconds during the run
results$elapsed_hours <- as.numeric(difftime(results$end_time, results$start_time, units = "hours"))
results$elapsed_hours[3] <- results$elapsed_hours[3] - (18 * 60 + 33) / 3600

model <- lm(elapsed_hours ~ n_gardens + I(n_gardens^2), data = results)

pred_gardens <- seq(min(results$n_gardens), max(results$n_gardens), length.out = 100)
pred_hours <- predict(model, newdata = data.frame(n_gardens = pred_gardens))

ggplot(results, aes(x = n_gardens, y = elapsed_hours)) +
  geom_point(size = 3, color = "blue") +
  geom_line(data = data.frame(n_gardens = pred_gardens, elapsed_hours = pred_hours),
            color = "red", size = 1) +
  labs(x = "Number of gardens and colonies", y = "Elapsed hours" ) +
  theme_minimal()
ggsave("model_run_time.png", bg="white", width = 2.5, height = 2.5, units = "in", dpi = 300)

hypo_double_batch_days <- (predict(model, newdata = data.frame(n_gardens = 173))) / 24
cross_validation_days <- (103 * results$elapsed_hours[5]) / 24
hypo_cross_validation_days <- (206 * predict(model, newdata = data.frame(n_gardens = 206))) / 24

{
cat("\n=== RESULTS ===\n")
cat("\nData:\n")
print(results[, c("n_gardens", "elapsed_hours")])
cat("\nModel coefficients:\n")
print(coef(model))
cat(sprintf("\nDays to if another batch the same: %.2f days\n", hypo_double_batch_days))
cat(sprintf("\nDays for full MCMC LOCO: %.2f days\n", cross_validation_days))
cat(sprintf("Days for full MCMC LOCO with doubled data: %.2f days\n", hypo_cross_validation_days))
}

