start.time <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

gardenNow <- read.csv(input_file) # load provided data

fixed_coefs <- c(
  Intercept = -0.1888526,
  Hedge = 1.2158663,
  Grass = -0.9773279,
  Tree = 1.4508558,
  Roof = 0.1177534,
  Artificial.Surface = -1.6754525,
  Bush = 1.4405024,
  DistHdg = -0.2054778,
  DistRf = -0.8187333
)

predictFunctionForUse <- function( newdata ) {
  distRfInverseScale <- 7.212275
  minTransformedDistRf <- -1.3929
  distHdgInverseScale <- 11.26151
  minTransformedDistHdg <- -1.148933
  newdata$DistRf <- ifelse( is.na(newdata$DistRf) | newdata$DistRf > 30, 30, newdata$DistRf)
  newdata$DistHdg <- ifelse( is.na(newdata$DistHdg) | newdata$DistHdg > 50, 50, newdata$DistHdg)
  newdata$DistRf <- newdata$DistRf / distRfInverseScale + minTransformedDistRf
  newdata$DistHdg <- newdata$DistHdg / distHdgInverseScale + minTransformedDistHdg
  lp <- fixed_coefs["Intercept"] +
        fixed_coefs["Hedge"] * newdata$Hedge +
        fixed_coefs["Grass"] * newdata$Grass +
        fixed_coefs["Tree"] * newdata$Tree +
        fixed_coefs["Roof"] * newdata$Roof +
        fixed_coefs["Artificial.Surface"] * newdata$Artificial.Surface +
        fixed_coefs["Bush"] * newdata$Bush +
        fixed_coefs["DistHdg"] * newdata$DistHdg +
        fixed_coefs["DistRf"] * newdata$DistRf
  plogis(lp)
}

scoresLP <- predictFunctionForUse(gardenNow)

scoresLP <- ifelse( gardenNow$Roof == 1, NA, scoresLP )

write.csv(scoresLP, output_file, row.names=FALSE )# Append the model coefficients to the output file

fixed_coefs_string <- paste(fixed_coefs, collapse = " ")
cat("above_model_predictions_below_model_coefficients\n", fixed_coefs_string, file = output_file, append = TRUE)

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
