source("../commonCode.R")
set.seed( 8647 )

caseControlOffset <- -2.056832

list_sorted_files <- function() {
  all_files <- list.files(path = "./cellingBatches", 
                          pattern = "^\\d{4}-\\d{2}-\\d{2}_\\d{6}.*\\.csv$", 
                          full.names = TRUE)
  file_timestamps <- sapply(all_files, function(file) {
    timestamp_str <- gsub("^.*/*(\\d{4}-\\d{2}-\\d{2}_\\d{6}).*\\.csv$", "\\1", file)
    as.POSIXct(timestamp_str, format = "%Y-%m-%d_%H%M%S")
  })
  sorted_files <- all_files[order(file_timestamps, decreasing = TRUE)]
  return(sorted_files)
}

process_files_in_order <- function() {
  gardensExcluded <- readLines("gardensExcluded.txt")
  files <- list_sorted_files()
  dF <- data.frame()
  seen_garden_ids <- character(0)
  for (file in files) {
    first_line <- readLines(file, n = 1, warn = FALSE)
    if (length(first_line) > 0 && grepl("no_gardens", first_line)) {
      cat(file, "contains no gardens, skipping\n")
      next
    }
    timestamp_str <- gsub("^.*/*(\\d{4}-\\d{2}-\\d{2}_\\d{6}).*\\.csv$", "\\1", file)
    batch_id <- gsub("^.*garden-celling-([a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12})\\.csv$", "\\1", file)
    temp_df <- read.csv(file)
    temp_df$timeStamp <- as.POSIXct(timestamp_str, format = "%Y-%m-%d_%H%M%S")
    temp_df$batch <- batch_id    
    new_gardens_df <- temp_df[
      !(temp_df$garden.guid %in% seen_garden_ids) &
      !(temp_df$garden.guid %in% gardensExcluded) &
      temp_df$country == "United Kingdom" &
      temp_df$realness == "real"
    , ]
    seen_garden_ids <- union(seen_garden_ids, temp_df$garden.guid)
    dF <- rbind(dF, new_gardens_df)
  }
  return(dF)
}

make_coefs_string <- function(coef_vector) {
  result <- "fixed_coefs <- c(\n"
  for (i in seq_along(coef_vector)) {
    name <- names(coef_vector)[i]
    value <- coef_vector[i]
    line <- sprintf("  %s = %f", name, value)
    if (i < length(coef_vector)) { # Add comma for all but the last item
      line <- paste0(line, ",")
    }
    result <- paste0(result, line, "\n")
  }
  result <- paste0(result, ")\n")
  return(result)
}

combine_scripts <- function(top_file, bottom_file, coefs_string, output_file) {
  top_content <- readLines(top_file, warn = FALSE)
  bottom_content <- readLines(bottom_file, warn = FALSE)
  final_content <- paste0(
    paste(top_content, collapse = "\n"),
    "\n",
    coefs_string,
    "\n",
    paste(bottom_content, collapse = "\n")
  )
  writeLines(final_content, output_file)
  cat("Successfully created", output_file, "\n")
}

dF <- process_files_in_order()
cat("REDS batches contain", nrow(dF), "rows and", ncol(dF), "columns\n")
cat("Number of unique batches:", length(unique(dF$batch)), "\n")
cat("Number of unique gardens:", length(unique(dF$garden.guid)), "\n")
cat("Gardens are:", unique(dF$garden.guid), "\n")
cat("Top left of data:\n")
print(head(dF[,1:4]))

dF <- cellingToAnalysisFormat(dF)
dFO <- readGlasgowData()
dFBoth <- merge_glasgow_with_REDS_data( dFO, dF, caseControlOffset )

updatedModel <- stan_glmer( formStandardUncorr, data = dFBoth, offset = prevalenceOffset, family=binomial, chains = 2 )
#updatedModel <- readRDS( "updatedModelForTesting.rds" )
filename <- paste0("./savedModels/", format(Sys.time(), "%Y-%m-%d_%H%M%S"), "_modelRun.rds")
saveRDS(updatedModel, filename)

# Get fixed effects
fe <- fixef(updatedModel)

# Create new vector with correct names
new_fixed_coefs <- c(
  Intercept = as.numeric(fe["(Intercept)"]) + caseControlOffset,
  Hedge = as.numeric(fe["hedge"]),
  Grass = as.numeric(fe["grass"]),
  Tree = as.numeric(fe["tree"]),
  Roof = as.numeric(fe["roof"]),
  Artificial.Surface = as.numeric(fe["artif"]),
  Bush = as.numeric(fe["bushAndCult"]),
  DistHdg = as.numeric(fe["anDistHdg"]),
  DistRf = as.numeric(fe["anDistRf"])
)

coefs_string <- make_coefs_string(new_fixed_coefs)
combine_scripts("modelScriptTop.R", "modelScriptBottom.R", coefs_string, "newlyUpdatedModel.R")
