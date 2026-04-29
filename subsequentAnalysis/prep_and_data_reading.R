# It's possible one or two of these libraries are no longer used in released code
library(rstanarm)
library(bayesplot)
library(pROC)
library(scoring)
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggtext)
library(INLA)
library(bayesplot)

old_GHSP_blue <- "#56B4E9"
old_GA_orange <- "#E69F00"
GHSP_blue <- "#4A9FD1"
GA_orange <- "#D18F00"

##############################################################################################
# DATA IMPORT AND PREPARATION FUNCTIONS
##############################################################################################

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

cellingToAnalysisFormat <- function(input, include_distTree = FALSE) {
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
  dF <- dF[ dF$boundary > 0,] # Drop all cells with no parts inside boundary
  dF <- dF[ dF$version == "existing", ] # Drop all planning data
  # Remove duplicates based on guid, cell.x, and cell.y (keeping first occurrence)
  # The only reason for duplicates to occur is that two different sparrow observation areas
  # both occur in the same cell - so for current purposes it is fine to just remove duplicates
  dup_key <- paste(dF$garden.guid, dF$cell.x, dF$cell.y, sep = "_")
  dup_indices <- duplicated(dup_key)
  dF <- dF[!dup_indices, ]
  dFS <- data.frame( id = sub( "-.*", "", dF$garden.guid ) ) # Make more convenient shorter garden ids
  dFS$cell_x <- dF$cell.x
  dFS$cell_y <- dF$cell.y
  dFS$boundary <- dF$boundary
  for (i in 1:length(habitatNamesInFile)) { # Add each variable divided by boundary to dFS
    dFS[[habitatNamesOutput[i]]] <- dF[[habitatNamesInFile[i]]] / dF$boundary
  }
  dFS$bushAndCult <- dFS$cult + dFS$bush # Count cultivated as bush
  # Add sparrows and distances
  dFS$sparrow <- ifelse( dF$sparrows.proportion > 0, 1, 0 )
  dFS$distRf <- dF$building.distance
  dFS$distHdg <- dF$hedge.distance
  # Add distTree if requested
  if (include_distTree) {
    dFS$distTree <- calculate_distTree(dF)
  }
  # Prepare distances for analysis (original two distances)
  if( !include_distTree ) {
    dFS$imputedDistHdg <- is.na(dFS$distHdg)
    dFS$imputedDistRf <- is.na(dFS$distRf)
    distances <- prepDistancesForAnalysis( dFS$distRf, dFS$distHdg, capForPrediction = FALSE )
    dFS$anDistRf <- distances$anDistRf
    dFS$anDistHdg <- distances$anDistHdg
  } else { # For distTree: just replace NA with max, no other transformation
    dFS$anDistTree <- ifelse(is.na(dFS$distTree), max(dFS$distTree, na.rm = TRUE), dFS$distTree)
    dFS$anDistHdg <- ifelse(is.na(dFS$distHdg), max(dFS$distHdg, na.rm = TRUE), dFS$distHdg)
    dFS$anDistRf <- ifelse(is.na(dFS$distRf), max(dFS$distRf, na.rm = TRUE), dFS$distRf)
    dFS$imputedDistTree <- is.na(dFS$distTree)
    dFS$imputedDistHdg <- is.na(dFS$distHdg)
    dFS$imputedDistRf <- is.na(dFS$distRf)
  }
  return( dFS )
}

readGlasgowData <- function() {
  # For convenience we redistribute the GHSP data file
  # This is allowed by the journal:
  # https://royalsociety.org/journals/permissions/ states that
  # "Supplementary material is made available under an open access CC-BY licence, meaning that others are free to share,
  # reuse and build upon the information, as long as they properly credit the original author."
  # We remind that the data comes from
  # Jason Matthiopoulos, Christopher Field, Ross MacLeod; Predicting population change from models based on
  # # habitat availability and utilization. Proc Biol Sci 1 April 2019; 286 (1901): 20182911.
  # https://doi.org/10.1098/rspb.2018.2911
  dFG <- read.csv("../redistributedData/rspb20182911_si_002.csv")
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

merge_glasgow_with_REDS_data <- function( glas, reds, caseControlOffset ) {
  if(nrow(reds) > 0) reds$batch <- "new"
  if(nrow(glas) > 0) glas$batch <- "old"
  if(nrow(glas) > 0) glas$prevalenceOffset <- caseControlOffset
  if(nrow(reds) > 0) reds$prevalenceOffset <- 0
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

list_sorted_files <- function() {
  all_files <- list.files(path = "../fullProofOfConcept/cellingBatches",
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
  gardensExcluded <- readLines("../fullProofOfConcept/gardensExcluded.txt")
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
    if(any(grepl("^79846938", temp_df$garden.guid))) {
      cat("Garden with guid starting 'adb5d46d' is from source:", file, "\n")
    }
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

output_coords <- function(dF, filename = "garden_coordinates.csv") {
  coords_df <- dF[!duplicated(dF$garden.guid), c("garden.guid", "epsg3857.x", "epsg3857.y")]
  write.csv(coords_df, filename, row.names = FALSE)
  return(coords_df)
}
