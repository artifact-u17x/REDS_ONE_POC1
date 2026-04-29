start.time <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

gardenNow <- read.csv(input_file) # load provided data

