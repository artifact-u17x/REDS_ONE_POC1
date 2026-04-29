setwd("C:/Users/benke/Documents/gardenREDSmodelUpdateCode/searchingForOtherREDS")

tools <- read.table("toolsDatabase.tsv",
                    sep = "\t",
                    header = TRUE,
                    quote = "\"",
                    stringsAsFactors = FALSE)

nrow(tools)

table(tools$source)
sum(table(tools$source))

tools <- tools[tools$tool_found=="YES",]

nrow(tools)

# A version without the long columns
tools_s <- subset(tools, select = -c(url, article_or_resources))

# And with abbreviated names (only if longer than 30 characters)
tools_s$name <- ifelse(nchar(tools_s$name) > 30, 
                       abbreviate(tools_s$name, minlength = 30),
                       tools_s$name)

# EDS "YES" is coded implicit or explicit, according to whether an EDS tool claims or does not claim to be EDS
# We remove this nuance here
tools_s$eds <- ifelse( substr(tools$eds,1,3)=="YES", "YES", tools_s$eds )

# What EDS is there?
tools_s[substr(tools$eds,1,3)=="YES",]

# What is there with model update from user data?
tools_s[tools$update=="YES",]

# There is nothing with both
tools_s[substr(tools$eds,1,3)=="YES" & tools$update=="YES",]

process_responses <- function(x) {
  # Count each response type
  n_yes <- sum(x == "YES")
  n_no <- sum(x == "NO")
  n_unknown <- sum(x == "UNKNOWN")
  n_known <- n_yes + n_no
  n_total <- length(x)

  # Calculate percentages
  pct_known <- (n_known / n_total) * 100
  pct_yes_of_known <- if (n_known > 0) (n_yes / n_known) * 100 else NA

  # Print results
  cat("% known:", round(pct_known, 0), "\n")
  cat("% of known that are YES:", round(pct_yes_of_known, 0), "\n")
}

process_responses(tools_s$networked)
process_responses(tools_s$model)
process_responses(tools_s$eds)
process_responses(tools_s$update)
process_responses(tools_s$environmental_conditions)
