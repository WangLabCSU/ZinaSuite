# Read the file
lines <- readLines("inst/shinyapp/modules/mod_quick_tcga.R")

# Find the UI function (lines 10-232)
start_line <- 10
end_line <- 232

# Count brackets in the UI function
open_count <- 0
close_count <- 0

for (i in start_line:end_line) {
  line <- lines[i]
  # Count non-string brackets
  # Remove strings first
  line_no_strings <- gsub('"[^"]*"', '', line)
  line_no_strings <- gsub("'[^']*'", '', line_no_strings)

  open_count <- open_count + length(gregexpr("\\(", line_no_strings)[[1]])
  close_count <- close_count + length(gregexpr("\\)", line_no_strings)[[1]])

  # Adjust for -1 values from gregexpr
  if (open_count < 0) open_count <- 0
}

# gregexpr returns -1 for no match, so we need to count properly
open_count <- sum(sapply(lines[start_line:end_line], function(line) {
  line_no_strings <- gsub('"[^"]*"', '', line)
  line_no_strings <- gsub("'[^']*'", '', line_no_strings)
  matches <- gregexpr("\\(", line_no_strings)[[1]]
  sum(matches > 0)
}))

close_count <- sum(sapply(lines[start_line:end_line], function(line) {
  line_no_strings <- gsub('"[^"]*"', '', line)
  line_no_strings <- gsub("'[^']*'", '', line_no_strings)
  matches <- gregexpr("\\)", line_no_strings)[[1]]
  sum(matches > 0)
}))

cat("Open brackets:", open_count, "\n")
cat("Close brackets:", close_count, "\n")
cat("Difference:", open_count - close_count, "\n")

# Try to find the issue by parsing line by line
cat("\nTrying to parse line by line...\n")
for (i in start_line:end_line) {
  snippet <- paste(lines[start_line:i], collapse = "\n")
  result <- tryCatch({
    parse(text = snippet)
    NULL
  }, error = function(e) {
    conditionMessage(e)
  })
  if (!is.null(result)) {
    cat("Error around line", i, ":", result, "\n")
    cat("Line content:", lines[i], "\n")
    break
  }
}
