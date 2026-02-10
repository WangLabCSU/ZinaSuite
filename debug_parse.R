# Read file
lines <- readLines("inst/shinyapp/modules/mod_quick_tcga.R")

# Try to parse the entire file
snippet <- paste(lines, collapse = "\n")
result <- tryCatch({
  parse(text = snippet)
  cat("File parsed successfully!\n")
  NULL
}, error = function(e) {
  msg <- conditionMessage(e)
  cat("Error:", msg, "\n")
  
  # Extract line number from error message
  if (grepl(":\\d+:", msg)) {
    line_num <- as.numeric(gsub(".*:(\\d+):.*", "\\1", msg))
    cat("Error around line", line_num, "\n")
    cat("\nContext:\n")
    start <- max(1, line_num-3)
    end <- min(length(lines), line_num+3)
    for (j in start:end) {
      cat(j, ":", lines[j], "\n")
    }
  }
  msg
})
