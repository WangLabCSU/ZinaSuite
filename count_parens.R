# Count parentheses in mod_quick_tcga.R
lines <- readLines("inst/shinyapp/modules/mod_quick_tcga.R")

# Count only significant parentheses (not in strings)
count_parens <- function(lines) {
  open <- 0
  close <- 0
  in_string <- FALSE
  string_char <- ""
  
  for (line in lines) {
    chars <- strsplit(line, "")[[1]]
    i <- 1
    while (i <= length(chars)) {
      char <- chars[i]
      
      # Handle escape sequences
      if (char == "\\" && i < length(chars)) {
        i <- i + 2
        next
      }
      
      # Handle string boundaries
      if (!in_string && (char == '"' || char == "'")) {
        in_string <- TRUE
        string_char <- char
      } else if (in_string && char == string_char) {
        in_string <- FALSE
      } else if (!in_string) {
        if (char == "(") open <- open + 1
        if (char == ")") close <- close + 1
      }
      i <- i + 1
    }
  }
  
  list(open = open, close = close, diff = open - close)
}

result <- count_parens(lines)
cat("Open parentheses:", result$open, "\n")
cat("Close parentheses:", result$close, "\n")
cat("Difference:", result$diff, "\n")

if (result$diff != 0) {
  cat("\nWARNING: Unbalanced parentheses!\n")
}
