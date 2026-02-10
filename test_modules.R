library(ZinaSuite)
library(shiny)

# Get all module files
module_files <- list.files(system.file("shinyapp/modules", package = "ZinaSuite"),
                           pattern = "mod_.*\\.R$", full.names = TRUE)

print(paste("Found", length(module_files), "module files"))

# Try to source each module
for (f in module_files) {
  tryCatch({
    source(f, local = TRUE)
    print(paste("Loaded:", basename(f)))
  }, error = function(e) {
    print(paste("Failed:", basename(f), "-", conditionMessage(e)))
  })
}

print("Module loading complete!")
