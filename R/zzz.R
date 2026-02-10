#' ZinaSuite Package Environment
#'
#' @keywords internal
.ZinaSuiteEnv <- new.env(parent = emptyenv())

#' @importFrom utils capture.output head packageVersion
#' @importFrom stats quantile median mad sd cor.test aov lm
#' @importFrom graphics par
#' @importFrom grDevices dev.off png pdf svg
.onLoad <- function(libname, pkgname) {
  # Initialize package environment
  assign(".zina_cache_manager", NULL, envir = .ZinaSuiteEnv)
  assign(".zina_async_engine", NULL, envir = .ZinaSuiteEnv)
}

.onUnload <- function(libpath) {
  # Clean up async engine if running
  if (exists(".zina_async_engine", envir = .ZinaSuiteEnv)) {
    engine <- get(".zina_async_engine", envir = .ZinaSuiteEnv)
    if (!is.null(engine) && engine$is_active()) {
      engine$stop()
    }
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("ZinaSuite loaded. Use query_molecule() to fetch data from UCSC Xena.")
}
