#' ZinaSuite: A Modern R Package for UCSC Xena Data Analysis
#'
#' @description
#' Provides a modernized framework for interactively exploring
#' UCSC Xena datasets. It enables downloading, analyzing and visualizing
#' datasets from UCSC Xena with enhanced performance through mirai-based
#' asynchronous computing and modular Shiny application design.
#'
#' @keywords internal
#' @importFrom stats aggregate aov complete.cases cor cor.test fisher.test
#' @importFrom stats kruskal.test kmeans lm median p.adjust pchisq prcomp quantile
#' @importFrom stats reorder residuals rnorm sd setNames t.test wilcox.test
#' @importFrom utils head read.csv read.delim
#' @importFrom methods new
#' @importFrom R6 R6Class
#' @importFrom digest digest
#' @importFrom lifecycle deprecated
#' @importFrom tibble tibble
#' @importFrom nanonext msleep
#' @importFrom rlang .data :=
#' @importFrom grDevices dev.off pdf png svg
#' @importFrom graphics par
#' @importFrom magrittr %>%
#'
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end

# Global environment for package
.ZinaSuiteEnv <- new.env(parent = emptyenv())

# Cache for data
.zina_data_cache <- new.env(parent = emptyenv())

# Async compute engine
.zina_async_engine <- NULL

#' Get or create global cache
#'
#' @return ZinaCache instance
#' @export
get_zina_cache <- function() {
  if (!exists(".zina_global_cache", envir = .ZinaSuiteEnv)) {
    cache <- ZinaCache$new()
    assign(".zina_global_cache", cache, envir = .ZinaSuiteEnv)
  }
  get(".zina_global_cache", envir = .ZinaSuiteEnv)
}

#' Clear global cache
#'
#' @export
clear_zina_cache <- function() {
  if (exists(".zina_global_cache", envir = .ZinaSuiteEnv)) {
    cache <- get(".zina_global_cache", envir = .ZinaSuiteEnv)
    cache$clear()
    rm(".zina_global_cache", envir = .ZinaSuiteEnv)
  }
  invisible(NULL)
}
