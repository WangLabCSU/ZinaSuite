#' Load Built-in Data
#'
#' @description
#' Load built-in datasets from the ZinaSuite package or download from remote
#' sources if not available locally.
#'
#' @param data_name Name of the dataset to load. Options include:
#'   - "tcga_gtex": TCGA and GTEx sample information
#'   - "tcga_clinical": TCGA clinical data
#'   - "tcga_surv": TCGA survival data
#'   - "pcawg_info": PCAWG sample information
#'   - "ccle_info": CCLE sample information
#'   - "tcga_TIL": Tumor-infiltrating lymphocyte data
#'   - "tcga_tmb": Tumor mutation burden data
#'   - "tcga_MSI": Microsatellite instability data
#'   - "tcga_stemness": Tumor stemness score data
#'   - "tcga_immune": Immune cell infiltration data
#' @param force_download Force re-download even if data exists locally
#' @return Data frame containing the requested data
#' @export
#'
#' @examples
#' \donttest{
#' # Load TCGA/GTEx sample information
#' sample_info <- load_data("tcga_gtex")
#' head(sample_info)
#'
#' # Load survival data
#' surv_data <- load_data("tcga_surv")
#' head(surv_data)
#'
#' # Load TIL data
#' til_data <- load_data("tcga_TIL")
#' head(til_data)
#' }
load_data <- function(data_name = c("tcga_gtex", "tcga_clinical", "tcga_clinical_fine",
                                     "tcga_surv", "tcga_purity", "tcga_subtypes",
                                     "TCGA.organ", "pcawg_info", "pcawg_info_fine",
                                     "pcawg_purity", "ccle_info", "ccle_info_fine",
                                     "ccle_absolute", "tcga_TIL", "tcga_tmb",
                                     "tcga_MSI", "tcga_stemness", "tcga_genome_instability"),
                      force_download = FALSE) {
  data_name <- match.arg(data_name)

  # Check if data is available in package
  if (exists(data_name, envir = .GlobalEnv)) {
    return(get(data_name, envir = .GlobalEnv))
  }

  # Try to load from package data
  data_path <- system.file("data", paste0(data_name, ".rda"), package = "ZinaSuite")

  if (file.exists(data_path) && !force_download) {
    env <- new.env()
    load(data_path, envir = env)
    return(env[[data_name]])
  }

  # Try to load from remote
  remote_data <- try_load_remote_data(data_name, force_download)
  if (!is.null(remote_data)) {
    return(remote_data)
  }

  stop("Data '", data_name, "' not found. Please ensure the package is properly installed.")
}

#' Try to Load Remote Data
#'
#' @description
#' Attempt to download data from remote sources.
#'
#' @param data_name Name of the dataset
#' @param force_download Force re-download
#' @return Data frame or NULL if not available
#' @keywords internal
try_load_remote_data <- function(data_name, force_download = FALSE) {
  # Define remote data URLs
  remote_urls <- list(
    tcga_TIL = "https://zenodo.org/record/xxx/files/tcga_TIL.rda",
    tcga_tmb = "https://zenodo.org/record/xxx/files/tcga_tmb.rda",
    tcga_MSI = "https://zenodo.org/record/xxx/files/tcga_MSI.rda",
    tcga_stemness = "https://zenodo.org/record/xxx/files/tcga_stemness.rda",
    tcga_immune = "https://zenodo.org/record/xxx/files/tcga_immune.rda"
  )

  if (!data_name %in% names(remote_urls)) {
    return(NULL)
  }

  # Check local cache
  cache_dir <- rappdirs::user_cache_dir("ZinaSuite")
  cache_file <- file.path(cache_dir, paste0(data_name, ".rda"))

  if (file.exists(cache_file) && !force_download) {
    env <- new.env()
    load(cache_file, envir = env)
    return(env[[data_name]])
  }

  # Download from remote
  url <- remote_urls[[data_name]]

  tryCatch({
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    utils::download.file(url, cache_file, mode = "wb", quiet = TRUE)

    env <- new.env()
    load(cache_file, envir = env)
    return(env[[data_name]])
  }, error = function(e) {
    warning("Failed to download remote data: ", conditionMessage(e))
    NULL
  })
}

#' Query Purity Data
#'
#' @description
#' Query tumor purity data for TCGA samples.
#'
#' @param source Data source (default: "tcga")
#' @return Named vector of purity values
#' @export
#'
#' @examples
#' \donttest{
#' # Query tumor purity
#' purity <- query_purity(source = "tcga")
#' summary(purity)
#' }
query_purity <- function(source = "tcga") {
  # Try to load from built-in data or remote
  tryCatch({
    purity_data <- load_data("tcga_purity")
    setNames(purity_data$purity, purity_data$sample)
  }, error = function(e) {
    warning("Purity data not available: ", conditionMessage(e))
    NULL
  })
}

#' Query TMB Data
#'
#' @description
#' Query tumor mutation burden (TMB) data.
#'
#' @param source Data source (default: "tcga")
#' @return Named vector of TMB values
#' @export
#'
#' @examples
#' \donttest{
#' # Query TMB
#' tmb <- query_tmb(source = "tcga")
#' summary(tmb)
#' }
query_tmb <- function(source = "tcga") {
  tryCatch({
    tmb_data <- load_data("tcga_tmb")
    setNames(tmb_data$tmb, tmb_data$sample)
  }, error = function(e) {
    warning("TMB data not available: ", conditionMessage(e))
    NULL
  })
}

#' Query MSI Data
#'
#' @description
#' Query microsatellite instability (MSI) data.
#'
#' @param source Data source (default: "tcga")
#' @return Named vector of MSI scores
#' @export
#'
#' @examples
#' \donttest{
#' # Query MSI
#' msi <- query_msi(source = "tcga")
#' table(msi)
#' }
query_msi <- function(source = "tcga") {
  tryCatch({
    msi_data <- load_data("tcga_MSI")
    setNames(msi_data$msi, msi_data$sample)
  }, error = function(e) {
    warning("MSI data not available: ", conditionMessage(e))
    NULL
  })
}

#' Query Stemness Data
#'
#' @description
#' Query tumor stemness score data.
#'
#' @param source Data source (default: "tcga")
#' @return Named vector of stemness scores
#' @export
#'
#' @examples
#' \donttest{
#' # Query stemness
#' stemness <- query_stemness(source = "tcga")
#' summary(stemness)
#' }
query_stemness <- function(source = "tcga") {
  tryCatch({
    stemness_data <- load_data("tcga_stemness")
    setNames(stemness_data$stemness, stemness_data$sample)
  }, error = function(e) {
    warning("Stemness data not available: ", conditionMessage(e))
    NULL
  })
}

#' Query TIL Data
#'
#' @description
#' Query tumor-infiltrating lymphocyte (TIL) data.
#'
#' @param cell_type Cell type to query (default: NULL for all)
#' @param source Data source (default: "tcga")
#' @return Data frame of TIL fractions
#' @export
#'
#' @examples
    #' \donttest{
    #' # Query all TIL data
    #' til <- query_til(source = "tcga")
    #' head(til)
    #' }
query_til <- function(cell_type = NULL, source = "tcga") {
  tryCatch({
    til_data <- load_data("tcga_TIL")

    if (!is.null(cell_type)) {
      if (cell_type %in% colnames(til_data)) {
        setNames(til_data[[cell_type]], til_data$sample)
      } else {
        warning("Cell type '", cell_type, "' not found in TIL data")
        NULL
      }
    } else {
      til_data
    }
  }, error = function(e) {
    warning("TIL data not available: ", conditionMessage(e))
    NULL
  })
}

#' Query Immune Infiltration Data
#'
#' @description
#' Query immune cell infiltration data.
#'
#' @param cell_type Immune cell type to query
#' @param source Data source (default: "tcga")
#' @return Named vector of infiltration scores
#' @export
#'
#' @examples
#' \donttest{
#' # Query CD8 T cell infiltration
#' cd8 <- query_immune_infiltration("CD8_T_cell", source = "tcga")
#' summary(cd8)
#' }
query_immune_infiltration <- function(cell_type = NULL, source = "tcga") {
  tryCatch({
    immune_data <- load_data("tcga_immune")

    if (!is.null(cell_type)) {
      if (cell_type %in% colnames(immune_data)) {
        setNames(immune_data[[cell_type]], immune_data$sample)
      } else {
        warning("Cell type '", cell_type, "' not found in immune data")
        NULL
      }
    } else {
      immune_data
    }
  }, error = function(e) {
    warning("Immune infiltration data not available: ", conditionMessage(e))
    NULL
  })
}
