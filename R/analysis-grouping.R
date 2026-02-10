#' Query TCGA Sample Groups
#'
#' @description
#' Group TCGA/PCAWG/CCLE samples by built-in phenotypes or custom criteria.
#' Supports filtering, merging operations, and custom group definitions.
#'
#' @param database Data source: "tcga", "pcawg", or "ccle"
#' @param cancer Cancer cohort(s) to include. If NULL, includes all cancers.
#' @param custom Custom phenotype data frame to merge with built-in data
#' @param group Target grouping variable name (e.g., "Gender", "Stage_ajcc")
#' @param filter_by List of filtering conditions. Each condition is a vector
#'   of format: c("variable", values, "operator"). Operators:
#'   - "+": include samples with these values
#'   - "-": exclude samples with these values
#'   - ">": greater than (numeric)
#'   - "<": less than (numeric)
#'   - "%>": greater than percentile (numeric)
#'   - "%<": less than percentile (numeric)
#' @param filter_id Vector of sample IDs to filter (direct inclusion)
#' @param merge_by List defining custom group merging. For categorical:
#'   list("Group1" = c("value1", "value2"), "Group2" = c("value3"))
#'   For numeric: list("Low" = c(0, 50), "High" = c(50, NA))
#' @param merge_quantile Whether to interpret numeric merge_by as quantiles
#' @param return_all Whether to return all phenotype columns
#'
#' @return A list with:
#'   - data: Data frame with Sample, Cancer, and group columns
#'   - stat: Summary statistics of the grouping
#' @export
#'
#' @importFrom stats quantile
#' @importFrom rlang .data
#'
#' @examples
#' \donttest{
#' # Basic grouping by gender across all cancers
#' result <- query_tcga_group(group = "Gender")
#' head(result$data)
#'
#' # Group by stage for specific cancer
#' result <- query_tcga_group(
#'   database = "tcga",
#'   cancer = "BRCA",
#'   group = "Stage_ajcc"
#' )
#'
#' # Filter and group
#' result <- query_tcga_group(
#'   cancer = "BRCA",
#'   group = "Stage_ajcc",
#'   filter_by = list(
#'     c("Code", c("TP"), "+"),
#'     c("Stage_ajcc", c(NA), "-")
#'   )
#' )
#'
#' # Numeric filtering with percentile
#' result <- query_tcga_group(
#'   cancer = "BRCA",
#'   group = "Age",
#'   filter_by = list(c("Age", c(0.5), "%>"))
#' )
#'
#' # Merge groups (categorical)
#' result <- query_tcga_group(
#'   cancer = "BRCA",
#'   group = "Stage_ajcc",
#'   merge_by = list(
#'     "Early" = c("Stage I"),
#'     "Late" = c("Stage II", "Stage III", "Stage IV")
#'   )
#' )
#'
#' # Merge groups (numeric with percentiles)
#' result <- query_tcga_group(
#'   cancer = "BRCA",
#'   group = "Age",
#'   merge_quantile = TRUE,
#'   merge_by = list(
#'     "Young" = c(0, 0.5),
#'     "Old" = c(0.5, 1)
#'   )
#' )
#' }
query_tcga_group <- function(database = c("tcga", "pcawg", "ccle"),
                             cancer = NULL,
                             custom = NULL,
                             group = "Gender",
                             filter_by = NULL,
                             filter_id = NULL,
                             merge_by = NULL,
                             merge_quantile = FALSE,
                             return_all = FALSE) {
  database <- match.arg(database)

  # Step 1: Load built-in phenotype data
  meta_data <- switch(database,
    "tcga" = load_data("tcga_clinical_fine"),
    "pcawg" = load_data("pcawg_info_fine"),
    "ccle" = load_data("ccle_info_fine")
  )

  # Standardize column names
  colnames(meta_data)[1:2] <- c("Sample", "Cancer")

  # Step 2: Merge custom phenotype data if provided
  if (!is.null(custom)) {
    dup_names <- intersect(colnames(meta_data)[-1:-2], colnames(custom)[-1])
    if (length(dup_names) > 0) {
      meta_data <- meta_data[, !colnames(meta_data) %in% dup_names]
    }
    colnames(custom)[1] <- "Sample"
    meta_data <- dplyr::inner_join(meta_data, custom, by = "Sample")
  }

  # Step 3: Filter by cancer type(s)
  if (is.null(cancer)) {
    cancer <- unique(meta_data$Cancer)
  }
  meta_data_sub <- meta_data[meta_data$Cancer %in% cancer, ]

  if (nrow(meta_data_sub) == 0) {
    stop("No samples found for specified cancer(s): ", paste(cancer, collapse = ", "))
  }

  # Check if grouping variable exists
  if (!group %in% colnames(meta_data_sub)) {
    available <- setdiff(colnames(meta_data_sub), c("Sample", "Cancer"))
    stop("Group variable '", group, "' not found. Available: ", paste(available, collapse = ", "))
  }

  # Step 4: Apply filters
  if (!is.null(filter_by)) {
    for (i in seq_along(filter_by)) {
      filter_var <- trimws(filter_by[[i]][1])
      filter_op <- trimws(filter_by[[i]][3])
      filter_vals <- strsplit(filter_by[[i]][2], "\\|")[[1]]
      filter_vals <- trimws(filter_vals)
      filter_vals[filter_vals == "NA"] <- NA

      if (!filter_var %in% colnames(meta_data_sub)) {
        warning("Filter variable '", filter_var, "' not found, skipping")
        next
      }

      meta_data_sub <- apply_filter(
        meta_data_sub, filter_var, filter_vals, filter_op
      )
    }
  }

  # Step 5: Filter by sample IDs
  if (!is.null(filter_id)) {
    meta_data_sub <- meta_data_sub[meta_data_sub$Sample %in% filter_id, ]
  }

  if (nrow(meta_data_sub) == 0) {
    warning("No samples remaining after filtering")
    return(list(data = data.frame(), stat = NULL))
  }

  # Step 6: Merge groups
  if (!is.null(merge_by)) {
    meta_data_sub <- merge_groups(
      meta_data_sub, group, merge_by, merge_quantile
    )
  }

  # Remove rows with NA in group variable
  meta_data_sub <- meta_data_sub[!is.na(meta_data_sub[[group]]), ]

  if (return_all) {
    return(list(data = meta_data_sub, stat = summary(meta_data_sub)))
  }

  # Return simplified result
  result_cols <- c("Sample", "Cancer", group)
  result_data <- meta_data_sub[, result_cols, drop = FALSE]
  result_data <- result_data[, !duplicated(colnames(result_data))]

  # Calculate statistics
  group_stats <- table(result_data[[group]], useNA = "ifany")

  list(
    data = result_data,
    stat = group_stats
  )
}

#' Apply Filter to Metadata
#'
#' @param data Data frame
#' @param var Variable name
#' @param values Filter values
#' @param operator Filter operator
#' @return Filtered data frame
#' @keywords internal
apply_filter <- function(data, var, values, operator) {
  switch(operator,
    "+" = {
      # Include values
      if (any(is.na(values))) {
        data[data[[var]] %in% values | is.na(data[[var]]), ]
      } else {
        data[data[[var]] %in% values, ]
      }
    },
    "-" = {
      # Exclude values
      if (any(is.na(values))) {
        data[!(data[[var]] %in% values | is.na(data[[var]])), ]
      } else {
        data[!data[[var]] %in% values, ]
      }
    },
    ">" = {
      # Greater than (absolute)
      thresh <- as.numeric(values[1])
      data[data[[var]] > thresh, ]
    },
    "<" = {
      # Less than (absolute)
      thresh <- as.numeric(values[1])
      data[data[[var]] < thresh, ]
    },
    "%>" = {
      # Greater than (percentile by cancer)
      thresh <- as.numeric(values[1])
      data <- data %>%
        dplyr::group_by(Cancer) %>%
        dplyr::filter(.data[[var]] > quantile(.data[[var]], thresh, na.rm = TRUE)) %>%
        dplyr::ungroup()
      as.data.frame(data)
    },
    "%<" = {
      # Less than (percentile by cancer)
      thresh <- as.numeric(values[1])
      data <- data %>%
        dplyr::group_by(Cancer) %>%
        dplyr::filter(.data[[var]] < quantile(.data[[var]], thresh, na.rm = TRUE)) %>%
        dplyr::ungroup()
      as.data.frame(data)
    },
    {
      warning("Unknown operator '", operator, "', skipping filter")
      data
    }
  )
}

#' Merge Groups
#'
#' @param data Data frame
#' @param group Group variable name
#' @param merge_by Merge definition
#' @param merge_quantile Whether to use quantiles
#' @return Data frame with merged groups
#' @keywords internal
merge_groups <- function(data, group, merge_by, merge_quantile) {
  var_data <- data[[group]]

  if (is.numeric(var_data)) {
    # Numeric variable
    data[[group]] <- merge_numeric_groups(
      var_data, data$Cancer, merge_by, merge_quantile
    )
  } else {
    # Categorical variable
    data[[group]] <- merge_categorical_groups(var_data, merge_by)
  }

  data
}

#' Merge Numeric Groups
#'
#' @param values Numeric values
#' @param cancers Cancer types
#' @param merge_by Merge definition
#' @param use_quantile Whether to use quantiles
#' @return Character vector of merged group labels
#' @keywords internal
merge_numeric_groups <- function(values, cancers, merge_by, use_quantile) {
  result <- rep(NA_character_, length(values))

  if (use_quantile) {
    # Process by cancer type
    for (cancer in unique(cancers)) {
      idx <- cancers == cancer
      cancer_values <- values[idx]

      for (group_name in names(merge_by)) {
        bounds <- merge_by[[group_name]]
        if (is.na(bounds[1])) bounds[1] <- 0
        if (is.na(bounds[2])) bounds[2] <- 1

        q_low <- quantile(cancer_values, bounds[1], na.rm = TRUE)
        q_high <- quantile(cancer_values, bounds[2], na.rm = TRUE)

        in_group <- cancer_values >= q_low & cancer_values <= q_high
        result[idx][in_group] <- group_name
      }
    }
  } else {
    # Absolute values
    for (group_name in names(merge_by)) {
      bounds <- merge_by[[group_name]]
      if (is.na(bounds[1])) bounds[1] <- min(values, na.rm = TRUE)
      if (is.na(bounds[2])) bounds[2] <- max(values, na.rm = TRUE)

      in_group <- values >= bounds[1] & values <= bounds[2]
      result[in_group] <- group_name
    }
  }

  result
}

#' Merge Categorical Groups
#'
#' @param values Character values
#' @param merge_by Merge definition
#' @return Character vector of merged group labels
#' @keywords internal
merge_categorical_groups <- function(values, merge_by) {
  result <- rep(NA_character_, length(values))

  for (group_name in names(merge_by)) {
    group_values <- merge_by[[group_name]]
    result[values %in% group_values] <- group_name
  }

  result
}

#' Get Available Grouping Variables
#'
#' @description
#' List all available grouping variables for a database.
#'
#' @param database Data source: "tcga", "pcawg", or "ccle"
#' @return Character vector of available grouping variables
#' @export
#'
#' @examples
#' \donttest{
#' # Get TCGA grouping variables
#' vars <- get_grouping_variables("tcga")
#' head(vars)
#'
#' # Get PCAWG variables
#' pcawg_vars <- get_grouping_variables("pcawg")
#' }
get_grouping_variables <- function(database = c("tcga", "pcawg", "ccle")) {
  database <- match.arg(database)

  meta_data <- switch(database,
    "tcga" = load_data("tcga_clinical_fine"),
    "pcawg" = load_data("pcawg_info_fine"),
    "ccle" = load_data("ccle_info_fine")
  )

  # Exclude Sample and Cancer columns
  setdiff(colnames(meta_data), c("Sample", "Cancer"))
}
