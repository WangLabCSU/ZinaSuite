#' CacheManager R6 Class
#'
#' @description
#' Manages caching of data to improve performance. Supports both memory and disk caching
#' with LRU (Least Recently Used) eviction policy.
#'
#' @examples
#' \dontrun{
#' # Create a cache manager
#' cm <- CacheManager$new(cache_dir = tempdir(), max_size = 1e9)
#'
#' # Store data
#' cm$set("gene_TP53", rnorm(100))
#'
#' # Retrieve data
#' data <- cm$get("gene_TP53")
#'
#' # Check if key exists
#' cm$has("gene_TP53")
#'
#' # Get cache info
#' cm$info()
#'
#' # Clear cache
#' cm$clear()
#' }
#'
#' @export
CacheManager <- R6::R6Class(
  "CacheManager",

  public = list(
    #' @description
    #' Initialize a new CacheManager instance
    #' @param cache_dir Directory for disk cache. If NULL, uses system cache dir.
    #' @param max_size Maximum cache size in bytes (default: 1GB)
    #' @param max_memory_items Maximum number of items to keep in memory
    initialize = function(cache_dir = NULL, max_size = 1e9, max_memory_items = 100) {
      private$cache_dir <- cache_dir %||% rappdirs::user_cache_dir("ZinaSuite")
      private$max_size <- max_size
      private$max_memory_items <- max_memory_items
      private$memory_cache <- new.env(parent = emptyenv())
      private$access_order <- character(0)
      private$ensure_dir()
    },

    #' @description
    #' Get data from cache
    #' @param key Cache key
    #' @return Cached data or NULL if not found
    get = function(key) {
      # Check memory cache first
      if (exists(key, envir = private$memory_cache, inherits = FALSE)) {
        private$update_access_order(key)
        return(get(key, envir = private$memory_cache, inherits = FALSE))
      }

      # Check disk cache
      path <- private$get_cache_path(key)
      if (file.exists(path)) {
        data <- readRDS(path)
        private$put_in_memory(key, data)
        return(data)
      }

      NULL
    },

    #' @description
    #' Store data in cache
    #' @param key Cache key
    #' @param value Data to cache
    #' @param persist Whether to persist to disk (default: TRUE)
    #' @return Invisible self
    set = function(key, value, persist = TRUE) {
      # Store in memory
      private$put_in_memory(key, value)

      # Store on disk if requested
      if (persist) {
        path <- private$get_cache_path(key)
        saveRDS(value, path, compress = "xz")
        private$cleanup_if_needed()
      }

      invisible(self)
    },

    #' @description
    #' Check if key exists in cache
    #' @param key Cache key
    #' @return TRUE if exists, FALSE otherwise
    has = function(key) {
      if (exists(key, envir = private$memory_cache, inherits = FALSE)) {
        return(TRUE)
      }
      file.exists(private$get_cache_path(key))
    },

    #' @description
    #' Remove item from cache
    #' @param key Cache key
    #' @return Invisible self
    remove = function(key) {
      # Remove from memory
      if (exists(key, envir = private$memory_cache, inherits = FALSE)) {
        rm(list = key, envir = private$memory_cache, inherits = FALSE)
      }
      private$access_order <- setdiff(private$access_order, key)

      # Remove from disk
      path <- private$get_cache_path(key)
      if (file.exists(path)) {
        file.remove(path)
      }

      invisible(self)
    },

    #' @description
    #' Clear all cached data
    #' @return Invisible self
    clear = function() {
      # Clear memory
      private$memory_cache <- new.env(parent = emptyenv())
      private$access_order <- character(0)

      # Clear disk
      if (dir.exists(private$cache_dir)) {
        files <- list.files(private$cache_dir, full.names = TRUE)
        unlink(files)
      }

      invisible(self)
    },

    #' @description
    #' Get cache statistics
    #' @return List with cache information
    info = function() {
      files <- list.files(private$cache_dir, full.names = TRUE)
      disk_size <- sum(file.size(files))

      list(
        cache_dir = private$cache_dir,
        memory_items = length(private$access_order),
        disk_files = length(files),
        disk_size = disk_size,
        max_size = private$max_size,
        utilization = disk_size / private$max_size
      )
    }
  ),

  private = list(
    cache_dir = NULL,
    max_size = NULL,
    max_memory_items = NULL,
    memory_cache = NULL,
    access_order = NULL,

    ensure_dir = function() {
      if (!dir.exists(private$cache_dir)) {
        dir.create(private$cache_dir, recursive = TRUE, showWarnings = FALSE)
      }
    },

    get_cache_path = function(key) {
      hash <- digest::digest(key, algo = "xxhash32")
      file.path(private$cache_dir, paste0(hash, ".rds"))
    },

    put_in_memory = function(key, value) {
      # Remove if exists to update position
      if (exists(key, envir = private$memory_cache, inherits = FALSE)) {
        private$access_order <- setdiff(private$access_order, key)
      }

      # Add to memory cache
      assign(key, value, envir = private$memory_cache)
      private$access_order <- c(private$access_order, key)

      # Evict oldest if over limit
      if (length(private$access_order) > private$max_memory_items) {
        oldest <- private$access_order[1]
        rm(list = oldest, envir = private$memory_cache, inherits = FALSE)
        private$access_order <- private$access_order[-1]
      }
    },

    update_access_order = function(key) {
      private$access_order <- c(setdiff(private$access_order, key), key)
    },

    cleanup_if_needed = function() {
      files <- list.files(private$cache_dir, full.names = TRUE)
      total_size <- sum(file.size(files))

      if (total_size > private$max_size) {
        # Sort by access time (oldest first)
        file_info <- file.info(files)
        files <- files[order(file_info$atime)]

        # Remove oldest files until under limit
        for (f in files) {
          if (total_size <= private$max_size * 0.8) break
          total_size <- total_size - file.size(f)
          file.remove(f)
        }
      }
    }
  )
)
