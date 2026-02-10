#' Data Cache Management
#'
#' @description
#' Advanced data caching system for ZinaSuite to improve performance
#' and reduce redundant data queries.
#'
#' @export
ZinaCache <- R6::R6Class(
  "ZinaCache",
  
  public = list(
    #' @field cache_dir Directory for disk cache
    cache_dir = NULL,
    
    #' @field max_memory_items Maximum items in memory cache
    max_memory_items = 100,
    
    #' @description
    #' Initialize cache manager
    #' @param cache_dir Cache directory path
    #' @param max_memory_items Maximum memory cache size
    initialize = function(cache_dir = NULL, max_memory_items = 100) {
      if (is.null(cache_dir)) {
        cache_dir <- file.path(tempdir(), "ZinaSuite_cache")
      }
      
      self$cache_dir <- cache_dir
      self$max_memory_items <- max_memory_items
      
      # Create cache directory
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Initialize memory cache
      private$memory_cache <- new.env(parent = emptyenv())
      private$access_times <- new.env(parent = emptyenv())
      
      invisible(self)
    },
    
    #' @description
    #' Store data in cache
    #' @param key Cache key
    #' @param data Data to cache
    #' @param use_disk Whether to cache to disk
    #' @param ttl Time to live in seconds (NULL = no expiration)
    set = function(key, data, use_disk = FALSE, ttl = NULL) {
      # Clean key for filename
      clean_key <- private$clean_key(key)
      
      # Store in memory
      private$memory_cache[[clean_key]] <- data
      private$access_times[[clean_key]] <- Sys.time()
      
      # Manage memory cache size
      private$prune_memory_cache()
      
      # Store to disk if requested
      if (use_disk) {
        file_path <- file.path(self$cache_dir, paste0(clean_key, ".rds"))
        metadata <- list(
          key = key,
          created = Sys.time(),
          ttl = ttl
        )
        saveRDS(list(data = data, metadata = metadata), file_path)
      }
      
      invisible(self)
    },
    
    #' @description
    #' Retrieve data from cache
    #' @param key Cache key
    #' @return Cached data or NULL if not found
    get = function(key) {
      clean_key <- private$clean_key(key)
      
      # Check memory cache first
      if (exists(clean_key, envir = private$memory_cache)) {
        private$access_times[[clean_key]] <- Sys.time()
        return(private$memory_cache[[clean_key]])
      }
      
      # Check disk cache
      file_path <- file.path(self$cache_dir, paste0(clean_key, ".rds"))
      if (file.exists(file_path)) {
        cached <- readRDS(file_path)
        
        # Check TTL
        if (!is.null(cached$metadata$ttl)) {
          age <- as.numeric(difftime(Sys.time(), cached$metadata$created, units = "secs"))
          if (age > cached$metadata$ttl) {
            # Expired
            file.remove(file_path)
            return(NULL)
          }
        }
        
        # Restore to memory cache
        private$memory_cache[[clean_key]] <- cached$data
        private$access_times[[clean_key]] <- Sys.time()
        private$prune_memory_cache()
        
        return(cached$data)
      }
      
      return(NULL)
    },
    
    #' @description
    #' Check if key exists in cache
    #' @param key Cache key
    #' @return TRUE if exists
    has = function(key) {
      clean_key <- private$clean_key(key)
      
      # Check memory
      if (exists(clean_key, envir = private$memory_cache)) {
        return(TRUE)
      }
      
      # Check disk
      file_path <- file.path(self$cache_dir, paste0(clean_key, ".rds"))
      if (file.exists(file_path)) {
        cached <- readRDS(file_path)
        if (!is.null(cached$metadata$ttl)) {
          age <- as.numeric(difftime(Sys.time(), cached$metadata$created, units = "secs"))
          if (age <= cached$metadata$ttl) {
            return(TRUE)
          }
        } else {
          return(TRUE)
        }
      }
      
      return(FALSE)
    },
    
    #' @description
    #' Remove item from cache
    #' @param key Cache key
    remove = function(key) {
      clean_key <- private$clean_key(key)
      
      # Remove from memory
      if (exists(clean_key, envir = private$memory_cache)) {
        rm(list = clean_key, envir = private$memory_cache)
        rm(list = clean_key, envir = private$access_times)
      }
      
      # Remove from disk
      file_path <- file.path(self$cache_dir, paste0(clean_key, ".rds"))
      if (file.exists(file_path)) {
        file.remove(file_path)
      }
      
      invisible(self)
    },
    
    #' @description
    #' Clear all cache
    clear = function() {
      # Clear memory cache
      private$memory_cache <- new.env(parent = emptyenv())
      private$access_times <- new.env(parent = emptyenv())
      
      # Clear disk cache
      files <- list.files(self$cache_dir, pattern = "\\.rds$", full.names = TRUE)
      unlink(files)
      
      message("Cache cleared")
      invisible(self)
    },
    
    #' @description
    #' Get cache statistics
    #' @return List with cache stats
    stats = function() {
      memory_items <- length(ls(private$memory_cache))
      
      disk_files <- list.files(self$cache_dir, pattern = "\\.rds$", full.names = TRUE)
      disk_items <- length(disk_files)
      disk_size <- sum(file.size(disk_files), na.rm = TRUE)
      
      list(
        memory_items = memory_items,
        disk_items = disk_items,
        disk_size_bytes = disk_size,
        disk_size_mb = round(disk_size / 1024 / 1024, 2),
        cache_dir = self$cache_dir
      )
    },
    
    #' @description
    #' List all cached keys
    #' @return Character vector of keys
    keys = function() {
      memory_keys <- ls(private$memory_cache)
      
      disk_files <- list.files(self$cache_dir, pattern = "\\.rds$")
      disk_keys <- sub("\\.rds$", "", disk_files)
      
      unique(c(memory_keys, disk_keys))
    }
  ),
  
  private = list(
    memory_cache = NULL,
    access_times = NULL,
    
    clean_key = function(key) {
      # Remove special characters that might cause issues
      gsub("[^a-zA-Z0-9_-]", "_", key)
    },
    
    prune_memory_cache = function() {
      keys <- ls(private$memory_cache)
      
      if (length(keys) > self$max_memory_items) {
        # Find least recently used
        times <- sapply(keys, function(k) private$access_times[[k]])
        lru_key <- keys[which.min(times)]
        
        # Remove it
        rm(list = lru_key, envir = private$memory_cache)
        rm(list = lru_key, envir = private$access_times)
      }
    }
  )
)

#' Get Global Cache Instance
#'
#' @description
#' Get or create the global ZinaCache instance
#' @return ZinaCache object
#' @export
get_zina_cache <- function() {
  if (!exists(".zina_global_cache", envir = .ZinaSuiteEnv)) {
    cache <- ZinaCache$new()
    assign(".zina_global_cache", cache, envir = .ZinaSuiteEnv)
  }
  get(".zina_global_cache", envir = .ZinaSuiteEnv)
}

#' Clear Global Cache
#'
#' @description
#' Clear the global ZinaCache instance
#' @export
clear_zina_cache <- function() {
  if (exists(".zina_global_cache", envir = .ZinaSuiteEnv)) {
    cache <- get(".zina_global_cache", envir = .ZinaSuiteEnv)
    cache$clear()
  }
}

#' Cache-aware Data Query
#'
#' @description
#' Query data with automatic caching
#' @param query_func Function to execute if cache miss
#' @param cache_key Cache key
#' @param use_cache Whether to use cache
#' @param ttl Cache TTL in seconds
#' @return Query result
#' @export
cached_query <- function(query_func, cache_key, use_cache = TRUE, ttl = 3600) {
  if (!use_cache) {
    return(query_func())
  }
  
  cache <- get_zina_cache()
  
  # Try to get from cache
  result <- cache$get(cache_key)
  if (!is.null(result)) {
    message("Cache hit: ", cache_key)
    return(result)
  }
  
  # Execute query
  message("Cache miss: ", cache_key)
  result <- query_func()
  
  # Store in cache
  cache$set(cache_key, result, use_disk = TRUE, ttl = ttl)
  
  return(result)
}
