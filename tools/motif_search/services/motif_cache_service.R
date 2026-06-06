# =====================================================================
# Motif Search Cache Service (Session-isolated & Bounded)
# =====================================================================

# Global fallback for non-session / testing contexts
.global_motif_cache_env <- new.env(parent = emptyenv())
.global_motif_cache_keys <- character()

motif_get_cache_store <- function() {
  session <- shiny::getDefaultReactiveDomain()
  if (!is.null(session)) {
    if (is.null(session$userData$motif_cache_env)) {
      session$userData$motif_cache_env <- new.env(parent = emptyenv())
      session$userData$motif_cache_keys <- character()
    }
    return(list(
      env = session$userData$motif_cache_env,
      get_keys = function() session$userData$motif_cache_keys,
      set_keys = function(k) { session$userData$motif_cache_keys <<- k }
    ))
  }
  list(
    env = .global_motif_cache_env,
    get_keys = function() .global_motif_cache_keys,
    set_keys = function(k) { .global_motif_cache_keys <<- k }
  )
}

motif_cache_key <- function(...) {
  parts <- vapply(list(...), function(x) paste(capture.output(str(x, max.level = 1)), collapse = ""), character(1))
  paste0("motif_", as.character(abs(sum(utf8ToInt(paste(parts, collapse = "::"))))))
}

motif_cache_get <- function(key) {
  store <- motif_get_cache_store()
  if (exists(key, envir = store$env, inherits = FALSE)) {
    return(get(key, envir = store$env))
  }
  NULL
}

motif_cache_set <- function(key, value, max_size = 100) {
  store <- motif_get_cache_store()
  
  # Manage bounded FIFO cache eviction
  keys <- store$get_keys()
  if (!exists(key, envir = store$env, inherits = FALSE)) {
    keys <- c(keys, key)
    if (length(keys) > max_size) {
      evict_key <- keys[1]
      if (exists(evict_key, envir = store$env, inherits = FALSE)) {
        rm(list = evict_key, envir = store$env)
      }
      keys <- keys[-1]
    }
    store$set_keys(keys)
  }
  
  assign(key, value, envir = store$env)
  value
}

motif_cached <- function(key, compute) {
  cached <- motif_cache_get(key)
  if (!is.null(cached)) return(cached)
  motif_cache_set(key, compute())
}

motif_cache_clear <- function() {
  store <- motif_get_cache_store()
  rm(list = ls(envir = store$env), envir = store$env)
  store$set_keys(character())
}
