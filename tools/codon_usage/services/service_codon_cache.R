# =====================================================================
# Codon Usage Cache Service (Session-isolated & Bounded)
# =====================================================================

# Global fallback for non-session / testing contexts
.global_codon_cache_env <- new.env(parent = emptyenv())
.global_codon_cache_keys <- character()

codon_get_cache_store <- function() {
  session <- shiny::getDefaultReactiveDomain()
  if (!is.null(session)) {
    if (is.null(session$userData$codon_cache_env)) {
      session$userData$codon_cache_env <- new.env(parent = emptyenv())
      session$userData$codon_cache_keys <- character()
    }
    return(list(
      env = session$userData$codon_cache_env,
      get_keys = function() session$userData$codon_cache_keys,
      set_keys = function(k) { session$userData$codon_cache_keys <<- k }
    ))
  }
  list(
    env = .global_codon_cache_env,
    get_keys = function() .global_codon_cache_keys,
    set_keys = function(k) { .global_codon_cache_keys <<- k }
  )
}

codon_cache_key <- function(sequence, host, genetic_code, window, step, rare_threshold) {
  paste(host, genetic_code, window, step, rare_threshold, digest::digest(sequence, algo = "md5"), sep = "::")
}

codon_cached_analysis <- function(sequence, host, genetic_code = "Standard", window = 30, step = 5, rare_threshold = 0.1, max_size = 100) {
  key <- codon_cache_key(sequence, host, genetic_code, window, step, rare_threshold)
  store <- codon_get_cache_store()
  
  if (exists(key, envir = store$env, inherits = FALSE)) {
    return(get(key, envir = store$env))
  }
  
  ref <- load_host_reference(host)
  analysis <- calculate_codon_metrics(sequence, ref, genetic_code = genetic_code)
  analysis$visualization <- build_codon_visualization_data(analysis, ref)
  analysis$sliding <- calculate_sliding_window(analysis$sequence, ref, window, step)
  
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
  
  assign(key, analysis, envir = store$env)
  analysis
}

codon_cache_clear <- function() {
  store <- codon_get_cache_store()
  rm(list = ls(envir = store$env), envir = store$env)
  store$set_keys(character())
}
