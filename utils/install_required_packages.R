# =====================================================================
# BioSeq Explorer package healing
# =====================================================================

bioseq_package_log <- function(message, level = "INFO") {
  dir.create("logs", showWarnings = FALSE, recursive = TRUE)
  line <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, message)
  try(cat(line, "\n", file = file.path("logs", "package_install.log"), append = TRUE), silent = TRUE)
  invisible(line)
}

bioseq_notify <- function(message, type = "message") {
  if (requireNamespace("shiny", quietly = TRUE) && shiny::isRunning()) {
    try(shiny::showNotification(message, type = type, duration = 8), silent = TRUE)
  }
  bioseq_package_log(message, if (identical(type, "error")) "ERROR" else "INFO")
}

# Helper to detect RSPM on Linux to avoid source compilation errors
get_linux_binary_repo <- function() {
  default_repo <- "https://cloud.r-project.org"
  if (Sys.info()["sysname"] != "Linux") return(default_repo)
  
  if (file.exists("/etc/os-release")) {
    lines <- readLines("/etc/os-release", warn = FALSE)
    id <- ""
    codename <- ""
    for (line in lines) {
      if (startsWith(line, "ID=")) {
        id <- gsub('^ID="|"$', "", substring(line, 4))
      }
      if (startsWith(line, "VERSION_CODENAME=")) {
        codename <- gsub('^VERSION_CODENAME="|"$', "", substring(line, 18))
      }
      if (startsWith(line, "UBUNTU_CODENAME=")) {
        codename <- gsub('^UBUNTU_CODENAME="|"$', "", substring(line, 17))
      }
    }
    
    if (id %in% c("ubuntu", "debian") && nzchar(codename)) {
      return(sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", codename))
    }
  }
  return(default_repo)
}

# Setup package installation configurations
configure_install_options <- function() {
  has_compiler <- TRUE
  if (Sys.info()["sysname"] %in% c("Windows", "Darwin") && nchar(Sys.which("make")) == 0) {
    has_compiler <- FALSE
  }
  
  pkg_type_val <- if (Sys.info()["sysname"] == "Linux") {
    "source"
  } else {
    if (has_compiler) "both" else "binary"
  }
  
  compile_from_source_val <- "interactive"
  if (Sys.info()["sysname"] %in% c("Windows", "Darwin")) {
    compile_from_source_val <- if (has_compiler) "always" else "never"
  }
  
  options(
    download.file.method = "libcurl",
    repos = c(CRAN = get_linux_binary_repo()),
    warn = 1,
    timeout = 600,
    install.packages.compile.from.source = compile_from_source_val,
    pkgType = pkg_type_val
  )
  
  pkg_type_val
}

ensure_cran_packages <- function(packages, retries = 2) {
  if (length(packages) == 0) return(invisible(TRUE))
  
  # Configure options dynamically
  pkg_type <- configure_install_options()

  results <- vapply(packages, function(pkg) {
    if (suppressWarnings(requireNamespace(pkg, quietly = TRUE))) return(TRUE)
    ok <- FALSE
    for (attempt in seq_len(retries)) {
      tryCatch({
        utils::install.packages(pkg, dependencies = TRUE, quiet = TRUE, type = pkg_type)
        ok <- suppressWarnings(requireNamespace(pkg, quietly = TRUE))
      }, error = function(e) {
        bioseq_package_log(sprintf("CRAN install failed for %s on attempt %s: %s", pkg, attempt, e$message), "WARN")
      })
      if (ok) break
    }
    if (!ok) bioseq_notify(sprintf("Package '%s' could not be installed. Related features will degrade gracefully.", pkg), "warning")
    ok
  }, logical(1))

  invisible(all(results))
}

ensure_bioc_packages <- function(packages, retries = 2) {
  if (length(packages) == 0) return(invisible(TRUE))
  if (!suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))) {
    ensure_cran_packages("BiocManager", retries = retries)
  }
  if (!suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))) {
    bioseq_notify("BiocManager is unavailable. Bioconductor packages could not be healed automatically.", "warning")
    return(invisible(FALSE))
  }

  # Configure options dynamically
  pkg_type <- configure_install_options()
  options(repos = BiocManager::repositories())

  results <- vapply(packages, function(pkg) {
    if (suppressWarnings(requireNamespace(pkg, quietly = TRUE))) return(TRUE)
    ok <- FALSE
    for (attempt in seq_len(retries)) {
      tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE, type = pkg_type)
        ok <- suppressWarnings(requireNamespace(pkg, quietly = TRUE))
      }, error = function(e) {
        bioseq_package_log(sprintf("Bioconductor install failed for %s on attempt %s: %s", pkg, attempt, e$message), "WARN")
      })
      if (ok) break
    }
    if (!ok) bioseq_notify(sprintf("Bioconductor package '%s' is unavailable. Local fallback code will be used where possible.", pkg), "warning")
    ok
  }, logical(1))

  invisible(all(results))
}

ensure_package_environment <- function() {
  cran <- c(
    "shiny", "bslib", "bsicons", "shinyjs", "xml2", "R6", "httr", "r3dmol",
    "seqinr", "data.table", "ggplot2", "plotly", "DT", "dplyr", "tidyr",
    "stringr", "scales", "echarts4r", "RColorBrewer", "patchwork",
    "shinycssloaders", "shinyWidgets", "memoise", "fs", "readr", "purrr",
    "jsonlite", "htmlwidgets"
  )
  
  needs_pwalign <- getRversion() >= "4.4.0"
  bioc <- c("Biostrings", "biomaRt", "AnnotationHub", "GenomicRanges", "DECIPHER")
  if (needs_pwalign) {
    bioc <- c(bioc, "pwalign")
  }

  bioseq_notify("Checking BioSeq Explorer package environment...", "message")
  ok_cran <- ensure_cran_packages(cran)
  ok_bioc <- ensure_bioc_packages(bioc)
  bioseq_package_log(sprintf("Package environment check complete. cran=%s bioc=%s", ok_cran, ok_bioc))
  invisible(ok_cran && ok_bioc)
}
