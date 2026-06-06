# =====================================================================
# FILE: requirements.R — Automatic Package Dependency Management
# =====================================================================
#
# PURPOSE:
#   This script automates the installation of all required R packages for
#   BioSeq-Explorer. It serves as a fallback when renv.lock is unavailable
#   or when users prefer manual installation.
#
#   FEATURES:
#   - Checks package installation status
#   - Compares versions against minimum requirements
#   - Auto-updates outdated packages
#   - Handles Bioconductor packages separately (different repo)
#   - Attempts retries on network failures
#   - Automatically configures Posit Package Manager (RSPM) on Linux for binary installs
#   - Soft-installs optional packages (memes, TFBSTools) without halting startup
#   - Gracefully handles pwalign differences between R < 4.4 and R >= 4.4
#   - Comprehensive diagnostic messages for Windows, macOS, and Linux
#
# =====================================================================

# ────────────────────────────────────────────────────────────────
# 1. DEFINE DEPENDENCIES AND COMPATIBILITY
# ────────────────────────────────────────────────────────────────

# pwalign is only available and required in R >= 4.4.0 (Bioc 3.19+)
needs_pwalign <- getRversion() >= "4.4.0"

# Critical packages that the app absolutely needs to boot
critical_pkgs <- c(
  "shiny", "bslib", "bsicons", "R6", "httr", "r3dmol", "echarts4r", "sass", 
  "shinyjs", "seqinr", "rentrez", "remotes", "pkgbuild", "BiocManager", "xml2",
  "DT", "plotly", "ggseqlogo",
  "BiocGenerics", "Biostrings", "DECIPHER", "universalmotif"
)
if (needs_pwalign) {
  critical_pkgs <- c(critical_pkgs, "pwalign")
}

# Optional packages that add advanced features (e.g. Motif Search discovery)
# but should not crash the app if they fail to compile/load
optional_pkgs <- c(
  "biomaRt", "AnnotationHub", "memes", "TFBSTools"
)

required_pkgs <- c(critical_pkgs, optional_pkgs)

required_min_versions <- list(
  shiny = "1.7.0", bslib = "0.5.0", bsicons = "0.1.0", R6 = "2.5.0", httr = "1.4.0",
  r3dmol = "0.1.0", echarts4r = "0.4.0", sass = "0.4.0", shinyjs = "2.1.0", seqinr = "4.2.0",
  rentrez = "1.2.0", remotes = "2.4.0", pkgbuild = "1.4.0", BiocManager = "1.30.22", xml2 = "1.3.0",
  DT = "0.20", plotly = "4.10.0", ggseqlogo = "0.1"
)
if (needs_pwalign) {
  required_min_versions$pwalign <- "1.0.0"
}

fast_check_packages_installed <- function() {
  for (pkg in critical_pkgs) {
    if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) return(FALSE)
  }
  for (pkg in names(required_min_versions)) {
    min_ver <- required_min_versions[[pkg]]
    installed <- tryCatch(utils::packageVersion(pkg), error = function(e) NULL)
    if (is.null(installed) || installed < package_version(min_ver)) return(FALSE)
  }
  TRUE
}

if (fast_check_packages_installed()) {
  cat("BioSeq Explorer — All critical package dependencies verified.\n\n")
} else {
  # ────────────────────────────────────────────────────────────────
  # 2. DETECT COMPILER AND SYSTEM PLATFORM
  # ────────────────────────────────────────────────────────────────
  
  # Check if a C++ compiler environment is connected (make/gcc/clang)
  has_compiler <- TRUE
  if (Sys.info()["sysname"] %in% c("Windows", "Darwin") && nchar(Sys.which("make")) == 0) {
    has_compiler <- FALSE
  }
  
  # ────────────────────────────────────────────────────────────────
  # 3. DETECT LINUX DISTRIBUTION & CONFIGURE BINARIES (RSPM)
  # ────────────────────────────────────────────────────────────────
  # Linux normally compiles everything from source, which takes hours and fails
  # if headers are missing. Posit Package Manager (RSPM) serves precompiled binaries
  # for specific Linux distributions (Ubuntu/Debian). We auto-detect this here.
  
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
      
      # Posit supports ubuntu and debian precompiled binary sets
      if (id %in% c("ubuntu", "debian") && nzchar(codename)) {
        bin_repo <- sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", codename)
        cat(sprintf("Linux system detected (%s %s). Configuring Posit Binary Repository:\n  %s\n\n", id, codename, bin_repo))
        return(bin_repo)
      }
    }
    return(default_repo)
  }
  
  # ────────────────────────────────────────────────────────────────
  # 4. PACKAGE INSTALLATION OPTIONS
  # ────────────────────────────────────────────────────────────────
  
  # Set package type: Linux only supports "source" (RSPM handles binary conversion under-the-hood)
  # Windows/macOS support "binary" or "both" (prefers source, falls back to binary)
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
  
  # ────────────────────────────────────────────────────────────────
  # 5. PRINT PLATFORM COMPILER WARNINGS
  # ────────────────────────────────────────────────────────────────
  if (Sys.info()["sysname"] == "Windows" && !has_compiler) {
    cat("\n========================================================================\n")
    cat("⚠ WARNING: Rtools is NOT connected to your Windows PATH!\n")
    cat("========================================================================\n")
    cat("BioSeq Explorer will install pre-compiled binary packages instead of compiling from source.\n")
    cat("If package installations fail, please consider installing Rtools.\n")
    cat("========================================================================\n\n")
  }
  
  if (Sys.info()["sysname"] == "Darwin" && !has_compiler) {
    cat("\n========================================================================\n")
    cat("⚠ WARNING: Apple Xcode Command Line Tools are NOT installed!\n")
    cat("========================================================================\n")
    cat("BioSeq Explorer will install pre-compiled binary packages instead of compiling from source.\n")
    cat("If installations fail, please run this command in terminal to install compilers:\n")
    cat("  xcode-select --install\n")
    cat("========================================================================\n\n")
  }
  
  # ────────────────────────────────────────────────────────────────
  # 6. PACKAGE MANIFEST
  # ────────────────────────────────────────────────────────────────
  required <- list(
    list(pkg = "shiny",          min = "1.7.0"),
    list(pkg = "bslib",          min = "0.5.0"),
    list(pkg = "bsicons",        min = "0.1.0"),
    list(pkg = "R6",             min = "2.5.0"),
    list(pkg = "httr",           min = "1.4.0"),
    list(pkg = "r3dmol",         min = "0.1.0"),
    list(pkg = "echarts4r",      min = "0.4.0"),
    list(pkg = "sass",           min = "0.4.0"),
    list(pkg = "shinyjs",        min = "2.1.0"),
    list(pkg = "seqinr",         min = "4.2.0"),
    list(pkg = "rentrez",        min = "1.2.0"),
    list(pkg = "remotes",        min = "2.4.0"),
    list(pkg = "pkgbuild",       min = "1.4.0"),
    list(pkg = "BiocManager",    min = "1.30.22"),
    list(pkg = "xml2",           min = "1.3.0"),
    list(pkg = "DT",             min = "0.20"),
    list(pkg = "plotly",         min = "4.10.0"),
    list(pkg = "ggseqlogo",      min = "0.1")
  )
  
  bioc_packages <- c(
    "BiocGenerics",
    "Biostrings",
    if (needs_pwalign) "pwalign" else NULL,
    "biomaRt",
    "AnnotationHub",
    "GenomicRanges",
    "DECIPHER",
    "universalmotif",
    "memes",
    "TFBSTools"
  )
  bioc_packages <- bioc_packages[!vapply(bioc_packages, is.null, logical(1))]
  
  # Global flag tracking if any CRAN/Bioc CRITICAL package failed to install
  all_ok <- TRUE
  
  # ────────────────────────────────────────────────────────────────
  # 7. INSTALLATION HELPER FUNCTION
  # ────────────────────────────────────────────────────────────────
  install_with_verify <- function(pkg, is_bioc = FALSE) {
    success <- FALSE
    attempts <- 1
    max_attempts <- 2
    
    while (!success && attempts <= max_attempts) {
      tryCatch({
        if (is_bioc) {
          # Bioconductor installation via BiocManager
          BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE, type = pkg_type_val)
        } else {
          # CRAN installation via install.packages()
          install.packages(pkg, dependencies = TRUE, quiet = TRUE, type = pkg_type_val)
        }
        
        # Physically test namespace loading to verify success
        if (suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
          success <- TRUE
        } else {
          stop("Package downloaded but failed to register (possible namespace load failure).")
        }
      }, error = function(e) {
        if (attempts == max_attempts) {
          cat(sprintf("  [✗ FAIL]   %-12s Error: %s\n", pkg, e$message))
          if (pkg %in% optional_pkgs) {
            cat(sprintf("  [INFO]     %-12s is optional. The app will launch without it.\n", pkg))
          } else {
            all_ok <<- FALSE
          }
        } else {
          cat(sprintf("  [RETRY]    %-12s Retrying...\n", pkg))
        }
      })
      attempts <- attempts + 1
    }
    return(success)
  }
  
  # ────────────────────────────────────────────────────────────────
  # 8. RUN INSTALLATIONS
  # ────────────────────────────────────────────────────────────────
  
  # Loop through CRAN packages
  for (item in required) {
    pkg <- item$pkg
    min_ver <- item$min
    
    installed <- tryCatch(utils::packageVersion(pkg), error = function(e) NULL)
    
    if (is.null(installed)) {
      cat(sprintf("  [MISSING]  %-12s → Installing...\n", pkg))
      if (install_with_verify(pkg, is_bioc = FALSE)) {
        cat(sprintf("  [✓ OK]     %-12s installed.\n", pkg))
      }
    } else if (installed < package_version(min_ver)) {
      cat(sprintf("  [UPDATE]   %-12s v%s → Updating to >= %s...\n", pkg, installed, min_ver))
      if (install_with_verify(pkg, is_bioc = FALSE)) {
        cat(sprintf("  [✓ OK]     %-12s updated.\n", pkg))
      }
    } else {
      cat(sprintf("  [✓ OK]     %-12s v%s\n", pkg, installed))
    }
  }
  
  # Loop through Bioconductor packages
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    # Keep repositories synced
    options(repos = BiocManager::repositories())
    
    for (pkg in bioc_packages) {
      installed <- tryCatch(utils::packageVersion(pkg), error = function(e) NULL)
      if (is.null(installed)) {
        cat(sprintf("  [BIOC]     %-12s → Installing...\n", pkg))
        if (install_with_verify(pkg, is_bioc = TRUE)) {
          cat(sprintf("  [✓ OK]     %-12s installed.\n", pkg))
        }
      } else {
        cat(sprintf("  [✓ OK]     %-12s v%s (Bioconductor)\n", pkg, installed))
      }
    }
  }
  
  cat(strrep("─", 50), "\n")
  
  # ────────────────────────────────────────────────────────────────
  # 9. COMPREHENSIVE TROUBLESHOOTING GUIDE FOR COMPILATION ERRORS
  # ────────────────────────────────────────────────────────────────
  if (!all_ok) {
    cat("\n========================================================================\n")
    cat("❌ CRITICAL ERROR: Some critical packages failed to install.\n")
    cat("========================================================================\n\n")
    
    # Platform-specific warning details
    sys_name <- Sys.info()["sysname"]
    if (sys_name == "Windows") {
      cat("🚨 DIAGNOSTIC: You are running Windows and may be missing Rtools compilers.\n")
      cat("✅ FIX: Add Rtools (make) to your Path variable or install Rtools44.\n")
      cat("Expected Path: C:\\rtools44\\usr\\bin\\make.exe\n\n")
    } else if (sys_name == "Linux") {
      cat("🚨 DIAGNOSTIC: You are running Linux and may be missing essential system libraries.\n")
      cat("✅ FIX: Run the following commands in your bash terminal to install dependency headers:\n")
      cat("  sudo apt-get update\n")
      cat("  sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libxt-dev make gcc g++ gfortran\n\n")
    } else if (sys_name == "Darwin") {
      cat("🚨 DIAGNOSTIC: You are running macOS and may be missing Xcode Command Line Tools.\n")
      cat("✅ FIX: Run the following command in terminal to install compilers:\n")
      cat("  xcode-select --install\n\n")
    }
    cat("========================================================================\n")
    cat("Proceeding anyway, hoping existing installs are sufficient.\n\n")
  } else {
    cat("All requirements satisfied. Starting BioSeq Explorer...\n\n")
  }
}
