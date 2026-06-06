# =====================================================================
# FILE: bootstrap.R — Environment Initialization & Startup
# =====================================================================
#
# PURPOSE:
#   This is the PRIMARY STARTUP ENTRY POINT for the BioSeq-Explorer application.
#   It orchestrates the full bootstrap sequence:
#   1. Checks for and installs renv (package manager)
#   2. Restores locked dependencies from renv.lock
#   3. Falls back to requirements.R if renv fails
#   4. Validates application structure
#   5. Launches the Shiny app in browser
#
# WHY RENV?
#   renv.lock is a "lockfile" that records EXACT package versions used in
#   development. This ensures reproducibility: any developer or user gets the
#   SAME package versions, preventing "works on my machine" bugs.
#   Think of it like npm/yarn for R.
#
# FALLBACK STRATEGY:
#   If renv.lock fails (corrupted, missing), we fall back to requirements.R
#   which performs manual package installation. This ensures the app never
#   fails to start due to environment issues.
#
# DEPENDENCIES:
#   - renv: Package version manager (auto-installed)
#   - renv.lock: Lockfile with exact package versions (required file)
#   - requirements.R: Fallback installation script
#   - app.R: The actual Shiny application to launch
#
# DEPENDENTS:
#   - Called by: User running source("bootstrap.R") in R console
#   - Calls: runApp() to launch the Shiny application
#
# OUTPUT:
#   Displays banner and progress messages, then launches Shiny app in browser.
#
# =====================================================================

# ────────────────────────────────────────────────────────────────
# STARTUP BANNER & TIME TRACKING
# ────────────────────────────────────────────────────────────────
# Display ASCII art branding and welcome message
# Record start time to measure total bootstrap duration

cat("\n")
cat("██████╗ ██╗ ██████╗ ███████╗███████╗ ██████╗ \n")
cat("██╔══██╗██║██╔═══██╗██╔════╝██╔════╝██╔═══██╗\n")
cat("██████╔╝██║██║   ██║███████╗█████╗  ██║   ██║\n")
cat("██╔══██╗██║██║   ██║╚════██║██╔══╝  ██║▄▄ ██║\n")
cat("██████╔╝██║╚██████╔╝███████║███████╗╚██████╔╝\n")
cat("╚═════╝ ╚═╝ ╚═════╝ ╚══════╝╚══════╝ ╚══▀▀═╝ \n")

cat("\n")
cat("BioSeq Explorer — Intelligent Bioinformatics Workspace\n")
cat("Version: 1.0.0\n")
cat("Environment Bootstrap Initializing...\n\n")

# Record exact start time for performance measurement
start_time <- Sys.time()

# ────────────────────────────────────────────────────────────────
# CONSOLE HELPER FUNCTIONS
# ────────────────────────────────────────────────────────────────
# These functions format console output with visual indicators
# for better user experience during bootstrap process

# Line separator for section headers (66 equal signs)
line_sep <- paste(rep("═", 66), collapse = "")

# Display info message with ℹ prefix
msg_info <- function(txt) {
  cat(sprintf("ℹ %s\n", txt))
}

# Display success message with ✓ checkmark
# Used after successful installation or validation
msg_ok <- function(txt) {
  cat(sprintf("✓ %s\n", txt))
}

# Display warning message with ⚠ symbol
# Used when something is suboptimal but not fatal
msg_warn <- function(txt) {
  cat(sprintf("⚠ %s\n", txt))
}

# Display failure message with ✖ symbol
# Used when critical error occurs
msg_fail <- function(txt) {
  cat(sprintf("✖ %s\n", txt))
}

# Display boxed section header
# Used to divide bootstrap into logical steps
section <- function(title) {
  cat("\n")
  cat("╔", line_sep, "╗\n", sep = "")
  cat(sprintf("║ %-66s ║\n", title))
  cat("╚", line_sep, "╝\n", sep = "")
}

# ────────────────────────────────────────────────────────────────
# STEP 1: ENSURE RENV PACKAGE MANAGER
# ────────────────────────────────────────────────────────────────
# renv is the package version manager. It's crucial for reproducibility.
# If renv is missing, we install it from CRAN.
# This is a one-time operation; subsequent runs just restore from lock.

section("STEP 1 — CHECKING PACKAGE MANAGER")

if (!requireNamespace("renv", quietly = TRUE)) {
  # renv not installed: show warning and install it
  msg_warn("renv not detected")
  msg_info("Installing renv from CRAN...")

  tryCatch({
    # Install renv from cloud.r-project.org (fast mirror)
    install.packages(
      "renv",
      repos = "https://cloud.r-project.org"
    )
    msg_ok("renv installed successfully")

  }, error = function(e) {
    # If installation fails, show error details and abort
    msg_fail("Failed to install renv")
    cat("\nError Details:\n")
    print(e)
    stop("Bootstrap aborted.")
  })
} else {
  # renv already installed: confirm and move on
  msg_ok("renv already installed")
}

# Load renv into memory for use in next step
library(renv)

# ────────────────────────────────────────────────────────────────
# STEP 2: RESTORE PROJECT ENVIRONMENT
# ────────────────────────────────────────────────────────────────
# renv.lock file contains a snapshot of all package versions at development time.
# renv::restore() installs all packages at their locked versions.
# This ensures "reproducible" environments: same packages, same versions, same behavior.

section("STEP 2 — RESTORING PROJECT ENVIRONMENT")

if (file.exists("renv.lock")) {
  # Lock file exists: parse it and restore all packages
  msg_info("renv.lock detected")
  msg_info("Restoring project dependencies...")

  restore_ok <- FALSE

  tryCatch({
    # Restore all packages from lock file
    # prompt = FALSE: Don't ask user confirmation, just restore
    renv::restore(prompt = FALSE)

    restore_ok <- TRUE
    msg_ok("All packages restored successfully")

  }, error = function(e) {
    # If renv::restore() fails (corrupted lock, network issues, etc.)
    # Fall back to manual installation via requirements.R
    msg_warn("renv restore failed")
    msg_info("Attempting fallback installation...")

    if (file.exists("requirements.R")) {
      tryCatch({
        # Fallback: manually install packages via requirements.R script
        source("requirements.R")

        restore_ok <- TRUE
        msg_ok("Fallback installation completed")

      }, error = function(e2) {
        # Even fallback failed: show error and stop
        msg_fail("Fallback installation failed")
        cat("\nError Details:\n")
        print(e2)
      })
    } else {
      # No fallback script available: abort
      msg_fail("requirements.R not found")
    }
  })

  # If either renv or fallback failed, abort startup
  if (!restore_ok) {
    stop("Dependency installation failed.")
  }

} else {
  # No renv.lock file: use requirements.R for manual installation
  msg_warn("renv.lock not found")

  if (file.exists("requirements.R")) {
    msg_info("Running requirements.R...")

    tryCatch({
      # Execute manual package installation script
      source("requirements.R")
      msg_ok("requirements.R completed")

    }, error = function(e) {
      # If requirements.R fails, show error and abort
      msg_fail("requirements.R failed")
      cat("\nError Details:\n")
      print(e)
      stop("Bootstrap aborted.")
    })
  } else {
    # Neither renv.lock nor requirements.R found: impossible to proceed
    msg_fail("No dependency configuration found")
    stop("Missing renv.lock and requirements.R")
  }
}

# ────────────────────────────────────────────────────────────────
# STEP 3: VALIDATE APPLICATION STRUCTURE
# ────────────────────────────────────────────────────────────────
# Ensure all required application files exist before launching.
# Missing files would cause cryptic startup errors, so we check upfront.

section("STEP 3 — VALIDATING APPLICATION")

# List of files that MUST exist for app to function
required_files <- c(
  "app.R"  # Entry point for Shiny app
)

# Check which required files are missing
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  # Some required files are missing: report them and abort
  msg_fail("Missing required application files:")

  for (f in missing_files) {
    cat("   • ", f, "\n", sep = "")
  }

  stop("Application validation failed.")
}

# All required files present: confirm and proceed
msg_ok("Application structure validated")

# ────────────────────────────────────────────────────────────────
# STEP 4: LAUNCH APPLICATION
# ────────────────────────────────────────────────────────────────
# All dependencies installed, all files validated.
# Now we launch the actual Shiny application.

section("STEP 4 — LAUNCHING APPLICATION")

# Calculate how long the bootstrap process took
# This is useful for identifying slow steps (slow package downloads, etc.)
elapsed <- round(
  as.numeric(difftime(Sys.time(), start_time, units = "secs")),
  2
)

msg_ok(paste("Bootstrap completed in", elapsed, "seconds"))

cat("\n")
cat("🚀 Launching BioSeq Explorer...\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# ── Launch the Shiny Application ──
# appDir = ".": Current directory contains app.R, ui.R, server.R
# launch.browser = TRUE: Automatically open browser to http://localhost:3838
tryCatch({
  shiny::runApp(
    appDir = ".",
    launch.browser = TRUE
  )
}, error = function(e) {
  # If Shiny app fails to run, show error details
  msg_fail("Application failed to launch")
  cat("\nError Details:\n")
  print(e)
})