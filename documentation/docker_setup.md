# Docker Container Setup & Deployment

This document describes how to containerize and deploy the **BioSeq-Explorer** R Shiny application using Docker and Docker Compose.

---

## 1. Containerization Architecture

To ensure the application runs reliably across different host environments, we package R, Shiny Server, required Linux system dependencies, and bioinformatics packages in a single container.

The architecture consists of:
- **Base Image**: `rocker/shiny:4.3.0` (provides R, Shiny Server, and the necessary daemon configurations).
- **System Packages**: Installs SSL, XML parsing, and network downloader libraries needed by the R dependencies.
- **R Package Restore**: Copying `renv.lock` or sourcing `requirements.R` to download, build, and register the identical package workspace.
- **Exposure**: Exposes port `3838` (the Shiny application port).

---

## 2. Docker Files Description

Two core files have been added to the root of the project:
1. **`Dockerfile`**: Sets up the operating system, installs system libraries, compiles R packages, and launches the app.
2. **`docker-compose.yml`**: Simplifies orchestration by mapping ports and mounting directories.

---

## 3. Quick Start Guide

### Step 1: Install Docker
Make sure Docker Desktop is installed and running on your host system:
- **Windows**: Install [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/) (WSL 2 backend recommended).
- **Linux/macOS**: Install Docker and Docker Compose.

### Step 2: Build the Container
Open your terminal in the project root directory and run:
```bash
docker compose build
```
*Note: The first build will take 15–20 minutes because heavy CRAN and Bioconductor packages (such as Biostrings and pwalign) need to download and compile. Subsequent builds are cached and take under 10 seconds.*

### Step 3: Run the Application
Start the container in the background:
```bash
docker compose up -d
```

### Step 4: Access the Workstation
Open your web browser and navigate to:
```text
http://localhost:3838
```

---

## 4. Container Management Commands

### View Logs
To view the live application logs and debug startup errors:
```bash
docker compose logs -f
```

### Stop the Container
To stop the application without deleting data:
```bash
docker compose down
```

### Force a Fresh Rebuild
If you update `requirements.R` or `renv.lock` and need to force R to reinstall packages:
```bash
docker compose build --no-cache
```

---

## 5. R Package Caching in Docker

The `Dockerfile` uses a separate build stage to cache compiled R packages on disk. If the container is destroyed or rebuilt, R will not have to recompile packages from source if they haven't changed.
- System dependencies (like `libxml2-dev` and `libssl-dev`) are cached at the operating system layer.
- The default CRAN repository is set to `https://cloud.r-project.org` to ensure maximum download speed and uptime.
- Libcurl is used for R package downloads inside the container to avoid SSL handshake certificate issues.
