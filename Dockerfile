# =====================================================================
# Dockerfile for BioSeq-Explorer (R Shiny Workstation)
# =====================================================================

# Use official rocker image which pre-installs R and Shiny Server daemon
FROM rocker/shiny:4.3.0

# Install system libraries required by Bioconductor & CRAN packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    zlib1g-dev \
    gfortran \
    make \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables for R package downloads
ENV R_REPOS="https://cloud.r-project.org"
ENV DOWNLOAD_METHOD="libcurl"

# Set working directory to default Shiny Server serving directory
WORKDIR /srv/shiny-server

# Remove default Shiny placeholder files
RUN rm -rf /srv/shiny-server/*

# Copy package lists first to leverage Docker layer caching
COPY requirements.R /srv/shiny-server/requirements.R
COPY utils/install_required_packages.R /srv/shiny-server/utils/install_required_packages.R

# Install CRAN and Bioconductor packages (takes 15-20 mins on first build)
RUN Rscript requirements.R

# Copy the rest of the application files
COPY . /srv/shiny-server/

# Set appropriate permissions for the shiny user
RUN chown -R shiny:shiny /srv/shiny-server

# Expose port 3838 (standard Shiny Server port)
EXPOSE 3838

# Launch the Shiny Server daemon
CMD ["/usr/bin/shiny-server"]
