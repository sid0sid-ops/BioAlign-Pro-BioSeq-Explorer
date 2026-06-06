# =====================================================================
# FILE: app.R — Application Entry Point
# =====================================================================
# Load all global configurations, libraries, and module sources
# This ensures every dependency is available before UI/server run
source("global.R")

# Define the user interface structure and layout
# Contains all HTML elements, Bootstrap components, and UI logic
source("ui.R")

# Define the server-side reactive logic and event handlers
# Handles all interactivity, data transformations, and state management
source("server.R")

# Create and launch the Shiny application
# ui: The user interface object created in ui.R
# server: The server function defined in server.R
# This combines the presentation layer with the logic layer
options(shiny.port = 3838, shiny.host = "127.0.0.1")
shinyApp(ui = ui, server = server)

