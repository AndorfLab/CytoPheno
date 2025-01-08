# Run Shiny app

library(shiny)
source("./R/ui.R")
source("./R/server.R")

shinyApp(
  ui = ui_app,
  server = server_app
)
