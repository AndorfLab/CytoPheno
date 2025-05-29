# Run Shiny app
#
# CytoPheno: Intakes post-clustered cytometry data and denotes marker description patterns and descriptive cell type names.
# Copyright (C) 2025 Amanda Tursi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

library(shiny)
source("./R/ui.R")
source("./R/server.R")

# Increase max file upload limit
options(shiny.maxRequestSize = 60 * 1024^2)

shinyApp(
  ui = ui_app,
  server = server_app,
  options = list(launch.browser = TRUE))
