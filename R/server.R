#' Shiny app server
#'
#' @name server_app
#' @param input Object that contains all the input data sent from the browser, read only
#' @param output Object named according to the output ID
#' @param session  Code run when session starts
#'
#' @import dplyr
#' @import ggplot2
#' @import ggridges
#' @import heatmaply
#' @import httr
#' @import jsonlite
#' @import magrittr
#' @import plotly
#' @import RColorBrewer
#' @import readxl
#' @import reshape2
#' @import scales
#' @import shiny
#' @import DT
#' @import shiny.destroy
#' @import shinybusy
#' @import shinyjs
#' @import shinyvalidate
#' @import stringi
#' @import stringr
#' @import tidyr
#' @import zip
#'
#' @importFrom magrittr %>%
#'
#' @return server
#'

# Cell.Naming: Intakes post-clustered cytometry data and denotes marker description patterns and descriptive cell type names.
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

source("endpoints.R")
source("protein_SPARQL.R")
source("integer_breaks.R")

# Increase max size for uploaded files
options(shiny.maxRequestSize = 60*1024^2)

server_app <- function(input, output, session) {
 
 url1 <-"https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/Kimmey_5hr_stim.csv"
 Kimmey_5hr_stim <- read.csv(url1, header = TRUE)  

 url2 <-"https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/Dusoswa_OMIP_54_markers.csv"
 Dusoswa_OMIP_54_markers <- read.csv(url2, header = TRUE) 

 url3 <- "https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/Lee_AML_cell_types_markers.csv"
 Lee_AML_cell_types_markers <- read.csv(url3, header = TRUE) 
 
  #################################################################
  ##                            Tab 1:                           ##
  ##                            Server                           ##
  ##                    Input expression data                    ##
  #################################################################
  ######~#~# Step 1 - Input expression csv, choose markers, change marker names if desired #~#~######

  # Side panel

  # Download example expression-cluster data
  output$download_example_expression_cluster <- shiny::downloadHandler(
    # File name
    filename = function() {
      c("Kimmey_5hr_stim.csv")
    },
    # Write to csv
    content = function(file) {
      utils::write.csv(Kimmey_5hr_stim, file, row.names = FALSE) 
    }
  )

  # Set seed if the user chooses to do so
  global_seed <- shiny::eventReactive(input$submit_tab1_step1,{
    shiny::req(input$seed_1 == 1)
    input$num_seed_1
  })

  shiny::observe({
    set.seed(global_seed())
  })

  # Reset everything
  shiny::observeEvent(input$reset_tab1_step1, {
    session$reload()
  })

  # Validation rules for text input
  val_cofactor <- shinyvalidate::InputValidator$new()
  val_cofactor$add_rule("type_cofactor", shinyvalidate::sv_numeric())
  val_cofactor$enable()

  val_seed <- shinyvalidate::InputValidator$new()
  val_seed$add_rule("num_seed_1", shinyvalidate::sv_numeric())
  val_seed$enable()

  val_seed <- shinyvalidate::InputValidator$new()
  val_seed$add_rule("num_random_1", shinyvalidate::sv_numeric())
  val_seed$enable()

  # Help buttons for side panel parameters

  shiny::observeEvent(input$help_input_type_1,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Upload Expression Data"),
      HTML("The user is required to upload an input file (CSV, XLS, or XLSX) that includes the cells (as rows) and markers (as columns) used in clustering. The first column, labeled Cluster, should contain the cluster number to which each cell was assigned.
      <br>
      <br>
      For example, if 40 markers were used to define 50,000 cells, the file should contain expression values for 50,000 rows and 41 columns (40 markers + the cluster column).
      <br>
      <br>"),
      tags$div("An example expression input file is available within the application (click 'Download Example' above the file upload box) and on the ", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/2.-Example-Data', "GitHub.", target="_blank")),
      easyClose = TRUE))
  })

  shiny::observeEvent(input$help_pregating,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Type Pre-Gated Markers"),
      HTML("Any markers that were pre-gated should be inputted here. These markers will be excluded from the Median Difference Equation and treated as homogeneous across all clusters.
      <br>
      <br>
      For example, if the data was pre-gated for CD3 positive expression, enter 'CD3+' in the text input box. CD3 will be interpreted as positive for all clusters.
      <br>
      <br>
      Any included pre-gated markers are excluded from Part 1 (Steps 1-3).
      <br>
      <br>
      If using internal references, CD45 will also be excluded from all analysis steps due to inconsistent usage in the CL."),
      easyClose = TRUE))
  })

  shiny::observeEvent(input$help_transform_option,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select Hyperbolic Arcsine (Arcsinh) Transformation"),
      HTML("The arcsinh transformation is widely used in cytometry to bring data closer to a normal distribution while handling zeros and negative values. It is recommended for use before downstream cytometry analysis.
      <br>
      <br>
      If you choose to apply this transformation, you must specify a numeric cofactor. Recommended values are:
      <br>
      <br>
      <li>5 for mass cytometry</li>
      <br>
      <li>15 for fluorescence/conventional flow cytometry</li>
      <br>
      <li>6000 for spectral flow cytometry</li>"),
      easyClose = TRUE))
  })

  shiny::observeEvent(input$help_random_type_1,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select Random Subsample of Data"),
      HTML("Subsampling can be used to reduce the dataset size by selecting a subset of cells. This is particularly useful if your expression data contains hundreds of thousands of cells, as very large files can slow down the tool.
      <br>
      <br>
      If you choose to subsample the data, specify the number of cells (rows) to be included in the subset."),
      easyClose = TRUE))
  })

  shiny::observeEvent(input$help_seed_type_1,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Set Seed"),
      HTML("A seed ensures that results are reproducible when the same dataset and parameters are inputted multiple times. This is especially relevant for subsampling, as the seed will ensure the same cells are selected each time, resulting in consistent outputs.
      <br>
      <br>
      If no seed is specified or if different seeds are used with each submission, slight variations in the results may occur due to the random nature of subsampling."),
      easyClose = TRUE))
  })

  # Main panel

  # Step 1, initial changes when the 'submit' button is hit
  shiny::observeEvent(input$submit_tab1_step1, {

    shinyjs::alert("This multi-step process may take a few minutes. Please have patience.")

    shiny.destroy::removeOutput("app_overview_1")
    shinyjs::show("overview_tab1_step1")
    shiny.destroy::removeOutput(id = "all_sidebar_tab1_step1")
    shinyjs::show(id = "all_sidebar_tab1_step2")

    if (is.null(uploaded_file()) == FALSE) {

      if (is.null(check_input_file_1()) == TRUE) {

        df_expr$data <- uploaded_expression_head()

        # Uploaded expression dataframe
        output$expression_data <- DT::renderDT({
          DT::datatable(
            df_expr$data,
            selection = 'multiple',
            editable = list(target = "cell", disable = list(columns = 0)),
            rownames = FALSE,
            escape = FALSE,
            options = list(pageLength = 75, processing=FALSE)
          )
        })
        shinyjs::show(id = "submit_tab1_step2")
        shinyjs::show("help_tab1_step1")
        shinyjs::show("delete_marker")
      }
    }
  })

  # Step 1, create editable dataframes
  df_expr <- shiny::reactiveValues(data = NULL)

  # Step 1, check if file type of uploaded input file is acceptable
  check_input_file_1 <- shiny::eventReactive(input$submit_tab1_step1, {

    # File inputted
    file_input <- input$upload_expression_csv

    if (endsWith(file_input$datapath, ".csv") | endsWith(file_input$datapath, ".xls") | endsWith(file_input$datapath, ".xlsx")) {
      return(NULL)
    }  else {
      return("a")
    }
  })

  # Step 1, outputs error stating that the input file is in the wrong format
  output$ui_error_check_input_file_1 <- shiny::renderUI({
    shiny::req(input$submit_tab1_step1)
    if (is.null(check_input_file_1()) == FALSE) {
      tags$div(id = "invalid_input_file_1", h4(
        "Error: Invalid uploaded input file format. Please upload a valid comma separated values (CSV) or Excel (XLSX or XLS) file."
      ))
    }})

  # Step 1, read in uploaded file
  uploaded_file <- shiny::eventReactive(input$submit_tab1_step1,{

    if (is.null(check_input_file_1()) == TRUE) {

      file_input <- input$upload_expression_csv

      if (endsWith(file_input$datapath, ".csv")) {

        read_file <- utils::read.csv(file_input$datapath, header = TRUE)

      } else if (endsWith(file_input$datapath, ".xls")) {

        read_file <- readxl::read_excel(file_input$datapath, col_names = TRUE)

      }  else if (endsWith(file_input$datapath, ".xlsx")) {

        read_file <- readxl::read_excel(file_input$datapath, col_names = TRUE)
      }

      if (input$random_1 != 1) {
        read_file <- read_file[sample(1:nrow(read_file), input$num_random_1, replace=FALSE),]
      }

      # Ensure first letter is capitalized and the rest in lower case
      colnames(read_file) <- toupper(colnames(read_file))

      # Check that 1st column name is correct. If not, return a warning and stop.
      if (colnames(read_file)[1] != "CLUSTER") {
        return(NULL)
      }
      return(read_file)
    }
  })

  # Step 1, uploaded input, error message if input has wrong column names
  output$ui_error_column_names_input_1 <- shiny::renderUI({
    shiny::req(input$submit_tab1_step1)

    if (is.null(check_input_file_1()) == TRUE) {
      if (is.null(uploaded_file()) == TRUE) {
        h4(
          "Error: Invalid uploaded input format. The first column should be titled 'Cluster' and contain the cluster each cell was assigned to."
        )
      }
    }
  })

  # Step 1, add in pre-gated markers
  added_MG_markers <- shiny::eventReactive(input$submit_tab1_step1,{

    P <- input$type_MG_markers

    if (!is.null(P)) {
      # Remove spaces
      P <- gsub(" ", "", P)

      if (P != "") {
        new_df <- utils::read.table(text = P, sep = ",")
        new_df <- as.data.frame(t(new_df))

        # Determine if positive or negative
        new_df$Sign[grepl("positive|\\+$", new_df$V1, ignore.case = TRUE)] <-
          "Positive"
        new_df$Sign[grepl("negative|\\-$", new_df$V1, ignore.case = TRUE)] <-
          "Negative"
        new_df$Sign[grepl("high|\\+\\+|\\+\\+\\+$", new_df$V1, ignore.case =
                            TRUE)] <- "High"
        new_df$Sign[grepl("low$", new_df$V1, ignore.case = TRUE)] <-
          "Low"

        # Remove the sign from the marker
        new_df$V1 <-
          gsub("\\+$", "", gsub("\\++$", "",  gsub("\\-$", "", gsub(
            "\\--$", "", gsub(
              "low$",
              "",
              ignore.case = TRUE,
              gsub(
                "high$",
                "",
                ignore.case = TRUE,
                gsub(
                  "positive$",
                  "",
                  ignore.case = TRUE,
                  gsub(
                    "negative$",
                    "",
                    as.character(new_df$V1),
                    ignore.case = TRUE
                  )
                )
              )
            )
          ))))


        get_clusters <- uploaded_file()

        names(get_clusters) <- tolower(names(get_clusters))

        unique_clusters <- unique(get_clusters$cluster)

        names(new_df) <- tolower(names(new_df))

        new_df <- merge(new_df, unique_clusters)

        # New column names
        colnames(new_df) <- c("Marker", "Sign", "Cluster")

        return(new_df)
      } else {
        P <- ""
        return(P)
      }
    } else {
      P <- ""
      return(P)
    }
  })

  # Show markers
  uploaded_expression_head <- shiny::eventReactive(input$submit_tab1_step1,{

    all_markers <- colnames(uploaded_file())

    indx <- grepl('cluster', all_markers, ignore.case=TRUE)

    all_markers <- all_markers[!indx]

    all_markers <- data.frame("Original marker names" = all_markers,
                              "New marker names" = all_markers, check.names = FALSE)

    all_markers <- all_markers[order(all_markers$'Original marker names'),]

    return(all_markers)
  })

  # Step 1, when help button is clicked show message
  shiny::observeEvent(input$help_tab1_step1,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 1<br>Select Markers Used in Analysis"),
      HTML("Before running any analysis, all markers present in the input data are displayed. The user has the option to edit marker names or remove markers from further analysis. It's recommended that only the markers used in clustering be included for cell type matching.
      <br>
      <br>
      <li>Edit Marker Names: Click on the marker in the 'New Marker Names' column to change the name.</li>
      <br>
      <li>Delete Markers: Select the row containing the marker(s) to delete. Multiple rows can be selected at once. Once the desired markers are highlighted, click 'Delete Marker(s)' at the bottom. Note: Deleting markers is irreversible unless you reset the tool via the sidebar.</li>
      <br>
      <br>"),
      tags$div("Additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/3.1.-Input-Expression-Data-and-Default-References', "here.", target="_blank")),
      easyClose = TRUE))
  })

  # Step 1, Option to delete markers
  shiny::observeEvent(input$delete_marker, {
    if (!is.null(input$expression_data_rows_selected)) {
      df_expr$data <<-
        df_expr$data[-as.numeric(input$expression_data_rows_selected), ]
    }
  })

  # Step 1, option to edit (add new marker names) in the table
  shiny::observeEvent(input$expression_data_cell_edit, {
    shiny::req(input$submit_tab1_step1)
    df_expr$data <<-
      DT::editData(df_expr$data,
                   input$expression_data_cell_edit,
                   'expression_data',
                   rownames = FALSE)
  })

  ######~#~# Step 2 - Median Difference Equation, change parameters if need be #~#~######

  # Step 2, side panel

  # Step 2, reset everything
  shiny::observeEvent(input$reset_tab1_step2, {
    session$reload()
  })

  # Step 2, help buttons for side panel parameters
  shiny::observeEvent(input$help_marker_diff_eq_option,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Specify Median Difference Equation Parameters"),
      HTML("After the median difference between the cluster of interest and the reference clusters is calculated, a marker is classified as positive, negative, or null for each cluster based on specific cutoff values.
      <br>
      <br>
      While default parameters are recommended, all cutoff values can be adjusted within the application. Changing parameters in this step will result in live updates to the results, allowing users to explore different thresholds."),
      easyClose = TRUE))
  })

  shiny::observeEvent(input$help_pos_neg,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Marker designated as positive or negative"),
      HTML("A marker is classified as positive, negative, or null for each cluster based on the following criteria:
      <br>
      <br>
      <li>Positive: Value is above a specified threshold (default = 0.25)</li>
      <br>
      <li>Negative: Value is below a specified threshold (default = -0.75)</li>
      <br>
      <li>Null: Values between the positive and negative cutoffs</li>
      "),
      easyClose = TRUE))
  })

  shiny::observeEvent(input$help_all_pos,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Marker designated as completely positive"),
      HTML("A marker is classified as positive across all clusters if:
      <br>
      <br>
      <li>The minimum median across all clusters exceeds a certain threshold (default = 0.2)</li>
      <br>
      <li>The maximum median across all clusters exceeds another value (default = 1)</li>
      <br>
      <li>The overall standard deviation is below a set value (default = 0.7)</li>
      "),
      easyClose = TRUE))
  })

  shiny::observeEvent(input$help_all_null,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Marker designated as completely null"),
      HTML("A marker is classified as null across all clusters if:
      <br>
      <br>
      <li>The minimum median across all clusters is below a certain threshold (default = 0.2)</li>
      <br>
      <li>The maximum median across clusters is below another value (default = 0.7)</li>
      <br>
      <li>The overall standard deviation is below a third threshold (default = 0.5)</li>
      "),
      easyClose = TRUE))
  })

  # Step 2, main panel

  # Step 2, run Median Difference Equation
  shiny::observeEvent(input$submit_tab1_step2, {

    shiny.destroy::removeOutput("help_tab1_step1")
    shiny.destroy::removeOutput("overview_tab1_step1")
    shiny.destroy::removeOutput("delete_marker")
    shiny.destroy::removeOutput("expression_data")
    shinyjs::show("overview_tab1_step2")
    shiny.destroy::removeOutput("submit_tab1_step2")
    shinyjs::show("submit_tab1_step2_5")
    shinyjs::show("marker_diff_eq_heatmap_binary_ex_plot2")

    # Step 2, output text telling the user what cutoff parameters were used
    output$cutoff_parameters <- shiny::renderUI({
      if (input$marker_diff_eq_option == 'default_marker_diff_eq') {

        # Marker designated as positive or negative
        type_positive_cutoff <- 0.25
        type_negative_cutoff <- -0.75

        # Marker designated as completely positive
        type_med_positive_cutoff <- 0.2
        type_max_med_positive_cutoff <- 1
        type_stdr_positive_cutoff <- 0.7

        # Marker designated as completely null
        type_med_negative_cutoff <- 0.2
        type_max_med_negative_cutoff <- 0.7
        type_stdr_negative_cutoff <- 0.5

      } else if (input$marker_diff_eq_option == 'choose_marker_diff_eq') {

        # Marker designated as positive or negative
        type_positive_cutoff <- input$type_positive_cutoff
        type_negative_cutoff <- input$type_negative_cutoff

        # Marker designated as completely positive
        type_med_positive_cutoff <- input$type_med_positive_cutoff
        type_max_med_positive_cutoff <- input$type_max_med_positive_cutoff
        type_stdr_positive_cutoff <- input$type_stdr_positive_cutoff

        # Marker designated as completely null
        type_max_med_negative_cutoff <- input$type_max_med_negative_cutoff
        type_med_negative_cutoff <- input$type_med_negative_cutoff
        type_stdr_negative_cutoff <- input$type_stdr_negative_cutoff
      }

      p(HTML(paste("<b>", "Marker designated as positive or negative:", "</b>",
                   "Positive cutoff = ",  type_positive_cutoff,
                   ", Negative cutoff = ", type_negative_cutoff,
                   "<br/> <b>", "Marker designated as completely positive:", "</b>",
                   "Min median >", type_med_positive_cutoff,
                   ", Max median >", type_max_med_positive_cutoff,
                   ", Standard deviation <", type_stdr_positive_cutoff,
                   "<br/> <b>", "Marker designated as completely null:", "</b>",
                   "Min median <", type_med_negative_cutoff,
                   ", Max median <", type_max_med_negative_cutoff,
                   ", Standard deviation <", type_stdr_negative_cutoff)), style = "font-size:20px; color: #b24bb4; text-align: center")
    })

    output$marker_diff_eq_heatmap_binary_ex_plot = plotly::renderPlotly({
      marker_diff_eq_heatmap_binary_ex()
    })
    output$marker_diff_eq_heatmap_raw_ex_plot = plotly::renderPlotly({
      marker_diff_eq_heatmap_raw_ex()
    })
    output$expr_heatmap_org_ex_plot = plotly::renderPlotly({
      expr_heatmap_org_ex()
    })
    output$expr_density_ex_plot = shiny::renderPlot({
      expr_density_ex()
    })
    output$expr_heatmap_max_ex_plot = plotly::renderPlotly({
      expr_heatmap_max_ex()
    })
    output$expr_heatmap_min_ex_plot = plotly::renderPlotly({
      expr_heatmap_min_ex()
    })
    output$expr_heatmap_std_ex_plot = plotly::renderPlotly({
      expr_heatmap_std_ex()
    })
  })

  # Step 2, Median Difference Equation
  raw_marker_diff_eq_results <- shiny::eventReactive(input$submit_tab1_step2,{

    marker_lookup <- df_expr$data

    org_df <- uploaded_file()

    indx <- grepl('cluster', colnames(org_df), ignore.case=TRUE)

    colnames(org_df)[indx] <- tolower(colnames(org_df[indx]))

    # Change marker names if applicable
    kept_markers <- org_df[, c(marker_lookup$`Original marker names`, "cluster")]

    change_names <- base::match(names(kept_markers), marker_lookup$`Original marker names`)

    names(kept_markers)[stats::na.omit(change_names)] = marker_lookup$`New marker names`[!is.na(change_names)]

    # Arcsinh transform the data if the user checked 'yes'
    if (input$transform_option == 'transform_yes') {
      cofactor <- input$type_cofactor
      kept_markers[, colnames(kept_markers)[colnames(kept_markers) != 'cluster']] <- asinh((kept_markers[, colnames(kept_markers)[colnames(kept_markers) != 'cluster']]) / cofactor)
    }

    # List of each cluster and the expression values
    each_cluster <- split(kept_markers, kept_markers$cluster)

    # Get vector of markers and the number of markers used
    markers_used <- names(each_cluster[[1]])
    markers_used <- markers_used[! markers_used %in% c("cluster", "Cluster", "population", "Population")]
    num_markers <- length(markers_used)

    # Empty vectors
    median_cluster_interest_vector <- c()
    median_cluster_reference_vector <- c()

    # Empty lists
    median_cluster_interest_list <- list()
    median_cluster_reference_list <- list()


    # Loop through each cluster
    for (cluster_num in 1:length(each_cluster)) {

      # Loop through each marker
      for (marker in 1:num_markers) {

        ## Cluster of interest ##

        # Specific cluster
        cluster_interest <- each_cluster[cluster_num]
        cluster_interest_name <- names(cluster_interest)

        # Specific marker in that cluster
        cluster_interest <- lapply(cluster_interest, "[", marker)

        # Put those marker values from that specific cluster into a vector
        cluster_interest <- unlist(cluster_interest, use.names=FALSE)

        # Get median
        median_cluster_interest <- abs(stats::median(cluster_interest))

        # Put into a vector
        median_cluster_interest_vector <- c(median_cluster_interest_vector, median_cluster_interest)

        ## Reference clusters (same marker, all other clusters) ##

        # Every other cluster except cluster of interest
        other_clusters <- each_cluster[-cluster_num]

        # Get marker of interest
        other_clusters <- lapply(other_clusters, "[", marker)

        # Put those marker values from other clusters into a vector
        cluster_reference <- unlist(other_clusters, use.names=FALSE)

        # Get median
        median_cluster_reference <- abs(stats::median(cluster_reference))

        # Put into a vector
        median_cluster_reference_vector <- c(median_cluster_reference_vector, median_cluster_reference)

      }

      # Put all median values (of interest) for each cluster into a list with all cluster IQR
      median_cluster_interest_list[[cluster_interest_name]] <- median_cluster_interest_vector

      # Empty the vector so it can start again for the next cluster
      median_cluster_interest_vector <- c()

      # Put all median values (reference) for each cluster into a list with all cluster IQR
      median_cluster_reference_list[[cluster_interest_name]] <- median_cluster_reference_vector

      # Empty the vector so it can start again for the next cluster
      median_cluster_reference_vector <- c()

    }

    ## Cluster of interest ##

    # Turn list of median (of interest) into dataframe
    median_cluster_interest_df <- as.data.frame(median_cluster_interest_list)

    # Tranpose
    median_cluster_interest_df <- t(median_cluster_interest_df)

    # Add column names
    colnames(median_cluster_interest_df) <- names(each_cluster[[1]])[1:num_markers]

    ## Reference clusters (same marker, all other clusters) ##

    # Turn list of median (reference) into dataframe
    median_cluster_reference_df <- as.data.frame(median_cluster_reference_list)

    # Tranpose
    median_cluster_reference_df <- t(median_cluster_reference_df)

    # Add column names
    colnames(median_cluster_reference_df) <- names(each_cluster[[1]])[1:num_markers]

    ## Get median difference and scale per marker ##

    # Median cluster of interest - median reference clusters
    median_difference <- median_cluster_interest_df - median_cluster_reference_df

    df_scaled <- as.data.frame(median_difference)

    df_scaled <- sapply(df_scaled, scales::rescale, to=c(-1,1))

    df_scaled <- data.frame(Cluster = names(each_cluster), df_scaled)

    return(df_scaled)
  })

  # Step 2, Median Difference Equation, cutoffs
  heatmap_marker_diff_eq_results <- reactive({
  shiny::req(input$submit_tab1_step2)
   
    if (input$marker_diff_eq_option == 'default_marker_diff_eq') {

      # Marker designated as positive or negative
      type_positive_cutoff <- 0.25
      type_negative_cutoff <- -0.75

      # Marker designated as completely positive
      type_med_positive_cutoff <- 0.2
      type_max_med_positive_cutoff <- 1
      type_stdr_positive_cutoff <- 0.7

      # Marker designated as completely null
      type_med_negative_cutoff <- 0.2
      type_max_med_negative_cutoff <- 0.7
      type_stdr_negative_cutoff <- 0.5

    } else if (input$marker_diff_eq_option == 'choose_marker_diff_eq') {

      # Marker designated as positive or negative
      type_positive_cutoff <- input$type_positive_cutoff
      type_negative_cutoff <- input$type_negative_cutoff

      # Marker designated as completely positive
      type_med_positive_cutoff <- input$type_med_positive_cutoff
      type_max_med_positive_cutoff <- input$type_max_med_positive_cutoff
      type_stdr_positive_cutoff <- input$type_stdr_positive_cutoff

      # Marker designated as completely null
      type_med_negative_cutoff <- input$type_med_negative_cutoff
      type_max_med_negative_cutoff <- input$type_max_med_negative_cutoff
      type_stdr_negative_cutoff <- input$type_stdr_negative_cutoff
    }

    marker_diff_eq_results_bi <- raw_marker_diff_eq_results()

    # Make sure everything is numeric
    marker_diff_eq_results_bi <- dplyr::mutate_all(marker_diff_eq_results_bi, function(x) as.numeric(as.character(x)))


    kept_markers <- new_names_expression()

    expr <- subset(kept_markers, select=-cluster)
    cell_clustering <- kept_markers$cluster

    # Calculate the median expression
    expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
      dplyr::group_by(cell_clustering) %>% dplyr::summarize_all(dplyr::funs(stats::median))


    expr_median <- subset(expr_median, select=-cell_clustering)

    ## Make everything null ##

    # If the max median per marker per cluster is under a certain value
    x <- names(which(apply(expr_median, 2, max) < type_max_med_negative_cutoff))

    # If the min median per marker per cluster is under a certain value
    y <- names(which(apply(expr_median, 2, min) < type_med_negative_cutoff))

    # If the standard deviation is under the a certain value per marker
    z <- names(which(apply(expr, 2, stats::sd) < type_stdr_negative_cutoff))

    # If both are true
    make_null <- intersect(x, y)
    make_null <- intersect(make_null, z)

    # Make everything null if median and standard deviation are under the threshold
    for (i in make_null) {
      marker_diff_eq_results_bi[i] <- 0
    }

    ## Make everything positive ##

    # If the max median per marker per cluster is over a certain value
    x2 <- names(which(apply(expr_median, 2, max) > type_max_med_positive_cutoff))

    # If the min median per marker per cluster is over a certain value
    y2 <- names(which(apply(expr_median, 2, min) > type_med_positive_cutoff))

    # If the standard deviation is under the a certain value per marker
    z2 <- names(which(apply(expr,2,stats::sd) < type_stdr_positive_cutoff))

    # If all are true
    make_pos <- intersect(x2, y2)
    make_pos <- intersect(make_pos, z2)

    # Make everything positive if median and standard deviation are under the threshold
    for (i in make_pos) {
      marker_diff_eq_results_bi[i] <- 1
    }

    marker_diff_eq_results_bi <- dplyr::mutate_all(marker_diff_eq_results_bi, function(x) as.numeric(as.character(x)))
    
    marker_diff_eq_results_bi[,2:length(marker_diff_eq_results_bi)][marker_diff_eq_results_bi[,2:length(marker_diff_eq_results_bi)] <= as.numeric(type_negative_cutoff)] <- "-"
    
    marker_diff_eq_results_bi[,2:length(marker_diff_eq_results_bi)][marker_diff_eq_results_bi[,2:length(marker_diff_eq_results_bi)] >= as.numeric(type_positive_cutoff)] <- "+"
    
    marker_diff_eq_results_bi[,2:length(marker_diff_eq_results_bi)][marker_diff_eq_results_bi[,2:length(marker_diff_eq_results_bi)] != "+" & marker_diff_eq_results_bi[,2:length(marker_diff_eq_results_bi)] != "-" ] <- ""
    
    marker_diff_eq_results_bi[is.na(marker_diff_eq_results_bi)] <- ""
    
    return(marker_diff_eq_results_bi)
  })

  # Step 2, Median Difference Equation heatmap (cutoffs)
  marker_diff_eq_heatmap_binary_ex <- reactive({
  shiny::req(input$submit_tab1_step2)
   
    # Delineate positives vs negatives
    marker_diff_eq_results_bi <- heatmap_marker_diff_eq_results()

    marker_diff_eq_results_bi[marker_diff_eq_results_bi == "+"] <- 1

    marker_diff_eq_results_bi[marker_diff_eq_results_bi == "-"] <- -1

    marker_diff_eq_results_bi[marker_diff_eq_results_bi == ""] <- 0

    breaks_list <- seq(-1, 1, by = .5)

    # Colors
    if (any(marker_diff_eq_results_bi == -1) & any(marker_diff_eq_results_bi == 1) & any(marker_diff_eq_results_bi == 0)) {
      color_heat <- c("#762A83", "white", "#1B7837")
    } else if (any(marker_diff_eq_results_bi == -1) & any(marker_diff_eq_results_bi == 1) & any(marker_diff_eq_results_bi != 0)) {
      color_heat <- c("#762A83", "#1B7837")
    } else if (any(marker_diff_eq_results_bi != -1) & any(marker_diff_eq_results_bi == 1) & any(marker_diff_eq_results_bi == 0)) {
      color_heat <- c("#1B7837", "white")
    } else if (any(marker_diff_eq_results_bi == -1) & any(marker_diff_eq_results_bi != 1) & any(marker_diff_eq_results_bi == 0)) {
      color_heat <- c("#762A83", "white")
    }

    # Put in numeric cluster order
    marker_diff_eq_results_bi <-
      marker_diff_eq_results_bi[order(marker_diff_eq_results_bi$Cluster), ]

    # Make cluster row names
    marker_diff_eq_results_bi2 <- marker_diff_eq_results_bi[,-1]

    # If there is only 1 marker
    if (ncol(marker_diff_eq_results_bi) == 2) {

      marker_name <- names(marker_diff_eq_results_bi)[2]

      marker_diff_eq_results_bi2 <- as.data.frame(marker_diff_eq_results_bi2)

      names(marker_diff_eq_results_bi2) <- marker_name
    }

    rownames(marker_diff_eq_results_bi2) <- paste("Cluster", marker_diff_eq_results_bi[,1])

    # Make sure everything is numeric
    marker_diff_eq_results_bi2 <- dplyr::mutate_all(marker_diff_eq_results_bi2, function(x) as.numeric(as.character(x)))


    binary_heatmap <- heatmaply::heatmaply(marker_diff_eq_results_bi2, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "#762A83", high = "#1B7837", midpoint = 0, limits = c(-1, 1)),  hide_colorbar = TRUE,
                                           fontsize_col = input$marker_diff_eq_heatmap_raw_column_size, fontsize_row = input$marker_diff_eq_heatmap_raw_row_size, main = "Median Difference Equation, binary results, green=positive & purple=negative", dend = "none", column_text_angle = 70) %>%
      plotly::layout(height = input$heatmap_height_size, titlefont = list(size=20), margin = list(l=50, r=50, b=100, t=100, pad=4))

    return(binary_heatmap)
  })

  # Step 2, Median Difference Equation heatmap (no cutoffs)
  marker_diff_eq_heatmap_raw_ex <- reactive({
    shiny::req(input$submit_tab1_step2)

    # Median Difference Equation: Raw Median Difference Equation scores
    marker_diff_eq_results <- raw_marker_diff_eq_results()

    # Make sure everything is numeric
    marker_diff_eq_results <- dplyr::mutate_all(marker_diff_eq_results, function(x) as.numeric(as.character(x)))

    breaks_list <- seq(min(marker_diff_eq_results), 1, by = .1)
    color_heat <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "PRGn"))(length(breaks_list))

    # Put in numeric cluster order
    marker_diff_eq_results <-
      marker_diff_eq_results[order(marker_diff_eq_results$Cluster), ]

    # Make cluster row names
    marker_diff_eq_results2 <- marker_diff_eq_results[,-1]

    # If there is only 1 marker
    if (ncol(marker_diff_eq_results) == 2) {

      marker_name <- names(marker_diff_eq_results)[2]

      marker_diff_eq_results2 <- as.data.frame(marker_diff_eq_results2)

      names(marker_diff_eq_results2) <- marker_name
    }

    rownames(marker_diff_eq_results2) <- paste("Cluster", marker_diff_eq_results[,1])


    cont_heatmap <- heatmaply::heatmaply(marker_diff_eq_results2,
                                         scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "#762A83", high = "#1B7837", midpoint = 0, limits = c(min(marker_diff_eq_results2), max(marker_diff_eq_results2))),
                                         fontsize_col = input$marker_diff_eq_heatmap_raw_column_size, fontsize_row = input$marker_diff_eq_heatmap_raw_row_size, main = "Median Difference Equation, continuous results",
                                         dend = "none", column_text_angle = 70) %>%
      plotly::layout(height = input$heatmap_height_size, titlefont = list(size=20), margin = list(l=50, r=50, b=100, t=100, pad=4))


    return(cont_heatmap)
  })

  # Step 2, get new names
  new_names_expression <- reactive({
    shiny::req(input$submit_tab1_step2)

    marker_lookup <- df_expr$data

    org_df <- uploaded_file()

    indx <- grepl('cluster', colnames(org_df), ignore.case=TRUE)

    colnames(org_df)[indx] <- tolower(colnames(org_df[indx]))

    # Change marker names if applicable
    kept_markers <- org_df[, c(marker_lookup$`Original marker names`, "cluster")]

    change_names <- base::match(names(kept_markers), marker_lookup$`Original marker names`)

    names(kept_markers)[stats::na.omit(change_names)] = marker_lookup$`New marker names`[!is.na(change_names)]


    return(kept_markers)
  })

  # Step 2, median expression heatmap
  expr_heatmap_org_ex <- reactive({
    shiny::req(input$submit_tab1_step2)

    kept_markers <- new_names_expression()

    expr <- subset(kept_markers, select=-cluster)
    cell_clustering <- kept_markers$cluster

    # Calculate the median expression
    expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
      dplyr::group_by(cell_clustering) %>% dplyr::summarize_all(dplyr::funs(stats::median))

    # Put in numeric cluster order
    expr_median <-
      expr_median[order(expr_median$cell_clustering), ]

    # Make cluster row names
    expr_median <- expr_median[,-1]
    rownames(expr_median) <- paste("Cluster", rownames(expr_median))

    # Colors for the heatmap
    color_heat <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "PRGn"))(100)

    expr_heatmap <- heatmaply::heatmaply(expr_median, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "#762A83", high = "#1B7837", midpoint = 0, limits = c(min(expr_median), max(expr_median))),
                                         fontsize_col = input$marker_diff_eq_heatmap_raw_column_size, fontsize_row = input$marker_diff_eq_heatmap_raw_row_size, main = "Median expression values", dend = "none", column_text_angle = 70) %>%
      plotly::layout(height = input$heatmap_height_size, titlefont = list(size=20), margin = list(l=50, r=50, b=100, t=100, pad=4))


    return(expr_heatmap)
  })

  # Step 2, density
  expr_density <- reactive({
    shiny::req(input$submit_tab1_step2)

    kept_markers <- new_names_expression()

    # Downsample number of cells for density plots so it is not time prohibitive
    if (nrow(kept_markers) > 50000) {
      kept_markers <- kept_markers %>% dplyr::slice_sample(n = 50000, replace = FALSE)
    }

    expr_1 <- subset(kept_markers, select=-cluster)

    cell_clustering <- kept_markers$cluster
    cell_clustering <- factor(cell_clustering)

    # Calculate the median expression
    expr_median <- data.frame(expr_1, cell_clustering = cell_clustering) %>%
      dplyr::group_by(cell_clustering) %>% dplyr::summarize_all(dplyr::funs(stats::median))

    ### Data organized per cluster
    ggd <- reshape2::melt(data.frame(Cluster = cell_clustering, expr_1),
                          id.vars = "Cluster", value.name = "Expression",
                          variable.name = "antigen")

    ggd$`All cells` <- "no"

    ### The reference data
    ggd_bg <- ggd
    ggd_bg$Cluster <- "All cells"
    ggd_bg$`All cells` <- "yes"

    ggd_plot <- rbind(ggd, ggd_bg)
    ggd_plot <- ggd_plot[!is.na(ggd_plot$antigen),]
    ggd_plot <- ggd_plot[!is.na(ggd_plot$Cluster),]
    ggd_plot$Cluster <- as.factor(ggd_plot$Cluster)
    ggd_plot$Cluster <- factor(paste("Cluster", ggd_plot$Cluster),
                               levels = c(rev(paste("Cluster",levels(cell_clustering))), "Cluster All cells"))

    levels(ggd_plot$Cluster) <- gsub("Cluster All cells", "All cells", levels(ggd_plot$Cluster), fixed=TRUE)

    density_plots <- ggplot2::ggplot() +
      ggridges::geom_density_ridges(data = ggd_plot, ggplot2::aes(x = Expression, y = Cluster,
                                                                  color = `All cells`, fill = `All cells`), alpha = 0.3) +
      ggplot2::scale_fill_manual(values = c("#762A83", "#1B7837")) +
      ggplot2::scale_color_manual(values = c("#762A83", "#1B7837"))

    return(density_plots)
  })

  # Step 2, density plot
  expr_density_ex <- reactive({
    shiny::req(input$submit_tab1_step2)

    density_plots <- expr_density()

    # density_plots <- ggplot2::ggplot() +
    #   ggridges::geom_density_ridges(data = ggd_plot, ggplot2::aes(x = Expression, y = Cluster,
    #                                            color = `All cells`, fill = `All cells`), alpha = 0.3) +
    #   ggplot2::scale_fill_manual(values = c("#762A83", "#1B7837")) +
    #   ggplot2::scale_color_manual(values = c("#762A83", "#1B7837"))
    #
    #
    density_plots <- density_plots +
      ggplot2::scale_x_continuous(breaks = integer_breaks(n = 4)) +
      ggplot2::facet_wrap( ~ antigen, scales = "free_x", nrow = input$density_facet) +
      ggridges::theme_ridges() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$density_row_size), #input$marker_diff_eq_heatmap_raw_row_size),  #17
                     axis.text.y = ggplot2::element_text(size = 17), #input$marker_diff_eq_heatmap_raw_col_size), #17
                     axis.title.x = ggplot2::element_text(hjust = .5, size = 25),
                     axis.title.y = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(size = 17),
                     legend.position = "none")


    return(density_plots)

  })

  # Step 2, max expression value per marker
  expr_max <- reactive({
    shiny::req(input$submit_tab1_step2)

    kept_markers <- new_names_expression()


    expr <- subset(kept_markers, select=-cluster)

    max_med <- apply(expr, 2, max, na.rm=TRUE)


    max_med <- t(as.data.frame(max_med))


    return(max_med)
  })

  # Step 2, max expression value per marker, heatmap
  expr_heatmap_max_ex <- reactive({
    shiny::req(input$submit_tab1_step2)

    max_med <- expr_max()

    color_heat <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "PRGn"))(100)


    expr_heatmap <- heatmaply::heatmaply(max_med, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "#762A83", high = "#1B7837", midpoint = 0, limits = c(min(max_med), max(max_med))),
                                         fontsize_col = input$marker_diff_eq_heatmap_raw_column_size, fontsize_row = input$marker_diff_eq_heatmap_raw_row_size, main = "Maximum median expression value", dend = "none", showticklabels = c(TRUE, FALSE), column_text_angle = 70) %>%
      plotly::layout(height = 425, titlefont = list(size=20), margin = list(l=50, r=50, b=100, t=100, pad=4))

    return(expr_heatmap)
  })

  # Step 2, min expression value per marker
  expr_min <-  reactive({
    shiny::req(input$submit_tab1_step2)

    kept_markers <- new_names_expression()

    expr <- subset(kept_markers, select=-cluster)

    min_med <- apply(expr, 2, min, na.rm=TRUE)

    min_med <- t(as.data.frame(min_med))

    return(min_med)
  })

  # Step 2, min expression value per marker, heatmap
  expr_heatmap_min_ex <- reactive({
    shiny::req(input$submit_tab1_step2)

    min_med <- expr_min()

    color_heat <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "PRGn"))(100)


    expr_heatmap <- heatmaply::heatmaply(min_med, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "#762A83", high = "#1B7837", midpoint = 0, limits = c(min(min_med), max(min_med))),
                                         fontsize_col = input$marker_diff_eq_heatmap_raw_column_size, fontsize_row = input$marker_diff_eq_heatmap_raw_row_size, main = "Minimum median expression value", dend = "none", showticklabels = c(TRUE, FALSE), column_text_angle = 70) %>%
      plotly::layout(height = 425, titlefont = list(size=20), margin = list(l=50, r=50, b=100, t=100, pad=4))

    return(expr_heatmap)

  })

  # Step 2, standard deviation expression value per marker
  expr_std <- reactive({
    shiny::req(input$submit_tab1_step2)

    kept_markers <- new_names_expression()

    expr <- subset(kept_markers, select=-cluster)

    sd_med <- apply(expr, 2, stats::sd, na.rm=TRUE)

    sd_med <- t(as.data.frame(sd_med))

    return(sd_med)
  })

  # Step 2, standard deviation expression value per marker, heatmap
  expr_heatmap_std_ex <- reactive({

    shiny::req(input$submit_tab1_step2)

    sd_med <- expr_std()

    color_heat <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "PRGn"))(100)


    expr_heatmap <- heatmaply::heatmaply(sd_med, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "#762A83", high = "#1B7837", midpoint = 0, limits = c(min(sd_med), max(sd_med))),
                                         fontsize_col = input$marker_diff_eq_heatmap_raw_column_size, fontsize_row = input$marker_diff_eq_heatmap_raw_row_size, main = "Standard Deviation of expression values", dend = "none", showticklabels = c(TRUE, FALSE), column_text_angle = 70) %>%
      plotly::layout(height = 425, titlefont = list(size=20), margin = list(l=50, r=50, b=100, t=100, pad=4))

    return(expr_heatmap)

  })

  output$marker_diff_eq_heatmap_binary_ex_plot2 = shiny::renderUI({
     #  shiny::req(input$submit_tab1_step2)
   
    heatmap_height = input$heatmap_height_size + 20

    density_height = input$density_height_size + 20

    shiny::tagList(
      plotly::plotlyOutput("marker_diff_eq_heatmap_binary_ex_plot", width = "100%", height = paste0(heatmap_height)),
      plotly::plotlyOutput("marker_diff_eq_heatmap_raw_ex_plot", width = "100%", height = paste0(heatmap_height)),
      plotly::plotlyOutput("expr_heatmap_org_ex_plot", width = "100%", height = paste0(heatmap_height)),
      plotOutput("expr_density_ex_plot", width = "100%", height = paste0(density_height)),
      plotly::plotlyOutput("expr_heatmap_max_ex_plot"),
      plotly::plotlyOutput("expr_heatmap_min_ex_plot"),
      plotly::plotlyOutput("expr_heatmap_std_ex_plot"))
  })

  # Step 2, when help button is clicked show message
  shiny::observeEvent(input$help_tab1_step2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 2<br>Adjust Median Difference Equation Parameters"),
      HTML("This step shows the results of the Median Difference Equation. Users can modify equation parameters via the sidebar, with changes reflected immediately in the results. After this step, parameters cannot be further adjusted, but individual overrides can still be made in Step 3.
      <br>
      <br>
      The median difference results are displayed as a data frame, where the rows correspond to the input clusters and the columns represent the markers. The values are color-coded: positives (green), negatives (purple), and nulls (white). A heatmap displays the continuous results, showing scaled median differences before categorical cutoffs are applied.
      <br>
      <br>
      Additional plots depict the underlying expression data. A heatmap shows the median expression values and a density plot shows the distribution of values for each cluster, as well as for all clusters combined. Other figures show the maximum and minimum median expression across clusters for each marker. This is included because these values are used when designating markers as entirely positive or null across all clusters. Likewise, the final figure shows standard deviation values for each marker, which is similarly used when designated markers as completely positive or null.
      <br>
      <br>
      The dimensions and label sizes of the plots can be adjusted via slider bars at the bottom. Heatmaps are interactive ??? hovering over any section reveals the row, column, and value, and users can zoom in on specific areas. These plots can also be downloaded as PNG files using the camera icon. Density plots are not interactive, but users can right-click to copy or save them.
      <br>
      <br>"),
      tags$div("The Median Difference Equation and additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/3.1.-Input-Expression-Data-and-Default-References', "here.", target="_blank")),
      easyClose = TRUE))
  })


  ######~#~# Step 3 - Override Median Difference Equation results if desired #~#~######

  # Step 3, side panel

  # Reset everything
  shiny::observeEvent(input$reset_tab1_step3, {
    session$reload()
  })

  # Download example marker-cell type reference data
  output$download_example_marker_cell_type_1 <- shiny::downloadHandler(
    # File name
    filename = function() {
      c("Lee_AML_cell_types_markers.csv")
    },
    # Write to csv
    content = function(file) {
      utils::write.csv(Lee_AML_cell_types_markers,  file, row.names = FALSE)
    }
  )

  # Step 3, help buttons for side panel parameters
  shiny::observeEvent(input$help_ref_type,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select Cell Type Reference"),
      HTML("The user can choose between two options for the descriptive cell type reference:
      <br>
      <br>
       <ul>
      <li>Default Reference: This option does not require the user to upload a reference file. It uses various built-in references, including the Cell Ontology, Protein Ontology, Wikidata, and curated datasets included with the application.</li>
      <br>
      <li>File Upload: This option allows the user to upload a custom reference file containing descriptive cell type names and their associated marker definitions. The file should be a CSV with two columns:</li>
      <br>
      <ul>
      <li>Name: The cell type identifier</li>
      <br>
      <li>Markers: The protein marker names with associated qualifiers ('positive', 'negative', 'low', 'high', '+', '++', '-'). Markers and qualifiers should be separated by commas (,).</li>
      <br>
      </ul>
      </ul>"),
      tags$div("An example marker-cell type reference file is available within the application (click 'Download Example' above the file upload box) and on the ", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/2.-Example-Data', "GitHub.", target="_blank")),
      easyClose = TRUE))
  })

  # Help box for species type
  shiny::observeEvent(input$help_species,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Specify Species"),
      HTML("Some proteins and cell type names are species-specific within ontologies. Therefore, when using 'Default References', specify the species of your sample(s)."),
      easyClose = TRUE))
  })

  # Help box for ontology types
  shiny::observeEvent(input$help_ontology_type_1,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select Cell Type Ontology(s)"),
      HTML("When using Default References, specify the cell type ontology. The primary source should be the 'Cell Ontology', which contains hundreds of cell types with marker definitions.
      <br>
      <br>
      The 'Provisional Cell Ontology' can also be optionally included, containing cell types that are currently only provisionally defined but may still be useful.
      <br>
      <br>
      At least one cell type ontology must be included."),
      easyClose = TRUE))
  })

  # Step 3, main panel

  # Step 3, run Median Difference Equation
  shiny::observeEvent(input$submit_tab1_step2_5, {

    df_marker_diff_eq$data <- heatmap_marker_diff_eq_results()

    shiny.destroy::removeOutput("help_tab1_step2")
    shiny.destroy::removeOutput("overview_tab1_step2")
    shinyjs::show("overview_tab1_step3")
    shinyjs::show("marker_diff_eq_results_editable")
    shiny.destroy::removeOutput("marker_diff_eq_heatmap_binary_ex_plot")
    shiny.destroy::removeOutput("all_sidebar_tab1_step2")
    shinyjs::show("all_sidebar_tab1_step3")
    shiny.destroy::removeOutput("expr_heatmap_max_ex_plot")
    shiny.destroy::removeOutput("expr_heatmap_min_ex_plot")
    shiny.destroy::removeOutput("expr_heatmap_std_ex_plot")

    # Median Difference Equation results
    output$marker_diff_eq_results_editable <- DT::renderDT({
      DT::datatable(
        df_marker_diff_eq$data,
        selection = 'none',
        editable = list(target = "cell", disable = list(columns = 0)), #TRUE,
        rownames = FALSE,
        escape = FALSE,
        options = list(pageLength = 75)
      )  %>%
        DT::formatStyle(names(df_marker_diff_eq$data), color = DT::styleEqual(c("-", "+"), c("white", "white")), fontWeight = 'bold', fontSize = '150%', textAlign = 'center', backgroundColor = DT::styleEqual(c("-", "+"), c("#762A83", "#1B7837")))
    }, server = FALSE)

    if (!is.null(heatmap_marker_diff_eq_results())) {
      shinyjs::show("download_marker_diff_eq_csv")
      shinyjs::show("help_tab1_step3")
    }

  })

  # Step 3, create editable dataframes
  df_marker_diff_eq <- shiny::reactiveValues(data = NULL)

  # Step 3, get max, min, std for table download
  max_min_std <- reactive({

    shiny::req(input$submit_tab1_step2_5)
    #  shiny::req(input$download_marker_diff_eq_csv)

    max <- expr_max()

    min <- expr_min()

    std <- expr_std()

    all <- rbind(max, min, std)

    Measurement <- c("Max. medium", "Min. medium", "Standard deviation")

    all <- cbind(Measurement, all)

    return(all)
  })

  # Step 3, parameter info for table download
  marker_diff_eq_parameter_df <- reactive({

    shiny::req(input$submit_tab1_step2_5)
    #  shiny::req(input$download_marker_diff_eq_csv)

    if (input$transform_option == 'transform_no') {
      transform_choice <- "False"
      cofactor_num <- NA
    } else if (input$transform_option == 'transform_yes') {
      transform_choice <- input$transform_option
      cofactor_num <- input$type_cofactor
    }

    # if (is.null(input$seed_1)) {
    #    seed_num <- "False"

    if (input$random_1 != 1) {
      downsample_num <- "False"
    } else {
      downsample_num <- input$num_random_1
    }

    if (input$seed_1 != 1) {
      seed_num <- "False"
    } else {
      seed_num <- input$num_seed_1
    }

    mg_string <- input$type_MG_markers
    mg_string <- gsub(" ", "", mg_string, fixed = TRUE)
    if (mg_string == "") {
      mg_markers <- "FALSE"
    } else {
      mg_markers <- mg_string
    }

    if (input$marker_diff_eq_option == 'default_marker_diff_eq') {

      # Marker designated as positive or negative
      type_positive_cutoff <- 0.25
      type_negative_cutoff <- -0.75

      # Marker designated as completely positive
      type_med_positive_cutoff <- 0.2
      type_max_med_positive_cutoff <- 1
      type_stdr_positive_cutoff <- 0.7

      # Marker designated as completely null
      type_med_negative_cutoff <- 0.2
      type_max_med_negative_cutoff <- 0.7
      type_stdr_negative_cutoff <- 0.5

    } else if (input$marker_diff_eq_option == 'choose_marker_diff_eq') {

      # Marker designated as positive or negative
      type_positive_cutoff <- input$type_positive_cutoff
      type_negative_cutoff <- input$type_negative_cutoff

      # Marker designated as completely positive
      type_med_positive_cutoff <- input$type_med_positive_cutoff
      type_max_med_positive_cutoff <- input$type_max_med_positive_cutoff
      type_stdr_positive_cutoff <- input$type_stdr_positive_cutoff

      # Marker designated as completely null
      type_med_negative_cutoff <- input$type_med_negative_cutoff
      type_max_med_negative_cutoff <- input$type_max_med_negative_cutoff
      type_stdr_negative_cutoff <- input$type_stdr_negative_cutoff
    }

    df <- data.frame("Parameter" = c("Date created", "Manual gating", "Transformation","Cofactor", "Downsample", "Seed", "Positive cutoff","Negative cutoff","Minimum medium all positive cutoff","Maximum medium all positive cutoff","Standard deviation all positive cutoff", "Minimum medium all null cutoff", "Maximum medium all null cutoff", "Standard deviation all null cutoff"),
                     "Value" = c(as.character(Sys.time()), mg_markers, transform_choice, cofactor_num, downsample_num, seed_num, type_positive_cutoff, type_negative_cutoff, type_med_positive_cutoff, type_max_med_positive_cutoff, type_stdr_positive_cutoff, type_med_negative_cutoff, type_max_med_negative_cutoff, type_stdr_negative_cutoff))
    return(df)
  })

  # Step 3, get marker difference equation results in string format for table download
  merged_marker_diff_eq_results <- reactive({

    shiny::req(input$submit_tab1_step2_5)

    marker_diff_eq_results_string <- df_marker_diff_eq$data

    marker_diff_eq_results_string[,2:length(marker_diff_eq_results_string)][marker_diff_eq_results_string[,2:length(marker_diff_eq_results_string)]==""]<-"xx"

    if (ncol(marker_diff_eq_results_string) == 2) {
      marker_name <- colnames(marker_diff_eq_results_string)[2]
      marker_diff_eq_results_string[,2] <- paste0(marker_name, marker_diff_eq_results_string[,2])
    } else {
      for (i in colnames(marker_diff_eq_results_string[,2:length(marker_diff_eq_results_string)])){
        marker_diff_eq_results_string[,i] <- paste0(i, marker_diff_eq_results_string[,i])
      }
    }

    for (i in 2:length(marker_diff_eq_results_string)) {
      marker_diff_eq_results_string[,i] <- gsub(".*xx", "", marker_diff_eq_results_string[,i])
    }

    marker_diff_eq_results_string <- marker_diff_eq_results_string %>%
      tidyr::unite("Marker definition", 2:length(marker_diff_eq_results_string), sep=" ")

    for (i in 2:length(marker_diff_eq_results_string)) {
      marker_diff_eq_results_string[,i] <- gsub("\\s+", " ", marker_diff_eq_results_string[,i])
      marker_diff_eq_results_string$`Marker definition` <- trimws(marker_diff_eq_results_string$`Marker definition`)
      marker_diff_eq_results_string[,i] <- gsub("\\s+", ", ", marker_diff_eq_results_string[,i])
    }

    return(marker_diff_eq_results_string)
  })

  # Step 3, when help button is clicked show message
  shiny::observeEvent(input$help_tab1_step3,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 3<br>Override Median Difference Equation Results"),
      HTML("In this step, users can override the results from the Median Difference Equation. The categorical results from Step 2 are shown in an editable data frame. Users can modify the designations for each marker-cluster combination by clicking the box and typing either positive ('+'), negative ('-'), or null (' ').
      <br>
      <br>
      A 'Download' button at the bottom of the page allows users to export the results as a ZIP file containing multiple CSVs. These files reflect the current editable data frame, so any changes will be immediately reflected in the CSV exports.
      <br>
      <br>"),
      tags$div("The Median Difference Equation and additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/3.1.-Input-Expression-Data-and-Default-References', "here.", target="_blank")),
      easyClose = TRUE))
  })

  # Step 3, edit table
  shiny::observeEvent(input$marker_diff_eq_results_editable_cell_edit, {
    shiny::req(input$submit_tab1_step2_5)
    df_marker_diff_eq$data <<-
      DT::editData(df_marker_diff_eq$data,
                   input$marker_diff_eq_results_editable_cell_edit,
                   'marker_diff_eq_results',
                   rownames = FALSE)
  })


  # Download marker definition csv files
  output$download_marker_diff_eq_csv <- shiny::downloadHandler(
    # Filename
    filename = function() {
      paste0(gsub(".csv|.xls|.xlsx", "", input$upload_expression_csv),
             "-marker-definitions-", Sys.Date(), ".zip")
    },

    content = function(file) {
      # definition of content to download
      to_dl <- list(
        # names to use in file names
        names = list(a = "Full_marker_definitions",
                     b = "Raw_marker_diff_eq_scores",
                     c = "Positive_negative_marker_diff_eq_scores",
                     d = "Max_min_std_per_marker",
                     e = "Set_parameters"),
        # data
        data = list(a = merged_marker_diff_eq_results(),
                    b = raw_marker_diff_eq_results(),
                    c = df_marker_diff_eq$data,
                    d = max_min_std(),
                    e = marker_diff_eq_parameter_df())
      )

      # temp dir for the csv's as we can only create
      # an archive from existent files and not data from R
      twd <- setwd(tempdir())
      on.exit(setwd(twd))
      files <- NULL

      # Loop on data to download and write individual csv's
      for (i in c("a", "b", "c", "d", "e")) {
        fileName <- paste0(to_dl[["names"]][[i]], ".csv") # csv file name
        utils::write.csv(to_dl[["data"]][[i]], fileName, row.names = FALSE) # write csv in temp dir

        files <- c(files, fileName) # store written file name
      }
      # Create zip file with all the csv files
      zip::zip(file, files)
    }
  )
  ######~#~# Step 4 - Get marker synonyms and matched PRO/GO terms, user can input new names if left unmatched  #~#~######

  # Step 4, side panel

  # Reset everything
  shiny::observeEvent(input$reset_tab1_step4, {
    session$reload()
  })


  #~#~#~#~#~# Uploaded reference & Internal references  #~#~#~#~#~#

  # Step 4, main panel

  # Step 4, change format of marker difference equation results
  reformatted_data_1 <- shiny::eventReactive(input$submit_tab1_step3, {

    # Marker difference equation results
    inputted_df <- merged_marker_diff_eq_results()

    # Change column names
    colnames(inputted_df) <- c("Cluster", "Markers")

    # Remove clusters that had no marker designations
    inputted_df <- inputted_df[!(is.na(inputted_df$Markers) | inputted_df$Markers==""), ]

    # Split marker column by , or ;
    new_df <-
      data.frame(Cluster = inputted_df$Cluster, do.call('rbind', strsplit(
        as.character(inputted_df$Markers), split = "[,;]+"
      )))

    # Remove spaces
    new_df <- as.data.frame(lapply(new_df, function(y) gsub(' ', '', y)))

    # Remove duplicates per row
    new_df <-
      t(apply(new_df, 1, function(x)
        replace(x, duplicated(x), "")))

    # Number of cell type names
    num_names <- nrow(new_df)

    # Change dataframe format
    new_df <- reshape2::melt(new_df, id = "Cluster")

    # Match numbers to names
    melt_key <- new_df[(1:num_names), , drop = FALSE]

    # Continue getting correct format
    new_df <- new_df[-(1:num_names), , drop = FALSE]
    new_df <- new_df[, -2]
    new_df <- new_df[!(new_df$value == ""),]

    # Determine if positive or negative
    new_df$Sign[grepl("positive|\\+$", new_df$value, ignore.case = TRUE)] <-
      "Positive"
    new_df$Sign[grepl("negative|\\-$", new_df$value, ignore.case = TRUE)] <-
      "Negative"
    new_df$Sign[grepl("high|\\+\\+|\\+\\+\\+$", new_df$value, ignore.case =
                        TRUE)] <- "High"
    new_df$Sign[grepl("low$", new_df$value, ignore.case = TRUE)] <-
      "Low"

    # Remove the sign from the marker
    new_df$value <-
      gsub("\\+$", "", gsub("\\++$", "",  gsub("\\-$", "", gsub(
        "\\--$", "", gsub(
          "low$",
          "",
          ignore.case = TRUE,
          gsub(
            "high$",
            "",
            ignore.case = TRUE,
            gsub(
              "positive$",
              "",
              ignore.case = TRUE,
              gsub(
                "negative$",
                "",
                as.character(new_df$value),
                ignore.case = TRUE
              )
            )
          )
        )
      ))))

    # Add names
    new_df <- merge(new_df, melt_key, by = "Var1", all = TRUE)

    # Remove old, unneeded column
    new_df <- subset(new_df, select = -c(Var1, Var2))

    # New column names
    colnames(new_df) <- c("Marker", "Sign", "Cluster")

    MG_df <- added_MG_markers()

    if (is.data.frame(MG_df) == TRUE) {
      new_df <- rbind(new_df, MG_df)
    }

    # Make everything capital
    new_df$Marker <- toupper(new_df$Marker)

    if (input$reference_type_1 == 'internal_ref_1') {
      # Remove CD45
      new_df <- new_df[new_df$Marker != "CD45", ]
    }


    # Order and sort by cluster
    new_df <- new_df[order(new_df$Cluster), ]

    if (sum(is.na(new_df$Sign)) > 0) {
      if (sum(is.na(new_df$Sign)) > 1) {
        no_sign <- unique(new_df$Marker[is.na(new_df$Sign)])

        return(no_sign)

      } else if (sum(is.na(new_df$Sign)) == 1) {
        no_sign <- new_df$Marker[is.na(new_df$Sign)]

        return(no_sign)
      }
    } else {
      # Make everything uppercase for the filter
      new_df$Marker <- toupper(new_df$Marker)

      return(new_df)
    }
  })

  # Step 4, outputs error if there is a marker designated twice or more in a cluster
  marker_sign_error_1 <- shiny::eventReactive(input$submit_tab1_step3, {

    df <- reformatted_data_1()

    same_marker_cluster <- data.frame(table(df$Cluster, df$Marker))

    same_marker_cluster <- same_marker_cluster[same_marker_cluster$Freq > 1,]

    if (nrow(same_marker_cluster) >= 1) {
      P <- "1"
      return(P)
    }
  })

  output$ui_marker_sign_error_1 <- shiny::renderUI({
    shiny::req(input$submit_tab1_step3)
    if (length(marker_sign_error_1()) >= 1) {
      tags$div(id = "id_marker_sign_error_1", h4(
        "Error: A marker is repeated at least twice within 1 or more clusters. Please ensure there are no repeated markers within each cluster (including specified manually gated markers)."
      ))}
  })

  #~#~#~#~#~#  Uploaded reference #~#~#~#~#~#

  # Step 4, uploaded reference, when help button is clicked show message
  shiny::observeEvent(input$help_uploaded_ref_1,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 4<br>Match Marker Patterns to Cell Types"),
      HTML("This step displays the results of the cell type matching process. The data frame shows:
      <br>
      <br>
      <li>Input clusters</li>
      <br>
      <li>Matched reference cell type names</li>
      <br>
      <li>Markers that are matched between the input clusters and the reference cell types</li>
      <br>
      <li>Markers present in the input clusters</li>
      <br>
      <li>Markers present in the matched reference cell types</li>
      <br>
      <li>Markers that are directly contradicted between the input clusters and the reference cell types</li>
      <br>
      <br>
      Three additional columns provide statistical metrics for evaluating the match:
      <br>
      <br>
      <li>Percentages of matched markers over the total markers in the reference cell type</li>
      <br>
      <li>Percentages of matched markers over the total markers in the input</li>
      <br>
      <li>The matching scores derived from Equation 3</li>
      <br>
      The 'Download' button at the bottom allows users to export a CSV file containing all relevant data, including additional columns detailing the number of markers in the input, reference cell types, matches, and contradictions.
      <br>
      <br>
     "),
      tags$div("Equation 3 and additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/3.2.-Input-Expression-Data-and-Uploaded-References', "here.", target="_blank")),
      easyClose = TRUE))
  })

  # Step 4, uploaded reference, check if file type of uploaded reference file is acceptable
  check_reference_file_1 <- shiny::eventReactive(input$submit_tab1_step3, {
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)

    # File inputted
    file_input <- input$upload_ref_csv_1

    if (endsWith(file_input$datapath, ".csv") | endsWith(file_input$datapath, ".xls") | endsWith(file_input$datapath, ".xlsx")) {
      return(NULL)
    }  else {
      return("a")
    }
  })

  # Step 4, uploaded reference, outputs error stating that the uploaded reference file is in the wrong format
  output$ui_error_check_reference_file_1 <- shiny::renderUI({
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)
  #  shiny::req(input$submit_tab1_step3)

    if (is.null(check_reference_file_1()) == FALSE) {
      tags$div(id = "invalid_ref_file_1", h4(
        "Error: Invalid uploaded reference format. Please upload a valid comma separated values (CSV) or Excel (XLSX or XLS) file."
      ))
    }})

  # Step 4, uploaded reference, change format of uploaded reference file
  reformatted_ref_1 <- shiny::eventReactive(input$submit_tab1_step3, {
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)

    if (is.null(check_reference_file_1()) == TRUE) {

      # File inputted
      file_input <- input$upload_ref_csv_1

      if (endsWith(file_input$datapath, ".csv")) {

        inputted_df <- utils::read.csv(file_input$datapath, header = TRUE)

      } else if (endsWith(file_input$datapath, ".xls")) {

        inputted_df <- readxl::read_excel(file_input$datapath, col_names = TRUE)

      }  else if (endsWith(file_input$datapath, ".xlsx")) {

        inputted_df <- readxl::read_excel(file_input$datapath, col_names = TRUE)
      }


      # Normalize all dashes
      inputted_df <- as.data.frame(lapply(inputted_df, function(y) gsub("\\p{Pd}", "-", y, perl=TRUE)))

      # Ensure first letter is capitalized and the rest in lower case
      colnames(inputted_df) <- tolower(colnames(inputted_df))
      substr(colnames(inputted_df), 1, 1) <- toupper(substr(colnames(inputted_df), 1, 1))

      # Check that column names are correct. If not, return a warning and stop.
      if (!all(c("Name", "Markers") == colnames(inputted_df))) {
        return(NULL)
      } else {
        # Split marker column by , or ;
        new_df <-
          data.frame(Name = inputted_df$Name, do.call('rbind', strsplit(
            as.character(inputted_df$Markers), split = "[,;]+"
          )))

        # Remove spaces
      #  new_df <- as.data.frame(lapply(new_df, function(y) gsub(' ', '', y)))
        new_df[2:ncol(new_df)] <- as.data.frame(lapply(new_df[2:ncol(new_df)], function(y) gsub(' ', '', y)))

        # Remove duplicates per row
        new_df <-
          t(apply(new_df, 1, function(x)
            replace(x, duplicated(x), "")))

        # Number of cell type names
        num_names <- nrow(new_df)

        # Change dataframe format
        new_df <- reshape2::melt(new_df, id = "Name")

        # Match numbers to names
        melt_key <- new_df[(1:num_names), , drop = FALSE]

        # Continue getting correct format
        new_df <- new_df[-(1:num_names), , drop = FALSE]
        new_df <- new_df[, -2]
        new_df <- new_df[!(new_df$value == ""),]

        # Determine if positive or negative
        new_df$Sign[grepl("positive|\\+$", new_df$value, ignore.case = TRUE)] <-
          "Positive"
        new_df$Sign[grepl("negative|\\-$", new_df$value, ignore.case = TRUE)] <-
          "Negative"
        new_df$Sign[grepl("high|\\+\\+|\\+\\+\\+$", new_df$value, ignore.case =
                            TRUE)] <- "High"
        new_df$Sign[grepl("low$", new_df$value, ignore.case = TRUE)] <-
          "Low"

        # Remove the sign from the marker
        new_df$value <-
          gsub("\\+$", "", gsub("\\++$", "",  gsub("\\-$", "", gsub(
            "\\--$", "", gsub(
              "low$",
              "",
              ignore.case = TRUE,
              gsub(
                "high$",
                "",
                ignore.case = TRUE,
                gsub(
                  "positive$",
                  "",
                  ignore.case = TRUE,
                  gsub(
                    "negative$",
                    "",
                    as.character(new_df$value),
                    ignore.case = TRUE
                  )
                )
              )
            )
          ))))

        # Add names
        new_df <- merge(new_df, melt_key, by = "Var1", all = TRUE)

        # Remove old, unneeded column
        new_df <- subset(new_df, select = -c(Var1, Var2))

        # New column names
        colnames(new_df) <- c("Marker", "Sign", "Name")

        if (sum(is.na(new_df$Sign)) > 0) {
          if (sum(is.na(new_df$Sign)) > 1) {
            no_sign <- unique(new_df$Marker[is.na(new_df$Sign)])

            return(no_sign)

          } else if (sum(is.na(new_df$Sign)) == 1) {
            no_sign <- new_df$Marker[is.na(new_df$Sign)]

            return(no_sign)
          }
        } else {

          # Make everything uppercase for the filter
          new_df$Marker <- toupper(new_df$Marker)

          return(new_df)
        }
      }
    }
  })

  # Step 4, uploaded reference, check if inputted marker is in uploaded reference
  marker_in_uploaded_reference_1 <- reactive({
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)

    input_df <- reformatted_data_1()
    ref_df <- reformatted_ref_1()

    different_markers <- setdiff(input_df$Marker, ref_df$Marker)

    if (length(different_markers) != 0) {
      return(different_markers)
    } else {
      return(NULL)
    }
  })

  # Step 4, uploaded reference, outputs warning if there is no sign for at least 1 marker in the uploaded reference
  output$ui_warning_no_sign_ref_1 <- shiny::renderUI({
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)

    if (is.data.frame(reformatted_ref_1()) == FALSE) {
      if (length(reformatted_ref_1()) == 1) {
        tags$div(id = "invalid_sign_1", h4(
          paste(
            "Error: Invalid uploaded reference format. The marker",
            reformatted_ref_1(),
            "does not have a valid sign."
          )
        ))
      } else if (length(reformatted_ref_1()) > 1) {
        no_sign <- reformatted_ref_1()
        no_sign <- data.frame(x = no_sign)
        no_sign <- paste(no_sign$x, collapse = ", ")
        tags$div(id = "invalid_sign_1", h4(
          paste(
            "Error: Invalid uploaded reference format. The markers",
            no_sign,
            "do not have valid signs."
          )
        ))
      }
    }
  })

  # Step 4, uploaded reference, error message if uploaded reference has wrong column names
  output$ui_error_column_names_ref_1 <- shiny::renderUI({
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)

    if (is.null(check_reference_file_1()) == TRUE) {
      if (is.null(reformatted_ref_1()) == TRUE) {
        h4(
          "Error: Invalid uploaded reference format. There should be 2 columns, the first labeled 'Name' and the second labeled 'Markers'."
        )
      }
    }
  })

  # Step 4, uploaded reference, get statement indicating if/which markers are not in uploaded reference
  output$ui_warning_unmatched_markers_uploaded_ref_1 <- shiny::renderUI({
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)

    if (is.data.frame(reformatted_data_1()) == TRUE &
        is.data.frame(reformatted_ref_1()) == TRUE) {
      if (length(marker_sign_error_1()) == 0) {

        new_df <- reformatted_data_1()
        ref_df <- reformatted_ref_1()

        different_markers <- setdiff(new_df$Marker, ref_df$Marker)

        if (length(different_markers) < length(unique(new_df$Marker))) {
          if (length(marker_in_uploaded_reference_1()) > 1) {
            marker_string <-
              paste(marker_in_uploaded_reference_1(), collapse = ", ")
            h4(
              paste(
                "Warning: The inputted markers",
                marker_string,
                "are not in the uploaded reference."
              )
            )
          } else if (length(marker_in_uploaded_reference_1()) == 1) {
            h4(
              paste(
                "Warning: The inputted marker",
                marker_in_uploaded_reference_1(),
                "is not in the uploaded reference."
              )
            )
          }
        } else {
          h4("Error: None of the inputted markers are in the uploaded reference.")
        }
      }
    }
  })

  # Step 4, uploaded reference, get results on how input and uploaded reference match
  uploaded_ref_matches_1 <- reactive({
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)
   
    # Make sure the uploaded reference and input are dataframes
    if (is.data.frame(reformatted_data_1()) == TRUE &
        is.data.frame(reformatted_ref_1()) == TRUE) {
      new_df <- reformatted_data_1()
      ref_df <- reformatted_ref_1()

      # Look at each cluster individually
      each_cluster <- split(new_df, new_df$Cluster)

      # Make empty, will be added to (+1) if there is a match
      ref_df[, 'Match'] <- 0

      # Empty list to put results in
      clusters_to_cell_types <- list()

      # Loop through each cluster of input
      for (x in 1:length(each_cluster)) {
        ref_df_new <- ref_df
        # Loop through each row of input and reference
        # First check to see if marker matches, then if the sign matches
        # If both match, add to the 'match' column
        for (a in 1:nrow(each_cluster[[x]])) {
          for (b in 1:nrow(ref_df_new)) {

            ## Matches ##

            if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                each_cluster[[x]]$Sign[a] == "Positive" &
                ref_df_new$Sign[b] == "Positive") {
              ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "Positive" &
                       ref_df_new$Sign[b] == "High") {
              ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "Positive" &
                       ref_df_new$Sign[b] == "Low") {
              ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "Negative" &
                       ref_df_new$Sign[b] == "Negative") {
              ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "Low" &
                       ref_df_new$Sign[b] == "Low") {
              ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "High" &
                       ref_df_new$Sign[b] == "High") {
              ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
            } else {
              ref_df_new$Match[b] <- ref_df_new$Match[b] + 0
            }

            # If the high/low positive box was clicked, then also consider 'low' and 'high' inputs as positive
            if (input$low_option_upload == "low_positive_upload") {
              if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                  each_cluster[[x]]$Sign[a] == "Low" &
                  ref_df_new$Sign[b] == "Positive") {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
              } else {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 0
              }}

            if (input$high_option_upload == "high_positive_upload") {
              if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                  each_cluster[[x]]$Sign[a] == "High" &
                  ref_df_new$Sign[b] == "Positive") {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
              } else {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 0
              }
            }

            # If the low negative box was clicked, then also consider 'low' input as negative
            if (input$low_option_upload == "low_negative_upload") {
              if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                  each_cluster[[x]]$Sign[a] == "Low" &
                  ref_df_new$Sign[b] == "Negative") {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
              } else {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 0
              }}}}

        # Remove unmatched
        ref_df_new <- subset(ref_df_new, Match != 0)


        # Sort, put into final format for export
        final_matched_output <- ref_df_new %>%
          dplyr::group_by(Name) %>%
          dplyr::summarise(
            `Matched inputted markers` = paste0(Marker, " (", Sign, ")", collapse = ", "),
            `num of inputs that match to cell type` = dplyr::n()
          )

        # Put into list
        cluster_name <- each_cluster[[x]]$Cluster[1]
        clusters_to_cell_types[[cluster_name]] <-
          as.data.frame(final_matched_output)
      }

      # Put all results in one data frame
      all_cluster_matches <-
        dplyr::bind_rows(clusters_to_cell_types, .id = "Cluster")


      ## Contradictions ##

      # Make empty, will be added to (+1) if there is a contradiction
      ref_df[, 'Contradiction'] <- 0

      # Empty list to put results in
      clusters_to_cell_types2 <- list()

      # Loop through each cluster of input
      for (x in 1:length(each_cluster)) {
        ref_df_new <- ref_df
        # Loop through each row of input and reference
        # First check to see if marker matches, then if the sign matches
        # If both match, add to the 'match' column
        for (a in 1:nrow(each_cluster[[x]])) {
          for (b in 1:nrow(ref_df_new)) {

            if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                each_cluster[[x]]$Sign[a] == "Positive" &
                ref_df_new$Sign[b] == "Negative") {
              ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "Negative" &
                       ref_df_new$Sign[b] == "Positive") {
              ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "High" &
                       ref_df_new$Sign[b] == "Negative") {
              ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "Negative" &
                       ref_df_new$Sign[b] == "High") {
              ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "Low" &
                       ref_df_new$Sign[b] == "High") {
              ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
            } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                       each_cluster[[x]]$Sign[a] == "High" &
                       ref_df_new$Sign[b] == "Low") {
              ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
            }

            # If the low negative box was clicked, then also consider 'low' input as negative
            if (input$low_option_upload == "low_negative_upload") {
              if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                  each_cluster[[x]]$Sign[a] == "Low" &
                  ref_df_new$Sign[b] == "Positive") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Positive" &
                         ref_df_new$Sign[b] == "Low") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              }}

            # If choosing low/positive, then any low signs are contradicted by negatives and highs
            if (input$low_option_upload == "low_positive_upload") {
              if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                  each_cluster[[x]]$Sign[a] == "Low" &
                  ref_df_new$Sign[b] == "Negative") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Negative" &
                         ref_df_new$Sign[b] == "Low") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              }}

            # If choosing high/high, then any high signs are contradicted by negatives, positives, and lows
            if (input$high_option_upload == "high_high_upload") {
              if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                  each_cluster[[x]]$Sign[a] == "High" &
                  ref_df_new$Sign[b] == "Positive") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Positive" &
                         ref_df_new$Sign[b] == "High") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              }}

            # If choosing low/low, then any low signs are contradicted by negatives, positives, and high
            if (input$low_option_upload == "low_low_upload") {
              if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                  each_cluster[[x]]$Sign[a] == "Low" &
                  ref_df_new$Sign[b] == "Positive") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Positive" &
                         ref_df_new$Sign[b] == "Low") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              }}}}

        # Remove unmatched
        ref_df_new2 <- subset(ref_df_new, Contradiction != 0)


        # Sort, put into final format for export
        final_matched_output2 <- ref_df_new2 %>%
          dplyr::group_by(Name) %>%
          dplyr::summarise(
            `Contradicted markers` = paste0(Marker, " (", Sign, ")", collapse = ", "),
            `num of contradictions` = dplyr::n()
          )

        # Put into list
        cluster_name2 <- each_cluster[[x]]$Cluster[1]
        clusters_to_cell_types2[[cluster_name2]] <-
          as.data.frame(final_matched_output2)
      }

      # Put all results in one data frame
      all_cluster_matches2 <-
        dplyr::bind_rows(clusters_to_cell_types2, .id = "Cluster")

      all_cluster_matches <-
        merge(all_cluster_matches, all_cluster_matches2, by = c("Cluster", "Name"), all.x = TRUE)

      # Remove protein identifier
      ref_df <- ref_df[!ref_df$Marker == "protein", ]

      # Get what and how many total markers in each cell type
      final_matched_output <- ref_df %>%
        dplyr::group_by(Name) %>%
        dplyr::summarise(
          `Cell type markers` = paste0(Marker, " (", Sign, ")", collapse = ", "),
          `total number of markers` = dplyr::n()
        )

      all_cluster_matches <-
        merge(all_cluster_matches, final_matched_output, by = "Name", all.x = TRUE)

      # Get how many total markers in each input
      #  num_markers_per_input <- as.data.frame(table(new_df$Cluster))
      # colnames(num_markers_per_input) <- c("Cluster", "input_num_total")

      # Get what and how many total markers in each input
      final_matched_output <- new_df %>%
        dplyr::group_by(Cluster) %>%
        dplyr::summarise(
          `Inputted markers` = paste0(Marker, " (", Sign, ")", collapse = ", "),
          input_num_total = dplyr::n()
        )

      all_cluster_matches <-
        merge(all_cluster_matches, final_matched_output, by = "Cluster")


      # Count how many markers are positive and negative
      all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_pos_input = stringr::str_count(`Inputted markers`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
      all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_neg_input = stringr::str_count(`Inputted markers`, '(Negative)'))

      all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_pos_cell_type = stringr::str_count(`Cell type markers`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
      all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_neg_cell_type = stringr::str_count(`Cell type markers`, '(Negative)'))

      all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_pos_match = stringr::str_count(`Matched inputted markers`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
      all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_neg_match = stringr::str_count(`Matched inputted markers`, '(Negative)'))

      all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_pos_contr = stringr::str_count(`Contradicted markers`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
      all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_neg_contr = stringr::str_count(`Contradicted markers`, 'Negative'))

      # Make empty values in number of contradictions NA values
      all_cluster_matches$`num of contradictions`[is.na(all_cluster_matches$`num of contradictions`)] <- 0
      all_cluster_matches$Num_pos_contr[is.na(all_cluster_matches$Num_pos_contr)] <- 0
      all_cluster_matches$Num_neg_contr[is.na(all_cluster_matches$Num_neg_contr)] <- 0

      # Percent match (match markers/markers in cell type)
      all_cluster_matches$`Percent match (match markers/markers in cell type)` <-
        (as.numeric(all_cluster_matches$`num of inputs that match to cell type`) /
           (as.numeric(all_cluster_matches$`total number of markers`)
           )) * 100

      all_cluster_matches$`Percent match (match markers/markers in cell type)` <- as.numeric(all_cluster_matches$`Percent match (match markers/markers in cell type)`)

      all_cluster_matches$`Percent match (match markers/markers in cell type)` <- round(all_cluster_matches$`Percent match (match markers/markers in cell type)`, 3)

      # Percent match (match markers/markers in input)
      all_cluster_matches$`Percent match (match markers/markers in input)` <-
        ((as.numeric(all_cluster_matches$`num of inputs that match to cell type`) /
            (as.numeric(all_cluster_matches$input_num_total))) * 100)

      all_cluster_matches$`Percent match (match markers/markers in input)` <- as.numeric(all_cluster_matches$`Percent match (match markers/markers in input)`)

      all_cluster_matches$`Percent match (match markers/markers in input)` <- round(all_cluster_matches$`Percent match (match markers/markers in input)`, 3)

      # Score
      all_cluster_matches$Score <-
        ((as.numeric(all_cluster_matches$Num_pos_match)) /
           (as.numeric(all_cluster_matches$Num_pos_cell_type) + as.numeric(all_cluster_matches$Num_pos_input) - as.numeric(all_cluster_matches$Num_pos_match))) - ((as.numeric(all_cluster_matches$Num_pos_contr) * 0.25) + (as.numeric(all_cluster_matches$Num_neg_contr) * 0.25))

      all_cluster_matches$Score <- as.numeric(all_cluster_matches$Score)

      all_cluster_matches$Score <- round(all_cluster_matches$Score, 3)

      # Make sure there are results
      if (nrow(all_cluster_matches) != 0) {

        # Get the columns for output, new order
        all_cluster_matches <-
          all_cluster_matches[, c("Cluster",
                                  "Name",
                                  "Matched inputted markers",
                                  "num of inputs that match to cell type",
                                  "Inputted markers",
                                  "input_num_total",
                                  "Cell type markers",
                                  "total number of markers",
                                  "Contradicted markers",
                                  "num of contradictions",
                                  "Percent match (match markers/markers in cell type)",
                                  "Percent match (match markers/markers in input)",
                                  "Score"
          )]

        # Final name
        colnames(all_cluster_matches) <-
          c(
            "Cluster",
            "Cell type name",
            "Matched inputted markers",
            "# matched inputted markers",
            "Inputted markers",
            "# inputted markers",
            "Cell type markers",
            "# cell type markers",
            "Contradictions",
            "# contradictions",
            "% (matched markers/markers in cell type)",
            "% (matched markers/markers in input)",
            "Score"
          )

        # Order and sort by cluster
        all_cluster_matches <- all_cluster_matches[with(all_cluster_matches, order(Cluster, -Score)), ]

        return(all_cluster_matches)
      }  else {
        # If no matches, return null
        return(NULL)
      }
    } else {
      # If no matches, return null
      return(NULL)
    }
  })

  # Step 4,  uploaded reference, if multiple clusters, change dataframe if filtered by cluster
  df_subsettable_uploaded_1 <- reactive({
  shiny::req(input$submit_tab1_step3)
   
    df <- uploaded_ref_results_1()

    if (!is.null(df) & length(unique(df$Cluster)) > 1) {
      df2 <- df[df$Cluster %in% input$show_specific_cluster_uploaded_1,]
      return(df2)
    } else {
      return(df)
    }
  })

  # Step 4, uploaded reference, get how many clusters are in final output for filtering purposes
  output$specific_cluster_uploaded_1 <- shiny::renderUI({
     shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)
   
    final_df <- uploaded_ref_results_1()
    if (!is.null(final_df) & length(unique(final_df$Cluster)) > 1) {
      shiny::selectInput(inputId = "show_specific_cluster_uploaded_1",
                         label = "Filter by cluster:  ",
                         selected = sort(unique(as.character(final_df$Cluster))),
                         choices = sort(unique(as.character(final_df$Cluster))),
                         multiple = TRUE)
    }
  })

  # Step 4, uploaded reference, changes when the 'submit' button is hit if using an uploaded reference
  shiny::observeEvent(input$submit_tab1_step3, {
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)

    shiny.destroy::removeOutput("help_tab1_step3")
    shiny.destroy::removeOutput("spaces1")
    shiny.destroy::removeOutput("spaces2")
    shiny.destroy::removeOutput("spaces3")
    shiny.destroy::removeOutput("overview_tab1_step3")
   shinyjs::show("overview_uploaded_ref_1")
   
    shiny.destroy::removeOutput("marker_diff_eq_heatmap_raw_row_size2")
    shiny.destroy::removeOutput("marker_diff_eq_heatmap_raw_column_size2")
    shiny.destroy::removeOutput("heatmap_height_size2")
    shiny.destroy::removeOutput("density_row_size2")
    shiny.destroy::removeOutput("density_height_size2")
    shiny.destroy::removeOutput("density_facet2")
    shiny.destroy::removeOutput("marker_diff_eq_results_editable")
    shiny.destroy::removeOutput("marker_diff_eq_heatmap_raw_ex_plot")
    shiny.destroy::removeOutput("expr_heatmap_org_ex_plot")
    shiny.destroy::removeOutput("expr_density_ex_plot")
    shiny.destroy::removeOutput("cutoff_parameters")
    shiny.destroy::removeOutput("download_marker_diff_eq_csv")

    # Step 4, uploaded reference, make datatable of matched cell types
    output$uploaded_ref_matches_1_dt <- DT::renderDT({
      if (!is.null(uploaded_ref_results_1()) & nrow(uploaded_ref_results_1()) == 0) {
      } else {
        DT::datatable(
          df_subsettable_uploaded_1(), #uploaded_ref_results_1(),
          selection = 'none',
          rownames = FALSE,
          escape = FALSE,
          options = list(pageLength = 25)
        )
      }
    }, server=FALSE)

    if (is.null(check_reference_file_1()) == TRUE) {

      if (is.data.frame(reformatted_data_1()) == TRUE &
          is.data.frame(reformatted_ref_1()) == TRUE) {

        if (!is.null(uploaded_ref_matches_1()) & nrow(uploaded_ref_matches_1()) > 0) {
          shinyjs::show("help_uploaded_ref_1")
          shiny.destroy::removeOutput("all_sidebar_tab1_step3")
          shinyjs::show("all_sidebar_tab1_step6")
          shinyjs::show("download_uploaded_ref_table_cell_types_1")
        } else {
          shinyjs::show("no_matched_cell_types_1")
        }
      }
    }
  })

  uploaded_ref_results_1 <- reactive({
    shiny::req(input$reference_type_1 == 'upload_ref_1')
    shiny::req(input$upload_ref_csv_1)
  # shiny::req(input$submit_tab1_step3)

    matches <- uploaded_ref_matches_1()

    # Remove number columns (still included in download)
    matches <-
      matches[,     c(
        "Cluster",
        "Cell type name",
        "Matched inputted markers",
        "Inputted markers",
        "Cell type markers",
        "Contradictions",
        "% (matched markers/markers in cell type)",
        "% (matched markers/markers in input)",
        "Score"
      )]

    # Make look nicer (round and add a % symbol)
    matches$`% (matched markers/markers in cell type)` <-
      paste(round(matches$`% (matched markers/markers in cell type)`, 1), "%")

    matches$`% (matched markers/markers in input)` <-
      paste(round(matches$`% (matched markers/markers in input)`, 1), "%")

    # Order by score
    matches <- matches[with(matches, order(Cluster, -Score)), ]

    return(matches)

  })

  # Step 4, uploaded reference, download final results
  output$download_uploaded_ref_table_cell_types_1 <- shiny::downloadHandler(
    # File name
    filename = function() {
      if (input$input_type == 'csv_upload') {
        paste0(gsub(".csv|.xls|.xlsx", "", input$csv_upload),
               "-matched-cell-types.csv")
      } else {
        paste0("typed-markers-matched-cell-types.csv")
      }
    },

    # Write to csv
    content = function(file) {
    utils::write.csv(uploaded_ref_matches_1(),  file, row.names = FALSE)
    }
  )

  #~#~#~#~#~# Internal references #~#~#~#~#~#

  # Step 4, when help button is clicked show message
  shiny::observeEvent(input$help_tab1_step4,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 4<br>Match Markers to PRO/GO Terms"),
      HTML("In this step, the inputted markers are processed through a multistep workflow to find matching ontology terms. If a marker cannot be automatically matched, the tool will return it to the user for manual name edits.
      <br>
      <br>
      A data frame is displayed, showing the original marker names and any suggested names. The last column, 'Type New Marker Name', is editable, allowing users to propose new names. These new names will be reprocessed to check for matches within the multistep workflow.
      <br>
      <br>"),
      tags$div("Additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/3.1.-Input-Expression-Data-and-Default-References', "here.", target="_blank")),
      easyClose = TRUE))
  })

  # Step 4, create editable dataframes
  df_1 <- shiny::reactiveValues(data = NULL)

  # Step 4, changes when the 'submit' button is hit if using the internal references
  shiny::observeEvent(input$submit_tab1_step3, {
    shiny::req(input$reference_type_1 == 'internal_ref_1')

    shiny.destroy::removeOutput("help_tab1_step3")
    shiny.destroy::removeOutput("spaces1")
    shiny.destroy::removeOutput("spaces2")
    shiny.destroy::removeOutput("spaces3")
    shiny.destroy::removeOutput("overview_tab1_step3")
    shiny.destroy::removeOutput("marker_diff_eq_results_editable")
    shiny.destroy::removeOutput("marker_diff_eq_heatmap_raw_ex_plot")
    shiny.destroy::removeOutput("expr_heatmap_org_ex_plot")
    shiny.destroy::removeOutput("expr_density_ex_plot")
    shiny.destroy::removeOutput("cutoff_parameters")
    shiny.destroy::removeOutput("marker_diff_eq_heatmap_raw_row_size2")
    shiny.destroy::removeOutput("marker_diff_eq_heatmap_raw_column_size2")
    shiny.destroy::removeOutput("heatmap_height_size2")
    shiny.destroy::removeOutput("density_row_size2")
    shiny.destroy::removeOutput("density_height_size2")
    shiny.destroy::removeOutput("density_facet2")
    shiny.destroy::removeOutput("download_marker_diff_eq_csv")

    shinyjs::show("overview_tab1_step4")

    if (length(marker_sign_error_1()) == 0) {

      df_1$data <- return_alt_markers_1()

      # Step 4, make datatable of unmatched markers to PRO terms
      output$table_new_format_1 <- DT::renderDT({
        DT::datatable(
          df_1$data,
          selection = 'none',
          editable = list(target = "cell", disable = list(columns = c(0, 1))),
          rownames = FALSE,
          escape = FALSE,
          options = list(pageLength = 25)
        )
      }, server = FALSE)

      if (is.data.frame(reformatted_data_1()) == TRUE) {
        shinyjs::show("help_tab1_step4")
        shiny.destroy::removeOutput("all_sidebar_tab1_step3")
        shinyjs::show("all_sidebar_tab1_step4")
      }
    }
  })

  # Step 4, match inputted markers to PRO terms
  PRO_results_1 <- reactive({
  shiny::req(input$submit_tab1_step3)
  shiny::req(input$reference_type_1 == 'internal_ref_1')

    if (is.data.frame(reformatted_data_1()) == TRUE) {

      ## Species, references, initial dataframe

      new_df <- reformatted_data_1()

      # Get NCBI taxon ID
      species_ID <- input$species_choice_1
      matched_row <-
        which(species_in_PRO == species_ID, arr.ind = TRUE)[1]
      NCBI_taxon_ID <-
        species_in_PRO[matched_row, ncol(species_in_PRO)]

      NCBI_taxon_ID_full <- species_in_PRO[matched_row, 1]
      NCBI_taxon_ID_short <- gsub(".*/","",NCBI_taxon_ID_full)

      # Get Wikidata taxon ID
      SPARQL_query <-

        paste0(
          "SELECT ?item ?itemLabel WHERE {

     ?item wdt:P2888 <", NCBI_taxon_ID_full, ">.

    SERVICE wikibase:label { bd:serviceParam wikibase:language 'en'. }
    }"
        )

      # Run SPARQL on the endpoint
      result <- httr::GET(
        url = wiki_endpoint,
        query = list(query = SPARQL_query),
        httr::user_agent(R.version.string))

      # Will show a warning/error if there is any
      httr::stop_for_status(result)

      # Get result in text JSON
      x <- httr::content(result, as = "text") #, encoding = "UTF-8")

      # Convert from JSON to a list
      df <- jsonlite::fromJSON(x, flatten = TRUE)

      # Extract the data frame
      df <- df$results$bindings

      if (length(df) != 0) {

        # Remove unneeded info
        df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", ':lang')))

        # Remove this part that was added onto the column names
        colnames(df) <- gsub(".value", "", colnames(df))
      }

      wikidata_taxon_ID <- df$item

      # List of specific protein suggestions
      marker_alt <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/specific_suggestions_V1.csv", header = TRUE)

      # List of GO complexes that are included in CL definitions
      GO_complexes <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/GO_complexes_V1.csv", header = TRUE)

      ## List of CD synonym (adapted from list published by HCDM)
      CD_syn <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/Curated_synonyms_V1.csv", header = TRUE)

      # Auto replacements
      replace_markers <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/replace_markers_V1.csv", header=TRUE)

      ## Lists of non-specific suggestions
      # List of non-protein words that were found in markers from the deveolpment data
      common_non_protein_words <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/common_non_protein_words_V1.csv", header = FALSE)

      # List of metal tags that were found in the development data
      commonly_used_metal_tags <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/commonly_used_metal_tags_V1.csv", header = FALSE)

      # List of fluorophores tags that were  found in the development data
      commonly_used_fluorophores <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/commonly_used_fluorophores_V1.csv", header = FALSE)

      # Make markers uppercase and remove dashes, underscores, periods, and spaces
      new_df$Marker2 <- toupper(new_df$Marker)
      new_df$Marker2 <- gsub('-', '', new_df$Marker2, fixed = TRUE)
      new_df$Marker2 <- gsub('_', '', new_df$Marker2, fixed = TRUE)
      new_df$Marker2 <- gsub(' ', '', new_df$Marker2, fixed = TRUE)
      new_df$Marker2 <- gsub('.', '', new_df$Marker2, fixed = TRUE)

      # Make empty dataframe
      auto_matches <- data.frame(Marker = unique(new_df$Marker2))

      auto_matches$PRO_term <- ""
      auto_matches$PRO_name <- ""
      auto_matches$Species <- ""
      auto_matches$Match_type <- ""
      auto_matches$Match_step <- ""
      auto_matches$Suggestions <- ""

      ## Automatic matching (CD3e, CD8a, TCR)

      # Automatically make any inputs of CD3e, CD8a, and TCR match to the wider complex

         # If TCR is also specified in the panel, don't auto match CD3E to it
    if (any(auto_matches=="TCRAB") | any(auto_matches=="TCRA/B" ) | 
        any(auto_matches=="TCRABCOMPLEX") | any(auto_matches=="TCRA/BCOMPLEX" ) |
        any(auto_matches=="ABTCR") | any(auto_matches=="ABTCRCOMPLEX" ) |
        any(auto_matches=="ALPHABETATCR") | any(auto_matches=="ALPHABETATCRCOMPLEX" ) |
        any(auto_matches=="TCRALPHABETA") | any(auto_matches=="TCRALPHABETACOMPLEX" ) |
        any(auto_matches=="TCRA") | any(auto_matches=="TCRB" ) |
        any(auto_matches=="TCRALPHA") | any(auto_matches=="TCRBETA" ) |
        any(auto_matches=="ALPHATCR") | any(auto_matches=="BETATCR" ) |
        any(auto_matches=="TCRGD") | any(auto_matches=="TCRGDCOMPLEX" ) |
        any(auto_matches=="TCRG/D") | any(auto_matches=="TCRG/DCOMPLEX" ) |
        any(auto_matches=="GDTCR") | any(auto_matches=="GDTCRCOMPLEX" ) |
        any(auto_matches=="GAMMADELTATCR") | any(auto_matches=="GAMMADELTATCRCOMPLEX" ) |
        any(auto_matches=="TCRGAMMADELTA") | any(auto_matches=="TCRGAMMADELTACOMPLEX" ) |
        any(auto_matches=="TCRG") | any(auto_matches=="TCRD" ) |
        any(auto_matches=="TCRGAMMA") | any(auto_matches=="TCRDELTA" ) |
        any(auto_matches=="GAMMATCR") | any(auto_matches=="DELTATCR" ) |
        any(auto_matches=="YTCR") | any(auto_matches=="TCRY" ) |
        any(auto_matches=="TCRYD") | any(auto_matches=="YDTCR" ) |
        any(auto_matches=="TCR") | any(auto_matches=="TCRCOMPLEX" ) |
        any(auto_matches=="TCELLRECEPTOR") | any(auto_matches=="TCELLRECEPTORCOMPLEX" )) {
      
      replace_markers <- replace_markers[!grepl("TCRAB", replace_markers$Replacement),]
      replace_markers <- replace_markers[!grepl("TCRGD", replace_markers$Replacement),]
      replace_markers <- replace_markers[!grepl("TCR", replace_markers$Replacement),]
    }
    
    # If CD8AB is also specified in the panel, don't auto match CD8 to it
    if (any(auto_matches=="CD8AB") | any(auto_matches=="CD8ALPHABETA")) {
      replace_markers <- replace_markers[!grepl("CD8ALPHABETA", replace_markers$Replacement),]
    }
     
      auto_matches <- merge(replace_markers, auto_matches,  by = "Marker", all.y = TRUE)

      auto_replaced <- auto_matches[!is.na(auto_matches$Replacement),]

      auto_matches$Replacement[is.na(auto_matches$Replacement)] <- auto_matches$Marker[is.na(auto_matches$Replacement)]

      auto_matches$Marker <- NULL

      colnames(auto_matches)[1] <- "Marker"

      ## Punctuation

      ## Check if these punctuation marks are in the input.

      Punctuation_matches <- auto_matches

      # Slash
      Punctuation_matches$slash <- ifelse(grepl("/", Punctuation_matches$Marker, useBytes = TRUE), "A slash '/' was detected in the marker name.", "")

      # Backslash
      Punctuation_matches$backslash <- ifelse(grepl("\\\\", Punctuation_matches$Marker, useBytes = TRUE), "A backslash '\' was detected in the marker name.", "")

      # Comma
      Punctuation_matches$comma <- ifelse(grepl(",", Punctuation_matches$Marker, useBytes = TRUE), "A comma ',' was detected in the marker name.", "")

      # Bracket
      Punctuation_matches$bracket <- ifelse(grepl("[\\[\\]", Punctuation_matches$Marker, useBytes = TRUE), "A bracket '[' was detected in the marker name.", "")

      # Parenthesis
      Punctuation_matches$parenthesis <- ifelse(grepl("[\\(\\)]", Punctuation_matches$Marker, useBytes = TRUE), "A paranthesis '(' was detected in the marker name.", "")

      # Paste together as suggestions
      Punctuation_matches$Suggestions <- paste(Punctuation_matches$slash, Punctuation_matches$backslash, Punctuation_matches$comma, Punctuation_matches$bracket, Punctuation_matches$parenthesis)
      Punctuation_matches$Suggestions <- gsub("^ *|(?<= ) | *$", "", Punctuation_matches$Suggestions, perl = TRUE)

      # Remove puncationation columns
      Punctuation_matches <- subset(Punctuation_matches, select = -c(slash, comma, bracket, parenthesis, backslash))

      # Indicate that marker was flagged for punctuation to change
      for(i in 1:nrow(Punctuation_matches)) {
        if(Punctuation_matches$Suggestions[i] != "") {
          Punctuation_matches$Match_step[i] <- "Punctuation suggestions"
        }
      }

      ## Non-specific Suggestions

      # If marker name partial matches to an entry on one of these lists, give specific suggestions on how to rewrite that name

      # Vector of unique marker names
      Non_protein_sugg_matches <- Punctuation_matches[Punctuation_matches$Match_step != "Punctuation suggestions",]

      if(length(Non_protein_sugg_matches != 0)) {

        # List of words that are commonly used in non-protein marker names/fcs headers
        Non_protein_sugg_matches$Common_non_protein_words <- ifelse(grepl(paste0(common_non_protein_words$V1, collapse = "|"), Non_protein_sugg_matches$Marker, useBytes = TRUE), "A non-protein word was detected in the marker name.", "")

        # List of metals
        Non_protein_sugg_matches$commonly_used_metal_tags <- ifelse(grepl(paste0(commonly_used_metal_tags$V1, collapse = "|"), Non_protein_sugg_matches$Marker, useBytes = TRUE), "A metal tag was detected in the marker name.", "")

        # List of fluorophores
        Non_protein_sugg_matches$commonly_used_fluorophores <- ifelse(grepl(paste0(commonly_used_fluorophores$V1, collapse = "|"), Non_protein_sugg_matches$Marker, useBytes = TRUE), "A fluorophore was detected in the marker name.", "")

        # Put all suggestions together
        Non_protein_sugg_matches$Suggestions <- paste(Non_protein_sugg_matches$Common_non_protein_words, Non_protein_sugg_matches$commonly_used_metal_tags, Non_protein_sugg_matches$commonly_used_fluorophores, Non_protein_sugg_matches$species_specific_letter)
        Non_protein_sugg_matches$Suggestions <- gsub("^ *|(?<= ) | *$", "", Non_protein_sugg_matches$Suggestions, perl = TRUE)

        for(i in 1:nrow(Non_protein_sugg_matches)) {
          if(Non_protein_sugg_matches$Suggestions[i] != "") {
            Non_protein_sugg_matches$Match_step[i] <- "Non-protein word suggestions"
          }
        }

        Non_protein_sugg_matches <- subset(Non_protein_sugg_matches, select = -c(Common_non_protein_words, commonly_used_metal_tags, commonly_used_fluorophores))
      }

      ## GO Complexes

      # Match entries to the GO complexes listed in the CL

      # Vector of unique marker names
      GO_matches <- Non_protein_sugg_matches[Non_protein_sugg_matches$Match_step != "Non-protein word suggestions",]

      unique_marker_names <- GO_matches$Marker

      all_ids <- list()

      for(marker2 in unique_marker_names) {

        matched_row <- GO_complexes[GO_complexes$Inputted_Name ==marker2,]

        if (nrow(matched_row) != 0) {

          GO_matches[GO_matches$Marker == marker2,]$Suggestions <- ""
          GO_matches[GO_matches$Marker == marker2,]$Match_step <- "GO exceptions"
          GO_matches[GO_matches$Marker == marker2,]$PRO_term <- matched_row$Matched_GO_term
          GO_matches[GO_matches$Marker == marker2,]$PRO_name <- matched_row$Matched_Descriptive_GO_name
        }
      }

      ## Direct SPAPRQL query to PRO

      # Directly query PRO

      Match_step <- "Directly to PRO"

      unique_marker_names <- GO_matches$Marker[GO_matches$Match_step != "GO exceptions"]

      # Run SPARQL
      all_ids <- PRO_SPARQL(unique_marker_names = unique_marker_names,
                            NCBI_taxon_ID = NCBI_taxon_ID_short,
                            onto_endpoint = onto_endpoint,
                            Match_step = Match_step)

      # Put all results in one data frame
      Direct_PRO_matches <- dplyr::bind_rows(all_ids, .id = "Marker")

      ## Specific Suggestions

      # If marker name matches to an entry on this list, give specific suggestions on how to rewrite that name

      # Vector of unique marker names
      Specific_protein_sugg_matches <- Direct_PRO_matches[Direct_PRO_matches$Match_step != "Directly to PRO",]

      unique_marker_names <- Specific_protein_sugg_matches$Marker

      all_ids <- list()

      if(length(unique_marker_names) != 0) {

        Specific_protein_sugg_matches$Suggestions <- ""

        for(marker2 in unique_marker_names) {

          matched_row <- marker_alt[marker_alt$CD_NAME == marker2,]

          if (nrow(matched_row) != 0) {

            remove_empties <- sapply(matched_row, function(x) all(is.na(x) | x == ""))
            matched_row <- matched_row[, !remove_empties]

            suggestions <- paste(matched_row[, -1], collapse = ", ")

            Specific_protein_sugg_matches[Specific_protein_sugg_matches$Marker == marker2,]$Suggestions <- paste("Try:", suggestions)
            Specific_protein_sugg_matches[Specific_protein_sugg_matches$Marker == marker2,]$Match_step <- "Specific protein suggestions"
            Specific_protein_sugg_matches[Specific_protein_sugg_matches$Marker == marker2,]$PRO_term <- ""
            Specific_protein_sugg_matches[Specific_protein_sugg_matches$Marker == marker2,]$PRO_name <- ""
          }
        }
      } else {
        Specific_protein_sugg_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
      }

      ## CD synonym list, then SPAPRQL query to PRO

      # Search CD synonym list, then if there are synonyms query them in PRO

      Match_step <- "CD synonym list"

      # Vector of unique marker names
      unique_marker_names <- Specific_protein_sugg_matches$Marker[Specific_protein_sugg_matches$Match_step != "Specific protein suggestions"]

      all_ids_final <- list()

      if(length(unique_marker_names) != 0) {

        for (marker in unique_marker_names) {

          matched_row <- CD_syn[which(CD_syn == marker, arr.ind=TRUE)[,1],]

          matched_row <- data.frame(x=unlist(matched_row))

          unique_syn_marker_names <- unique(stats::na.omit(matched_row))

          unique_syn_marker_names <- as.vector(unique_syn_marker_names[,1])

          # Run SPARQL
          all_ids <- PRO_SPARQL(unique_marker_names = unique_syn_marker_names,
                                NCBI_taxon_ID = NCBI_taxon_ID_short,
                                onto_endpoint = onto_endpoint,
                                Match_step = Match_step)

          # Put all results in one data frame
          all_ids <- dplyr::bind_rows(all_ids, .id = "Matched_synonym")

          all_ids <- stats::na.omit(all_ids)

          if (length(all_ids) != 0) {
            all_ids_final[[marker]] <- all_ids
          } else if (length(all_ids) == 0) {
            all_ids_final[[marker]] <- data.frame(Matched_synonym = "", PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
          }
        }

        # Put all results in one data frame
        all_ids_final2 <- dplyr::bind_rows(all_ids_final, .id = "Marker")

        # Remove duplicated entries (Marker, PRO term, and type of match are the same)
        all_ids_final2 <- all_ids_final2[!duplicated(all_ids_final2[,c("Marker","PRO_term", "Match_type")]),]

        # Get entries were inputted marker and matched species is the same

        # Subset into entries with and without duplicates
        x <-  all_ids_final2[duplicated(all_ids_final2[,c("Marker","Species")]) | duplicated(all_ids_final2[,c("Marker","Species")], fromLast = TRUE) ,]
        y <-  dplyr::setdiff(all_ids_final2, x)

        # Get vector of markers (duplicates)
        unique_marker_names <- unique(x$Marker)

        # If there is both an exact match and a secondary match, remove the secondary
        remove_secondary_per_marker <- data.frame()

        for(marker in unique_marker_names) {
          df_subsetted <- x[x$Marker == marker,]
          if(any(df_subsetted$Match_type == "Primary")) {
            df_per_marker <- df_subsetted[df_subsetted$Match_type != "Secondary",]
            df_per_marker <- df_per_marker[df_per_marker$Match_type != "",]
            remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
          } else if(any(df_subsetted$Match_type == "Secondary")) {
            df_per_marker <- df_subsetted[df_subsetted$Match_type != "",]
            remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
          }
        }

        # Final
        CD_syn_list_matches <- rbind(remove_secondary_per_marker, y)

        # Just columns we want
        CD_syn_list_matches <- unique(CD_syn_list_matches[ , c("Marker", "PRO_term", "PRO_name", "Species", "Match_type", "Match_step")])

      } else {
        CD_syn_list_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
      }

      ## Wikidata, then SPAPRQL query to PRO

      # Query Wikidata, then if there are synonyms query them in PRO

      Match_step <- "Wikidata"

      # Vector of unique marker names
      unique_marker_names <- CD_syn_list_matches$Marker[CD_syn_list_matches$Match_step != "CD synonym list"]

      # Remove duplicates introduces by synonyms from last step
      remove <- CD_syn_list_matches$Marker[CD_syn_list_matches$Match_step != ""]
      unique_marker_names <-  dplyr::setdiff(unique_marker_names, remove)

      all_ids_final <- list()

      if(length(unique_marker_names) != 0){
        if(unique_marker_names[1] != "") {

          # Loop through each and do the SPARQL query
          for(marker in unique_marker_names) {

            specific_marker <-  shQuote(marker)

            SPARQL_query <-
              paste0(
                "
SELECT DISTINCT ?item ?itemLabel ?altLabel WHERE {
  # Look either for non-species specific or add a specific species
  OPTIONAL {
    ?item wdt:P703 <", wikidata_taxon_ID, "> .

  }
  # Get Wikidata protein name and alternative names
  ?item rdfs:label ?itemLabel .
  ?item skos:altLabel ?altLabel.
  # Return results that are either proteins or protein-coding genes
  {?item wdt:P31 wd:Q8054}
  UNION
  {?item wdt:P279 wd:Q20747295}
  # Search for the marker input, either if it is the entry name or an alternative name
  FILTER (UCASE(REPLACE(str(?itemLabel),'[ -.]','')) = ", specific_marker, " || UCASE(REPLACE(str(?altLabel),'[ -.]','')) = ", specific_marker, ")
}
")

            # Run SPARQL on the endpoint
            result <- httr::GET(
              url = wiki_endpoint,
              query = list(query = SPARQL_query),
              httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result)

            # Get result in text JSON
            x <- httr::content(result, as = "text") #, encoding = "UTF-8")

            # Convert from JSON to a list
            df <- jsonlite::fromJSON(x, flatten = TRUE)

            # Extract the data frame
            df <- df$results$bindings

            if (length(df) != 0) {

              # Remove unneeded info
              df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", 'lang')))

              # Remove this part that was added onto the column names
              colnames(df) <- gsub(".value", "", colnames(df))

              # Make everything case, remove spaces, dashes
              df <- dplyr::mutate_all(df, .funs = toupper)
              df <- as.data.frame(lapply(df, function(y) gsub('-', '', y)))
              df <- as.data.frame(lapply(df, function(y) gsub('_', '', y)))
              df <- as.data.frame(lapply(df, function(y) gsub(' ', '', y)))
              df <- as.data.frame(lapply(df, function(y) gsub('\\.', '', y)))

              df <- df[(which(nchar(df$altLabel) > 2)),]
              df <- df[(which(nchar(df$itemLabel) > 2)),]

              # No marker can start with a number, so add X to any marker names that start with a number
              df <- as.data.frame(lapply(df, function(y) gsub('^(\\d)', 'X\\1', y)))

              # Vector of unique marker names
              unique_syn_names <- unique(c(df[, "altLabel"], df[, "itemLabel"]))

              unique_syn_names <- unique_syn_names[!grepl("[^\x01-\x7F]+", unique_syn_names)]

              unique_syn_names <- unique_syn_names[grep('[A-Za-z]',unique_syn_names)]

              unique_syn_names <- unique_syn_names[!grepl(paste0('^', marker, '$'), unique_syn_names, useBytes = TRUE)]

              # Run SPARQL
              all_ids <- PRO_SPARQL(unique_marker_names = unique_syn_names,
                                    NCBI_taxon_ID = NCBI_taxon_ID_short,
                                    onto_endpoint = onto_endpoint,
                                    Match_step = Match_step)

              # Put all results in one data frame
              all_ids <- dplyr::bind_rows(all_ids, .id = "Matched synonym")

              all_ids <- stats::na.omit(all_ids)

              if (length(all_ids) != 0) {
                all_ids_final[[marker]] <- all_ids
              } else if (length(all_ids) == 0) {
                all_ids_final[[marker]] <- data.frame(PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
              }
            } else if (length(df) == 0) {
              all_ids_final[[marker]] <- data.frame(PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
            }
          }

        # Put all results in one data frame
        all_ids_final2 <- dplyr::bind_rows(all_ids_final, .id = "Marker")

        # Remove duplicated entries (Marker, PRO term, and type of match are the same)
        all_ids_final2 <- all_ids_final2[!duplicated(all_ids_final2[,c("Marker","PRO_term", "Match_type")]),]

        # Get entries were inputted marker and matched species is the same

        # Subset into entries with and without duplicates
        x <-  all_ids_final2[duplicated(all_ids_final2[,c("Marker","Species")]) | duplicated(all_ids_final2[,c("Marker","Species")], fromLast = TRUE) ,]
        y <-  dplyr::setdiff(all_ids_final2, x)

        # Get vector of markers (duplicates)
        unique_marker_names <- unique(x$Marker)

        # If there is both an exact match and a secondary match, remove the secondary
        remove_secondary_per_marker <- data.frame()

        for(marker in unique_marker_names) {
          df_subsetted <- x[x$Marker == marker,]
          if(any(df_subsetted$Match_type == "Primary")) {
            df_per_marker <- df_subsetted[df_subsetted$Match_type != "Secondary",]
            df_per_marker <- df_per_marker[df_per_marker$Match_type != "",]
            remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
          } else if(any(df_subsetted$Match_type == "Secondary")) {
            df_per_marker <- df_subsetted[df_subsetted$Match_type != "",]
            remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
          }
        }

        # Final
        Wikidata_matches <- rbind(remove_secondary_per_marker, y)

        # Just columns we want
        Wikidata_matches <- unique(Wikidata_matches[ , c("Marker", "PRO_term", "PRO_name", "Species", "Match_type", "Match_step")])

        } else {
          Wikidata_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
        }
      } else {
        Wikidata_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
      }

    ## Final results

    # Put everything together and get the final list of matched/unmatched marker names

    # Add empty suggestion column to parts where suggestions were not relevant
    Direct_PRO_matches$Suggestions <- ""
    CD_syn_list_matches$Suggestions <- ""
    Wikidata_matches$Suggestions <- ""

    # Subset entries that were matched at each step
    Punctuation_matches2 <- Punctuation_matches[Punctuation_matches$Match_step == "Punctuation suggestions",]
    Non_protein_sugg_matches2 <- Non_protein_sugg_matches[Non_protein_sugg_matches$Match_step == "Non-protein word suggestions",]
    GO_matches2 <- GO_matches[GO_matches$Match_step == "GO exceptions",]
    Direct_PRO_matches2 <- Direct_PRO_matches[Direct_PRO_matches$Match_step == "Directly to PRO",]
    CD_syn_list_matches2 <- CD_syn_list_matches[CD_syn_list_matches$Match_step == "CD synonym list",]
    Wikidata_matches2 <- Wikidata_matches[Wikidata_matches$Match_step == "Wikidata",]
    Specific_protein_sugg_matches2 <- Specific_protein_sugg_matches[Specific_protein_sugg_matches$Match_step == "Specific protein suggestions",]

    # Put everything together into 1 dataframe
    all_ids <- do.call("rbind", list(Punctuation_matches2, Non_protein_sugg_matches2, GO_matches2,  Direct_PRO_matches2, CD_syn_list_matches2, Wikidata_matches2, Specific_protein_sugg_matches2))

    # Put original inputted names back in
    final_df <- merge(new_df, all_ids, by.x = "Marker2", by.y = "Marker", all.x = TRUE, all.y = TRUE)

    # Change names
    colnames(final_df)[1:2] <- c("New_Name", "Original_Name")

    if (nrow(auto_replaced) != 0) {
      # Put back original auto changed names
      replacements <- auto_replaced[as.character(auto_replaced$Marker) != as.character(auto_replaced$Replacement),]

      replacements <- replacements[,1:2]

      colnames(replacements) <- c("Marker2", "New_Name")

      replacements <- merge(replacements, new_df, by = "Marker2")

      replacements <- replacements[,c(2,3)]

      colnames(replacements) <- c("New_Name", "Original_Name")

      replacements <- merge(replacements, final_df, by = "New_Name")

      replacements <- replacements[,-3]

      colnames(replacements)[2] <- "Original_Name"

      final_df <- rbind(final_df, replacements)
    }

    final_df <- final_df[!is.na(final_df$Original_Name), ]

    # If there was no match or suggestion, indicate that
    final_df$Match_step[is.na(final_df$Match_step)] <- "Did not match"

    # Turn any remaining NA entries into blanks
    final_df[is.na(final_df)] <- ""

    # Change species URIs into human readable form
    final_df$Species <- gsub(NCBI_taxon_ID_full, input$species_choice_1, final_df$Species)
    keep <- c(input$species_choice_1, "Not species specific")

    # Remove IDs that aren't actual species
    remove_non_species <- dplyr::setdiff(final_df$Species, keep)
    remove_non_species <- remove_non_species[remove_non_species != ""]

    if (length(remove_non_species) != 0) {
      # Fill these in as species unspecific
      final_df$Species <- stringi::stri_replace_all_regex(final_df$Species,
                                                          pattern=remove_non_species,
                                                          replacement="Not species specific",
                                                          vectorize=FALSE)
    }

    # Remove duplicates
    final_df <- final_df[!duplicated(final_df), ]

    # Order rows by workflow step
    new_order <- c("Punctuation suggestions", "Non-protein word suggestions", "GO exceptions", "Directly to PRO", "Specific protein suggestions", "CD synonym list", "Wikidata")
    final_df <- final_df[order(base::match(final_df$Match_step, new_order)),]

    # Get new column names
    colnames(final_df) <-
      c("Revised Name",
        "Original Name",
        "Sign",
        "Cluster",
        "PRO/GO Term",
        "PRO/GO Name",
        "Species Specific",
        "Match type",
        "Workflow Step",
        "Suggestion")

    # Change column order
    final_df <- final_df[, c("Original Name",	"Revised Name", "Sign", "Cluster", "Workflow Step",	"PRO/GO Name",	"PRO/GO Term",	"Species Specific",	"Match type",	"Suggestion")]

    return(final_df)
    } else {
      return(NULL)
    }
  })

  # Step 4, if some markers were not matched to PRO terms, prompt the user to rename them
  return_alt_markers_1 <- shiny::eventReactive(input$submit_tab1_step3, {
   shiny::req(input$reference_type_1 == 'internal_ref_1')
   
    if (is.data.frame(reformatted_data_1()) == TRUE) {
      all_ids <- PRO_results_1()

      # Just markers that did not match
      all_ids <- all_ids[all_ids$`PRO/GO Term` == "",]

      # Get certain columns
      all_ids <- all_ids[c("Original Name", "Suggestion")]

      if (nrow(all_ids) != 0) {

        # If no suggestion, link to PRO
        all_ids$Suggestion[all_ids$Suggestion == ""]<- "<a href='https://proconsortium.org/' target='_blank'>Search PRO</a>"

        # New column
        all_ids$'Type new marker name' <- ""

        # Remove duplicates
        all_ids <- all_ids[!duplicated(all_ids), ]

        if (nrow(all_ids) != 0) {
          return(all_ids)
        }  else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    }
  })

  # Step 4, get statement indicating if/which markers were not able to be matched
  return_unmatched_markers_1 <- reactive({
    shiny::req(input$reference_type_1 == 'internal_ref_1')
    shiny::req(input$submit_tab1_step3)

    all_ids <- PRO_results_1()

    # Get markers that did not match
    unmatched_markers <- all_ids[all_ids$`PRO/GO Term` == "",]

    # Get just the original inputted marker name
    unmatched_markers <- unmatched_markers$`Original Name`

    # Just unique markers
    unmatched_markers <- unique(unmatched_markers)

    # Return nothing if all markers had located IDs
    if (length(unmatched_markers) == 0) {
      return(NULL)
    } else {
      # Return markers that were not able to be matched
      return(unmatched_markers)
    }
  })

  # Step 4, output message about if markers matched to PRO term
  output$ui_warning_unmatched_markers_1 <- shiny::renderUI({
    shiny::req(input$reference_type_1 == 'internal_ref_1')
    shiny::req(input$submit_tab1_step3)

    if (length(marker_sign_error_1()) == 0) {
      # Return nothing if all markers had located IDs
      if (is.null(return_unmatched_markers_1())) {
        h4(
          paste(
            "All inputted markers matched to at least 1 protein ID. Please continue to the next step."
          )
        )
      } else if (length(return_unmatched_markers_1()) > 1) {
        marker_string <- paste(return_unmatched_markers_1(), collapse = ", ")
        h4(
          paste(
            "Protein IDs for the inputted markers",
            marker_string,
            "were not found. Please input new marker names."
          )
        )
      } else if (length(return_unmatched_markers_1()) == 1) {
        h4(
          paste(
            "A Protein ID for the inputted marker",
            return_unmatched_markers_1(),
            "was not found. Please input a new marker name."
          )
        )
      }
    }
  })

  # Step 4, option to edit (add new marker names) in the table
  shiny::observeEvent(input$table_new_format_1_cell_edit, {
    df_1$data <<-
      DT::editData(df_1$data,
                   input$table_new_format_1_cell_edit,
                   'table_new_format',
                   rownames = FALSE)
  })

  ######~#~# Step 5 - Get marker synonyms of markers with an edited name, user can choose which matched marker terms to procees with #~#~######

  # Step 5, side panel

  # Reset everything
  shiny::observeEvent(input$reset_tab1_step5, {
    session$reload()
  })

  # Step 5, help button for side panel parameters
  shiny::observeEvent(input$help_top_N_matches,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Return Top Matches per Cluster"),
      HTML("Use this slider to specify the number of top matches (1-10) to return per cluster. Matches will be ordered by the modified Jaccard similarity score. If multiple cell types have the same score, they will all be included. For example, if you request the top 3 matches but 5 cell types share the third highest score, the tool will output 7 matches for that cluster.
      <br>
      It is recommended to output at least 5 matches per cluster, as some cell types may be closely related."),
      easyClose = TRUE))
  })

  # Step 5, changes when the 'save and continue' button is hit if using internal references
  # Make editable dataframe, show buttons once that dataframe is made
  shiny::observeEvent(input$submit_tab1_step4, {
    shiny.destroy::removeOutput("table_new_format_1")
    shiny.destroy::removeOutput("help_tab1_step4")
    shiny.destroy::removeOutput("ui_warning_unmatched_markers_1")
    shiny.destroy::removeOutput("overview_tab1_step4")
    shinyjs::show("overview_tab1_step5")
    df2_1$data <- return_protein_synonym_1()
    if (is.data.frame(downloadable_step_1_long()) == TRUE) {
      # Step 5, make datatable of markers and their matched PRO terms
      output$table_protein_synonym_1 <- DT::renderDT({
        if (is.null(df2_1$data)) {
        } else {
          DT::datatable(
            df2_1$data,
            rownames = FALSE,
            escape = FALSE,
            options = list(pageLength = 25)
          ) }
      }, server=FALSE)

      if (is.data.frame(return_protein_synonym_1()) == TRUE) {
        shinyjs::show("help_tab1_step5")
        shinyjs::show("delete_PRO_synonym_1")
        shinyjs::show("download_PRO_GO_matches_1")
        shiny.destroy::removeOutput("all_sidebar_tab1_step4")
        shinyjs::show("all_sidebar_tab1_step5")
      }
    }
  })

  # Step 5, create editable dataframes
  df2_1 <- shiny::reactiveValues(data = NULL)

  # OG Step 5, look up PRO IDs for the revised marker names
  new_PRO_matches_1 <- reactive({
  shiny::req(input$submit_tab1_step4)

    ######## Species, references, initial dataframe ########

    new_df <- df_1$data

    # If empty, just fill with inputted value
    new_df$`Type new marker name`[new_df$`Type new marker name` == ""] = new_df$`Original Name`[new_df$`Type new marker name` == ""]

    # Get NCBI taxon ID
    species_ID <- input$species_choice_1
    matched_row <-
      which(species_in_PRO == species_ID, arr.ind = TRUE)[1]
    NCBI_taxon_ID <-
      species_in_PRO[matched_row, ncol(species_in_PRO)]

    NCBI_taxon_ID_full <- species_in_PRO[matched_row, 1]
    NCBI_taxon_ID_short <- gsub(".*/","",NCBI_taxon_ID_full)

    # Get Wikidata taxon ID
    SPARQL_query <-

      paste0(
        "SELECT ?item ?itemLabel WHERE {

     ?item wdt:P2888 <", NCBI_taxon_ID_full, ">.

    SERVICE wikibase:label { bd:serviceParam wikibase:language 'en'. }
    }"
      )

    # Run SPARQL on the endpoint
    result <- httr::GET(
      url = wiki_endpoint,
      query = list(query = SPARQL_query),
      httr::user_agent(R.version.string))

    # Will show a warning/error if there is any
    httr::stop_for_status(result)

    # Get result in text JSON
    x <- httr::content(result, as = "text") #, encoding = "UTF-8")

    # Convert from JSON to a list
    df <- jsonlite::fromJSON(x, flatten = TRUE)

    # Extract the data frame
    df <- df$results$bindings

    if (length(df) != 0) {

      # Remove unneeded info
      df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", ':lang')))

      # Remove this part that was added onto the column names
      colnames(df) <- gsub(".value", "", colnames(df))
    }

    wikidata_taxon_ID <- df$item

    # List of specific protein suggestions
    marker_alt <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/specific_suggestions_V1.csv", header = TRUE)

    # List of GO complexes that are included in CL definitions
    GO_complexes <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/GO_complexes_V1.csv", header = TRUE)

    ## List of CD synonym (adapted from list published by HCDM)
    CD_syn <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/Curated_synonyms_V1.csv", header = TRUE)

    # Make markers uppercase and remove dashes, underscores, periods, and spaces
    new_df$Marker2 <- toupper(new_df$`Type new marker name`)
    new_df$Marker2 <- gsub('-', '', new_df$Marker2, fixed = TRUE)
    new_df$Marker2 <- gsub('_', '', new_df$Marker2, fixed = TRUE)
    new_df$Marker2 <- gsub(' ', '', new_df$Marker2, fixed = TRUE)
    new_df$Marker2 <- gsub('.', '', new_df$Marker2, fixed = TRUE)

    # Just keep old marker name and new marker name
    new_df <- new_df[,c("Original Name", "Marker2")]

    # Make empty dataframe
    auto_matches <- data.frame(Marker = unique(new_df$Marker2))

    auto_matches$PRO_term <- ""
    auto_matches$PRO_name <- ""
    auto_matches$Species <- ""
    auto_matches$Match_type <- ""
    auto_matches$Match_step <- ""
    auto_matches$Suggestions <- ""

    ######## Automatic matching (CD3e, CD8a, TCR) ########

    # Automatically make any inputs of CD3e, CD8a, and TCR match to the wider complex

    # Auto replacements
    replace_markers <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/replace_markers_V1.csv", header=TRUE)

            # If TCR is also specified in the panel, don't auto match CD3E to it
    if (any(auto_matches=="TCRAB") | any(auto_matches=="TCRA/B" ) | 
        any(auto_matches=="TCRABCOMPLEX") | any(auto_matches=="TCRA/BCOMPLEX" ) |
        any(auto_matches=="ABTCR") | any(auto_matches=="ABTCRCOMPLEX" ) |
        any(auto_matches=="ALPHABETATCR") | any(auto_matches=="ALPHABETATCRCOMPLEX" ) |
        any(auto_matches=="TCRALPHABETA") | any(auto_matches=="TCRALPHABETACOMPLEX" ) |
        any(auto_matches=="TCRA") | any(auto_matches=="TCRB" ) |
        any(auto_matches=="TCRALPHA") | any(auto_matches=="TCRBETA" ) |
        any(auto_matches=="ALPHATCR") | any(auto_matches=="BETATCR" ) |
        any(auto_matches=="TCRGD") | any(auto_matches=="TCRGDCOMPLEX" ) |
        any(auto_matches=="TCRG/D") | any(auto_matches=="TCRG/DCOMPLEX" ) |
        any(auto_matches=="GDTCR") | any(auto_matches=="GDTCRCOMPLEX" ) |
        any(auto_matches=="GAMMADELTATCR") | any(auto_matches=="GAMMADELTATCRCOMPLEX" ) |
        any(auto_matches=="TCRGAMMADELTA") | any(auto_matches=="TCRGAMMADELTACOMPLEX" ) |
        any(auto_matches=="TCRG") | any(auto_matches=="TCRD" ) |
        any(auto_matches=="TCRGAMMA") | any(auto_matches=="TCRDELTA" ) |
        any(auto_matches=="GAMMATCR") | any(auto_matches=="DELTATCR" ) |
        any(auto_matches=="YTCR") | any(auto_matches=="TCRY" ) |
        any(auto_matches=="TCRYD") | any(auto_matches=="YDTCR" ) |
        any(auto_matches=="TCR") | any(auto_matches=="TCRCOMPLEX" ) |
        any(auto_matches=="TCELLRECEPTOR") | any(auto_matches=="TCELLRECEPTORCOMPLEX" )) {
      
      replace_markers <- replace_markers[!grepl("TCRAB", replace_markers$Replacement),]
      replace_markers <- replace_markers[!grepl("TCRGD", replace_markers$Replacement),]
      replace_markers <- replace_markers[!grepl("TCR", replace_markers$Replacement),]
    }
    
    # If CD8AB is also specified in the panel, don't auto match CD8 to it
    if (any(auto_matches=="CD8AB") | any(auto_matches=="CD8ALPHABETA")) {
      replace_markers <- replace_markers[!grepl("CD8ALPHABETA", replace_markers$Replacement),]
    }
   
    auto_matches <- merge(replace_markers, auto_matches,  by = "Marker", all.y = TRUE)

    auto_replaced <- auto_matches[!is.na(auto_matches$Replacement),]

    auto_matches$Replacement[is.na(auto_matches$Replacement)] <- auto_matches$Marker[is.na(auto_matches$Replacement)]

    auto_matches$Marker <- NULL

    colnames(auto_matches)[1] <- "Marker"

    ######## GO Complexes ########

    # Match entries to the GO complexes listed in the CL

    # Vector of unique marker names
    GO_matches <- auto_matches

    unique_marker_names <- GO_matches$Marker

    all_ids <- list()

    for(marker2 in unique_marker_names) {

      matched_row <- GO_complexes[GO_complexes$Inputted_Name ==marker2,]

      if (nrow(matched_row) != 0) {

        GO_matches[GO_matches$Marker == marker2,]$Suggestions <- ""
        GO_matches[GO_matches$Marker == marker2,]$Match_step <- "GO exceptions"
        GO_matches[GO_matches$Marker == marker2,]$PRO_term <- matched_row$Matched_GO_term
        GO_matches[GO_matches$Marker == marker2,]$PRO_name <- matched_row$Matched_Descriptive_GO_name
      }
    }

    ######## Direct SPAPRQL query to PRO ########

    # Directly query PRO

    Match_step <- "Directly to PRO"

    unique_marker_names <- GO_matches$Marker[GO_matches$Match_step != "GO exceptions"]

    # Run SPARQL
    all_ids <- PRO_SPARQL(unique_marker_names = unique_marker_names,
                          NCBI_taxon_ID = NCBI_taxon_ID_short,
                          onto_endpoint = onto_endpoint,
                          Match_step = Match_step)

    # Put all results in one data frame
    Direct_PRO_matches <- dplyr::bind_rows(all_ids, .id = "Marker")

    ######## CD synonym list, then SPAPRQL query to PRO ########

    # Search CD synonym list, then if there are synonyms query them in PRO

    Match_step <- "CD synonym list"

    # Vector of unique marker names
    unique_marker_names <- Direct_PRO_matches$Marker[Direct_PRO_matches$Match_step != "Directly to PRO"] # Specific_protein_sugg_matches$Marker[Specific_protein_sugg_matches$Match_step == ""]

    all_ids_final <- list()

    if(length(unique_marker_names) != 0) {

      for (marker in unique_marker_names) {

        matched_row <- CD_syn[which(CD_syn == marker, arr.ind=TRUE)[,1],]

        matched_row <- data.frame(x=unlist(matched_row))

        unique_syn_marker_names <- unique(stats::na.omit(matched_row))

        unique_syn_marker_names <- as.vector(unique_syn_marker_names[,1])

        # Run SPARQL
        all_ids <- PRO_SPARQL(unique_marker_names = unique_syn_marker_names,
                              NCBI_taxon_ID = NCBI_taxon_ID_short,
                              onto_endpoint = onto_endpoint,
                              Match_step = Match_step)

        # Put all results in one data frame
        all_ids <- dplyr::bind_rows(all_ids, .id = "Matched_synonym")

        all_ids <- stats::na.omit(all_ids)

        if (length(all_ids) != 0) {
          all_ids_final[[marker]] <- all_ids
        } else if (length(all_ids) == 0) {
          all_ids_final[[marker]] <- data.frame(Matched_synonym = "", PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
        }
      }

      # Put all results in one data frame
      all_ids_final2 <- dplyr::bind_rows(all_ids_final, .id = "Marker")

      # Remove duplicated entries (Marker, PRO term, and type of match are the same)
      all_ids_final2 <- all_ids_final2[!duplicated(all_ids_final2[,c("Marker","PRO_term", "Match_type")]),]

      # Get entries were inputted marker and matched species is the same

      # Subset into entries with and without duplicates
      x <-  all_ids_final2[duplicated(all_ids_final2[,c("Marker","Species")]) | duplicated(all_ids_final2[,c("Marker","Species")], fromLast = TRUE) ,]
      y <-  dplyr::setdiff(all_ids_final2, x)

      # Get vector of markers (duplicates)
      unique_marker_names <- unique(x$Marker)

      # If there is both an exact match and a secondary match, remove the secondary
      remove_secondary_per_marker <- data.frame()

      for(marker in unique_marker_names) {
        df_subsetted <- x[x$Marker == marker,]
        if(any(df_subsetted$Match_type == "Primary")) {
          df_per_marker <- df_subsetted[df_subsetted$Match_type != "Secondary",]
          df_per_marker <- df_per_marker[df_per_marker$Match_type != "",]
          remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
        } else if(any(df_subsetted$Match_type == "Secondary")) {
          df_per_marker <- df_subsetted[df_subsetted$Match_type != "",]
          remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
        }
      }

      # Final
      CD_syn_list_matches <- rbind(remove_secondary_per_marker, y)

      # Just columns we want
      CD_syn_list_matches <- unique(CD_syn_list_matches[ , c("Marker", "PRO_term", "PRO_name", "Species", "Match_type", "Match_step")])

    } else {
      CD_syn_list_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
    }

    ######## Wikidata, then SPAPRQL query to PRO ########

    # Query Wikidata, then if there are synonyms query them in PRO

    Match_step <- "Wikidata"

    # Vector of unique marker names
    unique_marker_names <- CD_syn_list_matches$Marker[CD_syn_list_matches$Match_step != "CD synonym list"]

    # Remove duplicates introduces by synonyms from last step
    remove <- CD_syn_list_matches$Marker[CD_syn_list_matches$Match_step != ""]
    unique_marker_names <-  dplyr::setdiff(unique_marker_names, remove)

    all_ids_final <- list()

    if(length(unique_marker_names) != 0) {
      if(unique_marker_names[1] != "") {

        # Loop through each and do the SPARQL query
        for(marker in unique_marker_names) {

          specific_marker <-  shQuote(marker)

          SPARQL_query <-
            paste0(
              "
SELECT DISTINCT ?item ?itemLabel ?altLabel WHERE {
  # Look either for non-species specific or add a specific species
  OPTIONAL {
    ?item wdt:P703 <", wikidata_taxon_ID, "> .

  }
  # Get Wikidata protein name and alternative names
  ?item rdfs:label ?itemLabel .
  ?item skos:altLabel ?altLabel.
  # Return results that are either proteins or protein-coding genes
  {?item wdt:P31 wd:Q8054}
  UNION
  {?item wdt:P279 wd:Q20747295}
  # Search for the marker input, either if it is the entry name or an alternative name
  FILTER (UCASE(REPLACE(str(?itemLabel),'[ -.]','')) = ", specific_marker, " || UCASE(REPLACE(str(?altLabel),'[ -.]','')) = ", specific_marker, ")
}
")
          # Run SPARQL on the endpoint
          result <- httr::GET(
            url = wiki_endpoint,
            query = list(query = SPARQL_query),
            httr::user_agent(R.version.string))

          # Get result in text JSON
          x <- httr::content(result, as = "text") #, encoding = "UTF-8")

          # Convert from JSON to a list
          df <- jsonlite::fromJSON(x, flatten = TRUE)

          # Extract the data frame
          df <- df$results$bindings

          if (length(df) != 0) {

            # Remove unneeded info
            df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", 'lang')))

            # Remove this part that was added onto the column names
            colnames(df) <- gsub(".value", "", colnames(df))

            # Make everything case, remove spaces, dashes
            df <- dplyr::mutate_all(df, .funs = toupper)
            df <- as.data.frame(lapply(df, function(y) gsub('-', '', y)))
            df <- as.data.frame(lapply(df, function(y) gsub('_', '', y)))
            df <- as.data.frame(lapply(df, function(y) gsub(' ', '', y)))
            df <- as.data.frame(lapply(df, function(y) gsub('\\.', '', y)))

            df <- df[(which(nchar(df$altLabel) > 2)),]
            df <- df[(which(nchar(df$itemLabel) > 2)),]

            # No marker can start with a number, so add X to any marker names that start with a number
            df <- as.data.frame(lapply(df, function(y) gsub('^(\\d)', 'X\\1', y)))

            # Vector of unique marker names
            unique_syn_names <- unique(c(df[, "altLabel"], df[, "itemLabel"]))

            unique_syn_names <- unique_syn_names[!grepl("[^\x01-\x7F]+", unique_syn_names)]

            unique_syn_names <- unique_syn_names[grep('[A-Za-z]',unique_syn_names)]

            unique_syn_names <- unique_syn_names[!grepl(paste0('^', marker, '$'), unique_syn_names, useBytes = TRUE)]

            # Run SPARQL
            all_ids <- PRO_SPARQL(unique_marker_names = unique_syn_names,
                                  NCBI_taxon_ID = NCBI_taxon_ID_short,
                                  onto_endpoint = onto_endpoint,
                                  Match_step = Match_step)

            # Put all results in one data frame
            all_ids <- dplyr::bind_rows(all_ids, .id = "Matched synonym")

            all_ids <- stats::na.omit(all_ids)

            if (length(all_ids) != 0) {
              all_ids_final[[marker]] <- all_ids
            } else if (length(all_ids) == 0) {
              all_ids_final[[marker]] <- data.frame(PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
            }
          } else if (length(df) == 0) {
            all_ids_final[[marker]] <- data.frame(PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
          }
        }

      # Put all results in one data frame
      all_ids_final2 <- dplyr::bind_rows(all_ids_final, .id = "Marker")

      # Remove duplicated entries (Marker, PRO term, and type of match are the same)
      all_ids_final2 <- all_ids_final2[!duplicated(all_ids_final2[,c("Marker","PRO_term", "Match_type")]),]

      # Get entries were inputted marker and matched species is the same

      # Subset into entries with and without duplicates
      x <-  all_ids_final2[duplicated(all_ids_final2[,c("Marker","Species")]) | duplicated(all_ids_final2[,c("Marker","Species")], fromLast = TRUE) ,]
      y <-  dplyr::setdiff(all_ids_final2, x)

      # Get vector of markers (duplicates)
      unique_marker_names <- unique(x$Marker)

      # If there is both an exact match and a secondary match, remove the secondary
      remove_secondary_per_marker <- data.frame()

      for(marker in unique_marker_names) {
        df_subsetted <- x[x$Marker == marker,]
        if(any(df_subsetted$Match_type == "Primary")) {
          df_per_marker <- df_subsetted[df_subsetted$Match_type != "Secondary",]
          df_per_marker <- df_per_marker[df_per_marker$Match_type != "",]
          remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
        } else if(any(df_subsetted$Match_type == "Secondary")) {
          df_per_marker <- df_subsetted[df_subsetted$Match_type != "",]
          remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
        }
      }

      # Final
      Wikidata_matches <- rbind(remove_secondary_per_marker, y)

      # Just columns we want
      Wikidata_matches <- unique(Wikidata_matches[ , c("Marker", "PRO_term", "PRO_name", "Species", "Match_type", "Match_step")])

      } else {
        Wikidata_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
      }
    } else {
      Wikidata_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
    }

    ######## Final results ########

    # Put everything together and get the final list of matched/unmatched marker names

    # Add empty suggestion column to parts where suggestions were not relevant
    Direct_PRO_matches$Suggestions <- ""
    CD_syn_list_matches$Suggestions <- ""
    Wikidata_matches$Suggestions <- ""

    # Subset entries that were matched at each step
    GO_matches2 <- GO_matches[GO_matches$Match_step == "GO exceptions",]
    Direct_PRO_matches2 <- Direct_PRO_matches[Direct_PRO_matches$Match_step == "Directly to PRO",]
    CD_syn_list_matches2 <- CD_syn_list_matches[CD_syn_list_matches$Match_step == "CD synonym list",]
    Wikidata_matches2 <- Wikidata_matches[Wikidata_matches$Match_step == "Wikidata",]

    # Put everything together into 1 dataframe
    all_ids <- do.call("rbind", list(GO_matches2,  Direct_PRO_matches2, CD_syn_list_matches2, Wikidata_matches2))

    # Put original inputted names back in
    final_df <- merge(new_df, all_ids, by.x = "Marker2", by.y = "Marker", all.x = TRUE, all.y = TRUE)

    # Change names
    colnames(final_df)[1:2] <- c("New_Name", "Original_Name")

    if (nrow(auto_replaced) != 0) {
      # Put back original auto changed names
      replacements <- auto_replaced[as.character(auto_replaced$Marker) != as.character(auto_replaced$Replacement),]

      replacements <- replacements[,1:2]

      colnames(replacements) <- c("Marker2", "New_Name")

      replacements <- merge(replacements, new_df, by = "Marker2")

      replacements <- replacements[,c(2, 3)]

      replacements <- merge(replacements, final_df, by = "New_Name")

      replacements <- replacements[,-3]

      colnames(replacements)[2] <- "Original_Name"

      final_df <- rbind(final_df, replacements)
    }

    final_df <- final_df[!is.na(final_df$Original_Name), ]

    # If there was no match or suggestion, indicate that
    final_df$Match_step[is.na(final_df$Match_step)] <- "Did not match"

    # Turn any remaining NA entries into blanks
    final_df[is.na(final_df)] <- ""

    # Change species URIs into human readable form
    final_df$Species <- gsub(NCBI_taxon_ID_full, input$species_choice_1, final_df$Species)
    keep <- c(input$species_choice_1, "Not species specific")

    # Remove IDs that aren't actual species
    remove_non_species <- dplyr::setdiff(final_df$Species, keep)
    remove_non_species <- remove_non_species[remove_non_species != ""]

    if (length(remove_non_species) != 0) {
      # Fill these in as species unspecific
      final_df$Species <- stringi::stri_replace_all_regex(final_df$Species,
                                                          pattern=remove_non_species,
                                                          replacement="Not species specific",
                                                          vectorize=FALSE)
    }

    new_df <- reformatted_data_1()

    # Remove duplicates
    final_df <- final_df[!duplicated(final_df), ]

    # Get new column names
    colnames(final_df) <-
      c("Revised Name",
        "Original Name",
        "PRO/GO Term",
        "PRO/GO Name",
        "Species Specific",
        "Match type",
        "Workflow Step",
        "Suggestion")

    return(final_df)

  })


  # Step 5, see if inputed marker is in reference
  markers_in_ref <- shiny::eventReactive(input$submit_tab1_step4, {

    if (length(input$ontology_type_1) == 1) {
      if (input$ontology_type_1 == 'internal_CL_1') {
        ontology_URI <- "<http://purl.obolibrary.org/obo/merged/CL>"
      } else if (input$ontology_type_1 == 'internal_pCL_1') {
        ontology_URI <- "<http://purl.obolibrary.org/obo/merged/PCL>"
      }
    } else if (length(input$ontology_type_1) == 2) {
      ontology_URI <- "<http://purl.obolibrary.org/obo/merged/CL> FROM <http://purl.obolibrary.org/obo/merged/PCL>"
    }

    relationship_type1 <- "<http://purl.obolibrary.org/obo/RO_0002104>"
    relationship_type2 <- "<http://purl.obolibrary.org/obo/RO_0015015>"
    relationship_type3 <- "<http://purl.obolibrary.org/obo/RO_0015016>"
    relationship_type4 <- "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
    relationship_type5 <- "<http://purl.obolibrary.org/obo/CL_4030046>"
    relationship_type6 <- "<http://purl.obolibrary.org/obo/BFO_0000051>"

    # SPARQL query
    query <-

      paste0(
        "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?Related_protein ?Protein_name ?Relationship_type

  FROM ",
        ontology_URI,
        "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }

  FILTER (?Relationship_type =", relationship_type1, "|| ?Relationship_type =", relationship_type2, "|| ?Relationship_type =", relationship_type3, "|| ?Relationship_type =", relationship_type4, "|| ?Relationship_type =", relationship_type5, "|| ?Relationship_type =", relationship_type6, ")

    }
  "
      )

    # Run SPARQL on the endpoint
    result <- httr::POST(onto_endpoint,
                         body = list(query = query),
                         httr::user_agent(R.version.string))

    # Will show a warning/error if there is any
    httr::stop_for_status(result)

    # Get result in text JSON
    x <- httr::content(result, "text", encoding = "UTF-8")

    # Convert from JSON to a list
    df <- jsonlite::fromJSON(x, flatten = TRUE)

    # Extract the dataframe
    df <- df$results$bindings

    # Put the results in a list of data frames
    if (length(df) != 0) {

      # Remove unneeded info
      df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

      # Remove this part that was added onto the column names
      colnames(df) <- gsub(".value", "", colnames(df))
    }

    # Only get relevant relationship types for PRO/GO terms
    just_PRO <- subset(df, startsWith(as.character(Related_protein),"http://purl.obolibrary.org/obo/PR_"))
    just_GO <- subset(df, startsWith(as.character(Related_protein),"http://purl.obolibrary.org/obo/GO_") & Relationship_type != "http://purl.obolibrary.org/obo/BFO_0000051")

    # Combine
    df <- rbind(just_PRO, just_GO)

    # Only keep first 2 columns
    df <- df[ , c("Related_protein","Protein_name")]

    # Remove duplicates
    df <- df[!duplicated(df), ]

    return(df)
  })

  # Step 5, match inputted proteins to entries in the Protein Ontology, get PRO IDs
  #  Make PRO link clickable
  #  Merge in info about whether or not marker is in the reference(s)
  return_protein_synonym_1 <- shiny::eventReactive(input$submit_tab1_step4, {
    shiny::req(input$reference_type_1 == 'internal_ref_1')

    # Original results
    org_pro_results <- PRO_results_1()

    # Remove no matches found from original PRO query
    org_pro_results <- org_pro_results[org_pro_results$`PRO/GO Term` != "",]

    # Get certain columns
    org_pro_results <- org_pro_results[c("Original Name", "Revised Name", "PRO/GO Term")]

    # Remove duplicates
    org_pro_results <- org_pro_results[!duplicated(org_pro_results), ]

    if (!is.null(df_1$data)) {
      # Redo results
      new_pro_results <- new_PRO_matches_1()

      # Get certain columns
      new_pro_results <- new_pro_results[c("Original Name", "Revised Name", "PRO/GO Term")]

      # Remove duplicates
      new_pro_results <- new_pro_results[!duplicated(new_pro_results), ]

      # Bind all PRO results together
      all_pro_results <- rbind(org_pro_results, new_pro_results)

    } else {
      all_pro_results <- org_pro_results
    }

    # If there was no match or suggestion, indicate that
    all_pro_results$`PRO/GO Term`[all_pro_results$`PRO/GO Term` == ""] <- "Did not match"

    df_ref <- markers_in_ref()

    # Get inputted markers that are not in the reference
    not_in_reference <- setdiff(all_pro_results$`PRO/GO Term`, df_ref$Related_protein)

    if (length(not_in_reference) != 0) {
      df <- data.frame("PRO/GO Term" = not_in_reference, "Reference" = "Not in reference. This PRO term will be removed from subsequent analysis.", check.names = FALSE)

      all_pro_results <-  merge(all_pro_results,
                                df, by = "PRO/GO Term", all = TRUE)

      all_pro_results$Reference[is.na(all_pro_results$Reference)] <- "Included in reference."

    } else {
      all_pro_results$Reference <- "Included in reference."

    }

    # Change column order
    all_pro_results <- all_pro_results[, c("Original Name",	"Revised Name",	"PRO/GO Term", "Reference")]

    # Order by marker name
    all_pro_results <-
      all_pro_results[order(all_pro_results$`Original Name`, decreasing = FALSE), ]

    # Just get PRO ID to make link in R Shiny look better
    just_IDs <-
      gsub("http://purl.obolibrary.org/obo/",
           "",
           all_pro_results$`PRO/GO Term`)

    # Remove '<' and '>' from URI
    all_pro_results$`PRO/GO Term` <-
      gsub('<', '', all_pro_results$`PRO/GO Term`)
    all_pro_results$`PRO/GO Term` <-
      gsub('>', '', all_pro_results$`PRO/GO Term`)
    just_IDs <- gsub('>', '', just_IDs)
    just_IDs <- gsub('<', '', just_IDs)

    # Make it so the link is clickable
    for (a in 1:nrow(all_pro_results)) {
      if (all_pro_results$`PRO/GO Term`[a] != "Did not match") {
        all_pro_results$`PRO/GO Term`[a] <-
          paste0(
            "<a href='",
            all_pro_results$`PRO/GO Term`[a],
            "' target='_blank'>",
            just_IDs[a],
            "</a>"
          )
      }
    }
    return(all_pro_results)
  })

  # Step 5, remove link from PRO/GO term for download, 1st dataset
  downloadable_step_1 <- reactive({ # eventReactive(input$submit_tab1_step4, {
  shiny::req(input$submit_tab1_step4)
   
    remove_link <- df2_1$data

    if (all(remove_link$`PRO/GO Term` == "Did not match")) {
      return(remove_link)
    } else {
      # Remove syntax that made this linkable in R Shiny
      remove_link$`PRO/GO Term` <-
        gsub("' target='_blank'>[^>]+</a>",
             "",
             remove_link$`PRO/GO Term`)
      remove_link$`PRO/GO Term` <-
        gsub("<a href='", "", remove_link$`PRO/GO Term`)

      # Remove any NA entries
      remove_link <- remove_link[stats::complete.cases(remove_link), ]

      return(remove_link)
    }
  })

  # Step 5, 2nd dataset, all info
  downloadable_step_1_long <- reactive({
  shiny::req(input$submit_tab1_step4)
   
    markers <- downloadable_step_1()

    signs_clusters <- reformatted_data_1()

    all <- merge(signs_clusters, markers, by.x = "Marker", by.y = "Original Name")

    # Order by cluster
    all <- all[order(all$Cluster, decreasing = FALSE), ]

    return(all)

  })

  # Step 5, when help button is clicked show message
  shiny::observeEvent(input$help_tab1_step5,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 5<br>Confirm Matched Marker Terms"),
      HTML("This step displays the results of the marker standardization process. The data frame shows:
      <br>
      <br>
      <li>Original marker names</li>
      <br>
      <li>Revised marker names (if applicable)</li>
      <br>
      <li>Matched Protein Ontology (PRO) or Gene Ontology (GO) terms, which are clickable for direct access to the ontology page</li>
      <br>
      <li>Reference, which indicates whether the term is included in the Cell Ontology</li>
      <br>
      If markers match to multiple ontology terms, the tool will list all possible matches. This can occur if an inputted name is unspecific and listed as a synonym for multiple different proteins. Additionally, an inputted marker may match to both a species specific and species non-specific identifier. Any ontology term that begins with a 'P' (PR_P???) is species specific. Oftentimes specific-specific terms are not included within the Cell Ontology, but when they are it is suggested that they remain included.
      <br>
      <br>
      Markers not found in the Cell Ontology are excluded from subsequent steps. Additionally, any markers can be manually removed by selecting the row(s) and clicking 'Delete Marker Entry'. Deleting markers is irreversible.
      <br>
      <br>
      Users can download the matched terms in CSV format using the 'Download' button.
      <br>
      <br>"),
      tags$div("Additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/3.1.-Input-Expression-Data-and-Default-References', "here.", target="_blank")),
      easyClose = TRUE))
  })

  # Step 5, option to delete PRO entries before doing a cell type matching search
  shiny::observeEvent(input$delete_PRO_synonym_1, {
    if (!is.null(input$table_protein_synonym_1_rows_selected)) {
      df2_1$data <<-
        df2_1$data[-as.numeric(input$table_protein_synonym_1_rows_selected), ]
    }
  })

  # Step 5, download PRO/GO matches
  output$download_PRO_GO_matches_1 <- shiny::downloadHandler(

    # Filename
    filename = function() {
      paste0(gsub(".csv|.xls|.xlsx", "", input$upload_expression_csv),
             "-matched-PRO-GO-terms-", Sys.Date(), ".zip")
    },

    content = function(file) {
      # definition of content to download
      to_dl <- list(
        # names to use in file names
        names = list(a = "Individual_PRO_GO_terms",
                     b = "All_PRO_GO_terms"),
        # data
        data = list(a = downloadable_step_1(),
                    b = downloadable_step_1_long())
      )

      # temp dir for the csv's as we can only create
      # an archive from existent files and not data from R
      twd <- setwd(tempdir())
      on.exit(setwd(twd))
      files <- NULL

      # loop on data to download and write individual csv's
      for (i in c("a", "b")) {
        fileName <- paste0(to_dl[["names"]][[i]], ".csv") # csv file name
        utils::write.csv(to_dl[["data"]][[i]], fileName, row.names = FALSE) # write csv in temp dir

        files <- c(files, fileName) # store written file name
      }
      # create zip file with all the csv files
      zip::zip(file, files)
    }
  )

  ######~#~# Step 6 - Match Median Difference Equation results and PRO/GO terms with cell types, show user ranked results and CL definitions #~#~######

  # Step 6, side panel

  # Reset everything
  shiny::observeEvent(input$reset_tab1_step6, {
    session$reload()
  })

  # Step 6, main panel

  # Step 6, changes when the 'save and continue' button is hit if using internal references
  # Output final dataframe, show buttons once that dataframe is made
  shiny::observeEvent(input$submit_tab1_step5, {
    shiny.destroy::removeOutput("help_tab1_step5")
    shiny.destroy::removeOutput("overview_tab1_step5")
    shiny.destroy::removeOutput("table_protein_synonym_1")
    shiny.destroy::removeOutput("delete_PRO_synonym_1")
    shiny.destroy::removeOutput("download_PRO_GO_matches_1")
    shinyjs::show("overview_tab1_step6")

    if (!is.null(return_cell_types_linked_1())) {
      # Step 6, make datatable of matched cell types
      output$table_cell_types_1 <- DT::renderDT({
        DT::datatable(
          df_subsettable_1(),
          rownames = FALSE,
          escape = FALSE,
          options = list(pageLength = 25)
        )
      }, server=FALSE)

      # Step 6/Uploaded reference, download final results
      output$download_table_cell_types_1 <- shiny::downloadHandler(
        # File name
        filename = function() {
          paste0(gsub(".csv|.xls|.xlsx", "", input$upload_expression_csv),
                 "-matched-cell-types-", Sys.Date(), ".csv")
        },

        # Write to csv
        content = function(file) {
          utils::write.csv(return_cell_types_1(),  file, row.names = FALSE)
        }
      )
      shiny.destroy::removeOutput("all_sidebar_tab1_step5")
      shinyjs::show("all_sidebar_tab1_step6")
      shinyjs::show("help_tab1_step6")
      shinyjs::show("table_cell_types_1")
      shinyjs::show("download_table_cell_types_1")
    } else {
      shinyjs::show("no_matched_cell_types_1")
    }

  })

  # Step 6, Match markers (PRO terms) to cell types in CL/pCL
  # Get wikidata link for each cell type
  return_cell_types_1 <- shiny::eventReactive(input$submit_tab1_step5, {

    low_option_1 <- "low_positive_1"
    high_option_1 <- "high_positive_1"

    # Load in data from Step 6
    input_reformatted <- reformatted_data_1()

    # Load in data from Step 4
    all_PRO_matches <- df2_1$data

    all_PRO_matches <- all_PRO_matches[all_PRO_matches$Reference == "Included in reference.",]

    all_PRO_matches <- all_PRO_matches [ , !(names(all_PRO_matches ) %in% "Reference")]

    all_id_results <- input_reformatted

    colnames(all_id_results)[which(names(all_id_results) == "Marker")] <-
      "Original Name"

    all_id_results <-
      merge(all_id_results,
            all_PRO_matches,
            by.x = "Original Name",
            by.y = "Original Name")

    if (all(all_id_results$`PRO/GO Term` == "No match found")) {
      return(NULL)
    } else {

      # Remove syntax that made this linkable in R Shiny
      all_id_results$`PRO/GO Term` <-
        gsub("' target='_blank'>[^>]+</a>",
             "",
             all_id_results$`PRO/GO Term`)
      all_id_results$`PRO/GO Term` <-
        gsub("<a href='", "", all_id_results$`PRO/GO Term`)

      # Add back in < and >
      all_id_results$`PRO/GO Term` <-
        paste0("<", all_id_results$`PRO/GO Term`, ">")

      all_id_results2 <- all_id_results

      ## Get contradictions

      # Add empty columns to be filled
      all_id_results[, 'Relations_ontology_1'] <- NA
      all_id_results[, 'Relations_ontology_1'] <- NA
      all_id_results[, 'Relations_ontology_3'] <- NA
      all_id_results[, 'Relations_ontology_4'] <- NA
      all_id_results[, 'Relations_ontology_5'] <- NA

      #Loop through data, make positives, negatives, high, lows into ontology relations terminology
      for (a in 1:nrow(all_id_results)) {

        # Positive signs are always contradicted by negatives
        if (all_id_results$Sign[a] == "Positive") {
          all_id_results$Relations_ontology_1[a] <-
            "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
          all_id_results$Relations_ontology_2[a] <-
            "<http://purl.obolibrary.org/obo/CL_4030046>"

          all_id_results$Marker_query[a] <-
            paste(
              "?Related_protein =",
              all_id_results$`PRO/GO Term`[a],
              "&& ?Relationship_type = ",
              all_id_results$Relations_ontology_1[a],
              " || ?Related_protein = ",
              all_id_results$`PRO/GO Term`[a],
              "&& ?Relationship_type = ",
              all_id_results$Relations_ontology_2[a]
            )
        }

        # If choosing low/negative, contradictions of negatives are any positive except low
        #  if (input$low_option_1 == "low_negative_1") {
        if (low_option_1 == "low_negative_1") {

          if (all_id_results$Sign[a] == "Negative") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a]
              )
            # If not low/negative, the contradictions of positive is any positive, high, or low sign
          }} #else  if (input$low_option_1 != "low_negative_1") {

        else  if (low_option_1 != "low_negative_1") {

          if (all_id_results$Sign[a] == "Negative") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part
            all_id_results$Relations_ontology_4[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_4[a]
              )
          }}

        # If choosing high/positive, then any high signs are contradicted by negatives and lows
        #    if (input$high_option_1 == "high_positive_1") {
        if (high_option_1 == "high_positive_1") {
          if (all_id_results$Sign[a] == "High") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>" # lacks plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a]
              )
          }}

        # If choosing high/high, then any high signs are contradicted by negatives, positives, and lows
        #   if (input$high_option_1 == "high_high_1") {
        if (high_option_1 == "high_high_1") {
          if (all_id_results$Sign[a] == "High") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>" # lacks plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part
            all_id_results$Relations_ontology_4[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results$Relations_ontology_5[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_4[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_5[a]
              )
          }}


        # If choosing low/positive, then any low signs are contradicted by negatives and highs
        #  if (input$low_option_1 == "low_positive_1") {
        if (low_option_1 == "low_positive_1") {
          if (all_id_results$Sign[a] == "Low") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>" # lacks plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a]
              )
          }}

        # If choosing low/low, then any low signs are contradicted by negatives, positives, and high
        # if (input$low_option_1 == "low_low_1") {
        if (low_option_1 == "low_low_1") {
          if (all_id_results$Sign[a] == "Low") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>" # lacks plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
            all_id_results$Relations_ontology_4[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results$Relations_ontology_5[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_4[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_5[a]
              )
          }}


        # If choosing low/negative, then any low signs are contradicted by positives and high
        #  if (input$low_option_1 == "low_negative_1") {
        if (low_option_1 == "low_negative_1") {
          if (all_id_results$Sign[a] == "Low") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a]
              )
          }}
      }

      # Make sure cluster name/number is read in as string
      all_id_results$Cluster <- as.character(all_id_results$Cluster)

      # Empty list to put results in
      marker_in_cluster_matches <- list()

      # Look through each group/cluster category
      for (a_cluster in sort(unique(all_id_results$Cluster))) {
        # Just the specific grouping
        specific_cluster <-
          all_id_results[all_id_results$Cluster == a_cluster, ]

        # Remove proteins that didn't get a match
        specific_cluster <-
          specific_cluster[specific_cluster$`PRO/GO Term` != "No match found", ]

        for (b in sort(unique(specific_cluster$`Revised Name`))) {
          # Just the specific marker
          specific_marker <-
            specific_cluster[specific_cluster$`Revised Name` == b, ]

          # Filter for SPARQL query
          full_filter <-
            paste(specific_marker$Marker_query, collapse = "||")


          if (length(input$ontology_type_1) == 1) {
            if (input$ontology_type_1 == 'internal_CL_1') {
              ontology_URI <- "<http://purl.obolibrary.org/obo/merged/CL>"
              full_filter_1 <-
                paste(
                  "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/CL_\"  ) &&",
                  full_filter
                )
            } else if (input$ontology_type_1 == 'internal_pCL_1') {
              ontology_URI <- "<http://purl.obolibrary.org/obo/merged/PCL>"
              full_filter_1 <-
                paste(
                  "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/PCL_\"  ) &&",
                  full_filter
                )
            }

            # SPARQL query
            query1 <-

              paste0(
                "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                ontology_URI,
                "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }



  FILTER(",
                full_filter_1,
                ")
    }
  "
              )
            # Run SPARQL on the endpoint
            result1 <- httr::POST(onto_endpoint,
                                  body = list(query = query1),
                                  httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result1)

            # Get result in text JSON
            x1 <- httr::content(result1, "text", encoding = "UTF-8")

            # Convert from JSON to a list
            df <- jsonlite::fromJSON(x1, flatten = TRUE)

            # Extract the dataframe
            df <- df$results$bindings

            # Put the results in a list of data frames
            if (length(df) != 0) {

              # Remove unneeded info
              # df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype")))
              df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

              # Remove this part that was added onto the column names
              colnames(df) <- gsub(".value", "", colnames(df))

              # Make sure they all have these columns
              six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

              df[six_cols[!(six_cols %in% colnames(df))]] <- NA
            }
          } else if (length(input$ontology_type_1) == 2) {
            ontology_URI_1 <- "<http://purl.obolibrary.org/obo/merged/CL>"
            full_filter_1 <-
              paste(
                "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/CL_\"  ) &&",
                full_filter
              )

            # SPARQL query
            query1 <-

              paste0(
                "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                ontology_URI_1,
                "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }



  FILTER(",
                full_filter_1,
                ")
    }
  "
              )
            # Run SPARQL on the endpoint
            result1 <- httr::POST(onto_endpoint,
                                  body = list(query = query1),
                                  httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result1)

            # Get result in text JSON
            x1 <- httr::content(result1, "text", encoding = "UTF-8")

            # Convert from JSON to a list
            df1 <- jsonlite::fromJSON(x1, flatten = TRUE)

            # Extract the dataframe
            df1 <- df1$results$bindings

            # Put the results in a list of data frames
            if (length(df1) != 0) {

              # Remove unneeded info
              df1 <- df1 %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

              # Remove this part that was added onto the column names
              colnames(df1) <- gsub(".value", "", colnames(df1))

              # Make sure they all have these columns
              six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

              df1[six_cols[!(six_cols %in% colnames(df1))]] <- NA
            }

            ontology_URI_2 <-
              "<http://purl.obolibrary.org/obo/merged/PCL>"
            full_filter_2 <-
              paste(
                "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/PCL_\"  ) &&",
                full_filter
              )

            # SPARQL query
            query2 <-

              paste0(
                "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                ontology_URI_2,
                "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }

  FILTER(",
                full_filter_2,
                ")
    }
  "
              )
            # Run SPARQL on the endpoint
            result2 <- httr::POST(onto_endpoint,
                                  body = list(query = query2),
                                  httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result2)

            # Get result in text JSON
            x2 <- httr::content(result2, "text", encoding = "UTF-8")

            # Convert from JSON to a list
            df2 <- jsonlite::fromJSON(x2, flatten = TRUE)

            # Extract the dataframe
            df2 <- df2$results$bindings

            # Put the results in a list of data frames
            if (length(df2) != 0) {

              # Remove unneeded info
              df2 <- df2 %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

              # Remove this part that was added onto the column names
              colnames(df2) <- gsub(".value", "", colnames(df2))

              # Make sure they all have these columns
              six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

              df2[six_cols[!(six_cols %in% colnames(df2))]] <- NA
            }

            # Add results
            df <- rbind(df1, df2)
          }

          # Put the results in a list of data frames
          if (length(df) != 0) {
            # Add back in < and >
            df$CL_term <- paste0("<", df$CL_term, ">")
            for (c in 1:nrow(df)) {
              if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0002104") {
                df$Relationship_type[c] <- "Positive"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/BFO_0000051") {
                df$Relationship_type[c] <- "Positive"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015015") {
                df$Relationship_type[c] <- "High"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015016") {
                df$Relationship_type[c] <- "Low"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part") {
                df$Relationship_type[c] <- "Negative"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/CL_4030046") {
                df$Relationship_type[c] <- "Negative"
              }
            }

            if(!"Cell_comment" %in% colnames(df)) {
              df$Cell_comment = NA
            }

            # Collapse to get ensure each marker is only counted 1 time per cell cluster match
            # Some entries may appear more then once bc there is more then 1 marker relation for that protein or bc a protein has more then one PRO ID
            df <- df %>%
              dplyr::group_by(CL_term, Cell_type, Cell_comment) %>%
              dplyr::summarise(
                Relationship_type = paste(Relationship_type, collapse = ", "),
                Related_protein = paste(Related_protein, collapse = ", "),
                Protein_name = paste(Protein_name, collapse = ", ")
              )

            # Remove repeats that the last step may have formed
            df$Relationship_type <-
              sapply(df$Relationship_type, function(x)
                paste(unique(unlist(
                  strsplit(x, ", ")
                )), collapse = ", "))
            df$Related_protein <-
              sapply(df$Related_protein, function(x)
                paste(unique(unlist(
                  strsplit(x, ", ")
                )), collapse = ", "))
            df$Protein_name <-
              sapply(df$Protein_name, function(x)
                paste(unique(unlist(
                  strsplit(x, ", ")
                )), collapse = ", "))

            # Include cluster name/number
            df <- cbind(Cluster = a_cluster, df)

          } else if (length(df) == 0) {
            df <- data.frame(CL_term = "No match found")
          }

          # Name
          cluster_marker_name <- paste0(a_cluster, "::", b)

          # Put into list
          marker_in_cluster_matches[[cluster_marker_name]] <-
            as.data.frame(df)
        }
      }

      # Put all results in one data frame
      all_cluster_matches <-
        dplyr::bind_rows(marker_in_cluster_matches, .id = "Marker")

      # Remove proteins that didn't get a match in a particular cluster
      all_cluster_matches <-
        all_cluster_matches[all_cluster_matches$CL_term != "No match found", ]

      # Remove cluster from unique name
      all_cluster_matches$Marker <-
        gsub(".*::", "", all_cluster_matches$Marker)

      if(nrow(all_cluster_matches) >= 1) {

        # New column
        all_cluster_matches$`Matched inputted marker` <-
          paste0(all_cluster_matches$Marker,
                 " (",
                 all_cluster_matches$Relationship_type,
                 ")")

        all_cluster_matches$Cluster <- as.numeric(all_cluster_matches$Cluster)

        # Sort, put into final format for export
        final_matched_output <- all_cluster_matches %>%
          dplyr::group_by(Cluster,
                          CL_term,
                          Cell_type,
                          Cell_comment) %>%
          dplyr::summarise(
            `Matched inputted marker` = paste(`Matched inputted marker`, collapse = ", "),
            `num of inputs that match to cell type` =  dplyr::n()) %>%
          dplyr::arrange(Cluster, dplyr::desc(`num of inputs that match to cell type`))

        # Remove '<' and '>' from URI
        final_matched_output$CL_term <-
          gsub('<', '', final_matched_output$CL_term)
        final_matched_output$CL_term <-
          gsub('>', '', final_matched_output$CL_term)

        final_matched_output2 <- final_matched_output

        # Final column name order
        final_matched_output2 <- final_matched_output2[, c(
          "CL_term",
          "Cluster",
          "Matched inputted marker",
          "num of inputs that match to cell type"

        )]

        colnames(final_matched_output2)[3] <- "Contradiction(s)"

        colnames(final_matched_output2)[4] <- "num contradiction(s)"

        ## Get matches

        # Add empty columns to be filled
        all_id_results2[, 'Relations_ontology_1'] <- NA
        all_id_results2[, 'Relations_ontology_2'] <- NA
        all_id_results2[, 'Relations_ontology_3'] <- NA
        all_id_results2[, 'Relations_ontology_4'] <- NA

        #Loop through data, make positives, negatives, high, lows into ontology relations terminology
        for (a in 1:nrow(all_id_results2)) {
          if (all_id_results2$Sign[a] == "Positive") {
            all_id_results2$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results2$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
            all_id_results2$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part
            all_id_results2$Relations_ontology_4[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results2$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_3[a],
                " || ?Related_protein = ",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_4[a]
              )
          }

          if (all_id_results2$Sign[a] == "Negative") {
            all_id_results2$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results2$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>"

            all_id_results2$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_2[a]
              )
          }

          # if (input$high_option_1 == "high_positive_1") {
          if (high_option_1 == "high_positive_1") {
            if (all_id_results2$Sign[a] == "High") {
              all_id_results2$Relations_ontology_1[a] <-
                "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
              all_id_results2$Relations_ontology_2[a] <-
                "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part
              all_id_results2$Relations_ontology_3[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_2[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_3[a]
                )
            }}

          #  if (input$high_option_1 == "high_high_1") {
          if (high_option_1 == "high_high_1") {
            if (all_id_results2$Sign[a] == "High") {
              all_id_results2$Relations_ontology_1[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a]
                )
            }}

          # if (input$low_option_1 == "low_positive_1") {
          if (low_option_1 == "low_positive_1") {
            if (all_id_results2$Sign[a] == "Low") {
              all_id_results2$Relations_ontology_1[a] <-
                "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
              all_id_results2$Relations_ontology_2[a] <-
                "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part
              all_id_results2$Relations_ontology_3[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_2[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_3[a]
                )
            }}


          #if (input$low_option_1 == "low_low_1") {
          if (low_option_1 == "low_low_1") {
            if (all_id_results2$Sign[a] == "Low") {
              all_id_results2$Relations_ontology_1[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a]
                )
            }}

          #if (input$low_option_1 == "low_negative_1") {
          if (low_option_1 == "low_negative_1") {
            if (all_id_results2$Sign[a] == "Low") {
              all_id_results2$Relations_ontology_1[a] <-
                "http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part"
              all_id_results2$Relations_ontology_2[a] <-
                "<http://purl.obolibrary.org/obo/CL_4030046>"
              all_id_results2$Relations_ontology_3[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_2[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_3[a]
                )
            }}
        }

        # Make sure cluster name/number is read in as string
        all_id_results2$Cluster <- as.character(all_id_results2$Cluster)

        # Empty list to put results in
        marker_in_cluster_matches <- list()

        # Look through each group/cluster category
        for (a_cluster in sort(unique(all_id_results2$Cluster))) {
          # Just the specific grouping
          specific_cluster <-
            all_id_results2[all_id_results2$Cluster == a_cluster, ]

          # Remove proteins that didn't get a match
          specific_cluster <-
            specific_cluster[specific_cluster$`PRO/GO Term` != "No match found", ]

          for (b in sort(unique(specific_cluster$`Revised Name`))) {
            # Just the specific marker
            specific_marker <-
              specific_cluster[specific_cluster$`Revised Name` == b, ]

            # Filter for SPARQL query
            full_filter <-
              paste(specific_marker$Marker_query, collapse = "||")


            if (length(input$ontology_type_1) == 1) {
              if (input$ontology_type_1 == 'internal_CL_1') {
                ontology_URI <- "<http://purl.obolibrary.org/obo/merged/CL>"
                full_filter_1 <-
                  paste(
                    "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/CL_\"  ) &&",
                    full_filter
                  )
              } else if (input$ontology_type_1 == 'internal_pCL_1') {
                ontology_URI <- "<http://purl.obolibrary.org/obo/merged/PCL>"
                full_filter_1 <-
                  paste(
                    "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/PCL_\"  ) &&",
                    full_filter
                  )
              }

              # SPARQL query
              query1 <-

                paste0(
                  "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                  ontology_URI,
                  "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }

  FILTER(",
                  full_filter_1,
                  ")
    }
  "
                )
              # Run SPARQL on the endpoint
              result1 <- httr::POST(onto_endpoint,
                                    body = list(query = query1),
                                    httr::user_agent(R.version.string))

              # Will show a warning/error if there is any
              httr::stop_for_status(result1)

              # Get result in text JSON
              x1 <- httr::content(result1, "text", encoding = "UTF-8")

              # Convert from JSON to a list
              df <- jsonlite::fromJSON(x1, flatten = TRUE)

              # Extract the dataframe
              df <- df$results$bindings

              # Put the results in a list of data frames
              if (length(df) != 0) {

                # Remove unneeded info
                # df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype")))
                df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

                # Remove this part that was added onto the column names
                colnames(df) <- gsub(".value", "", colnames(df))

                # Make sure they all have these columns
                six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

                df[six_cols[!(six_cols %in% colnames(df))]] <- NA
              }
            } else if (length(input$ontology_type_1) == 2) {
              ontology_URI_1 <- "<http://purl.obolibrary.org/obo/merged/CL>"
              full_filter_1 <-
                paste(
                  "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/CL_\"  ) &&",
                  full_filter
                )

              # SPARQL query
              query1 <-

                paste0(
                  "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                  ontology_URI_1,
                  "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }

  FILTER(",
                  full_filter_1,
                  ")
    }
  "
                )
              # Run SPARQL on the endpoint
              result1 <- httr::POST(onto_endpoint,
                                    body = list(query = query1),
                                    httr::user_agent(R.version.string))

              # Will show a warning/error if there is any
              httr::stop_for_status(result1)

              # Get result in text JSON
              x1 <- httr::content(result1, "text", encoding = "UTF-8")

              # Convert from JSON to a list
              df1 <- jsonlite::fromJSON(x1, flatten = TRUE)

              # Extract the dataframe
              df1 <- df1$results$bindings

              # Put the results in a list of data frames
              if (length(df1) != 0) {

                # Remove unneeded info
                df1 <- df1 %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

                # Remove this part that was added onto the column names
                colnames(df1) <- gsub(".value", "", colnames(df1))

                # Make sure they all have these columns
                six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

                df1[six_cols[!(six_cols %in% colnames(df1))]] <- NA
              }

              ontology_URI_2 <-
                "<http://purl.obolibrary.org/obo/merged/PCL>"
              full_filter_2 <-
                paste(
                  "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/PCL_\"  ) &&",
                  full_filter
                )

              # SPARQL query
              query2 <-

                paste0(
                  "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                  ontology_URI_2,
                  "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }

  FILTER(",
                  full_filter_2,
                  ")
    }
  "
                )
              # Run SPARQL on the endpoint
              result2 <- httr::POST(onto_endpoint,
                                    body = list(query = query2),
                                    httr::user_agent(R.version.string))

              # Will show a warning/error if there is any
              httr::stop_for_status(result2)

              # Get result in text JSON
              x2 <- httr::content(result2, "text", encoding = "UTF-8")

              # Convert from JSON to a list
              df2 <- jsonlite::fromJSON(x2, flatten = TRUE)

              # Extract the dataframe
              df2 <- df2$results$bindings

              # Put the results in a list of data frames
              if (length(df2) != 0) {

                # Remove unneeded info
                df2 <- df2 %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

                # Remove this part that was added onto the column names
                colnames(df2) <- gsub(".value", "", colnames(df2))

                # Make sure they all have these columns
                six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

                df2[six_cols[!(six_cols %in% colnames(df2))]] <- NA
              }

              # Add results
              df <- rbind(df1, df2)
            }

            # Put the results in a list of data frames
            if (length(df) != 0) {
              # Add back in < and >
              df$CL_term <- paste0("<", df$CL_term, ">")
              for (c in 1:nrow(df)) {
                if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0002104") {
                  df$Relationship_type[c] <- "Positive"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/BFO_0000051") {
                  df$Relationship_type[c] <- "Positive"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015015") {
                  df$Relationship_type[c] <- "High"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015016") {
                  df$Relationship_type[c] <- "Low"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part") {
                  df$Relationship_type[c] <- "Negative"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/CL_4030046") {
                  df$Relationship_type[c] <- "Negative"
                }
              }

              if(!"Cell_comment" %in% colnames(df)) {
                df$Cell_comment = NA
              }

              # Collapse to get ensure each marker is only counted 1 time per cell cluster match
              # Some entries may appear more then once bc there is more then 1 marker relation for that protein or bc a protein has more then one PRO ID
              df <- df %>%
                dplyr::group_by(CL_term, Cell_type, Cell_comment) %>%
                dplyr::summarise(
                  Relationship_type = paste(Relationship_type, collapse = ", "),
                  Related_protein = paste(Related_protein, collapse = ", "),
                  Protein_name = paste(Protein_name, collapse = ", ")
                )

              # Remove repeats that the last step may have formed
              df$Relationship_type <-
                sapply(df$Relationship_type, function(x)
                  paste(unique(unlist(
                    strsplit(x, ", ")
                  )), collapse = ", "))
              df$Related_protein <-
                sapply(df$Related_protein, function(x)
                  paste(unique(unlist(
                    strsplit(x, ", ")
                  )), collapse = ", "))
              df$Protein_name <-
                sapply(df$Protein_name, function(x)
                  paste(unique(unlist(
                    strsplit(x, ", ")
                  )), collapse = ", "))

              # Include cluster name/number
              df <- cbind(Cluster = a_cluster, df)

            } else if (length(df) == 0) {
              df <- data.frame(CL_term = "No match found")
            }

            # Name
            cluster_marker_name <- paste0(a_cluster, "::", b)

            # Put into list
            marker_in_cluster_matches[[cluster_marker_name]] <-
              as.data.frame(df)
          }
        }

        # Put all results in one data frame
        all_cluster_matches <-
          dplyr::bind_rows(marker_in_cluster_matches, .id = "Marker")

        # Remove proteins that didn't get a match in a particular cluster
        all_cluster_matches <-
          all_cluster_matches[all_cluster_matches$CL_term != "No match found", ]

        # Remove cluster from unique name
        all_cluster_matches$Marker <-
          gsub(".*::", "", all_cluster_matches$Marker)

        if(nrow(all_cluster_matches) >= 1) {

          # New column
          all_cluster_matches$`Matched inputted marker` <-
            paste0(all_cluster_matches$Marker,
                   " (",
                   all_cluster_matches$Relationship_type,
                   ")")

          # Sort, put into final format for export
          final_matched_output <- all_cluster_matches %>%
            dplyr::group_by(Cluster,
                            CL_term,
                            Cell_type,
                            Cell_comment) %>%
            dplyr::summarise(
              `Matched inputted marker` = paste(`Matched inputted marker`, collapse = ", "),
              `num of inputs that match to cell type` =  dplyr::n()
            ) %>%
            dplyr::arrange(Cluster, dplyr::desc(`num of inputs that match to cell type`))


          unq_cell_types <- final_matched_output[!duplicated(final_matched_output[,"CL_term"]),]

          # Put all cell types in another SPARQL query to get full set of proteins
          unq_cell_types$Cell_type_query <-
            paste(
              "strstarts(str(?Related_protein), \"http://purl.obolibrary.org/obo/PR_\") && str(?CL_term) =",
              unq_cell_types$CL_term, "|| strstarts(str(?Related_protein), \"http://purl.obolibrary.org/obo/GO_\") && str(?CL_term) =",
              unq_cell_types$CL_term
            )

          # Split up into smaller chunks to the SPARQL query can handle the filter
          # Number of cell types per query = 100
          chunk <- 70 #100
          total_num_rows <- nrow(unq_cell_types)
          how_to_split  <- rep(1:ceiling(total_num_rows/chunk),each=chunk)[1:total_num_rows]
          final_matched_output_list <- split(unq_cell_types,how_to_split)

          # Empty dataframe to append to
          df_all_pro <- data.frame()

          for(i in 1:length(final_matched_output_list)){

            each_part <- final_matched_output_list[[i]]

            full_filter <-
              paste(each_part$Cell_type_query, collapse = " || ")

            query <-

              paste0(
                "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Relationship_type ?Protein_name

  FROM <http://purl.obolibrary.org/obo/merged/CL>

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  ?Class_of owl:someValuesFrom ?Related_protein .

  ?Related_protein rdfs:label ?Protein_name .

  ?Class_of owl:onProperty ?Relationship_type .

  FILTER(",
                full_filter,
                ")

    }
  "
              )

            # Run SPARQL on the endpoint
            result <- httr::POST(onto_endpoint,
                                 body = list(query = query),
                                 httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result)

            # Get result in text JSON
            x <- httr::content(result, "text", encoding = "UTF-8")

            # Convert from JSON to a list
            df_all_pro_part <- jsonlite::fromJSON(x, flatten = TRUE)

            # Extract the dataframe
            df_all_pro_part <- df_all_pro_part$results$bindings

            # Remove unneeded info
            df_all_pro_part <- df_all_pro_part %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

            # Remove this part that was added onto the column names
            colnames(df_all_pro_part) <- gsub(".value", "", colnames(df_all_pro_part))

            # Make sure they all have these columns
            three_cols <- c("CL_term", "Relationship_type", "Protein_name")
            df_all_pro_part[three_cols[!(three_cols %in% colnames(df_all_pro_part))]] <- NA


            df_all_pro <- rbind(df_all_pro, df_all_pro_part)
          }

          # Change relations into human readable
          for (c in 1:nrow(df_all_pro)) {
            if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0002104") {
              df_all_pro$Relationship_type[c] <- "Positive"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/BFO_0000051") {
              df_all_pro$Relationship_type[c] <- "Positive"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015015") {
              df_all_pro$Relationship_type[c] <- "High"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015016") {
              df_all_pro$Relationship_type[c] <- "Low"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part") {
              df_all_pro$Relationship_type[c] <- "Negative"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/CL_4030046") {
              df_all_pro$Relationship_type[c] <- "Negative"
            } else {
              df_all_pro$Relationship_type[c] <- ""
            }
          }

          # Remove proteins that didn't get a match in a particular cluster
          df_all_pro <- df_all_pro[df_all_pro$Relationship_type != "", ]

          # Remove any duplicate rows
          df_all_pro <- df_all_pro[!duplicated(df_all_pro),]


          # If protein as listed as having 2 different relationship types
          df_all_pro <- df_all_pro %>%
            dplyr::group_by(CL_term, Protein_name) %>%
            dplyr::summarise(
              Relationship_type = paste(Relationship_type, collapse = ", ")
            )



          # Remove protein identifier
          df_all_pro <- df_all_pro[!df_all_pro$Protein_name == "protein", ]


          df_all_pro$x <-
            paste0(df_all_pro$Protein_name,
                   " (",
                   df_all_pro$Relationship_type,
                   ")")

          # If a relationship type is negative and positive remove entirely
          # Shouldn't be happening unless mistake in CL
          df_all_pro <- df_all_pro[!grepl("Negative, Positive", df_all_pro$Relationship_type),]
          df_all_pro <- df_all_pro[!grepl("Positive, Negative", df_all_pro$Relationship_type),]


          # Sort, put into final format for export
          df_all_pro <- df_all_pro %>%
            dplyr::group_by(CL_term) %>%
            dplyr::summarise(
              `Full marker description` = paste(x, collapse = ", "),
              `total number of markers` = dplyr::n()
            )

          # Remove '<' and '>' from URI
          final_matched_output$CL_term <-
            gsub('<', '', final_matched_output$CL_term)
          final_matched_output$CL_term <-
            gsub('>', '', final_matched_output$CL_term)

          final_matched_output <-
            merge(final_matched_output,
                  df_all_pro,
                  by.x = "CL_term",
                  by.y = "CL_term",
                  all.x = TRUE)

          # Get wikidata link
          for (i in 1:nrow(final_matched_output)) {

            CL_term <- final_matched_output[i,"CL_term"]

            # Make CL term to be in the proper format
            CL_term <- gsub(".*/","",CL_term)
            CL_term <- shQuote(CL_term)

            SPARQL_query <-

              paste0(
                "SELECT ?item ?itemLabel WHERE {
  ?item wdt:P7963 ", CL_term, ".
  ?item wdt:P31 wd:Q189118

  SERVICE wikibase:label { bd:serviceParam wikibase:language 'en'. }
}"
              )

            # Run SPARQL on the endpoint
            result <- httr::GET(
              url = wiki_endpoint,
              query = list(query = SPARQL_query),
              httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result)

            # Get result in text JSON
            x <- httr::content(result, as = "text") #, encoding = "UTF-8")

            # Convert from JSON to a list
            df <- jsonlite::fromJSON(x, flatten = TRUE)

            # Extract the data frame
            df <- df$results$bindings

            # Put the results in a list of data frames
            if (length(df) != 0) {

              # Remove unneeded info
              df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

              # Remove this part that was added onto the column names
              colnames(df) <- gsub(".value", "", colnames(df))

              if (length(df) == 1) {
                final_matched_output[i, "Wikidata"] <- df$item
              } else if (length(df) > 1) {
                final_matched_output[i, "Wikidata"] <- paste0(df$item, collapse = ", ")
              }
            } else if (length(df) == 0) {
              final_matched_output[i, "Wikidata"] <- ""
            }
          }

          # All input markers
          input_markers <- reformatted_data_1()

          input_markers$`All inputted markers` <-
            paste0(input_markers$Marker,
                   " (",
                   input_markers$Sign,
                   ")")

          input_markers <- input_markers %>%
            dplyr::group_by(Cluster) %>%
            dplyr::summarise(
              `Full input marker description` = paste(`All inputted markers`, collapse = ", "),
              `input_num_total` = dplyr::n()
            )

          final_matched_output <-
            merge(final_matched_output,
                  input_markers,
                  by.x = "Cluster",
                  by.y = "Cluster")

          final_matched_output3 <-
            merge(x = final_matched_output,
                  y = final_matched_output2,
                  by = c("CL_term","Cluster"), all.x = TRUE)



          final_matched_output3$`num contradiction(s)`[is.na(final_matched_output3$`num contradiction(s)`)] <- 0

          # Count how many markers are positive and negative
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_pos_input = stringr::str_count(`Full input marker description`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_neg_input = stringr::str_count(`Full input marker description`, '(Negative)'))

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_pos_cell_type = stringr::str_count(`Full marker description`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_neg_cell_type = stringr::str_count(`Full marker description`, '(Negative)'))

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_pos_match = stringr::str_count(`Matched inputted marker`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_neg_match = stringr::str_count(`Matched inputted marker`, '(Negative)'))

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_pos_contr = stringr::str_count(`Contradiction(s)`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_neg_contr = stringr::str_count(`Contradiction(s)`, 'Negative'))

          # TCR

          # Positive
          TCR_pos <- c('TCR \\(POSITIVE\\)', 'TCRAB \\(POSITIVE\\)', 'TCRA \\(POSITIVE\\)', 'TCRB \\(POSITIVE\\)', 'TCRGD \\(POSITIVE\\)', 'TCRG \\(POSITIVE\\)', 'TCRD \\(POSITIVE\\)', 'CD3E \\(POSITIVE\\)',
                       'TCR \\(HIGH\\)', 'TCRAB \\(HIGH\\)', 'TCRA \\(HIGH\\)', 'TCRB \\(HIGH\\)', 'TCRGD \\(HIGH\\)', 'TCRG \\(HIGH\\)', 'TCRD \\(HIGH\\)', 'CD3E \\(HIGH\\)',
                       'TCR \\(LOW\\)', 'TCRAB \\(LOW\\)', 'TCRA \\(LOW\\)', 'TCRB \\(LOW\\)', 'TCRGD \\(LOW\\)', 'TCRG \\(LOW\\)', 'TCRD \\(LOW\\)', 'CD3E \\(LOW\\)',
                       'TCR \\(LOW, POSITIVE\\)', 'TCRAB \\(LOW, POSITIVE\\)', 'TCRA \\(LOW, POSITIVE\\)', 'TCRB \\(LOW, POSITIVE\\)', 'TCRGD \\(LOW, POSITIVE\\)', 'TCRG \\(LOW, POSITIVE\\)', 'TCRD \\(LOW, POSITIVE\\)', 'CD3E \\(LOW, POSITIVE\\)',
                       'TCR \\(POSITIVE, LOW\\)', 'TCRAB \\(POSITIVE, LOW\\)', 'TCRA \\(POSITIVE, LOW\\)', 'TCRB \\(POSITIVE, LOW\\)', 'TCRGD \\(POSITIVE, LOW\\)', 'TCRG \\(POSITIVE, LOW\\)', 'TCRD \\(POSITIVE, LOW\\)', 'CD3E \\(POSITIVE, LOW\\)',
                       'TCR \\(LOW, HIGH\\)', 'TCRAB \\(LOW, HIGH\\)', 'TCRA \\(LOW, HIGH\\)', 'TCRB \\(LOW, HIGH\\)', 'TCRGD \\(LOW, HIGH\\)', 'TCRG \\(LOW, HIGH\\)', 'TCRD \\(LOW, HIGH\\)', 'CD3E \\(LOW, HIGH\\)',
                       'TCR \\(HIGH, POSITIVE\\)', 'TCRAB \\(HIGH, POSITIVE\\)', 'TCRA \\(HIGH, POSITIVE\\)', 'TCRB \\(HIGH, POSITIVE\\)', 'TCRGD \\(HIGH, POSITIVE\\)', 'TCRG \\(HIGH, POSITIVE\\)', 'TCRD \\(HIGH, POSITIVE\\)', 'CD3E \\(HIGH, POSITIVE\\)',
                       'TCR \\(POSITIVE, HIGH\\)', 'TCRAB \\(POSITIVE, HIGH\\)', 'TCRA \\(POSITIVE, HIGH\\)', 'TCRB \\(POSITIVE, HIGH\\)', 'TCRGD \\(POSITIVE, HIGH\\)', 'TCRG \\(POSITIVE, HIGH\\)', 'TCRD \\(POSITIVE, HIGH\\)', 'CD3E \\(POSITIVE, HIGH\\)',
                       'TCR \\(HIGH, LOW\\)', 'TCRAB \\(HIGH, LOW\\)', 'TCRA \\(HIGH, LOW\\)', 'TCRB \\(HIGH, LOW\\)', 'TCRGD \\(HIGH, LOW\\)', 'TCRG \\(HIGH, LOW\\)', 'TCRD \\(HIGH, LOW\\)', 'CD3E \\(HIGH, LOW\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_input_pos = stringr::str_count(toupper(`Full marker description`), paste(TCR_pos, collapse='|')))

          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 1] <- 0
          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 2] <- 1
          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 3] <- 2
          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 4] <- 3
          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_match_pos = stringr::str_count(toupper(`Matched inputted marker`), paste(TCR_pos, collapse='|')))

          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 1] <- 0
          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 2] <- 1
          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 3] <- 2
          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 4] <- 3
          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_contradictions_pos = stringr::str_count(toupper(`Contradiction(s)`), paste(TCR_pos, collapse='|')))

          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 1] <- 0
          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 2] <- 1
          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 3] <- 2
          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 4] <- 3
          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 5] <- 4
          final_matched_output3$TCR_contradictions_pos[is.na(final_matched_output3$TCR_contradictions_pos)] <- 0

          TCR_pos <- c('T cell receptor complex \\(Positive\\)', 'alpha-beta T cell receptor complex \\(Positive\\)', 'gamma-delta T cell receptor complex \\(Positive\\)', 'CD3 epsilon \\(Positive\\)',
                       'T cell receptor complex \\(High\\)', 'alpha-beta T cell receptor complex \\(High\\)', 'gamma-delta T cell receptor complex \\(High\\)', 'CD3 epsilon \\(High\\)',
                       'T cell receptor complex \\(Low\\)', 'alpha-beta T cell receptor complex \\(Low\\)', 'gamma-delta T cell receptor complex \\(Low\\)', 'CD3 epsilon \\(Low\\)',
                       'T cell receptor complex \\(Low, Positive\\)', 'alpha-beta T cell receptor complex \\(Low, Positive\\)', 'gamma-delta T cell receptor complex \\(Low, Positive\\)', 'CD3 epsilon \\(Low, Positive\\)',
                       'T cell receptor complex \\(Positive, Low\\)', 'alpha-beta T cell receptor complex \\(Positive, Low\\)', 'gamma-delta T cell receptor complex \\(Positive, Low\\)', 'CD3 epsilon \\(Positive, Low\\)',
                       'T cell receptor complex \\(Low, High\\)', 'alpha-beta T cell receptor complex \\(Low, High\\)', 'gamma-delta T cell receptor complex \\(Low, High\\)', 'CD3 epsilon \\(Low, High\\)',
                       'T cell receptor complex \\(High, Positive\\)', 'alpha-beta T cell receptor complex \\(High, Positive\\)', 'gamma-delta T cell receptor complex \\(High, Positive\\)', 'CD3 epsilon \\(High, Positive\\)',
                       'T cell receptor complex \\(Positive, High\\)', 'alpha-beta T cell receptor complex \\(Positive, High\\)', 'gamma-delta T cell receptor complex \\(Positive, High\\)', 'CD3 epsilon \\(Positive, High\\)',
                       'T cell receptor complex \\(High, Low\\)', 'alpha-beta T cell receptor complex \\(High, Low\\)', 'gamma-delta T cell receptor complex \\(High, Low\\)', 'CD3 epsilon \\(High, Low\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_cell_type_pos = stringr::str_count(`Full marker description`, paste(TCR_pos, collapse='|')))

          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 1] <- 0
          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 2] <- 1
          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 3] <- 2
          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 4] <- 3
          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 5] <- 4

          # Negative

          TCR_negs <- c('TCR \\(NEGATIVE\\)', 'TCRAB \\(NEGATIVE\\)', 'TCRA \\(NEGATIVE\\)', 'TCRB \\(NEGATIVE\\)', 'TCRGD \\(NEGATIVE\\)', 'TCRG \\(NEGATIVE\\)', 'TCRD \\(NEGATIVE\\)', 'CD3E \\(NEGATIVE\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_input_neg = stringr::str_count(toupper(`Full marker description`), paste(TCR_negs, collapse='|')))

          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg_pos == 1] <- 0
          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg == 2] <- 1
          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg == 3] <- 2
          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg == 4] <- 3
          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_match_neg = stringr::str_count(toupper(`Matched inputted marker`), paste(TCR_negs, collapse='|')))

          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 1] <- 0
          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 2] <- 1
          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 3] <- 2
          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 4] <- 3
          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_contradictions_neg = stringr::str_count(toupper(`Contradiction(s)`), paste(TCR_negs, collapse='|')))

          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 1] <- 0
          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 2] <- 1
          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 3] <- 2
          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 4] <- 3
          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 5] <- 4
          final_matched_output3$TCR_contradictions_neg[is.na(final_matched_output3$TCR_contradictions_neg)] <- 0

          TCR_negs <- c('T cell receptor complex \\(Negative\\)', 'alpha-beta T cell receptor complex \\(Negative\\)', 'gamma-delta T cell receptor complex \\(Negative\\)', 'CD3 epsilon \\(Negative\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_cell_type_neg = stringr::str_count(`Full marker description`, paste(TCR_negs, collapse='|')))

          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 1] <- 0
          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 2] <- 1
          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 3] <- 2
          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 4] <- 3
          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 5] <- 4

          final_matched_output3$`total number of markers` <- as.numeric(final_matched_output3$`total number of markers`)
          final_matched_output3$input_num_total <- as.numeric(final_matched_output3$input_num_total)
          final_matched_output3$`num of inputs that match to cell type` <- as.numeric(final_matched_output3$`num of inputs that match to cell type`)
          final_matched_output3$`num contradiction(s)` <- as.numeric(final_matched_output3$`num contradiction(s)`)

          final_matched_output3$`total number of markers` <- final_matched_output3$`total number of markers` - final_matched_output3$TCR_cell_type_pos - final_matched_output3$TCR_cell_type_neg

          final_matched_output3$input_num_total <- final_matched_output3$input_num_total - final_matched_output3$TCR_input_pos - final_matched_output3$TCR_input_neg

          final_matched_output3$`num of inputs that match to cell type` <- final_matched_output3$`num of inputs that match to cell type` - final_matched_output3$TCR_match_neg - final_matched_output3$TCR_match_pos

          final_matched_output3$`num contradiction(s)` <- final_matched_output3$`num contradiction(s)` - final_matched_output3$TCR_contradictions_neg - final_matched_output3$TCR_contradictions_pos

          final_matched_output3$Num_pos_input <- final_matched_output3$Num_pos_input - final_matched_output3$TCR_input_pos
          final_matched_output3$Num_neg_input <- final_matched_output3$Num_neg_input - final_matched_output3$TCR_input_neg

          final_matched_output3$Num_pos_cell_type <- final_matched_output3$Num_pos_cell_type - final_matched_output3$TCR_cell_type_pos
          final_matched_output3$Num_neg_cell_type <- final_matched_output3$Num_neg_cell_type -final_matched_output3$TCR_cell_type_neg

          final_matched_output3$Num_pos_match <- final_matched_output3$Num_pos_match - final_matched_output3$TCR_match_pos
          final_matched_output3$Num_neg_match <- final_matched_output3$Num_neg_match - final_matched_output3$TCR_match_neg

          final_matched_output3$Num_pos_contr <- final_matched_output3$Num_pos_contr - final_matched_output3$TCR_contradictions_pos
          final_matched_output3$Num_neg_contr <- final_matched_output3$Num_neg_contr -final_matched_output3$TCR_contradictions_neg

          # CD8
          # Positive
          CD8_pos <- c('CD8A \\(POSITIVE\\)', 'CD8ALPHABETA \\(POSITIVE\\)',
                       'CD8A \\(HIGH\\)', 'CD8ALPHABETA \\(HIGH\\)',
                       'CD8A \\(LOW\\)', 'CD8ALPHABETA \\(LOW\\)',
                       'CD8A \\(LOW, POSITIVE\\)', 'CD8ALPHABETA \\(LOW, POSITIVE\\)',
                       'CD8A \\(POSITIVE, LOW\\)', 'CD8ALPHABETA \\(POSITIVE, LOW\\)',
                       'CD8A \\(LOW, HIGH\\)', 'CD8ALPHABETA \\(LOW, HIGH\\)',
                       'CD8A \\(HIGH, POSITIVE\\)', 'CD8ALPHABETA \\(HIGH, POSITIVE\\)',
                       'CD8A \\(POSITIVE, HIGH\\)', 'CD8ALPHABETA \\(POSITIVE, HIGH\\)',
                       'CD8A \\(HIGH, LOW\\)', 'CD8ALPHABETA \\(HIGH, LOW\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_input_pos = stringr::str_count(toupper(`Full marker description`), paste(CD8_pos, collapse='|')))

          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 1] <- 0
          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 2] <- 1
          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 3] <- 2
          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 4] <- 3
          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_match_pos = stringr::str_count(toupper(`Matched inputted marker`), paste(CD8_pos, collapse='|')))

          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 1] <- 0
          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 2] <- 1
          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 3] <- 2
          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 4] <- 3
          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_contradictions_pos = stringr::str_count(toupper(`Contradiction(s)`), paste(CD8_pos, collapse='|')))

          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 1] <- 0
          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 2] <- 1
          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 3] <- 2
          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 4] <- 3
          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 5] <- 4
          final_matched_output3$CD8_contradictions_pos[is.na(final_matched_output3$CD8_contradictions_pos)] <- 0

          CD8_pos <- c('T-cell surface glycoprotein CD8 alpha chain \\(Positive\\)', 'T cell receptor co-receptor CD8 \\(Positive\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(High\\)', 'T cell receptor co-receptor CD8 \\(High\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Low\\)', 'T cell receptor co-receptor CD8 \\(Low\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Low, Positive\\)', 'T cell receptor co-receptor CD8 \\(Low, Positive\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Positive, Low\\)', 'T cell receptor co-receptor CD8 \\(Positive, Low\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Low, High\\)', 'T cell receptor co-receptor CD8 \\(Low, High\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(High, Positive\\)', 'T cell receptor co-receptor CD8 \\(High, Positive\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Positive, High\\)', 'T cell receptor co-receptor CD8 \\(Positive, High\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(High, Low\\)', 'T cell receptor co-receptor CD8 \\(High, Low\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_cell_type_pos = stringr::str_count(`Full marker description`, paste(CD8_pos, collapse='|')))
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 1] <- 0
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 2] <- 1
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 3] <- 2
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 4] <- 3
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 5] <- 4

          # Negative

          CD8_negs <- c('CD8A \\(NEGATIVE\\)', 'CD8ALPHABETA \\(NEGATIVE\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_input_neg = stringr::str_count(toupper(`Full marker description`), paste(CD8_negs, collapse='|')))

          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg_pos == 1] <- 0
          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg == 2] <- 1
          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg == 3] <- 2
          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg == 4] <- 3
          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_match_neg = stringr::str_count(toupper(`Matched inputted marker`), paste(CD8_negs, collapse='|')))

          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 1] <- 0
          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 2] <- 1
          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 3] <- 2
          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 4] <- 3
          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_contradictions_neg = stringr::str_count(toupper(`Contradiction(s)`), paste(CD8_negs, collapse='|')))

          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 1] <- 0
          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 2] <- 1
          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 3] <- 2
          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 4] <- 3
          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 5] <- 4
          final_matched_output3$CD8_contradictions_neg[is.na(final_matched_output3$CD8_contradictions_neg)] <- 0

          CD8_negs <- c('T-cell surface glycoprotein CD8 alpha chain \\(Negative\\)', 'T cell receptor co-receptor CD8 \\(Negative\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_cell_type_neg = stringr::str_count(`Full marker description`, paste(CD8_negs, collapse='|')))

          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 1] <- 0
          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 2] <- 1
          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 3] <- 2
          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 4] <- 3
          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 5] <- 4

          final_matched_output3$`total number of markers` <- final_matched_output3$`total number of markers` - final_matched_output3$CD8_cell_type_pos - final_matched_output3$CD8_cell_type_neg

          final_matched_output3$input_num_total <- final_matched_output3$input_num_total - final_matched_output3$CD8_input_pos - final_matched_output3$CD8_input_neg

          final_matched_output3$`num of inputs that match to cell type` <- final_matched_output3$`num of inputs that match to cell type` - final_matched_output3$CD8_match_neg - final_matched_output3$CD8_match_pos


          final_matched_output3$Num_pos_input <- final_matched_output3$Num_pos_input - final_matched_output3$CD8_input_pos
          final_matched_output3$Num_neg_input <- final_matched_output3$Num_neg_input - final_matched_output3$CD8_input_neg

          final_matched_output3$Num_pos_cell_type <- final_matched_output3$Num_pos_cell_type - final_matched_output3$CD8_cell_type_pos
          final_matched_output3$Num_neg_cell_type <- final_matched_output3$Num_neg_cell_type -final_matched_output3$CD8_cell_type_neg

          final_matched_output3$Num_pos_match <- final_matched_output3$Num_pos_match - final_matched_output3$CD8_match_pos
          final_matched_output3$Num_neg_match <- final_matched_output3$Num_neg_match - final_matched_output3$CD8_match_neg

          final_matched_output3$Num_pos_contr[is.na(final_matched_output3$Num_pos_contr)] <- 0
          final_matched_output3$Num_neg_contr[is.na(final_matched_output3$Num_neg_contr)] <- 0

          final_matched_output3$Num_pos_contr <- final_matched_output3$Num_pos_contr - final_matched_output3$CD8_contradictions_pos
          final_matched_output3$Num_neg_contr <- final_matched_output3$Num_neg_contr -final_matched_output3$CD8_contradictions_neg

          # Score
          final_matched_output3$Score <-
            ((as.numeric(final_matched_output3$Num_pos_match)) /
               (as.numeric(final_matched_output3$Num_pos_cell_type) + as.numeric(final_matched_output3$Num_pos_input) - as.numeric(final_matched_output3$Num_pos_match))) - ((as.numeric(final_matched_output3$Num_pos_contr) * 0.25) + (as.numeric(final_matched_output3$Num_neg_contr) * 0.25))

          final_matched_output3$Score <- as.numeric(final_matched_output3$Score)

          final_matched_output3$Score <- round(final_matched_output3$Score, 3)

          # Percent match (match markers/markers in cell type)
          final_matched_output3$`Percent match (match markers/markers in cell type)` <-
            (as.numeric(final_matched_output3$`num of inputs that match to cell type`) /
               (as.numeric(final_matched_output3$`total number of markers`)
               )) * 100

          final_matched_output3$`Percent match (match markers/markers in cell type)` <- as.numeric(final_matched_output3$`Percent match (match markers/markers in cell type)`)

          final_matched_output3$`Percent match (match markers/markers in cell type)` <- round(final_matched_output3$`Percent match (match markers/markers in cell type)`, 3)

          # Percent match (match markers/markers in input)
          final_matched_output3$`Percent match (match markers/markers in input)` <-
            ((as.numeric(final_matched_output3$`num of inputs that match to cell type`) /
                (as.numeric(final_matched_output3$input_num_total))) * 100)

          final_matched_output3$`Percent match (match markers/markers in input)` <- as.numeric(final_matched_output3$`Percent match (match markers/markers in input)`)

          final_matched_output3$`Percent match (match markers/markers in input)` <- round(final_matched_output3$`Percent match (match markers/markers in input)`, 3)

          # Get the columns for output, new order
          final_matched_output4 <- final_matched_output3[, c(
            "Cluster",
            "CL_term",
            "Cell_type",
            "Cell_comment",
            "Wikidata",
            "Matched inputted marker",
            "num of inputs that match to cell type",
            "Full input marker description",
            "input_num_total",
            "Full marker description",
            "total number of markers",
            "Contradiction(s)",
            "num contradiction(s)",
            "Percent match (match markers/markers in cell type)",
            "Percent match (match markers/markers in input)",
            "Score")]

          # New column names
          colnames(final_matched_output4)<-
            c(
              "Cluster",
              "Cell ontology ID",
              "Cell type name",
              "Cell type description",
              "Wikidata ID",
              "Matched markers",
              "# matched markers",
              "Inputted markers (also in reference)",
              "# inputted markers (also in reference)",
              "Cell type markers",
              "# cell type markers",
              "Contradictions",
              "# contradictions",
              "% (matched markers/markers in cell type)",
              "% (matched markers/markers in input)",
              "Score")

          # Order by score, get top N scores per cluster
          final_matched_output4 <-  dplyr::tbl_df(final_matched_output4) %>%
            dplyr::group_by(Cluster) %>%
            dplyr::filter(as.integer(ordered(-Score)) %in% 1:input$num_top_matches_1)

          final_matched_output4 <- final_matched_output4[with(final_matched_output4, order(Cluster, -Score)), ]

          return(final_matched_output4)
        } else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    }
  })

  # Step 6, make cell ontology term linkable in RShiny
  return_cell_types_linked_1 <- shiny::eventReactive(input$submit_tab1_step5, {

    final_matched_output <- return_cell_types_1()

    if (!is.null(final_matched_output)) {

      final_matched_output <- final_matched_output[, c(
        "Cluster",
        "Cell ontology ID",
        "Cell type name",
        "Cell type description",
        "Wikidata ID",
        "Matched markers",
        "Inputted markers (also in reference)",
        "Cell type markers",
        "Contradictions",
        "% (matched markers/markers in cell type)",
        "% (matched markers/markers in input)",
        "Score")]


      # Make look nicer (round and add a % symbol)
      final_matched_output$`% (matched markers/markers in cell type)` <-
        paste(round(final_matched_output$`% (matched markers/markers in cell type)`, 1), "%")

      final_matched_output$`% (matched markers/markers in input)` <-
        paste(round(final_matched_output$`% (matched markers/markers in input)`, 1), "%")

      # Just get cell ID to make link in R Shiny look better
      just_IDs <-
        gsub(
          "http://purl.obolibrary.org/obo/",
          "",
          final_matched_output$`Cell ontology ID`
        )

      just_IDs_wiki <-
        gsub(
          "http://www.wikidata.org/entity/",
          "",
          final_matched_output$`Wikidata ID`
        )

      # Remove '<' and '>' from URI
      final_matched_output$`Cell ontology ID` <-
        gsub('<', '', final_matched_output$`Cell ontology ID`)
      final_matched_output$`Cell ontology ID` <-
        gsub('>', '', final_matched_output$`Cell ontology ID`)
      just_IDs <- gsub('>', '', just_IDs)
      just_IDs <- gsub('<', '', just_IDs)

      # Make it so the link is clickable
      for (a in 1:nrow(final_matched_output)) {
        if (final_matched_output$`Cell ontology ID`[a] != "No match found") {
          final_matched_output$`Cell ontology ID`[a] <-
            paste0(
              "<a href='",
              final_matched_output$`Cell ontology ID`[a],
              "' target='_blank'>",
              just_IDs[a],
              "</a>"
            )
        }

        if (final_matched_output$`Wikidata ID`[a] != "") {
          if (lengths(regmatches(final_matched_output$`Wikidata ID`[a], gregexpr(",", final_matched_output$`Wikidata ID`[a]))) == 0) {
            final_matched_output$`Wikidata ID`[a] <- paste0(
              "<a href='",
              final_matched_output$`Wikidata ID`[a],
              "' target='_blank'>",
              just_IDs_wiki[a],
              "</a>"
            )
          } else {
            final_matched_output$`Wikidata ID`[a] <- paste0(
              "<a href='",
              strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][1],
              "' target='_blank'>",
              strsplit(just_IDs_wiki[a], ",")[[1]][1],
              "</a> , <a href='",
              strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][2],
              "' target='_blank'>",
              strsplit(just_IDs_wiki[a], ",")[[1]][2],
              "</a>"
            )
            if (lengths(regmatches(final_matched_output$`Wikidata ID`[a], gregexpr(",", final_matched_output$`Wikidata ID`[a]))) == 2) {
              final_matched_output$`Wikidata ID`[a] <- paste0(final_matched_output$`Wikidata ID`[a], "<a href='",
                                                              strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][3],
                                                              "' target='_blank'>",
                                                              strsplit(just_IDs_wiki[a], ",")[[1]][3],
                                                              "</a>")
            } else if (lengths(regmatches(final_matched_output$`Wikidata ID`[a], gregexpr(",", final_matched_output$`Wikidata ID`[a]))) == 3) {
              final_matched_output$`Wikidata ID`[a] <- paste0(
                "<a href='",
                strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][3],
                "' target='_blank'>",
                strsplit(just_IDs_wiki[a], ",")[[1]][3],
                "</a> , <a href='",
                strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][4],
                "' target='_blank'>",
                strsplit(just_IDs_wiki[a], ",")[[1]][4],
                "</a>"
              )
            }
          }
        }
      }

    #  final_matched_output$Cluster <- as.numeric(final_matched_output$Cluster)

      return(final_matched_output)
    } else {
      return(NULL)
    }
  })

  # Step 6, if multiple clusters, change dataframe if filtered by cluster
  df_subsettable_1 <- reactive({
   shiny::req(input$submit_tab1_step5)
    df <- return_cell_types_linked_1()
    if (!is.null(df) & length(unique(df$Cluster)) > 1) {
      df2 <- df[df$Cluster %in% input$show_specific_cluster_1,]
      return(df2)
    } else {
      return(df)
    }
  })

  # Step 6, get how many clusters are in final output for filtering purposes
  output$specific_cluster_1 <- shiny::renderUI({
   shiny::req(input$submit_tab1_step5)
    final_df <- return_cell_types_linked_1()
    if (!is.null(final_df) & length(unique(final_df$Cluster)) > 1) {
      shiny::selectInput(inputId = "show_specific_cluster_1",
                         label = "Filter by cluster:  ",
                         selected = sort(unique(as.character(final_df$Cluster))),
                         choices = sort(unique(as.character(final_df$Cluster))),
                         multiple = TRUE)
    }
  })

  # Step 6, when help button is clicked show message
  shiny::observeEvent(input$help_tab1_step6,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 6<br>Match Marker Patterns to Cell Types"),
      HTML("This step displays the results of the cell type matching process. The data frame shows:
      <br>
      <br>
      <li>Input clusters</li>
      <br>
      <li>Cell Ontology (CL) identification terms</li>
      <br>
      <li>Matched CL cell type names</li>
      <br>
      <li>Descriptions of the matched cell types</li>
      <br>
      <li>Wikidata identification terms</li>
      <br>
      <li>Markers that are matched between the input clusters and the CL cell types</li>
      <br>
      <li>Markers present in the input clusters</li>
      <br>
      <li>Markers present in the matched CL cell types</li>
      <br>
      <li>Markers that are directly contradicted between the input clusters and the CL cell types</li>
      <br>
      <br>
      Both CL and Wikidata terms are clickable and lead to their respective web pages.
      <br>
      <br>
      Three additional columns provide statistical metrics for evaluating the match:
      <br>
      <br>
      <li>Percentages of matched markers over the total markers in the CL cell type</li>
      <br>
      <li>Percentages of matched markers over the total markers in the input</li>
      <br>
      <li>The matching scores derived from Equation 3</li>
      <br>
      The 'Download' button at the bottom allows users to export a CSV file containing all relevant data, including additional columns detailing the number of markers in the input, CL cell types, matches, and contradictions.
      <br>
      <br>"),
      tags$div("Equation 3 and additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/3.1.-Input-Expression-Data-and-Default-References', "here.", target="_blank")),
      easyClose = TRUE))
  })


  ##################################################################
  ##                            Tab 2:                            ##
  ##                            Server                            ##
  ##                   Input marker descriptors                   ##
  ##################################################################
  ######~#~# Step 1 - Get marker synonyms and matched PRO/GO terms, user can input new names if left unmatched  #~#~######

  # Step 1, side panel

  # Reset everything
  shiny::observeEvent(input$reset_tab2_step1, {
    session$reload()
  })

  # Download example marker-cluster input data
  output$download_example <- shiny::downloadHandler(
    # File name
    filename = function() {
      c("Dusoswa_OMIP_54_markers.csv")
    },
    # Write to csv
    content = function(file) {
      utils::write.csv(Dusoswa_OMIP_54_markers,  file, row.names = FALSE)
    }
  )

  # Download example marker-cell type reference data
  output$download_example_marker_cell_type_2 <- shiny::downloadHandler(
    # File name
    filename = function() {
      c("Lee_AML_cell_types_markers.csv")
    },
    # Write to csv
    content = function(file) {
      utils::write.csv(Lee_AML_cell_types_markers,  file, row.names = FALSE)
    }
  )

  # Side panel help buttons

  # Help buttons for side panel parameters
  shiny::observeEvent(input$help_pregating_2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Type Pre-Gated Markers"),
      HTML("Any markers that were pre-gated should be inputted here. These markers will be treated as homogeneous across all clusters.
      <br>
      <br>
      For example, if the data was pre-gated for CD3 positive expression, enter 'CD3+' in the text input box. CD3 will be interpreted as positive for all clusters.
      <br>
      <br>
      If using internal references, CD45 will also be excluded from all analysis steps due to inconsistent usage in the CL."),
      easyClose = TRUE))
  })

  # Help box for input type
  shiny::observeEvent(input$help_input_type_2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Input Marker Definitions"),
      HTML("Upload Marker Definitions
      <br>
      <br>
      <ul>
      <li>The input file (CSV, XLS, or XLSX) should contain two columns: Cluster and Markers.</li>
      <ul>
      <br>
      <li>The Cluster column should contain the cluster number.</li>
      <br>
      <li>The Markers column should list protein marker names, with qualifiers such as positive, negative, low, high, +, ++, and -.</li>
      <br>
      <li>Each marker and qualifier should be separated by a comma (,).</li>
      <br>
      </ul>
      </ul>
      Type Marker Definitions
      <br>
      <br>
      <ul>
      <li>Instead of uploading a file containing the marker descriptions, markers and their qualifiers can be directly typed into a box within the application. This can only be used when querying a single marker definition/cluster.</li>
      <br>
      <ul>
      <li>The input should contain protein marker names, with qualifiers such as positive, negative, low, high, +, ++, and -.</li>
      <br>
      <li>Each marker and qualifier should be separated by a comma (,).</li>
      </ul>
      </ul>
       <br>"),
      tags$div(" An example marker definition file is available within the application (click 'Download Example' above the file upload box) and on the ", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/2.-Example-Data', "GitHub.", target="_blank")),
      easyClose = TRUE))
  })

  # Help box for reference type
  shiny::observeEvent(input$help_ref_type_2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select Cell Type Reference"),
      HTML("The user can choose between two options for the descriptive cell type reference:
      <br>
      <br>
       <ul>
      <li>Default Reference: This option does not require the user to upload a reference file. It uses various built-in references, including the Cell Ontology, Protein Ontology, Wikidata, and curated datasets included with the application.</li>
      <br>
      <li>File Upload: This option allows the user to upload a custom reference file containing descriptive cell type names and their associated marker definitions. The file should be a CSV with two columns:</li>
      <br>
      <ul>
      <li>Name: The cell type identifier</li>
      <br>
      <li>Markers: The protein marker names with associated qualifiers ('positive', 'negative', 'low', 'high', '+', '++', '-'). Markers and qualifiers should be separated by commas (,).</li>
      <br>
      </ul>
      </ul>"),
      tags$div("An example marker-cell type reference file is available within the application (click 'Download Example' above the file upload box) and on the ", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/2.-Example-Data', "GitHub.", target="_blank")),
      easyClose = TRUE))
  })

  # Help box for species type
  shiny::observeEvent(input$help_species_2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Specify Species"),
      HTML("Some proteins and cell type names are species-specific within ontologies. Therefore, when using 'Default References', specify the species of your sample(s)."),
      easyClose = TRUE))
  })

  # Help box for ontology types
  shiny::observeEvent(input$help_ontology_type_2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select Cell Type Ontology(s)"),
      HTML("When using Default References, specify the cell type ontology. The primary source should be the 'Cell Ontology', which contains hundreds of cell types with marker definitions.
      <br>
      <br>
      The 'Provisional Cell Ontology' can also be optionally included, containing cell types that are currently only provisionally defined but may still be useful.
      <br>
      <br>
      At least one cell type ontology must be included."),
      easyClose = TRUE))
  })

  # Help box for low marker sign
  shiny::observeEvent(input$help_low_upload,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select the Interpretation of 'Low'"),
      HTML("Specify how you would like the 'low' modifier to be interpreted. This is used to determine how inputted markers will be matched to markers within the reference cell types. Low can be interpreted as either:
      <br>
      <br>
      <ul>
      <li>'low'</li>
      <br>
      <li>'low' and 'positive'</li>
      <br>
      <li>'low' and 'negative'</li>"),
      easyClose = TRUE))
  })

  # Help box for high marker sign
  shiny::observeEvent(input$help_high_upload,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select the Interpretation of 'High'"),
      HTML("Specify how you would like the 'high' modifier to be interpreted. This is used to determine how inputted markers will be matched to markers within the reference cell types. High can be interpreted as either:
      <br>
      <br>
      <ul>
      <li>'high'</li>
      <br>
      <li>'high' and 'positive'</li>"),
      easyClose = TRUE))
  })

  # Step 1, main panel

  # Step 1, change format of uploaded input file or inputted text
  reformatted_data_2 <-  shiny::eventReactive(input$submit_tab2_step1, {

    if (is.null(check_input_file_2()) == TRUE) {

      # Upload file
      if (input$input_type == 'csv_upload') {

        file_input <- input$csv_upload

        if (endsWith(file_input$datapath, ".csv")) {

          inputted_df <- utils::read.csv(file_input$datapath, header = TRUE)

        } else if (endsWith(file_input$datapath, ".xls")) {

          inputted_df <- readxl::read_excel(file_input$datapath, col_names = TRUE)

        }  else if (endsWith(file_input$datapath, ".xlsx")) {

          inputted_df <- readxl::read_excel(file_input$datapath, col_names = TRUE)
        }

      } else if (input$input_type == 'type_input') {
        type_input <- input$type_input

        # Make dataframe
        inputted_df <- data.frame("Cluster" = "1", "Markers" = type_input)
      }

      # Ensure first letter is capitalized and the rest in lower case
      colnames(inputted_df) <- tolower(colnames(inputted_df))
      substr(colnames(inputted_df), 1, 1) <- toupper(substr(colnames(inputted_df), 1, 1))

      # Check that column names are correct. If not, return a warning and stop.
      if (!all(c("Cluster", "Markers") == colnames(inputted_df))) {
        return(NULL)
      }

      # Change column names
      colnames(inputted_df) <- c("Cluster", "Markers")

      # Split marker column by , or ;
      new_df <-
        data.frame(Cluster = inputted_df$Cluster, do.call('rbind', strsplit(
          as.character(inputted_df$Markers), split = "[,;]+"
        )))

      # Remove spaces
      new_df <- as.data.frame(lapply(new_df, function(y) gsub(' ', '', y)))

    # Remove duplicates per row
      new_df <-
        t(apply(new_df, 1, function(x)
          replace(x, duplicated(x), "")))

      # Number of cell type names
      num_names <- nrow(new_df)

      # Change dataframe format
      new_df <- reshape2::melt(new_df, id = "Cluster")

      # Match numbers to names
      melt_key <- new_df[(1:num_names), , drop = FALSE]

      # Continue getting correct format
      new_df <- new_df[-(1:num_names), , drop = FALSE]
      new_df <- new_df[, -2]
      new_df <- new_df[!(new_df$value == ""),]

      # Determine if positive or negative
      new_df$Sign[grepl("positive|\\+$", new_df$value, ignore.case = TRUE, useBytes = TRUE)] <-
        "Positive"
      new_df$Sign[grepl("negative|\\-$", new_df$value, ignore.case = TRUE, useBytes = TRUE)] <-
        "Negative"
      new_df$Sign[grepl("high|\\+\\+|\\+\\+\\+$", new_df$value, ignore.case =
                          TRUE, useBytes = TRUE)] <- "High"
      new_df$Sign[grepl("low$", new_df$value, ignore.case = TRUE, useBytes = TRUE)] <-
        "Low"

      # Remove the sign from the marker
      new_df$value <-
        gsub("\\+$", "", gsub("\\++$", "",  gsub("\\-$", "", gsub(
          "\\--$", "", gsub(
            "low$",
            "",
            ignore.case = TRUE,
            gsub(
              "high$",
              "",
              ignore.case = TRUE,
              gsub(
                "positive$",
                "",
                ignore.case = TRUE,
                gsub(
                  "negative$",
                  "",
                  as.character(new_df$value),
                  ignore.case = TRUE
                )
              )
            )
          )
        ))))

      # Add names
      new_df <- merge(new_df, melt_key, by = "Var1", all = TRUE)

      # Remove old, unneeded column
      new_df <- subset(new_df, select = -c(Var1, Var2))

      # New column names
      colnames(new_df) <- c("Marker", "Sign", "Cluster")

      MG_df <- added_MG_markers_2()

      if (is.data.frame(MG_df) == TRUE) {
        new_df <- rbind(new_df, MG_df)
      }

      # Make everything capital
      new_df$Marker <- toupper(new_df$Marker)

      if (input$reference_type_2 == 'internal_ref_2') {
        # Remove CD45
        new_df <- new_df[new_df$Marker != "CD45", ]
      }

      # Order and sort by cluster
      new_df <- new_df[order(new_df$Cluster), ]

      if (sum(is.na(new_df$Sign)) > 0) {
        if (sum(is.na(new_df$Sign)) > 1) {
          no_sign <- unique(new_df$Marker[is.na(new_df$Sign)])

          return(no_sign)

        } else if (sum(is.na(new_df$Sign)) == 1) {
          no_sign <- new_df$Marker[is.na(new_df$Sign)]

          return(no_sign)
        }
      } else {

        return(new_df)
      }
    }
  })

  # Step 1, uploaded input, error message if input has wrong column names
  output$ui_error_column_names_input_2 <- shiny::renderUI({
    shiny::req(input$input_type == 'csv_upload')
    shiny::req(input$submit_tab2_step1)

    if (is.null(check_input_file_2()) == TRUE) {
      if (is.null(reformatted_data_2()) == TRUE) {
        h4(
          "Error: Invalid uploaded input format. There should be 2 columns, the first labeled 'Cluster' and the second labeled 'Markers'."
        )
      }
    }
  })

  # Step 1, outputs warning if there is no sign for at least 1 marker in the input
  output$ui_warning_no_sign_2 <- shiny::renderUI({
    shiny::req(input$submit_tab2_step1)
    if (is.data.frame(reformatted_data_2()) == FALSE) {
      if (length(reformatted_data_2()) == 1) {
        tags$div(id = "invalid_sign_2", h4(
          paste(
            "Error: Invalid input format. The marker",
            reformatted_data_2(),
            "does not have a valid sign."
          )
        ))
      } else if (length(reformatted_data_2()) > 1) {
        no_sign <- reformatted_data_2()
        no_sign <- data.frame(x = no_sign)
        no_sign <- paste(no_sign$x, collapse = ", ")
        tags$div(id = "invalid_sign_2", h4(
          paste(
            "Error: Invalid input format. The markers",
            no_sign,
            "do not have valid signs."
          )
        ))
      }
    }
  })

  added_MG_markers_2 <- shiny::eventReactive(input$submit_tab2_step1, {

    P <- input$type_MG_markers_2

    if (!is.null(P)) {
      # Remove spaces
      P <- gsub(" ", "", P)

      if (P != "") {
        new_df <- utils::read.table(text = P, sep = ",")
        new_df <- as.data.frame(t(new_df))

        # Determine if positive or negative
        new_df$Sign[grepl("positive|\\+$", new_df$V1, ignore.case = TRUE)] <-
          "Positive"
        new_df$Sign[grepl("negative|\\-$", new_df$V1, ignore.case = TRUE)] <-
          "Negative"
        new_df$Sign[grepl("high|\\+\\+|\\+\\+\\+$", new_df$V1, ignore.case =
                            TRUE)] <- "High"
        new_df$Sign[grepl("low$", new_df$V1, ignore.case = TRUE)] <-
          "Low"

        # Remove the sign from the marker
        new_df$V1 <-
          gsub("\\+$", "", gsub("\\++$", "",  gsub("\\-$", "", gsub(
            "\\--$", "", gsub(
              "low$",
              "",
              ignore.case = TRUE,
              gsub(
                "high$",
                "",
                ignore.case = TRUE,
                gsub(
                  "positive$",
                  "",
                  ignore.case = TRUE,
                  gsub(
                    "negative$",
                    "",
                    as.character(new_df$V1),
                    ignore.case = TRUE
                  )
                )
              )
            )
          ))))

        # Uploaded file
        if (input$input_type == 'csv_upload') {

          file_input <- input$csv_upload

          if (endsWith(file_input$datapath, ".csv")) {

            inputted_df <- utils::read.csv(file_input$datapath, header = TRUE)

          } else if (endsWith(file_input$datapath, ".xls")) {

            inputted_df <- readxl::read_excel(file_input$datapath, col_names = TRUE)

          }  else if (endsWith(file_input$datapath, ".xlsx")) {

            inputted_df <- readxl::read_excel(file_input$datapath, col_names = TRUE)
          }

        } else if (input$input_type == 'type_input') {
          type_input <- input$type_input

          # Make dataframe
          inputted_df <- data.frame("Cluster" = "1", "Markers" = type_input)
        }

        # Change column names
        colnames(inputted_df) <- c("Cluster", "Markers")

        names(inputted_df) <- tolower(names(inputted_df))

        unique_clusters <- unique(inputted_df$cluster)

        names(new_df) <- tolower(names(new_df))

        new_df <- merge(new_df, unique_clusters)

        # New column names
        colnames(new_df) <- c("Marker", "Sign", "Cluster")

        return(new_df)
      } else {
        P <- ""
        return(P)
      }
    } else {
      P <- ""
      return(P)
    }
  })

  # Step 1, check if file type of uploaded input file is acceptable
  check_input_file_2 <- shiny::eventReactive(input$submit_tab2_step1, {

    if (input$input_type == 'csv_upload') {

      # File inputted
      file_input <- input$csv_upload

      if (endsWith(file_input$datapath, ".csv") | endsWith(file_input$datapath, ".xls") | endsWith(file_input$datapath, ".xlsx")) {
        return(NULL)
      } else {
        return("a")
      }
    } else {
      return(NULL)
    }
  })

  # Step 1, outputs error stating that the input file is in the wrong format
  output$ui_error_check_input_file_2 <- shiny::renderUI({
    shiny::req(input$submit_tab2_step1)
    if (is.null(check_input_file_2()) == FALSE) {
      tags$div(id = "invalid_input_file_2", h4(
        "Error: Invalid uploaded input file format. Please upload a valid comma separated values (CSV) or Excel (XLSX or XLS) file."
      ))
    }})

  # # Step 1, outputs error if there is a marker designated twice or more in a cluster
  marker_sign_error_2 <- shiny::eventReactive(input$submit_tab2_step1, {

    df <- reformatted_data_2()

    same_marker_cluster <- data.frame(table(df$Cluster, df$Marker))

    same_marker_cluster <- same_marker_cluster[same_marker_cluster$Freq > 1,]

    if (nrow(same_marker_cluster) >= 1) {
      P <- "1"
      return(P)
    }
  })

  output$ui_marker_sign_error_2 <- shiny::renderUI({
    shiny::req(input$submit_tab2_step1)
    if (length(marker_sign_error_2()) >= 1) {
      tags$div(id = "id_marker_sign_error_2", h4(
        "Error: A marker is repeated at least twice within 1 or more clusters. Please ensure there are no repeated markers within each cluster (including specified manually gated markers)."
      ))}
  })

  #~#~#~#~#~#  Uploaded reference #~#~#~#~#~#

  # Step 1, uploaded reference, when help button is clicked show message
  shiny::observeEvent(input$help_uploaded_ref_2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Match Marker Patterns to Cell Types"),
      HTML("This step displays the results of the cell type matching process. The data frame shows:
      <br>
      <br>
      <li>Input clusters</li>
      <br>
      <li>Matched reference cell type names</li>
      <br>
      <li>Markers that are matched between the input clusters and the reference cell types</li>
      <br>
      <li>Markers present in the input clusters</li>
      <br>
      <li>Markers present in the matched reference cell types</li>
      <br>
      <li>Markers that are directly contradicted between the input clusters and the reference cell types</li>
      <br>
      <br>
      Three additional columns provide statistical metrics for evaluating the match:
      <br>
      <br>
      <li>Percentages of matched markers over the total markers in the reference cell type</li>
      <br>
      <li>Percentages of matched markers over the total markers in the input</li>
      <br>
      <li>The matching scores derived from Equation 3</li>
      <br>
      The 'Download' button at the bottom allows users to export a CSV file containing all relevant data, including additional columns detailing the number of markers in the input, reference cell types, matches, and contradictions.
      <br>
      <br>
     "),
      tags$div("Equation 3 and additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/4.2.-Input-Marker-Descriptors-and-Uploaded-References', "here.", target="_blank")),
      easyClose = TRUE))
  })

  # Step 1, uploaded reference, check if file type of uploaded reference file is acceptable
  check_reference_file_2 <- shiny::eventReactive(input$submit_tab2_step1, {
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)

    # File inputted
    file_input <- input$upload_ref_csv_2

    if (endsWith(file_input$datapath, ".csv") | endsWith(file_input$datapath, ".xls") | endsWith(file_input$datapath, ".xlsx")) {
      return(NULL)
    }  else {
      return("a")
    }
  })

  # Step 1, uploaded reference, outputs error stating that the uploaded reference file is in the wrong format
  output$ui_error_check_reference_file_2 <- shiny::renderUI({
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)
    shiny::req(input$submit_tab2_step1)
    if (is.null(check_reference_file_2()) == FALSE) {
      tags$div(id = "invalid_ref_file_2", h4(
        "Error: Invalid uploaded reference format. Please upload a valid comma separated values (CSV) or Excel (XLSX or XLS) file."
      ))
    }})

  # Step 1, uploaded reference, change format of uploaded reference file
  reformatted_ref_2 <- shiny::eventReactive(input$submit_tab2_step1, {
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)

    if (is.null(check_reference_file_2()) == TRUE) {

      # File inputted
      file_input <- input$upload_ref_csv_2

      if (endsWith(file_input$datapath, ".csv")) {

        inputted_df <- utils::read.csv(file_input$datapath, header = TRUE)

      } else if (endsWith(file_input$datapath, ".xls")) {

        inputted_df <- readxl::read_excel(file_input$datapath, col_names = TRUE)

      }  else if (endsWith(file_input$datapath, ".xlsx")) {

        inputted_df <- readxl::read_excel(file_input$datapath, col_names = TRUE)
      }

      # Normalize all dashes
      inputted_df <- as.data.frame(lapply(inputted_df, function(y) gsub("\\p{Pd}", "-", y, perl=TRUE)))

      # Ensure first letter is capitalized and the rest in lower case
      colnames(inputted_df) <- tolower(colnames(inputted_df))
      substr(colnames(inputted_df), 1, 1) <- toupper(substr(colnames(inputted_df), 1, 1))

      # Check that column names are correct. If not, return a warning and stop.
      if (!all(c("Name", "Markers") == colnames(inputted_df))) {
        return(NULL)
      } else {
        # Split marker column by , or ;
        new_df <-
          data.frame(Name = inputted_df$Name, do.call('rbind', strsplit(
            as.character(inputted_df$Markers), split = "[,;]+"
          )))

        # Remove spaces
       # new_df <- as.data.frame(lapply(new_df, function(y) gsub(' ', '', y)))
        new_df[2:ncol(new_df)] <- as.data.frame(lapply(new_df[2:ncol(new_df)], function(y) gsub(' ', '', y)))

        # Remove duplicates per row
        new_df <-
          t(apply(new_df, 1, function(x)
            replace(x, duplicated(x), "")))

        # Number of cell type names
        num_names <- nrow(new_df)

        # Change dataframe format
        new_df <- reshape2::melt(new_df, id = "Name")

        # Match numbers to names
        melt_key <- new_df[(1:num_names), , drop = FALSE]

        # Continue getting correct format
        new_df <- new_df[-(1:num_names), , drop = FALSE]
        new_df <- new_df[, -2]
        new_df <- new_df[!(new_df$value == ""),]

        # Determine if positive or negative
        new_df$Sign[grepl("positive|\\+$", new_df$value, ignore.case = TRUE)] <-
          "Positive"
        new_df$Sign[grepl("negative|\\-$", new_df$value, ignore.case = TRUE)] <-
          "Negative"
        new_df$Sign[grepl("high|\\+\\+|\\+\\+\\+$", new_df$value, ignore.case =
                            TRUE)] <- "High"
        new_df$Sign[grepl("low$", new_df$value, ignore.case = TRUE)] <-
          "Low"

        # Remove the sign from the marker
        new_df$value <-
          gsub("\\+$", "", gsub("\\++$", "",  gsub("\\-$", "", gsub(
            "\\--$", "", gsub(
              "low$",
              "",
              ignore.case = TRUE,
              gsub(
                "high$",
                "",
                ignore.case = TRUE,
                gsub(
                  "positive$",
                  "",
                  ignore.case = TRUE,
                  gsub(
                    "negative$",
                    "",
                    as.character(new_df$value),
                    ignore.case = TRUE
                  )
                )
              )
            )
          ))))

        # Add names
        new_df <- merge(new_df, melt_key, by = "Var1", all = TRUE)

        # Remove old, unneeded column
        new_df <- subset(new_df, select = -c(Var1, Var2))

        # New column names
        colnames(new_df) <- c("Marker", "Sign", "Name")

        if (sum(is.na(new_df$Sign)) > 0) {
          if (sum(is.na(new_df$Sign)) > 1) {
            no_sign <- unique(new_df$Marker[is.na(new_df$Sign)])

            return(no_sign)

          } else if (sum(is.na(new_df$Sign)) == 1) {
            no_sign <- new_df$Marker[is.na(new_df$Sign)]

            return(no_sign)
          }
        } else {

          # Make everything uppercase for the filter
          new_df$Marker <- toupper(new_df$Marker)

          return(new_df)
        }
      }
    }
  })

  # Step 1, uploaded reference, outputs warning if there is no sign for at least 1 marker in the uploaded reference
  output$ui_warning_no_sign_ref_2 <- shiny::renderUI({
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)
    shiny::req(input$submit_tab2_step1)

    if (is.data.frame(reformatted_ref_2()) == FALSE) {
      if (length(reformatted_ref_2()) == 1) {
        tags$div(id = "invalid_sign_2", h4(
          paste(
            "Error: Invalid uploaded reference format. The marker",
            reformatted_ref_2(),
            "does not have a valid sign."
          )
        ))
      } else if (length(reformatted_ref_2()) > 1) {
        no_sign <- reformatted_ref_2()
        no_sign <- data.frame(x = no_sign)
        no_sign <- paste(no_sign$x, collapse = ", ")
        tags$div(id = "invalid_sign_2", h4(
          paste(
            "Error: Invalid uploaded reference format. The markers",
            no_sign,
            "do not have valid signs."
          )
        ))
      }
    }
  })

  # Step 1, uploaded reference, error message if uploaded reference has wrong column names
  output$ui_error_column_names_ref_2 <- shiny::renderUI({
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)
    shiny::req(input$submit_tab2_step1)

    if (is.null(check_reference_file_2()) == TRUE) {
      if (is.null(reformatted_ref_2()) == TRUE) {
        h4(
          "Error: Invalid uploaded reference format. There should be 2 columns, the first labeled 'Name' and the second labeled 'Markers'."
        )
      }
    }
  })

  # Step 1, uploaded reference, get warning statement indicating if/which markers are not in uploaded reference
  output$ui_warning_unmatched_markers_uploaded_ref_2 <- shiny::renderUI({
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)
    shiny::req(input$submit_tab2_step1)

    if (is.data.frame(reformatted_data_2()) == TRUE &
        is.data.frame(reformatted_ref_2()) == TRUE) {
      if (length(marker_sign_error_2()) == 0) {
        new_df <- reformatted_data_2()
        ref_df <- reformatted_ref_2()

        different_markers <- setdiff(new_df$Marker, ref_df$Marker)

        if (length(different_markers) < length(unique(new_df$Marker))) {
          if (length(different_markers) > 1) {
            marker_string <-
              paste(different_markers, collapse = ", ")
            h4(
              paste(
                "Warning: The inputted markers",
                marker_string,
                "are not in the uploaded reference."
              )
            )
          } else if (length(different_markers) == 1) {
            h4(
              paste(
                "Warning: The inputted marker",
                different_markers,
                "is not in the uploaded reference."
              )
            )
          }
        } else {
          h4("Error: None of the inputted markers are in the uploaded reference.")
        }
      }
    }
  })

  # Step 1, uploaded reference, get results on how input and uploaded reference match
  uploaded_ref_matches_2 <- shiny::eventReactive(input$submit_tab2_step1, {
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)

    # Make sure the uploaded reference and input are dataframes
    if (is.data.frame(reformatted_data_2()) == TRUE &
        is.data.frame(reformatted_ref_2()) == TRUE) {
      new_df <- reformatted_data_2()
      ref_df <- reformatted_ref_2()

      different_markers <- setdiff(new_df$Marker, ref_df$Marker)

      if (length(different_markers) < length(unique(new_df$Marker))) {

        # Look at each cluster individually
        each_cluster <- split(new_df, new_df$Cluster)

        # Make empty, will be added to (+1) if there is a match
        ref_df[, 'Match'] <- 0

        # Empty list to put results in
        clusters_to_cell_types <- list()

        # Loop through each cluster of input
        for (x in 1:length(each_cluster)) {
          ref_df_new <- ref_df
          # Loop through each row of input and reference
          # First check to see if marker matches, then if the sign matches
          # If both match, add to the 'match' column
          for (a in 1:nrow(each_cluster[[x]])) {
            for (b in 1:nrow(ref_df_new)) {

              ## Matches ##

              if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                  each_cluster[[x]]$Sign[a] == "Positive" &
                  ref_df_new$Sign[b] == "Positive") {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Positive" &
                         ref_df_new$Sign[b] == "High") {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Positive" &
                         ref_df_new$Sign[b] == "Low") {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Negative" &
                         ref_df_new$Sign[b] == "Negative") {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Low" &
                         ref_df_new$Sign[b] == "Low") {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "High" &
                         ref_df_new$Sign[b] == "High") {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
              } else {
                ref_df_new$Match[b] <- ref_df_new$Match[b] + 0
              }

              # If the high/low positive box was clicked, then also consider 'low' and 'high' inputs as positive
              if (input$low_option_upload == "low_positive_upload") {
                if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                    each_cluster[[x]]$Sign[a] == "Low" &
                    ref_df_new$Sign[b] == "Positive") {
                  ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
                } else {
                  ref_df_new$Match[b] <- ref_df_new$Match[b] + 0
                }}

              if (input$high_option_upload == "high_positive_upload") {
                if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                    each_cluster[[x]]$Sign[a] == "High" &
                    ref_df_new$Sign[b] == "Positive") {
                  ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
                } else {
                  ref_df_new$Match[b] <- ref_df_new$Match[b] + 0
                }
              }

              # If the low negative box was clicked, then also consider 'low' input as negative
              if (input$low_option_upload == "low_negative_upload") {
                if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                    each_cluster[[x]]$Sign[a] == "Low" &
                    ref_df_new$Sign[b] == "Negative") {
                  ref_df_new$Match[b] <- ref_df_new$Match[b] + 1
                } else {
                  ref_df_new$Match[b] <- ref_df_new$Match[b] + 0
                }}}}

          # Remove unmatched
          ref_df_new <- subset(ref_df_new, Match != 0)

          # Sort, put into final format for export
          final_matched_output <- ref_df_new %>%
            dplyr::group_by(Name) %>%
            dplyr::summarise(
              `Matched inputted markers` = paste0(Marker, " (", Sign, ")", collapse = ", "),
              `num of inputs that match to cell type` = dplyr::n()
            )

          # Put into list
          cluster_name <- each_cluster[[x]]$Cluster[1]
          clusters_to_cell_types[[cluster_name]] <-
            as.data.frame(final_matched_output)
        }

        # Put all results in one data frame
        all_cluster_matches <-
          dplyr::bind_rows(clusters_to_cell_types, .id = "Cluster")

        ## Contradictions ##

        # Make empty, will be added to (+1) if there is a contradiction
        ref_df[, 'Contradiction'] <- 0

        # Empty list to put results in
        clusters_to_cell_types2 <- list()

        # Loop through each cluster of input
        for (x in 1:length(each_cluster)) {
          ref_df_new <- ref_df
          # Loop through each row of input and reference
          # First check to see if marker matches, then if the sign matches
          # If both match, add to the 'match' column
          for (a in 1:nrow(each_cluster[[x]])) {
            for (b in 1:nrow(ref_df_new)) {

              if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                  each_cluster[[x]]$Sign[a] == "Positive" &
                  ref_df_new$Sign[b] == "Negative") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Negative" &
                         ref_df_new$Sign[b] == "Positive") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "High" &
                         ref_df_new$Sign[b] == "Negative") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Negative" &
                         ref_df_new$Sign[b] == "High") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "Low" &
                         ref_df_new$Sign[b] == "High") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                         each_cluster[[x]]$Sign[a] == "High" &
                         ref_df_new$Sign[b] == "Low") {
                ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
              }

              # If the low negative box was clicked, then also consider 'low' input as negative
              if (input$low_option_upload == "low_negative_upload") {
                if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                    each_cluster[[x]]$Sign[a] == "Low" &
                    ref_df_new$Sign[b] == "Positive") {
                  ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
                } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                           each_cluster[[x]]$Sign[a] == "Positive" &
                           ref_df_new$Sign[b] == "Low") {
                  ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
                }}


              # If choosing low/positive, then any low signs are contradicted by negatives and highs
              if (input$low_option_upload == "low_positive_upload") {
                if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                    each_cluster[[x]]$Sign[a] == "Low" &
                    ref_df_new$Sign[b] == "Negative") {
                  ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
                } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                           each_cluster[[x]]$Sign[a] == "Negative" &
                           ref_df_new$Sign[b] == "Low") {
                  ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
                }}

              # If choosing high/high, then any high signs are contradicted by negatives, positives, and lows
              if (input$high_option_upload == "high_high_upload") {
                if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                    each_cluster[[x]]$Sign[a] == "High" &
                    ref_df_new$Sign[b] == "Positive") {
                  ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
                } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                           each_cluster[[x]]$Sign[a] == "Positive" &
                           ref_df_new$Sign[b] == "High") {
                  ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
                }}

              # If choosing low/low, then any low signs are contradicted by negatives, positives, and high
              if (input$low_option_upload == "low_low_upload") {
                if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                    each_cluster[[x]]$Sign[a] == "Low" &
                    ref_df_new$Sign[b] == "Positive") {
                  ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
                } else if (each_cluster[[x]]$Marker[a] == ref_df_new$Marker[b] &
                           each_cluster[[x]]$Sign[a] == "Positive" &
                           ref_df_new$Sign[b] == "Low") {
                  ref_df_new$Contradiction[b] <- ref_df_new$Contradiction[b] + 1
                }}}}

          # Remove unmatched
          ref_df_new2 <- subset(ref_df_new, Contradiction != 0)

          # Sort, put into final format for export
          final_matched_output2 <- ref_df_new2 %>%
            dplyr::group_by(Name) %>%
            dplyr::summarise(
              `Contradicted markers` = paste0(Marker, " (", Sign, ")", collapse = ", "),
              `num of contradictions` = dplyr::n()
            )

          # Put into list
          cluster_name2 <- each_cluster[[x]]$Cluster[1]
          clusters_to_cell_types2[[cluster_name2]] <-
            as.data.frame(final_matched_output2)
        }

        # Put all results in one data frame
        all_cluster_matches2 <-
          dplyr::bind_rows(clusters_to_cell_types2, .id = "Cluster")

        all_cluster_matches <-
          merge(all_cluster_matches, all_cluster_matches2, by = c("Cluster", "Name"), all.x = TRUE)

        # Get what and how many total markers in each cell type
        final_matched_output <- ref_df %>%
          dplyr::group_by(Name) %>%
          dplyr::summarise(
            `Cell type markers` = paste0(Marker, " (", Sign, ")", collapse = ", "),
            `total number of markers` = dplyr::n()
          )

        all_cluster_matches <-
          merge(all_cluster_matches, final_matched_output, by = "Name", all.x = TRUE)

        # Get what and how many total markers in each input
        final_matched_output <- new_df %>%
          dplyr::group_by(Cluster) %>%
          dplyr::summarise(
            `Inputted markers` = paste0(Marker, " (", Sign, ")", collapse = ", "),
            input_num_total = dplyr::n()
          )

        all_cluster_matches <-
          merge(all_cluster_matches, final_matched_output, by = "Cluster")

        # Count how many markers are positive and negative
        all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_pos_input = stringr::str_count(`Inputted markers`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
        all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_neg_input = stringr::str_count(`Inputted markers`, '(Negative)'))

        all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_pos_cell_type = stringr::str_count(`Cell type markers`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
        all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_neg_cell_type = stringr::str_count(`Cell type markers`, '(Negative)'))

        all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_pos_match = stringr::str_count(`Matched inputted markers`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
        all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_neg_match = stringr::str_count(`Matched inputted markers`, '(Negative)'))

        all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_pos_contr = stringr::str_count(`Contradicted markers`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
        all_cluster_matches <- all_cluster_matches %>% dplyr::mutate(Num_neg_contr = stringr::str_count(`Contradicted markers`, 'Negative'))

        # Make empty values in number of contradictions NA values
        all_cluster_matches$`num of contradictions`[is.na(all_cluster_matches$`num of contradictions`)] <- 0
        all_cluster_matches$Num_pos_contr[is.na(all_cluster_matches$Num_pos_contr)] <- 0
        all_cluster_matches$Num_neg_contr[is.na(all_cluster_matches$Num_neg_contr)] <- 0

        # Percent match (match markers/markers in cell type)
        all_cluster_matches$`Percent match (match markers/markers in cell type)` <-
          (as.numeric(all_cluster_matches$`num of inputs that match to cell type`) /
             (as.numeric(all_cluster_matches$`total number of markers`)
             )) * 100

        all_cluster_matches$`Percent match (match markers/markers in cell type)` <- as.numeric(all_cluster_matches$`Percent match (match markers/markers in cell type)`)

        all_cluster_matches$`Percent match (match markers/markers in cell type)` <- round(all_cluster_matches$`Percent match (match markers/markers in cell type)`, 3)

        # Percent match (match markers/markers in input)
        all_cluster_matches$`Percent match (match markers/markers in input)` <-
          ((as.numeric(all_cluster_matches$`num of inputs that match to cell type`) /
              (as.numeric(all_cluster_matches$input_num_total))) * 100)

        all_cluster_matches$`Percent match (match markers/markers in input)` <- as.numeric(all_cluster_matches$`Percent match (match markers/markers in input)`)

        all_cluster_matches$`Percent match (match markers/markers in input)` <- round(all_cluster_matches$`Percent match (match markers/markers in input)`, 3)

        # Score
        all_cluster_matches$Score <-
          ((as.numeric(all_cluster_matches$Num_pos_match)) /
             (as.numeric(all_cluster_matches$Num_pos_cell_type) + as.numeric(all_cluster_matches$Num_pos_input) - as.numeric(all_cluster_matches$Num_pos_match))) - ((as.numeric(all_cluster_matches$Num_pos_contr) * 0.25) + (as.numeric(all_cluster_matches$Num_neg_contr) * 0.25))

        all_cluster_matches$Score <- as.numeric(all_cluster_matches$Score)

        all_cluster_matches$Score <- round(all_cluster_matches$Score, 3)

        # Make sure there are results
        if (nrow(all_cluster_matches) != 0) {

          # Get the columns for output, new order
          all_cluster_matches <-
            all_cluster_matches[, c("Cluster",
                                    "Name",
                                    "Matched inputted markers",
                                    "num of inputs that match to cell type",
                                    "Inputted markers",
                                    "input_num_total",
                                    "Cell type markers",
                                    "total number of markers",
                                    "Contradicted markers",
                                    "num of contradictions",
                                    "Percent match (match markers/markers in cell type)",
                                    "Percent match (match markers/markers in input)",
                                    "Score"
            )]

          # Final name
          colnames(all_cluster_matches) <-
            c(
              "Cluster",
              "Cell type name",
              "Matched inputted markers",
              "# matched inputted markers",
              "Inputted markers",
              "# inputted markers",
              "Cell type markers",
              "# cell type markers",
              "Contradictions",
              "# contradictions",
              "% (matched markers/markers in cell type)",
              "% (matched markers/markers in input)",
              "Score"
            )

          # Order by score
          all_cluster_matches <-
            all_cluster_matches[order(all_cluster_matches$Score, decreasing = TRUE), ]

          return(all_cluster_matches)
        } else {
          return("p")
        }
      }  else {
        return("p")
      }
    }
  })

  # Step 1, uploaded reference, changes when the 'submit' button is hit if using an uploaded reference
  shiny::observeEvent(input$submit_tab2_step1, {
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)

    shiny.destroy::removeOutput("app_overview_2")
    shinyjs::show("overview_uploaded_ref_2")

    if (is.null(check_input_file_2)) {
      if (is.data.frame(reformatted_data_2()) == TRUE &
          is.data.frame(reformatted_ref_2()) == TRUE) {
        if (is.data.frame(uploaded_ref_matches_2()) == FALSE) {
          shinyjs::show("no_matched_cell_types_2")
        } else if (nrow(uploaded_ref_matches_2()) == 0) {
          shinyjs::show("no_matched_cell_types_2")
        }
        else {
          shinyjs::show("help_uploaded_ref_2")
          shinyjs::show("download_uploaded_ref_table_cell_types_2")
        }
      }
    }
  })

  uploaded_ref_results_2 <- shiny::eventReactive(input$submit_tab2_step1, {

    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)

    matches <- uploaded_ref_matches_2()

    if (is.data.frame(matches) == TRUE) {

      # Remove number columns (still included in download)
      matches <-
        matches[,     c(
          "Cluster",
          "Cell type name",
          "Matched inputted markers",
          "Inputted markers",
          "Cell type markers",
          "Contradictions",
          "% (matched markers/markers in cell type)",
          "% (matched markers/markers in input)",
          "Score"
        )]

      # Make look nicer (round and add a % symbol)
      matches$`% (matched markers/markers in cell type)` <-
        paste(round(matches$`% (matched markers/markers in cell type)`, 1), "%")

      matches$`% (matched markers/markers in input)` <-
        paste(round(matches$`% (matched markers/markers in input)`, 1), "%")

      # Order by score
      matches <- matches[with(matches, order(Cluster, -Score)), ]

      return(matches)
    } else {
      P <- "1"
      return(P)
    }
  })

  # Step 1, uploaded reference, make datatable of matched cell types
  output$uploaded_ref_matches_2_dt <- DT::renderDT({
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)
    shiny::req(input$submit_tab2_step1)

    df <- df_subsettable_uploaded_2()

    if (is.data.frame(df) == FALSE) {
    } else if (nrow(df) == 0) {
    }
    else {
      DT::datatable(
        df,
        selection = 'none',
        rownames = FALSE,
        escape = FALSE,
        options = list(pageLength = 25)
      )
    }
  }, server=FALSE)

  # Step 1, uploaded reference, if multiple clusters, change dataframe if filtered by cluster
  df_subsettable_uploaded_2 <- reactive({
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)
    shiny::req(input$submit_tab2_step1)

    df <- uploaded_ref_results_2()
    if (!is.null(df) & length(unique(df$Cluster)) > 1) {
    #  df$Cluster <- as.numeric(df$Cluster)
      df2 <- df[df$Cluster %in% input$show_specific_cluster_uploaded_2,]
      return(df2)
    } else {
      return(df)
    }
  })

  # Step 1, uploaded reference, get how many clusters are in final output for filtering purposes
  output$specific_cluster_uploaded_2 <- shiny::renderUI({
    shiny::req(input$reference_type_2 == 'upload_ref_2')
    shiny::req(input$upload_ref_csv_2)
    shiny::req(input$submit_tab2_step1)

    final_df <- uploaded_ref_results_2()
    b <- sort(unique(as.character(final_df$Cluster)))
    if (!is.null(final_df) & length(unique(final_df$Cluster)) > 1) {
      shiny::selectInput(inputId = "show_specific_cluster_uploaded_2",
                         label = "Filter by cluster:  ",
                         selected = sort(unique(as.character(final_df$Cluster))),
                         choices = sort(unique(as.character(final_df$Cluster))),
                         multiple = TRUE)
    }
  })

  # Step 1, uploaded reference, download final results
  output$download_uploaded_ref_table_cell_types_2 <- shiny::downloadHandler(
    # File name
    filename = function() {
      if (input$input_type == 'csv_upload') {
        paste0(gsub(".csv|.xls|.xlsx", "", input$csv_upload),
               "-matched-cell-types.csv")
      } else {
        paste0("typed-markers-matched-cell-types.csv")
      }
    },

    # Write to csv
    content = function(file) {
      utils::write.csv(uploaded_ref_matches_2(),  file, row.names = FALSE)
    }
  )
  #   #~#~#~#~#~# Internal references #~#~#~#~#~#

  # Step 1, when help button is clicked show message
  shiny::observeEvent(input$help_tab2_step1,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 1<br>Match Markers to PRO/GO Terms"),
      HTML("In this step, the inputted markers are processed through a multistep workflow to find matching ontology terms. If a marker cannot be automatically matched, the tool will return it to the user for manual name edits.
      <br>
      <br>
      A data frame is displayed, showing the original marker names and any suggested names. The last column, 'Type New Marker Name', is editable, allowing users to propose new names. These new names will be reprocessed to check for matches within the multistep workflow.
      <br>
      <br>"),
      tags$div("Additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/4.1.-Input-Marker-Descriptors-and-Default-References', "here.", target="_blank")),
      easyClose = TRUE))
  })

  # Step 1, create editable dataframes
  df_1_2 <- shiny::reactiveValues(data = NULL)

  # Step 1, changes when the 'submit' button is hit if using the internal references
  shiny::observeEvent(input$submit_tab2_step1, {
    shiny::req(input$reference_type_2 == 'internal_ref_2')

    shiny.destroy::removeOutput("app_overview_2")

    shinyjs::show("overview_tab2_step1")

    if (is.null(check_input_file_2())) {

      if (is.data.frame(reformatted_data_2()) == TRUE) {

        if (length(marker_sign_error_2()) == 0) {

          df_1_2$data <- return_alt_markers_2()

          # Step 1, make datatable of unmatched markers to PRO terms
          output$table_new_format_2 <- DT::renderDT({

            if (is.null(df_1_2$data)) {
            } else {
              DT::datatable(
                df_1_2$data,
                selection = 'none',
                editable = list(target = "cell", disable = list(columns = c(0, 1))),
                rownames = FALSE,
                escape = FALSE,
                options = list(pageLength = 25)
              )
            }
          })


          # if (is.data.frame(reformatted_data_2()) == TRUE) {
          #  if (!is.null(return_unmatched_markers_2())) {
          shinyjs::show("help_tab2_step1")
          # }
          shiny.destroy::removeOutput("all_sidebar_tab2_step1")
          shinyjs::show("all_sidebar_tab2_step2")
        }
      }
    }
  })

  # Step 1, if some markers were not matched to PRO terms, prompt the user to rename them
  return_alt_markers_2 <- shiny::eventReactive(input$submit_tab2_step1, {
    if (is.data.frame(reformatted_data_2()) == TRUE) {
      all_ids <- PRO_results_2()

      # Just markers that did not match
      all_ids <- all_ids[all_ids$`PRO/GO Term` == "",]

      # Get certain columns
      all_ids <- all_ids[c("Original Name", "Suggestion")]

      if (nrow(all_ids) != 0) {

        # If no suggestion, link to PRO
        all_ids$Suggestion[all_ids$Suggestion == ""]<- "<a href='https://proconsortium.org/' target='_blank'>Search PRO</a>"

        # New column
        all_ids$'Type new marker name' <- ""

        # Remove duplicates
        all_ids <- all_ids[!duplicated(all_ids), ]

        if (nrow(all_ids) != 0) {
          return(all_ids)
        }  else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    }
  })

  #   # Step 1, match inputted markers to PRO terms
  PRO_results_2 <- shiny::eventReactive(input$submit_tab2_step1, {

    if (is.data.frame(reformatted_data_2()) == TRUE) {

      ## Species, references, initial dataframe

      new_df <- reformatted_data_2()

      # Get NCBI taxon ID
      species_ID <- input$species_choice_2
      matched_row <-
        which(species_in_PRO == species_ID, arr.ind = TRUE)[1]
      NCBI_taxon_ID <-
        species_in_PRO[matched_row, ncol(species_in_PRO)]

      NCBI_taxon_ID_full <- species_in_PRO[matched_row, 1]
      NCBI_taxon_ID_short <- gsub(".*/","",NCBI_taxon_ID_full)

      # Get Wikidata taxon ID
      SPARQL_query <-

        paste0(
          "SELECT ?item ?itemLabel WHERE {

     ?item wdt:P2888 <", NCBI_taxon_ID_full, ">.

    SERVICE wikibase:label { bd:serviceParam wikibase:language 'en'. }
    }"
        )

      # Run SPARQL on the endpoint
      result <- httr::GET(
        url = wiki_endpoint,
        query = list(query = SPARQL_query),
        httr::user_agent(R.version.string))

      # Will show a warning/error if there is any
      httr::stop_for_status(result)

      # Get result in text JSON
      x <- httr::content(result, as = "text") #, encoding = "UTF-8")

      # Convert from JSON to a list
      df <- jsonlite::fromJSON(x, flatten = TRUE)

      # Extract the data frame
      df <- df$results$bindings

      if (length(df) != 0) {

        # Remove unneeded info
        df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", ':lang')))

        # Remove this part that was added onto the column names
        colnames(df) <- gsub(".value", "", colnames(df))
      }

      wikidata_taxon_ID <- df$item

      # List of specific protein suggestions
      marker_alt <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/specific_suggestions_V1.csv", header = TRUE)

      # List of GO complexes that are included in CL definitions
      GO_complexes <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/GO_complexes_V1.csv", header = TRUE)

      ## List of CD synonym (adapted from list published by HCDM)
      CD_syn <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/Curated_synonyms_V1.csv", header = TRUE)

      # Auto replacements
      replace_markers <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/replace_markers_V1.csv", header=TRUE)

      ## Lists of non-specific suggestions
      # List of non-protein words that were found in markers from the deveolpment data
      common_non_protein_words <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/common_non_protein_words_V1.csv", header = FALSE)

      # List of metal tags that were found in the development data
      commonly_used_metal_tags <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/commonly_used_metal_tags_V1.csv", header = FALSE)

      # List of fluorophores tags that were  found in the development data
      commonly_used_fluorophores <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/commonly_used_fluorophores_V1.csv", header = FALSE)

      # Make markers uppercase and remove dashes, underscores, periods, and spaces
      new_df$Marker2 <- toupper(new_df$Marker)
      new_df$Marker2 <- gsub('-', '', new_df$Marker2, fixed = TRUE)
      new_df$Marker2 <- gsub('_', '', new_df$Marker2, fixed = TRUE)
      new_df$Marker2 <- gsub(' ', '', new_df$Marker2, fixed = TRUE)
      new_df$Marker2 <- gsub('.', '', new_df$Marker2, fixed = TRUE)

      # Make empty dataframe
      auto_matches <- data.frame(Marker = unique(new_df$Marker2))

      auto_matches$PRO_term <- ""
      auto_matches$PRO_name <- ""
      auto_matches$Species <- ""
      auto_matches$Match_type <- ""
      auto_matches$Match_step <- ""
      auto_matches$Suggestions <- ""

      ## Automatic matching (CD3e, CD8a, TCR)

      # Automatically make any inputs of CD3e, CD8a, and TCR match to the wider complex

          # If TCR is also specified in the panel, don't auto match CD3E to it
    if (any(auto_matches=="TCRAB") | any(auto_matches=="TCRA/B" ) | 
        any(auto_matches=="TCRABCOMPLEX") | any(auto_matches=="TCRA/BCOMPLEX" ) |
        any(auto_matches=="ABTCR") | any(auto_matches=="ABTCRCOMPLEX" ) |
        any(auto_matches=="ALPHABETATCR") | any(auto_matches=="ALPHABETATCRCOMPLEX" ) |
        any(auto_matches=="TCRALPHABETA") | any(auto_matches=="TCRALPHABETACOMPLEX" ) |
        any(auto_matches=="TCRA") | any(auto_matches=="TCRB" ) |
        any(auto_matches=="TCRALPHA") | any(auto_matches=="TCRBETA" ) |
        any(auto_matches=="ALPHATCR") | any(auto_matches=="BETATCR" ) |
        any(auto_matches=="TCRGD") | any(auto_matches=="TCRGDCOMPLEX" ) |
        any(auto_matches=="TCRG/D") | any(auto_matches=="TCRG/DCOMPLEX" ) |
        any(auto_matches=="GDTCR") | any(auto_matches=="GDTCRCOMPLEX" ) |
        any(auto_matches=="GAMMADELTATCR") | any(auto_matches=="GAMMADELTATCRCOMPLEX" ) |
        any(auto_matches=="TCRGAMMADELTA") | any(auto_matches=="TCRGAMMADELTACOMPLEX" ) |
        any(auto_matches=="TCRG") | any(auto_matches=="TCRD" ) |
        any(auto_matches=="TCRGAMMA") | any(auto_matches=="TCRDELTA" ) |
        any(auto_matches=="GAMMATCR") | any(auto_matches=="DELTATCR" ) |
        any(auto_matches=="YTCR") | any(auto_matches=="TCRY" ) |
        any(auto_matches=="TCRYD") | any(auto_matches=="YDTCR" ) |
        any(auto_matches=="TCR") | any(auto_matches=="TCRCOMPLEX" ) |
        any(auto_matches=="TCELLRECEPTOR") | any(auto_matches=="TCELLRECEPTORCOMPLEX" )) {
      
      replace_markers <- replace_markers[!grepl("TCRAB", replace_markers$Replacement),]
      replace_markers <- replace_markers[!grepl("TCRGD", replace_markers$Replacement),]
      replace_markers <- replace_markers[!grepl("TCR", replace_markers$Replacement),]
    }
    
    # If CD8AB is also specified in the panel, don't auto match CD8 to it
    if (any(auto_matches=="CD8AB") | any(auto_matches=="CD8ALPHABETA")) {
      replace_markers <- replace_markers[!grepl("CD8ALPHABETA", replace_markers$Replacement),]
    }

      auto_matches <- merge(replace_markers, auto_matches,  by = "Marker", all.y = TRUE)

      auto_replaced <- auto_matches[!is.na(auto_matches$Replacement),]

      auto_matches$Replacement[is.na(auto_matches$Replacement)] <- auto_matches$Marker[is.na(auto_matches$Replacement)]

      auto_matches$Marker <- NULL

      colnames(auto_matches)[1] <- "Marker"

      ## Punctuation

      ## Check if these punctuation marks are in the input.

      Punctuation_matches <- auto_matches

      # Slash
      Punctuation_matches$slash <- ifelse(grepl("/", Punctuation_matches$Marker, useBytes = TRUE), "A slash '/' was detected in the marker name.", "")

      # Backslash
      Punctuation_matches$backslash <- ifelse(grepl("\\\\", Punctuation_matches$Marker, useBytes = TRUE), "A backslash '\' was detected in the marker name.", "")

      # Comma
      Punctuation_matches$comma <- ifelse(grepl(",", Punctuation_matches$Marker, useBytes = TRUE), "A comma ',' was detected in the marker name.", "")

      # Bracket
      Punctuation_matches$bracket <- ifelse(grepl("[\\[\\]", Punctuation_matches$Marker, useBytes = TRUE), "A bracket '[' was detected in the marker name.", "")

      # Parenthesis
      Punctuation_matches$parenthesis <- ifelse(grepl("[\\(\\)]", Punctuation_matches$Marker, useBytes = TRUE), "A paranthesis '(' was detected in the marker name.", "")

      # Paste together as suggestions
      Punctuation_matches$Suggestions <- paste(Punctuation_matches$slash, Punctuation_matches$backslash, Punctuation_matches$comma, Punctuation_matches$bracket, Punctuation_matches$parenthesis)
      Punctuation_matches$Suggestions <- gsub("^ *|(?<= ) | *$", "", Punctuation_matches$Suggestions, perl = TRUE)

      # Remove puncationation columns
      Punctuation_matches <- subset(Punctuation_matches, select = -c(slash, comma, bracket, parenthesis, backslash))

      # Indicate that marker was flagged for punctuation to change
      for(i in 1:nrow(Punctuation_matches)) {
        if(Punctuation_matches$Suggestions[i] != "") {
          Punctuation_matches$Match_step[i] <- "Punctuation suggestions"
        }
      }

      ## Non-specific Suggestions

      # If marker name partial matches to an entry on one of these lists, give specific suggestions on how to rewrite that name

      # Vector of unique marker names
      Non_protein_sugg_matches <- Punctuation_matches[Punctuation_matches$Match_step != "Punctuation suggestions",]

      if(length(Non_protein_sugg_matches != 0)) {

        # List of words that are commonly used in non-protein marker names/fcs headers
        Non_protein_sugg_matches$Common_non_protein_words <- ifelse(grepl(paste0(common_non_protein_words$V1, collapse = "|"), Non_protein_sugg_matches$Marker, useBytes = TRUE), "A non-protein word was detected in the marker name.", "")

        # List of metals
        Non_protein_sugg_matches$commonly_used_metal_tags <- ifelse(grepl(paste0(commonly_used_metal_tags$V1, collapse = "|"), Non_protein_sugg_matches$Marker, useBytes = TRUE), "A metal tag was detected in the marker name.", "")

        # List of fluorophores
        Non_protein_sugg_matches$commonly_used_fluorophores <- ifelse(grepl(paste0(commonly_used_fluorophores$V1, collapse = "|"), Non_protein_sugg_matches$Marker, useBytes = TRUE), "A fluorophore was detected in the marker name.", "")

        # Put all suggestions together
        Non_protein_sugg_matches$Suggestions <- paste(Non_protein_sugg_matches$Common_non_protein_words, Non_protein_sugg_matches$commonly_used_metal_tags, Non_protein_sugg_matches$commonly_used_fluorophores, Non_protein_sugg_matches$species_specific_letter)
        Non_protein_sugg_matches$Suggestions <- gsub("^ *|(?<= ) | *$", "", Non_protein_sugg_matches$Suggestions, perl = TRUE)

        for(i in 1:nrow(Non_protein_sugg_matches)) {
          if(Non_protein_sugg_matches$Suggestions[i] != "") {
            Non_protein_sugg_matches$Match_step[i] <- "Non-protein word suggestions"
          }
        }

        Non_protein_sugg_matches <- subset(Non_protein_sugg_matches, select = -c(Common_non_protein_words, commonly_used_metal_tags, commonly_used_fluorophores))
      }

      ## GO Complexes

      # Match entries to the GO complexes listed in the CL

      # Vector of unique marker names
      GO_matches <- Non_protein_sugg_matches[Non_protein_sugg_matches$Match_step != "Non-protein word suggestions",]

      unique_marker_names <- GO_matches$Marker

      all_ids <- list()

      for(marker2 in unique_marker_names) {

        matched_row <- GO_complexes[GO_complexes$Inputted_Name ==marker2,]

        if (nrow(matched_row) != 0) {

          GO_matches[GO_matches$Marker == marker2,]$Suggestions <- ""
          GO_matches[GO_matches$Marker == marker2,]$Match_step <- "GO exceptions"
          GO_matches[GO_matches$Marker == marker2,]$PRO_term <- matched_row$Matched_GO_term
          GO_matches[GO_matches$Marker == marker2,]$PRO_name <- matched_row$Matched_Descriptive_GO_name
        }
      }

      ## Direct SPAPRQL query to PRO

      # Directly query PRO

      Match_step <- "Directly to PRO"

      unique_marker_names <- GO_matches$Marker[GO_matches$Match_step != "GO exceptions"]

      # Run SPARQL
      all_ids <- PRO_SPARQL(unique_marker_names = unique_marker_names,
                            NCBI_taxon_ID = NCBI_taxon_ID_short,
                            onto_endpoint = onto_endpoint,
                            Match_step = Match_step)

      # Put all results in one data frame
      Direct_PRO_matches <- dplyr::bind_rows(all_ids, .id = "Marker")

      ## Specific Suggestions

      # If marker name matches to an entry on this list, give specific suggestions on how to rewrite that name

      # Vector of unique marker names
      Specific_protein_sugg_matches <- Direct_PRO_matches[Direct_PRO_matches$Match_step != "Directly to PRO",]

      unique_marker_names <- Specific_protein_sugg_matches$Marker

      all_ids <- list()

      if(length(unique_marker_names) != 0) {

        Specific_protein_sugg_matches$Suggestions <- ""

        for(marker2 in unique_marker_names) {

          matched_row <- marker_alt[marker_alt$CD_NAME == marker2,]

          if (nrow(matched_row) != 0) {

            remove_empties <- sapply(matched_row, function(x) all(is.na(x) | x == ""))
            matched_row <- matched_row[, !remove_empties]

            suggestions <- paste(matched_row[, -1], collapse = ", ")

            Specific_protein_sugg_matches[Specific_protein_sugg_matches$Marker == marker2,]$Suggestions <- paste("Try:", suggestions)
            Specific_protein_sugg_matches[Specific_protein_sugg_matches$Marker == marker2,]$Match_step <- "Specific protein suggestions"
            Specific_protein_sugg_matches[Specific_protein_sugg_matches$Marker == marker2,]$PRO_term <- ""
            Specific_protein_sugg_matches[Specific_protein_sugg_matches$Marker == marker2,]$PRO_name <- ""
          }
        }
      } else {
        Specific_protein_sugg_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
      }


      ## CD synonym list, then SPAPRQL query to PRO

      # Search CD synonym list, then if there are synonyms query them in PRO

      Match_step <- "CD synonym list"

      # Vector of unique marker names
      unique_marker_names <- Specific_protein_sugg_matches$Marker[Specific_protein_sugg_matches$Match_step != "Specific protein suggestions"]

      all_ids_final <- list()

      if(length(unique_marker_names) != 0) {

        for (marker in unique_marker_names) {

          matched_row <- CD_syn[which(CD_syn == marker, arr.ind=TRUE)[,1],]

          matched_row <- data.frame(x=unlist(matched_row))

          unique_syn_marker_names <- unique(stats::na.omit(matched_row))

          unique_syn_marker_names <- as.vector(unique_syn_marker_names[,1])

          # Run SPARQL
          all_ids <- PRO_SPARQL(unique_marker_names = unique_syn_marker_names,
                                NCBI_taxon_ID = NCBI_taxon_ID_short,
                                onto_endpoint = onto_endpoint,
                                Match_step = Match_step)

          # Put all results in one data frame
          all_ids <- dplyr::bind_rows(all_ids, .id = "Matched_synonym")

          all_ids <- stats::na.omit(all_ids)

          if (length(all_ids) != 0) {
            all_ids_final[[marker]] <- all_ids
          } else if (length(all_ids) == 0) {
            all_ids_final[[marker]] <- data.frame(Matched_synonym = "", PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
          }
        }

        # Put all results in one data frame
        all_ids_final2 <- dplyr::bind_rows(all_ids_final, .id = "Marker")

        # Remove duplicated entries (Marker, PRO term, and type of match are the same)
        all_ids_final2 <- all_ids_final2[!duplicated(all_ids_final2[,c("Marker","PRO_term", "Match_type")]),]

        # Get entries were inputted marker and matched species is the same

        # Subset into entries with and without duplicates
        x <-  all_ids_final2[duplicated(all_ids_final2[,c("Marker","Species")]) | duplicated(all_ids_final2[,c("Marker","Species")], fromLast = TRUE) ,]
        y <-  dplyr::setdiff(all_ids_final2, x)

        # Get vector of markers (duplicates)
        unique_marker_names <- unique(x$Marker)

        # If there is both an exact match and a secondary match, remove the secondary
        remove_secondary_per_marker <- data.frame()

        for(marker in unique_marker_names) {
          df_subsetted <- x[x$Marker == marker,]
          if(any(df_subsetted$Match_type == "Primary")) {
            df_per_marker <- df_subsetted[df_subsetted$Match_type != "Secondary",]
            df_per_marker <- df_per_marker[df_per_marker$Match_type != "",]
            remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
          } else if(any(df_subsetted$Match_type == "Secondary")) {
            df_per_marker <- df_subsetted[df_subsetted$Match_type != "",]
            remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
          }
        }

        # Final
        CD_syn_list_matches <- rbind(remove_secondary_per_marker, y)

        # Just columns we want
        CD_syn_list_matches <- unique(CD_syn_list_matches[ , c("Marker", "PRO_term", "PRO_name", "Species", "Match_type", "Match_step")])

      } else {
        CD_syn_list_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
      }

      ## Wikidata, then SPAPRQL query to PRO

      # Query Wikidata, then if there are synonyms query them in PRO

      Match_step <- "Wikidata"

      # Vector of unique marker names
      unique_marker_names <- CD_syn_list_matches$Marker[CD_syn_list_matches$Match_step != "CD synonym list"]

      # Remove duplicates introduces by synonyms from last step
      remove <- CD_syn_list_matches$Marker[CD_syn_list_matches$Match_step != ""]
      unique_marker_names <-  dplyr::setdiff(unique_marker_names, remove)

      all_ids_final <- list()

      if(length(unique_marker_names) != 0){
        if(unique_marker_names[1] != "") {

          # Loop through each and do the SPARQL query
          for(marker in unique_marker_names) {

            specific_marker <-  shQuote(marker)

            SPARQL_query <-
              paste0(
                "
SELECT DISTINCT ?item ?itemLabel ?altLabel WHERE {
  # Look either for non-species specific or add a specific species
  OPTIONAL {
    ?item wdt:P703 <", wikidata_taxon_ID, "> .

  }
  # Get Wikidata protein name and alternative names
  ?item rdfs:label ?itemLabel .
  ?item skos:altLabel ?altLabel.
  # Return results that are either proteins or protein-coding genes
  {?item wdt:P31 wd:Q8054}
  UNION
  {?item wdt:P279 wd:Q20747295}
  # Search for the marker input, either if it is the entry name or an alternative name
  FILTER (UCASE(REPLACE(str(?itemLabel),'[ -.]','')) = ", specific_marker, " || UCASE(REPLACE(str(?altLabel),'[ -.]','')) = ", specific_marker, ")
}
")

            # Run SPARQL on the endpoint
            result <- httr::GET(
              url = wiki_endpoint,
              query = list(query = SPARQL_query),
              httr::user_agent(R.version.string))


            # Will show a warning/error if there is any
            httr::stop_for_status(result)

            # Get result in text JSON
            x <- httr::content(result, as = "text") #, encoding = "UTF-8")

            # Convert from JSON to a list
            df <- jsonlite::fromJSON(x, flatten = TRUE)

            # Extract the data frame
            df <- df$results$bindings

            if (length(df) != 0) {

              # Remove unneeded info
              df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", 'lang')))

              # Remove this part that was added onto the column names
              colnames(df) <- gsub(".value", "", colnames(df))

              # Make everything case, remove spaces, dashes
              df <- dplyr::mutate_all(df, .funs = toupper)
              df <- as.data.frame(lapply(df, function(y) gsub('-', '', y)))
              df <- as.data.frame(lapply(df, function(y) gsub('_', '', y)))
              df <- as.data.frame(lapply(df, function(y) gsub(' ', '', y)))
              df <- as.data.frame(lapply(df, function(y) gsub('\\.', '', y)))

              df <- df[(which(nchar(df$altLabel) > 2)),]
              df <- df[(which(nchar(df$itemLabel) > 2)),]

              # No marker can start with a number, so add X to any marker names that start with a number
              df <- as.data.frame(lapply(df, function(y) gsub('^(\\d)', 'X\\1', y)))

              # Vector of unique marker names
              unique_syn_names <- unique(c(df[, "altLabel"], df[, "itemLabel"]))

              unique_syn_names <- unique_syn_names[!grepl("[^\x01-\x7F]+", unique_syn_names)]
                                         
              unique_syn_names <- unique_syn_names[grep('[A-Za-z]',unique_syn_names)]

              unique_syn_names <- unique_syn_names[!grepl(paste0('^', marker, '$'), unique_syn_names, useBytes = TRUE)]

              # Run SPARQL
              all_ids <- PRO_SPARQL(unique_marker_names = unique_syn_names,
                                    NCBI_taxon_ID = NCBI_taxon_ID_short,
                                    onto_endpoint = onto_endpoint,
                                    Match_step = Match_step)

              # Put all results in one data frame
              all_ids <- dplyr::bind_rows(all_ids, .id = "Matched synonym")

              all_ids <- stats::na.omit(all_ids)

              if (length(all_ids) != 0) {
                all_ids_final[[marker]] <- all_ids
              } else if (length(all_ids) == 0) {
                all_ids_final[[marker]] <- data.frame(PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
              }
            } else if (length(df) == 0) {
              all_ids_final[[marker]] <- data.frame(PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
            }
          }

        # Put all results in one data frame
        all_ids_final2 <- dplyr::bind_rows(all_ids_final, .id = "Marker")

        # Remove duplicated entries (Marker, PRO term, and type of match are the same)
        all_ids_final2 <- all_ids_final2[!duplicated(all_ids_final2[,c("Marker","PRO_term", "Match_type")]),]

        # Get entries were inputted marker and matched species is the same

        # Subset into entries with and without duplicates
        x <-  all_ids_final2[duplicated(all_ids_final2[,c("Marker","Species")]) | duplicated(all_ids_final2[,c("Marker","Species")], fromLast = TRUE) ,]
        y <-  dplyr::setdiff(all_ids_final2, x)

        # Get vector of markers (duplicates)
        unique_marker_names <- unique(x$Marker)

        # If there is both an exact match and a secondary match, remove the secondary
        remove_secondary_per_marker <- data.frame()

        for(marker in unique_marker_names) {
          df_subsetted <- x[x$Marker == marker,]
          if(any(df_subsetted$Match_type == "Primary")) {
            df_per_marker <- df_subsetted[df_subsetted$Match_type != "Secondary",]
            df_per_marker <- df_per_marker[df_per_marker$Match_type != "",]
            remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
          } else if(any(df_subsetted$Match_type == "Secondary")) {
            df_per_marker <- df_subsetted[df_subsetted$Match_type != "",]
            remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
          }
        }

        # Final
        Wikidata_matches <- rbind(remove_secondary_per_marker, y)

        # Just columns we want
        Wikidata_matches <- unique(Wikidata_matches[ , c("Marker", "PRO_term", "PRO_name", "Species", "Match_type", "Match_step")])

        } else {
          Wikidata_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
        }
      } else {
        Wikidata_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
      }

    ## Final results

    # Put everything together and get the final list of matched/unmatched marker names

    # Add empty suggestion column to parts where suggestions were not relevant
    Direct_PRO_matches$Suggestions <- ""
    CD_syn_list_matches$Suggestions <- ""
    Wikidata_matches$Suggestions <- ""

    # Subset entries that were matched at each step
    Punctuation_matches2 <- Punctuation_matches[Punctuation_matches$Match_step == "Punctuation suggestions",]
    Non_protein_sugg_matches2 <- Non_protein_sugg_matches[Non_protein_sugg_matches$Match_step == "Non-protein word suggestions",]
    GO_matches2 <- GO_matches[GO_matches$Match_step == "GO exceptions",]
    Direct_PRO_matches2 <- Direct_PRO_matches[Direct_PRO_matches$Match_step == "Directly to PRO",]
    CD_syn_list_matches2 <- CD_syn_list_matches[CD_syn_list_matches$Match_step == "CD synonym list",]
    Wikidata_matches2 <- Wikidata_matches[Wikidata_matches$Match_step == "Wikidata",]
    Specific_protein_sugg_matches2 <- Specific_protein_sugg_matches[Specific_protein_sugg_matches$Match_step == "Specific protein suggestions",]

    # Put everything together into 1 dataframe
    all_ids <- do.call("rbind", list(Punctuation_matches2, Non_protein_sugg_matches2, GO_matches2,  Direct_PRO_matches2, CD_syn_list_matches2, Wikidata_matches2, Specific_protein_sugg_matches2))

    # Put original inputted names back in
    final_df <- merge(new_df, all_ids, by.x = "Marker2", by.y = "Marker", all.x = TRUE, all.y = TRUE)

    # Change names
    colnames(final_df)[1:2] <- c("New_Name", "Original_Name")

    if (nrow(auto_replaced) != 0) {
      # Put back original auto changed names
      replacements <- auto_replaced[as.character(auto_replaced$Marker) != as.character(auto_replaced$Replacement),]

      replacements <- replacements[,1:2]

      colnames(replacements) <- c("Marker2", "New_Name")

      replacements <- merge(replacements, new_df, by = "Marker2")

      replacements <- replacements[,c(2,3)]

      colnames(replacements) <- c("New_Name", "Original_Name")

      replacements <- merge(replacements, final_df, by = "New_Name")

      replacements <- replacements[,-3]

      colnames(replacements)[2] <- "Original_Name"

      final_df <- rbind(final_df, replacements)
    }

    final_df <- final_df[!is.na(final_df$Original_Name), ]

    # If there was no match or suggestion, indicate that
    final_df$Match_step[is.na(final_df$Match_step)] <- "Did not match"

    # Turn any remaining NA entries into blanks
    final_df[is.na(final_df)] <- ""

    # Change species URIs into human readable form
    final_df$Species <- gsub(NCBI_taxon_ID_full, input$species_choice_2, final_df$Species)
    keep <- c(input$species_choice_2, "Not species specific")

    # Remove IDs that aren't actual species
    remove_non_species <- dplyr::setdiff(final_df$Species, keep)
    remove_non_species <- remove_non_species[remove_non_species != ""]

    if (length(remove_non_species) != 0) {
      # Fill these in as species unspecific
      final_df$Species <- stringi::stri_replace_all_regex(final_df$Species,
                                                          pattern=remove_non_species,
                                                          replacement="Not species specific",
                                                          vectorize=FALSE)
    }

    # Remove duplicates
    final_df <- final_df[!duplicated(final_df), ]

    # Order rows by workflow step
    new_order <- c("Punctuation suggestions", "Non-protein word suggestions", "GO exceptions", "Directly to PRO", "Specific protein suggestions", "CD synonym list", "Wikidata")
    final_df <- final_df[order(base::match(final_df$Match_step, new_order)),]

    # Get new column names
    colnames(final_df) <-
      c("Revised Name",
        "Original Name",
        "Sign",
        "Cluster",
        "PRO/GO Term",
        "PRO/GO Name",
        "Species Specific",
        "Match type",
        "Workflow Step",
        "Suggestion")

    # Change column order
    final_df <- final_df[, c("Original Name",	"Revised Name", "Sign", "Cluster", "Workflow Step",	"PRO/GO Name",	"PRO/GO Term",	"Species Specific",	"Match type",	"Suggestion")]

    return(final_df)
    } else {
      return(NULL)
    }
  })

  # Step 1, get statement indicating if/which markers were not able to be matched
  return_unmatched_markers_2 <- reactive({
    shiny::req(input$reference_type_2 == 'internal_ref_2')
    shiny::req(input$submit_tab2_step1)

    if (is.data.frame(reformatted_data_2()) == TRUE) {
      all_ids <- PRO_results_2()

      # Get markers that did not match
      unmatched_markers <- all_ids[all_ids$`PRO/GO Term` == "",]

      # Get just the original inputted marker name
      unmatched_markers <- unmatched_markers$`Original Name`

      # Just unique markers
      unmatched_markers <- unique(unmatched_markers)

      # Return nothing if all markers had located IDs
      if (length(unmatched_markers) == 0) {
        return(NULL)
      } else {
        # Return markers that were not able to be matched
        return(unmatched_markers)
      }
    }
  })

  # Step 1, output message about if markers matched to PRO term
  output$ui_warning_unmatched_markers_2 <- shiny::renderUI({
    shiny::req(input$reference_type_2 == 'internal_ref_2')
    shiny::req(input$submit_tab2_step1)

    if (is.data.frame(reformatted_data_2()) == TRUE) {
      if (length(marker_sign_error_2()) == 0) {
        # Return nothing if all markers had located IDs
        if (is.null(return_unmatched_markers_2())) {
          h4(
            paste(
              "All inputted markers matched to at least 1 protein ID. Please continue to the next step."
            )
          )
        } else if (length(return_unmatched_markers_2()) > 1) {
          marker_string <- paste(return_unmatched_markers_2(), collapse = ", ")
          h4(
            paste(
              "Protein IDs for the inputted markers",
              marker_string,
              "were not found. Please input new marker names."
            )
          )
        } else if (length(return_unmatched_markers_2()) == 1) {
          h4(
            paste(
              "A Protein ID for the inputted marker",
              return_unmatched_markers_2(),
              "was not found. Please input a new marker name."
            )
          )
        }
      }
    }
  })

  # Step 1, option to edit (add new marker names) in the table
  shiny::observeEvent(input$table_new_format_2_cell_edit, {
    #  shiny::req(input$reference_type_2 == 'internal_ref_2')
    df_1_2$data <<-
      DT::editData(df_1_2$data,
                   input$table_new_format_2_cell_edit,
                   'table_new_format_2',
                   rownames = FALSE)
  })

  ######~#~# Step 2 - Get marker synonyms of markers with an edited name, user can choose which matched marker terms to procees with #~#~######

  # Step 2, side panel

  # Reset everything
  shiny::observeEvent(input$reset_tab2_step2, {
    session$reload()
  })

  shiny::observeEvent(input$help_low_2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select the Interpretation of 'Low'"),
      HTML("Specify how you would like the 'low' modifier to be interpreted. This is used to determine how inputted markers will be matched to markers within the reference cell types. Low can be interpreted as either:
      <br>
      <br>
      <ul>
      <li>'low'</li>
      <br>
      <li>'low' and 'positive'</li>
      <br>
      <li>'low' and 'negative'</li>"),
      easyClose = TRUE))
  })

  shiny::observeEvent(input$help_high_2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Select the Interpretation of 'High'"),
      HTML("Specify how you would like the 'high' modifier to be interpreted. This is used to determine how inputted markers will be matched to markers within the reference cell types. High can be interpreted as either:
      <br>
      <br>
      <ul>
      <li>'high'</li>
      <br>
      <li>'high' and 'positive'</li>"),
      easyClose = TRUE))
  })


  # Step 2, main panel

  # Step 2, create editable dataframes
  df2_2 <- shiny::reactiveValues(data = NULL)

  # Step 2, when help button is clicked show message
  shiny::observeEvent(input$help_tab2_step2,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 2<br>Confirm Matched Marker Terms"),
      HTML("This step displays the results of the marker standardization process. The data frame shows:
      <br>
      <br>
      <li>Original marker names</li>
      <br>
      <li>Revised marker names (if applicable)</li>
      <br>
      <li>Matched Protein Ontology (PRO) or Gene Ontology (GO) terms, which are clickable for direct access to the ontology page</li>
      <br>
      <li>Reference, which indicates whether the term is included in the Cell Ontology</li>
      <br>
      If markers match to multiple ontology terms, the tool will list all possible matches. This can occur if an inputted name is unspecific and listed as a synonym for multiple different proteins. Additionally, an inputted marker may match to both a species specific and species non-specific identifier. Any ontology term that begins with a 'P' (PR_P???) is species specific. Oftentimes specific-specific terms are not included within the Cell Ontology, but when they are it is suggested that they remain included.
      <br>
      <br>
      Markers not found in the Cell Ontology are excluded from subsequent steps. Additionally, any markers can be manually removed by selecting the row(s) and clicking 'Delete Marker Entry'. Deleting markers is irreversible.
      <br>
      <br>
      Users can download the matched terms in CSV format using the 'Download' button.
      <br>
      <br>"),
      tags$div("Additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/4.1.-Input-Marker-Descriptors-and-Default-References', "here.", target="_blank")),
      easyClose = TRUE))
  })

  # Step 2, changes when the 'save and continue' button is hit if using internal references
  # Make editable dataframe, show buttons once that dataframe is made
  shiny::observeEvent(input$submit_tab2_step2, {

    shiny.destroy::removeOutput("help_tab2_step1")
    shiny.destroy::removeOutput("ui_warning_unmatched_markers_2")
    shiny.destroy::removeOutput("overview_tab2_step1")
    shinyjs::show("overview_tab2_step2")
    shiny.destroy::removeOutput("table_new_format_2")

    df2_2$data <- return_protein_synonym_2()

    # Step 2, make datatable of markers and their matched PRO terms
    output$table_protein_synonym_2 <- DT::renderDT({
      if (is.null(df2_2$data)) {
      } else {
        DT::datatable(
          df2_2$data,
          rownames = FALSE,
          escape = FALSE,
          options = list(pageLength = 25)
        )
      }
    }, server=FALSE)

    if (is.data.frame(return_protein_synonym_2()) == TRUE) {
      shinyjs::show("help_tab2_step2")
      shinyjs::show("delete_PRO_synonym_2")
      shinyjs::show("download_PRO_GO_matches_2")
      shiny.destroy::removeOutput("all_sidebar_tab2_step2")
      shinyjs::show("all_sidebar_tab2_step3")
    }
  })

  # Step 2, match inputted proteins to entries in the Protein Ontology, get PRO IDs
  # Make PRO link clickable
  # Merge in info about whether or not marker is in the reference(s)
  return_protein_synonym_2 <- shiny::eventReactive(input$submit_tab2_step2, {

    #   shiny::req(input$reference_type_2 == 'internal_ref_2')
    # Original results
    org_pro_results <- PRO_results_2()

    # Remove no matches found from original PRO query
    org_pro_results <- org_pro_results[org_pro_results$`PRO/GO Term` != "",]

    # Get certain columns
    org_pro_results <- org_pro_results[c("Original Name", "Revised Name", "PRO/GO Term")]

    # Remove duplicates
    org_pro_results <- org_pro_results[!duplicated(org_pro_results), ]

    if (!is.null(df_1_2$data)) {
      # Redo results
      new_pro_results <- new_PRO_matches_2()

      # Get certain columns
      new_pro_results <- new_pro_results[c("Original Name", "Revised Name", "PRO/GO Term")]

      # Remove duplicates
      new_pro_results <- new_pro_results[!duplicated(new_pro_results), ]

      # Bind all PRO results together
      all_pro_results <- rbind(org_pro_results, new_pro_results)

    } else {
      all_pro_results <- org_pro_results
    }

    # If there was no match or suggestion, indicate that
    all_pro_results$`PRO/GO Term`[all_pro_results$`PRO/GO Term` == ""] <- "Did not match"

    df_ref <- markers_in_ref_2()

    # Get inputted markers that are not in the reference
    not_in_reference <- setdiff(all_pro_results$`PRO/GO Term`, df_ref$Related_protein)

    if (length(not_in_reference) != 0) {
      df <- data.frame("PRO/GO Term" = not_in_reference, "Reference" = "Not in reference. This PRO term will be removed from subsequent analysis.", check.names = FALSE)

      all_pro_results <-  merge(all_pro_results,
                                df, by = "PRO/GO Term", all = TRUE)

      all_pro_results$Reference[is.na(all_pro_results$Reference)] <- "Included in reference."

    } else {
      all_pro_results$Reference <- "Included in reference."

    }

    # Change column order
    all_pro_results <- all_pro_results[, c("Original Name",	"Revised Name",	"PRO/GO Term", "Reference")]

    # Order by marker name
    all_pro_results <-
      all_pro_results[order(all_pro_results$`Original Name`, decreasing = FALSE), ]

    # Just get PRO ID to make link in R Shiny look better
    just_IDs <-
      gsub("http://purl.obolibrary.org/obo/",
           "",
           all_pro_results$`PRO/GO Term`)

    # Remove '<' and '>' from URI
    all_pro_results$`PRO/GO Term` <-
      gsub('<', '', all_pro_results$`PRO/GO Term`)
    all_pro_results$`PRO/GO Term` <-
      gsub('>', '', all_pro_results$`PRO/GO Term`)
    just_IDs <- gsub('>', '', just_IDs)
    just_IDs <- gsub('<', '', just_IDs)

    # Make it so the link is clickable
    for (a in 1:nrow(all_pro_results)) {
      if (all_pro_results$`PRO/GO Term`[a] != "Did not match") {
        all_pro_results$`PRO/GO Term`[a] <-
          paste0(
            "<a href='",
            all_pro_results$`PRO/GO Term`[a],
            "' target='_blank'>",
            just_IDs[a],
            "</a>"
          )
      }
    }

    return(all_pro_results)
  })

  # Step 2, see if inputed marker is in reference
  markers_in_ref_2 <- reactive({
    shiny::req(input$submit_tab2_step2)

    if (length(input$ontology_type_2) == 1) {
      if (input$ontology_type_2 == 'internal_CL_2') {
        ontology_URI <- "<http://purl.obolibrary.org/obo/merged/CL>"
      } else if (input$ontology_type_2 == 'internal_pCL_2') {
        ontology_URI <- "<http://purl.obolibrary.org/obo/merged/PCL>"
      }
    } else if (length(input$ontology_type_2) == 2) {
      ontology_URI <- "<http://purl.obolibrary.org/obo/merged/CL> FROM <http://purl.obolibrary.org/obo/merged/PCL>"
    }

    relationship_type1 <- "<http://purl.obolibrary.org/obo/RO_0002104>"
    relationship_type2 <- "<http://purl.obolibrary.org/obo/RO_0015015>"
    relationship_type3 <- "<http://purl.obolibrary.org/obo/RO_0015016>"
    relationship_type4 <- "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
    relationship_type5 <- "<http://purl.obolibrary.org/obo/CL_4030046>"
    relationship_type6 <- "<http://purl.obolibrary.org/obo/BFO_0000051>"

    # SPARQL query
    query <-

      paste0(
        "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?Related_protein ?Protein_name ?Relationship_type

  FROM ",
        ontology_URI,
        "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }

  FILTER (?Relationship_type =", relationship_type1, "|| ?Relationship_type =", relationship_type2, "|| ?Relationship_type =", relationship_type3, "|| ?Relationship_type =", relationship_type4, "|| ?Relationship_type =", relationship_type5, "|| ?Relationship_type =", relationship_type6, ")

    }
  "
      )

    # Run SPARQL on the endpoint
    result <- httr::POST(onto_endpoint,
                         body = list(query = query),
                         httr::user_agent(R.version.string))

    # Will show a warning/error if there is any
    httr::stop_for_status(result)

    # Get result in text JSON
    x <- httr::content(result, "text", encoding = "UTF-8")

    # Convert from JSON to a list
    df <- jsonlite::fromJSON(x, flatten = TRUE)

    # Extract the dataframe
    df <- df$results$bindings

    # Put the results in a list of data frames
    if (length(df) != 0) {

      # Remove unneeded info
      df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

      # Remove this part that was added onto the column names
      colnames(df) <- gsub(".value", "", colnames(df))
    }

    # Only get relevant relationship types for PRO/GO terms
    just_PRO <- subset(df, startsWith(as.character(Related_protein),"http://purl.obolibrary.org/obo/PR_"))
    just_GO <- subset(df, startsWith(as.character(Related_protein),"http://purl.obolibrary.org/obo/GO_") & Relationship_type != "http://purl.obolibrary.org/obo/BFO_0000051")

    # Combine
    df <- rbind(just_PRO, just_GO)

    # Remove protein identifier
    df <- df[!df$Related_protein == "http://purl.obolibrary.org/obo/PR_000000001", ]

    # Only keep first 2 columns
    df <- df[ , c("Related_protein","Protein_name")]

    # Remove duplicates
    df <- df[!duplicated(df), ]

    return(df)
  })

  # Step 2, look up PRO IDs for the revised marker names
  new_PRO_matches_2 <- shiny::eventReactive(input$submit_tab2_step2, {

    ######## Species, references, initial dataframe ########

    new_df <- df_1_2$data

    # If empty, just fill with inputted value
    new_df$`Type new marker name`[new_df$`Type new marker name` == ""] = new_df$`Original Name`[new_df$`Type new marker name` == ""]

    # Get NCBI taxon ID
    species_ID <- input$species_choice_2
    matched_row <-
      which(species_in_PRO == species_ID, arr.ind = TRUE)[1]
    NCBI_taxon_ID <-
      species_in_PRO[matched_row, ncol(species_in_PRO)]

    NCBI_taxon_ID_full <- species_in_PRO[matched_row, 1]
    NCBI_taxon_ID_short <- gsub(".*/","",NCBI_taxon_ID_full)

    # Get Wikidata taxon ID
    SPARQL_query <-

      paste0(
        "SELECT ?item ?itemLabel WHERE {

     ?item wdt:P2888 <", NCBI_taxon_ID_full, ">.

    SERVICE wikibase:label { bd:serviceParam wikibase:language 'en'. }
    }"
      )

    # Run SPARQL on the endpoint
    result <- httr::GET(
      url = wiki_endpoint,
      query = list(query = SPARQL_query),
      httr::user_agent(R.version.string))

    # Will show a warning/error if there is any
    httr::stop_for_status(result)

    # Get result in text JSON
    x <- httr::content(result, as = "text") #, encoding = "UTF-8")

    # Convert from JSON to a list
    df <- jsonlite::fromJSON(x, flatten = TRUE)

    # Extract the data frame
    df <- df$results$bindings

    if (length(df) != 0) {

      # Remove unneeded info
      df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", ':lang')))

      # Remove this part that was added onto the column names
      colnames(df) <- gsub(".value", "", colnames(df))
    }

    wikidata_taxon_ID <- df$item

    # List of specific protein suggestions
    marker_alt <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/specific_suggestions_V1.csv", header = TRUE)

    # List of GO complexes that are included in CL definitions
    GO_complexes <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/GO_complexes_V1.csv", header = TRUE)

    ## List of CD synonym (adapted from list published by HCDM)
    CD_syn <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/Curated_synonyms_V1.csv", header = TRUE)

    # Make markers uppercase and remove dashes, underscores, periods, and spaces
    new_df$Marker2 <- toupper(new_df$`Type new marker name`)
    new_df$Marker2 <- gsub('-', '', new_df$Marker2, fixed = TRUE)
    new_df$Marker2 <- gsub('_', '', new_df$Marker2, fixed = TRUE)
    new_df$Marker2 <- gsub(' ', '', new_df$Marker2, fixed = TRUE)
    new_df$Marker2 <- gsub('.', '', new_df$Marker2, fixed = TRUE)

    # Just keep old marker name and new marker name
    new_df <- new_df[,c("Original Name", "Marker2")]

    # Make empty dataframe
    auto_matches <- data.frame(Marker = unique(new_df$Marker2))

    auto_matches$PRO_term <- ""
    auto_matches$PRO_name <- ""
    auto_matches$Species <- ""
    auto_matches$Match_type <- ""
    auto_matches$Match_step <- ""
    auto_matches$Suggestions <- ""

    ######## Automatic matching (CD3e, CD8a, TCR) ########

    # Automatically make any inputs of CD3e, CD8a, and TCR match to the wider complex

    # Auto replacements
    replace_markers <- utils::read.csv("https://raw.githubusercontent.com/AmandaRT18/Cell.Naming/refs/heads/main/data/replace_markers_V1.csv", header=TRUE)

            # If TCR is also specified in the panel, don't auto match CD3E to it
    if (any(auto_matches=="TCRAB") | any(auto_matches=="TCRA/B" ) | 
        any(auto_matches=="TCRABCOMPLEX") | any(auto_matches=="TCRA/BCOMPLEX" ) |
        any(auto_matches=="ABTCR") | any(auto_matches=="ABTCRCOMPLEX" ) |
        any(auto_matches=="ALPHABETATCR") | any(auto_matches=="ALPHABETATCRCOMPLEX" ) |
        any(auto_matches=="TCRALPHABETA") | any(auto_matches=="TCRALPHABETACOMPLEX" ) |
        any(auto_matches=="TCRA") | any(auto_matches=="TCRB" ) |
        any(auto_matches=="TCRALPHA") | any(auto_matches=="TCRBETA" ) |
        any(auto_matches=="ALPHATCR") | any(auto_matches=="BETATCR" ) |
        any(auto_matches=="TCRGD") | any(auto_matches=="TCRGDCOMPLEX" ) |
        any(auto_matches=="TCRG/D") | any(auto_matches=="TCRG/DCOMPLEX" ) |
        any(auto_matches=="GDTCR") | any(auto_matches=="GDTCRCOMPLEX" ) |
        any(auto_matches=="GAMMADELTATCR") | any(auto_matches=="GAMMADELTATCRCOMPLEX" ) |
        any(auto_matches=="TCRGAMMADELTA") | any(auto_matches=="TCRGAMMADELTACOMPLEX" ) |
        any(auto_matches=="TCRG") | any(auto_matches=="TCRD" ) |
        any(auto_matches=="TCRGAMMA") | any(auto_matches=="TCRDELTA" ) |
        any(auto_matches=="GAMMATCR") | any(auto_matches=="DELTATCR" ) |
        any(auto_matches=="YTCR") | any(auto_matches=="TCRY" ) |
        any(auto_matches=="TCRYD") | any(auto_matches=="YDTCR" ) |
        any(auto_matches=="TCR") | any(auto_matches=="TCRCOMPLEX" ) |
        any(auto_matches=="TCELLRECEPTOR") | any(auto_matches=="TCELLRECEPTORCOMPLEX" )) {
      
      replace_markers <- replace_markers[!grepl("TCRAB", replace_markers$Replacement),]
      replace_markers <- replace_markers[!grepl("TCRGD", replace_markers$Replacement),]
      replace_markers <- replace_markers[!grepl("TCR", replace_markers$Replacement),]
    }
    
    # If CD8AB is also specified in the panel, don't auto match CD8 to it
    if (any(auto_matches=="CD8AB") | any(auto_matches=="CD8ALPHABETA")) {
      replace_markers <- replace_markers[!grepl("CD8ALPHABETA", replace_markers$Replacement),]
    }
   
    auto_matches <- merge(replace_markers, auto_matches,  by = "Marker", all.y = TRUE)

    auto_replaced <- auto_matches[!is.na(auto_matches$Replacement),]

    auto_matches$Replacement[is.na(auto_matches$Replacement)] <- auto_matches$Marker[is.na(auto_matches$Replacement)]

    auto_matches$Marker <- NULL

    colnames(auto_matches)[1] <- "Marker"

    ######## GO Complexes ########

    # Match entries to the GO complexes listed in the CL

    # Vector of unique marker names
    GO_matches <- auto_matches

    unique_marker_names <- GO_matches$Marker

    all_ids <- list()

    for(marker2 in unique_marker_names) {

      matched_row <- GO_complexes[GO_complexes$Inputted_Name ==marker2,]

      if (nrow(matched_row) != 0) {

        GO_matches[GO_matches$Marker == marker2,]$Suggestions <- ""
        GO_matches[GO_matches$Marker == marker2,]$Match_step <- "GO exceptions"
        GO_matches[GO_matches$Marker == marker2,]$PRO_term <- matched_row$Matched_GO_term
        GO_matches[GO_matches$Marker == marker2,]$PRO_name <- matched_row$Matched_Descriptive_GO_name
      }
    }

    ######## Direct SPAPRQL query to PRO ########

    # Directly query PRO

    Match_step <- "Directly to PRO"

    unique_marker_names <- GO_matches$Marker[GO_matches$Match_step != "GO exceptions"]

    # Run SPARQL
    all_ids <- PRO_SPARQL(unique_marker_names = unique_marker_names,
                          NCBI_taxon_ID = NCBI_taxon_ID_short,
                          onto_endpoint = onto_endpoint,
                          Match_step = Match_step)

    # Put all results in one data frame
    Direct_PRO_matches <- dplyr::bind_rows(all_ids, .id = "Marker")

    ######## CD synonym list, then SPAPRQL query to PRO ########

    # Search CD synonym list, then if there are synonyms query them in PRO

    Match_step <- "CD synonym list"

    # Vector of unique marker names
    unique_marker_names <- Direct_PRO_matches$Marker[Direct_PRO_matches$Match_step != "Directly to PRO"] # Specific_protein_sugg_matches$Marker[Specific_protein_sugg_matches$Match_step == ""]

    all_ids_final <- list()

    if(length(unique_marker_names) != 0) {

      for (marker in unique_marker_names) {

        matched_row <- CD_syn[which(CD_syn == marker, arr.ind=TRUE)[,1],]

        matched_row <- data.frame(x=unlist(matched_row))

        unique_syn_marker_names <- unique(stats::na.omit(matched_row))

        unique_syn_marker_names <- as.vector(unique_syn_marker_names[,1])

        # Run SPARQL
        all_ids <- PRO_SPARQL(unique_marker_names = unique_syn_marker_names,
                              NCBI_taxon_ID = NCBI_taxon_ID_short,
                              onto_endpoint = onto_endpoint,
                              Match_step = Match_step)

        # Put all results in one data frame
        all_ids <- dplyr::bind_rows(all_ids, .id = "Matched_synonym")

        all_ids <- stats::na.omit(all_ids)

        if (length(all_ids) != 0) {
          all_ids_final[[marker]] <- all_ids
        } else if (length(all_ids) == 0) {
          all_ids_final[[marker]] <- data.frame(Matched_synonym = "", PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
        }
      }

      # Put all results in one data frame
      all_ids_final2 <- dplyr::bind_rows(all_ids_final, .id = "Marker")

      # Remove duplicated entries (Marker, PRO term, and type of match are the same)
      all_ids_final2 <- all_ids_final2[!duplicated(all_ids_final2[,c("Marker","PRO_term", "Match_type")]),]

      # Get entries were inputted marker and matched species is the same

      # Subset into entries with and without duplicates
      x <-  all_ids_final2[duplicated(all_ids_final2[,c("Marker","Species")]) | duplicated(all_ids_final2[,c("Marker","Species")], fromLast = TRUE) ,]
      y <-  dplyr::setdiff(all_ids_final2, x)

      # Get vector of markers (duplicates)
      unique_marker_names <- unique(x$Marker)

      # If there is both an exact match and a secondary match, remove the secondary
      remove_secondary_per_marker <- data.frame()

      for(marker in unique_marker_names) {
        df_subsetted <- x[x$Marker == marker,]
        if(any(df_subsetted$Match_type == "Primary")) {
          df_per_marker <- df_subsetted[df_subsetted$Match_type != "Secondary",]
          df_per_marker <- df_per_marker[df_per_marker$Match_type != "",]
          remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
        } else if(any(df_subsetted$Match_type == "Secondary")) {
          df_per_marker <- df_subsetted[df_subsetted$Match_type != "",]
          remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
        }
      }

      # Final
      CD_syn_list_matches <- rbind(remove_secondary_per_marker, y)

      # Just columns we want
      CD_syn_list_matches <- unique(CD_syn_list_matches[ , c("Marker", "PRO_term", "PRO_name", "Species", "Match_type", "Match_step")])

    } else {
      CD_syn_list_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
    }

    ######## Wikidata, then SPAPRQL query to PRO ########

    # Query Wikidata, then if there are synonyms query them in PRO

    Match_step <- "Wikidata"

    # Vector of unique marker names
    unique_marker_names <- CD_syn_list_matches$Marker[CD_syn_list_matches$Match_step != "CD synonym list"]

    # Remove duplicates introduces by synonyms from last step
    remove <- CD_syn_list_matches$Marker[CD_syn_list_matches$Match_step != ""]
    unique_marker_names <-  dplyr::setdiff(unique_marker_names, remove)

    all_ids_final <- list()

    if(length(unique_marker_names) != 0) {
      if(unique_marker_names[1] != "") {

        # Loop through each and do the SPARQL query
        for(marker in unique_marker_names) {

          specific_marker <-  shQuote(marker)

          SPARQL_query <-
            paste0(
              "
SELECT DISTINCT ?item ?itemLabel ?altLabel WHERE {
  # Look either for non-species specific or add a specific species
  OPTIONAL {
    ?item wdt:P703 <", wikidata_taxon_ID, "> .

  }
  # Get Wikidata protein name and alternative names
  ?item rdfs:label ?itemLabel .
  ?item skos:altLabel ?altLabel.
  # Return results that are either proteins or protein-coding genes
  {?item wdt:P31 wd:Q8054}
  UNION
  {?item wdt:P279 wd:Q20747295}
  # Search for the marker input, either if it is the entry name or an alternative name
  FILTER (UCASE(REPLACE(str(?itemLabel),'[ -.]','')) = ", specific_marker, " || UCASE(REPLACE(str(?altLabel),'[ -.]','')) = ", specific_marker, ")
}
")
          # Run SPARQL on the endpoint
          result <- httr::GET(
            url = wiki_endpoint,
            query = list(query = SPARQL_query),
            httr::user_agent(R.version.string))

          # Get result in text JSON
          x <- httr::content(result, as = "text") #, encoding = "UTF-8")

          # Convert from JSON to a list
          df <- jsonlite::fromJSON(x, flatten = TRUE)

          # Extract the data frame
          df <- df$results$bindings

          if (length(df) != 0) {

            # Remove unneeded info
            df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", 'lang')))

            # Remove this part that was added onto the column names
            colnames(df) <- gsub(".value", "", colnames(df))

            # Make everything case, remove spaces, dashes
            df <- dplyr::mutate_all(df, .funs = toupper)
            df <- as.data.frame(lapply(df, function(y) gsub('-', '', y)))
            df <- as.data.frame(lapply(df, function(y) gsub('_', '', y)))
            df <- as.data.frame(lapply(df, function(y) gsub(' ', '', y)))
            df <- as.data.frame(lapply(df, function(y) gsub('\\.', '', y)))

            df <- df[(which(nchar(df$altLabel) > 2)),]
            df <- df[(which(nchar(df$itemLabel) > 2)),]

            # No marker can start with a number, so add X to any marker names that start with a number
            df <- as.data.frame(lapply(df, function(y) gsub('^(\\d)', 'X\\1', y)))

            # Vector of unique marker names
            unique_syn_names <- unique(c(df[, "altLabel"], df[, "itemLabel"]))

             unique_syn_names <- unique_syn_names[!grepl("[^\x01-\x7F]+", unique_syn_names)]
                                       
            unique_syn_names <- unique_syn_names[grep('[A-Za-z]',unique_syn_names)]

            unique_syn_names <- unique_syn_names[!grepl(paste0('^', marker, '$'), unique_syn_names, useBytes = TRUE)]

            # Run SPARQL
            all_ids <- PRO_SPARQL(unique_marker_names = unique_syn_names,
                                  NCBI_taxon_ID = NCBI_taxon_ID_short,
                                  onto_endpoint = onto_endpoint,
                                  Match_step = Match_step)

            # Put all results in one data frame
            all_ids <- dplyr::bind_rows(all_ids, .id = "Matched synonym")

            all_ids <- stats::na.omit(all_ids)

            if (length(all_ids) != 0) {
              all_ids_final[[marker]] <- all_ids
            } else if (length(all_ids) == 0) {
              all_ids_final[[marker]] <- data.frame(PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
            }
          } else if (length(df) == 0) {
            all_ids_final[[marker]] <- data.frame(PRO_term = "", PRO_name = "", Species = "", Match_step = "", Match_type = "")
          }
        }

      # Put all results in one data frame
      all_ids_final2 <- dplyr::bind_rows(all_ids_final, .id = "Marker")

      # Remove duplicated entries (Marker, PRO term, and type of match are the same)
      all_ids_final2 <- all_ids_final2[!duplicated(all_ids_final2[,c("Marker","PRO_term", "Match_type")]),]

      # Get entries were inputted marker and matched species is the same

      # Subset into entries with and without duplicates
      x <-  all_ids_final2[duplicated(all_ids_final2[,c("Marker","Species")]) | duplicated(all_ids_final2[,c("Marker","Species")], fromLast = TRUE) ,]
      y <-  dplyr::setdiff(all_ids_final2, x)

      # Get vector of markers (duplicates)
      unique_marker_names <- unique(x$Marker)

      # If there is both an exact match and a secondary match, remove the secondary
      remove_secondary_per_marker <- data.frame()

      for(marker in unique_marker_names) {
        df_subsetted <- x[x$Marker == marker,]
        if(any(df_subsetted$Match_type == "Primary")) {
          df_per_marker <- df_subsetted[df_subsetted$Match_type != "Secondary",]
          df_per_marker <- df_per_marker[df_per_marker$Match_type != "",]
          remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
        } else if(any(df_subsetted$Match_type == "Secondary")) {
          df_per_marker <- df_subsetted[df_subsetted$Match_type != "",]
          remove_secondary_per_marker <- rbind(remove_secondary_per_marker, df_per_marker)
        }
      }

      # Final
      Wikidata_matches <- rbind(remove_secondary_per_marker, y)

      # Just columns we want
      Wikidata_matches <- unique(Wikidata_matches[ , c("Marker", "PRO_term", "PRO_name", "Species", "Match_type", "Match_step")])

      } else {
        Wikidata_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
      }
    } else {
      Wikidata_matches <- data.frame(Marker = "", PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
    }

    ######## Final results ########

    # Put everything together and get the final list of matched/unmatched marker names

    # Add empty suggestion column to parts where suggestions were not relevant
    Direct_PRO_matches$Suggestions <- ""
    CD_syn_list_matches$Suggestions <- ""
    Wikidata_matches$Suggestions <- ""

    # Subset entries that were matched at each step
    GO_matches2 <- GO_matches[GO_matches$Match_step == "GO exceptions",]
    Direct_PRO_matches2 <- Direct_PRO_matches[Direct_PRO_matches$Match_step == "Directly to PRO",]
    CD_syn_list_matches2 <- CD_syn_list_matches[CD_syn_list_matches$Match_step == "CD synonym list",]
    Wikidata_matches2 <- Wikidata_matches[Wikidata_matches$Match_step == "Wikidata",]

    # Put everything together into 1 dataframe
    all_ids <- do.call("rbind", list(GO_matches2,  Direct_PRO_matches2, CD_syn_list_matches2, Wikidata_matches2))

    # Put original inputted names back in
    final_df <- merge(new_df, all_ids, by.x = "Marker2", by.y = "Marker", all.x = TRUE, all.y = TRUE)

    # Change names
    colnames(final_df)[1:2] <- c("New_Name", "Original_Name")

    if (nrow(auto_replaced) != 0) {
      # Put back original auto changed names
      replacements <- auto_replaced[as.character(auto_replaced$Marker) != as.character(auto_replaced$Replacement),]

      replacements <- replacements[,1:2]

      colnames(replacements) <- c("Marker2", "New_Name")

      replacements <- merge(replacements, new_df, by = "Marker2")

      replacements <- replacements[,c(2, 3)]

      replacements <- merge(replacements, final_df, by = "New_Name")

      replacements <- replacements[,-3]

      colnames(replacements)[2] <- "Original_Name"

      final_df <- rbind(final_df, replacements)
    }

    final_df <- final_df[!is.na(final_df$Original_Name), ]

    # If there was no match or suggestion, indicate that
    final_df$Match_step[is.na(final_df$Match_step)] <- "Did not match"

    # Turn any remaining NA entries into blanks
    final_df[is.na(final_df)] <- ""

    # Change species URIs into human readable form
    final_df$Species <- gsub(NCBI_taxon_ID_full, input$species_choice_2, final_df$Species)
    keep <- c(input$species_choice_2, "Not species specific")

    # Remove IDs that aren't actual species
    remove_non_species <- dplyr::setdiff(final_df$Species, keep)
    remove_non_species <- remove_non_species[remove_non_species != ""]

    if (length(remove_non_species) != 0) {
      # Fill these in as species unspecific
      final_df$Species <- stringi::stri_replace_all_regex(final_df$Species,
                                                          pattern=remove_non_species,
                                                          replacement="Not species specific",
                                                          vectorize=FALSE)
    }

    new_df <- reformatted_data_2()

    # Remove duplicates
    final_df <- final_df[!duplicated(final_df), ]

    # Get new column names
    colnames(final_df) <-
      c("Revised Name",
        "Original Name",
        "PRO/GO Term",
        "PRO/GO Name",
        "Species Specific",
        "Match type",
        "Workflow Step",
        "Suggestion")

    return(final_df)

  })

  # Step 2, option to delete PRO entries before doing a cell type matching search
  shiny::observeEvent(input$delete_PRO_synonym_2, {
    if (!is.null(input$table_protein_synonym_2_rows_selected)) {
      df2_2$data <<-
        df2_2$data[-as.numeric(input$table_protein_synonym_2_rows_selected), ]
    }
  })

  # Step 2, remove link from PRO/GO term for download, 1st dataset
  downloadable_step_2 <- reactive({
  shiny::req(input$submit_tab2_step2)
   
    remove_link <- df2_2$data

    if (all(remove_link$`PRO/GO Term` == "Did not match")) {
      return(remove_link)
    } else {
      # Remove syntax that made this linkable in R Shiny
      remove_link$`PRO/GO Term` <-
        gsub("' target='_blank'>[^>]+</a>",
             "",
             remove_link$`PRO/GO Term`)
      remove_link$`PRO/GO Term` <-
        gsub("<a href='", "", remove_link$`PRO/GO Term`)

      # Remove any NA entries
      remove_link <- remove_link[stats::complete.cases(remove_link), ]

      return(remove_link)
    }
  })

  # Step 2, 2nd dataset, all info
  downloadable_step_2_long <-  reactive({
  shiny::req(input$submit_tab2_step2)
   
    markers <- downloadable_step_2()

    signs_clusters <- reformatted_data_2()

    all <- merge(signs_clusters, markers, by.x = "Marker", by.y = "Original Name")

    # Order by cluster
    all <- all[order(all$Cluster, decreasing = FALSE), ]

    return(all)

  })

  # Step 2, download PRO/GO matches
  output$download_PRO_GO_matches_2 <- shiny::downloadHandler(

    # Filename
    filename = function() {
      if (input$input_type == 'csv_upload') {
        paste0(gsub(".csv|.xls|.xlsx", "", input$csv_upload),
               "-matched-PRO-GO-terms-", Sys.Date(), ".zip")
      } else {
        paste0("typed-markers-matched-PRO-GO-terms-", Sys.Date(), ".zip")
      }
    },

    content = function(file) {
      # definition of content to download
      to_dl <- list(
        # names to use in file names
        names = list(a = "Individual_PRO_GO_terms",
                     b = "All_PRO_GO_terms"),
        # data
        data = list(a = downloadable_step_2(),
                    b = downloadable_step_2_long())
      )

      # temp dir for the csv's as we can only create
      # an archive from existent files and not data from R
      twd <- setwd(tempdir())
      on.exit(setwd(twd))
      files <- NULL

      # loop on data to download and write individual csv's
      for (i in c("a", "b")) {
        fileName <- paste0(to_dl[["names"]][[i]], ".csv") # csv file name
        utils::write.csv(to_dl[["data"]][[i]], fileName, row.names = FALSE) # write csv in temp dir

        files <- c(files, fileName) # store written file name
      }
      # create zip file with all the csv files
      zip::zip(file, files)
    }
  )

  ######~#~# Step 3 - Match Median Difference Equation results and PRO/GO terms with cell types, show user ranked results and CL definitions #~#~######

  # Step 3, side panel

  # Reset everything
  shiny::observeEvent(input$reset_tab2_step3, {
    session$reload()
  })

  # Step 3, main panel

  # Step 3, when help button is clicked show message
  shiny::observeEvent(input$help_tab2_step3,{
    shiny::showModal(shiny::modalDialog(
      title = HTML("Step 3<br>Match Marker Patterns to Cell Types"),
      HTML("This step displays the results of the cell type matching process. The data frame shows:
      <br>
      <br>
      <li>Input clusters</li>
      <br>
      <li>Cell Ontology (CL) identification terms</li>
      <br>
      <li>Descriptions of the matched cell types</li>
      <br>
      <li>Wikidata identification terms</li>
      <br>
      <li>Markers that are matched between the input clusters and the CL cell types</li>
      <br>
      <li>Markers present in the input clusters</li>
      <br>
      <li>Markers present in the matched CL cell types</li>
      <br>
      <li>Markers that are directly contradicted between the input clusters and the CL cell types</li>
      <br>
      <br>
      Both CL and Wikidata terms are clickable and lead to their respective web pages.
      <br>
      <br>
      Three additional columns provide statistical metrics for evaluating the match:
      <br>
      <br>
      <li>Percentages of matched markers over the total markers in the CL cell type</li>
      <br>
      <li>Percentages of matched markers over the total markers in the input</li>
      <br>
      <li>The matching scores derived from Equation 3</li>
      <br>
      The 'Download' button at the bottom allows users to export a CSV file containing all relevant data, including additional columns detailing the number of markers in the input, CL cell types, matches, and contradictions.
      <br>
      <br>"),
      tags$div("Equation 3 and additional information can be found", tags$a(href = 'https://github.com/AmandaRT18/Cell.Naming/wiki/4.1.-Input-Marker-Descriptors-and-Default-References', "here.", target="_blank")),
      easyClose = TRUE))
  })

  # Step 3, Match markers (PRO terms) to cell types in CL/pCL
  # Get wikidata link for each cell type
  return_cell_types_2 <- shiny::eventReactive(input$submit_tab2_step3, {
   browser()
    # Load in data from Step 3
    input_reformatted <- reformatted_data_2()

    # Load in data from Step 4
    all_PRO_matches <- df2_2$data

    all_PRO_matches <- all_PRO_matches[all_PRO_matches$Reference == "Included in reference.",]

    all_PRO_matches <- all_PRO_matches [ , !(names(all_PRO_matches ) %in% "Reference")]

    all_id_results <- input_reformatted

    colnames(all_id_results)[which(names(all_id_results) == "Marker")] <-
      "Original Name"

    all_id_results <-
      merge(all_id_results,
            all_PRO_matches,
            by.x = "Original Name",
            by.y = "Original Name")

    if (all(all_id_results$`PRO/GO Term` == "No match found")) {
      return(NULL)
    } else {

      # Remove syntax that made this linkable in R Shiny
      all_id_results$`PRO/GO Term` <-
        gsub("' target='_blank'>[^>]+</a>",
             "",
             all_id_results$`PRO/GO Term`)
      all_id_results$`PRO/GO Term` <-
        gsub("<a href='", "", all_id_results$`PRO/GO Term`)

      # Add back in < and >
      all_id_results$`PRO/GO Term` <-
        paste0("<", all_id_results$`PRO/GO Term`, ">")

      all_id_results2 <- all_id_results

      ## Get contradictions

      # Add empty columns to be filled
      all_id_results[, 'Relations_ontology_1'] <- NA
      all_id_results[, 'Relations_ontology_2'] <- NA
      all_id_results[, 'Relations_ontology_3'] <- NA
      all_id_results[, 'Relations_ontology_4'] <- NA
      all_id_results[, 'Relations_ontology_5'] <- NA

      #Loop through data, make positives, negatives, high, lows into ontology relations terminology
      for (a in 1:nrow(all_id_results)) {

        # Positive signs are always contradicted by negatives
        if (all_id_results$Sign[a] == "Positive") {
          all_id_results$Relations_ontology_1[a] <-
            "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
          all_id_results$Relations_ontology_2[a] <-
            "<http://purl.obolibrary.org/obo/CL_4030046>"

          all_id_results$Marker_query[a] <-
            paste(
              "?Related_protein =",
              all_id_results$`PRO/GO Term`[a],
              "&& ?Relationship_type = ",
              all_id_results$Relations_ontology_1[a],
              " || ?Related_protein = ",
              all_id_results$`PRO/GO Term`[a],
              "&& ?Relationship_type = ",
              all_id_results$Relations_ontology_2[a]
            )
        }

        # If choosing low/negative, contradictions of negatives are any positive except low
        if (input$low_option_2 == "low_negative_2") {

          if (all_id_results$Sign[a] == "Negative") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a]
              )
            # If not low/negative, the contradictions of positive is any positive, high, or low sign
          }} else  if (input$low_option_2 != "low_negative_2") {

            if (all_id_results$Sign[a] == "Negative") {
              all_id_results$Relations_ontology_1[a] <-
                "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
              all_id_results$Relations_ontology_2[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
              all_id_results$Relations_ontology_3[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part
              all_id_results$Relations_ontology_4[a] <-
                "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

              all_id_results$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results$Relations_ontology_1[a],
                  " || ?Related_protein = ",
                  all_id_results$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results$Relations_ontology_2[a],
                  " || ?Related_protein = ",
                  all_id_results$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results$Relations_ontology_3[a],
                  " || ?Related_protein = ",
                  all_id_results$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results$Relations_ontology_4[a]
                )
            }}

        # If choosing high/positive, then any high signs are contradicted by negatives and lows
        if (input$high_option_2 == "high_positive_2") {
          if (all_id_results$Sign[a] == "High") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>" # lacks plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a]
              )
          }}

        # If choosing high/high, then any high signs are contradicted by negatives, positives, and lows
        if (input$high_option_2 == "high_high_2") {
          if (all_id_results$Sign[a] == "High") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>" # lacks plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part
            all_id_results$Relations_ontology_4[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results$Relations_ontology_5[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_4[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_5[a]
              )
          }}


        # If choosing low/positive, then any low signs are contradicted by negatives and highs
        if (input$low_option_2 == "low_positive_2") {
          if (all_id_results$Sign[a] == "Low") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>" # lacks plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a]
              )
          }}

        # If choosing low/low, then any low signs are contradicted by negatives, positives, and high
        if (input$low_option_2 == "low_low_2") {
          if (all_id_results$Sign[a] == "Low") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>" # lacks plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
            all_id_results$Relations_ontology_4[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results$Relations_ontology_5[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_4[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_5[a]
              )
          }}


        # If choosing low/negative, then any low signs are contradicted by positives and high
        if (input$low_option_2 == "low_negative_2") {
          if (all_id_results$Sign[a] == "Low") {
            all_id_results$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
            all_id_results$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results$Relations_ontology_3[a]
              )
          }}
      }

      # Make sure cluster name/number is read in as string
      all_id_results$Cluster <- as.character(all_id_results$Cluster)

      # Empty list to put results in
      marker_in_cluster_matches <- list()

      # Look through each group/cluster category
      for (a_cluster in sort(unique(all_id_results$Cluster))) {
        # Just the specific grouping
        specific_cluster <-
          all_id_results[all_id_results$Cluster == a_cluster, ]

        # Remove proteins that didn't get a match
        specific_cluster <-
          specific_cluster[specific_cluster$`PRO/GO Term` != "No match found", ]

        for (b in sort(unique(specific_cluster$`Revised Name`))) {
          # Just the specific marker
          specific_marker <-
            specific_cluster[specific_cluster$`Revised Name` == b, ]

          # Filter for SPARQL query
          full_filter <-
            paste(specific_marker$Marker_query, collapse = "||")

          if (length(input$ontology_type_2) == 1) {
            if (input$ontology_type_2 == 'internal_CL_2') {
              ontology_URI <- "<http://purl.obolibrary.org/obo/merged/CL>"
              full_filter_1 <-
                paste(
                  "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/CL_\"  ) &&",
                  full_filter
                )
            } else if (input$ontology_type_2 == 'internal_pCL_2') {
              ontology_URI <- "<http://purl.obolibrary.org/obo/merged/PCL>"
              full_filter_1 <-
                paste(
                  "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/PCL_\"  ) &&",
                  full_filter
                )
            }

            # SPARQL query
            query1 <-

              paste0(
                "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                ontology_URI,
                "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }



  FILTER(",
                full_filter_1,
                ")
    }
  "
              )
            # Run SPARQL on the endpoint
            result1 <- httr::POST(onto_endpoint,
                                  body = list(query = query1),
                                  httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result1)

            # Get result in text JSON
            x1 <- httr::content(result1, "text", encoding = "UTF-8")

            # Convert from JSON to a list
            df <- jsonlite::fromJSON(x1, flatten = TRUE)

            # Extract the dataframe
            df <- df$results$bindings

            # Put the results in a list of data frames
            if (length(df) != 0) {

              # Remove unneeded info
              # df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype")))
              df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

              # Remove this part that was added onto the column names
              colnames(df) <- gsub(".value", "", colnames(df))

              # Make sure they all have these columns
              six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

              df[six_cols[!(six_cols %in% colnames(df))]] <- NA
            }
          } else if (length(input$ontology_type_2) == 2) {
            ontology_URI_1 <- "<http://purl.obolibrary.org/obo/merged/CL>"
            full_filter_1 <-
              paste(
                "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/CL_\"  ) &&",
                full_filter
              )

            # SPARQL query
            query1 <-

              paste0(
                "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                ontology_URI_1,
                "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }



  FILTER(",
                full_filter_1,
                ")
    }
  "
              )
            # Run SPARQL on the endpoint
            result1 <- httr::POST(onto_endpoint,
                                  body = list(query = query1),
                                  httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result1)

            # Get result in text JSON
            x1 <- httr::content(result1, "text", encoding = "UTF-8")

            # Convert from JSON to a list
            df1 <- jsonlite::fromJSON(x1, flatten = TRUE)

            # Extract the dataframe
            df1 <- df1$results$bindings

            # Put the results in a list of data frames
            if (length(df1) != 0) {

              # Remove unneeded info
              df1 <- df1 %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

              # Remove this part that was added onto the column names
              colnames(df1) <- gsub(".value", "", colnames(df1))

              # Make sure they all have these columns
              six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

              df1[six_cols[!(six_cols %in% colnames(df1))]] <- NA
            }

            ontology_URI_2 <-
              "<http://purl.obolibrary.org/obo/merged/PCL>"
            full_filter_2 <-
              paste(
                "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/PCL_\"  ) &&",
                full_filter
              )

            # SPARQL query
            query2 <-

              paste0(
                "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                ontology_URI_2,
                "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }

  FILTER(",
                full_filter_2,
                ")
    }
  "
              )
            # Run SPARQL on the endpoint
            result2 <- httr::POST(onto_endpoint,
                                  body = list(query = query2),
                                  httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result2)

            # Get result in text JSON
            x2 <- httr::content(result2, "text", encoding = "UTF-8")

            # Convert from JSON to a list
            df2 <- jsonlite::fromJSON(x2, flatten = TRUE)

            # Extract the dataframe
            df2 <- df2$results$bindings

            # Put the results in a list of data frames
            if (length(df2) != 0) {

              # Remove unneeded info
              df2 <- df2 %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

              # Remove this part that was added onto the column names
              colnames(df2) <- gsub(".value", "", colnames(df2))

              # Make sure they all have these columns
              six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

              df2[six_cols[!(six_cols %in% colnames(df2))]] <- NA
            }

            # Add results
            df <- rbind(df1, df2)
          }

          # Put the results in a list of data frames
          if (length(df) != 0) {
            # Add back in < and >
            df$CL_term <- paste0("<", df$CL_term, ">")
            for (c in 1:nrow(df)) {
              if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0002104") {
                df$Relationship_type[c] <- "Positive"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/BFO_0000051") {
                df$Relationship_type[c] <- "Positive"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015015") {
                df$Relationship_type[c] <- "High"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015016") {
                df$Relationship_type[c] <- "Low"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part") {
                df$Relationship_type[c] <- "Negative"
              } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/CL_4030046") {
                df$Relationship_type[c] <- "Negative"
              }
            }

            if(!"Cell_comment" %in% colnames(df)) {
              df$Cell_comment = NA
            }

            # Collapse to get ensure each marker is only counted 1 time per cell cluster match
            # Some entries may appear more then once bc there is more then 1 marker relation for that protein or bc a protein has more then one PRO ID
            df <- df %>%
              dplyr::group_by(CL_term, Cell_type, Cell_comment) %>%
              dplyr::summarise(
                Relationship_type = paste(Relationship_type, collapse = ", "),
                Related_protein = paste(Related_protein, collapse = ", "),
                Protein_name = paste(Protein_name, collapse = ", ")
              )

            # Remove repeats that the last step may have formed
            df$Relationship_type <-
              sapply(df$Relationship_type, function(x)
                paste(unique(unlist(
                  strsplit(x, ", ")
                )), collapse = ", "))
            df$Related_protein <-
              sapply(df$Related_protein, function(x)
                paste(unique(unlist(
                  strsplit(x, ", ")
                )), collapse = ", "))
            df$Protein_name <-
              sapply(df$Protein_name, function(x)
                paste(unique(unlist(
                  strsplit(x, ", ")
                )), collapse = ", "))

            # Include cluster name/number
            df <- cbind(Cluster = a_cluster, df)

          } else if (length(df) == 0) {
            df <- data.frame(CL_term = "No match found")
          }

          # Name
          cluster_marker_name <- paste0(a_cluster, "::", b)

          # Put into list
          marker_in_cluster_matches[[cluster_marker_name]] <-
            as.data.frame(df)
        }
      }

      # Put all results in one data frame
      all_cluster_matches <-
        dplyr::bind_rows(marker_in_cluster_matches, .id = "Marker")

      # Remove proteins that didn't get a match in a particular cluster
      all_cluster_matches <-
        all_cluster_matches[all_cluster_matches$CL_term != "No match found", ]

      # Remove cluster from unique name
      all_cluster_matches$Marker <-
        gsub(".*::", "", all_cluster_matches$Marker)

      if(nrow(all_cluster_matches) >= 1) {

        # New column
        all_cluster_matches$`Matched inputted marker` <-
          paste0(all_cluster_matches$Marker,
                 " (",
                 all_cluster_matches$Relationship_type,
                 ")")

        all_cluster_matches$Cluster <- as.numeric(all_cluster_matches$Cluster)

        # Sort, put into final format for export
        final_matched_output <- all_cluster_matches %>%
          dplyr::group_by(Cluster,
                          CL_term,
                          Cell_type,
                          Cell_comment) %>%
          dplyr::summarise(
            `Matched inputted marker` = paste(`Matched inputted marker`, collapse = ", "),
            `num of inputs that match to cell type` =  dplyr::n()) %>%
          dplyr::arrange(Cluster, dplyr::desc(`num of inputs that match to cell type`))

        # Remove '<' and '>' from URI
        final_matched_output$CL_term <-
          gsub('<', '', final_matched_output$CL_term)
        final_matched_output$CL_term <-
          gsub('>', '', final_matched_output$CL_term)

        final_matched_output2 <- final_matched_output

        # Final column name order
        final_matched_output2 <- final_matched_output2[, c(
          "CL_term",
          "Cluster",
          "Matched inputted marker",
          "num of inputs that match to cell type"

        )]

        colnames(final_matched_output2)[3] <- "Contradiction(s)"

        colnames(final_matched_output2)[4] <- "num contradiction(s)"

        ## Get matches

        # Add empty columns to be filled
        all_id_results2[, 'Relations_ontology_1'] <- NA
        all_id_results2[, 'Relations_ontology_2'] <- NA
        all_id_results2[, 'Relations_ontology_3'] <- NA
        all_id_results2[, 'Relations_ontology_4'] <- NA

        #Loop through data, make positives, negatives, high, lows into ontology relations terminology
        for (a in 1:nrow(all_id_results2)) {
          if (all_id_results2$Sign[a] == "Positive") {
            all_id_results2$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
            all_id_results2$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part
            all_id_results2$Relations_ontology_3[a] <-
              "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part
            all_id_results2$Relations_ontology_4[a] <-
              "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part

            all_id_results2$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_2[a],
                " || ?Related_protein = ",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_3[a],
                " || ?Related_protein = ",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_4[a]
              )
          }

          if (all_id_results2$Sign[a] == "Negative") {
            all_id_results2$Relations_ontology_1[a] <-
              "<http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>"
            all_id_results2$Relations_ontology_2[a] <-
              "<http://purl.obolibrary.org/obo/CL_4030046>"

            all_id_results2$Marker_query[a] <-
              paste(
                "?Related_protein =",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_1[a],
                " || ?Related_protein = ",
                all_id_results2$`PRO/GO Term`[a],
                "&& ?Relationship_type = ",
                all_id_results2$Relations_ontology_2[a]
              )
          }

          if (input$high_option_2 == "high_positive_2") {
            if (all_id_results2$Sign[a] == "High") {
              all_id_results2$Relations_ontology_1[a] <-
                "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
              all_id_results2$Relations_ontology_2[a] <-
                "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part
              all_id_results2$Relations_ontology_3[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_2[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_3[a]
                )
            }}

          if (input$high_option_2 == "high_high_2") {
            if (all_id_results2$Sign[a] == "High") {
              all_id_results2$Relations_ontology_1[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015015>" # has high plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a]
                )
            }}


          if (input$low_option_2 == "low_positive_2") {
            if (all_id_results2$Sign[a] == "Low") {
              all_id_results2$Relations_ontology_1[a] <-
                "<http://purl.obolibrary.org/obo/RO_0002104>" # has plasma membrane part
              all_id_results2$Relations_ontology_2[a] <-
                "<http://purl.obolibrary.org/obo/BFO_0000051>" # has part
              all_id_results2$Relations_ontology_3[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_2[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_3[a]
                )
            }}


          if (input$low_option_2 == "low_low_2") {
            if (all_id_results2$Sign[a] == "Low") {
              all_id_results2$Relations_ontology_1[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a]
                )
            }}

          if (input$low_option_2 == "low_negative_2") {
            if (all_id_results2$Sign[a] == "Low") {
              all_id_results2$Relations_ontology_1[a] <-
                "http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part"
              all_id_results2$Relations_ontology_2[a] <-
                "<http://purl.obolibrary.org/obo/CL_4030046>"
              all_id_results2$Relations_ontology_3[a] <-
                "<http://purl.obolibrary.org/obo/RO_0015016>" # has low plasma membrane part

              all_id_results2$Marker_query[a] <-
                paste(
                  "?Related_protein =",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_1[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_2[a],
                  " || ?Related_protein = ",
                  all_id_results2$`PRO/GO Term`[a],
                  "&& ?Relationship_type = ",
                  all_id_results2$Relations_ontology_3[a]
                )
            }}
        }

        # Make sure cluster name/number is read in as string
        all_id_results2$Cluster <- as.character(all_id_results2$Cluster)

        # Empty list to put results in
        marker_in_cluster_matches <- list()

        # Look through each group/cluster category
        for (a_cluster in sort(unique(all_id_results2$Cluster))) {
          # Just the specific grouping
          specific_cluster <-
            all_id_results2[all_id_results2$Cluster == a_cluster, ]

          # Remove proteins that didn't get a match
          specific_cluster <-
            specific_cluster[specific_cluster$`PRO/GO Term` != "No match found", ]

          for (b in sort(unique(specific_cluster$`Revised Name`))) {
            # Just the specific marker
            specific_marker <-
              specific_cluster[specific_cluster$`Revised Name` == b, ]

            # Filter for SPARQL query
            full_filter <-
              paste(specific_marker$Marker_query, collapse = "||")


            if (length(input$ontology_type_2) == 1) {
              if (input$ontology_type_2 == 'internal_CL_2') {
                ontology_URI <- "<http://purl.obolibrary.org/obo/merged/CL>"
                full_filter_1 <-
                  paste(
                    "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/CL_\"  ) &&",
                    full_filter
                  )
              } else if (input$ontology_type_2 == 'internal_pCL_2') {
                ontology_URI <- "<http://purl.obolibrary.org/obo/merged/PCL>"
                full_filter_1 <-
                  paste(
                    "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/PCL_\"  ) &&",
                    full_filter
                  )
              }

              # SPARQL query
              query1 <-

                paste0(
                  "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                  ontology_URI,
                  "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }



  FILTER(",
                  full_filter_1,
                  ")
    }
  "
                )
              # Run SPARQL on the endpoint
              result1 <- httr::POST(onto_endpoint,
                                    body = list(query = query1),
                                    httr::user_agent(R.version.string))

              # Will show a warning/error if there is any
              httr::stop_for_status(result1)

              # Get result in text JSON
              x1 <- httr::content(result1, "text", encoding = "UTF-8")

              # Convert from JSON to a list
              df <- jsonlite::fromJSON(x1, flatten = TRUE)

              # Extract the dataframe
              df <- df$results$bindings

              # Put the results in a list of data frames
              if (length(df) != 0) {

                # Remove unneeded info
                # df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype")))
                df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

                # Remove this part that was added onto the column names
                colnames(df) <- gsub(".value", "", colnames(df))

                # Make sure they all have these columns
                six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

                df[six_cols[!(six_cols %in% colnames(df))]] <- NA
              }
            } else if (length(input$ontology_type_2) == 2) {
              ontology_URI_1 <- "<http://purl.obolibrary.org/obo/merged/CL>"
              full_filter_1 <-
                paste(
                  "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/CL_\"  ) &&",
                  full_filter
                )

              # SPARQL query
              query1 <-

                paste0(
                  "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                  ontology_URI_1,
                  "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }



  FILTER(",
                  full_filter_1,
                  ")
    }
  "
                )
              # Run SPARQL on the endpoint
              result1 <- httr::POST(onto_endpoint,
                                    body = list(query = query1),
                                    httr::user_agent(R.version.string))

              # Will show a warning/error if there is any
              httr::stop_for_status(result1)

              # Get result in text JSON
              x1 <- httr::content(result1, "text", encoding = "UTF-8")

              # Convert from JSON to a list
              df1 <- jsonlite::fromJSON(x1, flatten = TRUE)

              # Extract the dataframe
              df1 <- df1$results$bindings

              # Put the results in a list of data frames
              if (length(df1) != 0) {

                # Remove unneeded info
                df1 <- df1 %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

                # Remove this part that was added onto the column names
                colnames(df1) <- gsub(".value", "", colnames(df1))

                # Make sure they all have these columns
                six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

                df1[six_cols[!(six_cols %in% colnames(df1))]] <- NA
              }

              ontology_URI_2 <-
                "<http://purl.obolibrary.org/obo/merged/PCL>"
              full_filter_2 <-
                paste(
                  "strstarts(str(?CL_term), \"http://purl.obolibrary.org/obo/PCL_\"  ) &&",
                  full_filter
                )

              # SPARQL query
              query2 <-

                paste0(
                  "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Cell_type ?Cell_comment ?Relationship_type ?Related_protein ?Protein_name

  FROM ",
                  ontology_URI_2,
                  "

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term rdf:type ?type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  OPTIONAL {
    ?CL_term rdfs:comment ?Cell_comment .
    }

  OPTIONAL {
    ?Class_of owl:onProperty ?Relationship_type .
    }

  OPTIONAL {
    ?Class_of owl:someValuesFrom ?Related_protein .
    ?Related_protein rdfs:label ?Protein_name
  }

  FILTER(",
                  full_filter_2,
                  ")
    }
  "
                )
              # Run SPARQL on the endpoint
              result2 <- httr::POST(onto_endpoint,
                                    body = list(query = query2),
                                    httr::user_agent(R.version.string))

              # Will show a warning/error if there is any
              httr::stop_for_status(result2)

              # Get result in text JSON
              x2 <- httr::content(result2, "text", encoding = "UTF-8")

              # Convert from JSON to a list
              df2 <- jsonlite::fromJSON(x2, flatten = TRUE)

              # Extract the dataframe
              df2 <- df2$results$bindings

              # Put the results in a list of data frames
              if (length(df2) != 0) {

                # Remove unneeded info
                df2 <- df2 %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

                # Remove this part that was added onto the column names
                colnames(df2) <- gsub(".value", "", colnames(df2))

                # Make sure they all have these columns
                six_cols <- c("CL_term","Cell_type", "Cell_comment",	"Relationship_type",	"Related_protein",	"Protein_name")

                df2[six_cols[!(six_cols %in% colnames(df2))]] <- NA
              }

              # Add results
              df <- rbind(df1, df2)
            }


            # Put the results in a list of data frames
            if (length(df) != 0) {
              # Add back in < and >
              df$CL_term <- paste0("<", df$CL_term, ">")
              for (c in 1:nrow(df)) {
                if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0002104") {
                  df$Relationship_type[c] <- "Positive"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/BFO_0000051") {
                  df$Relationship_type[c] <- "Positive"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015015") {
                  df$Relationship_type[c] <- "High"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015016") {
                  df$Relationship_type[c] <- "Low"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part") {
                  df$Relationship_type[c] <- "Negative"
                } else if (df$Relationship_type[c] == "http://purl.obolibrary.org/obo/CL_4030046") {
                  df$Relationship_type[c] <- "Negative"
                }
              }

              if(!"Cell_comment" %in% colnames(df)) {
                df$Cell_comment = NA
              }

              # Collapse to get ensure each marker is only counted 1 time per cell cluster match
              # Some entries may appear more then once bc there is more then 1 marker relation for that protein or bc a protein has more then one PRO ID
              df <- df %>%
                dplyr::group_by(CL_term, Cell_type, Cell_comment) %>%
                dplyr::summarise(
                  Relationship_type = paste(Relationship_type, collapse = ", "),
                  Related_protein = paste(Related_protein, collapse = ", "),
                  Protein_name = paste(Protein_name, collapse = ", ")
                )

              # Remove repeats that the last step may have formed
              df$Relationship_type <-
                sapply(df$Relationship_type, function(x)
                  paste(unique(unlist(
                    strsplit(x, ", ")
                  )), collapse = ", "))
              df$Related_protein <-
                sapply(df$Related_protein, function(x)
                  paste(unique(unlist(
                    strsplit(x, ", ")
                  )), collapse = ", "))
              df$Protein_name <-
                sapply(df$Protein_name, function(x)
                  paste(unique(unlist(
                    strsplit(x, ", ")
                  )), collapse = ", "))

              # Include cluster name/number
              df <- cbind(Cluster = a_cluster, df)

            } else if (length(df) == 0) {
              df <- data.frame(CL_term = "No match found")
            }

            # Name
            cluster_marker_name <- paste0(a_cluster, "::", b)

            # Put into list
            marker_in_cluster_matches[[cluster_marker_name]] <-
              as.data.frame(df)
          }
        }

        # Put all results in one data frame
        all_cluster_matches <-
          dplyr::bind_rows(marker_in_cluster_matches, .id = "Marker")

        # Remove proteins that didn't get a match in a particular cluster
        all_cluster_matches <-
          all_cluster_matches[all_cluster_matches$CL_term != "No match found", ]

        # Remove cluster from unique name
        all_cluster_matches$Marker <-
          gsub(".*::", "", all_cluster_matches$Marker)

        if(nrow(all_cluster_matches) >= 1) {

          # New column
          all_cluster_matches$`Matched inputted marker` <-
            paste0(all_cluster_matches$Marker,
                   " (",
                   all_cluster_matches$Relationship_type,
                   ")")


          # Sort, put into final format for export
          final_matched_output <- all_cluster_matches %>%
            dplyr::group_by(Cluster,
                            CL_term,
                            Cell_type,
                            Cell_comment) %>%
            dplyr::summarise(
              `Matched inputted marker` = paste(`Matched inputted marker`, collapse = ", "),
              `num of inputs that match to cell type` =  dplyr::n()
            ) %>%
            dplyr::arrange(Cluster, dplyr::desc(`num of inputs that match to cell type`))


          unq_cell_types <- final_matched_output[!duplicated(final_matched_output[,"CL_term"]),]

          # Put all cell types in another SPARQL query to get full set of proteins
          unq_cell_types$Cell_type_query <-
            paste(
              "strstarts(str(?Related_protein), \"http://purl.obolibrary.org/obo/PR_\") && str(?CL_term) =",
              unq_cell_types$CL_term, "|| strstarts(str(?Related_protein), \"http://purl.obolibrary.org/obo/GO_\") && str(?CL_term) =",
              unq_cell_types$CL_term
            )

          # Split up into smaller chunks to the SPARQL query can handle the filter
          # Number of cell types per query = 100
          chunk <- 70 #100
          total_num_rows <- nrow(unq_cell_types)
          how_to_split  <- rep(1:ceiling(total_num_rows/chunk),each=chunk)[1:total_num_rows]
          final_matched_output_list <- split(unq_cell_types,how_to_split)

          # Empty dataframe to append to
          df_all_pro <- data.frame()

          for(i in 1:length(final_matched_output_list)){

            each_part <- final_matched_output_list[[i]]

            full_filter <-
              paste(each_part$Cell_type_query, collapse = " || ")

            query <-

              paste0(
                "
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

  SELECT DISTINCT ?CL_term ?Relationship_type ?Protein_name

  FROM <http://purl.obolibrary.org/obo/merged/CL>

  WHERE
  {

  ?CL_term rdfs:label ?Cell_type .

  ?CL_term  rdfs:subClassOf* ?Class_of .

  ?Class_of owl:someValuesFrom ?Related_protein .

  ?Related_protein rdfs:label ?Protein_name .

  ?Class_of owl:onProperty ?Relationship_type .

  FILTER(",
                full_filter,
                ")

    }
  "
              )

            # Run SPARQL on the endpoint
            result <- httr::POST(onto_endpoint,
                                 body = list(query = query),
                                 httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result)

            # Get result in text JSON
            x <- httr::content(result, "text", encoding = "UTF-8")

            # Convert from JSON to a list
            df_all_pro_part <- jsonlite::fromJSON(x, flatten = TRUE)

            # Extract the dataframe
            df_all_pro_part <- df_all_pro_part$results$bindings

            # Remove unneeded info
            df_all_pro_part <- df_all_pro_part %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

            # Remove this part that was added onto the column names
            colnames(df_all_pro_part) <- gsub(".value", "", colnames(df_all_pro_part))

            # Make sure they all have these columns
            three_cols <- c("CL_term", "Relationship_type", "Protein_name")
            df_all_pro_part[three_cols[!(three_cols %in% colnames(df_all_pro_part))]] <- NA


            df_all_pro <- rbind(df_all_pro, df_all_pro_part)
          }

          # Change relations into human readable
          for (c in 1:nrow(df_all_pro)) {
            if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0002104") {
              df_all_pro$Relationship_type[c] <- "Positive"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/BFO_0000051") {
              df_all_pro$Relationship_type[c] <- "Positive"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015015") {
              df_all_pro$Relationship_type[c] <- "High"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/RO_0015016") {
              df_all_pro$Relationship_type[c] <- "Low"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part") {
              df_all_pro$Relationship_type[c] <- "Negative"
            } else if (df_all_pro$Relationship_type[c] == "http://purl.obolibrary.org/obo/CL_4030046") {
              df_all_pro$Relationship_type[c] <- "Negative"
            } else {
              df_all_pro$Relationship_type[c] <- ""
            }
          }

          # Remove proteins that didn't get a match in a particular cluster
          df_all_pro <- df_all_pro[df_all_pro$Relationship_type != "", ]

          # Remove any duplicate rows
          df_all_pro <- df_all_pro[!duplicated(df_all_pro),]

          # If protein as listed as having 2 different relationship types
          df_all_pro <- df_all_pro %>%
            dplyr::group_by(CL_term, Protein_name) %>%
            dplyr::summarise(
              Relationship_type = paste(Relationship_type, collapse = ", ")
            )

          df_all_pro$x <-
            paste0(df_all_pro$Protein_name,
                   " (",
                   df_all_pro$Relationship_type,
                   ")")

          # If a relationship type is negative and positive remove entirely
          # Shouldn't be happening unless mistake in CL
          df_all_pro <- df_all_pro[!grepl("Negative, Positive", df_all_pro$Relationship_type),]
          df_all_pro <- df_all_pro[!grepl("Positive, Negative", df_all_pro$Relationship_type),]

          # Sort, put into final format for export
          df_all_pro <- df_all_pro %>%
            dplyr::group_by(CL_term) %>%
            dplyr::summarise(
              `Full marker description` = paste(x, collapse = ", "),
              `total number of markers` = dplyr::n()
            )

          # Remove '<' and '>' from URI
          final_matched_output$CL_term <-
            gsub('<', '', final_matched_output$CL_term)
          final_matched_output$CL_term <-
            gsub('>', '', final_matched_output$CL_term)

          final_matched_output <-
            merge(final_matched_output,
                  df_all_pro,
                  by.x = "CL_term",
                  by.y = "CL_term",
                  all.x = TRUE)

          # Get wikidata link
          for (i in 1:nrow(final_matched_output)) {

            CL_term <- final_matched_output[i,"CL_term"]

            # Make CL term to be in the proper format
            CL_term <- gsub(".*/","",CL_term)
            CL_term <- shQuote(CL_term)

            SPARQL_query <-

              paste0(
                "SELECT ?item ?itemLabel WHERE {
  ?item wdt:P7963 ", CL_term, ".
  ?item wdt:P31 wd:Q189118

  SERVICE wikibase:label { bd:serviceParam wikibase:language 'en'. }
}"
              )

            # Run SPARQL on the endpoint
            result <- httr::GET(
              url = wiki_endpoint,
              query = list(query = SPARQL_query),
              httr::user_agent(R.version.string))

            # Will show a warning/error if there is any
            httr::stop_for_status(result)

            # Get result in text JSON
            x <- httr::content(result, as = "text") #, encoding = "UTF-8")

            # Convert from JSON to a list
            df <- jsonlite::fromJSON(x, flatten = TRUE)

            # Extract the data frame
            df <- df$results$bindings

            # Put the results in a list of data frames
            if (length(df) != 0) {

              # Remove unneeded info
              df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

              # Remove this part that was added onto the column names
              colnames(df) <- gsub(".value", "", colnames(df))

              if (length(df) == 1) {
                final_matched_output[i, "Wikidata"] <- df$item
              } else if (length(df) > 1) {
                final_matched_output[i, "Wikidata"] <- paste0(df$item, collapse = ", ")
              }
            } else if (length(df) == 0) {
              final_matched_output[i, "Wikidata"] <- ""
            }
          }

          # All input markers
          input_markers <- reformatted_data_2()

          input_markers$`All inputted markers` <-
            paste0(input_markers$Marker,
                   " (",
                   input_markers$Sign,
                   ")")

          input_markers <- input_markers %>%
            dplyr::group_by(Cluster) %>%
            dplyr::summarise(
              `Full input marker description` = paste(`All inputted markers`, collapse = ", "),
              `input_num_total` = dplyr::n()
            )

          final_matched_output <-
            merge(final_matched_output,
                  input_markers,
                  by.x = "Cluster",
                  by.y = "Cluster")

          final_matched_output3 <-
            merge(x = final_matched_output,
                  y = final_matched_output2,
                  by = c("CL_term","Cluster"), all.x = TRUE)

          final_matched_output3$`num contradiction(s)`[is.na(final_matched_output3$`num contradiction(s)`)] <- 0

          # Count how many markers are positive and negative
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_pos_input = stringr::str_count(`Full marker description`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_neg_input = stringr::str_count(`Full marker description`, '(Negative)'))

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_pos_cell_type = stringr::str_count(`Full marker description`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_neg_cell_type = stringr::str_count(`Full marker description`, '(Negative)'))

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_pos_match = stringr::str_count(`Matched inputted marker`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_neg_match = stringr::str_count(`Matched inputted marker`, '(Negative)'))

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_pos_contr = stringr::str_count(`Contradiction(s)`, '\\(Positive\\)|\\(High\\)|\\(Low\\)|\\(Low, Positive\\)|\\(High, Positive\\)|\\(Positive, Low\\)|\\(Positive, High\\)|\\(Low, High\\)|\\(High, Low\\)'))
          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(Num_neg_contr = stringr::str_count(`Contradiction(s)`, 'Negative'))

          # TCR

          # Positive
          TCR_pos <- c('TCR \\(POSITIVE\\)', 'TCRAB \\(POSITIVE\\)', 'TCRA \\(POSITIVE\\)', 'TCRB \\(POSITIVE\\)', 'TCRGD \\(POSITIVE\\)', 'TCRG \\(POSITIVE\\)', 'TCRD \\(POSITIVE\\)', 'CD3E \\(POSITIVE\\)',
                       'TCR \\(HIGH\\)', 'TCRAB \\(HIGH\\)', 'TCRA \\(HIGH\\)', 'TCRB \\(HIGH\\)', 'TCRGD \\(HIGH\\)', 'TCRG \\(HIGH\\)', 'TCRD \\(HIGH\\)', 'CD3E \\(HIGH\\)',
                       'TCR \\(LOW\\)', 'TCRAB \\(LOW\\)', 'TCRA \\(LOW\\)', 'TCRB \\(LOW\\)', 'TCRGD \\(LOW\\)', 'TCRG \\(LOW\\)', 'TCRD \\(LOW\\)', 'CD3E \\(LOW\\)',
                       'TCR \\(LOW, POSITIVE\\)', 'TCRAB \\(LOW, POSITIVE\\)', 'TCRA \\(LOW, POSITIVE\\)', 'TCRB \\(LOW, POSITIVE\\)', 'TCRGD \\(LOW, POSITIVE\\)', 'TCRG \\(LOW, POSITIVE\\)', 'TCRD \\(LOW, POSITIVE\\)', 'CD3E \\(LOW, POSITIVE\\)',
                       'TCR \\(POSITIVE, LOW\\)', 'TCRAB \\(POSITIVE, LOW\\)', 'TCRA \\(POSITIVE, LOW\\)', 'TCRB \\(POSITIVE, LOW\\)', 'TCRGD \\(POSITIVE, LOW\\)', 'TCRG \\(POSITIVE, LOW\\)', 'TCRD \\(POSITIVE, LOW\\)', 'CD3E \\(POSITIVE, LOW\\)',
                       'TCR \\(LOW, HIGH\\)', 'TCRAB \\(LOW, HIGH\\)', 'TCRA \\(LOW, HIGH\\)', 'TCRB \\(LOW, HIGH\\)', 'TCRGD \\(LOW, HIGH\\)', 'TCRG \\(LOW, HIGH\\)', 'TCRD \\(LOW, HIGH\\)', 'CD3E \\(LOW, HIGH\\)',
                       'TCR \\(HIGH, POSITIVE\\)', 'TCRAB \\(HIGH, POSITIVE\\)', 'TCRA \\(HIGH, POSITIVE\\)', 'TCRB \\(HIGH, POSITIVE\\)', 'TCRGD \\(HIGH, POSITIVE\\)', 'TCRG \\(HIGH, POSITIVE\\)', 'TCRD \\(HIGH, POSITIVE\\)', 'CD3E \\(HIGH, POSITIVE\\)',
                       'TCR \\(POSITIVE, HIGH\\)', 'TCRAB \\(POSITIVE, HIGH\\)', 'TCRA \\(POSITIVE, HIGH\\)', 'TCRB \\(POSITIVE, HIGH\\)', 'TCRGD \\(POSITIVE, HIGH\\)', 'TCRG \\(POSITIVE, HIGH\\)', 'TCRD \\(POSITIVE, HIGH\\)', 'CD3E \\(POSITIVE, HIGH\\)',
                       'TCR \\(HIGH, LOW\\)', 'TCRAB \\(HIGH, LOW\\)', 'TCRA \\(HIGH, LOW\\)', 'TCRB \\(HIGH, LOW\\)', 'TCRGD \\(HIGH, LOW\\)', 'TCRG \\(HIGH, LOW\\)', 'TCRD \\(HIGH, LOW\\)', 'CD3E \\(HIGH, LOW\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_input_pos = stringr::str_count(toupper(`Full marker description`), paste(TCR_pos, collapse='|')))

          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 1] <- 0
          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 2] <- 1
          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 3] <- 2
          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 4] <- 3
          final_matched_output3$TCR_input_pos[final_matched_output3$TCR_input_pos == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_match_pos = stringr::str_count(toupper(`Matched inputted marker`), paste(TCR_pos, collapse='|')))

          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 1] <- 0
          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 2] <- 1
          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 3] <- 2
          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 4] <- 3
          final_matched_output3$TCR_match_pos[final_matched_output3$TCR_match_pos == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_contradictions_pos = stringr::str_count(toupper(`Contradiction(s)`), paste(TCR_pos, collapse='|')))

          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 1] <- 0
          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 2] <- 1
          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 3] <- 2
          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 4] <- 3
          final_matched_output3$TCR_contradictions_pos[final_matched_output3$TCR_contradictions_pos == 5] <- 4
          final_matched_output3$TCR_contradictions_pos[is.na(final_matched_output3$TCR_contradictions_pos)] <- 0

          TCR_pos <- c('T cell receptor complex \\(Positive\\)', 'alpha-beta T cell receptor complex \\(Positive\\)', 'gamma-delta T cell receptor complex \\(Positive\\)', 'CD3 epsilon \\(Positive\\)',
                       'T cell receptor complex \\(High\\)', 'alpha-beta T cell receptor complex \\(High\\)', 'gamma-delta T cell receptor complex \\(High\\)', 'CD3 epsilon \\(High\\)',
                       'T cell receptor complex \\(Low\\)', 'alpha-beta T cell receptor complex \\(Low\\)', 'gamma-delta T cell receptor complex \\(Low\\)', 'CD3 epsilon \\(Low\\)',
                       'T cell receptor complex \\(Low, Positive\\)', 'alpha-beta T cell receptor complex \\(Low, Positive\\)', 'gamma-delta T cell receptor complex \\(Low, Positive\\)', 'CD3 epsilon \\(Low, Positive\\)',
                       'T cell receptor complex \\(Positive, Low\\)', 'alpha-beta T cell receptor complex \\(Positive, Low\\)', 'gamma-delta T cell receptor complex \\(Positive, Low\\)', 'CD3 epsilon \\(Positive, Low\\)',
                       'T cell receptor complex \\(Low, High\\)', 'alpha-beta T cell receptor complex \\(Low, High\\)', 'gamma-delta T cell receptor complex \\(Low, High\\)', 'CD3 epsilon \\(Low, High\\)',
                       'T cell receptor complex \\(High, Positive\\)', 'alpha-beta T cell receptor complex \\(High, Positive\\)', 'gamma-delta T cell receptor complex \\(High, Positive\\)', 'CD3 epsilon \\(High, Positive\\)',
                       'T cell receptor complex \\(Positive, High\\)', 'alpha-beta T cell receptor complex \\(Positive, High\\)', 'gamma-delta T cell receptor complex \\(Positive, High\\)', 'CD3 epsilon \\(Positive, High\\)',
                       'T cell receptor complex \\(High, Low\\)', 'alpha-beta T cell receptor complex \\(High, Low\\)', 'gamma-delta T cell receptor complex \\(High, Low\\)', 'CD3 epsilon \\(High, Low\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_cell_type_pos = stringr::str_count(`Full marker description`, paste(TCR_pos, collapse='|')))

          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 1] <- 0
          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 2] <- 1
          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 3] <- 2
          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 4] <- 3
          final_matched_output3$TCR_cell_type_pos[final_matched_output3$TCR_cell_type_pos == 5] <- 4

          # Negative

          TCR_negs <- c('TCR \\(NEGATIVE\\)', 'TCRAB \\(NEGATIVE\\)', 'TCRA \\(NEGATIVE\\)', 'TCRB \\(NEGATIVE\\)', 'TCRGD \\(NEGATIVE\\)', 'TCRG \\(NEGATIVE\\)', 'TCRD \\(NEGATIVE\\)', 'CD3E \\(NEGATIVE\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_input_neg = stringr::str_count(toupper(`Full marker description`), paste(TCR_negs, collapse='|')))

          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg_pos == 1] <- 0
          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg == 2] <- 1
          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg == 3] <- 2
          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg == 4] <- 3
          final_matched_output3$TCR_input_neg[final_matched_output3$TCR_input_neg == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_match_neg = stringr::str_count(toupper(`Matched inputted marker`), paste(TCR_negs, collapse='|')))

          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 1] <- 0
          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 2] <- 1
          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 3] <- 2
          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 4] <- 3
          final_matched_output3$TCR_match_neg[final_matched_output3$TCR_match_neg == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_contradictions_neg = stringr::str_count(toupper(`Contradiction(s)`), paste(TCR_negs, collapse='|')))

          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 1] <- 0
          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 2] <- 1
          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 3] <- 2
          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 4] <- 3
          final_matched_output3$TCR_contradictions_neg[final_matched_output3$TCR_contradictions_neg == 5] <- 4
          final_matched_output3$TCR_contradictions_neg[is.na(final_matched_output3$TCR_contradictions_neg)] <- 0

          TCR_negs <- c('T cell receptor complex \\(Negative\\)', 'alpha-beta T cell receptor complex \\(Negative\\)', 'gamma-delta T cell receptor complex \\(Negative\\)', 'CD3 epsilon \\(Negative\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(TCR_cell_type_neg = stringr::str_count(`Full marker description`, paste(TCR_negs, collapse='|')))

          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 1] <- 0
          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 2] <- 1
          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 3] <- 2
          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 4] <- 3
          final_matched_output3$TCR_cell_type_neg[final_matched_output3$TCR_cell_type_neg == 5] <- 4

          final_matched_output3$`total number of markers` <- as.numeric(final_matched_output3$`total number of markers`)
          final_matched_output3$input_num_total <- as.numeric(final_matched_output3$input_num_total)
          final_matched_output3$`num of inputs that match to cell type` <- as.numeric(final_matched_output3$`num of inputs that match to cell type`)
          final_matched_output3$`num contradiction(s)` <- as.numeric(final_matched_output3$`num contradiction(s)`)

          final_matched_output3$`total number of markers` <- final_matched_output3$`total number of markers` - final_matched_output3$TCR_cell_type_pos - final_matched_output3$TCR_cell_type_neg

          final_matched_output3$input_num_total <- final_matched_output3$input_num_total - final_matched_output3$TCR_input_pos - final_matched_output3$TCR_input_neg

          final_matched_output3$`num of inputs that match to cell type` <- final_matched_output3$`num of inputs that match to cell type` - final_matched_output3$TCR_match_neg - final_matched_output3$TCR_match_pos

          final_matched_output3$`num contradiction(s)` <- final_matched_output3$`num contradiction(s)` - final_matched_output3$TCR_contradictions_neg - final_matched_output3$TCR_contradictions_pos

          final_matched_output3$Num_pos_input <- final_matched_output3$Num_pos_input - final_matched_output3$TCR_input_pos
          final_matched_output3$Num_neg_input <- final_matched_output3$Num_neg_input - final_matched_output3$TCR_input_neg

          final_matched_output3$Num_pos_cell_type <- final_matched_output3$Num_pos_cell_type - final_matched_output3$TCR_cell_type_pos
          final_matched_output3$Num_neg_cell_type <- final_matched_output3$Num_neg_cell_type -final_matched_output3$TCR_cell_type_neg

          final_matched_output3$Num_pos_match <- final_matched_output3$Num_pos_match - final_matched_output3$TCR_match_pos
          final_matched_output3$Num_neg_match <- final_matched_output3$Num_neg_match - final_matched_output3$TCR_match_neg

          final_matched_output3$Num_pos_contr <- final_matched_output3$Num_pos_contr - final_matched_output3$TCR_contradictions_pos
          final_matched_output3$Num_neg_contr <- final_matched_output3$Num_neg_contr -final_matched_output3$TCR_contradictions_neg

          # CD8
          # Positive
          CD8_pos <- c('CD8A \\(POSITIVE\\)', 'CD8ALPHABETA \\(POSITIVE\\)',
                       'CD8A \\(HIGH\\)', 'CD8ALPHABETA \\(HIGH\\)',
                       'CD8A \\(LOW\\)', 'CD8ALPHABETA \\(LOW\\)',
                       'CD8A \\(LOW, POSITIVE\\)', 'CD8ALPHABETA \\(LOW, POSITIVE\\)',
                       'CD8A \\(POSITIVE, LOW\\)', 'CD8ALPHABETA \\(POSITIVE, LOW\\)',
                       'CD8A \\(LOW, HIGH\\)', 'CD8ALPHABETA \\(LOW, HIGH\\)',
                       'CD8A \\(HIGH, POSITIVE\\)', 'CD8ALPHABETA \\(HIGH, POSITIVE\\)',
                       'CD8A \\(POSITIVE, HIGH\\)', 'CD8ALPHABETA \\(POSITIVE, HIGH\\)',
                       'CD8A \\(HIGH, LOW\\)', 'CD8ALPHABETA \\(HIGH, LOW\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_input_pos = stringr::str_count(toupper(`Full marker description`), paste(CD8_pos, collapse='|')))

          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 1] <- 0
          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 2] <- 1
          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 3] <- 2
          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 4] <- 3
          final_matched_output3$CD8_input_pos[final_matched_output3$CD8_input_pos == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_match_pos = stringr::str_count(toupper(`Matched inputted marker`), paste(CD8_pos, collapse='|')))

          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 1] <- 0
          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 2] <- 1
          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 3] <- 2
          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 4] <- 3
          final_matched_output3$CD8_match_pos[final_matched_output3$CD8_match_pos == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_contradictions_pos = stringr::str_count(toupper(`Contradiction(s)`), paste(CD8_pos, collapse='|')))

          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 1] <- 0
          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 2] <- 1
          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 3] <- 2
          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 4] <- 3
          final_matched_output3$CD8_contradictions_pos[final_matched_output3$CD8_contradictions_pos == 5] <- 4
          final_matched_output3$CD8_contradictions_pos[is.na(final_matched_output3$CD8_contradictions_pos)] <- 0

          CD8_pos <- c('T-cell surface glycoprotein CD8 alpha chain \\(Positive\\)', 'T cell receptor co-receptor CD8 \\(Positive\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(High\\)', 'T cell receptor co-receptor CD8 \\(High\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Low\\)', 'T cell receptor co-receptor CD8 \\(Low\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Low, Positive\\)', 'T cell receptor co-receptor CD8 \\(Low, Positive\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Positive, Low\\)', 'T cell receptor co-receptor CD8 \\(Positive, Low\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Low, High\\)', 'T cell receptor co-receptor CD8 \\(Low, High\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(High, Positive\\)', 'T cell receptor co-receptor CD8 \\(High, Positive\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(Positive, High\\)', 'T cell receptor co-receptor CD8 \\(Positive, High\\)',
                       'T-cell surface glycoprotein CD8 alpha chain \\(High, Low\\)', 'T cell receptor co-receptor CD8 \\(High, Low\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_cell_type_pos = stringr::str_count(`Full marker description`, paste(CD8_pos, collapse='|')))
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 1] <- 0
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 2] <- 1
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 3] <- 2
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 4] <- 3
          final_matched_output3$CD8_cell_type_pos[final_matched_output3$CD8_cell_type_pos == 5] <- 4

          # Negative

          CD8_negs <- c('CD8A \\(NEGATIVE\\)', 'CD8ALPHABETA \\(NEGATIVE\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_input_neg = stringr::str_count(toupper(`Full marker description`), paste(CD8_negs, collapse='|')))

          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg_pos == 1] <- 0
          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg == 2] <- 1
          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg == 3] <- 2
          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg == 4] <- 3
          final_matched_output3$CD8_input_neg[final_matched_output3$CD8_input_neg == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_match_neg = stringr::str_count(toupper(`Matched inputted marker`), paste(CD8_negs, collapse='|')))

          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 1] <- 0
          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 2] <- 1
          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 3] <- 2
          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 4] <- 3
          final_matched_output3$CD8_match_neg[final_matched_output3$CD8_match_neg == 5] <- 4

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_contradictions_neg = stringr::str_count(toupper(`Contradiction(s)`), paste(CD8_negs, collapse='|')))

          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 1] <- 0
          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 2] <- 1
          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 3] <- 2
          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 4] <- 3
          final_matched_output3$CD8_contradictions_neg[final_matched_output3$CD8_contradictions_neg == 5] <- 4
          final_matched_output3$CD8_contradictions_neg[is.na(final_matched_output3$CD8_contradictions_neg)] <- 0

          CD8_negs <- c('T-cell surface glycoprotein CD8 alpha chain \\(Negative\\)', 'T cell receptor co-receptor CD8 \\(Negative\\)')

          final_matched_output3 <- final_matched_output3 %>% dplyr::mutate(CD8_cell_type_neg = stringr::str_count(`Full marker description`, paste(CD8_negs, collapse='|')))

          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 1] <- 0
          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 2] <- 1
          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 3] <- 2
          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 4] <- 3
          final_matched_output3$CD8_cell_type_neg[final_matched_output3$CD8_cell_type_neg == 5] <- 4

          final_matched_output3$`total number of markers` <- final_matched_output3$`total number of markers` - final_matched_output3$CD8_cell_type_pos - final_matched_output3$CD8_cell_type_neg

          final_matched_output3$input_num_total <- final_matched_output3$input_num_total - final_matched_output3$CD8_input_pos - final_matched_output3$CD8_input_neg

          final_matched_output3$`num of inputs that match to cell type` <- final_matched_output3$`num of inputs that match to cell type` - final_matched_output3$CD8_match_neg - final_matched_output3$CD8_match_pos


          final_matched_output3$Num_pos_input <- final_matched_output3$Num_pos_input - final_matched_output3$CD8_input_pos
          final_matched_output3$Num_neg_input <- final_matched_output3$Num_neg_input - final_matched_output3$CD8_input_neg

          final_matched_output3$Num_pos_cell_type <- final_matched_output3$Num_pos_cell_type - final_matched_output3$CD8_cell_type_pos
          final_matched_output3$Num_neg_cell_type <- final_matched_output3$Num_neg_cell_type -final_matched_output3$CD8_cell_type_neg

          final_matched_output3$Num_pos_match <- final_matched_output3$Num_pos_match - final_matched_output3$CD8_match_pos
          final_matched_output3$Num_neg_match <- final_matched_output3$Num_neg_match - final_matched_output3$CD8_match_neg

          final_matched_output3$Num_pos_contr[is.na(final_matched_output3$Num_pos_contr)] <- 0
          final_matched_output3$Num_neg_contr[is.na(final_matched_output3$Num_neg_contr)] <- 0

          final_matched_output3$Num_pos_contr <- final_matched_output3$Num_pos_contr - final_matched_output3$CD8_contradictions_pos
          final_matched_output3$Num_neg_contr <- final_matched_output3$Num_neg_contr -final_matched_output3$CD8_contradictions_neg

          # Score
          final_matched_output3$Score <-
            ((as.numeric(final_matched_output3$Num_pos_match)) /
               (as.numeric(final_matched_output3$Num_pos_cell_type) + as.numeric(final_matched_output3$Num_pos_input) - as.numeric(final_matched_output3$Num_pos_match))) - ((as.numeric(final_matched_output3$Num_pos_contr) * 0.25) + (as.numeric(final_matched_output3$Num_neg_contr) * 0.25))

          final_matched_output3$Score <- as.numeric(final_matched_output3$Score)

          final_matched_output3$Score <- round(final_matched_output3$Score, 3)

          # Percent match (match markers/markers in cell type)
          final_matched_output3$`Percent match (match markers/markers in cell type)` <-
            (as.numeric(final_matched_output3$`num of inputs that match to cell type`) /
               (as.numeric(final_matched_output3$`total number of markers`)
               )) * 100

          final_matched_output3$`Percent match (match markers/markers in cell type)` <- as.numeric(final_matched_output3$`Percent match (match markers/markers in cell type)`)

          final_matched_output3$`Percent match (match markers/markers in cell type)` <- round(final_matched_output3$`Percent match (match markers/markers in cell type)`, 3)

          # Percent match (match markers/markers in input)
          final_matched_output3$`Percent match (match markers/markers in input)` <-
            ((as.numeric(final_matched_output3$`num of inputs that match to cell type`) /
                (as.numeric(final_matched_output3$input_num_total))) * 100)

          final_matched_output3$`Percent match (match markers/markers in input)` <- as.numeric(final_matched_output3$`Percent match (match markers/markers in input)`)

          final_matched_output3$`Percent match (match markers/markers in input)` <- round(final_matched_output3$`Percent match (match markers/markers in input)`, 3)

          # Get the columns for output, new order
          final_matched_output4 <- final_matched_output3[, c(
            "Cluster",
            "CL_term",
            "Cell_type",
            "Cell_comment",
            "Wikidata",
            "Matched inputted marker",
            "num of inputs that match to cell type",
            "Full input marker description",
            "input_num_total",
            "Full marker description",
            "total number of markers",
            "Contradiction(s)",
            "num contradiction(s)",
            "Percent match (match markers/markers in cell type)",
            "Percent match (match markers/markers in input)",
            "Score")]

          # New column names
          colnames(final_matched_output4)<-
            c(
              "Cluster",
              "Cell ontology ID",
              "Cell type name",
              "Cell type description",
              "Wikidata ID",
              "Matched markers",
              "# matched markers",
              "Inputted markers (also in reference)",
              "# inputted markers (also in reference)",
              "Cell type markers",
              "# cell type markers",
              "Contradictions",
              "# contradictions",
              "% (matched markers/markers in cell type)",
              "% (matched markers/markers in input)",
              "Score")

          # Order by score, get top N scores per cluster
          final_matched_output4 <-  dplyr::tbl_df(final_matched_output4) %>%
            dplyr::group_by(Cluster) %>%
            dplyr::filter(as.integer(ordered(-Score)) %in% 1:input$num_top_matches_2)

          final_matched_output4 <- final_matched_output4[with(final_matched_output4, order(Cluster, -Score)), ]

          return(final_matched_output4)
        } else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    }
  })

  # Step 3, make cell ontology term linkable in RShiny
  return_cell_types_linked_2 <- shiny::eventReactive(input$submit_tab2_step3, {

    final_matched_output <- return_cell_types_2()

    if (!is.null(final_matched_output)) {

      final_matched_output <- final_matched_output[, c(
        "Cluster",
        "Cell ontology ID",
        "Cell type name",
        "Cell type description",
        "Wikidata ID",
        "Matched markers",
        "Inputted markers (also in reference)",
        "Cell type markers",
        "Contradictions",
        "% (matched markers/markers in cell type)",
        "% (matched markers/markers in input)",
        "Score")]

      # Make look nicer (round and add a % symbol)
      final_matched_output$`% (matched markers/markers in cell type)` <-
        paste(round(final_matched_output$`% (matched markers/markers in cell type)`, 1), "%")

      final_matched_output$`% (matched markers/markers in input)` <-
        paste(round(final_matched_output$`% (matched markers/markers in input)`, 1), "%")

      # Just get cell ID to make link in R Shiny look better
      just_IDs <-
        gsub(
          "http://purl.obolibrary.org/obo/",
          "",
          final_matched_output$`Cell ontology ID`
        )

      just_IDs_wiki <-
        gsub(
          "http://www.wikidata.org/entity/",
          "",
          final_matched_output$`Wikidata ID`
        )

      # Remove '<' and '>' from URI
      final_matched_output$`Cell ontology ID` <-
        gsub('<', '', final_matched_output$`Cell ontology ID`)
      final_matched_output$`Cell ontology ID` <-
        gsub('>', '', final_matched_output$`Cell ontology ID`)
      just_IDs <- gsub('>', '', just_IDs)
      just_IDs <- gsub('<', '', just_IDs)

      # Make it so the link is clickable
      for (a in 1:nrow(final_matched_output)) {
        if (final_matched_output$`Cell ontology ID`[a] != "No match found") {
          final_matched_output$`Cell ontology ID`[a] <-
            paste0(
              "<a href='",
              final_matched_output$`Cell ontology ID`[a],
              "' target='_blank'>",
              just_IDs[a],
              "</a>"
            )
        }

        if (final_matched_output$`Wikidata ID`[a] != "") {
          if (lengths(regmatches(final_matched_output$`Wikidata ID`[a], gregexpr(",", final_matched_output$`Wikidata ID`[a]))) == 0) {
            final_matched_output$`Wikidata ID`[a] <- paste0(
              "<a href='",
              final_matched_output$`Wikidata ID`[a],
              "' target='_blank'>",
              just_IDs_wiki[a],
              "</a>"
            )
          } else {
            final_matched_output$`Wikidata ID`[a] <- paste0(
              "<a href='",
              strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][1],
              "' target='_blank'>",
              strsplit(just_IDs_wiki[a], ",")[[1]][1],
              "</a> , <a href='",
              strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][2],
              "' target='_blank'>",
              strsplit(just_IDs_wiki[a], ",")[[1]][2],
              "</a>"
            )
            if (lengths(regmatches(final_matched_output$`Wikidata ID`[a], gregexpr(",", final_matched_output$`Wikidata ID`[a]))) == 2) {
              final_matched_output$`Wikidata ID`[a] <- paste0(final_matched_output$`Wikidata ID`[a], "<a href='",
                                                              strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][3],
                                                              "' target='_blank'>",
                                                              strsplit(just_IDs_wiki[a], ",")[[1]][3],
                                                              "</a>")
            } else if (lengths(regmatches(final_matched_output$`Wikidata ID`[a], gregexpr(",", final_matched_output$`Wikidata ID`[a]))) == 3) {
              final_matched_output$`Wikidata ID`[a] <- paste0(
                "<a href='",
                strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][3],
                "' target='_blank'>",
                strsplit(just_IDs_wiki[a], ",")[[1]][3],
                "</a> , <a href='",
                strsplit(final_matched_output$`Wikidata ID`[a], ",")[[1]][4],
                "' target='_blank'>",
                strsplit(just_IDs_wiki[a], ",")[[1]][4],
                "</a>"
              )
            }
          }
        }
      }

   #   final_matched_output$Cluster <- as.numeric(final_matched_output$Cluster)

      return(final_matched_output)
    } else {
      return(NULL)
    }
  })

  # Step 3, get how many clusters are in final output for filtering purposes
  output$specific_cluster_2 <- shiny::renderUI({
    shiny::req(input$submit_tab2_step3)

    final_df <- return_cell_types_linked_2()
    b <- sort(unique(as.character(final_df$Cluster)))
    if (!is.null(final_df) & length(unique(final_df$Cluster)) > 1) {
      shiny::selectInput(inputId = "show_specific_cluster_2",
                         label = "Filter by cluster:  ",
                         selected = sort(unique(as.character(final_df$Cluster))),
                         choices = sort(unique(as.character(final_df$Cluster))),
                         multiple = TRUE)
    }
  })

  # Step 3, if multiple clusters, change dataframe if filtered by cluster
  df_subsettable_2 <- reactive({
  shiny::req(input$submit_tab2_step3)
   
    df <- return_cell_types_linked_2()

    if (!is.null(df) & length(unique(df$Cluster)) > 1) {
      df2 <- df[df$Cluster %in% input$show_specific_cluster_2,]
      return(df2)
    } else {
      return(df)
    }
  })

  # Step 3, make datatable of matched cell types
  output$table_cell_types_2 <- DT::renderDT({
     req(input$submit_tab2_step3)

    df <- df_subsettable_2()

    if (is.data.frame(df) == FALSE) {
    } else if (nrow(df) == 0) {
    } else {
      DT::datatable(
        df,
        selection = 'none',
        rownames = FALSE,
        escape = FALSE,
        options = list(pageLength = 25)
      )
    }
  }, server=FALSE)

  # Step 3, changes when the 'save and continue' button is hit if using internal references
  # Output final dataframe, show buttons once that dataframe is made
  shiny::observeEvent(input$submit_tab2_step3, {

    shiny.destroy::removeOutput("help_tab2_step2")
    shiny.destroy::removeOutput("overview_tab2_step2")
    shiny.destroy::removeOutput("table_protein_synonym_2")
    shiny.destroy::removeOutput("delete_PRO_synonym_2")
    shiny.destroy::removeOutput("download_PRO_GO_matches_2")
    shinyjs::show("overview_tab2_step3")

    if (!is.null(return_cell_types_linked_2())) {
      shiny.destroy::removeOutput("all_sidebar_tab2_step3")
      shinyjs::show("all_sidebar_tab2_step4")
      shinyjs::show("help_tab2_step3")
      shinyjs::show("download_table_cell_types_2")
    } else {
      shinyjs::show("no_matched_cell_types_2")
    }
  })

  # Step 3/Uploaded reference, download final results
  output$download_table_cell_types_2 <- shiny::downloadHandler(
    # File name
    filename = function() {
      if (input$input_type == 'csv_upload') {
        paste0(gsub(".csv|.xls|.xlsx", "", input$csv_upload),
               "-matched-cell-types-", Sys.Date(), ".csv")
      } else {
        paste0("typed_markers_matched_cell_types-", Sys.Date(), ".csv")
      }
    },

    # Write to csv
    content = function(file) {
      utils::write.csv(return_cell_types_2(),  file, row.names = FALSE)
    }
  )

  # Reset everything
  shiny::observeEvent(input$reset_tab2_step4, {
    session$reload()
  })
}


