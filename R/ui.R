#' Shiny app user interface
#'
#' @name ui_app
#' @param req Required
#'
#' @import dplyr
#' @import DT
#' @import httr
#' @import magrittr
#' @import reshape2
#' @import scales
#' @import shiny
#' @import shinybusy
#' @import shinyjs
#' @import shinyvalidate
#' @import stringi
#' @import stringr
#' @import tidyr
#' @importFrom magrittr %>%
#'
#'
source("./R/get-input_choices.R")
source("./R/get-species.R")
source("./R/endpoints.R")

#source("get-input_choices.R")
#source("get-species.R")
#source("endpoints.R")

ui_app <- function(req){
  fluidPage(

    # Use Shinyjs package
    shinyjs::useShinyjs(),

    # CSS, design of app
    tags$head(tags$style(
      #      type="text/css",
      #      ".shiny-output-error {
      #      visibility: hidden;
      #      }",
      #      ".shiny-output-error:before {
      #      visibility: hidden;
      #      }"
      #    ),
      HTML(
        '
        .tabbable > .nav > li > a {
        background-color: #f1efef;
        color:black;
        }
        .tabbable > .nav > li[class=active] > a {
           background-color: #f6eb84;
           color:black;
        }
       #inital_alert {
         background-color:#b24bb4;
       }
       #sidebar_tab1_step1, #sidebar_tab1_step2,
       #sidebar_tab1_step3, #sidebar_tab1_step4,
       #sidebar_tab1_step5, #sidebar_tab1_step6,
       #sidebar_tab2_step1, #sidebar_tab2_step2,
       #sidebar_tab2_step3, #sidebar_tab2_step4 {
         background-color: #18a1b2;
         color: white;
       }
       .progress-bar {
         background-color: #f6eb84;
         color: black;
       }
       .js-irs-0 .irs-single,
       .js-irs-0 .irs-bar-edge,
       .js-irs-0 .irs-bar {
         background: #f6eb84;
         color: black;
       }
       .js-irs-1 .irs-single,
       .js-irs-1 .irs-bar-edge,
       .js-irs-1 .irs-bar {
         background: #a4cca3;
         color: black;
       }
        .js-irs-2 .irs-single,
       .js-irs-2 .irs-bar-edge,
       .js-irs-2 .irs-bar {
         background: #a4cca3;
         color: black;
       }
        .js-irs-3 .irs-single,
       .js-irs-3 .irs-bar-edge,
       .js-irs-3 .irs-bar {
         background: #a4cca3;
         color: black;
       }
        .js-irs-4 .irs-single,
       .js-irs-4 .irs-bar-edge,
       .js-irs-4 .irs-bar {
         background: #a4cca3;
         color: black;
       }
        .js-irs-5 .irs-single,
       .js-irs-5 .irs-bar-edge,
       .js-irs-5 .irs-bar {
         background: #a4cca3;
         color: black;
       }
        .js-irs-6 .irs-single,
       .js-irs-6 .irs-bar-edge,
       .js-irs-6 .irs-bar {
         background: #a4cca3;
         color: black;
       }
        .js-irs-7 .irs-single,
       .js-irs-7 .irs-bar-edge,
       .js-irs-7 .irs-bar {
         background: #a4cca3;
         color: black;
       }
        .js-irs-8 .irs-single,
       .js-irs-8 .irs-bar-edge,
       .js-irs-8 .irs-bar {
         background: #f6eb84;
         color: black;
       }
       body {
         background-color: white;
         }
       #download_example, #download_example_expression_cluster,
       #download_example_marker_cell_type, #download_example_marker_cell_type_1,
       #download_example_marker_cell_type_2 {
         color: #f6eb84;
       }
       #submit_tab1_step1, #submit_tab1_step2,
      #submit_tab1_step2_5, #submit_tab1_step3,
      #submit_tab1_step4, #submit_tab1_step5,
      #reset_tab1_step1, #reset_tab1_step2,
       #reset_tab1_step3, #reset_tab1_step4,
       #reset_tab1_step5, #reset_tab1_step6,
       #submit_tab2_step1, #submit_tab2_step2,
       #submit_tab2_step3, #reset_tab2_step1,
       #reset_tab2_step2, #reset_tab2_step3,
       #reset_tab2_step4 {
         background-color: #f6eb84;
         color: black;
         font-size: 12px;
       }
       #download_PRO_GO_matches_1,
       #download_table_cell_types_1, #download_table_cell_types,
       #delete_PRO_synonym_1, #delete_marker,
       #download_marker_diff_eq_csv, #download_uploaded_ref_table_cell_types_1,
       #download_uploaded_ref_table_cell_types_2, #delete_PRO_synonym_2,
       #download_PRO_GO_matches_2, #download_table_cell_types_2 {
         background-color: #f6eb84;
         color: black;
       }
       #choose_cluster label{
       display: table-cell;
       text-align: middle;
       vertical-align: middle;
       }
       #choose_cluster .form-group {
       display: table-row;
       }
       #choose_cluster {
       display: inline-block;
       float: right;
       }
       #help_tab1_step1, #help_tab1_step2,
       #help_tab1_step3, #help_tab1_step4,
       #help_tab1_step5, #help_tab1_step6,
       #help_tab2_step1, #help_tab2_step2,
       #help_tab2_step3, #help_uploaded_ref_1,
       #help_uploaded_ref_2 {
       float: right;
       }
       #help_pregating, #help_pregating_2, #help_transform_option,
       #help_seed_type_1, #help_random_type_1,
       #help_input_type_1, #help_input_type_2,
       #help_marker_diff_eq_option, #help_pos_neg,
       #help_all_pos, #help_all_null,
       #help_ref_type, #help_ref_type_2,
       #help_species,#help_species_2,
       #help_ontology_type_1, #help_ontology_type_2,
       #help_ontology_type_2,
       #help_low_upload, #help_high_upload,
       #help_low_2, #help_high_2 {
       padding: 3.5px;
       font-size: 110%;
       font-weight: bold;
       margin-left: 7px;
       }
       h1 {
         color: #b24bb4;
         text-align: center;
         font-size: 55px;
       }
       h2 {
         color: #18a1b2;
         text-align: center;
         font-size: 30px;
       }
       h3 {
         color: #b24bb4;
         text-align: center;
         font-size: 25px;
         font-style: italic;
       }
       h4 {
         color: #b24bb4;
         font-size: 20px;
         text-align: center;
       }
        h5 {
         color: #f1efef;
         text-align: center;
         font-size: 17px;
         font-weight: bold;
        }
        h6 {
         color: #f1efef;
         text-align: center;
         font-size: 15px;
         font-weight: bold;
        }
       p {
         color: black;
          font-size: 15px;
          text-align: left;
         }
      .shiny-output-error {
      visibility: hidden;
      }
      .shiny-output-error:before {
      visibility: hidden;
      }'
      )
    )),

    # Title of app
    h1(strong("Descriptive Cell Type Naming")),

    shiny::tabsetPanel(

      #################################################################
      ##                            Tab 1:                           ##
      ##                        User interface                       ##
      ##                    Input expression data                    ##
      #################################################################

      shiny::tabPanel("Input expression data",
                      br(),
                      div(id = "all_sidebar_tab1_step1", shiny::sidebarPanel(id = "sidebar_tab1_step1",


                                                                             h5('Expression data'),
                                                                             shiny::fileInput(
                                                                               inputId = "upload_expression_csv",
                                                                               label = HTML(
                                                                                 paste(
                                                                                   "Upload post-clustering expression",  shiny::actionButton("help_input_type_1", "?"),"<br /> CSV or XLSX file <br />",
                                                                                   shiny::downloadLink("download_example_expression_cluster", label = "(download example)"),
                                                                                   sep = " "
                                                                                 )
                                                                               ),
                                                                               multiple = FALSE,
                                                                               # accept = c('text/csv',
                                                                               #            'text/comma-separated-values',
                                                                               #            '.csv'),
                                                                               width = NULL,
                                                                               buttonLabel = "Browse...",
                                                                               placeholder = "No file selected"
                                                                             ),
                                                                             shiny::textInput(
                                                                               inputId = "type_MG_markers",
                                                                               label = HTML(
                                                                                 paste("Type any marker(s) used for pre-gating,",  shiny::actionButton("help_pregating", "?"), "<br /> using commas to seperate <br /> (e.g. CD8+, CD4-)")),
                                                                               value = "",
                                                                               width = NULL
                                                                             ),
                                                                             # Choose if you want to  transform the data
                                                                             shiny::radioButtons(
                                                                               "transform_option",
                                                                               label = HTML(paste("Do you want to arcsinh transform the data?", shiny::actionButton("help_transform_option", "?"))),
                                                                               choices = input_choices_transform,
                                                                               selected = "transform_no"
                                                                             ),

                                                                             # If choosing transform,
                                                                             shiny::conditionalPanel(
                                                                               condition = "input.transform_option == 'transform_yes'",
                                                                               # Set cofactor
                                                                               shiny::numericInput(
                                                                                 inputId = "type_cofactor",
                                                                                 label = HTML(
                                                                                   "Set the cofactor"
                                                                                 ),
                                                                                 value = "",
                                                                                 width = NULL
                                                                               )
                                                                             ),

                                                                             # Checkbox, check if you want to set a seed
                                                                             shiny::checkboxInput(
                                                                               inputId = "random_1",
                                                                               label = HTML(paste("Randomly downsample expression data", shiny::actionButton("help_random_type_1", "?"))),
                                                                               value = FALSE),

                                                                             # If downsampling
                                                                             shiny::conditionalPanel(
                                                                               condition = "input.random_1 == 1",
                                                                               shiny::numericInput(
                                                                                 inputId = "num_random_1",
                                                                                 label = HTML(
                                                                                   "Number of cells"
                                                                                 ),
                                                                                 value = 50000,
                                                                                 width = NULL
                                                                               )
                                                                             ),

                                                                             # Checkbox, check if you want to set a seed
                                                                             shiny::checkboxInput(
                                                                               inputId = "seed_1",
                                                                               label = HTML(paste("Set seed (to ensure reproducible results)", shiny::actionButton("help_seed_type_1", "?"))),
                                                                               value = FALSE),

                                                                             # If choosing set seed
                                                                             shiny::conditionalPanel(
                                                                               condition = "input.seed_1 == 1",
                                                                               shiny::numericInput(
                                                                                 inputId = "num_seed_1",
                                                                                 label = HTML(
                                                                                   "Seed number"
                                                                                 ),
                                                                                 value = 123,
                                                                                 width = NULL
                                                                               )
                                                                             ),

                                                                             # Submit button
                                                                             shiny::actionButton("submit_tab1_step1", label = "Submit"),

                                                                             # Reset button
                                                                             shiny::actionButton("reset_tab1_step1", label = "Reset")

                      )),

                      shinyjs::hidden(div(id = "all_sidebar_tab1_step2", shiny::sidebarPanel(id = "sidebar_tab1_step2",

                                                                                             h5('Median Difference Equation'),

                                                                                             # Choose if you want to use suggested Median Difference Equation parameters
                                                                                             shiny::radioButtons(
                                                                                               "marker_diff_eq_option",
                                                                                               label = HTML(paste("Do you want to use default Median Difference Equation parameters?", shiny::actionButton("help_marker_diff_eq_option", "?"))),
                                                                                               choices = input_choices_marker_diff_eq,
                                                                                               selected = "default_marker_diff_eq"
                                                                                             ),

                                                                                             # If choosing set Median Difference Equation parameters,
                                                                                             shiny::conditionalPanel(
                                                                                               condition = "input.marker_diff_eq_option == 'choose_marker_diff_eq'",

                                                                                               h6(HTML(paste('Marker designated as positive or negative',
                                                                                                             shiny::actionButton("help_pos_neg", "?")))),

                                                                                               # Set Positive cutoff
                                                                                               shiny::numericInput(
                                                                                                 inputId = "type_positive_cutoff",
                                                                                                 label = "Set the positive cutoff",
                                                                                                 value = 0.25,
                                                                                                 width = NULL
                                                                                               ),

                                                                                               # Set Negative cutoff
                                                                                               shiny::numericInput(
                                                                                                 inputId = "type_negative_cutoff",
                                                                                                 label = "Set the negative cutoff",
                                                                                                 value = -0.75,
                                                                                                 width = NULL
                                                                                               ),

                                                                                               h6(HTML(paste('Marker designated as completely positive', shiny::actionButton("help_all_pos", "?")))),

                                                                                               # Set all median positive cutoff
                                                                                               shiny::numericInput(
                                                                                                 inputId = "type_med_positive_cutoff",
                                                                                                 label = HTML(paste("Set minimum median cutoff  <br /> (minimum must be above this value)  <br />")),
                                                                                                 value = 0.2,
                                                                                                 width = NULL
                                                                                               ),

                                                                                               # Set max median positive cutoff
                                                                                               shiny::numericInput(
                                                                                                 inputId = "type_max_med_positive_cutoff",
                                                                                                 label = HTML(paste("Set maximum median cutoff  <br /> (maximum must be above this value)  <br />")),
                                                                                                 value = 1,
                                                                                                 width = NULL
                                                                                               ),

                                                                                               # Set standard deviation positive cutoff
                                                                                               shiny::numericInput(
                                                                                                 inputId = "type_stdr_positive_cutoff",
                                                                                                 label = HTML(paste("Set standard deviation cutoff  <br /> (standard deviation must be below this value)  <br />")),
                                                                                                 value = 0.7,
                                                                                                 width = NULL
                                                                                               ),

                                                                                               h6(HTML(paste('Marker designated as completely null', shiny::actionButton("help_all_null", "?")))),

                                                                                               # Set all median null cutoff
                                                                                               shiny::numericInput(
                                                                                                 inputId = "type_med_negative_cutoff",
                                                                                                 label = HTML(paste("Set minimum median cutoff <br />  (minimum must be below this value) <br /> ")),
                                                                                                 value = 0.2,
                                                                                                 width = NULL
                                                                                               ),

                                                                                               # Set maximum median null cutoff
                                                                                               shiny::numericInput(
                                                                                                 inputId = "type_max_med_negative_cutoff",
                                                                                                 label = HTML(paste("Set maximum median cutoff <br /> (maximum must be below this value) <br /> ")),
                                                                                                 value = 0.7,
                                                                                                 width = NULL
                                                                                               ),

                                                                                               # Set standard deviation null cutoff
                                                                                               shiny::numericInput(
                                                                                                 inputId = "type_stdr_negative_cutoff",
                                                                                                 label = HTML(paste("Set standard deviation cutoff <br /> (standard deviation must be below this value) <br />")),
                                                                                                 value = 0.5,
                                                                                                 width = NULL
                                                                                               ),
                                                                                             ),

                                                                                             # Submit button
                                                                                             shinyjs::hidden(shiny::actionButton("submit_tab1_step2", label = "Save and continue")),

                                                                                             shinyjs::hidden(shiny::actionButton("submit_tab1_step2_5", label = "Save and continue")),

                                                                                             # Reset button
                                                                                             shiny::actionButton("reset_tab1_step2", label = "Reset to beginning")
                      ))),

                      shinyjs::hidden(div(id = "all_sidebar_tab1_step3", shiny::sidebarPanel(id = "sidebar_tab1_step3",

                                                                                             h5('Cell type matching'),

                                                                                             # Choose whether to use internal or uploaded reference
                                                                                             shiny::radioButtons(
                                                                                               "reference_type_1",
                                                                                               label = HTML(paste("Reference type", shiny::actionButton("help_ref_type", "?"))),
                                                                                               choices = input_choices_ref_1,
                                                                                               selected = "internal_ref_1"
                                                                                             ),

                                                                                             # If choosing internal reference, specify the species
                                                                                             shiny::conditionalPanel(
                                                                                               condition = "input.reference_type_1 == 'internal_ref_1'",
                                                                                               shiny::selectInput(
                                                                                                 "species_choice_1",
                                                                                                 label = HTML(paste("Select the species", shiny::actionButton("help_species", "?"))),
                                                                                                 choices = species_in_PRO$full_name,
                                                                                                 multiple = FALSE
                                                                                               )
                                                                                             ),

                                                                                             # If choosing internal reference, specify if you want to use CL, pCL, or both
                                                                                             shiny::conditionalPanel(
                                                                                               condition = "input.reference_type_1 == 'internal_ref_1'",
                                                                                               shiny::checkboxGroupInput(
                                                                                                 "ontology_type_1",
                                                                                                 label = HTML(paste("Select cell type ontology(s)", shiny::actionButton("help_ontology_type_1", "?"))),
                                                                                                 choices = input_choices_onto_1,
                                                                                                 selected = c("internal_CL_1", "internal_pCL_1")
                                                                                               )
                                                                                             ),

                                                                                             # If choosing uploaded reference, give the user a file upload box
                                                                                             shiny::conditionalPanel(
                                                                                               condition = "input.reference_type_1 == 'upload_ref_1'",
                                                                                               shiny::fileInput(
                                                                                                 inputId = "upload_ref_csv_1",
                                                                                                 label = HTML(
                                                                                                   paste(
                                                                                                     "Upload reference marker-cell type <br /> CSV or XLSX file <br />",
                                                                                                     shiny::downloadLink("download_example_marker_cell_type_1", label = "(download example)"),
                                                                                                     sep = " "
                                                                                                   )
                                                                                                 ),
                                                                                                 multiple = FALSE,
                                                                                                 # accept = c('text/csv',
                                                                                                 #            'text/comma-separated-values',
                                                                                                 #            '.csv'),
                                                                                                 width = NULL,
                                                                                                 buttonLabel = "Browse...",
                                                                                                 placeholder = "No file selected"
                                                                                               )
                                                                                             ),


                                                                                             # Submit button
                                                                                             shiny::actionButton("submit_tab1_step3", label = "Save and continue"),

                                                                                             # Reset button
                                                                                             shiny::actionButton("reset_tab1_step3", label = "Reset to beginning")
                      ))),

                      shinyjs::hidden(div(id = "all_sidebar_tab1_step4", shiny::sidebarPanel(id = "sidebar_tab1_step4",


                                                                                             h5('Cell type matching'),

                                                                                             # Submit button
                                                                                             shiny::actionButton("submit_tab1_step4", label = "Save and continue"),

                                                                                             # Reset button
                                                                                             shiny::actionButton("reset_tab1_step4", label = "Reset to beginning")

                      ))),

                      shinyjs::hidden(div(id = "all_sidebar_tab1_step5", shiny::sidebarPanel(id = "sidebar_tab1_step5",


                                                                                             h5('Cell type matching'),

                                                                                             # Choose how many results (top N matches) you want to return per cluster
                                                                                             shiny::sliderInput(
                                                                                               "num_top_matches_1",
                                                                                               label = "Return top N matches per cluster (including ties)",
                                                                                               min = 1,
                                                                                               max = 5,
                                                                                               value = 3
                                                                                             ),

                                                                                             # Submit button
                                                                                             shiny::actionButton("submit_tab1_step5", label = "Save and continue"),

                                                                                             # Reset button
                                                                                             shiny::actionButton("reset_tab1_step5", label = "Reset to beginning")

                      ))),


                      shinyjs::hidden(div(id = "all_sidebar_tab1_step6", shiny::sidebarPanel(id = "sidebar_tab1_step6",

                                                                                             # Reset button
                                                                                             shiny::actionButton("reset_tab1_step6", label = "Reset to beginning")

                      ))),

                      shiny::mainPanel(align="center",

                                       # Loading spinner
                                       shinybusy::add_busy_spinner(spin = "fading-circle", height = "140px",
                                                                   width = "140px", position = "top-left",  margins = c("50vh", "60vw"), color = "#b24bb4"),


                                       # Instructions/information about the app
                                       tags$div(
                                         id = "app_overview_1",

                                         h2("Overview"),
                                         p("The full application takes post-clustering flow or mass cytometry expression data as the input, returns marker definitions per cluster (Part 1), standardizes marker names (Part 2), and finally matches to cell type names (Part 3)."),
                                         p("If the user prefers to directly input marker descriptions (e.g. CD4+, CD8-) and skip Part 1, they can do so through the ", strong("Input Marker Descriptors"), "tab"),
                                         p("More information on this tool, including examples files and a comprehensive user guide, can be found on the", a(href = 'https://github.com/AmandaRT18/work-in-progress', target="_blank", " GitHub", .noWS = "outside"),".", .noWS = c("after-begin", "before-end")),


                                         img(src='Overview-input-expression-values.png', align = "center", width = "100%")
                                       ),

                                       #~#~#~#~#~#~#~#~#~# #~#~#~#~#~#~#~#~# #~#~#~#~#~#~#~#~#~#
                                       #~#~#~#~#~#~#~#~#~# Input expression #~#~#~#~#~#~#~#~#~#
                                       #~#~#~#~#~#~#~#~#~# #~#~#~#~#~#~#~#~# #~#~#~#~#~#~#~#~#~#

                                       #~#~# #~#~#~# #~#~#
                                       #~#~# Step 1 #~#~#
                                       #~#~# #~#~#~# #~#~#

                                       # Step 1, uploaded expression data, help button
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_tab1_step1", "Help"
                                       )),

                                       # Step 1, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_tab1_step1",
                                         h2("Part 1: Get Marker Definitions", style = 'color:black; font-size:35px'),
                                         h2("Step 1: Select Markers Used in Analysis"),
                                         br()
                                       )),

                                       # Step 1, error message if the uploaded input is not a csv, xlsx, or xls file
                                       shiny::uiOutput("ui_error_check_input_file_1"), #works

                                       # Uploaded reference, error message if the column names for the uploaded query file are incorrect
                                       shiny::uiOutput("ui_error_column_names_input_1"), #works

                                       # Step 1, show marker names
                                       DT::DTOutput("expression_data"),

                                       # Step 1, 'delete' button
                                       shinyjs::hidden(shiny::actionButton(
                                         "delete_marker", "Delete marker(s)"
                                       )),

                                       #~#~#~#~#~#~#~#~#~# #~#~#~#~#~#~#~#~# #~#~#~#~#~#~#~#~#~#
                                       #~#~#~#~#~# Median Difference Equation #~#~#~#~#~#
                                       #~#~#~#~#~#~#~#~#~# #~#~#~#~#~#~#~#~# #~#~#~#~#~#~#~#~#~#

                                       #~#~# #~#~#~# #~#~#
                                       #~#~# Step 2/3 #~#~#
                                       #~#~# #~#~#~# #~#~#

                                       # Step 2, Median Difference Equation help button
                                       shiny::conditionalPanel(
                                         condition = "output.marker_diff_eq_heatmap_binary_ex_plot",

                                         shiny::actionButton(
                                           "help_tab1_step2", "Help"
                                         )
                                       ),

                                       # Step 2, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_tab1_step2",
                                         h2("Part 1: Get Marker Definitions", style = 'color:black; font-size:35px'),
                                         h2("Step 2: Adjust Median Difference Equation Parameters"),
                                         br()
                                       )),

                                       # Step 3, help
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_tab1_step3", "Help"
                                       )),

                                       # Step 3, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_tab1_step3",
                                         h2("Part 1: Get Marker Definitions", style = 'color:black; font-size:35px'),
                                         h2("Step 3: Override Median Difference Equation Results"),
                                         br()
                                       )),

                                       # # Step 3, dataframe

                                       shiny::conditionalPanel(
                                         condition = "output.marker_diff_eq_heatmap_binary_ex_plot2",

                                         shiny::tagList(shiny::uiOutput("cutoff_parameters"),

                                                        # Step 3, dataframe
                                                        shinyjs::hidden(DT::DTOutput("marker_diff_eq_results_editable")),

                                                        shiny::uiOutput("marker_diff_eq_heatmap_binary_ex_plot2"),

                                                        tags$div(id = "spaces1",
                                                                 br(),
                                                                 br()),

                                                        tags$div(id ="marker_diff_eq_heatmap_raw_row_size2" , style="display: inline-block;",

                                                                 shiny::sliderInput(
                                                                   "marker_diff_eq_heatmap_raw_row_size",
                                                                   label = "Change text size of heatmap rows",
                                                                   min = 1,
                                                                   max = 30,
                                                                   value = 15
                                                                 )),

                                                        tags$div(id = "marker_diff_eq_heatmap_raw_column_size2", style="display: inline-block;",
                                                                 shiny::sliderInput(
                                                                   "marker_diff_eq_heatmap_raw_column_size",
                                                                   label = "Change text size of heatmap columns",
                                                                   min = 1,
                                                                   max = 30,
                                                                   value = 15
                                                                 )),

                                                        tags$div(id = "heatmap_height_size2", style="display: inline-block;",
                                                                 shiny::sliderInput(
                                                                   "heatmap_height_size",
                                                                   label = "Change height of heatmap plots",
                                                                   min = 100,
                                                                   max = 1200,
                                                                   value = 800
                                                                 )),

                                                        tags$div(id = "spaces2",
                                                                 br()
                                                        ),

                                                        tags$div(id = "density_row_size2", style="display: inline-block;",
                                                                 shiny::sliderInput(
                                                                   "density_row_size",
                                                                   label = "Change text size of density plot (x-axis)",
                                                                   min = 1,
                                                                   max = 30,
                                                                   value = 17
                                                                 )),

                                                        tags$div(id = "density_height_size2", style="display: inline-block;",
                                                                 shiny::sliderInput(
                                                                   "density_height_size",
                                                                   label = "Change height of density plot",
                                                                   min = 200,
                                                                   max = 1500,
                                                                   value = 1000
                                                                 )),

                                                        tags$div(id = "density_facet2", style="display: inline-block;",
                                                                 shiny::sliderInput(
                                                                   "density_facet",
                                                                   label = "Change number of splits (rows) in density plot",
                                                                   min = 1,
                                                                   max = 5,
                                                                   value = 3
                                                                 )),

                                                        tags$div(id = "spaces3",
                                                                 br()
                                                        ))),

                                       shinyjs::hidden(shiny::downloadButton("download_marker_diff_eq_csv", label = "Download marker definitions (.csv)", suspendWhenHidden = FALSE)),

                                       #~#~# #~#~#~# #~#~#
                                       #~#~# Step 4 #~#~#
                                       #~#~# #~#~#~# #~#~#

                                       #~#~#~#~#~#  Uploaded reference #~#~#~#~#~#

                                       # Uploaded reference, help
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_uploaded_ref_1", "Help"
                                       )),

                                       # Uploaded reference, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_uploaded_ref_1",
                                         h2("Part 2: Get Cell Type Names", style = 'color:black; font-size:35px'),
                                         h2("Step 4: Match Marker Patterns to Cell Types"),
                                         br()
                                       )),

                                       # Uploaded reference, error message if the uploaded reference is not a csv, xlsx, or xls file
                                       shiny::uiOutput("ui_error_check_reference_file_1"), #works

                                       # Uploaded reference, error message if the column names for the uploaded reference are incorrect
                                       shiny::uiOutput("ui_error_column_names_ref_1"), #works

                                       # Uploaded reference, error message if there is a marker with no sign in the uploaded reference
                                       shiny::uiOutput("ui_warning_no_sign_ref_1"), #works

                                       # Uploaded reference, get statement indicating if/which markers are not in uploaded reference
                                       shiny::uiOutput("ui_warning_unmatched_markers_uploaded_ref_1"), #works

                                       # Uploaded reference,, filter dataframe by cluster
                                       tags$div(id = "choose_cluster_uploaded_1",
                                                shiny::uiOutput("specific_cluster_uploaded_1")), #works

                                       #  Uploaded reference, results table
                                       DT::DTOutput("uploaded_ref_matches_1_dt"), #works

                                       # Uploaded reference, download final dataframe
                                       shinyjs::hidden(shiny::downloadButton("download_uploaded_ref_table_cell_types_1", "Download")),

                                       #~#~#~#~#~#  Internal reference #~#~#~#~#~#

                                       # Step 4, help
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_tab1_step4", "Help"
                                       )),

                                       # Step 4, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_tab1_step4",
                                         h2("Part 2: Get Standardized Marker Names", style = 'color:black; font-size:35px'),
                                         h2("Step 4: Match Markers to PRO/GO Terms"),
                                         br()
                                       )),

                                       # Step 4, error message if there is a marker duplicated within a cluster
                                       shiny::uiOutput("ui_marker_sign_error_1"),

                                       # Step 4, message unmatched markers in PRO
                                       shiny::uiOutput("ui_warning_unmatched_markers_1"),

                                       # Step 4, dataframe
                                       DT::DTOutput("table_new_format_1"),

                                       #~#~# #~#~#~# #~#~#
                                       #~#~# Step 5 #~#~#
                                       #~#~# #~#~#~# #~#~#

                                       # Step 5, help
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_tab1_step5", "Help"
                                       )),

                                       # Step 5, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_tab1_step5",
                                         h2("Part 2: Get Standardized Marker Names", style = 'color:black; font-size:35px'),
                                         h2("Step 5: Confirm Matched Marker Terms"),
                                         br()
                                       )),

                                       DT::DTOutput("table_protein_synonym_1"),

                                       # Step 5, 'delete' button
                                       shinyjs::hidden(shiny::actionButton(
                                         "delete_PRO_synonym_1", "Delete marker entry"
                                       )),

                                       # Step 5, download dataframe showing how inputted markers matched to PRO/GO terms
                                       shinyjs::hidden(shiny::downloadButton("download_PRO_GO_matches_1", "Download")),


                                       #~#~# #~#~#~# #~#~#
                                       #~#~# Step 6 #~#~#
                                       #~#~# #~#~#~# #~#~#

                                       # Step 6, help
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_tab1_step6", "Help"
                                       )),

                                       # Step 6, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_tab1_step6",
                                         h2("Part 3: Get Cell Type Names", style = 'color:black; font-size:35px'),
                                         h2("Step 6: Match Marker Patterns to Cell Types"),
                                         br()
                                       )),

                                       # Step 6, filter dataframe by cluster
                                       tags$div(id = "choose_cluster_1",
                                                shiny::uiOutput("specific_cluster_1")),

                                       DT::DTOutput("table_cell_types_1"),

                                       # Step 6, message if there are no cell type matches
                                       shinyjs::hidden(tags$div(
                                         id = "no_matched_cell_types_1",
                                         h4("No cell types matched with your query")
                                       )),

                                       # Step 6/uploaded reference, download final dataframe
                                       shinyjs::hidden(shiny::downloadButton("download_table_cell_types_1", "Download"))
                      )

      ),

      ##################################################################
      ##                            Tab 2:                            ##
      ##                        User interface                        ##
      ##                   Input marker descriptors                   ##
      ##################################################################

      shiny::tabPanel("Input marker descriptors",


                      div(id = "all_sidebar_tab2_step1", shiny::sidebarPanel(id = "sidebar_tab2_step1",

                                                                             h5('Upload Marker Definitions'),

                                                                             # Choose input type (type in box or uploaded file)
                                                                             shiny::radioButtons(
                                                                               "input_type",
                                                                               label = HTML(
                                                                                 paste("Input type", shiny::actionButton("help_input_type_2", "?"))),
                                                                               choices = input_choices,
                                                                               selected = "type_input"
                                                                             ),

                                                                             # If choosing type input, give text input box
                                                                             shiny::conditionalPanel(
                                                                               condition = "input.input_type == 'type_input'",
                                                                               shiny::textInput(
                                                                                 inputId = "type_input",
                                                                                 label = HTML(
                                                                                   "Type each marker name and sign, <br /> using commas to seperate <br /> (e.g. CCR3+, CD203c++, CD19-)"
                                                                                 ),
                                                                                 value = "",
                                                                                 width = NULL
                                                                               )
                                                                             ),

                                                                             # If choosing file upload, give the user a file upload box
                                                                             shiny::conditionalPanel(
                                                                               condition = "input.input_type == 'csv_upload'",
                                                                               shiny::fileInput(
                                                                                 inputId = "csv_upload",
                                                                                 label = HTML(paste(
                                                                                   "Upload input CSV or XLSX file <br />",
                                                                                   shiny::downloadLink("download_example", label = "(download example)"),
                                                                                   sep = " "
                                                                                 )),
                                                                                 multiple = FALSE,
                                                                                 # accept = c('text/csv',
                                                                                 #            'text/comma-separated-values',
                                                                                 #            '.csv',
                                                                                 #            ".xlsx",
                                                                                 #            ".xls"),
                                                                                 width = NULL,
                                                                                 buttonLabel = "Browse...",
                                                                                 placeholder = "No file selected"
                                                                               )
                                                                             ),

                                                                             shiny::textInput(
                                                                               inputId = "type_MG_markers_2",
                                                                               label = HTML(
                                                                                 paste("Type any marker(s) used for pre-gating,",  shiny::actionButton("help_pregating_2", "?"), "<br /> using commas to seperate <br /> (e.g. CD8+, CD4-)")),
                                                                               value = "",
                                                                               width = NULL
                                                                             ),

                                                                             h5('Cell Type Reference'),

                                                                             # Choose whether to use internal or uploaded reference
                                                                             shiny::radioButtons(
                                                                               "reference_type_2",
                                                                               label = HTML(paste("Reference type", shiny::actionButton("help_ref_type_2", "?"))),
                                                                               choices = input_choices_ref_2,
                                                                               selected = "internal_ref_2"
                                                                             ),

                                                                             # If choosing internal reference, specify the species
                                                                             shiny::conditionalPanel(
                                                                               condition = "input.reference_type_2 == 'internal_ref_2'",
                                                                               shiny::selectInput(
                                                                                 "species_choice_2",
                                                                                 label = HTML(paste("Select the species", shiny::actionButton("help_species_2", "?"))),
                                                                                 choices = species_in_PRO$full_name,
                                                                                 multiple = FALSE
                                                                               )
                                                                             ),

                                                                             # If choosing internal reference, specify if you want to use CL, pCL, or both
                                                                             shiny::conditionalPanel(
                                                                               condition = "input.reference_type_2 == 'internal_ref_2'",
                                                                               shiny::checkboxGroupInput(
                                                                                 "ontology_type_2",
                                                                                 label = HTML(paste("Select cell type ontology(s)", shiny::actionButton("help_ontology_type_2", "?"))),
                                                                                 choices = input_choices_onto_2,
                                                                                 selected = c("internal_CL_2", "internal_pCL_2")
                                                                               )
                                                                             ),

                                                                             # If choosing uploaded reference, give the user a file upload box
                                                                             shiny::conditionalPanel(
                                                                               condition = "input.reference_type_2 == 'upload_ref_2'",
                                                                               shiny::fileInput(
                                                                                 inputId = "upload_ref_csv_2",
                                                                                 label = HTML(
                                                                                   paste(
                                                                                     "Upload reference marker-cell type <br /> CSV or XLSX file <br />",
                                                                                     shiny::downloadLink("download_example_marker_cell_type_2", label = "(download example)"),
                                                                                     sep = " "
                                                                                   )
                                                                                 ),
                                                                                 multiple = FALSE,
                                                                                 # accept = c('text/csv',
                                                                                 #            'text/comma-separated-values',
                                                                                 #            '.csv',
                                                                                 #            ".xlsx",
                                                                                 #            ".xls"),
                                                                                 width = NULL,
                                                                                 buttonLabel = "Browse...",
                                                                                 placeholder = "No file selected"
                                                                               ),

                                                                               # Checkbox, check if positive marker signs should be considered if using 'high' or 'low'
                                                                               shiny::radioButtons(
                                                                                 "low_option_upload",
                                                                                 label = HTML(paste("Select the interpretation of 'low'", shiny::actionButton("help_low_upload", "?"))),
                                                                                 choices = input_choices_low_upload,
                                                                                 selected = "low_positive_upload"
                                                                               ),

                                                                               shiny::radioButtons(
                                                                                 "high_option_upload",
                                                                                 label = HTML(paste("Select the interpretation of 'high'", shiny::actionButton("help_high_upload", "?"))),
                                                                                 choices = input_choices_high_upload,
                                                                                 selected = "high_positive_upload"
                                                                               ),
                                                                             ),

                                                                             # Submit button
                                                                             shiny::actionButton("submit_tab2_step1", label = "Save and continue"),

                                                                             # Reset button
                                                                             shiny::actionButton("reset_tab2_step1", label = "Reset")
                      )),

                      shinyjs::hidden(div(id = "all_sidebar_tab2_step2", shiny::sidebarPanel(id = "sidebar_tab2_step2",

                                                                                             h5('Cell type matching'),

                                                                                             # Submit button
                                                                                             shiny::actionButton("submit_tab2_step2", label = "Save and continue"),

                                                                                             # Reset button
                                                                                             shiny::actionButton("reset_tab2_step2", label = "Reset to beginning")

                      ))),

                      shinyjs::hidden(div(id = "all_sidebar_tab2_step3", shiny::sidebarPanel(id = "sidebar_tab2_step3",


                                                                                             h5('Cell type matching'),

                                                                                             # Choose how many results (top N matches) you want to return per cluster
                                                                                             shiny::sliderInput(
                                                                                               "num_top_matches_2",
                                                                                               label = "Return top N matches per cluster (including ties)",
                                                                                               min = 1,
                                                                                               max = 5,
                                                                                               value = 3
                                                                                             ),

                                                                                             # Checkbox, check if positive marker signs should be considered if using 'high' or 'low'
                                                                                             shiny::radioButtons(
                                                                                               "low_option_2",
                                                                                               label = HTML(paste("Select the interpretation of 'low'", shiny::actionButton("help_low_2", "?"))),
                                                                                               choices = input_choices_low_2,
                                                                                               selected = "low_positive_2"
                                                                                             ),

                                                                                             shiny::radioButtons(
                                                                                               "high_option_2",
                                                                                               label = HTML(paste("Select the interpretation of 'high'", shiny::actionButton("help_high_2", "?"))),
                                                                                               choices = input_choices_high_2,
                                                                                               selected = "high_positive_2"
                                                                                             ),

                                                                                             # Submit button
                                                                                             shiny::actionButton("submit_tab2_step3", label = "Save and continue"),

                                                                                             # Reset button
                                                                                             shiny::actionButton("reset_tab2_step3", label = "Reset to beginning")

                      ))),

                      shinyjs::hidden(div(id = "all_sidebar_tab2_step4", shiny::sidebarPanel(id = "sidebar_tab2_step4",


                                                                                             h5('Cell type matching'),

                                                                                             # Reset button
                                                                                             shiny::actionButton("reset_tab2_step4", label = "Reset to beginning")

                      ))),




                      shiny::mainPanel(align="center",

                                       # Loading spinner
                                       shinybusy::add_busy_spinner(spin = "fading-circle", height = "140px",
                                                                   width = "140px", position = "top-left",  margins = c("50vh", "60vw"), color = "#b24bb4"),


                                       # Instructions/information about the app
                                       tags$div(
                                         id = "app_overview_2",

                                         h2("Overview"),
                                         p("This part of the application uses marker descriptions (e.g. CD4+, CD8-) as the input, standardizes marker names (Part 1) and matches to cell type names (Part 2)."),
                                         p("If the user prefers to directly input post-clustered expression data and add an additional Part, they can do so through the ", strong("Input Expression Data"), "tab"),
                                         p("More information on this tool, including examples files and a comprehensive user guide, can be found on the", a(href = 'https://github.com/AmandaRT18/work-in-progress', target="_blank", " GitHub", .noWS = "outside"),".", .noWS = c("after-begin", "before-end")),

                                         img(src='Overview-input-marker-definitions.png', align = "center", width = "100%")

                                       ),

                                       #~#~#~#~#~#  Uploaded reference #~#~#~#~#~#

                                       # Uploaded reference, help
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_uploaded_ref_2", "Help"
                                       )),

                                       # Uploaded reference, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_uploaded_ref_2",
                                         h2("Match Marker Patterns to Cell Types"),
                                         br()
                                       )),

                                       # Uploaded reference, error message if the uploaded reference is not a csv, xlsx, or xls file
                                       shiny::uiOutput("ui_error_check_reference_file_2"), #works

                                       # Uploaded reference, error message if the column names for the uploaded reference are incorrect
                                       shiny::uiOutput("ui_error_column_names_ref_2"), #works

                                       # Uploaded reference, error message if there is a marker with no sign in the uploaded reference
                                       shiny::uiOutput("ui_warning_no_sign_ref_2"), #works

                                       # Uploaded reference, error message if there is a marker with no sign in the inputted marker pattern
                                       shiny::uiOutput("ui_warning_no_sign_2"), #works

                                       # Uploaded reference, get statement indicating if/which markers are not in uploaded reference
                                       shiny::uiOutput("ui_warning_unmatched_markers_uploaded_ref_2"), #works

                                       # Uploaded reference,, filter dataframe by cluster
                                       tags$div(id = "choose_cluster_uploaded_2",
                                                shiny::uiOutput("specific_cluster_uploaded_2")),

                                       #  Uploaded reference, results table
                                       DT::DTOutput("uploaded_ref_matches_2_dt"), #works

                                       # Uploaded reference, download final dataframe
                                       shinyjs::hidden(shiny::downloadButton("download_uploaded_ref_table_cell_types_2", "Download")), #works

                                       #~#~#~#~#~#  Internal reference #~#~#~#~#~#

                                       # Step 1, help
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_tab2_step1", "Help"
                                       )),

                                       # Step 1, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_tab2_step1",
                                         h2("Part 1: Get Standardized Marker Names", style = 'color:black; font-size:35px'),
                                         h2("Step 1: Match Markers to PRO/GO Terms"),
                                         br()
                                       )),

                                       # Step 1, error message if the uploaded input is not a csv, xlsx, or xls file
                                       shiny::uiOutput("ui_error_check_input_file_2"), #works?

                                       # Uploaded reference, error message if the column names for the uploaded query file are incorrect
                                       shiny::uiOutput("ui_error_column_names_input_2"), #works

                                       # Step 1, error message if there is a marker duplicated within a cluster
                                       shiny::uiOutput("ui_marker_sign_error_2"), #works

                                       # Step 1, error message if there is a marker with no sign in the inputted marker pattern
                                       shiny::uiOutput("ui_warning_no_sign_2"), # works

                                       # Step 1, message unmatched markers in PRO
                                       shiny::uiOutput("ui_warning_unmatched_markers_2"),

                                       DT::DTOutput("table_new_format_2"),

                                       # Step 2, help
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_tab2_step2", "Help"
                                       )),

                                       # Step 2, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_tab2_step2",
                                         h2("Part 1: Get Standardized Marker Names", style = 'color:black; font-size:35px'),
                                         h2("Step 2: Confirm Matched Marker Terms"),
                                         br()
                                       )),

                                       # Step 2, dataframe
                                       DT::DTOutput("table_protein_synonym_2"),

                                       # Step 2, 'delete' button
                                       shinyjs::hidden(shiny::actionButton(
                                         "delete_PRO_synonym_2", "Delete marker entry"
                                       )),

                                       # Step 2, download dataframe showing how inputted markers matched to PRO/GO terms
                                       shinyjs::hidden(shiny::downloadButton("download_PRO_GO_matches_2", "Download")),

                                       # Step 3, help
                                       shinyjs::hidden(shiny::actionButton(
                                         "help_tab2_step3", "Help"
                                       )),

                                       # Step 3, title
                                       shinyjs::hidden(tags$div(
                                         id = "overview_tab2_step3",
                                         h2("Part 2: Get Cell Type Names", style = 'color:black; font-size:35px'),
                                         h2("Step 3: Match Marker Patterns to Cell Types"),
                                         br()
                                       )),

                                       # Step 3, filter dataframe by cluster
                                       tags$div(id = "choose_cluster_2",
                                                shiny::uiOutput("specific_cluster_2")),

                                       # Step 3, final dataframe
                                       # shiny::conditionalPanel(
                                       #   condition = "input.submit_tab2_step3",
                                       #   DT::DTOutput("table_cell_types_2")),

                                       DT::DTOutput("table_cell_types_2"),

                                       # Step 3, message if there are no cell type matches
                                       shinyjs::hidden(tags$div(
                                         id = "no_matched_cell_types_2",
                                         h4("No cell types matched with your query")
                                       )),

                                       # Step 3/uploaded reference, download final dataframe
                                       shinyjs::hidden(shiny::downloadButton("download_table_cell_types_2", "Download"))
                      )
      )
    )
  )
}
