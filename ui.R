library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(tibble)
library(magrittr)
library(shinyWidgets)
library(ggprism)
library(shinyBS)
library(readxl)
library(writexl)
library(data.table)
library(scales)

ui <- fluidPage(id = "page",
                
                fluidRow(
                  # App title
                  column(width = 4, 
                         titlePanel("shinyPCR")),
                  
                  
                  column(width = 5, div(style = "padding: 15px 0px; margin:00%"),
                         uiOutput(outputId = "go_button")),
                         
                  
                  column(width = 2, div(style = "padding: 0px 0px; margin:0%"),
                         userFileUI("datafile")),
                  
                  bsButton("file_help", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                  
                  bsTooltip(id = "file_help", 
                            "Select one or multiple .xls files (directly exported from QuantStudio Flex)", 
                            # or a single .csv file with 3 columns, titled Sample Name, Target Name and Ct mean", 
                            trigger = "hover")
                ),
                
                
                # Sidebar layout with input and output definitions
                fluidRow(
                  
                  #Sidebar panel for inputs
                  column(width = 2,
                         
                         # Input : sample selector
                         uiOutput(outputId = "sample_selector"),
                         
                         # Input : GOI selector
                         uiOutput(outputId = "gene_selector"),
                         
                         # Input: housekeeping gene selector
                         uiOutput(outputId = "housekeeping_selector"),
                         
                         # Are we comparing groups?
                         switchInput(
                           inputId = "groups_switch",
                           label = "Groups", 
                           value = FALSE
                         ),
                         
                         # Input : number of groups
                         conditionalPanel(condition = "input.groups_switch == TRUE",
                                          uiOutput(outputId = "num_groups")),
                         
                         # Input : if group analysis is chosen, define members of groups 
                         uiOutput(outputId = "groups_members"),
                         
                         # Input : if group analysis is chosen, define labels of groups and draw means
                         uiOutput(outputId = "label_groups")
                  ),
                  
                  column(width = 2,
                         
                         # Input: Delta or Delta Delta Ct?
                         radioButtons(inputId = "delta", 
                                      label = "Analysis", 
                                      choices = c("Delta Ct" = "delta1", "Delta Delta Ct" = "delta2")),
                         
                         # Should we use a logarithmic scale?
                         checkboxInput(inputId = "log_scale",
                                       label = "Log scale",
                                       value = FALSE),
                         
                         # If Groups : Should we display means or error bars?
                         uiOutput(outputId = "groups_means"),
                         
                         # Input : If delta delta ct, define control samples for each sample
                         conditionalPanel(condition = "input.delta == 'delta2'",
                                          uiOutput(outputId = "controls")),
                         
                         
                  ),
                  
                  # Main panel for displaying outputs
                  column(width = 7, 
                         tags$head(tags$script(
                           'var dimension = [0, 0];
              $(document).on("shiny:connected", function(e) {
                  dimension[0] = window.innerWidth;
                  dimension[1] = window.innerHeight;
                  Shiny.onInputChange("dimension", dimension);
              });
              $(window).resize(function(e) {
                  dimension[0] = window.innerWidth;
                  dimension[1] = window.innerHeight;
                  Shiny.onInputChange("dimension", dimension);
              });'
                         )),
                         tabsetPanel(
                           # Plot output
                           tabPanel("Plot", 
                                    verbatimTextOutput("info"),
                                    uiOutput("plot_dl"),
                                    plotOutput("plot_out", 
                                               width = "100%",
                                               hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce"))
                                   ),

                           # DeltaCt +/- DeltaDeltaCt Data table
                           tabPanel("Data table", dataTableOutput("table_out"),
                                    uiOutput("data_dl"))
                           
                         )
                  )
                )
)
