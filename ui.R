library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(tibble)
library(magrittr)
library(shinyWidgets)
library(ggprism)

ui <- fluidPage(id = "page",
                
                fluidRow(
                  # App title
                  column(width = 4, 
                         titlePanel("qPCR shiny v1.0")),
                  
                  
                  column(width = 3, div(style = "padding: 15px 0px; margin:00%"),
                         uiOutput(outputId = "go_button")),
                         
                  
                  column(width = 3, div(style = "padding: 0px 0px; margin:0%"),
                         fileInput(inputId = "file", "label" = "Upload dataset",
                                   multiple = FALSE, 
                                   accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv"), 
                                   width = NULL))
                ),
                
                
                # Sidebar layout with input and output definitions
                fluidRow(
                  
                  #Sidebar panel for inputs
                  column(width = 2,
                         # Are we comparing groups?
                         switchInput(
                           inputId = "groups_switch",
                           label = "Groups", 
                           value = FALSE
                         ),
                         
                         # Input : number of groups
                         conditionalPanel(condition = "input.groups_switch == TRUE",
                                          uiOutput(outputId = "num_groups")),
                         
                         # Input : sample selector
                         uiOutput(outputId = "sample_selector"),
                         
                         # Input : GOI selector
                         uiOutput(outputId = "gene_selector"),
                         
                         # Input: housekeeping gene selector
                         uiOutput(outputId = "housekeeping_selector"),
                         
                         # Input : if group analysis is chosen, define members of groups 
                         uiOutput(outputId = "groups_members")
                  ),
                  
                  column(width = 2,
                         
                         # Input: Delta or Delta Delta Ct?
                         radioButtons(inputId = "delta", 
                                      label = "Analysis", 
                                      choices = c("Delta Ct" = "delta1", "Delta Delta Ct" = "delta2")),
                         
                         # Input : If delta delta ct, define control samples for each sample
                         conditionalPanel(condition = "input.delta == 'delta2'",
                                          uiOutput(outputId = "controls")),
                         
                         # Input : if group analysis is chosen, choose if you want means computed
                         conditionalPanel(condition = "input.groups_switch == TRUE",
                                          uiOutput(outputId = "groups_means"))
                         
                  ),
                  
                  # Main panel for displaying outputs
                  column(width = 6, 
                         tabsetPanel(
                           # Plot output
                           tabPanel("Plot", 
                                    plotOutput("plot_out", 
                                               width = "100%",
                                               hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce")),
                                    verbatimTextOutput("info"),
                                    uiOutput("plot_dl")),
                           # Input : if group analysis is chosen, define names of groups 
                           # uiOutput(outputId = "groups_labels")),
                           # DeltaCt +/- DeltaDeltaCt Data table
                           tabPanel("Data table", dataTableOutput("table_out"),
                                    uiOutput("data_dl"))
                         )
                  )
                )
)
