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

# Module UI function
userFileUI <- function(id, label = "Upload dataset") {
  # `NS(id)` returns a namespace function, which was save as `ns` and will
  # invoke later.

  ns <- NS(id)
  
  
  fileInput(ns("file"), label,
            multiple = TRUE, 
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv", 
              "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", 	
              "application/vnd.ms-excel"))
}

# Module server function
userFileServer <- function(id) {
  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      # The selected file, if any
      userFile <- reactive({
        # If no file is selected, don't do anything
        validate(need(input$file, message = FALSE))
        input$file
      })
      

      # The user's data is parsed into a data frame
      dataframe <- reactive({
        # Either a single CSV with three columns (prepared by the user)
        if (any(userFile()$type == c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv"))) 
          {
          read_csv2(userFile()$datapath,
                    col_types = c(col_character(),
                                  col_character(),
                                  col_number()))
          }
        
        # Or one or multiple xls files (created by QuantStudio(TM) 6 Flex)
        else 
          {
          rbindlist(lapply(userFile()$datapath,read_excel, sheet = "Results", range = "A41:AE137"), fill = TRUE) %>%
            select(`Sample Name`, `Target Name`, `Ct Mean`) %>% 
            distinct() %>% 
            na.omit() %>% 
            arrange(`Sample Name`)
            }
      
      })
  
      # 
      
  # samples <- dataframe() %>% select('Sample Name') %>% distinct()
  # genes <- dataframe() %>% select('Target Name') %>% distinct()

  # We can run observers in here if we want to
  observe({
    msg <- sprintf("Uploaded : %s", userFile()$name)
    cat(msg, "\n")
  })
  
  # Return the reactive that yields the data frame
  return(dataframe)
    }
  )
}
