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


server <- function(input, output, session) {
  
  datafile <- userFileServer("datafile")
  
  samples <- reactive({
    datafile() %>% 
      select('Sample Name') %>% 
      distinct()
  })
  
  genes <- reactive({
    datafile() %>% 
      select('Target Name') %>% 
      distinct()
  })
  
  # Mandatory dynamic inputs
  
  ### 1. Samples
  output$sample_selector <- renderUI({
    pickerInput(inputId = "selected_samples",
                label = "Samples",
                choices = as.character(samples()$`Sample Name`),
                options = list(`actions-box` = TRUE), 
                multiple = TRUE)
  })
  
  ### 2. Genes of interest
  output$gene_selector <- renderUI({
    pickerInput(inputId = "selected_genes",
                label = "Genes of interest",
                choices = as.character(genes()$`Target Name`),
                options = list(`actions-box` = TRUE), 
                multiple = TRUE)
  })
  
  ### 3. Dynamic housekeeping gene selector
  output$housekeeping_selector <- renderUI({
    pickerInput(inputId = "housekeeping_gene", 
                label = "Housekeeping gene", 
                choices = as.character(genes()$`Target Name`))
  })
  
  # Optional dynamic inputs
  
  # If DDCt is selected
  # User defines controls
  output$controls <- renderUI({
    if(input$delta == 'delta2'){
      numSamples <- as.integer(length(input$selected_samples))
      lapply(1:numSamples, function(i) {
        pickerInput(inputId = paste("ctrl",input$selected_samples[i]),
                    label = paste("Control for ", input$selected_samples[i]),
                    choices = input$selected_samples)
      })
    }
  })
  
  # If Group analysis is selected
  # User defines number of groups 
  output$num_groups <- renderUI({
    if(input$groups_switch == TRUE){
      numericInput("num_groups", "Number of groups", value = 1, min = 1, max = 10)
    }
  })
  
  # If Group analysis is selected
  # User defines members of each group
  output$groups_members <- renderUI({
    if(input$groups_switch == TRUE & !is.null(input$num_groups)){
      lapply(1:input$num_groups, function(i) {
        pickerInput(inputId = paste("group",i),
                    label = paste("Samples in group ", i),
                    choices = input$selected_samples,
                    options = list(`actions-box` = TRUE), 
                    multiple = TRUE)
      })
    }
  })
  
  # If Group analysis is selected
  # User defines labels of each group
  output$label_groups <- renderUI({
    if(input$groups_switch == TRUE & !is.null(input$num_groups)){
      lapply(1:input$num_groups, function(i) {
        textInput(inputId = paste("label_group",i),
                  label = paste("Label for group", i),
                  placeholder = paste(i),
                  value = i)
      })
    }
  })
  
  # If Group analysis is selected
  # Do we compute mean of groups?
  output$groups_means <- renderUI({
    if(input$groups_switch == TRUE){
      checkboxInput(inputId = "means_check",
                    label = "Display means",
                    value = TRUE)
    }
  })
  
  # # User defines labels of each group
  # output$groups_labels <- renderUI({
  #   if(input$groups_switch == TRUE & !is.null(input$num_groups)){
  #     lapply(1:input$num_groups, function(i) {
  #       textInput(inputId = paste("grouplabel",i),
  #                 label = paste("Group ", i, " label"),
  #                 placeholder = i, value = i)
  #     })
  #   }
  # })
  
  # Run analysis
  output$go_button <- renderUI({
    actionButton("go", "Run", icon = icon("angle-double-right"))
  })
  
  output$plot_dl <- renderUI({
    downloadButton("downloadPlot", "Download plot")
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { 
      paste(Sys.Date(), "_plot", '.png', sep='') 
    },
    content = function(file) {
      ggsave(plot = pcr_plot(), device = "png", file)
    }
  )
  
  output$data_dl <- renderUI({
    req(input$go)
    downloadButton("downloadData", "Download data")
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_data",".xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(dct(), file)
    }
  )
  
  
  
  # MAIN DATA WRANGLING FUNCTION
  # 0. For each sample, calculates DCt / 2^-DCt 
  # 1. If Delta Delta Ct is required : compute 2^-DDCt
  # 2. IF group analysis is selected : compute group data
  dct <- eventReactive(input$go, {
    
    # 0. ##########################
    
    # Initialize dataframe
    df <- data.frame(character(),character(),double(),double(), double())
    
    # Loop through genes in every sample
    for (sample in input$selected_samples) {
      for (gene in input$selected_genes) {
        # Compute deltaCt
        ct_g <- datafile() %>% filter(`Target Name` == gene & `Sample Name` == sample) %>% use_series(`Ct`)
        ct_r <- datafile() %>% filter(`Target Name` == input$housekeeping_gene & `Sample Name` == sample) %>% use_series(`Ct`)
        dct <- ct_g - ct_r
        edct <- 2^(-1*(dct))
        df <- rbind(df, c(sample,gene,as.double(ct_g)[1],dct[1],edct[1]))
      }
    }
    
    names(df) = c('Sample', 'Gene', 'Ct', 'DeltaCt', 'DeltaCt_exp')
    
    # 1. ##########################
    
    # If DeltaDeltaCt required, add new columns
    # The first added column contains controls for each sample (for verification)
    # The second column is the actual DeltaDeltaCt
    if(input$delta == "delta2") {
      
      # This is an intermediate dataframe
      # for housing the Sample, Gene, Control, DeltaDeltaCt and DeltaDeltaCt_exp values
      # Before merging with main dataframe
      ctrl_df <- data.frame(character(),character(), character(), double())
      
      # Loop through genes in every sample
      for (sample in input$selected_samples) {
        ctrl <- input[[paste("ctrl", sample)]]
        for (gene in input$selected_genes) {
          sampledct <- as.numeric(df %>% filter(Sample == sample, Gene == gene) %>% use_series("DeltaCt"))
          ctrldct <- as.numeric(df %>% filter(Sample == ctrl, Gene == gene) %>% use_series("DeltaCt"))
          # Compute DeltaDeltaCt
          ddct <- as.numeric(sampledct-ctrldct)
          # Compute 2^(-DeltaDeltaCt)
          ddct_exp <- 2^(-1*ddct)
          # Add to the intermediate dataframe
          ctrl_df <- rbind(ctrl_df, c(sample, gene, ctrl, ddct, ddct_exp))
        }
      }
      
      names(ctrl_df) = c('Sample', 'Gene', 'Control', 'DeltaDeltaCt', 'DeltaDeltaCt_exp')
      
      # Merge results to the main dataframe
      df$Control <- ctrl_df$Control
      df$DeltaDeltaCt <- ctrl_df$DeltaDeltaCt
      df$DeltaDeltaCt_exp <- ctrl_df$DeltaDeltaCt_exp
    }
    
    # 2. ##########################
    
    if(input$groups_switch == TRUE) {
      
      df[,"Group"] = NA
      
      for(i in 1:input$num_groups) 
        
      {
        for(s in input[[paste("group", i)]])
        {
          df <- df %>% mutate(Group = replace(Group, Sample == s, input[[paste("label_group", i)]]))
          
        }
      }
    }
    
    df <- df %>% na.omit()
    return(df)
  })
  
  # This function draws the actual plots
  pcr_plot <- reactive({
    req(input$go)
    
    d <- dct()
    ly <- ""
    lx <- ""
    
    if(input$delta == "delta2") {
      # Remove controls (will be replaced by horizontal line at 1.0)
      d <- d %>% filter(DeltaDeltaCt != 0)
      ly <- "Relative expression"
      
      if(input$groups_switch == TRUE){
        lx <- "Group"
        
        pcr_plot <- ggplot(data = d) +
          geom_point(aes(x = as.factor(Group), y = as.numeric(DeltaDeltaCt_exp))) + 
          geom_hline(yintercept = 1, linetype = "dashed") 
        
        if(input$means_check == TRUE) {
          d_means <- d %>% 
            group_by(Gene, Group) %>%
            summarise(mDeltaDeltaCt_exp = mean(as.numeric(DeltaDeltaCt_exp)), sd = sd(DeltaDeltaCt_exp))

          pcr_plot <- pcr_plot + 
            geom_point(data = d_means, 
                       aes(x = as.factor(Group), 
                           y = mDeltaDeltaCt_exp), 
                       color = "gray", shape = 17, size = 2,
                       show.legend = FALSE) +
            geom_errorbar(data = d_means, 
                          aes(x = as.factor(Group), 
                              y = mDeltaDeltaCt_exp, 
                              ymin = mDeltaDeltaCt_exp - sd, 
                              ymax = mDeltaDeltaCt_exp + sd),
                          size = 0.3, color = "gray", width = 0.05)
        }
        
      }
      
      else {
        lx <- "Sample"
        pcr_plot <- ggplot(data = d) +
          geom_point(aes(x = Sample, y = as.numeric(DeltaDeltaCt_exp))) + 
          geom_hline(yintercept = 1, linetype = "dashed")
        
      }
    }
    
    else {
      ly <- "Expression"
      if(input$groups_switch == TRUE){
        
        lx <- "Group"
        pcr_plot <- ggplot(data = d) +
          geom_point(aes(x = as.factor(Group), y = as.numeric(DeltaCt_exp)))
        
        if(input$means_check == TRUE) {
          d_means <- d %>% 
            group_by(Group, Gene) %>%
            summarise(mDeltaCt_exp = mean(as.numeric(DeltaCt_exp)), sd = sd(DeltaCt_exp))
          pcr_plot <- pcr_plot + 
            geom_point(data = d_means,
                       aes(x = as.factor(Group),
                           y = mDeltaCt_exp),
                       color = "gray", shape = 17, size = 2,
                       show.legend = FALSE) +
            geom_errorbar(data = d_means, 
                          aes(x = as.factor(Group), 
                              y = mDeltaCt_exp, 
                              ymin = mDeltaCt_exp - sd, 
                              ymax = mDeltaCt_exp + sd),
                          size = 0.3, color = "gray", width = 0.05)
        }
        
      }
      else{
        lx <- "Sample"
        pcr_plot <- ggplot(data = d) +
          geom_point(aes(x = Sample, y = as.numeric(DeltaCt_exp)))
      }
    }
    
    if(input$log_scale) {
      pcr_plot <- pcr_plot + 
        scale_y_log10() + 
        annotation_logticks(base = 10, sides = "l")  
    }
    
    pcr_plot <-  pcr_plot + theme_prism() + xlab(lx) + ylab(ly) +
    facet_wrap(vars(as.factor(Gene)), ncol = 2, scales = "free") +
    theme(aspect.ratio = 1)
    
    
    return(pcr_plot)
  })
  
  
  # These functions output the data 
  output$table_out <- renderDataTable(dct())
  
  output$plot_out <- renderPlot({pcr_plot()}, height = function(){
    req(input$go)
    isolate(400*ceiling(length(input$selected_genes)*0.5))
    }
    )
  
  output$groups_out <- renderDataTable(groups())
  
  output$info <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
    }
    xy_range_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
             " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
    }
    
    paste0(
      "Pointer coordinates: ", xy_str(input$plot_hover)
    )
  })
  
}
