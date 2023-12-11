#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(bslib)
library(tidyverse)
library(ggplot2)
library(tools)
library(colourpicker) # you might need to install this
library(DT)
library(matrixStats)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  titlePanel("BF591 Final Project"),
  HTML(paste(rep("<p>On this webiste you can run analysis related to Huntington's Disease</p>"))),
  HTML(paste(rep("<p>Developed by Christian Landaverde</p>"))),
    tabsetPanel(type = "tabs",
      tabPanel("Samples",
               HTML(paste(rep('<h4>Sample Information Exploration</h4>'))),
                sidebarLayout(
                  sidebarPanel(
                    fileInput("samplesFile", paste0("Load sample information matrix CSV"), accept = c('text/csv', 'text/comma-separated-values',
                                                                                                       'text/tab-separated-values','text/plain', 'csv','tsv')),
                ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("Summary", DTOutput("Ssummary")),
                    tabPanel("Metadata", DTOutput("Stable")),
                    tabPanel("Plots", plotOutput("Shist"), plotOutput("Shist2"), plotOutput("Shist3"))
                    )
                  )
                )
        ),
        tabPanel("Counts",
                 sidebarLayout(
                   sidebarPanel(
                     fileInput("countsFile", paste0("Load normalized counts matrix CSV"), accept = c('text/csv', 'text/comma-separated-values',
                                                                                                         'text/tab-separated-values','text/plain', 'csv','tsv')),
                     sliderInput("countsVariance", min = 0, max = 100, value = 80, label = "Select the minimum percentile  variance of genes"),
                     sliderInput("countsNonzero", min = 0, max = 69, value = 5, label = "Select the minimum number of non-zero samples"),
                     actionButton("countsSubmit", "display", icon("star"))
                 ),
                 mainPanel(
                   tabsetPanel(
                     
                     tabPanel("Summary", DTOutput("Csummary")),
                     tabPanel("Diagnostic Plots", plotOutput("CDiagplot1"), plotOutput("CDiagplot2")),
                     tabPanel("Heatmap"),
                     tabPanel("PCA")
                    )
                  )
                )
        ),
        tabPanel("Diff. Expr.",
                 sidebarLayout(
                   sidebarPanel(
                     fileInput("DiffFile", paste0("Load differential expression results"), accept = c('text/csv', 'text/comma-separated-values',
                                                                                                  'text/tab-separated-values','text/plain', 'csv','tsv')),
                     HTML(paste(rep('<p>A volcano plot can be generated with <b>"log2 fold"</b> change on the x-axis and <b>"p-adjusted"</b> on the y-axis.</p>'), collapse = "")),
                     radioButtons("DiffButtonx", "Choose the column for the x-axis",
                                  choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                     radioButtons("DiffButtony", "Choose the column for the y-axis",
                                  choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                     colourInput("DiffBase_col", "Base point color", "#22577A"),
                     colourInput("DiffHighlight_col", "Highlight point color", "#FFCF56"),
                     sliderInput("DiffPadj_coloring", min = -50, max = 0, value = -10,
                                 label = "Select magnitude of the p adjusted coloring:"),
                     actionButton("DiffSubmit", "display", icon("star"))
                     ),
                   
                   mainPanel(
                     tabsetPanel(
                       type = "tabs",
                       tabPanel("Table", DTOutput("Dtable")),
                       tabPanel("Volcano plot", plotOutput("Dvolcano")),
                       tabPanel("Plot table", DTOutput("DPlottable"))
                   )
                 )

                 )
        ),
        tabPanel("Indiv. Gene Expr.")
      )
    )


# Define server logic required to draw a histogram
options(shiny.maxRequestSize=30*1024^2)
server <- function(input, output) {
  
  #SAMPLES
  # loads in the metadata for use in the Samples tab
  sample_load_data <- reactive({
    get_file <- input$samplesFile$datapath
    if(is.null(get_file)){ return() }
    validate(
      need(file_ext(input$samplesFile$name) %in% c(
        'text/csv',
        'text/comma-separated-values',
        'text/tab-separated-values',
        'text/plain',
        'csv',
        'tsv'
      ), "Wrong File Format try again!"))
    read.table(get_file, header = TRUE, sep = ",")
  })
  # outputs a table of the metadata into the Samples > Metadata tab 
  output$Stable <- renderDT({
    #browser()
    file <- sample_load_data()
  })
  # outputs a summary table of the metadata into the Samples > Summary tab
  output$Ssummary <- renderDT({
    sum_df <-data.frame(
      Columns = names(sample_load_data()),
      Type = sapply(sample_load_data(), class),
      Mean_SD_or_Distinct_Values = sapply(sample_load_data(), function(col) {
        if (is.numeric(col)) {
          paste(round(mean(col), 2), "(+/", round(sd(col), 2), ")", sep="")
        } else {
          paste(length(unique(col)), " distinct values", sep ="")
        }
        
      })
    )
  })
  # outputs plots of the metadata into the Samples > Plots tab
  output$Shist <- renderPlot({
    file <- sample_load_data()
    ggplot(file, aes(x=Age.of.Death)) + 
      geom_histogram(color = "black", fill="orange", bins = 10)
  })
  output$Shist2 <- renderPlot({
    file <- sample_load_data()
    ggplot(file, aes(x=PMI)) + 
      geom_histogram(color = "black", fill="blue", bins = 10)
  })
  output$Shist3 <- renderPlot({
    file <- sample_load_data()
    ggplot(file, aes(x=RIN)) + 
      geom_histogram(color = "black", fill="purple", bins = 10)
  })
    
    
  
  #COUNTS
  counts_load_data <- eventReactive({input$countsFile}, {
    #browser()
    get_file <- input$countsFile$datapath
    if(is.null(get_file)){ return() }
    validate(
      need(file_ext(input$countsFile$name) %in% c(
        'text/csv',
        'text/comma-separated-values',
        'text/tab-separated-values',
        'text/plain',
        'csv',
        'tsv'
      ), "Wrong File Format try again!"))
    read.table(get_file, header = TRUE, sep = ",")
  })
  
  
  # filtering variance and non zero genes
  filter_zero_var_genes <- function(verse_counts, button_variance, button_nonzero, num_samples, num_genes) {
    # Compute variance for each gene
    variance <- apply(verse_counts[, -1], 1, var)
    
    # Filter genes based on variance and non-zero criteria
    filter_variance <- variance > button_variance
    nonzero <- rowSums(verse_counts[, -1] > 0) 
    filter_nonzero <- nonzero > button_nonzero
    filtered_genes <- which(filter_variance & filter_nonzero)
    
    # Apply the filters to get the final counts matrix
    filtered_counts <- verse_counts[filtered_genes, ]
    percentage_passing <- round((nrow(filtered_counts) / num_genes) * 100, 2)
    not_passing <- num_genes - nrow(filtered_counts)
    percentage_not_passing <- round((not_passing / num_genes) * 100, 2)
    
    # Create a summary dataframe
    c_sum_df <- data.frame(
      Measure = c("Number of samples", "Number of genes", "No. genes passing filter", 
                  "Percentage (%) passing filter", "No. genes not passing filter", "Percentage (%) not passing filter"),
      Value = c(num_samples, num_genes, nrow(filtered_counts), percentage_passing, not_passing, percentage_not_passing)
    )
    return(c_sum_df)
  }
  
  # plotting median vs variance
  plot_variance_vs_median <- function(data, button_variance, scale_y_axis=TRUE, title="") {
    #browser()
    matrix_data <- matrix(data[-1,-1])
    median_val <- rowMedians(matrix_data)
    calc_variance <- apply(matrix_data, 1, var)
    filter_variance <- calc_variance > button_variance
    
    df_data <- data.frame(Median_count = median_val, Variance = filter_variance)
    
    df_data$Rank <- rank(median_val)
    
    plot_val <- ggplot(df_data, aes(x = Rank, y = Variance)) +
      geom_point() +
      labs(x='Rank(Mean)', y='Variance') +
      geom_abline() +
      geom_smooth()
    ggtitle(title)
    if (scale_y_axis) {
      plot_val <- plot_val + scale_y_log10()
    }
    
    return(plot_val)
    
  }
  
  output$Csummary <- renderDT({
    input$countsSubmit
    c_file <<- counts_load_data()
    num_samples <-  ncol(c_file)
    num_genes <- nrow(c_file)
    isolate({filter_zero_var_genes(c_file, input$countsVariance, input$countsNonzero, num_samples, num_genes)})
  })
  
  
  output$CDiagplot1 <- renderPlot({
    input$countsSubmit
    isolate({plot_variance_vs_median(c_file, input$countsVariance)})
    
  })
  
  output$CDiagplot2 <- renderPlot({
    
    
    
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #DIFF EXP   
  # # loads in the diff expression data for use in the Diff. Expr. tab
  diff_load_data <- reactive({
    d_get_file <- input$DiffFile$datapath
    if(is.null(d_get_file)){ return() }
    validate(
      need(file_ext(input$DiffFile$name) %in% c(
        'text/csv',
        'text/comma-separated-values',
        'text/tab-separated-values',
        'text/plain',
        'csv',
        'tsv'
      ), "Wrong File Format try again!"))
    read.table(d_get_file, header = TRUE, sep = ",")
  })

  
  output$Dtable <- renderDT({
    d_file3 <<- diff_load_data()
  })
  
  
  #browser()
  # ggplot object of a volcano plot
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    ggplot_out <- ggplot(dataf, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), color = (padj < 1 * 10^slider))) +
      geom_point() +
      scale_color_manual(values = c(color1, color2)) +
      labs(title = "Volcano Plot",
           x = x_name,
           y = paste0("-log10(",y_name, ")"),
           color = "p-adjusted") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    return(ggplot_out)
  }
  
  #' Draw and filter table
  #' @return Data frame filtered to p-adjusted values that are less than 
  #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
  #' displayed.
  draw_table <- function(dataf, slider) {
    filtered_df <- dataf %>%
      filter(padj < 1 * 10^slider)
    formatted_df <- filtered_df %>%
      mutate(across(c("pvalue", "padj"), ~ formatC(., format = "e", digits = 20)))
    
    return(formatted_df)
  }
  
  # outputs volcano plot of the diff expr results into the Diff. Expr > Volcano plot tab
  output$Dvolcano <- renderPlot({
    input$DiffSubmit
    isolate({volcano_plot(d_file3, input$DiffButtonx, input$DiffButtony, input$DiffPadj_coloring, input$DiffBase_col, input$DiffHighlight_col)})
  })
  # outputs a filtered table of the results into the Diff. Expr > Plots table tab
  output$DPlottable <- renderDT({
    #browser()
    input$DiffSubmit
    isolate({draw_table(d_file3, input$DiffPadj_coloring)})
    
  })
  
  
  
  
  
  
  
    
}

# Run the application 
shinyApp(ui = ui, server = server)
