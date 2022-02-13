#'
#' @rdname shiny_ancom
#' 
#' @import shiny
#' @export
#' 

shiny_ancomUI <- fluidPage(
  headerPanel("Analysis of Composition of Microbiomes (ANCOM) v1.1-3"),
  sidebarLayout(
    sidebarPanel(     
      ##
      ## Data input
      ##
      fileInput('file1', 'Choose data file',
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv' )
      ),
      helpText("Accepted file formats: csv, txt"),
      textInput( 'file_out' , "Path of output file (without file name)" ),
      ##
      ## Main controls
      ##
      hr(),
      #helpText("Box edges flash red if selected value is out of bounds."),
      #hr(),
      selectInput(inputId="datafmt", width= '90%' ,
                  label  = "Format of input dataset:",
                  choices = c("Subjects on rows, OTUs on columns" = "wide",
                               "OTUs on rows, Subjects on columns" = "tall")
      ),
      #checkboxInput(inputId = "wexact",
      #              label = "Use exact p-values",
      #              value=FALSE
      #),
      checkboxInput(inputId = "repeated",
                    label = "Repeated measures",
                    value=FALSE
      ),
      checkboxInput(inputId = "adjust",
                    label = "Correct for multiple testing",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.adjust == true ",
                       numericInput(inputId = "fdr",
                                    label = "Enter desired FDR (between 0 and 1):",
                                    min=0 , max=1 , value=0.05, step=0.005
                       )
      ),
      conditionalPanel(condition = "input.adjust == false ",
                       numericInput(inputId = "alpha",
                                    label = "Enter significance level (between 0 and 1):",
                                    min=0 , max=1 , value=0.05, step=0.005
                       )
      ),
      #numericInput(inputId = "ncores",
      #             label = "Number of cores (to run in parallel)",
      #             min=1 , max=1000 , value=1, step=1
      #),
      checkboxInput(inputId = "fixPlot",
                    label = "Adjust figure:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.fixPlot == true ",
                       numericInput(inputId = "pltWidth",
                                    label = "Width of plot (%)",
                                    min=10 , max=100 , value=100, step=1
                       ),
                       numericInput(inputId = "pltHeight",
                                    label = "Height of plot (pixel)",
                                    min=100 , max=10000 , value=400, step=1
                       ),
                       sliderInput(inputId = "ncols",
                                   label = "Number of columns for the plot",
                                   min=1 , max=10 , value=3, step=1
                       )
      ),
      
      ##
      ## Action buttons
      ##
      hr(),
      actionButton(inputId = "compute1",
                   label = "Run ANCOM"
      ),
      actionButton(inputId = "compute2",
                   label = "Update Plot"
      ),
      helpText("'Update Plot' must be clicked before the plot will appear.")
    ),
    mainPanel(
      #uiOutput("plot.ui")
      #tabsetPanel(type = "tabs", 
      #            tabPanel("Summary", verbatimTextOutput("summary") ),
      #            tabPanel("Boxplots of Abundances", plotOutput("plot") ),
      #            tabPanel("Plot B", uiOutput("plot.ui") ),
      #            tabPanel("Etc"     , verbatimTextOutput("other") )
      #)
      
      uiOutput("theTabset")
      
    )
  ))












