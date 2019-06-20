#'Histogram of Data
#'
#'@description Shiny app to go through GEOquery datasets
#'
#'@usage Shiny()
#'
#'@author Nicholas Hutson
#'
#'@examples Shiny()
#'
#'@export
#'@import "DT"
#'@import "shiny"

library(DT)
library(shiny)
library(geneSummary)

# Define UI for app that allows a user to get data tables for gene expression, gene annotations, or sample annotations ----
ui <- fluidPage(

  titlePanel("Data"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      selectInput(inputId = "tableType", label = "Type of data table:", choices = c("Gene Expression" ,"Gene Annotations","Sample Annotations"), selected = "Gene Expression"),
      textInput(inputId = "DataID", label = "GEO accession ID:"),
      conditionalPanel(
        condition = "input.tableType != 'Sample Annotations'",
        selectizeInput(inputId = "geneSymbol", label = "Gene Symbol",choices = NULL)
      ),
      conditionalPanel(
        condition = "input.tableType == 'Sample Annotations'",
        selectizeInput(inputId = "sampleName", label = "Sample Name",choices = NULL)
      ),
      actionButton("submit", label = "Submit"),
      br(),
      br(),
      strong("Graphs for GSE43452:"),
      br(),
      actionButton("graph1", label = "TP53 Concentration"),
      br(),
      br(),
      actionButton("graph2", label = "TP53 in DBTRG and U87"),
      br(),
      br(),
      actionButton("graph3", label = "Top Variance Expressions"),
      br(),
      br(),
      actionButton("reset", label = "Reset")
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),

      conditionalPanel(
        condition = "output.table",
        conditionalPanel(
          condition = "input.tableType == 'Gene Expression'",
          titlePanel("Gene Expression Data")
        ),
        conditionalPanel(
          condition = "input.tableType == 'Gene Annotations'",
          titlePanel("Gene Annotations")
        ),
        conditionalPanel(
          condition = "input.tableType == 'Sample Annotations'",
          titlePanel("Sample Annotations")
        )#,
        #page selector, have it default to 1 so the first element of the list will be ouput by render graph at first

      ),

      uiOutput("length"),

      conditionalPanel(
        condition = "output.table",
        actionButton("save", label = "Save Data")
      ),

      #Hide errors
      #tags$style(type="text/css",
      #           ".shiny-output-error { visibility: hidden; }",
      #           ".shiny-output-error:before { visibility: hidden; }"),
      # Output: Data table ----
      DT::dataTableOutput("table"),

      plotOutput("graph")

    )
  )
)
# Define server logic required to draw a histogram ----
server <- function(input, output) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot

  data <- reactive({getGEO(input$DataID)})
  observe({
    input$DataID
    counter$countervalue <- 0
    setwd("/home/bios/Documents/geneSummary/Shiny") #find way to do this on any computer, not just this file path
  })

  observeEvent(input$submit, {
    ##create function to "render" data table before it is displayed to pass the length of the list back to the UI
    output$length <- renderUI({

      x <- data()

      y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

      z <- strsplit(input$sampleName,", ",fixed = TRUE)[[1]]

      #browser()

      if(identical(y, character(0))){
        if(identical(z, character(0))){
          tbl <- switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x))
        }else{
          tbl <- switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x,z))
        }
      }else if(identical(z, character(0))){
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x))
      }else{
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x,z))
      }

      #browser()

      updateSelectizeInput("geneSymbol",choices = tbl[[2]]$Symbol,server = TRUE) #may also have to create a separate function to update this based on the page selected
      updateSelectizeInput("sampleName",choices = tbl[[2]]$title,server = TRUE) #not realistic for over 50 options

      selectInput(inputId = "page", label = "Dataset:",choices = tbl[[1]])
    })

    output$table <- DT::renderDataTable({

      x <- data()

      y <- input$geneSymbol

      z <- input$sampleName

      #browser()

      if(identical(y, character(0))){
        if(identical(z, character(0))){
          tbl <- switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x))
        }else{
          tbl <- switch(input$tableType, "Gene Expression" = extExp(x), "Gene Annotations" = extGene(x),"Sample Annotations" = extSample(x,z))
        }
      }else if(identical(z, character(0))){
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x))
      }else{
        tbl <- switch(input$tableType, "Gene Expression" = extExp(x,y), "Gene Annotations" = extGene(x,y),"Sample Annotations" = extSample(x,z))
      }

      #browser()
      if(!is.data.frame(tbl[[2]])){
          n <- input$page #being delayed so is passing a NULL
          n <- match(n, tbl[[1]])
          #browser()
          if(is.null(n)){
            n <- 1
          }
          #browser()
          tbl <- tbl[[2]][[n]]
      }else{
        tbl <- tbl[[2]]
      }

      #browser()

      tbl

    }, rownames = TRUE)
  })

  counter <- reactiveValues(countervalue = 0, countervalueG = 0,countervalueE = 0, countervalueA = 0)

  observeEvent(input$save, {

    if(counter$countervalue == 0){
      counter$countervalue <- counter$countervalue + 1
      dir.create(paste("data.",input$page,sep = ""))
      setwd(paste("data.",input$page,sep = ""))
    }

    x <- data()

    y <- strsplit(input$geneSymbol,", ",fixed = TRUE)[[1]]

    z <- strsplit(input$sampleName,", ",fixed = TRUE)[[1]]

    if(identical(y, character(0))){
      if(identical(z, character(0))){
        y <- NA
        z <- NA
      }else{
        y <- NA
      }
    }else if(identical(z, character(0))){
      z <- NA
    }

    if(input$tableType == "Gene Expression"){
      tbl <- extExp(x,y)
      #browser()
      saveRDS(tbl, file = paste("geneExpression.",paste(y,collapse = "."),".rds", sep="")) #makes 2 names if geneSymbol is larger than 1
    }else if(input$tableType == "Gene Annotations"){
      tbl <- extGene(x,y)
      saveRDS(tbl, file = paste("geneAnnotation.",paste(y,collapse = "."),".rds", sep=""))
    }else{
      counter$countervalueA <- counter$countervalueA + 1
      #browser()
      tbl <- extSample(x,z)
      saveRDS(tbl, file = paste("sampleAnnotation",counter$countervalueE,".rds", sep=""))
    }
  })

  ##ONLY WORKS FOR GSE43452
  #use reactive values to define graph and just call the variable in render plot
  v <- reactiveValues(data = NULL)
  dataDef <- reactive({getGEO("GSE43452")})

  observeEvent(input$graph1, {
    x <- dataDef()

    D1a <- extGene(x)
    D2a <- extExp(x)
    D2b <- extSample(x)

    v$graph <- bar(D1a[[2]],D2a[[2]],D2b[[2]])
  })
  observeEvent(input$graph2, {
    x <- dataDef()

    D1a <- extGene(x)
    D2a <- extExp(x)
    D2b <- extSample(x)

    v$graph <- box(D1a[[2]],D2a[[2]])
  })
  observeEvent(input$graph3, {
    x <- dataDef()

    D1a <- extGene(x)
    D2a <- extExp(x)
    D2b <- extSample(x)

    v$graph <- hist(D1a[[2]],D2a[[2]])
  })

  observeEvent(input$reset, {
    v$graph <- NULL
  })

  output$graph <- renderPlot({

    v$graph

  })

}

shinyApp(ui = ui, server = server)
