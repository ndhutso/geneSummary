# Define UI for app that allows a user to get data tables for gene expression, gene annotations, or sample annotations ----
#GEOquery, stringr
ui <- fluidPage(

  titlePanel("Data"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      selectInput(inputId = "tableType", label = "Type of data table:", choices = c("Gene Expression" ,"Gene Annotations","Sample Annotations"), selected = "Gene Expression"),
      conditionalPanel(
        condition = "input.tableType == 'Gene Expression'",
        checkboxInput("long", "Long format data", value = FALSE)
      ),
      textInput(inputId = "DataID", label = "GEO accession ID:"),
      actionButton("submit", label = "Submit"),
      br(),
      br(),
      actionButton("add_btn", "Add Filter"),
      actionButton("rm_btn", "Remove Filter"),
      br(),
      br(),
      textOutput("counter2"),
      br(),
      uiOutput("textbox_ui1")
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
        )
      ),

      uiOutput("length"),

      conditionalPanel(
        condition = "output.table",
        actionButton("save", label = "Save Data")
      ),
      br(),

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
# Define server logic
server <- function(input, output) {

  #START: mess with adding filters

  # Track the number of input boxes to render
  counter2 <- reactiveValues(n = 0)

  # Track all user inputs
  AllInputs <- reactive({
    x <- reactiveValuesToList(input)
  })

  observeEvent(input$add_btn, {counter2$n <- counter2$n + 1})
  observeEvent(input$rm_btn, {
    if (counter2$n > 0) counter2$n <- counter2$n - 1
  })

  output$counter2 <- renderPrint(print(counter2$n))

  textboxes1 <- reactive({

    n <- counter2$n

    if (n > 0) {
      lapply(seq_len(n), function(i) {
        list(
          selectInput(inputId = paste0("textin", i), #have to pull choices from below
                    label = paste0("Textbox", i),
                    choices = NULL),
          selectizeInput(inputId = paste0("selctin", i),
                         label = paste0("Select", i),
                         choices = NULL)
        )
      })
    }

  })

  output$textbox_ui1 <- renderUI({ textboxes1() })

  #END

  data <- reactive({getGEO(input$DataID)})

  tab <- reactive({
    x <- data()

    y <- character(0)

    z <- character(0)

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
    tbl
  })

  observeEvent(input$submit, {
    ##create function to "render" data table before it is displayed to pass the length of the list back to the UI
    output$length <- renderUI({

        tbl <- tab()

        #browser()
        selectInput(inputId = "page", label = "Dataset:",choices = tbl[[1]])
    })

    output$table <- DT::renderDataTable({

      tbl <- tab()

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

      tbl

    }, rownames = FALSE)
  })

  observe({
    input$DataID
    counter$countervalue <- 1 #sets to one whenever DataID is changed
    #setwd(counter$choice) #have to pass choice to here
  })

  counter <- reactiveValues(countervalue = 0,choice = getwd())

  observeEvent(input$save, {

    #browser()

    if(counter$countervalue == 1){
      if (exists('utils::choose.dir')) {
        counter$choice <- choose.dir(caption = "Select folder")
      } else {
        counter$choice <- tk_choose.dir(caption = "Select folder")
      }
      setwd(counter$choice)
        #only call choose.dir if input$DataID is changed after save is hit

      d <- paste("data.",input$page,sep = "")
      dir.create(d)
      setwd(d)
    }else if(counter$countervalue > 1){
      setwd(counter$choice)
      d <- paste("data.",input$page,sep = "")
      dir.create(d)
      setwd(d)
    }

    counter$countervalue <- counter$countervalue + 1

    #browser()

    x <- data()

    y <- character(0)

    z <- character(0)

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

    if(is.na(y)){
      if(is.na(z)){
        if(input$tableType == "Gene Expression"){
          tbl <- extExp(x,y,long = input$long)
          #browser()
          if(input$long){
            saveRDS(tbl, file = paste("geneExpression.long.rds", sep=""))
          }else{
            saveRDS(tbl, file = paste("geneExpression.rds", sep=""))
          }
        }else if(input$tableType == "Gene Annotations"){
          tbl <- extGene(x,y)
          saveRDS(tbl, file = paste("geneAnnotation.rds", sep=""))
        }else{
          #browser()
          tbl <- extSample(x,z)
          saveRDS(tbl, file = paste("sampleAnnotation.rds", sep=""))
        }
      }else{
        if(input$tableType == "Gene Expression"){
          tbl <- extExp(x,y,long = input$long)
          #browser()
          if(input$long){
            saveRDS(tbl, file = paste("geneExpression.long.rds", sep=""))
          }else{
            saveRDS(tbl, file = paste("geneExpression.rds", sep=""))
          }
        }else if(input$tableType == "Gene Annotations"){
          tbl <- extGene(x,y)
          saveRDS(tbl, file = paste("geneAnnotation.rds", sep=""))
        }else{
          #browser()
          tbl <- extSample(x,z)
          saveRDS(tbl, file = paste("sampleAnnotation.",paste(z,collapse = "."),".rds", sep=""))
        }
      }
    }else if(is.na(z)){
      if(input$tableType == "Gene Expression"){
        tbl <- extExp(x,y,long = input$long)
        #browser()
        if(input$long){
          saveRDS(tbl, file = paste("geneExpression.",paste(y,collapse = "."),".long.rds", sep=""))
        }else{
          saveRDS(tbl, file = paste("geneExpression.",paste(y,collapse = "."),".rds", sep=""))
        }
      }else if(input$tableType == "Gene Annotations"){
        tbl <- extGene(x,y)
        saveRDS(tbl, file = paste("geneAnnotation.",paste(y,collapse = "."),".rds", sep=""))
      }else{
        #browser()
        tbl <- extSample(x,z)
        saveRDS(tbl, file = paste("sampleAnnotation.rds", sep=""))
      }
    }else{
      if(input$tableType == "Gene Expression"){
        tbl <- extExp(x,y,long = input$long)
        #browser()
        if(input$long){
          saveRDS(tbl, file = paste("geneExpression.",paste(y,collapse = "."),".long.rds", sep=""))
        }else{
          saveRDS(tbl, file = paste("geneExpression.",paste(y,collapse = "."),".rds", sep=""))
        }
      }else if(input$tableType == "Gene Annotations"){
        tbl <- extGene(x,y)
        saveRDS(tbl, file = paste("geneAnnotation.",paste(y,collapse = "."),".rds", sep=""))
      }else{
        #browser()
        tbl <- extSample(x,z)
        saveRDS(tbl, file = paste("sampleAnnotation.",paste(z,collapse = "."),".rds", sep=""))
      }
    }
  })
}

shinyApp(ui = ui, server = server)
