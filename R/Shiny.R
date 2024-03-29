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
      selectInput(inputId = "importType", label = "Data type:", choices = c("GEO data", "RSE data"), selected = "GEO data"),
      conditionalPanel(
        condition = "input.importType == 'GEO data'",
        textInput(inputId = "DataID1", label = "GEO accession ID:")
      ),
      conditionalPanel(
        condition = "input.importType == 'RSE data'",
        textInput(inputId = "DataID2", label = "RSE accession ID:")
      ),
      actionButton("submit", label = "Submit"),
      br(),
      br(),
      conditionalPanel(
        condition = "output.n",
        actionButton("add_btn", "Add Filter")
      ),
      conditionalPanel(
        condition = "output.table",
        actionButton("rm_btn", "Remove Filter"),
        br(),
        br()
      ),
      uiOutput("textbox_ui1") #displays the filter list output
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      #prevents background warnings from being output for the user to see
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

      conditionalPanel(
        condition = "input.importType == 'GEO data'",
        uiOutput("length")),

      conditionalPanel(
        condition = "output.table",
        actionButton("save", label = "Save Data")
      ),
      br(),

      # Output: Data table ----
      DT::dataTableOutput("table")

    )
  )
)
# Define server logic
server <- function(input, output, session) {

  data <- reactive({
    if(input$importType == "GEO data"){
      #browser()
      getGEO(input$DataID1)
    }else{
      #browser()
      x <- grep("[A-Z]{3}[0-9]{6}", input$DataID2) #outputs "''" at first
      if(identical(x, integer(0))){
        rse_gene <- SummarizedExperiment()
      }else{
        rse_gene <- NULL
        download_study(input$DataID2) #throws an error that stops all code if the study ID doesn't exist
        y <- file.path(input$DataID2, 'rse_gene.Rdata')
        load(y)
      }
      rse_gene
    }
    })

  tab <- reactive({ #pre-rendered data and table so that the different types of tables can be accessed quickly without having to re-load the data from online

    if(input$importType == "GEO data"){
      #as summarized experiment, gene expressions = exprs(object)
      x <- data() #need to edit for longer expression sets
      #browser()
      tbl <- switch(input$tableType, "Gene Expression" = extExpGEO(x,long = input$long), "Gene Annotations" = extGeneGEO(x),"Sample Annotations" = extSampleGEO(x))
    }else{
      #browser()
      rse_gene <- data()
      #counts(rse_gene) = count of genes at each place
      #rowData = gene annotations
      #colData = sample annotations
      tbl <- switch(input$tableType, "Gene Expression" = extExpRSE(rse_gene,long = input$long), "Gene Annotations" = rowData(rse_gene),"Sample Annotations" = colData(rse_gene))
      tbl <- data.frame(tbl)
    }

    tbl
  })

  # Track all user inputs so that they are not reset when the number of filters is changed (they used to get re-rendered blank)
  AllInputs <- reactive({
    reactiveValuesToList(input)
  })

  #track the number of input boxes to render
  counter2 <- reactiveValues(n = 0)

  observeEvent(input$submit, {

    observeEvent(input$add_btn, {
      counter2$n <- counter2$n + 1
      Sys.sleep(1)
      #browser()
      })

    observeEvent(input$rm_btn, {
      if (counter2$n > 0) counter2$n <- counter2$n - 1
    })

    output$counter2 <- renderPrint(print(counter2$n))

    #used to remove add button when max filters have been reached
    output$n <- reactive({
      tbl <- tab()
      if(input$importType == "GEO data"){
        tbl <- tbl[[2]]
      }
      if(counter2$n >= length(colnames(tbl))){
        counter2$n <- length(colnames(tbl))
        FALSE
      }else{
        TRUE
      }
      #browser()
    })
    outputOptions(output, 'n', suspendWhenHidden = FALSE)

    #create UI inputs
    textboxes1 <- reactive({
      n <- counter2$n
      tbl <- tab()
      #if these 2 values are changed, the filters are re-rendered because they are reactive values

      if(input$importType == "GEO data"){
        tbl <- tbl[[2]]
      }

      if (n > 0) {
          isolate({lapply(seq_len(n), function(i) {
            #browser()
            list( #putting these in a list makes sure they come out in the correct order: select1, text1, select2, text2, etc.
              selectInput(inputId = paste0("selctin", i),
                          label = paste0(i, ". Filter by:"),
                          choices = colnames(tbl),
                          selected = AllInputs()[[paste0("selctin", i)]]),
              textInput(inputId = paste0("textin", i),
                          label = paste0("Value:"),
                          value = AllInputs()[[paste0("textin", i)]])
            )
          })
        })
      }
    })

    #output UI list
    output$textbox_ui1 <- renderUI({ textboxes1() })

    #create "page" input
    output$length <- renderUI({
        tbl <- tab()
        if(input$importType == "GEO data"){
          selectInput(inputId = "page", label = "Dataset:",choices = tbl[[1]])
        }else{
          NULL
        }
    })

    #create data table
    output$table <- DT::renderDataTable({

      tbl <- tab()
      num <- counter2$n
      if(input$importType == "GEO data"){
        if(num > 0){
          #searches for selct with a number on the end and gets all the inputs from inputId's like this
          #problem passing values from removed filters, might have to change rmv filter function to remove the value from the list
          y <- sapply(grep(pattern = "selctin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]])
          z <- as.character(sapply(grep(pattern = "textin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]]))
          z <- filterTbl(tbl[[2]], input$tableType, input$long, y, z)
        }else{
          z <- 0
        }

        #code to change "page" or dataset
        if(!is.data.frame(tbl[[2]])){
          n <- input$page
          n <- match(n, tbl[[1]]) #tbl[[1]] is a vector of data set names passed by the extraction functions
          #browser()
          if(is.null(n) | identical(n, integer(0))){
            n <- 1
          }
          #browser()
          tbl <- tbl[[2]][[n]]
        }else{
          tbl <- tbl[[2]]
        }

        #this narrows down the table if the filters exist
        if(num==0 | identical(z,integer(0)) | identical(z,0))  {
          tbl
        }else{
          tbl[z,]
        }
      }else{
        #browser()
        if(num > 0){
          #searches for selct with a number on the end and gets all the inputs from inputId's like this
          #problem passing values from removed filters, might have to change rmv filter function to remove the value from the list
          y <- sapply(grep(pattern = "selctin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]])
          z <- as.character(sapply(grep(pattern = "textin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]]))
          z <- filterTbl(tbl, input$tableType, input$long, y, z)
        }else{
          z <- 0
        }
        #this narrows down the table if the filters exist
        if(num==0 | identical(z,integer(0)) | identical(z,0))  {
          tbl
        }else{
          tbl[z,]
        }
      }
    }, rownames = FALSE)
  })

  observe({
    input$DataID1
    input$DataID2
    counter$countervalue <- 1 #sets to one whenever DataID is changed so that the user is asked which directory they want to save their data in later in the code, whenever this happens
  })

  counter <- reactiveValues(countervalue = 0,choice = getwd())

  #SAVE function
  observeEvent(input$save, {

    #browser()

    #set correct directory
    if(input$importType == "GEO data"){
      if(counter$countervalue == 1){
        if (exists('utils::choose.dir')) {
          counter$choice <- choose.dir(caption = "Select folder") #this is only available on Windows
        } else {
          counter$choice <- tk_choose.dir(caption = "Select folder")
        }
        setwd(counter$choice)
        #only call choose.dir if input$DataID is changed after save is hit
        d <- paste("data.",input$page,sep = "") #sets name of the new folder to data. the "page" or dataset selected
        dir.create(d)
        setwd(d)
      }else if(counter$countervalue > 1){ #if the data set hasn't been changed, the directory stays the same
        setwd(counter$choice)
        d <- paste("data.",input$page,sep = "")
        dir.create(d)
        setwd(d)
      }
    }else{
      if(counter$countervalue == 1){
        if (exists('utils::choose.dir')) {
          counter$choice <- choose.dir(caption = "Select folder") #this is only available on Windows
        } else {
          counter$choice <- tk_choose.dir(caption = "Select folder")
        }
        setwd(counter$choice)
        #only call choose.dir if input$DataID is changed after save is hit
        d <- paste("data.",input$DataID2,sep = "") #sets name of the new folder to data. the "page" or dataset selected
        dir.create(d)
        setwd(d)
      }else if(counter$countervalue > 1){ #if the data set hasn't been changed, the directory stays the same
        setwd(counter$choice)
        d <- paste("data.",input$DataID2,sep = "")
        dir.create(d)
        setwd(d)
      }
    }

    counter$countervalue <- counter$countervalue + 1

    #browser()

    #similar code to above, but has to keep the inputs from the text inputs, so that they can be used in the save names
    tbl <- tab()
    num <- counter2$n
    if(input$importType == "GEO data"){
      if(num > 0){
        #searches for selct with a number on the end and gets all the inputs from inputId's like this
        #problem passing values from removed filters, might have to change rmv filter function to remove the value from the list
        y <- sapply(grep(pattern = "selctin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]])
        z <- as.character(sapply(grep(pattern = "textin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]]))
        z <- filterTbl(tbl[[2]], input$tableType, input$long, y, z)
      }else{
        z <- 0
      }

      #code to change "page" or dataset
      if(!is.data.frame(tbl[[2]])){
        n <- input$page
        n <- match(n, tbl[[1]]) #tbl[[1]] is a vector of data set names passed by the extraction functions
        #browser()
        if(is.null(n)){
          n <- 1
        }
        #browser()
        tbl <- tbl[[2]][[n]]
      }else{
        tbl <- tbl[[2]]
      }

      #this narrows down the table if the filters exist
      if(num==0 | identical(z,integer(0)) | identical(z,0))  {
        tbl <- tbl
      }else{
        tbl <- tbl[z,]
      }
    }else{
      #browser()
      if(num > 0){
        #searches for selct with a number on the end and gets all the inputs from inputId's like this
        #problem passing values from removed filters, might have to change rmv filter function to remove the value from the list
        y <- sapply(grep(pattern = "selctin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]])
        z <- as.character(sapply(grep(pattern = "textin+[[:digit:]]", x = names(input), value = TRUE), function(x) input[[x]]))
        z <- filterTbl(tbl, input$tableType, input$long, y, z)
      }else{
        z <- 0
      }
      #this narrows down the table if the filters exist
      if(num==0 | identical(z,integer(0)) | identical(z,0))  {
        tbl <- tbl
      }else{
        tbl <- tbl[z,]
      }
    }

    #browser()

    if(identical(z,0) | identical(z,integer(0)))  {
      if(input$tableType == "Gene Expression"){
        #browser()
        if(input$long){
          saveRDS(tbl, file = paste("geneExpression.long.rds", sep="")) #adds a .long when the long data checkbox has been marked
        }else{
          saveRDS(tbl, file = paste("geneExpression.rds", sep=""))
        }
      }else if(input$tableType == "Gene Annotations"){
        saveRDS(tbl, file = paste("geneAnnotation.rds", sep=""))
      }else{
        #browser()
        saveRDS(tbl, file = paste("sampleAnnotation.rds", sep=""))
      }
      #the next save functions take the vector of user text inputs and add them to the end of the file name
    }else{
      if(input$tableType == "Gene Expression"){
        #browser()
        if(input$long){
          saveRDS(tbl, file = paste("geneExpression.long.",paste(z,collapse = "."),".rds", sep=""))
        }else{
          saveRDS(tbl, file = paste("geneExpression.",paste(z,collapse = "."),".rds", sep=""))
        }
      }else if(input$tableType == "Gene Annotations"){
        saveRDS(tbl, file = paste("geneAnnotation.",paste(z,collapse = "."),".rds", sep=""))
      }else{
        #browser()
        saveRDS(tbl, file = paste("sampleAnnotation.",paste(z,collapse = "."),".rds", sep=""))
      }
    }
  })
}
shinyApp(ui = ui, server = server)
