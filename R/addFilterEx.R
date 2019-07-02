#to use this, I have to change it so that a filter by dropdown renders and then a selectize input comes under it. extraction functions will have to take a
#vector or list of parameters, increment through it, and make sure all parameters match

library(shiny)

ui <- shinyUI(fluidPage(

  sidebarPanel(

    actionButton("add_btn", "Add Textbox"),
    actionButton("rm_btn", "Remove Textbox"),
    textOutput("counter")

  ),

  mainPanel(uiOutput("textbox_ui"))

))

server <- shinyServer(function(input, output, session) {

  # Track the number of input boxes to render
  counter <- reactiveValues(n = 0)

  # Track all user inputs
  AllInputs <- reactive({
    x <- reactiveValuesToList(input)
  })

  observeEvent(input$add_btn, {counter$n <- counter$n + 1})
  observeEvent(input$rm_btn, {
    if (counter$n > 0) counter$n <- counter$n - 1
  })

  output$counter <- renderPrint(print(counter$n))

  textboxes <- reactive({

    n <- counter$n

    if (n > 0) {
      isolate({
        lapply(seq_len(n), function(i) {
          textInput(inputId = paste0("textin", i),
                    label = paste0("Textbox", i),
                    value = AllInputs()[[paste0("textin", i)]])
        })
      })
    }

  })

  output$textbox_ui <- renderUI({ textboxes() })

})

shinyApp(ui, server)
