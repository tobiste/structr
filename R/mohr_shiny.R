library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  headerPanel(
    'Mohr circle'
  ),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      h4("Stress magnitudes"),
       sliderInput(
        inputId = "S1",
        label = "Sigma 1 (MPa)",
        min = -100,
        max = 1200,
        value = 1025
      ),
      sliderInput(
        inputId = "S2",
        label = "Sigma 2 (MPa)",
        min = -100,
        max = 1200,
        value = 400
      ),
      sliderInput(
        inputId = "S3",
        label = "Sigma 3 (MPa)",
        min = -100,
        max = 1200,
        value = 250
      ),
      
      fluidRow(
        h4("Coulomb criteria"),
        sliderInput(
        inputId = "coulomb1",
        label = "Cohesion",
        min = -500,
        max = 500,
        value = 70,
        round = FALSE, 
        step = 1
      ),
      sliderInput(
        inputId = "coulomb2",
        label = "Slope",
        min = 0,
        max = 2,
        value = .6,
        round = FALSE,
        step = 0.01
      )
      ),
      
      sliderInput(
        inputId = "sliding",
        label = "Sliding criteria",
        min = 0,
        max = 2,
        value = .81,
        round = FALSE,
        step = 0.01
      )
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Interactive map ----
      plotOutput(outputId = "mohr", height = "600px")
    )
  )
)

server <- function(input, output) {
  s1 <- reactive({input$S1})
  s2 <- reactive({input$S2})
  s3 <- reactive({input$S3})
  
  coulomb1 <- reactive({input$coulomb1})
  coulomb2 <- reactive({input$coulomb2})
  sliding <- reactive({input$sliding})
  
  
  
  output$mohr <- renderPlot({
    s1 <- s1()
    s2 <- s2()
    s3 <- s3()
    
    if(s2 > s1){
      sx <- s1
      s1 <- s2
      s2 <- sx
    }
    if(s3 > s2){
      sx <- s2
      s2 <- s3
      s3 <- sx
    }
    if(s1 < s3){
      sx1 <- s1
      sx2 <- s2
      sx3 <- s3
      s3 <- sx1
      s2 <- sx3
      s1 <- sx2
    }
    # if(s3 < s2){
    #   sx <- s2
    #   s2 <- s3
    #   s3 <- sx
    # }
    
    coulomb <- c(coulomb1(), coulomb2())
    sliding <- sliding()
    ggMohr(s1, s2, s3, coulomb, sliding)
  })
   
}

shinyApp(ui, server)

