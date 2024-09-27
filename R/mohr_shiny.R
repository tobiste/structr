library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  headerPanel(
    "Mohr circle"
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
      sliderInput(
        inputId = "pf",
        label = "Pore fluid pressure (MPa)",
        min = 0,
        max = 1200,
        value = 0
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
  s1 <- reactive({
    input$S1
  })
  s2 <- reactive({
    input$S2
  })
  s3 <- reactive({
    input$S3
  })

  coulomb1 <- reactive({
    input$coulomb1
  })
  coulomb2 <- reactive({
    input$coulomb2
  })
  sliding <- reactive({
    input$sliding
  })

  PF <- reactive({
    input$pf
  })



  output$mohr <- renderPlot({
    s1 <- s1()
    s2 <- s2()
    s3 <- s3()
    pf <- PF()

    if (s2 > s1) {
      sx <- s1
      s1 <- s2
      s2 <- sx
    }
    if (s3 > s2) {
      sx <- s2
      s2 <- s3
      s3 <- sx
    }
    if (s1 < s3) {
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


    if (pf != 0) {
      s1_pf <- s1 - pf
      s2_pf <- s2 - pf
      s3_pf <- s3 - pf




      circle13.r <- diff_stress(s1, s3) / 2
      circle13.m <- mean_stress(s1, s3)

      circle12.r <- diff_stress(s1, s2) / 2
      circle12.m <- mean_stress(s1, s2)

      circle23.r <- diff_stress(s2, s3) / 2
      circle23.m <- mean_stress(s2, s3)


      pf_circle13.r <- diff_stress(s1_pf, s3_pf) / 2
      pf_circle13.m <- mean_stress(s1_pf, s3_pf)

      pf_circle12.r <- diff_stress(s1_pf, s2_pf) / 2
      pf_circle12.m <- mean_stress(s1_pf, s2_pf)

      pf_circle23.r <- diff_stress(s2_pf, s3_pf) / 2
      pf_circle23.m <- mean_stress(s2_pf, s3_pf)


      if (!is.null(coulomb)) {
        theta.f <- theta(coulomb[2]) # (90 + tectonicr:::atand(coulomb[2]))/2
      } else {
        theta.f <- 0
      }

      sigma_s <- shear_stress(s1_pf, s3_pf, theta.f / 2)
      sigma_n <- normal_stress(s1_pf, s3_pf, theta.f / 2)


      ggplot2::ggplot() +
        ggforce::geom_circle(aes(x0 = circle13.m, y0 = 0, r = circle13.r), fill = "grey", alpha = .3) +
        ggforce::geom_circle(aes(x0 = circle23.m, y0 = 0, r = circle23.r), fill = "white", alpha = .2) +
        ggforce::geom_circle(aes(x0 = circle12.m, y0 = 0, r = circle12.r), fill = "white", alpha = .2) +
        ggforce::geom_circle(aes(x0 = pf_circle13.m, y0 = 0, r = pf_circle13.r), fill = "slategrey", alpha = .5) +
        ggforce::geom_circle(aes(x0 = pf_circle23.m, y0 = 0, r = pf_circle23.r), fill = "white", alpha = .3) +
        ggforce::geom_circle(aes(x0 = pf_circle12.m, y0 = 0, r = pf_circle12.r), fill = "white", alpha = .3) +
        ggplot2::geom_abline(intercept = 0, slope = sliding, lty = 2) +
        ggplot2::geom_abline(intercept = coulomb[1], slope = coulomb[2], lty = 1) +
        ggplot2::geom_point(aes(x = circle13.m, 0), alpha = .2) +
        ggplot2::geom_point(aes(x = pf_circle13.m, 0)) +
        ggplot2::geom_line(aes(x = c(pf_circle13.m, sigma_n), y = c(0, sigma_s)), lty = 3) +
        ggplot2::geom_hline(yintercept = 0, alpha = .2) +
        ggplot2::geom_vline(xintercept = 0, alpha = .2) +
        ggplot2::geom_text(aes(x = (s1_pf+s3_pf)/2, y = 0), label = expression(sigma["m"]), vjust = -.5, hjust = -1) +
        ggplot2::geom_text(aes(x = s3_pf, y = 0), label = expression(sigma[3]), vjust = -.5, hjust = -1) +
        ggplot2::geom_text(aes(x = s2_pf, y = 0), label = expression(sigma[2]), vjust = -.5, hjust = -1) +
        ggplot2::geom_text(aes(x = s1_pf, y = 0), label = expression(sigma[1]), vjust = -.5, hjust = -1) +
        ggplot2::coord_fixed() +
        ggplot2::labs(x = bquote(sigma[n] ~ (.(units))), y = bquote(sigma[s] ~ (.(units))), caption = bquote(theta[f] == .(round(theta.f, 2)) ~ alpha[f] == .(round(90 - theta.f, 2)))) +
        ggplot2::theme_classic()
    } else {
      ggMohr(s1, s2, s3, coulomb, sliding)
    }
  })
}

shinyApp(ui, server)
