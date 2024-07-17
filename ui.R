library(shiny)

fluidPage(
  titlePanel("Two Sample T-Test P-Value Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("n", "Sample Size (n)", min = 3, max = 20, value = 3),
      sliderInput("cv", "Coefficient of variation (%)", min = 0, max = 200, value = 50),
      sliderInput("log2_fold_change", "Log2 Fold Change", min = 0, max = 3, value = 1, step = 0.1),
      sliderInput("number", "Number", min = 1, max = 10000, value = 1),
      sliderInput("rank", "Rank", min = 1, max = 1000, value = 1)
    ),
    
    mainPanel(
      h3("Result"),
      verbatimTextOutput("p_value_text"),
      plotOutput("volcano_plot")
    )
  )
)
