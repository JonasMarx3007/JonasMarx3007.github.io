library(ggplot2)
library(dplyr)

two_sample_t_test_p_value <- function(n, cv, log2_fold_change) {
  mean_y <- 1
  mean_x <- mean_y * 2^log2_fold_change
  
  sd_x <- (cv / 100) * mean_x
  sd_y <- (cv / 100) * mean_y
  var_x <- sd_x^2
  var_y <- sd_y^2
  
  pooled_var <- ((n - 1) * var_x + (n - 1) * var_y) / (2 * (n - 1))
  SE <- sqrt(pooled_var * (2 / n))
  t_value <- (mean_x - mean_y) / SE
  df <- 2 * (n - 1)
  p_value <- 2 * pt(-abs(t_value), df)
  
  return(p_value)
}

adjust_BH <- function(p_value, number, rank) {
  adj_pval <- p_value * number / rank
  adj_pval <- min(1, adj_pval)
  return(adj_pval)
}

create_volcano_plot <- function(number, cv, n) {
  set.seed(42)

  means1 <- rnorm(number, mean = 1, sd = (cv / 100))
  means2 <- rnorm(number, mean = 1, sd = (cv / 100))
  sds1 <- abs(rnorm(number, mean = 0.2, sd = (cv / 100)))
  sds2 <- abs(rnorm(number, mean = 0.2, sd = (cv / 100)))
  
  random_values1 <- abs(rnorm(number, mean = means1, sd = sds1))
  random_values2 <- abs(rnorm(number, mean = means2, sd = sds2))
  
  random_data <- data.frame(random_values1, random_values2)
  
  random_data <- random_data %>%
    mutate(log2foldchange = log2(random_values1) - log2(random_values2),
           pvalue = mapply(function(x, y) two_sample_t_test_p_value(n = n, cv = runif(1, cv-(1/2*cv), cv+(1/2*cv)), log2_fold_change = log2(x) - log2(y)),
                           random_values1, random_values2))
  
  random_data$adjpval = p.adjust(random_data$pvalue, method = "BH")
  
  random_data <- random_data %>%
    mutate(log10pval = -log10(adjpval),
           significance = case_when(
             adjpval < 0.05 & log2foldchange > 1 ~ "Upregulated",
             adjpval < 0.05 & log2foldchange < -1 ~ "Downregulated",
             TRUE ~ "Not Significant"
           ))
  
  ggplot(random_data, aes(x = log2foldchange, y = log10pval, color = significance)) + 
    geom_point(alpha = 0.5) + 
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
    theme_minimal() + 
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")
}

shinyServer(function(input, output) {
  
  output$p_value_text <- renderPrint({
    p_value <- two_sample_t_test_p_value(input$n, input$cv, input$log2_fold_change)
    
    adj_p_value <- adjust_BH(p_value, input$number, input$rank)
    
    cat("Original P-Value:", p_value, "\n")
    cat("Adjusted P-Value:", adj_p_value, "\n\n\n")
    cat("-log10 P-Value:", -log10(p_value), "\n")
    cat("-log10 adj. P-Value:", -log10(adj_p_value), "\n")
  })
  
  output$volcano_plot <- renderPlot({
    create_volcano_plot(input$number, input$cv, input$n)
  })
  
})

