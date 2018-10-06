library(tidyverse)
library(ggplot2)
library(shiny)
library(skewR)
library(truncnorm)
library(mvtnorm)
library(sn)

# Define UI for data upload app ----
ui <- fluidPage(
    
    # App title ----
    titlePanel("Bayesian Estimation of Skewness for Skew Normal Data"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Input: Select a file ----
            fileInput("file1", "Choose CSV File",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            
            # Horizontal line ----
            tags$hr(),
            
            # Input: Checkbox if file has header ----
            checkboxInput("header", "Header", TRUE),
            
            # Input: Select separator ----
            radioButtons("sep", "Separator",
                         choices = c(Comma = ",",
                                     Semicolon = ";",
                                     Tab = "\t"),
                         selected = ","),
            
            
            # Horizontal line ----
            tags$hr()
            
            
            
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Data file ----
            plotOutput("post"),
            textOutput("CI")
            
        )
        
    )
)

# Define server logic to read selected file ----
server <- function(input, output) {
    
    output$post <- renderPlot({
        
        req(input$file1)
        
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        tryCatch(
            {
                df <- read.csv(input$file1$datapath,
                               header = input$header,
                               sep = input$sep)
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        y <- df[,2]
        withProgress(message = "Running Gibbs sampler",
                     expr = {alphas <- skewR::sample_skew_posterior(y)})
        ggplot() + 
            geom_histogram(aes(x = alphas),alpha = 0.8,fill = "cadetblue") + 
            ggtitle("Posterior sample of skewness parameter") + 
            geom_vline(xintercept = quantile(alphas,0.025), color ="grey") + 
            geom_vline(xintercept = quantile(alphas,0.975), color = "grey") + 
            theme_minimal() + 
            labs(subtitle = paste("With 95% Credible Interval between",round(quantile(alphas,0.025),2),"and",round(quantile(alphas,0.975),2)))
        
    })
    
    
    
}

# Create Shiny app ----
shinyApp(ui, server)

