# Load packages ----
library(stringr)
# Load data ----
vocab = readLines("data/SLA/vocab.sla.txt")
title = readLines("data/SLA/title.sla.txt")
doi = readLines("data/SLA/doi.sla.txt")


# User interface ----
ui <- fluidPage(
  titlePanel("topicView"),

  sidebarLayout(
    sidebarPanel(
      helpText("Show topics computed from ebpmf-wbg model"),

      selectInput("n_topic",
                  label = "choose model (number of topics)",
                  choices = c(10, 50, 100),
                  selected = 10),

      selectInput("k",
                  label = "which topic to show",
                  choices = 1:100,
                  selected = 1)
      ),
    mainPanel(
      tabsetPanel(
        tabPanel("top words", htmlOutput("top_words")),
        tabPanel("top document", htmlOutput("top_doc"))
      )
    )

  )
)



server <- function(input, output) {
  ## load model
    output$top_words <- renderText({
      K = as.integer(input$n_topic)
      k = as.integer(input$k)
      model = readRDS(sprintf("output/SLA/v0.4.5/sla_ebpmf_wbg_initLF50_K%d_maxiter5000.Rds", K))
      L = model$qg$qls_mean
      F = model$qg$qfs_mean
      topic_k = sort(F[,k], decreasing = TRUE, index.return=TRUE)
      str_out = "word & weights(sum to 1)"
      for(i in 1:10){
        idx = topic_k$ix[i]
        word_weight = sprintf("%s:   %.3f", vocab[idx], F[idx,k]/sum(F[,k]))
        str_out = paste(str_out,word_weight, sep = "<br/>")
      }
      HTML(str_out)
    })

    output$top_doc <- renderText({
      K = as.integer(input$n_topic)
      k = as.integer(input$k)
      model = readRDS(sprintf("output/SLA/v0.4.5/sla_ebpmf_wbg_initLF50_K%d_maxiter5000.Rds", K))
      L = model$qg$qls_mean
      F = model$qg$qfs_mean
      doc_weight_k = sort(L[,k], decreasing = TRUE, index.return=TRUE)

      str_out = "paper title & DOI"
      for(i in 1:5){
        idx = doc_weight_k$ix[i]
        title_doi = sprintf("%d. %s\n   %s", i, title[idx], doi[idx])
        str_out = paste(str_out,title_doi, sep = "<br/>")
      }
      HTML(str_out)
    })
}

# Run app ----
shinyApp(ui, server)
