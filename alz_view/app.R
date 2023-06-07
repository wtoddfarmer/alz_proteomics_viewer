library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(org.Hs.eg.db)
library(GO.db)

johnson_summary <- readRDS("summary.RDS")
johnson_expression <- readRDS("expression.RDS")
johnson_sample_meta <- readRDS("sample_meta.RDS")

getGeneSummary <- function(gene_symbol){
  
  gene_summary <- johnson_summary %>%
    filter(Symbol == gene_symbol) %>%
    drop_na(ANOVA_PValue) %>% arrange(ANOVA_PValue)
  
  return(gene_summary)
}

getPlotData <- function(gene){
  
  gene_summary <- getGeneSummary(gene)
  
  gene_data <- johnson_expression %>% 
    filter(Symbol == gene, UniprotAccession %in% gene_summary$UniprotAccession) %>% dplyr::select(-Symbol)
  
  plot_data <- gene_data %>%
    pivot_longer(cols = starts_with(c("banner", "rosmap")), names_to = "case", values_to = c("l")) %>%
    mutate(l=as.numeric(l)) %>%
    mutate(UniprotAccession = factor(UniprotAccession, levels = gene_summary$UniprotAccession)) %>%
    left_join(johnson_sample_meta, by = "case") 
  
  return(plot_data)
}

multiPlot <- function(gene){
  
  gene_summary <- getGeneSummary(gene)
  
  
  p <- ggplot(getPlotData(gene)) +
    geom_violin(aes(y=l, x=diagnosis, color=diagnosis)) +
    geom_jitter(aes(y=l, x=diagnosis, color=diagnosis), alpha=0.5, width = 0.2) + 
    facet_wrap(~UniprotAccession,  ncol=2, scales = "free") +
    ylab("Normalized Expression")
  return(p)
}

getMultiGeneSummary <- function(gene_symbol_list){
  
  multi_gene_summary <- johnson_summary %>%
    filter(Symbol %in% gene_symbol_list) %>%
    drop_na(ANOVA_PValue) %>% arrange(ANOVA_PValue)
  
  return(multi_gene_summary)
}

getGoGenes <- function(term){
  go_res <- AnnotationDbi::select(org.Hs.eg.db, term, keytype = "GO", columns = "SYMBOL")
  genes <- unique(go_res$SYMBOL)
  return(genes)
}

getGoDes <- function(term){
  go_des <- AnnotationDbi::select(GO.db, term, keytype = "GOID", columns = "TERM")
  return(go_des$TERM)
  
}

ui <- fluidPage(
  
  
  titlePanel("Human Proteomics Data Viewer"),
  tags$h4("This app shows the expression of the gene products in the contex of Alzheimer's disease. The data used is from the publication", 
          tags$a(href="https://www.nature.com/articles/s41593-021-00999-y", "Johnson et. al, 2022."), "The data consists of the quantification of over 13k peptides from the dorsolateral prefrontal cortex of 516 human subjects."),
  tabsetPanel(
    tabPanel("Single Gene Search",
             fluidPage(
               sidebarLayout(
                 sidebarPanel(
                   tags$h4("Check to see if a single gene is changed in AD"),
                   textInput("GOI", "Input a Human gene symbol without quotes (i.e. PARK7) :", value = "MTOR", width = NULL,
                             placeholder = NULL),
                   tags$h4("If you get an error then the symbol you entered is either incorrect or the products of the gene were not captured in the study."),
                   width = 3,
                   fluid = FALSE
                 ),
                 
                 # Show a plot of the generated distribution
                 mainPanel(
                   plotOutput("violin"),
                   tableOutput("table"),
                   width = 9,
                   fluid = TRUE
                 )
                 
               )
             )
    ), 
    tabPanel("GO Term Query",
             fluidPage(
               sidebarLayout(
                 sidebarPanel(
                   tags$h4("Check to see if members of a GO term are changed in AD. Useful to see which of several related proteins change."),
                   textInput("GO", "Input a GO Term :", value = "GO:0005245", width = NULL,
                             placeholder = NULL),
                   tags$h4("A good place to explore GO terms and their relatives is:",  tags$a(href="http://www.informatics.jax.org/vocab/gene_ontology/GO:0005245", "MGI Gene Ontology Browser")),
                   width = 3
                 ),
                 mainPanel(
                   tags$h1(textOutput("textGO")),
                   tableOutput("table2"),
                   width = 9
                   
                 )
               )
             )
    )
  )
)



server <- function(input, output) {
  
  output$table <- renderTable(getGeneSummary(input$GOI), spacing = "xs",digits = -2)
  output$violin <- renderPlot(multiPlot(input$GOI))
  output$textGO <- renderText(getGoDes(input$GO))
  output$table2 <- renderTable(getMultiGeneSummary(getGoGenes(input$GO)), spacing = "xs",digits = -2, striped = TRUE)
}

# Run the application 
shinyApp(ui = ui, server = server)
