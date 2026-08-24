library(shiny)
library(httr)
library(jsonlite)
library(visNetwork)
library(dplyr)
library(igraph)  
library(tidygraph)
library(shinycssloaders)
library(shinyjs)
library(htmltools)

# Função para chamar a API, para teste local descomente a linha com "localhost"
api_base_url <- "http://tissueppi-api:8003"
#api_base_url <- "http://localhost:8003"

addResourcePath("web_data", "/app/data")

sqlite_path <- file.path(getwd(), "data", "interactions.sqlite")

# Função para buscar dados da API
fetch_from_api <- function(endpoint, params = NULL) {
  url <- paste0(api_base_url, endpoint)
  response <- tryCatch({
    if (is.null(params)) {
      GET(url)
    } else {
      GET(url, query = params)
    }
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(response) || status_code(response) != 200) {
    return(NULL)
  }
  
  content(response, "parsed", simplifyDataFrame = TRUE)
}

# Carregar os códigos dos tecidos
tissue_input <- fetch_from_api("/tissue_list")

# Carregar os códigos das proteínas
protein_input <- fetch_from_api("/protein_list")

ui_page <- reactiveVal("home")

#ui <- uiOutput("main_ui")
ui <- tagList(
  useShinyjs(), # Inicializa o shinyjs
  tags$head(
    tags$style(HTML("
  .btn-download-modern {
    background-color: #d97706;
    color: white;
    font-weight: 600;
    font-size: 1rem;
    padding: 12px 24px;
    border: none;
    border-radius: 8px;
    cursor: pointer;
    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
    transition: all 0.3s ease;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    text-decoration: none;
  }

  /* Efeito ao passar o mouse (Hover) */
  .btn-download-modern:hover {
    background-color: #b45309;
    box-shadow: 0 6px 12px -2px rgba(0, 0, 0, 0.15);
    transform: translateY(-1px);
  }
  /* Efeito ao clicar (Active) */
  .btn-download-modern:active {
    transform: translateY(1px);
    box-shadow: 0 2px 4px -1px rgba(0, 0, 0, 0.1);
  }
      .spinner { border: 3px solid rgba(255,255,255,0.4); border-top: 3px solid #ffffff; 
                 border-radius: 50%; width: 20px; height: 20px; 
                 animation: spin 0.8s linear infinite;
                 display: inline-block; 
                 vertical-align: middle; margin-right: 10px; }
      @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }

    "))
  ),
  uiOutput("main_ui")
)
server <- function(input, output, session) {

  # Criar o grafo para exibição no visNetwork
  createGraph <- function(arestas, metricas) {
     grafo <- graph_from_data_frame(arestas, directed = FALSE)
     metricas$color <- color_map(metricas$clustering*100)
     metricas <- metricas %>% mutate(size = (scale_size(metricas$degree)))
     nodesGraph <- data.frame(proteinName = V(grafo)$name, stringsAsFactors = FALSE)
     nodesGraph <- nodesGraph %>% left_join(metricas, by = c("proteinName" = "geneSymbol"))
  
  # Converter o igraph para um objeto tbl_graph
     g_tidy <- as_tbl_graph(grafo)
  
  # Fazer o join direto como se fosse uma tabela
     g_tidy <- g_tidy %>%
       activate(nodes) %>%
       left_join(nodesGraph, by = c("name" = "proteinName"))
     grafo <- as.igraph(g_tidy)
     coords <- layout_with_fr(grafo) * 100
     V(grafo)$x <- coords[, 1]
     V(grafo)$y <- coords[, 2]
     return(grafo)
  }
  # Calcular a escala para o tamaho do nó.
  scale_size <- function(x) {
    if (all(x == 0)) return(rep(5, length(x)))
    scaled <- log10(x) * 25
    scaled
  }
  # Definir a cor baseada no coeficiente de clusterização.
  color_map <- function(values, initialColor = "#1f77b4", finalColor = "#d62728") {
    
    # 1. Cria a função de rampa de cores
    graphPalette <- colorRampPalette(c(initialColor, finalColor))
    
    # 2. Garante que os valores estejam dentro do intervalo [0, 100]
    # Isso evita erros se algum valor estiver fora do esperado
    limitValues <- pmax(pmin(values, 100), 0)
    
    # 3. Mapeia os valores para a paleta de 100 níveis
    # +1 pois a paleta é indexada de 1 a 100
    indices <- round(limitValues) + 1
    
    # Retorna o vetor de cores correspondente
    return(graphPalette(101)[indices])
  }
  output$protein_selector_ui <- renderUI({
    selectInput("protein", "", 
                   choices = protein_input$geneSymbol,
                    selected = "FUCA2")
  })
  
  output$tissue_selector_ui <- renderUI({
    selectInput("tissue", "",
                   choices = setNames(tissue_input$idTissue, tissue_input$descTissue),
                   selected = "BLDDER")
  })
  
  output$multiple_protein_selector_ui <- renderUI({
    selectizeInput("protein_multi", "", 
                   choices = protein_input$geneSymbol, multiple = TRUE)
  })
  
  # Dados das proteínas  
  output$protein_degree_ui <- renderText({
    datalist_rv()$proteinMetrics$degree
  })
  
  output$protein_betweenness_ui <- renderText({
    datalist_rv()$proteinMetrics$betweenness
  })
  
  output$protein_closeness_ui <- renderText({
    datalist_rv()$proteinMetrics$closeness
  })
  
  output$protein_count_ui <- renderText({
    datalist_rv()$proteinMetrics$degree
  })
  
  output$protein_name_ui <- renderText({
    datalist_rv()$proteinMetrics$geneSymbol
  })
  
  output$protein_geneid_ui <- renderText({
    datalist_rv()$proteinMetrics$idGene
  })
  
  output$protein_interactions_ui <- renderText({
    nrow(datalist_rv()$all_interactions)
  })

  output$protein_clustering_ui <- renderText({
    datalist_rv()$proteinMetrics$clustering
  })

  # Dados dos tecidos
  output$tissue_name_ui <- renderText({
    datalist_rv()$tissueMetrics$descTissue
  })
  
  output$protein_diameter_ui <- renderText({
    datalist_rv()$tissueMetrics$diameter
  })
  
  output$tissue_interactions_ui <- renderText({
    datalist_rv()$tissueMetrics$totalEdges
  })
  
  output$tissue_distance_ui <- renderText({
    datalist_rv()$tissueMetrics$distancia
  })
  
  output$tissue_proteins_ui <- renderText({
    datalist_rv()$tissueMetrics$totalNodes
  })
  
  output$tissue_diameter_ui <- renderText({
    datalist_rv()$tissueMetrics$diameter
  })
  
  output$tissue_density_ui <- renderText({
    datalist_rv()$tissueMetrics$density
  })

  output$tissue_clustering_ui <- renderText({
    datalist_rv()$tissueMetrics$clustering
  })

    output$tissue_degree_ui <- renderText({
    datalist_rv()$tissueMetrics$avgDegree
  })
  
  datalist_rv <- reactiveVal(NULL)
  graph_edges_rv <- reactiveVal(NULL)  

  output$download_edges <- downloadHandler(

    filename = function() {
      paste0("protein_", datalist_rv()$proteinMetrics$geneSymbol, "_", datalist_rv()$tissueMetrics$idTissue, "_edgelist_", Sys.Date(), ".csv")
    },
    content = function(file) {
      edges <- graph_edges_rv()
      if (!is.data.frame(edges) || nrow(edges) == 0) {
        div(
            style = "padding: 8px 14px;
                   color: #b71c1c;
                   background-color: #fdecea;
                   border: 1px solid #f5c6cb;
                   border-radius: 6px;
                   font-size: 16px;
                   font-family: Arial, sans-serif;",
            "No interactions found !"
        )
      } else {
        write.csv(edges, file, row.names = FALSE)
      }
    }
  )

  output$download_database <- downloadHandler(
    
    filename = function() {
      paste0("interactions.sqlite")
    },
    content = function(file) {
      runjs("$('#download_database').prop('disabled', true).addClass('loading');")
      runjs("$('#download_database').prepend('<div class=\"spinner\"></div>');")
      file.copy("/app/data/interactions.sqlite", file)
      runjs("
            setTimeout(function() {
              $('.spinner').remove();
              $('#download_database').prop('disabled', false).removeClass('loading');
            }, 10000);")
    }        
  )
  
  output$download_tissue_edges <- downloadHandler(

    filename = function() {
      paste0("tissue_",datalist_rv()$tissueMetrics$idTissue, "_edgelist_", Sys.Date(), ".csv")
    },
    content = function(file) {
      edges <- datalist_rv()$all_tissue_edges
      if (!is.data.frame(edges) || nrow(edges) == 0) {
        div(
          style = "padding: 8px 14px;
                   color: #b71c1c;
                   background-color: #fdecea;
                   border: 1px solid #f5c6cb;
                   border-radius: 6px;
                   font-size: 16px;
                   font-family: Arial, sans-serif;",
          "No interactions found !"
        )
      } else {
        write.csv(edges, file, row.names = FALSE)
      }
    }
  )

  output$download_network <- downloadHandler(

    filename = function() {
      paste0("network_",datalist_rv()$proteinMetrics$geneSymbol, "_", datalist_rv()$tissueMetrics$idTissue, "_", Sys.Date(), ".html")
    },
      content = function(file) {
        edges <- graph_edges_rv()
        if (!is.data.frame(edges) || nrow(edges) == 0) {
        div(
            style = "padding: 8px 14px;
                   color: #b71c1c;
                   background-color: #fdecea;
                   border: 1px solid #f5c6cb;
                   border-radius: 6px;
                   font-size: 16px;
                   font-family: Arial, sans-serif;",
            "No interactions found !"
            
        )
      } else {
           graphTitle <- paste0(datalist_rv()$proteinMetrics$geneSymbol, " protein network ")
           graphSubTitle <- paste0("Tissue ", datalist_rv()$tissueMetrics$idTissue)
           nodesMetrics <- datalist_rv()$allProteinsMetrics
           grafo <- createGraph(edges, nodesMetrics)
           grafoVis <- toVisNetworkData(grafo)
           vis <- visNetwork(nodes = grafoVis$nodes, edges = grafoVis$edges,
                      main = graphTitle,
                      submain = graphSubTitle) %>% 
           visNetwork::visOptions(highlightNearest = FALSE, nodesIdSelection = FALSE) %>%
           visIgraphLayout() %>% 
           visInteraction(
             dragNodes = FALSE,       # Impede arrastar os nós
             dragView = FALSE,        # Impede arrastar a tela de fundo
             zoomView = FALSE,        # Desativa o zoom (roda do mouse)
             selectable = FALSE,      # Impede selecionar nós/arestas com o clique
             hover = FALSE            # Desativa efeitos ao passar o mouse por cima
           ) %>% 
             visPhysics(enabled = FALSE) # Desativa a física
           visNetwork::visSave(vis, file = file)
      }

    }
  )
  
  # Renderizar o HTML
  output$main_ui <- renderUI({
    nav_buttons <- list(
      home = actionButton("home", "TissuePPI"),
      go_home = actionButton("go_home", "Home"),
      go_about = actionButton("go_about", "About"),
      go_search = actionButton("single_search", "Search"),
      go_help = actionButton("go_help", "Help"),
      go_downloadDB = actionButton("go_downloadDB", "Download DB"),
      multiple_search = actionButton("multiple_search", "Multiple Proteins"),
      single_search = actionButton("single_search", "Single Protein"),
      search_single_protein = actionButton("search_single_protein", "Search"),
      search_multiple_protein = actionButton("search_multiple_protein", "Search")
    )
    

    output$graph_ui_single <- renderUI({
      req(input$search_single_protein)
      edges <- datalist_rv()$all_interactions
      
      if (!is.data.frame(edges) || is.null(nrow(edges))) {
        div(
          style = "padding: 8px 14px;
                   color: #b71c1c;
                   background-color: #fdecea;
                   border: 1px solid #f5c6cb;
                   border-radius: 6px;
                   font-size: 16px;
                   font-family: Arial, sans-serif;",
          "No interactions found !"

        )
      } else {
        nodesMetrics <- datalist_rv()$allProteinsMetrics
        grafo <- createGraph(edges, nodesMetrics)
       V(grafo)$title <- paste0(
          "Proteín: ", V(grafo)$name, "<br>",
          "Degree: ", V(grafo)$degree, "<br>",
          "Betweenness: ", V(grafo)$betweenness, "<br>",
          "Coef. Clustering: ", V(grafo)$clustering
        )
          visIgraph(grafo) %>%
          visEdges(arrows = "none") %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

              }
    })
    # Cria a legenda para o grafo, o visNetwork não cria legendas em gradiente 
    output$graph_ui_legend <- renderUI({
      nodesMetrics <- datalist_rv()$allProteinsMetrics
      minDegree <- min(nodesMetrics$degree,na.rm = TRUE)
      avgDegree <- mean(nodesMetrics$degree,na.rm = TRUE)
      maxDegree <- max(nodesMetrics$degree,na.rm = TRUE)
      minDegree <- round(minDegree,1)
      avgDegree <- round(avgDegree,1)
      maxDegree <- round(maxDegree,1)

      tags$div(
      style = "display: flex; justify-content: center; gap: 40px; background: #ffffff; padding: 15px; border-top: 1px solid #ddd; font-family: Arial, sans-serif; font-size: 12px; color: #333;",
      
      # Escala Gradiente de Cor
      tags$div(
        style = "display: flex; flex-direction: column; align-items: center;",
        tags$div("Node color: Clustering", style = "font-weight: bold; margin-bottom: 5px;"),
        tags$div(style = "width: 180px; height: 12px; background: linear-gradient(to right, #1f77b4, #d62728); border-radius: 3px; margin-bottom: 3px;"),
        tags$div(style = "display: flex; justify-content: space-between; width: 180px; font-size: 10px; color: #666;",
                 tags$span("Low"), tags$span("High")
        )
      ),
      
      # Escala em Círculos para o tamanho baseado no grau
      tags$div(
        style = "display: flex; flex-direction: column; align-items: center;",
        tags$div("Node size: Degree", style = "font-weight: bold; margin-bottom: 5px;"),
        tags$div(
          style = "display: flex; align-items: flex-end; gap: 15px;",
          tags$div(style = "display: flex; flex-direction: column; align-items: center;",
                   tags$div(style = "width: 8px; height: 8px; background: #888; border-radius: 50%; margin-bottom: 2px;"),
                   tags$span(minDegree, style = "font-size: 10px; color: #666;")
          ),
          tags$div(style = "display: flex; flex-direction: column; align-items: center;",
                   tags$div(style = "width: 15px; height: 15px; background: #888; border-radius: 50%; margin-bottom: 2px;"),
                   tags$span(avgDegree, style = "font-size: 10px; color: #666;")
          ),
          tags$div(style = "display: flex; flex-direction: column; align-items: center;",
                   tags$div(style = "width: 22px; height: 22px; background: #888; border-radius: 50%; margin-bottom: 2px;"),
                   tags$span(maxDegree, style = "font-size: 10px; color: #666;")
          )
        )
      )
    )
    })
    

# Alterar para usar a função createGraph  (Esta função ainda não foi implementada)      
    output$graph_ui_multiple <- renderUI({
      req(input$search_multiple_protein)
      
      edges <- graph_edges_rv()

      if (!is.data.frame(edges) || nrow(edges) == 0) {
        div(
          style = "padding: 8px 14px;
                   color: #b71c1c;
                   background-color: #fdecea;
                   border: 1px solid #f5c6cb;
                   border-radius: 6px;
                   font-size: 16px;
                   font-family: Arial, sans-serif;",
          "No interactions found !"
          
        )
      } else {
        nodes <- unique(c(edges$from, edges$to))
        nodes_df <- data.frame(id = nodes, label = nodes)
        
        visNetwork::visNetwork(nodes_df, edges, width = "100%", main = " ") %>%
          visNetwork::visEdges(arrows = "none") %>%
          visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
          visPhysics(stabilization = TRUE) %>%
          visNetwork::visLayout(randomSeed = 50)

              }
    })
    
    switch(
      ui_page(),
      home = htmlTemplate("www/home.html", !!!nav_buttons),
      downloadDB = htmlTemplate("www/download.html", !!!nav_buttons,
         go_downloadDB = actionButton("go_downloadDB", "Download DB")
      ),
      about = htmlTemplate("www/about.html", !!!nav_buttons),
      single_search = htmlTemplate(
        "www/single_search.html",
        !!!nav_buttons,
        multiple_search = nav_buttons$multiple_search,
        tissue_selector = uiOutput("tissue_selector_ui"),
        protein_selector = uiOutput("protein_selector_ui"),
        protein_count = textOutput("protein_count_ui"),
        protein_degree = textOutput("protein_degree_ui"),
        protein_name = textOutput("protein_name_ui"),
        protein_interactions = textOutput("protein_interactions_ui"),
        protein_closeness = textOutput("protein_closeness_ui"),
        protein_betweenness = textOutput("protein_betweenness_ui"),
        protein_diameter_ui = textOutput("protein_diameter_ui"),
        protein_geneid = textOutput("protein_geneid_ui"),
        protein_clustering = textOutput("protein_clustering_ui"),
        tissue_name = textOutput("tissue_name_ui"),
        tissue_interactions = textOutput("tissue_interactions_ui"),
        tissue_proteins = textOutput("tissue_proteins_ui"),
        tissue_diameter = textOutput("tissue_diameter_ui"),
        tissue_density = textOutput("tissue_density_ui"),
        tissue_clustering = textOutput("tissue_clustering_ui"),
        tissue_degree = textOutput("tissue_degree_ui"),
        graph_interactions = uiOutput("graph_ui_single"),
        graph_legend = uiOutput("graph_ui_legend"),
        downloadHTML = downloadButton("download_network", "Network (HTML)"),
        downloadEDGES = downloadButton("download_edges", "Edge List (CSV)"),
        downloadTISSUE = downloadButton("download_tissue_edges", "Tissue's Edge (CSV)"),
        
      ),
      help = htmlTemplate("www/help.html", !!!nav_buttons),
      multiple_search = htmlTemplate(
        "www/multiple_search.html",
        !!!nav_buttons,
        multiple_protein_selector = uiOutput("multiple_protein_selector_ui"),
        tissue_selector = uiOutput("tissue_selector_ui"),
        tissue_name = textOutput("tissue_name_ui"),
        tissue_interactions = textOutput("tissue_interactions_ui"),
        tissue_proteins = textOutput("tissue_proteins_ui"),
        tissue_diameter = textOutput("tissue_diameter_ui"),
        tissue_density = textOutput("tissue_density_ui"),
        tissue_clustering = textOutput("tissue_clustering_ui"),
        tissue_degree = textOutput("tissue_degree_ui"),
        graph_interactions = uiOutput("graph_ui_multiple"),
        downloadHTML = downloadButton("download_network", "Download Network (HTML)"),
        downloadEDGES = downloadButton("download_edges", "Download Edge List (CSV)")
      )
    )
  })
  
  # Observadores de navegação
  observeEvent(input$home, ui_page("home"))
  observeEvent(input$go_home, ui_page("home"))
  observeEvent(input$go_about, ui_page("about"))
  observeEvent(input$go_downloadDB, ui_page("downloadDB"))
  observeEvent(input$go_search, ui_page("search"))
  observeEvent(input$go_help, ui_page("help"))
  observeEvent(input$single_search, ui_page("single_search"))
  observeEvent(input$multiple_search, ui_page("multiple_search"))
  
  # Observador para pesquisa individual
  observeEvent(input$search_single_protein, {
    req(input$protein, input$tissue)
    response <- POST(
        url <- paste0(api_base_url, "/single_query"),
        body = list(
        prot = input$protein,
        tissue = input$tissue,
        layer = 1
        
      ),
      encode = "json"
    )
    
    if (http_type(response) != "application/json") {
      showNotification("Error accessing the API", type = "error")
      return()
    }
    
    data <- content(response, as = "parsed", simplifyDataFrame = TRUE)
    if (!is.null(data$error)) {
      showNotification(paste("Error:", data$error), type = "error")
      return()
    }
    
    datalist_rv(data)
    graph_edges_rv(data$all_interactions)
  }, ignoreInit = FALSE)
  
  # Observadores para pesquisa múltipla
  observeEvent(input$search_multiple_protein, {

    req(input$protein_multi, input$tissue)
    
    response <- POST(
      url <- paste0(api_base_url, "/multiple_query"),
      body = list(
        genes = paste(input$protein_multi, collapse = ","),
        tissue = input$tissue
      ),
      encode = "json"
    )
    
    if (http_type(response) != "application/json") {
      showNotification("Error accessing the API", type = "error")
      return()
    }
    
    data <- content(response, as = "parsed", simplifyDataFrame = TRUE)
    
    if (!is.null(data$error)) {
      showNotification(paste("Error:", data$error), type = "error")
      return()
    }
    datalist_rv(data)
    graph_edges_rv(data$multi_interactions)
    
    if (is.data.frame(graph_edges_rv())) {
      edges <- graph_edges_rv() %>%
        rename(from = geneSymbol1, to = geneSymbol2) %>%
        mutate(pair = ifelse(from < to, paste(from, to, sep = "_"), paste(to, from, sep = "_"))) %>%
        distinct(pair, .keep_all = TRUE) %>%
        select(from, to)
        graph_edges_rv(edges)
    } else {
            graph_edges_rv(data.frame(from = character(), to = character()))
    }
   })
}

shinyApp(ui, server)