# plumber.R
#library(vroom)
library(dplyr)
library(magrittr)
library(DBI)
library(RSQLite)
library(glue)
library(igraph)
library(pool)

# 1. Pool criado para gerenciar as conexões com o banco (inicia quando a API sobe no servidor)
#pool <- dbPool(
#  drv = RSQLite::SQLite(),
#  dbname = "/app/data/interactions.sqlite",
#  onCreate = function(con) {
#    dbExecute(con, "PRAGMA journal_mode = WAL;")     # Permite leitura e escrita simultânea
#    dbExecute(con, "PRAGMA synchronous = NORMAL;")   # Reduz a frequência de gravação no disco
#    dbExecute(con, "PRAGMA cache_size = -131072;")   # 128MB
#    dbExecute(con, "PRAGMA mmap_size = 1073741824;") # 1GB
#    dbExecute(con, "PRAGMA busy_timeout = 5000;")    # Evita travamento da aplicação em caso de tentativa de escrita simultânea
#  }
#)

# Garante que o pool feche quando o servidor Plumber for encerrado
#*@plumber
function(pr) {
  pr$registerHook("exit", function() {
    poolClose(pool)
  })
}
db_name <- "/app/data/interactions.sqlite"

database_connect <- function() {
  con <- dbConnect(RSQLite::SQLite(), db_name)
  
  # Pragmas essenciais para ambientes web concorrentes
  dbExecute(con, "PRAGMA journal_mode = WAL;")
  dbExecute(con, "PRAGMA synchronous = NORMAL;")
  dbExecute(con, "PRAGMA busy_timeout = 5000;") # Espera até 5s se o DB estiver ocupado
  dbExecute(con, "PRAGMA foreign_keys = ON;")
  dbExecute(con, "PRAGMA cache_size = -131072;")   # 128MB
  dbExecute(con, "PRAGMA mmap_size = 3221225472;") # 3GB
  
  return(con)
}

extract_subgraph_by_degree <- function(edges, seed_protein, layer) {
  if (nrow(edges) == 0 || !seed_protein %in% unique(c(edges$from, edges$to))) {
    return(data.frame(from = character(), to = character()))
  }
  
  # Criar o grafo não direcionado
  g <- graph_from_data_frame(edges, directed = FALSE)
  # Remover arestas duplicadas e loops
  g <- simplify(g)
  
  # Encontrar todos os vértices até o grau especificado
#  if (degree > 0) {
    # Obter distâncias do vértice inicial
#    distances <- distances(g, v = seed_protein)
    
    # Filtrar vértices dentro do grau especificado
#    valid_vertices <- names(which(distances[1, ] <= degree & distances[1, ] > 0))
#    valid_vertices <- c(seed_protein, valid_vertices)
  if (seed_protein %in% V(g)$name) {
    
    # Descobrir o ID numérico ou nome do vértice 'palio'
    v_protein <- seed_protein
    
    # Obter todos os nós vizinhos até a distância 2
    valid_vertices <- neighborhood(g, order = layer, nodes = v_protein)[[1]]    

    # Extrair o subgrafo induzido
    if (length(valid_vertices) > 1) {
      subg <- induced_subgraph(g, valid_vertices)
      subg_edges <- as_data_frame(subg)
      return(subg_edges)
    } else {
      return(data.frame(from = character(), to = character()))
    }
  } else {
    # Grau 0: apenas o nó inicial
    return(data.frame(from = character(), to = character()))
  }
}

#* Gerar a lista de códigos de proteinas
#* @get /protein_list
function() {
  con <- database_connect()
  on.exit(dbDisconnect(con), add = TRUE) # Garante o fechamento seguro da conexão
  query <- glue_sql("SELECT geneSymbol, idProtein FROM proteins ORDER BY geneSymbol", .con = con)
  protein_list <- dbGetQuery(con, query)
  return(protein_list)
}

#* Gerar a lista de códigos dos tecidos
#* @get /tissue_list
function() {
  con <- database_connect()
  on.exit(dbDisconnect(con), add = TRUE) # Garante o fechamento seguro da conexão
  query <- glue_sql("SELECT idTissue,descTissue FROM tissues ORDER BY descTissue", .con = con)
  tissue_list <- dbGetQuery(con, query)
  
  return(tissue_list)
}

#* Baixar o banco todo ou a rede de um tecido
#* @get /download_edges
function(sqltxt) {
  con <- database_connect()
  on.exit(dbDisconnect(con), add = TRUE) # Garante o fechamento seguro da conexão
  query <- glue_sql({sqltxt}, .con = con)
  edges_list <- dbGetQuery(con, query)
  return(edge_list)
}

#* Receber uma lista de proteínas e tecido, e retornar subconjunto para criar o grafo
#* @param genes A lista de genes (separados por virgula)
#* @param tissue O tecido
#* @post /multiple_query
function(genes, tissue) {
  con <- database_connect()
  on.exit(dbDisconnect(con), add = TRUE) # Garante o fechamento seguro da conexão
  genes_list <- unlist(strsplit(genes, ",")) 
  genes_list <- trimws(genes_list)
    query <- glue_sql(
    "SELECT geneSymbol1, geneSymbol2, idTissue
     FROM interactions
     WHERE idTissue = {tissue} 
     AND (geneSymbol1 IN ({genes_list*}) AND geneSymbol2 IN ({genes_list*}))",
    .con = con
  )

  multi_interactions <- dbGetQuery(con, query)

  query <- glue_sql("SELECT m.*, t.* FROM networkMetrics m, tissues t 
                    WHERE m.idTissue = {tissue} AND t.idTissue = {tissue}", .con = con)
  
  calculated_tissue_metrics <- dbGetQuery(con, query)
  
  return(list(multi_interactions = multi_interactions, tissueMetrics = calculated_tissue_metrics))
}

#* Buscar todas as interações de um tecido e retornar o dados para processamento no cliente
#* @param prot Proteína inicial (para cálculo de métricas)
#* @param layer Grau máximo de interação (1-4)
#* @param tissue Nome do tecido (ex: "cortex")
#* @post /single_query
function(prot, tissue, layer) {
  con <- database_connect()
  on.exit(dbDisconnect(con), add = TRUE) # Garante o fechamento seguro da conexão
  layer <- as.integer(layer)
  # Buscar todas as interações do tecido
  query <- glue_sql("
    SELECT geneSymbol1, geneSymbol2, idTissue
    FROM interactions
    WHERE idTissue = {tissue}", .con = con)

  all_interactions <- dbGetQuery(con, query)
  
  # Processar o grafo com igraph
  if (is.data.frame(all_interactions) && nrow(all_interactions) > 0) {
      edges <- all_interactions %>%
      rename(from = geneSymbol1, to = geneSymbol2) %>%
      mutate(pair = ifelse(from < to, paste(from, to, sep = "_"), paste(to, from, sep = "_"))) %>%
      distinct(pair, .keep_all = TRUE) %>%
      select(from, to)
  
  # Extrair o subgrafo baseado no grau de interação
      subgraph_edges <- extract_subgraph_by_degree(edges, prot, 
                                                     as.integer(layer))
      all_edges <- edges
  } else {
        subgraph_edges <- data.frame(from = character(), to = character())
        all_edges <- data.frame(from = character(), to = character())
        
  }
  
  # Obter os nomes das proteínas
  vProteins <- unique(c(subgraph_edges$from, subgraph_edges$to))

  # Buscar as metricas calculadas paras as proteinas
  query <- glue_sql("
    SELECT m.*, p.*
    FROM proteinMetrics m, proteins p
    WHERE p.geneSymbol IN ({vProteins*}) 
      AND m.idTissue = {tissue} 
      AND m.idProtein = p.idProtein", .con = con)
  
  all_proteins_metrics <- dbGetQuery(con, query)
  
  # Buscar as métricas calculadas para uma proteínas especifica
  
  query <- glue_sql("
    SELECT m.*, p.*
    FROM proteinMetrics m, proteins p
    WHERE p.geneSymbol = {prot} 
      AND m.idTissue = {tissue} 
      AND m.idProtein = p.idProtein
  ", .con = con)
  
  selected_proteins_metrics <- dbGetQuery(con, query)

  # Buscar as métricas calculadas para o tecido
  
  query <- glue_sql("SELECT m.*, t.* FROM networkMetrics m, tissues t 
                    WHERE t.idTissue = {tissue}
                    AND m.idTissue = {tissue}", .con = con)
  calculated_tissue_metrics <- dbGetQuery(con, query)
  
  return(list(
    all_interactions = subgraph_edges,
    all_tissue_edges = all_edges,
    proteinMetrics = selected_proteins_metrics,
    allProteinsMetrics = all_proteins_metrics,
    tissueMetrics = calculated_tissue_metrics
  ))
}
