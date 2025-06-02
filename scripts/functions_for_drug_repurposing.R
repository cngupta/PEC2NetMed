
# Function to detect delimiter
detect_delimiter = function(filepath) {
  # Read the first line of the file
  first_line = readLines(filepath, n = 1)

  # Check for common delimiters
  if (grepl("\t", first_line)) {
    return("\t")
  } else if (grepl(",", first_line)) {
    return(",")
  } else {
    stop("Unknown delimiter")
  }
}

# Function to read file with detected delimiter
read_delimited_file = function(filepath) {
  # Detect delimiter
  sep = detect_delimiter(filepath)

  # Read the file with the correct delimiter
  return(read.table(filepath, header = TRUE, sep = sep))
}

read_grn = function(file_path) {
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    net = read.table(file_path, header = TRUE, sep = ",")
  } else if (grepl("\\.txt$", file_path, ignore.case = TRUE) || grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    net = read.table(file_path, header = TRUE, sep = "\t")
  } else {
    stop("Unsupported file type. Please provide a CSV, TSV, or TXT file.")
  }
  return(net)
}


# Function to process GRN files and construct coregulatory graphs
process_grn_files = function(grn_path)
{
  tsv_files = list.files(path = grn_path, pattern = "\\.csv$", full.names = TRUE)
  graphs = list()
  for (path in tsv_files)
  {
    message("Processing file: ", path)
    net = read_delimited_file(path)
    net = net %>%
      mutate(TF = gsub("\\(\\+\\)|\\(-\\)", "", TF)) %>%
      filter(link != 'distal') %>%
      select(TF, target) %>%
      mutate(score = 1) %>%
      find_target_pairs_matrix() # Assuming this function exists and is appropriate
    inet=graph_from_adjacency_matrix(net,mode =c("undirected"),weighted = TRUE,diag=F)
    df=as_data_frame(inet)
    df=df[df$weight >0.5,]
    df$weight=1
    inet=graph_from_data_frame(df)
    g = simplify(inet, edge.attr.comb = "sum") # Simplify to combine multiple edges
    node = names(V(g))
    degree = degree(g, v = V(g))
    g.list=list(graph=g,node=node,degree=degree)
    graphs[[basename(path)]] = g.list
    message("Completed processing for ", path)
  }
  return(graphs)
}

# for abalation
process_grn_files_abalation = function(grn_path,topn) # select top n % edges
{
  tsv_files = list.files(path = grn_path, pattern = "\\.csv$", full.names = TRUE)
  graphs = list()
  for (path in tsv_files)
  {
    message("Processing file: ", path)
    net = read_delimited_file(path)
    net = net[order(-net$importance), ]
    net = net[1:ceiling(nrow(net) * topn), ]

    net = net %>%
      mutate(TF = gsub("\\(\\+\\)|\\(-\\)", "", TF)) %>%
      filter(link != 'distal') %>%
      select(TF, target) %>%
      mutate(score = 1) %>%
      find_target_pairs_matrix() # Assuming this function exists and is appropriate

    inet=graph_from_adjacency_matrix(net,mode =c("undirected"),weighted = TRUE,diag=F)
    df=as_data_frame(inet)
    df=df[df$weight >0.5,]
    df$weight=1
    inet=graph_from_data_frame(df)
    g = simplify(inet, edge.attr.comb = "sum") # Simplify to combine multiple edges
    node = names(V(g))
    degree = degree(g, v = V(g))
    g.list=list(graph=g,node=node,degree=degree)
    graphs[[basename(path)]] = g.list
    message("Completed processing for ", path)
  }
  return(graphs)
}

#Following functions are from the SaveRunner package

computeMinimum = function(distance_matrix,dim=1,perc_thr=5)
{
  # dim = 1 -> row mean (ie, from DrugTarget to DiseaseGene),
  # dim = 2 -> col mean (ie, from DiseaseGene to DrugTarget)
  minimum = apply(distance_matrix,dim,min,na.rm = T)
  count = length(minimum[is.infinite(minimum)])
  perc = ( count / length(minimum) ) * 100
  if( perc > perc_thr )
  {
    minimum = minimum
  }
  else
  {
    minimum = minimum[is.finite(minimum)]
  }
}

computeDegreeDistribution = function(list,graph_info)
{
  graph = graph_info$graph
  node = graph_info$node
  from = which(node %in% list)
  d = degree(graph, v = V(graph)[from])
  t = table(d)
  degree_sorted = as.numeric(names(t))
  freq = as.numeric(t)
  df = data.frame(degree=degree_sorted,frequency=freq)
  return(df)
}

computeRandomProximity = function(targets_degree_distribution,genes_degree_distribution,graph_info,iter=100)
{
  random_proximity = NULL
  for (i in 1:iter)
  {
    random_node1 = selectRandomNodes(targets_degree_distribution,graph_info)
    random_node2 = selectRandomNodes(genes_degree_distribution,graph_info)
    tmp = computeProximity(random_node1,random_node2,graph_info)
    random_proximity = c(random_proximity,tmp)
  }
  random_proximity[is.infinite(random_proximity)] = NaN
  return(random_proximity)
}

selectRandomNodes = function(degree_distribution,graph_info)
{
  node = graph_info$node
  degree = graph_info$degree
  d = degree_distribution$degree
  freq = degree_distribution$freq
  random_node = unlist(lapply(d, function(x){
    sample(node[degree == x], freq[d == x])
  }))
  return(random_node)
}

computeProximity = function(list1,list2,graph_info)
{
  graph = graph_info$graph
  node = graph_info$node
  from = which(node %in% list1)
  to = which(node %in% list2)
  distance_matrix = distances(graph, v = V(graph)[from], to = V(graph)[to])
  minimum = computeMinimum(distance_matrix)
  proximity = mean(minimum, na.rm = T)
  return(proximity)
}
computeStatistics = function(distribution,observation,plot=F,filename=NULL)
{
  m = mean(distribution, na.rm = T)
  sd = sd(distribution, na.rm = T)
  z = (observation - m) / sd
  pval = pnorm(z, lower.tail = T) # H1 = be less than expected by chance
  # ######################################
  # the xpd parameter is "A logical value or NA.
  # If FALSE, all plotting is clipped to the plot region
  # if TRUE, all plotting is clipped to the figure region
  # if NA, all plotting is clipped to the device region"

  if(plot){

    par(xpd=F)

    xmin = min((m-3*sd),observation)
    xmax = max((m+3*sd),observation)

    fun = function(x) dnorm(x,m,sd)

    pdf(file = filename)

    plot = curve(fun, xlim = c(xmin,xmax),
                  xlab = "Network proximity",
                  ylab = "Probability density",
                  n = 1000)

    polygon(plot, col = "grey")

    abline(v = observation, lty = 2, lwd = 2, col = "red")

    legend("topright",
           paste("Observation =", format(observation, digits = 3), "\n",
                 "p-value =", format(pval, digits=1)),
           text.col=2,
           text.font = 2,
           bty="o",
           xjust = 0.5,
           yjust = 0.5)

    dev.off()

  }
  return(list(pval,z))
}
