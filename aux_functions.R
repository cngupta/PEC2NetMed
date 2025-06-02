

find_module_GO_enrichment=function(module_df,geneset_db,universe)
{
  GO_tbl=data.frame(label=character(),pval=numeric(),fdr=numeric(),
    signature=numeric(),geneset=numeric(),overlap=numeric(),background=numeric(),hits=character(),module=character())
  module_df=module_df[as.character(module_df$module) != "0",]
  modnames=unique(factor(module_df$module))
  for (j in 1:length(modnames))
  {
    mod.genes=module_df[module_df$module %in% modnames[j],]$gene
    hyp_obj = hypeR(mod.genes, geneset_db,background=universe,fdr=0.05)
    hyp_df =  hyp_obj$data
    message(nrow(hyp_df))
    if(nrow(hyp_df) > 0)
    {
			hyp_df$module=modnames[j]
      GO_tbl=rbind(GO_tbl, hyp_df)
    }
  }
  return(GO_tbl)
}


get_centrality=function(net, type, tag)
{
  colnames(net)=c("TF","TG","scores","link")
	net.igraph=graph_from_data_frame(net, directed = TRUE, vertices = NULL)
	if(type=="pr")
	{
		cent=page_rank(net.igraph)
		cent=as.data.frame(cent$vector)
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="closeness")
	{
		cent=closeness(net.igraph)
		cent=as.data.frame(cent)
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="betweenness")
	{
		cent=betweenness(net.igraph, directed=TRUE)
		cent=as.data.frame(cent)
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="hub_score")
	{
		cent=hub_score(net.igraph)
		cent=as.data.frame(cent$vector)
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="degree_in")
	{
		net.igraph=graph_from_data_frame(net[net$TG %in% net$TF,], directed = TRUE, vertices = NULL)
		cent=	as.data.frame(igraph::degree(net.igraph,mode="in"))
		colnames(cent)=c("scores")
		cent=cent[order(-cent$scores), , drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="degree_out")
	{
		cent=	as.data.frame(igraph::degree(net.igraph,mode="out"))
		colnames(cent)=c("scores")
		#cent=cent[order(-cent$scores),drop = FALSE]
		colnames(cent)=paste(tag,type,sep=".")
		cent$gene=rownames(cent)
		rownames(cent)=c()
		return(cent)
	}
	else if (type=="components")
	{
	 #no of connected components
		nCC=components(net.igraph)$no
		return(nCC)
	}
}


find_target_pairs_matrix=function(net) #network
{
	colnames(net)=c("TF","target","score")
	m=acast(net, TF~target, value.var="score")
	m=t(m)
	m[is.na(m)]=0 #set NA =0
	#find cardinalities
	# Find that paper and add reference
	i12 = m %*% t(m)
	s = diag(i12) %*% matrix(1, ncol = length(diag(i12)))
	u12 = s + t(s) - i12
	jacc= i12/u12
	jacc
}

get_coregnet_graph=function(net,th) #network, JI threshold, tag
{
	colnames(net)=c("TF","target","score")
	m=acast(net, TF~target, value.var="score")
	m=t(m)
	m[is.na(m)]=0 #set NA =0
	#find cardinalities
	# Find that paper and add reference
	i12 = m %*% t(m)
	s = diag(i12) %*% matrix(1, ncol = length(diag(i12)))
	u12 = s + t(s) - i12
	jacc= i12/u12
	grt=graph_from_adjacency_matrix(jacc,"undirected", weighted=T, diag=FALSE)
	tg=as.data.frame(get.edgelist(grt))
	tg$jaccard=E(grt)$weight
	tg=tg[tg$jaccard > th,]
	colnames(tg)=c("g1","g2","jaccard")
	tg
}

detect_modules = function(matrix, msize)
{
 #ref: http://pklab.med.harvard.edu/scw2014/WGCNA.html	#ref:
	dissMatrix = 1 - matrix
	# Call the hierarchical clustering function
	geneTree = flashClust(as.dist(dissMatrix),method="average");
	minModuleSize = msize;
	# Module identification using dynamic tree cut:
	#dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize)
	dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissMatrix, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize =msize)
	#dynamicColors = labels2colors(dynamicMods)
	modules=cbind(as.data.frame(dynamicMods),rownames(matrix))
	modules=modules[,c(2,1)]
	colnames(modules)=c("gene_name","module")
  modules$module=paste("M",modules$module,sep="_")
	modules
}
