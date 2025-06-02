#read GRNs, find coreg matrix, find modules

rm(list=ls())
set.seed(123)
source('load_libraries.R')
source('aux_functions.R')
library(clusterProfiler)
#set min mod size
msize=10

library(aricode)
library(ggnewscale)
library(org.Hs.eg.db)


#funciton for NMI
get_nmi=function(dis.df,ctrl.df)
{
  colnames(dis.df)=c("gene_name","disease")
  colnames(ctrl.df)=c("gene_name","ctrl")
  common=intersect(dis.df$gene_name,ctrl.df$gene_name)
  ref=ctrl.df[ctrl.df$gene_name %in% common,]
  test=dis.df[dis.df$gene_name %in% common,]
  df=merge(test,ref)
  nmi=NMI(df$ctrl, df$disease)
  return(nmi)
}

# Function to perform enrichment for each cluster
enrich_cluster_hyper = function(cluster_genes, disease_genes, universe) {
  # Calculate the overlap
  overlap = length(intersect(cluster_genes, disease_genes))
  # Total genes in the universe
  total_genes = universe
  # Genes in the disease gene set
  disease_set_size = length(disease_genes)
  # Genes in the cluster
  cluster_size = length(cluster_genes)

  # Perform hypergeometric test
  p_value = phyper(overlap - 1, disease_set_size, total_genes - disease_set_size, cluster_size, lower.tail = FALSE)

  # Return the result
  return(data.frame(cluster_size = cluster_size, overlap = overlap, p_value = p_value))
}


get_gmt_for_magma=function(mod_df)
{
  df=mod_df
  df$entrez_id = mapIds(org.Hs.eg.db,
                       keys = df$gene_name,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")

  df = df %>% filter(!is.na(entrez_id))

  gmt_list = df %>%
   group_by(module) %>%
   summarise(genes = paste(entrez_id, collapse = "\t")) %>%
   ungroup() %>%
   mutate(description = module) %>%
   dplyr::select(module, description, genes)

  return(gmt_list)
}

get_gmt_for_gene_symbols=function(mod_df)
{
  df=mod_df
  gmt_list = df %>%
   group_by(module) %>%
   summarise(genes = paste(gene_name, collapse = "\t")) %>%
   ungroup() %>%
   mutate(description = module) %>%
   dplyr::select(module, description, genes)

  return(gmt_list)
}


#Read GO data
#data=GSA.read.gmt('data/Human_GO_bp_with_GO_iea_symbol.gmt')
#genesets=data$genesets
#names(genesets)=data$geneset.descriptions

# Read autism gene data
autism_genes = read.csv("data/SFARI-Gene_genes_01-23-2023release_03-02-2023export.csv", sep = ",") %>%
  filter(gene.score == 1 | syndromic == 1) %>%
  .$gene.symbol %>%
  unique()

dis=read.table("data/all_gene_disease_associations.tsv",header=T,sep="\t",quote="")
scz= dis[dis$diseaseName %like% "Schizophrenia" & (dis$source %like% "CTD_human" | dis$source %like% "GWASCAT"),"geneSymbol"]
bp= dis[dis$diseaseName %like% "Bipolar Disorder" & (dis$source %like% "CTD_human" | dis$source %like% "GWASCAT"),"geneSymbol"]

cohorts=c("UCLA_ASD","CMC","Urban_DLPFC")

grn_path = "data/PEC2_NetMed_GRNs"

for(cohort in cohorts)
{
  phenos = list.dirs(paste(grn_path,cohort,sep="/"))
  print(phenos)
  for(i in 1:length(phenos))
  {
    if(basename(phenos[i]) != cohort)
    {
      pheno = basename(phenos[i])
      print(pheno)
      # Directory path containing the GRN files
      pheno_grn_path = paste(grn_path,cohort,sep="/")
      pheno_grn_path = paste(pheno_grn_path,pheno,sep="/")
      # Loop through GRN files in the directory
      for (filename in list.files(pheno_grn_path))
      {
          print(filename)
          ptm = proc.time()
          message(paste("working on",filename,sep="::"))
          tag =  sub("^(.*?)_.*$", "\\1", filename)
          net = read.table(paste(pheno_grn_path,filename,sep="/"),header=T,sep="\t")
          net=net %>% filter(link != 'distal')
          net=net[,c(1,2)]
          net$score=1

          mat=find_target_pairs_matrix(net)
          m=detect_modules(mat,msize)
          #enrich_df=find_module_GO_enrichment(m,genesets,7200)
          name=paste(cohort,pheno,sep="::")
          name=paste(name,tag,sep="::")
          #assign(name,enrich_df)
          name=paste(name,"modules",sep="::")
          assign(name,m)
      }
    }
  }
}

# write GO to tbl
#pattern=paste0("^(", paste(cohorts, collapse = "|"), ")")
#list=grep(pattern, ls(), value = TRUE)
#modified_dfs = list()
#for(item in list)
#{#
#  df=get(item)
#  if(nrow(df)>0)
#  {
#    df$celltype=item
#    modified_dfs[[item]] = df

#  }
#}
#combined_df = do.call(rbind, modified_dfs)
#write.table(combined_df,file="tables/WGCNA_modules_GO_enrichment_merged.csv",quote=F)


#write modules to tbl
list=grep("*::modules", ls(), value = TRUE)
modified_dfs = list()
for(item in list)
{
  df=get(item)
  if(nrow(df)>0)
  {
    df$celltype=gsub("::modules","",item)
    modified_dfs[[item]] = df
  }
}
combined_df = do.call(rbind, modified_dfs)
combined_df$moduleID=paste(combined_df$module,combined_df$Celltype,sep="::")
write.table(combined_df,file="tables/WGCNA_modules_merged.csv",quote=F)

#get NMI for ctrl vs disease
celltypes=c("astro","endo","micro","oligo","opc","vlmc","excitatory","inhibitory")
cohorts=c("Urban_DLPFC","UCLA_ASD","CMC")
nmi_df=data.frame(celltype=character(),cohort=character(),disorder=character(),nmi=numeric())
nmod_df=data.frame(celltype=character(),cohort=character(),disorder=character(),nmod=numeric())
ndismod_df=data.frame(celltype=character(),cohort=character(),disorder=character(),ndismod=numeric())
enrichment_results_tbl=data.frame(module=character(),cluster_size=numeric(),
  overlap=numeric(),p_value=numeric(),cohort=character(),disorder=character(),celltype=character())


for (cohort in cohorts)
{
  if (cohort == "UCLA_ASD")
  {
    pheno=c("ASD")
  }
  if (cohort == "CMC")
  {
    pheno=c("SCZ")
  }
  if (cohort == "Urban_DLPFC")
  {
    pheno=c("BPD")
  }
  for (celltype in celltypes)
  {
    #
    name=paste(cohort,pheno,sep="::")
    name=paste(name,celltype,sep="::")
    name=paste(name,"::modules",sep="")
    pheno_mod=get(name)

    # test for disease genes
    if(pheno == "ASD")
    {
      disease_genes=autism_genes
    }
    if(pheno == "BPD")
    {
      disease_genes=bp
    }
    if(pheno == "SCZ")
    {
      disease_genes=scz
    }
    universe_genes = unique(pheno_mod$gene_name)
    enrichment_results = lapply(split(pheno_mod$gene_name, pheno_mod$module), enrich_cluster_hyper, disease_genes, 7000)
    enrichment_results_df = bind_rows(enrichment_results, .id = "module")
    enrichment_results_df$cohort=cohort
    enrichment_results_df$disorder=pheno
    enrichment_results_df$celltype=celltype
    enrichment_results_tbl=rbind(enrichment_results_tbl,enrichment_results_df)
    n = nrow(enrichment_results_df %>% filter(p_value <= 0.05))

    df=data.frame(celltype=celltype,cohort=cohort,disorder=pheno,ndismod=n)
    ndismod_df=rbind(ndismod_df,df)

    # NMI
    name=paste(cohort,"ctrl",sep="::")
    name=paste(name,celltype,sep="::")
    name=paste(name,"::modules",sep="")
    ctrl_mod=get(name)
    nmi=get_nmi(pheno_mod,ctrl_mod)
    df=data.frame(celltype=celltype,cohort=cohort,disorder=pheno,nmi=nmi)
    nmi_df=rbind(nmi_df,df)
    df=data.frame(celltype=celltype,cohort=cohort,disorder=pheno,nmod=length(unique(pheno_mod$module)))
    nmod_df=rbind(nmod_df,df)

    #GMT with gene symbols
    pheno_gmt=get_gmt_for_gene_symbols(pheno_mod)
    name=paste(cohort,pheno,sep="::")
    name=paste(name,celltype,sep="::")
    name=paste(name,".gmt",sep="")
    name=paste("tables/",name,sep="/")
    write.table(pheno_gmt, file = name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    ctrl_gmt=get_gmt_for_gene_symbols(ctrl_mod)
    name=paste(cohort,"ctrl",sep="::")
    name=paste(name,celltype,sep="::")
    name=paste(name,".gmt",sep="")
    name=paste("tables/",name,sep="/")
    write.table(ctrl_gmt, file = name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    #GMT for magma
    pheno_gmt=get_gmt_for_magma(pheno_mod)
    name=paste(cohort,pheno,sep="::")
    name=paste(name,celltype,sep="::")
    name=paste(name,".gmt",sep="")
    name=paste("tables/",name,sep="/")
    write.table(pheno_gmt, file = name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    ctrl_gmt=get_gmt_for_magma(ctrl_mod)
    name=paste(cohort,"ctrl",sep="::")
    name=paste(name,celltype,sep="::")
    name=paste(name,".gmt",sep="")
    name=paste("tables/",name,sep="/")
    write.table(ctrl_gmt, file = name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}
main.nimi_df=nmi_df
write.table(enrichment_results_tbl, file = "tables/risk_gene_enrichment_results.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

mat=t(acast(main.nimi_df[,c(1,3,4)], celltype ~ disorder, value.var = "nmi"))
scale_vector = function(x) {
  (x - min(x)) / (max(x) - min(x))
}
mat = t(apply(mat, 1, scale_vector))
nmi_df=melt(mat)
colnames(nmi_df)=c("disorder","celltype","nmi")

mat=t(acast(nmod_df[,c(1,3,4)], celltype ~ disorder, value.var = "nmod"))
nmod_df=melt(mat)
colnames(nmod_df)=c("disorder","celltype","nmod")

ndismod_df[is.na(ndismod_df)]=0
mat=t(acast(ndismod_df[,c(1,3,4)], celltype ~ disorder, value.var = "ndismod"))
ndismod_df=melt(mat)
colnames(ndismod_df)=c("disorder","celltype","ndismod")

df=merge(nmi_df,nmod_df,by=c("celltype","disorder"))
df=merge(df,ndismod_df,by=c("celltype","disorder"))

df$disorder = factor(df$disorder, levels = c("BPD", "SCZ", "ASD"))

p=ggplot(df, aes(y = celltype, x = disorder)) +
  geom_tile(aes(fill = nmi), color = "black") +
  scale_fill_gradient2(low = "#1B191999", mid = "#80818099", high = "white", midpoint = median(df$nmi),
                       breaks = c(0, 0.4, 1),
                       labels = c("Low", "", "High"), name = "NMI \n(ctrl vs dis.)") +
  new_scale_fill() +
  geom_point(aes(size = ndismod), fill = "grey", color = "black", shape = 22) +
  scale_size(name = "no. disease \n modules") +
  theme_minimal(base_size = 10) +
  labs(title = "", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = unit(0.25, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.title.position = "top",
    legend.position = "left",
    legend.text = element_text(size = 8),  # Adjust legend text size
    legend.title = element_text(size = 8)
  )

ggsave(p,file="3A.pdf",width=2.75, height=3.3,units=c("in"))
