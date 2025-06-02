rm(list=ls())
source('load_libraries.R')
source('aux_functions.R')

cohorts=c("CMC","UCLA_ASD","Urban_DLPFC")
centrality="betweenness"

grn_path = "data/PEC2_NetMed_GRNs/"

# Read autism gene data
asd.risk = read.csv("data/SFARI-Gene_genes_01-23-2023release_03-02-2023export.csv", sep = ",") %>%
  filter(gene.score == 1 | syndromic == 1) %>%
  .$gene.symbol %>%
  unique()

#from dis gene net db
dis=read.table("data/all_gene_disease_associations.tsv",header=T,sep="\t",quote="")
scz.risk= dis[dis$diseaseName %like% "Schizophrenia" & (dis$source %like% "CTD_human" ),"geneSymbol"]
bp.risk= dis[dis$diseaseName %like% "Bipolar Disorder" & (dis$source %like% "CTD_human" ),"geneSymbol"]

all_dis_genes=union(asd.risk,union(scz.risk,bp.risk))

outdegree_df_all = data.frame(Network = character(), gene = character(), degree = numeric(),celltype=character())
outdegree_df_top = data.frame(Network = character(), gene = character(), degree = numeric(),celltype=character())

calculate_delta=function(dis,ctrl)
{
  #
  all_tfs=union(dis$gene,ctrl$gene)

  delta = data.frame(gene = character(), delta = numeric(),celltype=character())
  for(gene in all_tfs)
  {
    for(cell in unique(dis$celltype))
    {
      dis.degree=dis[dis$gene == gene & dis$celltype == cell,"degree"]
      ctrl.degree=ctrl[ctrl$gene == gene & ctrl$celltype == cell,"degree"]
      diff=log2(dis.degree/ctrl.degree)
      if (length(diff)>0)
      {
        df= data.frame(gene = gene, delta = diff,celltype=cell)
        delta=rbind(delta,df)
      }
    }
  }
  df=delta[abs(delta$delta)>=0,]
  selgenes=unique(df$gene)
  df=delta[delta$gene %in% selgenes ,]
  mat1=acast(df, celltype ~ gene, value.var="delta", fun.aggregate=sum)
  return(mat1)
}


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
      print(pheno_grn_path)
      # Loop through GRN files in the directory
      for (filename in list.files(pheno_grn_path))
      {
          print(filename)
          ptm = proc.time()
          message(paste("working on",filename,sep="::"))
          tag =  sub("^(.*?)_.*$", "\\1", filename)
          network_data = read.table(paste(pheno_grn_path,filename,sep="/"),header=T,sep="\t")
          df=get_centrality(network_data,centrality,paste(cohort,pheno,sep="_"))
          df$Network=paste(cohort,pheno,sep="_")
          colnames(df)=c("degree","gene","Network")
          df=df[,c("Network","gene","degree")]
          df$celltype=tag
          df=df[order(-df$degree),]
          df=df[df$degree > 0,]
          outdegree_df_all=rbind(outdegree_df_all,df)
          df=df[1:(nrow(df)*0.1),]
          outdegree_df_top=rbind(outdegree_df_top,df)
      }
    }
  }
}

# Calculate delta
dis_df=outdegree_df_all[outdegree_df_all$Network %in% c("UCLA_ASD_ASD"),]
ctrl_df=outdegree_df_all[outdegree_df_all$Network %in% c("UCLA_ASD_ctrl"),]
mat1=calculate_delta(dis_df,ctrl_df)
write.csv(mat1,"tables/ASD_delta_influence_gt3.csv",quote=F)


dis_df=outdegree_df_all[outdegree_df_all$Network %in% c("CMC_SCZ"),]
ctrl_df=outdegree_df_all[outdegree_df_all$Network %in% c("CMC_ctrl"),]
mat2=calculate_delta(dis_df,ctrl_df)
write.csv(mat2,"tables/CMC_SCZ_delta_influence_gt3.csv",quote=F)

dis_df=outdegree_df_all[outdegree_df_all$Network %in% c("Urban_DLPFC_BPD"),]
ctrl_df=outdegree_df_all[outdegree_df_all$Network %in% c("Urban_DLPFC_ctrl"),]
mat3=calculate_delta(dis_df,ctrl_df)
write.csv(mat3,"tables/Urban_DLPFC_BPD_delta_influence_gt3.csv",quote=F)

all_columns = Reduce(union, list(colnames(mat1), colnames(mat2), colnames(mat3)))

add_missing_columns = function(matrix, all_columns) {
  # Determine missing columns and create them with NA values
  missing_columns = setdiff(all_columns, colnames(matrix))
  if (length(missing_columns) > 0) {
    new_cols = matrix(NA, nrow=nrow(matrix), ncol=length(missing_columns))
    colnames(new_cols) = missing_columns

    matrix = cbind(matrix, new_cols)
  }
  matrix = matrix[, all_columns, drop = FALSE]
  return(matrix)
}

mat1.n = add_missing_columns(mat1, all_columns)
mat2.n = add_missing_columns(mat2, all_columns)
mat3.n = add_missing_columns(mat3, all_columns)

mat1.n[is.na(mat1.n)]=0
mat2.n[is.na(mat2.n)]=0
mat3.n[is.na(mat3.n)]=0


de=read.table("data/ASD_DEGcombined.csv",sep=",",header=T)
filt=de[(de$gene %in% colnames(mat1)),]
mat3 = acast(filt, cell_type ~ gene, value.var="log2FoldChange", fun.aggregate=sum)

missing_genes=setdiff(colnames(mat1),colnames(mat3))
# add missing genes
for(gene in missing_genes)
{
  mat3 = cbind(mat3, setNames(rep(0, nrow(mat3)), gene))
  colnames(mat3)[ncol(mat3)] = gene
}

# for main fig
rownames.cex = ifelse(rownames(t(mat3.n)) %in% all_dis_genes, 0.7, 0.01)
rownames.col = ifelse(rownames(t(mat3.n)) %in% bp.risk, "#e18727ff",
  ifelse(rownames(t(mat3.n)) %in% asd.risk,"#ee4c97ff","#0072b5ff"))

circlize_plot = function(){
circos.par(start.degree = 90, gap.degree = 15)
circos.heatmap.initialize(t(mat3.n))
circos.heatmap(t(mat3.n),col = col_fun1,rownames.side = "outside",rownames.cex=0.5,
rownames.col="black",cluster = FALSE,
  bg.border = "#e18727ff", bg.lwd = 3, bg.lty = 2,cell.lwd=0.1)
#for labeling cell types
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
         if(CELL_META$sector.numeric.index == 1) { # the first sector
             cn = colnames(t(mat1))
             n = length(cn)
             circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                 1:n + 1, cn,
                 cex = 1, adj = c(0, 1.8), facing = "inside")
         }
     }, bg.border = NA)
set_track_gap(cm_h(0.25)) # 0.5cm
circos.heatmap(t(mat2.n),col = col_fun3,bg.border = "#0072b5ff", bg.lwd = 3, bg.lty = 2,cluster = FALSE,
cell.lwd=0.1)
circos.heatmap(t(mat1.n),col = col_fun3,bg.border = "#ee4c97ff", bg.lwd = 3, bg.lty = 2,cluster = TRUE,
cell.lwd=0.1)
circos.clear()
}

pdf("2D.pdf",height=12,width=12)
col_fun1 = colorRamp2(c(-1, 0, 1), c("#ee0000ff", "white", "#3b4992ff"))
col_fun2 = colorRamp2(c(-1, 0, 1), c("#ee0000ff", "white", "#3b4992ff"))
col_fun3 = colorRamp2(c(-1, 0, 1), c("#ee0000ff", "white", "#3b4992ff"))

lgd_lines = Legend(at = c("BPD", "SCZ", "ASD"), type = "lines",
    legend_gp = gpar(col = c("#e18727ff","#0072b5ff","#ee4c97ff"), lwd = 2,fontsize = 24), title_position = "topleft",
    title = "Tracks",title_gp = gpar(fontsize = 18))

# continuous
lgd_links = Legend(at = c( -1, 0, 1), col_fun = col_fun1,
    title_position = "topleft", title = "DI score",title_gp = gpar(fontsize = 18,
    labels_gp = gpar(fontsize = 18))) #title = expression(Delta ~ "TF influence"))

lgd_list = packLegend(lgd_lines, lgd_links,direction = "horizontal")
circlize_plot()
draw(lgd_list, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

dev.off()
