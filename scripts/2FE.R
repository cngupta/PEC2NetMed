#rm(list=ls())
source('load_libraries.R')
library(GSEABase)

# for ex neur
Exsub=c('L2.3.IT', 'L4.IT', 'L5.IT', 'L5.ET', 'L5.6.NP','L6b','L6.IT','L6.CT','L6.IT.Car3')

# DI scores
#scz=read.csv("tables/Urban_DLPFC_BPD_delta_influence_gt3.csv")
scz=read.csv("tables/CMC_SCZ_delta_influence_gt3.csv")
#scz=read.csv("tables/ASD_delta_influence_gt3.csv")


scz=melt(scz)
scz=scz %>% filter (X == "excitatory")
colnames(scz)=c("celltype","gene","DI")


# deg FC
#degdf=read.csv("data/ASD_DEGcombined.csv")
degdf=read.csv("data/Schizophrenia_DEGcombined.csv")

degscz <- degdf %>% dplyr::filter(cell_type %in% Exsub) %>% dplyr::select(gene,log2FoldChange,cell_type)
#get average FC across ex
degscz <- degscz %>% group_by(gene) %>%
  summarise( average_log2FoldChange = mean(log2FoldChange, na.rm = TRUE))


# BBB+ drugs
b3db = read.csv("data/B3DB_classification.tsv", sep = "\t")
bbb=b3db[b3db[,7] == "BBB+",]

# LINCS 1000 from mayan lab
#coef = fread("data/cp_mean_coeff_mat.tsv", sep = "\t", header = TRUE)
metadata = read.delim("data/LINCS_small_molecules.tsv", header = TRUE, sep = "\t", fill = TRUE, quote = "")

# common bet bbb and lincs
bbb_lincs=intersect(coef$V1,bbb$compound_name)
coef.bbb=as.data.frame(coef[coef$V1 %in% bbb_lincs,])
rownames(coef.bbb)=coef.bbb$V1
coef.bbb$V1=NULL


# merge FC and DI
scz_di_fc <- scz %>% left_join(degscz, by = "gene")

#merge DI FC and lincs
coef.bbb=coef.bbb[,colnames(coef.bbb) %in% scz_di_fc$gene]
coef.bbb=melt(as.matrix(coef.bbb))
colnames(coef.bbb)=c("compound_name","gene","DrugDE")
scz_di_fc = scz_di_fc %>% left_join(coef.bbb,by="gene")


# filter drugs and TFs
scz_di_fc=scz_di_fc %>%  filter(average_log2FoldChange * DrugDE < 0 , abs(average_log2FoldChange) > 0 | abs(DrugDE) >0)
scz_di_fc=scz_di_fc %>%  filter(abs(DI) >2)
scz_di_fc <- scz_di_fc %>% filter(complete.cases(.))

scz_di_fc$celltype=NULL


####################
#### plot a heatmap
mat=acast(scz_di_fc,gene ~ compound_name,value.var="DrugDE")
mat[is.na(mat)]=0

# add dis fc as row annotations
lfc_data <- scz_di_fc %>%
  group_by(gene) %>%
  summarize(avg_log2FoldChange = mean(average_log2FoldChange))

#
DI <- scz_di_fc %>%
  group_by(gene) %>%
  summarize(DI = mean(DI))

#write.table(scz_di_fc,file="tables/dataS3.ASD.txt",sep="\t",col.names=T,row.names=F,quote=F)
#write.table(scz_di_fc,file="tables/dataS3.SCZ.txt",sep="\t",col.names=T,row.names=F,quote=F)
#write.table(scz_di_fc,file="tables/dataS3.BPD.txt",sep="\t",col.names=T,row.names=F,quote=F)


col_fun = colorRamp2(c(min(scale(mat)), 0, max(scale(mat))), c("#ee0000ff", "white", "#3b4992ff"))

row_annotations <- rowAnnotation(
  BPD = anno_barplot(lfc_data$avg_log2FoldChange, gp = gpar(fill = "grey")),
  DI = anno_barplot(DI$DI, gp = gpar(fill = "grey")),
  gap = unit(2, "mm"),
  width = unit(2.5, "cm"),
  annotation_name_gp = gpar(fontsize = 7)
)


h=Heatmap(
  scale(mat),
  name = "Drug DE",
  col=col_fun,
  border=TRUE,
  right_annotation = row_annotations,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(
    legend_height = unit(2, "cm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 7),
    legend_width = unit(2, "cm"),
    legend_direction = "horizontal"
  )
)


#pdf("2F.pdf",height=2,width=2.75)
#draw(h,heatmap_legend_side = "bottom")
#dev.off()

#pdf("2E.pdf",height=2,width=2.75)
#draw(h,heatmap_legend_side = "bottom")
#dev.off()

#pdf("S4E.pdf",height=2,width=2.75)
#draw(h,heatmap_legend_side = "bottom")
#dev.off()
