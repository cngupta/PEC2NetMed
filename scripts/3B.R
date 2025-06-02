# run magma before this and store all .out files to data/
rm(list=ls())
source('load_libraries.R')
source('aux_functions.R')


safri=read.table("data/UCLA_ASD_ASD_micro_safri_enrichment_results.txt",header=T)
deg1=read.table("data/gandal_asd.fc-UCLA_ASD_ASD_micro.zscore.mat",header=T)
deg2=read.table("data/41586_2022_5377_MOESM5_ESM.LFC-UCLA_ASD_ASD_micro.zscore.mat",header=T)


# Function to process MAGMA files
process_magma = function(file_path, variable_name) {
  magma = read.table(file_path, header = TRUE)
  magma = magma %>% select(VARIABLE, P)
  colnames(magma) = c("Module", variable_name)
  return(magma)
}

# outputs from magma runs 
magma_asd = process_magma("data/asd_UCLA_ASD_ASD_micro.gmt.gsa.out", "ASD_GWAS_P")
magma_scz = process_magma("data/scz_UCLA_ASD_ASD_micro.gmt.gsa.out", "SCZ_GWAS_P")
magma_bpd = process_magma("data/bpd_UCLA_ASD_ASD_micro.gmt.gsa.out", "BPD_GWAS_P")
magma_epilepsy = process_magma("data/epilepsy_UCLA_ASD_ASD_micro.gmt.gsa.out", "Epilepsy_GWAS_P")


#magma=magma %>% select(VARIABLE,P)
#colnames(magma)=c("Module","GWAS_P")

deg1=deg1 %>% select(gs.desc,ASD)
colnames(deg1)=c("Module","ASD_FC")

deg2$gs.id=NULL
deg2$gs.id=NULL
deg2$gs.ngenes=NULL
colnames(deg2)[1]="Module"

library(ggnewscale)

safri=safri %>% select(module,p_value)
colnames(safri)=c("Module","HyperP")

#df=merge(deg1,deg2)
#df=merge(df,merge(safri,magma))
df = merge(deg1, deg2)
df = merge(df, safri)
df = merge(df, magma_asd)
df = merge(df, magma_scz)
df = merge(df, magma_bpd)
df = merge(df, magma_epilepsy)


rownames(df)=df$Module
rownames(df)=gsub("_"," ",rownames(df))

plotdf=df %>% select(Module,WholeCortex_ASD_logFC,SexF_logFC,Age_logFC,HyperP,ASD_GWAS_P, SCZ_GWAS_P, BPD_GWAS_P, Epilepsy_GWAS_P)
plotdf$Module=NULL

plotdf=t(plotdf)
colnames(plotdf)=gsub("_"," ",colnames(plotdf))
rownames(plotdf)=gsub("_logFC"," ",rownames(plotdf))
rownames(plotdf)=gsub("_ASD"," ",rownames(plotdf))


# Create labels for x-axis
module_labels = colnames(plotdf)
labels = rep("", length(module_labels))


# revision 1
labels[module_labels == "M 38"] = "M 38"
labels[module_labels == "M 6"] = "M 6"
labels[module_labels == "M 52"] = "M 52"


col_fun = colorRamp2(c(-2, 0, 2), c("#ee0000ff", "white", "#3b4992ff"))

#c("#e18727ff","#0072b5ff","#ee4c97ff")
column_ha = HeatmapAnnotation(
  SAFRI = anno_lines(-log10(plotdf["HyperP", ]), gp = gpar(col = "grey"), add_points = TRUE,pt_gp = gpar(col = "#ee4c97ff")),
  ASD_GWAS = anno_lines(-log10(plotdf["ASD_GWAS_P", ]), gp = gpar(col = "#ee4c97ff"), add_points = TRUE,pt_gp = gpar(col = "#ee4c97ff")),
  SCZ_GWAS = anno_lines(-log10(plotdf["SCZ_GWAS_P", ]), gp = gpar(col = "#0072b5ff"), add_points = TRUE,pt_gp = gpar(col = "#0072b5ff")),
  BPD_GWAS = anno_lines(-log10(plotdf["BPD_GWAS_P", ]), gp = gpar(col = "#e18727ff"), add_points = TRUE,pt_gp = gpar(col =  "#e18727ff")),
  Epilepsy_GWAS = anno_lines(-log10(plotdf["Epilepsy_GWAS_P", ]), gp = gpar(col = "grey"), add_points = TRUE,pt_gp = gpar(col = "grey")),
  border = c(SAFRI = FALSE),
  gap = unit(1, "mm"),
  annotation_name_gp = gpar(fontsize = 7)
)


h=Heatmap(
  plotdf[1:3, ],
  col = col_fun,
  top_annotation = column_ha,
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_labels =  labels,  #module_labels, # Add labels for specific modules
  border = TRUE,
  height = unit(0.5, "in"),
  heatmap_legend_param = list(
    title = "Module Expression\n(Z score)",
    legend_height = unit(2, "cm"),
    legend_height = unit(4, "cm"),
    legend_direction = "horizontal",
    legend_position = "top"
  )
)



#revision 1
pdf("3B.pdf",height=4,width=12)
draw(h)
dev.off()
