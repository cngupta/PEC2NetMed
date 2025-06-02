rm(list=ls())
source('load_libraries.R')
source('scripts/aux_functions.R')

# Read autism gene data; repeat for other two
autism_genes = read.csv("data/SFARI-Gene_genes_01-23-2023release_03-02-2023export.csv", sep = ",") %>%
  filter(gene.score == 1 | syndromic == 1) %>%
  .$gene.symbol %>%
  unique()

# read disease genes
# created from disgenenet + a few additional genes from PEC2 cornerstone working group 
dis=read.table("disease_genes_disgennet.txt")
bp=dis[dis$V2=="BP",][,1]
scz=dis[dis$V2=="scz",][,1]


get_known_tf_plot=function(matrix,dis_genes,disease_name)
{
  matrix=matrix[,colnames(matrix) %in% dis_genes]
  n_cols = ncol(matrix)
  title_text = bquote(Delta ~ "influence of known" ~ .(disease_name) ~ "TFs " ~ n == .(n_cols))
  df=data.frame(celltype=rownames(matrix),mean_delta_influence=rowSums(matrix))
  df$pos_neg = ifelse(df$mean_delta_influence >= 0, "Positive", "Negative")
  p= ggplot(df, aes(x = celltype, y = mean_delta_influence, fill = pos_neg)) +
   geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Negative" = "#ee0000ff", "Positive" = "#3b4992ff")) +
    labs(title = "", y = title_text, x = "") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill=NA, size=0.25),
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8, family = "Helvetica"),
          axis.title = element_text(size = 7, family = "Helvetica"),
          plot.title = element_text(size = 7, family = "Helvetica", hjust = 0.5)) +
          coord_flip()

  return(list(plot=p,tfs=colnames(matrix)))
}

# tables created in 2D.R
mat=read.csv("tables/ASD_delta_influence_gt3.csv", row.names=1)
p.asd=get_known_tf_plot(mat,autism_genes,"ASD")

mat=read.csv("tables/CMC_SCZ_delta_influence_gt3.csv", row.names=1)
p.scz=get_known_tf_plot(mat,scz, "SCZ")

mat=read.csv("tables/Urban_DLPFC_BPD_delta_influence_gt3.csv", row.names=1)
p.bp=get_known_tf_plot(mat,bp, "BPD")

ggsave(p.asd$plot,file="S4Da.png",device = "png", width = 2.5, height = 1.75, units = "in")
#ggsave(p.scz$plot,file="S4Db.png",device = "png", width = 2.5, height = 1.75, units = "in")
#ggsave(p.bp$plot,file="S4c.png",device = "png", width = 2.5, height = 1.75, units = "in")
