rm(list=ls())
set.seed(123)
source('load_libraries.R')
source('aux_functions.R')
library(dplyr)
library(ggplot2)
library(ggalluvial)


#Read and filter GO data
data=GSA.read.gmt('data/Human_GO_bp_with_GO_iea_symbol.gmt')
genesets=data$genesets
names(genesets)=data$geneset.descriptions
geneset_sizes = sapply(genesets, length)
genesets = genesets[geneset_sizes < 300 & geneset_sizes > 30]

# read module data
modules=read.table("tables/WGCNA_modules_merged.csv",header=T)

# Count the number of genes within each moduleID
gene_count_per_moduleID = modules %>%
  group_by(moduleID) %>%
  summarise(gene_count = n(), .groups = 'drop')
gene_count_per_moduleID = gene_count_per_moduleID %>%
  separate(moduleID, into = c("module","Cohort", "Phenotype", "Celltype"), sep = "::")
average_module_size = gene_count_per_moduleID %>%
  group_by(Phenotype, Celltype, module) %>%
  summarise(average_module_size = mean(gene_count), .groups = 'drop')


modules = modules %>%
  separate(celltype, into = c("Cohort", "Phenotype", "Celltype"), sep = "::")
unique_module_count = modules %>%
  group_by(Phenotype, Celltype) %>%
  summarise(unique_module_count = n_distinct(module), .groups = 'drop')

# combine
combined_df = average_module_size %>%
  inner_join(unique_module_count, by = c("Phenotype", "Celltype")) %>%
  pivot_longer(cols = c("average_module_size", "unique_module_count"),
               names_to = "Metric",
               values_to = "Value")

p=ggplot(combined_df, aes(x = Celltype, y = Value, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Metric ~ ., scales = "free_y") +
  labs(title = "Avg. module size & total modules",
       x = "",
       y = "Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")+
  theme_minimal(base_size = 10) +
  labs(title = "", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = unit(0.25, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.title.position = "top",
    legend.position = "top",
    legend.text = element_text(size = 8),  # Adjust legend text size
    legend.title = element_text(size = 8)
  )

ggsave(p,file="S6A.pdf",width=4, height=6,units=c("in"))

modules=modules %>% dplyr::filter(module != "M_0")
moduleID=unique(modules$moduleID)

enrich.tbl=data.frame("label"=NULL,pval=NULL,fdr=NULL,signature=NULL,geneset=NULL,
overlap=NULL,background=NULL, "hits"=NULL,"moduleID"=NULL)

i=1
total_mods=length(unique(moduleID))
for (module in moduleID)
{
  i=i+1
  print(paste("module progress",(i/total_mods)*100,sep="::"))
  genes=modules[moduleID %in% module,"gene_name"]
  if (length(genes) > 10 & length(genes) < 100)
  {
    hyp_obj = hypeR(genes, genesets, fdr=0.1,background=length(unique(modules$gene_name)))
    hyp_df =  hyp_obj$data
    if(nrow(hyp_df) > 0)
    {
      print(nrow(hyp_df))
      hyp_df$moduleID=module
      enrich.tbl=rbind(enrich.tbl,hyp_df)
    }
  }
  print(nrow(enrich.tbl))
}

df=enrich.tbl %>% separate(moduleID, into = c("module","Cohort", "Phenotype", "Celltype"), sep = "::")
write.table(df,file="tables/modules_GO_enrichment.txt",sep="\t",quote=F,col.names=TRUE,row.names=F)


term_count_per_dis = df %>%
  group_by(Phenotype) %>%
  summarise(label = n(), .groups = 'drop')

p=ggplot(term_count_per_dis, aes(x = Phenotype, y = label, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "",
       y = "# of GO BP \n terms recovered") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")+
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = unit(0.25, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.title.position = "top",
    legend.position = "none",
    legend.text = element_text(size = 8),  # Adjust legend text size
    legend.title = element_text(size = 8)
  )
#ggsave(p,file="S6B.pdf",width=1.5, height=1.5,units=c("in"))


term_count_per_celltype = df %>%
  group_by(Celltype) %>%
  summarise(label = n(), .groups = 'drop')

p=ggplot(term_count_per_celltype, aes(x = Celltype, y = label, fill = Celltype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "",
       y = "# of GO BP \n terms recovered") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")+
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = unit(0.25, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.title.position = "top",
    legend.position = "none",
    legend.text = element_text(size = 8),  # Adjust legend text size
    legend.title = element_text(size = 8)
  )
#ggsave(p,file="S6C.pdf",width=2.5, height=1.5,units=c("in"))

df=read.table("tables/modules_GO_enrichment.txt",header=T,sep="\t")
df=df %>% filter(fdr <= 0.01) %>% dplyr::select(label,Celltype,Phenotype)

sankey_data = df %>%
  group_by(label, Celltype,Phenotype) %>%
  summarize(count = n(), .groups = 'drop')

# Create the Sankey plot
p=ggplot(sankey_data,
       aes(axis1 = Celltype, axis2 = label, ,axis3=Phenotype, y = count)) +
  geom_alluvium(aes(fill = Celltype), width = 0.1) +
  scale_fill_brewer(palette = "Set1")+
  geom_stratum(width = 0.2, fill = "grey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(str_wrap(stratum, width = 75))), size = 3, angle = 0) +  # Wrap text labels
  scale_x_discrete(limits = c("Celltype", "GO","Disorder"), expand = c(0.15, 0.15)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  labs(x = "", y = "", title = "Module GO enrichment")+
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1))


pdf("S6D.pdf",width=6, height=8)
p
dev.off()
