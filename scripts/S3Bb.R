rm(list=ls())
source('load_libraries.R')
source('aux_functions.R')

compare_bins_overlap = function(data1, data2, bins = 10)
{
  # Add bin column to each dataset based on the importance
  data1 = data1 %>%
   mutate(importance = as.numeric(importance)) %>%
   arrange(desc(importance)) %>%
  mutate(bin = ntile(-(importance), bins))

  data2 = data2 %>%
   mutate(importance = as.numeric(importance)) %>%
   arrange(desc(importance)) %>%
  mutate(bin = ntile(-(importance), bins))

  overlap_results = list()

  for (i in 1:bins)
  {
    bin1_links = data1 %>%
      filter(bin == i) %>%
      pull(edge)

    bin2_links = data2 %>%
      filter(bin == i) %>%
      pull(edge)

    overlap = intersect(bin1_links, bin2_links)

    overlap_results[[i]] = list(
      bin = i,
      overlap_count = length(overlap),
      bin1_size = length(bin1_links),
      bin2_size = length(bin2_links),
      overlap_percentage = (length(overlap) / min(length(bin1_links), length(bin2_links))) * 100
    )
  }
  overlap_df = do.call(rbind, lapply(overlap_results, as.data.frame))

  return(overlap_df)
}


cohorts=c("CMC","Urban_DLPFC")
celltypes=c("endo","astro","opc","micro","oligo","vlmc","inhibitory","excitatory")

grn_path = "data/PEC2_NetMed_GRNs/major_class"
pheno="SCZ"

all_overlap_results = data.frame()


for(cell in celltypes)
{
  # Directory path containing the GRN files
  file=paste(cell,"_merged.csv",sep="")
  pheno_grn_path_1 = paste(grn_path,"CMC",sep="/")
  pheno_grn_path_1 = paste(pheno_grn_path_1,pheno,sep="/")
  pheno_grn_path_1 = paste(pheno_grn_path_1,file,sep="/")

  pheno_grn_path_2 = paste(grn_path,"Urban_DLPFC",sep="/")
  pheno_grn_path_2 = paste(pheno_grn_path_2,pheno,sep="/")
  pheno_grn_path_2 = paste(pheno_grn_path_2,file,sep="/")

  network_data1 = read.table(pheno_grn_path_1,header=T,sep="\t")
  network_data1$edge=paste(network_data1$TF,network_data1$target,sep="::")
  network_data1=network_data1 %>% filter(!(link == 'distal'))
  network_data2 = read.table(pheno_grn_path_2,header=T,sep="\t")
  network_data2$edge=paste(network_data2$TF,network_data2$target,sep="::")
  network_data2=network_data2 %>% filter(!(link == 'distal'))

  overlap_comparison = compare_bins_overlap(network_data1, network_data2, bins = 10)
  overlap_comparison$cell=cell

  all_overlap_results = rbind(all_overlap_results, overlap_comparison)
}

p=ggplot(all_overlap_results, aes(x = as.factor(bin), y = overlap_percentage, color = cell, group = cell)) +
  geom_line(size = 0.5) +
  geom_point(size = 1.5) +
  labs(
    title = "",
    x = "SCZ bins",
    y = "Overlap Percentage (%)",
    color = "Cell Type"
  ) +
  theme_bw() +
  theme(
  axis.title.x = element_text(size = 8),
  axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.25, "cm")
  )

ggsave(p,file="S3B.pdf",width=4,height=2.25)
