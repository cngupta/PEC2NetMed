#read GRNs, find coreg matrix, find modules using WGCNA

rm(list=ls())
set.seed(123)
source('load_libraries.R')
source('aux_functions.R')

detect_delimiter = function(filepath) {
  first_line = readLines(filepath, n = 1)

  if (grepl("\t", first_line)) {
    return("\t")
  } else if (grepl(",", first_line)) {
    return(",")
  } else {
    stop("Unknown delimiter")
  }
}

read_delimited_file = function(filepath) {
  # Detect delimiter
  sep = detect_delimiter(filepath)

  return(read.table(filepath, header = TRUE, sep = sep))
}

phenos=c("ASD","BPD","SCZ")

grn_path = "data/PEC2_NetMed_GRNs/subclass/"

all_link_counts = data.frame()


for(pheno in phenos)
{
  pheno_grn_path = paste(grn_path,pheno,sep="")
  # Loop through GRN files in the directory
  for (filename in list.files(pheno_grn_path))
  {
      print(filename)
      ptm = proc.time()
      message(paste("working on",filename,sep="::"))
      tag =  sub("^(.*?)_.*$", "\\1", filename)
      grn_data = read_delimited_file(paste(pheno_grn_path,filename,sep="/"))
      link_counts = grn_data %>%
        group_by(link) %>%
        summarise(count = n()) %>%
        mutate(cell_type = tag) %>%
        mutate(disorder = pheno)
    all_link_counts = rbind(all_link_counts, link_counts)
  }
}

celltypes=c('Lamp5', 'Pvalb', 'Sncg', 'Sst', 'Lamp5.Lhx6', 'Vip','Pax6','Chandelier',
'L2.3.IT', 'L4.IT', 'L5.IT', 'L5.ET', 'L5.6.NP','L6b','L6.IT','L6.CT','L6.IT.Car3',
"endo","astro","opc","micro","oligo","vlmc")

all_link_counts$cell_type = factor(all_link_counts$cell_type, levels = celltypes)

p=ggplot(all_link_counts, aes(x = cell_type, y = count, fill = link)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
    facet_wrap(~disorder) + # Rotate the "x-axis" labels after flipping
  labs(x = "Cell Type", y = "# Links", fill = "Edge") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_fill_manual(values = c("proximal" = "#DF8F44FF", "distal" = "#00A1D5FF")) +
  theme_bw() +
  theme(
  axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),  # Rotate y-axis labels before the flip
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8), # Adjusts the title size in the legend
    legend.key.size = unit(0.5, "cm")
  )

ggsave(p,file="S2B.pdf",width=4,height=4)
