rm(list=ls())
set.seed(123)
source('load_libraries.R')
source('aux_functions.R')

detect_delimiter = function(filepath) {
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

read_delimited_file = function(filepath) {
  # Detect delimiter
  sep = detect_delimiter(filepath)
  return(read.table(filepath, header = TRUE, sep = sep))
}

phenos=c("ASD","BPD","SCZ")

grn_path = "data/PEC2_NetMed_GRNs/subclass/"

all_link_counts = data.frame()

node.tbl=data.frame(disorder=character(),celltype=character(),nTFs=numeric(),ngenes=numeric())

for(pheno in phenos)
{
  pheno_grn_path = paste(grn_path,pheno,sep="")
  for (filename in list.files(pheno_grn_path))
  {
      print(filename)
      ptm = proc.time()
      message(paste("working on",filename,sep="::"))
      tag =  sub("^(.*?)_.*$", "\\1", filename)
      grn_data = read_delimited_file(paste(pheno_grn_path,filename,sep="/"))
      df=data.frame(disorder=pheno,celltype=tag,nTFs=length(unique(grn_data$TF)),ngenes=length(grn_data$target))
      node.tbl=rbind(node.tbl,df)
  }
}

celltypes=c('Lamp5', 'Pvalb', 'Sncg', 'Sst', 'Lamp5.Lhx6', 'Vip','Pax6','Chandelier',
'L2.3.IT', 'L4.IT', 'L5.IT', 'L5.ET', 'L5.6.NP','L6b','L6.IT','L6.CT','L6.IT.Car3',
"endo","astro","opc","micro","oligo","vlmc")

node.tbl$celltype = factor(node.tbl$celltype, levels = celltypes)


plot_data = node.tbl %>%
  gather(key = "Metric", value = "Count", nTFs, ngenes)

# Calculate the average ngenes count for each disorder
average_ngenes_per_disorder = plot_data %>%
  filter(Metric == "ngenes") %>%
  group_by(disorder) %>%
  summarize(avg_ngenes = mean(Count))

p=ggplot(plot_data, aes(x = celltype, y = Count, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~ disorder, ncol = 1) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.box.spacing = unit(0.15, "cm")
  ) +
  labs(y = "Count") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +  # Set y-axis breaks to 3
  scale_fill_manual(
    values = c("nTFs" = "#1f77b4", "ngenes" = "#ff7f0e"),
    labels = c("# TFs", "# genes")
  ) +
  geom_hline(data = average_ngenes_per_disorder, aes(yintercept = avg_ngenes),
             linetype = "dashed", color = "red", size = 0.8) +
  geom_text(
    data = average_ngenes_per_disorder,
    aes(x = 1, y = avg_ngenes, label = paste0("Avg: ", round(avg_ngenes, 1))),
    inherit.aes = FALSE,
    color = "black", size = 3,
    hjust = +0.1, vjust = -0.5
  )


ggsave(p,file="S2C.pdf",width=4,height=4)
