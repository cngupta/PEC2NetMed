
rm(list=ls())
set.seed(123)
source('load_libraries.R')
source('aux_functions.R')
library(biomaRt)


# Function to detect delimiter
detect_delimiter = function(filepath) {
  # Read the first line of the file
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

# Function to read file with detected delimiter
read_delimited_file = function(filepath) {
  sep = detect_delimiter(filepath)

  return(read.table(filepath, header = TRUE, sep = sep))
}

# read gold standard (.txt files in chip seq from gtrd database)
dir_path = "data/gtrddb"  # path to GRTD files.
files = list.files(dir_path, pattern = "\\.txt$", full.names = TRUE)
combined = do.call(rbind, lapply(files, function(f) {
  df = read.table(f, header = TRUE, sep = "\t")  # adjust sep if needed
  if(nrow(df>2))
  {
    df$TF = tools::file_path_sans_ext(basename(f))  # add file name (without .txt) as 'TF'
  }
  return(df)
}))
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
uniprot_ids = unique(combined$TF)
#Query for mapping UniProt -> Gene Symbol
mapping = getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = uniprot_ids,
  mart = mart
)
combined = merge(combined, mapping, by.x = "TF", by.y = "uniprotswissprot", all.x = TRUE)
combined=combined[,c("hgnc_symbol","Gene.symbol")]
gold_standard = combined[!is.na(combined$hgnc_symbol), ]
colnames(gold_standard)=c("TF","target")
gold_standard=unique(gold_standard)
#####

phenos=c("ASD","BPD","SCZ")
grn_path = "data/PEC2_NetMed_GRNs/subclass/"
all_link_counts = data.frame()
f1.tbl=data.frame(disorder=character(),celltype=character(),f1=numeric(),precison=numeric(),recall=numeric())

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

      ref=gold_standard[gold_standard$TF %in% grn_data$TF,]
      pred=grn_data[grn_data$TF %in% ref$TF,]

      pred_edges = paste(pred$TF, pred$target, sep = "->")
      ref_edges  = paste(ref$TF, ref$target, sep = "->")
      tp = sum(pred_edges %in% ref_edges)
      fp = sum(!(pred_edges %in% ref_edges))
      fn = sum(!(ref_edges %in% pred_edges))
      precision = tp / (tp + fp)
      recall    = tp / (tp + fn)
      f1        = if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
      df=data.frame(disorder=pheno,celltype=filename,f1=f1,precison=precision,recall=recall)
      f1.tbl=rbind(f1.tbl,df)
  }
}
f1.tbl$celltype=gsub("_merged.csv","",f1.tbl$celltype)
library(ggplot2)
library(dplyr)

library(ggplot2)
library(dplyr)

# Calculate median precision
median_precision = median(f1.tbl$precison, na.rm = TRUE)

p=ggplot(f1.tbl, aes(x = celltype, y = precison, fill = disorder)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = median_precision, linetype = "dashed", color = "red", size = 0.8) +
  annotate("text", x = 1, y = median_precision,
           label = sprintf("Median = %.3f", median_precision),
           vjust = -1, hjust = 0, color = "black", size = 3.5) +
  labs(x = "Cell Type", y = "Precision in ChIP TF targets", title = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

ggsave(p,file="S3B.pdf", height=3.5,width=3.5)
