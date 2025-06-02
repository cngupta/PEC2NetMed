rm(list=ls())
source('load_libraries.R')
source('aux_functions.R')

# Read in databases and preprocess data
db_complete = read.csv("data/drugbank.tsv", sep = "\t") # download drugbank
dgib=read.csv("data/interactions.tsv",sep="\t") # from DGIB
dgib$drug_claim_name=tolower(dgib$drug_claim_name)

# find compounds that are BBB+
b3db = read.csv("data/B3DB_classification.tsv", sep = "\t")
bbb=b3db[b3db[,7] == "BBB+",]
bbb$compound_name=tolower(bbb$compound_name)

# Intersect based on common drugs bet db and bbb+
common_drugs = intersect(bbb$IUPAC_name, db_complete$IUPAC)
db_n_b3db_dbids = db_complete$drugbank_id[db_complete$IUPAC %in% common_drugs]

# Intersect based on common drugs bet dgib and bbb+
common_drugs = intersect(bbb$compound_name, dgib$drug_claim_name)
dgib_n_b3db = dgib %>% filter(drug_claim_name %in% common_drugs) %>%
  select(drug_claim_name,gene_name)
colnames(dgib_n_b3db)=c("Drug","Target")


# Read and preprocess drug-target interactions, select drugs and targets with BBB+
drug_target = read.csv("data/pharmacologically_active.csv", sep = ",")
drug_target=drug_target[,c("Gene.Name","Drug.IDs")]
colnames(drug_target)=c("Target","drugbank_id")
drug_target=drug_target[,c("drugbank_id","Target")]
drug_target=separate_rows(drug_target, drugbank_id, Target, sep=";",convert = TRUE)
drug_target$drugbank_id=gsub(" ","",drug_target$drugbank_id)
drug_target=drug_target[drug_target$drugbank_id %in% db_n_b3db_dbids,]

db_names=db_complete %>% select(drugbank_id,name)
db_names$name=tolower(db_names$name)
drug_target=left_join(drug_target,db_names)
drug_target=drug_target %>% select(name,Target)
colnames(drug_target)=c("Drug","Target")

# merge two sources
full_drug_target=rbind(drug_target,dgib_n_b3db)

# Function to calculate enrichment
calculate_enrichment = function(drug_distal, total_distal_links, total_drug_targets, total_genes_in_grn) {
  observed = length(unique(drug_distal$target))
  expected = (total_drug_targets / total_genes_in_grn) * total_distal_links
  enrichment_score = observed / expected
  return(enrichment_score)
}


# Function to detect delimiter
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

# Function to read file with detected delimiter
read_delimited_file = function(filepath) {
  sep = detect_delimiter(filepath)

  return(read.table(filepath, header = TRUE, sep = sep))
}

phenos=c("ASD","BPD","SCZ")

grn_path = "data/PEC2_NetMed_GRNs/subclass/"

drug_ep = data.frame(pheno=character(), celltype=character(), ndrugs=numeric(), total_distal=numeric(), p_value=numeric(), enrichment=numeric())

for(pheno in phenos)
{
  pheno_grn_path = paste(grn_path,pheno,sep="")
  for (filename in list.files(pheno_grn_path))
  {
      print(filename)
      ptm = proc.time()
      message(paste("working on",filename,sep="::"))
      tag =  sub("^(.*?)_.*$", "\\1", filename)

      network_data = read_delimited_file(paste(pheno_grn_path,filename,sep="/"))
      drug_target_filtered=full_drug_target %>% filter(Target %in% network_data$target )
      drug_distal=network_data %>% dplyr::filter(link %in% "distal") %>%
        dplyr::filter(target %in% drug_target$Target)

      ndr=nrow(drug_target %>% filter(Target %in% drug_distal$target) %>% dplyr::select(Drug) %>% unique())

      # Enrichment test
      total_genes_in_grn = length(unique(network_data$target))
      total_drug_targets_in_grn = length(unique(drug_target_filtered$Target))
      total_distal_links = length(unique(network_data %>% filter(link == "distal") %>% pull(target)))

      # Normalization: Calculate enrichment score
      enrichment_score = calculate_enrichment(drug_distal, total_distal_links, total_drug_targets_in_grn, total_genes_in_grn)

      # Constructing the contingency table for Fisher's Exact Test
      n_distal_drug_targets = length(unique(drug_distal$target))
      n_distal_non_drug_targets = total_distal_links - n_distal_drug_targets
      n_non_distal_drug_targets = total_drug_targets_in_grn - n_distal_drug_targets
      n_non_distal_non_drug_targets = total_genes_in_grn - total_drug_targets_in_grn - n_distal_non_drug_targets

      contingency_table = matrix(c(n_distal_drug_targets, n_distal_non_drug_targets,
                                   n_non_distal_drug_targets, n_non_distal_non_drug_targets),
                                 nrow=2, byrow=TRUE)

      fisher_result = fisher.test(contingency_table)
      p_value = fisher_result$p.value

      df=data.frame(pheno=pheno,celltype=tag,ndrugs=ndr,total_distal=nrow(network_data),p_value=p_value,enrichment=enrichment_score)
      drug_ep=rbind(drug_ep,df)

  }
}

celltypes=c('Lamp5', 'Pvalb', 'Sncg', 'Sst', 'Lamp5.Lhx6', 'Vip','Pax6','Chandelier',
'L2.3.IT', 'L4.IT', 'L5.IT', 'L5.ET', 'L5.6.NP','L6b','L6.IT','L6.CT','L6.IT.Car3',
"endo","astro","opc","micro","oligo","vlmc")

drug_ep$celltype = factor(drug_ep$celltype, levels = celltypes)


p1=ggplot(drug_ep, aes(x = -log10(p_value) , y =celltype, fill = pheno)) +
  geom_bar(stat = "identity", color = "white", alpha = 0.8, position = position_dodge()) +
  labs(y = "", x = "-log10(p-value)") +
  coord_flip() +
  scale_fill_manual(values =  c("SCZ" = "#0072b5ff", "BPD" = "#e18727ff", "ASD" = "#ee4c97ff")) +
  theme_bw() +
  theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 8, family = "Helvetica"),
    legend.position = c(0.85, 0.75),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.box.spacing = unit(0.1, "cm"),
    legend.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.15, "cm"),
    legend.title.align = 0.15
  )+ coord_flip()

ggsave(p1, file="S2D.pdf",width=4,height=2,units=c("in"))
