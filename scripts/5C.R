# run Fig5.R before this
# rm(list=ls())


source('load_libraries.R')
source('functions_for_network_analysis.R')
library(GSEABase)

library(dplyr)
library(purrr)
library(fgsea)
library(tibble)
library(rlang)
library(dplyr)
library(stringr)
library(ggnewscale)


get_fc_and_ks_by_celltype = function(degdf, q_geneset) {
  results_list = list()

  cell_types = unique(degdf$cell_type)

  for (set_name in names(q_geneset)) {
    gene_ids = geneIds(q_geneset[[set_name]])

    for (ct in cell_types) {
      df_ct = degdf %>%
        filter(cell_type == ct, !is.na(log2FoldChange))

      if (nrow(df_ct) < 5) next

      df_ct = df_ct %>%
        arrange(desc(log2FoldChange)) %>%
        mutate(in_gene_set = gene %in% gene_ids)

      mean_fc = df_ct %>%
        filter(in_gene_set) %>%
        summarize(mean_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)) %>%
        pull(mean_log2FoldChange)

      if (grepl(" up$", set_name)) {
        direction = ifelse(mean_fc > 0, "Mimicker", "Reverser")
      } else if (grepl(" down$", set_name)) {
        direction = ifelse(mean_fc < 0, "Mimicker", "Reverser")
      } else {
        direction = NA
      }

      if (sum(df_ct$in_gene_set) > 1) {
        ks = ks.test(
          x = df_ct$log2FoldChange[df_ct$in_gene_set],
          y = df_ct$log2FoldChange[!df_ct$in_gene_set],
          alternative = "two.sided"
        )

        results_list[[paste(set_name, ct, sep = "::")]] = tibble(
          set_name = set_name,
          cell_type = ct,
          mean_log2FoldChange = mean_fc,
          ks_statistic = ks$statistic,
          p_value = ks$p.value,
          direction = direction
        )
      }
    }
  }

  return(bind_rows(results_list))
}

#LINCS
gene_sets = getGmt("data/l1000_cp.gmt")
metadata = read.delim("data/LINCS_small_molecules.tsv", header = TRUE, sep = "\t", fill = TRUE, quote = "")

#DEGs
degasd=read.csv("ASD_DEGcombined.csv")
degscz=read.csv("Schizophrenia_DEGcombined.csv")
degbpd=read.csv("Bipolar_DEGcombined.csv")

ASD=read.csv("tables/proximity_results.signif.subclass.ASD.csv",sep= "\t")
ASD$dis="ASD"
ASD$cell=gsub("_merged.csv","",ASD$cell)

SCZ=read.csv("tables/proximity_results.signif.subclass.SCZ.csv",sep="\t")
SCZ$dis="SCZ"
SCZ$cell=gsub("_merged.csv","",SCZ$cell)

BPD=read.csv("tables/proximity_results.signif.subclass.BPD.csv",sep="\t")
BPD$dis="BPD"
BPD$cell=gsub("_merged.csv","",BPD$cell)


proximity_results.signif=rbind(ASD,rbind(SCZ,BPD))
proximity_results.signif=proximity_results.signif %>% filter(proximity != Inf)

df=proximity_results.signif %>% filter(z<=-3, p<=0.001)
db_complete = read.csv("drugbank.tsv", sep = "\t")
db_complete$name=tolower(db_complete$name)
combined_df = merge(df, db_complete, by.x = "drug", by.y = "name", all = F)
colnames(combined_df)[6]="Disorder"
proximity_results.signif.merged=combined_df
proximity_results.signif.merged$cell_salt=paste(proximity_results.signif.merged$cell,proximity_results.signif.merged$drug,sep="::")

matched_drugs = combined_df %>% left_join(metadata, by = c("inchikey" = "inchi_key"),relationship = "many-to-many")

unique_drug_names = unique(matched_drugs$drug)

subset_gene_sets = lapply(unique_drug_names, function(drug) {
  grep(drug, names(gene_sets), value = TRUE, ignore.case = TRUE)
})

subset_gene_sets_names = unique(unlist(subset_gene_sets))

subset_gene_sets_final = gene_sets[subset_gene_sets_names]
gene_set_names = names(subset_gene_sets_final)

#CENTRAL nervous system cell line NEU
neu_gene_sets_names = grep("NEU", gene_set_names, value = TRUE, ignore.case = TRUE)
neu_gene_sets = subset_gene_sets_final[neu_gene_sets_names]

#if using fgsea
neu_gene_sets_list = geneIds(neu_gene_sets)

ASD_effects=get_fc_and_ks_by_celltype(degasd,neu_gene_sets)
SCZ_effects=get_fc_and_ks_by_celltype(degscz,neu_gene_sets)
BPD_effects=get_fc_and_ks_by_celltype(degbpd,neu_gene_sets)

ASD_effects$Disorder="ASD"
SCZ_effects$Disorder="SCZ"
BPD_effects$Disorder="BPD"

combined_df=rbind(ASD_effects,rbind(SCZ_effects,BPD_effects))

#
combined_df = combined_df %>%
  mutate(salt_name = gsub("^[^_]*_[^_]*_[^_]*_[^_]*_", "", set_name))
combined_df = combined_df %>%
  mutate(salt_name = gsub("_.*", "", salt_name))

combined_df$cell_salt=paste(combined_df$cell_type,combined_df$salt_name,sep="::")

report_df = combined_df %>% mutate(padj = p.adjust(p_value, method = "BH"))
combined_df=report_df

# Plot top drugs data using ggplot2 with the extracted salt names and borders
df=proximity_results.signif %>% filter(z< -8, p<=0.001)
#df=proximity_results.signif
top_drugs=unique(df$drug)
plotdf=combined_df  %>% filter(salt_name %in% top_drugs) %>% as.data.frame()
plotdf=plotdf %>% dplyr::filter(!grepl("_6H_", set_name))
plotdf=plotdf %>% dplyr::filter(!grepl("NMH002_NEU_24H_O06_riluzole", set_name))
plotdf=plotdf %>% dplyr::filter(!grepl("CPC016_NEU_24H_F03_riluzole", set_name))

top_drugs_with_effects=unique(plotdf$salt_name)
top_drugs_with_effects_metadata=db_complete %>% filter(name %in% top_drugs) %>% dplyr::select(name, description)

# a drug is selected only if its a reverser or a mimiker in both up and down gene sets.
plotdf_consistent = plotdf %>%
  mutate(
    set_type = case_when(
      str_detect(set_name, " up$") ~ "up",
      str_detect(set_name, " down$") ~ "down",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(cell_type, Disorder, salt_name) %>%
  filter(
    n_distinct(set_type) == 2,
    dplyr::n_distinct(.data$direction) == 1
  ) %>%
  ungroup()

plotdf_consistent = plotdf_consistent %>%
  mutate(
    logp = -log10(padj)
  )

p1=ggplot() +
  # Reversers first
  geom_tile(
    data = plotdf_consistent %>% filter(direction == "Reverser"),
    aes(x = 1, y = Disorder, fill = logp),
    color = "black"
  ) +
  scale_fill_gradient(
    name = "Reverser -log10(p)",
    low = "#CFE2F3", high = "#0B5394", limits = c(0, max(plotdf_consistent$logp))
  ) +

  new_scale_fill() +

  # Mimickers next
  geom_tile(
    data = plotdf_consistent %>% filter(direction == "Mimicker"),
    aes(x = 1, y = Disorder, fill = logp),
    color = "black"
  ) +
  scale_fill_gradient(
    name = "Mimicker -log10(p)",
    low = "#EA9999", high = "#990000", limits = c(0, max(plotdf_consistent$logp))
  ) +

  facet_grid(salt_name ~ cell_type, switch = "both") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.2, "lines"),
    panel.border = element_rect(fill = NA),
    axis.text = element_text(size = 7),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    strip.text.y.left = element_text(angle = 0),
    strip.text.x = element_text(angle = 90, vjust = 0.1, hjust = 0.1),
    legend.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )+
  theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 6),
  legend.key.height = unit(0.25, "cm"),
  legend.key.width = unit(0.5, "cm"),
  legend.spacing.x = unit(0.2, "cm")
)


ggsave(p1,file="5C.pdf", height=4,width=5.5)


p2= ggplot(proximity_results.signif %>% filter(drug == 'asenapine'), aes(x = cell, y = z, fill = dis)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(proximity, 2)), vjust = -0.5) +
  labs(title = "Asenapine",
       x = "",
       y = "Proximity (Z)",
       fill= "Disease") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p2,file="S11C.pdf", height=3.5,width=3.5)


p4= ggplot(proximity_results.signif %>% filter(drug == 'iloperidone'), aes(x = cell, y = z, fill = dis)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(proximity, 2)), vjust = -0.5) +
  labs(title = "Iloperidone",
       x = "",
       y = "Proximity (Z)",
       fill= "Disease") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p4,file="S11D.pdf", height=3.5,width=3.5)
