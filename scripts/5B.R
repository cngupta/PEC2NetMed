# run Fig5.R before this
# rm(list=ls())


source('load_libraries.R')
source('functions_for_network_analysis.R')
source('functions_for_drug_repurposing.R')
library(ggalluvial)

metadata = read.delim("data/LINCS_small_molecules.tsv", header = TRUE, sep = "\t", fill = TRUE, quote = "")

ASD=read.csv("tables/proximity_results.signif.subclass.ASD.csv",sep= "\t")
ASD$dis="ASD"
SCZ=read.csv("tables/proximity_results.signif.subclass.SCZ.csv",sep="\t")
SCZ$dis="SCZ"
BPD=read.csv("tables/proximity_results.signif.subclass.BPD.csv",sep="\t")
BPD$dis="BPD"

proximity_results.signif=rbind(ASD,rbind(SCZ,BPD))
proximity_results.signif$cell=gsub("_merged.csv","",proximity_results.signif$cell)
colnames(proximity_results.signif)[6]="Disorder"

celltypes=c('Lamp5', 'Pvalb', 'Sncg', 'Sst', 'Lamp5.Lhx6', 'Vip','Pax6','Chandelier',
  'L2.3.IT', 'L4.IT', 'L5.IT', 'L5.ET', 'L5.6.NP','L6b','L6.IT','L6.CT','L6.IT.Car3',
  "endo","astro","opc","micro","oligo","vlmc")

proximity_results.signif$cell = factor(proximity_results.signif$cell, levels = celltypes)


proximity_results.signif=proximity_results.signif %>%
  filter(proximity != Inf) %>%
  filter(p <= 0.001)  %>%
  filter(z <= -3)


drug_counts=proximity_results.signif %>%
   group_by(cell, Disorder) %>%
   summarise(num_drugs = n())
p1=ggplot(drug_counts, aes(x = cell, y = num_drugs)) +
   geom_bar(stat = "identity") +
   facet_wrap(~ Disorder) +
   labs(title = "Z <= -3", y = "# of molecules", x = "")+
   theme_bw()+
   theme(
   axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
     legend.position = "top",
     legend.text = element_text(size = 8),
     legend.title = element_text(size = 8),
     legend.key.size = unit(0.5, "cm")
   )

p1=p1+coord_flip()
ggsave(p1,file="S11A.pdf",height=3.5,width=3.25,units=c("in"))

p2=ggplot(drug_counts, aes(x = cell, y = num_drugs)) +
   geom_bar(stat = "identity") +
   labs(title = "Z <= -3", y = "# of molecules", x = "")+
   theme_bw()+
   theme(
   axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
     legend.position = "top",
     legend.text = element_text(size = 8),
     legend.title = element_text(size = 8),
     legend.key.size = unit(0.5, "cm")
   )

p2=p2+coord_flip()
ggsave(p2,file="S11B.pdf",height=3.5,width=2,units=c("in"))



df=proximity_results.signif %>% filter(z<=-2)
db_complete = read.csv("~/work/DrugBank2022/complete/drugbank.tsv", sep = "\t")
db_complete$name=tolower(db_complete$name)
combined_df = merge(df, db_complete, by.x = "drug", by.y = "name", all = F)
combined_df$cell=gsub("_merged.csv","",combined_df$cell)
colnames(combined_df)[6]="Disorder"
df=combined_df

plotdf = df %>%
  dplyr::filter(z < -3) %>%
  dplyr::filter(if_all(where(is.numeric), is.finite))

matched_drugs = plotdf %>% left_join(metadata, by = c("inchikey" = "inchi_key"),relationship = "many-to-many") %>%
  dplyr::filter(moa != 'NA') %>% dplyr::filter(moa != '-')

#remove long name, take only first moa
matched_drugs$moa = sub(",.*", "", matched_drugs$moa)

# for rev 1
df_for_cytoscape=matched_drugs[,c("cell","Disorder","drug","moa","z")]

# Keep top 10% drugs per disorder based on absolute z
top_drugs_df = df_for_cytoscape %>%
  group_by(Disorder, drug) %>%
  summarise(mean_abs_z = mean(abs(z)), .groups = "drop") %>%
  group_by(Disorder) %>%
  slice_max(mean_abs_z, prop = 0.1) %>%
  select(Disorder, drug)

filtered_df = df_for_cytoscape %>%
  semi_join(top_drugs_df, by = c("Disorder", "drug"))

# rev 1
# flow transparency colored by z
p2=ggplot(data = filtered_df,
       aes(axis1 = Disorder, axis2 = cell, axis3 = drug, axis4 = moa, y = 1)) +
  geom_alluvium(aes(fill = Disorder, alpha = z), width = 0.1) +
  scale_alpha_continuous(range = c(1, 0.2)) +
geom_stratum(width = 0.3, size = 0.1) +
geom_text(stat = "stratum", aes(label = after_stat(str_wrap(stratum, width = 22))),
          size = 2.5, angle = 0, hjust = 0.5) +
geom_text(stat = "stratum", aes(label = after_stat(str_wrap(stratum, width = 22))),
          size = 2.5, angle = 0, hjust = 1,
          position = position_nudge(x = -0.1),
          data = function(df) df[df$axis == 4, ]) +
scale_x_discrete(expand = c(0.1, 0.1)) +
theme_minimal() +
scale_fill_manual(values = c("SCZ" = "#12909F", "BPD" = "#F7BC25", "ASD" = "#75405A"),
guide = guide_legend(direction="vertical")) +
labs(fill = "Disorder", alpha = "Proximity \n Z-score")+
theme(
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  text = element_text(size = 8),
  legend.position = "right"
)

ggsave(p2,file="5B.pdf",height=8,width=4.25,units=c("in"))
