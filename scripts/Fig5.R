# Clear workspace and load necessary libraries
rm(list = ls())

start_time = Sys.time()

source('load_libraries.R') # Assume this loads required libraries like dplyr, tidyr, igraph, etc.
source('functions_for_network_analysis.R')
source('functions_for_drug_repurposing.R')

InSub=c('Lamp5', 'Pvalb', 'Sncg', 'Sst','Sst.Chodl', 'Lamp5.Lhx6', 'Vip','Pax6','Chandelier')
ExSub=c('L2.3.IT', 'L4.IT', 'L5.IT', 'L5.ET', 'L5.6.NP','L6b','L6.IT','L6.CT','L6.IT.Car3')


# Read in databases and preprocess data
db_complete = read.csv("data/drugbank.tsv", sep = "\t")
dgib=read.csv("data/interactions.tsv",sep="\t")
dgib=dgib %>% filter(approved == TRUE)
dgib=dgib %>% filter(anti_neoplastic == FALSE)
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
dgib_n_b3db = dgib %>% dplyr::filter(drug_claim_name %in% common_drugs) %>%
  dplyr::select(drug_claim_name,gene_name)
colnames(dgib_n_b3db)=c("Drug","Target")


# Read and preprocess drug-target interactions, select drugs and targets with BBB+
drug_target = read.csv("data/pharmacologically_active.csv", sep = ",")
drug_target=drug_target[,c("Gene.Name","Drug.IDs")]
colnames(drug_target)=c("Target","drugbank_id")
drug_target=drug_target[,c("drugbank_id","Target")]
drug_target=separate_rows(drug_target, drugbank_id, Target, sep=";",convert = TRUE)
drug_target$drugbank_id=gsub(" ","",drug_target$drugbank_id)
drug_target=drug_target[drug_target$drugbank_id %in% db_n_b3db_dbids,]

db_names=db_complete %>% dplyr::select(drugbank_id,name)
db_names$name=tolower(db_names$name)
drug_target=left_join(drug_target,db_names)
drug_target=drug_target %>% dplyr::select(name,Target)
colnames(drug_target)=c("Drug","Target")

# merge two sources
drug_target=rbind(drug_target,dgib_n_b3db)
write.table(drug_target,file="tables/drug_target_network.txt",sep="\t",quote=F,row.names=F,col.names=T)

# baits = module + risk genes
baits_path="data/GNN_labels"

network_proximity_tests = function(graphs, dt_df, disease_name, N)
{

  results = data.frame()
  for (graph_name in names(graphs))
  {
    cellname=gsub("_merged.csv","",graph_name)
    if(cellname %in% ExSub)
    {
      cellname="excitatory"
    }
    if(cellname %in% InSub)
    {
      cellname="inhibitory"
    }
    # cell type disease gene from labels
    label_file_name=paste(baits_path,paste(cellname,paste(disease_name,".txt",sep=""),sep="_"),sep="/")
    disease_labels=read.csv(label_file_name)
    disease_labels=disease_labels %>% filter(label == "positive") %>% select(gene)

    g.info = graphs[[graph_name]]
    dt_df_sub = dt_df[dt_df$Target %in% g.info$node,]
    drugs = unique(dt_df_sub$Drug)
    k=0
    for (drug in drugs)
    {
      k=k+1
      t=length(drugs)
      msg1=paste("Processed ", paste((k/t)*100,"% of drugs",sep=" "),sep="")
      msg2=paste(" for celltype ", cellname ,sep="")
      msg=paste(msg1,msg2)
      message(msg)
      graph = g.info$graph
      node = g.info$node
      disease_genes = intersect(disease_labels$gene, node)
      drug_targets = intersect(dt_df_sub[dt_df_sub$Drug == drug,]$Target, node)
      DG_degree_distribution=computeDegreeDistribution(disease_genes,g.info)
      DT_degree_distribution = computeDegreeDistribution(drug_targets,g.info)
      if (length(disease_genes) > 0 && length(drug_targets) > 0)
      {
        distance_matrix = distances(graph, v = V(graph)[drug_targets], to = V(graph)[disease_genes])
        minimum = computeMinimum(distance_matrix)
        proximity = mean(minimum, na.rm = T)
        random_distribution_proximity = computeRandomProximity(DT_degree_distribution,DG_degree_distribution,g.info)
        pval = computeStatistics(random_distribution_proximity,proximity)
        results = rbind(results, data.frame(
          cell = graph_name,
          drug = drug,
          proximity = proximity,
          z = pval[[2]],
          p = pval[[1]]
        ))
      }
      message(" ")
      message(paste(drug, pval[[1]],sep="::"))
    }
    message("Finished graph: ", graph_name)
  }
  return(results)
}



# Process GRNs and get the list of graphs
#asd_graphs = process_grn_files("data/PEC2_NetMed_GRNs/subclass/ASD")
#scz_Urban_graphs = process_grn_files("data/PEC2_NetMed_GRNs/subclass/SCZ")
bp_Urban_graphs = process_grn_files("data/PEC2_NetMed_GRNs/subclass/BPD")

# Run network proximity tests
#proximity_results_asd = network_proximity_tests(asd_graphs, drug_target, "UCLA_ASD_ASD",100)
#proximity_results_scz = network_proximity_tests(scz_Urban_graphs, drug_target, "CMC_SCZ",100)
proximity_results_bpd = network_proximity_tests(bp_Urban_graphs, drug_target, "Urban_DLPFC_BPD",100)


#write.table(proximity_results_asd %>% filter(p <=0.1,), file="tables/proximity_results.signif.subclass.ASD.csv",quote=F,sep="\t",row.names=F)
#write.table(proximity_results_scz %>% filter(p <=0.1,), file="tables/proximity_results.signif.subclass.SCZ.csv",quote=F,sep="\t",row.names=F)
write.table(proximity_results_bpd %>% filter(p <=0.1,), file="tables/proximity_results.signif.subclass.BPD.csv",quote=F,sep="\t",row.names=F)

#ASD: 5.5 hours
#SCZ:6.6 hours
end_time = Sys.time()
duration = end_time - start_time
print(paste("Total time taken:", duration))
