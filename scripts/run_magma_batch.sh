#!/bin/bash
# please follow magma instructions to create the correct dir structure

traits=("scz" "bpd" "asd" "epilepsy")
celltypes=("astro" "endo" "micro" "oligo" "opc" "vlmc" "excitatory" "inhibitory") # for major cell type class

celltype_dir="./"

for celltype_file in ${celltype_dir}*; do
    celltype=$(basename ${celltype_file})

    for trait in "${traits[@]}"; do
        geneset="${celltype_file}"
        # Perform gene-set analysis with MAGMA
        ~/tools/magma/magma \
            --gene-results ${trait}.genes.raw \
            --set-annot ${geneset} \
            --out ${trait}_${celltype}
    done
done
