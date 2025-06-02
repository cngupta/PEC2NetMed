#!/bin/bash

traits=("scz" "bpd" "asd" "epilepsy")

for trait in "${traits[@]}"; do
    snploc="./gwas_files/${trait}.magma.input.snp.chr.loc.txt"
      ncbi38="./magma_files/geneloc/NCBI38.gene.loc"

    # Annotate SNPs
    ~/tools/magma/magma --annotate \
        --snp-loc ${snploc} \
        --gene-loc ${ncbi38} \
        --out ${trait}

    # Set sample size based on trait
    if [ "$trait" == "scz" ]; then
        N=76755
    elif [ "$trait" == "asd" ]; then
        N=46350
    elif [ "$trait" == "bpd" ]; then
        N=41917
    elif [ "$trait" == "epilepsy" ]; then
        N=29944
    else
        echo "Unknown trait: $trait"
        exit 1
    fi

    ref=./magma_files/ref/g1000_eur/g1000_eur
    ~/tools/magma/magma \
        --bfile $ref \
        --pval ./gwas_files/$trait.magma.input.snp.p.txt N=$N \
        --gene-annot ${trait}.genes.annot \
        --out ${trait}
done
