#!/bin/bash

#all gexprs
DIRECTORY="data"

for file in "$DIRECTORY"/*
do
  echo "Processing $file"
  base=$(basename "$file" .mat)
  newfile="${file%/*}/$base""_grnboost.tsv"
	echo $newfile
  python ~/anaconda3/envs/pyscenic/bin/arboreto_with_multiprocessing.py "$file" allTFs_hg38.txt --method grnboost2  --output $newfile --seed 123
done
