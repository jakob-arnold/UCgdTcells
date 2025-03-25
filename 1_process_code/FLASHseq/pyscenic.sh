#!/bin/bash

# Files
loom=/home/jakob/Desktop/UC_Project/Paper/UCgdTcells/data_processed/FLASHseq/gex.loom

tf=/home/jakob/reference_human/pySCENIC/TF_list/allTFs_hg38.txt
feather=/home/jakob/reference_human/pySCENIC/feather
motif=/home/jakob/reference_human/pySCENIC/Motif2TF/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

nCores=32

###

cd '/home/jakob/Desktop/UC_Project/Paper/UCgdTcells/data_processed/FLASHseq'

mkdir results

source activate pyscenic

# Step 1
 pyscenic grn \
	 $loom $tf \
	 --num_workers $nCores \
	 -o results/grn.csv


/home/jakob/Desktop/UC_Project/Paper/UCgdTcells/data_processed/FLASHseq/step2_3.sh
