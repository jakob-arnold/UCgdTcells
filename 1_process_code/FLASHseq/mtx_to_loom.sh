#!/bin/bash

# In house funtion to convert .mtx file to .loom

source activate pyscenic

python /home/jakob/pySCENIC/scripts/save_mtx_to_loom.py \
	--barcodes /home/jakob/Desktop/UC_Project/Paper/UCgdTcells/data_processed/FLASHseq/cells.tsv \
	--features /home/jakob/Desktop/UC_Project/Paper/UCgdTcells/data_processed/FLASHseq/genes.tsv \
	--matrix /home/jakob/Desktop/UC_Project/Paper/UCgdTcells/data_processed/FLASHseq/gex.mtx \
	--output /home/jakob/Desktop/UC_Project/Paper/UCgdTcells/data_processed/FLASHseq/gex.loom
	
