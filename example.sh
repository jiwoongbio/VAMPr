#!/bin/env bash

# 1. Download an example input genome FASTA file
wget --no-check-certificate https://cdc.biohpc.swmed.edu/VAMPr/example/Klebsiella_pneumoniae_MGH_43.fasta

# 2. Download VAMPr protein database
perl VAMP_database.pl

# 3. Identify AMR genotypes
perl VAMP.pl Klebsiella_pneumoniae_MGH_43.fasta > Klebsiella_pneumoniae_MGH_43.VAMP.txt

# 4. Download VAMPr prediction model
wget --no-check-certificate https://cdc.biohpc.swmed.edu/VAMPr/VAMP_model/573.amikacin.VAMP_model.RData
wget --no-check-certificate https://cdc.biohpc.swmed.edu/VAMPr/VAMP_model/573.cefepime.VAMP_model.RData

# 5. Predict AMR phenotypes
perl VAMP_predict.pl 573.amikacin.VAMP_model.RData amikacin=Klebsiella_pneumoniae_MGH_43.VAMP.txt
perl VAMP_predict.pl 573.cefepime.VAMP_model.RData cefepime=Klebsiella_pneumoniae_MGH_43.VAMP.txt
