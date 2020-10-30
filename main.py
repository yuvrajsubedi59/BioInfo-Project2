# Intro to BioInformatics
# Project 2

import os
import numpy

# Read Covid file
print("Reading the Cov2 SARS fasta file...\n")
with open('CoV.fasta') as fileCov:
    commentsCov = fileCov.readline()
    SARS_COV2_genome = fileCov.read()

fileCov.close()

# validate the file and extract nsp1 gene
SARS_COV2_genome = SARS_COV2_genome.replace('\n', '')
Cov_nsp1 = SARS_COV2_genome[266:805]
Cov_nsp1 = Cov_nsp1.upper()

# Read Mers file
print("Reading the MERS fasta file...\n")
with open('MERS.fasta') as fileMers:
    commentsMers = fileMers.readline()
    Mers_genome = fileMers.read()
fileMers.close()

# validate the file and extract nsp1 gene
Mers_genome = Mers_genome.replace('\n', '')
Mers_nsp1 = Mers_genome[279:857]
Mers_nsp1 = Mers_nsp1.upper()
