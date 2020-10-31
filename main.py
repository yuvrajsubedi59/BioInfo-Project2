# Intro to BioInformatics
# Project 2
# By            and Fred Fikter

import os
import numpy

def main():
    #cov_nsp1 = get_cov()
    #mers_nsp1 = get_mers()
    #short sequences for testing
    cov_nsp1 = "GCATGCAT"
    mers_nsp1 = "CATTCAT"
    print(cov_nsp1)
    print(mers_nsp1)
    gene_matx=get_alignment_matrix(cov_nsp1, mers_nsp1)
    for x in gene_matx:
        print(x)

#function to get the covid gene
#parameters: none
#returns: a string containg the cov 2 gene
def get_cov():
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
    return Cov_nsp1

#function to get the covid gene
#parameters: none
#returns: a string containg the cov 2 gene
def get_mers():
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
    return Mers_nsp1

#get_alignment_matrix completes the first step of needleman-wunsch
#   and creates the alignment grid for the two genes
#parameters: the genes to compare
#returns: a matrix with the comparisons of the
def get_alignment_matrix(gene1, gene2):
    #alignment scores
    gap = -2
    mismatch = -1
    match = 1
    #get length plus one for size of matrix
    g1_len=len(gene1)+1
    g2_len=len(gene1)+1
    # initialize gene_matx to 0
    gene_matx=[[ 0 for y in range(g2_len) ]  for x in range(g1_len)]

    i = 1
    # set the scores of the row that lines up with empty
    while i < g1_len:
        gene_matx[0][i] = i * gap
        i+=1
    i = 1
    # set the scores of the column that lines up with empty
    while i < g2_len:
        gene_matx[i][0] = i * gap
        i+=1
    #get the score for the rest of the spots
    i=1
    while i < g1_len:
        j=1
        while j<g2_len-1:
            #diagonel score
            diag = gene_matx[i-1][j-1] +  (match,mismatch) [gene1[i-1] == gene2[j-1]]
            #vertical score
            vert = gene_matx[i][j-1] + gap
            #horisontal score
            hori = gene_matx[i-1][j] + gap
            #get max score
            gene_matx[i][j] = max(diag, vert, hori)
            j+=1
        i+=1

    return gene_matx

#run main when script is run
if __name__ == "__main__":
    main()
