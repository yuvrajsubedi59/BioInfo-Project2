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
    for x in gene_matx[0]:
        print(x)
    for x in gene_matx[1]:
        print(x)
    trace_back(gene_matx,cov_nsp1, mers_nsp1)

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
#returns: a tuple of matrixes with the comparisons and direction
def get_alignment_matrix(gene1, gene2):
    #alignment scores
    gap = -1
    mismatch = -1
    match = 1
    #get length plus one for size of matrix
    g1_len=len(gene1)+1
    g2_len=len(gene2)+1
    # initialize gene_matx to 0
    gene_matx=[[ 0 for x in range(g2_len) ]  for y in range(g1_len)]
    gene_matx_dir=[[ "" for x in range(g2_len) ]  for y in range(g1_len)]
    x = 1

    # set the scores of the row that lines up with empty
    while x < g2_len:
        gene_matx[0][x] = x * gap
        gene_matx_dir[0][x] = "h"
        x+=1
    y = 1
    # set the scores of the column that lines up with empty
    while y < g1_len:
        gene_matx[y][0] = y * gap
        gene_matx_dir[y][0] = "v"
        y+=1
    #get the score for the rest of the spots
    y=1
    while y < g1_len:
        x=1
        while x < g2_len:
            #print("i:",i,", j:",j)
            #diagonel score
            diag = gene_matx[y-1][x-1] +  (match,mismatch) [gene1[y-1] == gene2[x-1]]
            #vertical score
            vert = gene_matx[y][x-1] + gap
            #horisontal score
            hori = gene_matx[y-1][x] + gap
            #get max score
            gene_matx[y][x] = max(diag, vert, hori)
            #store direction
            if diag == gene_matx[y][x]:
                gene_matx_dir[y][x] = "d"
            elif x>y:
                gene_matx_dir[y][x] = ("h","v") [hori>vert]
            else:
                gene_matx_dir[y][x] = ("v","h") [vert>hori]

            x+=1
        y+=1
    return (gene_matx,gene_matx_dir)
#trace_back traces the alignment in the gene_matx
#parameters: gene_matx: the tuple of gene matrixes
#           gene1,gene2 the genes we are compairing
#               must be in same order they were given to
# Returns: tuple containing aligned genes
def trace_back(gene_matxs,gene1,gene2):
    gene_matx = gene_matxs[0]
    gene_matx_dir = gene_matxs[1]
    for x in gene_matx:
        print(x)
    for x in gene_matx_dir:
        print(x)
    x_len = len(gene_matx[0])
    y_len = len(gene_matx)
    print('x length:',x_len,", y length:" ,y_len)
    x = x_len-1
    y = y_len-1
    #x = len(gene2)
    #y = len(gene1)
    print("x:",x,", y:",y)
    g1_align = ""
    g2_align = ""
    while y > 0 and x > 0:
        print("x:",x,", y:",y, ", dir:", gene_matx_dir[y][x] )
        if gene_matx_dir[y][x] == "d":
            g1_align = gene1[y-1] + g1_align
            g2_align = gene2[x-1] + g2_align
            x -= 1
            y -= 1
        elif gene_matx_dir[y][x] == "h":
                g1_align = "_" + g1_align
                g2_align = gene2[x-1] + g2_align
                x -= 1
        elif gene_matx_dir[y][x] == "v":
                g1_align = gene1[y-1] + g1_align
                g2_align = "_" + g2_align
                y -= 1
        else:
            x -= 1
        print("x:",x,", y:",y)
    print(g1_align)
    print(g2_align)
    return (g1_align,g2_align)


#run main when script is run
if __name__ == "__main__":
    main()
