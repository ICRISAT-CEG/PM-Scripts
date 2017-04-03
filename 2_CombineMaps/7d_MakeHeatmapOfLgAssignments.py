__author__ = 'jgw87'
"""
Take the list of LG assignments from 1c_ and turn it into a heatmap
    -i Input LG assignment file from 1c_
    -o Output file prefix (for both .txt and .png)
    -l Linker file showing equivalencies in LGs across my maps
"""

import argparse
import matplotlib
#import matplotlib.colors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

sys.path.append('/home/jgw87/Software/Python')
import JasonUtils


# sys.argv[1:] = ["-i","1c_contig_lg_assignments_som.txt",
#                 "-o","1d_som_matrix",
#                 "-l","0_my_lingake_group_equivalencies_MANUAL.txt"]


def main():
    infile, outprefix = parse_args()
    matrix, consensus = parse_lg_mappings(infile)
    output_matrix(matrix, consensus, outprefix)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outprefix")
    args = parser.parse_args()
    return args.infile, args.outprefix


class map_matrix:
    matrix = None
    contigs = None
    lg = None


    def __init__(self, nrow, ncol, mycontigs, mylg):
        self.matrix = np.zeros((nrow, ncol), dtype=int)
        self.contigs = mycontigs
        self.lg = mylg

    def add_contig(self, contig, group):
        self.matrix[self.contigs == contig, self.lg == group] += 1

    def normalize(self):
        for i in range(len(self.matrix)):
            total = np.sum(self.matrix[1,:])
            if total > 0:
                self.matrix[i,:] = self.matrix[i,:] / total


#Take the listings of LGs and tabulate them
def parse_lg_mappings(infile):
    print("Loading LG assignments")
    data = pd.read_csv(infile, sep="\t")
    data =data.sort_index(axis=0, by=["consensus_lg", "contig"])

    #Create matrix
    contigs = data["contig"].unique()
    lg = np.sort(data["lg"].unique())
    lg = lg[lg != "unknown"]
    #print(lg)
    nrow = len(contigs)  #Rows and columns of matrix. Rows = contigs from consensus map, cols = my LG
    ncol = len(lg)
    matrix = map_matrix(nrow=nrow, ncol=ncol, mycontigs=contigs, mylg=lg)
    print("Created matrix with",nrow,"contigs among",ncol,"linkage groups")

    #Fill matrix and create key of consensus marker locations
    consensus = dict()
    for i in range(len(data)):
        mylg = data["lg"].iloc[i]
        mycontig = data["contig"].iloc[i]
        myconsensus = data["consensus_lg"].iloc[i]
        myrank = data["rank"].iloc[i]

        #Record consensus LG
        if mycontig in consensus and consensus[mycontig] != myconsensus:
            print("WARNING! Mismatch in consensus lg for",mycontig,":",consensus[mycontig], myconsensus)
        consensus[mycontig] = myconsensus

        #Skip if no mapped LG found
        if mylg == "unknown":
            continue
        #mylg = linker[infile][mylg]


        if myrank == "1" or myrank.startswith("1 of "):  #Only take the best ranking
            matrix.add_contig(contig = mycontig, group = mylg)
    matrix.normalize()
    return matrix, consensus


def output_matrix(matrix, consensus, outprefix):
    #Output as text
    outfile = outprefix + ".txt"
    OUT = open(outfile, "w")
    header = ["contig", "consensusLG"] + list(matrix.lg)
    OUT.writelines("\t".join(header) + "\n")
    for i in range(len(matrix.contigs)):
        mycontig = matrix.contigs[i]
        line = [mycontig, consensus[mycontig]] + list(matrix.matrix[i,:].astype(dtype=str))
        OUT.writelines("\t".join(line) + "\n")
    OUT.close()

    #Output as heatmap
    fig = plt.figure()
    ax = fig.add_subplot(111, title="Heatmap of LG assignments", xlabel="My map LG", ylabel="Consensus Contigs")
    ax.pcolormesh(matrix.matrix, cmap=plt.cm.Blues)
    ax.invert_yaxis()
    ##Add colorbars based on linkage group
    colors = get_consensus_colors(matrix, consensus)
    ax.scatter(np.zeros(len(matrix.contigs)), range(len(matrix.contigs)), linewidths=0, c=list(colors), cmap=plt.cm.rainbow_r)
    fig.savefig(filename = outprefix + ".png", transparent=True)
    plt.close()


def get_consensus_colors(matrix, consensus):
    #First fill with just the linkage group while assembling unique ones into a set
    groups = set()
    colors = np.ndarray(len(matrix.contigs), dtype=object)
    for i in range(len(matrix.contigs)):
        colors[i] = consensus[matrix.contigs[i]]
        groups.add(colors[i])
    groups = sorted(groups)

    #Then put in actual values for each
    for i in range(len(groups)):
        colors[colors == groups[i]] = i
    return colors


#
##Load a file showing which LGs match to each across populations. Returns a nested dict() of file->map_lg = standard_lg
#def make_linker(linkfile):
#    print("Creating linker structure from", linkfile, "; assuming first column is default")
#    linker = dict()
#    LINK = open(linkfile, "r")
#    #Parse header
#    header = LINK.readline().strip().split("\t")  # Separate into entries
#    for h in header:
#        linker[h] = dict()
#    #Parse linkage groups
#    for line in LINK:
#        lg = line.strip().split("\t")
#        for i in range(len(lg)):
#            mymap = header[i]
#            mylg = lg[i]
#            reflg = lg[0]
#            linker[mymap][mylg] = reflg
#    LINK.close()
#    return linker


if __name__ == "__main__":
    main()