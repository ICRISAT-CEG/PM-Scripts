__author__ = 'jgw87'
"""
Take a genetic map and renumber its linkage groups as needed
    -i input scaffold map
    -o Output scaffold map
    -l Linker file showing which linkage groups correspond across the different maps. Should be two columns; First is the
        reference, and second is what will be renumbered to match
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

# sys.argv[1:] = ["-i","1b_841_core_scaffold_bins.txt",
#                 "-o","1e_841_core_scaffold_bins_standardized.txt",
#                 "-l","1d_841_linkage_correspondance.txt"]

def main():
    infile, outfile, linkfile = parse_args()
    linker = make_linker(linkfile)
    renumber_scaffolds(infile, linker, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-l", "--linkfile")
    args = parser.parse_args()
    return args.infile, args.outfile, args.linkfile


#Load a file showing which LGs match to each across populations
def make_linker(linkfile):
    print("Creating linker structure from",linkfile,"; assuming first column is reference, so the second will be renumbered to match first")
    linker = dict()
    LINK = open(linkfile, "r")
    #Parse header
    header = LINK.readline()
    #Parse linkage groups
    for line in LINK:
        lg = line.strip().split("\t")
        for i in range(len(lg)):
            mylg = lg[1]
            reflg = lg[0]
            linker[mylg] = reflg
    LINK.close()
    return linker


#Creates a nested dictionary of scaffold -> lg -> count, where count will be used to determine the consensus
def renumber_scaffolds(infile, linker,outfile):
    scaffolds=dict()
    print("Loading scaffold data from",infile)

    #Find appropriate columns in header
    IN = open(infile, "r")
    OUT = open(outfile, "w")
    header = IN.readline()
    OUT.write(header)

    #Process header to find correct columns
    header = header.strip().split("\t")
    header=np.array(header, dtype=str)
    lgID = np.where(header=="lg")[0] # [0] b/c np.where returns tuple
    if len(lgID) > 1 :
        print("WARNING! More than one header column found for linkage groups (",lgID)
    else:   # Reduce to simple ints
        lgID = lgID[0]
    #print("LG at column", lgID, "; scaffold at column",scaffoldID)

    #Parse data line by line
    for line in IN:
        data = line.strip().split("\t")
        data = np.array(data, dtype=str)
        mylg = data[lgID]
        new_lg = linker[mylg]
        data[lgID] = new_lg
        OUT.write("\t".join(data) + "\n")
        #print("LG",lg,"in",infile,"becomes",new_lg)

    IN.close()
    OUT.close()

#
# def determine_consensus_linkage_group(lg, linker, infile):
#     return linker[infile][lg]
#
#
# def consolidate_scaffolds(scaffolds):
#     print("Determining consensus linkage group by majority rules")
#     consensus = dict()
#     for scaffold in scaffolds:
#         lg = scaffolds[scaffold]
#         #print("\t",scaffold, ":",lg)
#         for mylg in lg:
#             if lg[mylg] > len(lg) / 2:
#                 consensus[scaffold] = mylg
#                 #print("\t\tAssigned to",mylg)
#     print("\t",len(consensus),"of",len(scaffolds),"assigned (",str(len(consensus)/len(scaffolds)),"%)")
#     return consensus
#
# def output_consensus(consensus, outfile):
#     print("Outputting results to", outfile)
#     #Build lists
#     scaffolds, lg = list(), list()
#     for scaff in consensus:
#         scaffolds.append(scaff)
#         lg.append(consensus[scaff])
#     output = pd.DataFrame()
#     output["scaffold"] = scaffolds
#     output["lg"] = lg
#     output = output.sort(columns = ["lg", "scaffold"])
#     output.to_csv(outfile, sep="\t", index=False)




if __name__ == "__main__":
    main()