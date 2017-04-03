__author__ = 'jgw87'
"""
Take the scaffold assignments from 1b_ and put with linkage groups I generated
    -i Primer scaffold assignments
    -l My linkage groups (should be something like "6d_[popname]_scaffold_bins.txt"
    -o  Output file
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

# sys.argv[1:] = ["-i", "1b_primer_scaffold_assignments.txt",
#                   "-l", "/home/jgw87/Working_Files/GBS/Analysis/PearlMillet/20140324_AlignToRefseqV0.3/MakeNewMaps/SomData/6d_som_scaffold_bins.txt",
#                   "-o", "1c_primer_lg_assignments.txt",
#                   "-c", "0_Rajaman_consensus_map_ssr_assignments.txt"]

contigID, scaffoldID = 0, 1 #Columns in primer input file with needed info
ssrID, consID   = 0, 1  #Columns in consensus map file with name of SSR and assigned linkage group

def main():
    infile, linkfile, outfile, consensusfile = parse_args()
    scaffolds = load_linkage_groups(linkfile)
    contigs = assign_linkage_groups(infile, scaffolds)
    contigs = connect_to_consensus(contigs, consensusfile)
    contigs = contigs.sort(columns=["rank", "lg","bin"])
    contigs.to_csv(outfile, sep="\t", index=False)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-l", "--linkfile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-c", "--consensusfile")
    args = parser.parse_args()
    return args.infile, args.linkfile, args.outfile, args.consensusfile


def load_linkage_groups(linkfile):
    IN = open(linkfile, "r")
    #Parse header and find appropriate columns
    header = np.array(IN.readline().strip().split("\t"))   #Clear header
    binID = np.where((header == "binval") | (header == "bin") | (header =="pos"))[0]
    lgID = np.where((header == "lg")| (header == "linkage_group"))[0]
    scaffoldID = np.where(header == "scaffold")[0]
    print(binID, lgID, scaffoldID)

    #Go through and save data
    scaffolds = dict()
    for line in IN:
        data=line.strip().split("\t")
        bin, lg, scaffold = data[binID], data[lgID], data[scaffoldID]
        if scaffold in scaffolds:
            print("WARNING!",scaffold,"already loaded into disctionary of scaffolds!")
        scaffolds[scaffold] = (lg, bin)

    print("Loaded data for",len(scaffolds),"scaffolds")
    return scaffolds


def assign_linkage_groups(infile, scaffolds):
    #Set up variables
    bins = list()
    lg = list()
    contigs = pd.read_csv(infile, sep="\t")
    #contigs=contigs[["contig","scaffold"] ]

    #Go through and create new data
    for i in range(len(contigs)):
        scaffold = contigs["scaffold"].iloc[i]
        mylg, mybin = "unknown", "unknown"
        if scaffold in scaffolds:   # Reassign if known
            mylg, mybin = scaffolds[scaffold]
        bins.append(mybin)
        lg.append(mylg)

    #Add to data
    contigs["lg"] = lg
    contigs["bin"] = bins
    contigs = contigs.sort(columns = ["lg", "bin"])
    return contigs


def connect_to_consensus(contigs, consensusfile):
    IN = open(consensusfile, "r")
    IN.readline()   #Clear header
    contigs["consensus_lg"] = "unknown"
    print(contigs["contig"])
    for line in IN:
        data= line.strip().split("\t")
        ssr, lg = data[ssrID], data[consID]
        ssr = ssr.upper().lstrip("X")   # Format same as primers
        print("Looking for",ssr, lg)
        if np.any(contigs["contig"] == ssr):
            i = contigs["contig"] == ssr
            if np.sum(i) != 1:
                print("WARNING! More than one contig set found for",ssr,":",np.sum(i))
            contigs["consensus_lg"].loc[i] = lg

    return contigs



if __name__ == "__main__":
    main()