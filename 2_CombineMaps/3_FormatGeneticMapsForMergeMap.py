__author__ = 'jgw87'
"""
Take a genetic map and reformat it for MergeMap
    -i Scaffold assignments for a specific map
    -o Output file name with consolidated scaffold locations
    -l Linker file showing which linkage groups correspond across the different maps. Each column should be labeled with
       the name of the corresponding file from -i
    -g Linkage group number (optional; if given, will only output the specified group)
    -b Bin divisor - how much to divide the bin value by to recreate the centimorgan range
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

# sys.argv[1:] = ["-i","0_som_scaffolds.txt",
#                 "-o","2_som_mergemap.txt",
#                 "-l","1a_lingake_group_equivalencies_MANUAL.txt"]

def main():
    infile, outfile, linkfile, group, bindiv = parse_args()
    linker = make_linker(linkfile)
    format_map(infile, linker, outfile, group, bindiv)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-l", "--linkfile")
    parser.add_argument("-g", "--group", required=False, type=int)
    parser.add_argument("-b", "--bindiv", required=False, type=int)
    args = parser.parse_args()
    return args.infile, args.outfile, args.linkfile, args.group, args.bindiv


#Load a file showing which LGs match to each across populations
def make_linker(linkfile):
    print("Creating linker structure from",linkfile,"; assuming first column is default")
    linker = dict()
    LINK = open(linkfile, "r")
    #Parse header
    header = LINK.readline().strip().split("\t")    # Separate into entries
    for h in header:
        linker[h] = dict()
    #Parse linkage groups
    for line in LINK:
        lg = line.strip().split("\t")
        for i in range(len(lg)):
            mymap = header[i]
            mylg = int(lg[i])
            reflg = int(lg[0])
            linker[mymap][mylg] = reflg
    LINK.close()
    return linker


def format_map(infile, linker, outfile, group, bindiv):
    print("Converting",infile,"to MergeMap format in",outfile)
    OUT = open(outfile, "w")
    data =pd.read_csv(infile, index_col=False, sep="\t")
    #new_lg = data["lg"]     #May need to fix this later; for now assuming MergeMap handles it correctly
    new_lg = match_lg(data, linker, infile)
    data["lg"] = new_lg

    #Go through and output one LG at a time
    lgs = np.unique(new_lg)
    if group is not None:   # If group given, output just that one
        lgs = [group]
    for lg in sorted(lgs):
        sub = data.loc[data["lg"] == lg]
        print("LG",lg,"has",len(sub),"scaffolds")
        output_lg("lg" + str(int(lg)), sub, OUT, bindiv)
    OUT.close()


#Figure out which consensus linkage groups
def match_lg(data, linker, infile):
    new_lg = np.zeros(len(data))
    lg_list = np.unique(data["lg"])
    for lg in lg_list:
        current = np.array(data["lg"] == lg)
        new = linker[infile][lg]
        new_lg[current] = new
    return new_lg


def output_lg(name, data, OUT, bindiv):
    data = data.sort(columns=["lg", "bin"])
    first_binval = np.min(data["bin"])

    OUT.writelines("group " + name + "\n")
    OUT.writelines(";BEGINOFGROUP\n")
    for i in range(len(data)):
        scaffold = data["scaffold"].iloc[i]
        binval = data["bin"].iloc[i] - first_binval
        binval /= bindiv
        OUT.writelines(scaffold + "\t" + str(binval) + "\n")
    OUT.writelines(";ENDOFGROUP\n\n")


# #Creates a nested dictionary of scaffold -> lg -> count, where count will be used to determine the consensus
# def load_scaffolds(infiles, linker):
#     scaffolds=dict()
#     for infile in infiles:
#         print("Loading scaffold data from",infile)
#
#         #Find appropriate columns in header
#         IN = open(infile, "r")
#         header = IN.readline().strip().split("\t")
#         header=np.array(header, dtype=str)
#         lgID, scaffoldID = np.where(header=="lg")[0], np.where(header=="scaffold")[0] # [0] b/c np.where returns tuple
#         if len(lgID) > 1 or len(scaffoldID) > 1:
#             print("WARNING! More than one header column found for linkage groups (",lgID,"or scaffolds",scaffoldID)
#         else:   # Reduce to simple ints
#             lgID = lgID[0]
#             scaffoldID = scaffoldID[0]
#         #print("LG at column", lgID, "; scaffold at column",scaffoldID)
#
#         #Parse data line by line
#         for line in IN:
#             data = line.strip().split("\t")
#             data = np.array(data, dtype=str)
#             lg, scaffold = data[lgID], data[scaffoldID]
#             new_lg = determine_consensus_linkage_group(lg, linker, infile)
#             #print("LG",lg,"in",infile,"becomes",new_lg)
#             if scaffold not in scaffolds:   # Initialize scaffold
#                 scaffolds[scaffold] = dict()
#             if new_lg not in scaffolds[scaffold]:   #Initialize scaffold->lg
#                 scaffolds[scaffold][new_lg] = 0
#             scaffolds[scaffold][new_lg] += 1
#         IN.close()
#     return scaffolds
#
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