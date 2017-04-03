__author__ = 'jgw87'
"""
Take my scaffold assignments and filter them by the consensus to come out of step 3_
    -i input (raw) scaffold assignments
    -l linker file
    -c consensus scaffold assignments
    -o output scaffold assignments
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys


# sys.argv[1:] = ["-i","0_som_scaffolds.txt",
#                 "-l","1a_lingake_group_equivalencies_MANUAL.txt",
#                 "-c","3_consensus_linkage_groups.txt",
#                 "-o","3a_som_scaffolds_filtered.txt"]


def main():
    infile, linkfile, consensusfile, outfile = parse_args()
    linker = make_linker(linkfile)
    consensus = load_consensus(consensusfile)
    filter_scaffolds(infile, linker, consensus, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-l", "--linkfile")
    parser.add_argument("-c", "--consensusfile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    return args.infile, args.linkfile, args.consensusfile, args.outfile


def load_consensus(consensusfile):
    print("Loading consensus linkage groups from",consensusfile)
    consensus=dict()
    CON = open(consensusfile, "r")
    header = CON.readline().strip().split("\t")
    scaffoldID, lgID = get_headerIDs(header, ["scaffold", "lg"])
    for line in CON:
        data = line.strip().split("\t")
        scaffold, lg = data[scaffoldID], data[lgID]
        if scaffold in consensus:
            print("\tWARNING!",scaffold,"already added to consensus list!")
        consensus[scaffold] = lg
    print("\t",len(consensus),"scaffolds loaded")
    return consensus


#Go through and remove scaffolds that are (1) not in the consensus, or (2) on the wrong linkage group
def filter_scaffolds(infile, linker, consensus, outfile):
    print("Filtering scaffolds from",infile)
    #Set up things
    IN=open(infile, "r")
    OUT = open(outfile, "w")
    header= IN.readline()
    header_data = header.strip().split("\t")
    scaffoldID, lgID = get_headerIDs(header_data, ["scaffold", "lg"])
    OUT.writelines(header)

    #Go through each scaffold and filter
    notfound=0
    wrongplace=0
    right=0
    for line in IN:
        data = np.array(line.strip().split("\t"), dtype=str)    #Split line into an array of strings
        scaffold, lg = data[scaffoldID], data[lgID]
        if scaffold not in consensus:   #Skip if not part of consensus
            notfound+=1
            continue
        new_lg = linker[infile][lg]
        if new_lg == consensus[scaffold]:   #If on correct LG, then write out
            right+=1
            OUT.writelines(line)
        else: wrongplace+=1
    OUT.close()
    print("\tFound",right,"scaffolds to pass on;",notfound,"not found and",wrongplace,"were in the wrong place")


def get_headerIDs(header, names):
    header= np.array(header)
    names = np.array(names)
    ids = np.zeros(len(names), dtype=int)
    for i in range(len(names)):
        myID = np.where(header==names[i])
        if len(myID) == 0:
            print("WARNING! Header name",names[i],"not found in header data:",header)
        if len(myID) > 1 or len(myID[0]) > 1:
            print("Warning! Trying to find header",names[i],"results in multiple matches from np.where:",myID)
        ids[i] = int(myID[0][0])
    return ids

#Load a file showing which LGs match to each across populations; copied from step 3_
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
            mylg = lg[i]
            reflg = lg[0]
            linker[mymap][mylg] = reflg
    LINK.close()
    return linker

if __name__ == "__main__":
    main()