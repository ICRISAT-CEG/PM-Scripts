__author__ = 'jgw87'
"""
Take the distance matrices for each scaffold combination and determine the best anchor scaffold for the unmapped ones
    -i Comma-separated list of input files (distance matrices)
    -m List of core map scaffolds (from 5e_)
    -o Final output file
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys
from collections import Counter
sys.path.append('/home/jgw87/Software/Python')
import JasonUtils

# sys.argv[1:]=["-i", "bak/6b_segregating_841_scaffold_distances.txt,bak/6b_segregating_boubacar_scaffold_distances.txt,bak/6b_segregating_som_scaffold_distances.txt",
#           "-m","5e_core_map_scaffolds.txt",
#           "-o","6c_extended_map_scaffolds",
#           "-c", "0.4"]


def main():
    infiles, mapfile, outfile, cutoff = parse_args()
    map = load_map(mapfile)
    scaffolds, matrices = load_matrices(infiles)
    anchors = find_anchors(scaffolds, matrices, map, cutoff)
    output_extended_map(anchors, map, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles")
    parser.add_argument("-m", "--mapfile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-c", "--cutoff", type=float)
    args = parser.parse_args()
    return args.infiles.split(","), args.mapfile, args.outfile, args.cutoff


def load_map(mapfile):
    print("Loading map from",mapfile)
    map=dict()
    MAP=open(mapfile, "r")
    scaffoldID, lgID, posID = JasonUtils.get_headerIDs_from_file(MAP, ["scaffold","lg","pos"])
    for line in MAP:
        data = line.strip().split("\t")
        scaffold, lg, pos = data[scaffoldID], data[lgID], data[posID]
        if scaffold == "unknown": continue  #Skip unknown ones
        if scaffold in map:
            print("WARNING! Scaffold",scaffold, "already loaded and will be overwritten!")
        map[scaffold] = dict()
        map[scaffold]["lg"] = lg
        map[scaffold]["pos"] = pos
    MAP.close()
    print("\tLoaded",len(map),"scaffolds")
    return map

def load_matrices(infiles):
    scaffolds = list()
    matrices = list()
    for infile in infiles:
        print("Loading",infile)
        #data = pd.read_csv(infile, index_col=0, sep="\t", header=None, nrows=100) #For debugging
        #rsq = data.as_matrix()[:,:100]  #For debugging
        data = pd.read_csv(infile, index_col=0, sep="\t", header=None)
        rsq = data.as_matrix()
        names = np.array(data.index, dtype=object)
        matrices.append(rsq)
        scaffolds.append(names)
    print("Loaded",len(matrices),"distance matrices")
    return scaffolds, matrices

def find_anchors(scaffolds, matrices, map, cutoff):
    print("Finding anchors")
    anchors = dict()
    tomap = get_anchors_to_map(scaffolds, map)
    for i in range(len(scaffolds)):
        myscaff = scaffolds[i]
        mymatrix = matrices[i]

        #Set self-LD to 0 to avoid only associating with itself
        for i in range(len(mymatrix)):
            mymatrix[i,i] = 0
        #Set LD to 0 for non-anchor scaffolds
        is_anchor = np.array([s in map for s in myscaff])
        mymatrix[:,~is_anchor] = 0


        totest = np.where([s in tomap for s in myscaff])[0]
        for t in totest:
            myname = myscaff[t]
            mymax = np.nanmax(mymatrix[t,:])
            if mymax < cutoff:  #Skip ones with too low LD
                continue
            max_i = np.where(mymatrix[t,:] == mymax)[0]
            myanchor = list(np.array(myscaff)[max_i])
            if myname not in anchors:
                anchors[myname] = list()
            anchors[myname] += myanchor
    print(anchors)
    return(anchors)


def get_anchors_to_map(scaffolds, map):
    tomap = set()
    for myscaff in scaffolds:
        for s in myscaff:
            if s not in map:
                tomap.add(s)
    print("\tFound",len(tomap),"scaffolds to be mapped")
    return list(tomap)


def output_extended_map(anchors, map, outfile):
    print("Outputting best anchors for",len(anchors), "scaffolds")
    OUT=open(outfile, "w")
    OUT.writelines( "\t".join(["scaffold","lg","pos", "best_anchor", "all_anchors"]) + "\n")
    for scaffold in anchors:

        all_anchors = set(anchors[scaffold])
        best_anchor = Counter(anchors[scaffold]).most_common(1)[0][0]
        print("All:",anchors[scaffold])
        print("\tBest:",best_anchor)
        mylg = map[best_anchor]["lg"]
        mypos = map[best_anchor]["pos"]
        OUT.writelines("\t".join([scaffold, mylg, mypos, best_anchor, ",".join(all_anchors)]) + "\n")

    OUT.close




if __name__ == "__main__":
    main()