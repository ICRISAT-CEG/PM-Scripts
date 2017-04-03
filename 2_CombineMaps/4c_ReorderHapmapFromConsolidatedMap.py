__author__ = 'jgw87'
"""
Take a hapmap and the list of SNP-scaffold assignments from 5a_ and reorder the hapmap based on the new map order
    -i Input hapmap
    -s SNP-scaffold assignment (columns of snp name, scaffold, linkage group, and position)
    -o Output hapmap file
    -m Multiplier for position (optional, default 100)
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys
import gzip
sys.path.append('/home/jgw87/Software/Python')
import JasonUtils

# sys.argv[1:] = ["-i","5b_core_scaffold_snps.hmp.txt",
#                 "-s","5a_consolidated_map_snps.txt",
#                 "-o","5c_consolidated_map_reordered.hmp.txt",
#                 "-m","100"]


def main():
    infile, snpfile, outfile, multiplier = parse_args()
    snpkey = load_snp_key(snpfile, multiplier)
    hmp = reorder_hapmap(infile, snpkey)  #Returns a pd.DataFrame
    hmp.to_csv(outfile, sep="\t", index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-s", "--snpfile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-m", "--multiplier", required=False, default=1, type=int)
    args = parser.parse_args()
    return args.infile, args.snpfile, args.outfile, args.multiplier

def load_snp_key(snpfile, multiplier):
    print("Loading SNP key to order hapmap from",snpfile)
    SNPS = open(snpfile, "r")

    #Parse file
    snps = dict()
    snpID, scaffoldID, lgID, posID = JasonUtils.get_headerIDs_from_file(SNPS, ["snp","scaffold","lg","pos"])
    for line in SNPS:
        data = line.strip().split("\t")
        snp, scaffold, lg, pos = data[snpID], data[scaffoldID], data[lgID], data[posID]
        snps[snp] = dict()
        snps[snp]["lg"] = lg
        snps[snp]["scaffold"] = scaffold
        snps[snp]["pos"] = str(int(float(pos) * multiplier))    #Multiply, round, and then turn to string
    SNPS.close()
    print("\tLoaded",len(snps),"snps")
    return snps

def reorder_hapmap(infile, snpkey):
    print("Reordering hapmap")
    print("\tLoading file...",infile)
    hmp=pd.read_csv(infile, sep="\t")
    #hmp=pd.read_csv(infile, sep="\t", nrows=1000)   #For debugging

    #Extract the columns I'll be working with
    snps = np.array(hmp["rs#"])
    scaffolds = np.empty(shape=len(snps), dtype=object)
    lgs = np.array(hmp["chrom"])
    pos = np.array(hmp["pos"])

    #Go through and update everything
    print("\tUpdating information...")
    #print(snps)
    for i in range(len(snps)):
        if snps[i] not in snpkey:
            print("\tWarning!",snps[i],"not in key!")
            continue
        scaffolds[i] = snpkey[snps[i]]["scaffold"]
        lgs[i] = snpkey[snps[i]]["lg"]
        pos[i] = snpkey[snps[i]]["pos"]

    #Put data back in hapmap and sort
    print("\tSorting...")
    hmp["rs#"] = snps   #Shouldn't change anything, but just to be safe
    hmp["alleles"] = scaffolds
    hmp["chrom"] = lgs
    hmp["pos"] = pos
    hmp = hmp.sort(columns=["chrom","pos","rs#"])

    return hmp


if __name__ == "__main__":
    main()