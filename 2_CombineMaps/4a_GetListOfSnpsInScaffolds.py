__author__ = 'jgw87'
"""
Take the consolidated map and list of snp-scaffold assignments and get out the list of SNPs on the core map and their
positions
    -m Consolidated map file
    -s SNP-scaffold file from 5_ (2 columns, with SNP name and scaffold ID)
    -o Output file
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys
sys.path.append('/home/jgw87/Software/Python')
import JasonUtils

# sys.argv[1:] = ["-m","5_core_map_combined.txt",
#                 "-s","999_test_scaffolds.txt",
#                 "-o","5a_consolidated_map_snps.txt"]

def main():
    mapfile, snpfile, outfile = parse_args()
    scaffolds = load_map(mapfile)
    output_snps(scaffolds, snpfile, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mapfile")
    parser.add_argument("-s", "--snpfile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    return args.mapfile, args.snpfile, args.outfile


def load_map(mapfile):
    print("Loading map file from",mapfile)
    master = dict()
    MAP = open(mapfile, "r")
    header = MAP.readline().strip().split("\t")
    lgID, scaffoldID, posID = JasonUtils.get_headerIDs(header, ["lg","scaffolds","pos"])
    for line in MAP:
        data = line.strip().split("\t")
        lg, pos, scaffolds = data[lgID], data[posID], data[scaffoldID],
        scaffolds = scaffolds.split(',')
        for scaffold in scaffolds:
            master[scaffold] = lg + "\t" + pos
    MAP.close()
    print("\tFound",len(master),"target scaffolds")
    return master


def output_snps(scaffolds, snpfile, outfile):
    print("Filtering and outputting snps on target scaffolds")
    SNP = open(snpfile, "r")
    OUT = open(outfile, "w")
    OUT.writelines("snp\tscaffold\tlg\tpos\n")

    #Parse header
    header = SNP.readline().strip().split("\t")
    snpID, scaffoldID = JasonUtils.get_headerIDs(header, ["snp", "scaffold"])

    #Go through
    n=0
    total=0
    for line in SNP:
        n+=1
        if n % 50000 == 0:
            print("Processed",n,"snps")
        #if n > 100: break   #For debugging
        data = line.strip().split("\t")
        snp, scaffold = data[snpID], data[scaffoldID]
        if scaffold in scaffolds:
            total+=1
            OUT.writelines(snp + "\t" + scaffold + "\t" + scaffolds[scaffold] + "\n")
    SNP.close()
    OUT.close()
    print("Found",total,"of",n,"snps in target scaffolds")

if __name__ == "__main__":
    main()