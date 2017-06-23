__author__ = 'jgw87'
"""
Take inputs of binned scaffolds from 6d_ and calculate the total lengths included in them
    -i Unambiguously binned scaffolds
    -a Ambiguous scaffolds
    -l File of scaffold names and lengths
    -o Output file
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

# sys.argv[1:] = ["-i", "SomData/6d_som_scaffold_bins.txt",
#                 "-a", "SomData/6d_som_scaffold_bins_ambiguous.txt",
#                 "-l", "/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v0.3/0_scaffold_lengths.txt",
#                 "-o", "SomData/6e_lg_totals.txt",]

def main():
    infile, ambigfile, lengthfile, outfile = parse_args()
    lengths = get_length_key(lengthfile)
    print("Adding counts")
    counts = add_good_counts(infile, lengths)
    counts = add_ambiguous_counts(ambigfile, lengths, counts)
    counts = add_unmapped_counts(lengths, counts)
    output_counts(counts, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-a", "--ambigfile")
    parser.add_argument("-l", "--lengthfile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    return args.infile, args.ambigfile, args.lengthfile, args.outfile


def get_length_key(lengthfile):
    print("Loading sequence length key")
    IN = open(lengthfile, "r")
    IN.readline()
    key=dict()
    for line in IN:
        seq, length = line.strip().split("\t")
        if seq in key:
            print("\tWARNING! More than one sequence named",seq)
        key[seq] = int(length)
    return key


def add_good_counts(infile, lengths):
    #Use header to identify appropriate columns
    IN = open(infile, "r")
    header = np.array(IN.readline().strip().split("\t"))
    lgID = np.where(header == "lg")[0]
    scaffoldID = np.where(header == "scaffold")[0]

    counts = dict()
    for line in IN:
        data = line.strip().split("\t")
        lg, scaffold = data[lgID], data[scaffoldID]
        if scaffold == "unknown": continue	#Rare, but happens
        if lg not in counts:
            counts[lg] = 0
        counts[lg] += int(lengths[scaffold])
        lengths[scaffold] = 0   #Set to zero to indicate that has already been added
    IN.close()
    return counts

def add_ambiguous_counts(ambigfile, lengths, counts):
    #Use header to identify appropriate columns
    IN = open(ambigfile, "r")
    header = np.array(IN.readline().strip().split("\t"))
    scaffoldID = np.where(header == "scaffold")[0]

    #Load all unique scaffolds
    scaffolds = set()
    for line in IN:
        data = line.strip().split("\t")
        if data[scaffoldID] != "unknown":
            scaffolds.add(data[scaffoldID])

    #Add up counts
    counts["ambiguous"] = 0
    for scaf in scaffolds:
        counts["ambiguous"] += lengths[scaf]
        lengths[scaf] = 0   #Set to zero to indicate that has already been added
    IN.close()
    return counts


def add_unmapped_counts(lengths, counts):
    count=0
    for scaffold in lengths:
        count+= lengths[scaffold]   # Already found ones are 0, so don't contribute
    counts["unknown"] = count
    return counts


def output_counts(counts, outfile):
    OUT = open(outfile, "w")
    OUT.writelines("lg\ttotal_nt\n")
    for lg in sorted(list(counts.keys())):
        OUT.writelines(str(lg) + "\t" + str(counts[lg]) + "\n")
    OUT.close()

if __name__ == "__main__":
    main()