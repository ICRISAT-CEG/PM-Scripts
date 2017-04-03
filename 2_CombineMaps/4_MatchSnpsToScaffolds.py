__author__ = 'jgw87'
"""
Take a hapmap file and determine which scaffold each SNP belongs on
    -i Input hapmap file
    -c Chromosome prefix
    -l Linker file
    -o Output file
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

#sys.argv[1:] = ["-i", "../06_HapMap/refseq_unfiltered.hmp.txt",
#                "-c", "chr",
#                "-l", "/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v0.3/0_pseudochromosome_scaffold_key.txt",
#                "-o", "999_test_scaffolds.txt"]

nameID, lgID, binID = 0, 2, 3  # Columns in hapmap with needed data


def main():
    print("Matching SNPs to scaffolds")
    infile, linkfile, outfile, chrom_prefix = parse_args()
    linker = pd.read_csv(linkfile, sep="\t")
    map = map_scaffolds(infile, linker, chrom_prefix)
    output_map(map, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--chrom_prefix")
    parser.add_argument("-i", "--infile")
    parser.add_argument("-l", "--linkfile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    return args.infile, args.linkfile, args.outfile, args.chrom_prefix


def map_scaffolds(infile, linker, chrom_prefix):
    """
    Take a hapmap and linker and return a nested dictionary of linkage groups -> bins -> scaffolds
    @param linker: pd.DataFrame
    """
    print("Parsing SNPs from", infile)
    print("NOTE: NEED TO RECODE THE SEARCH IN NUMPY TO GO 40-50x FASTER"); map = list()
    IN = open(infile, "r")
    IN.readline()  # Clear out header
    n = 0
    for line in IN:
        n += 1
        if n % 1000 == 0:
            print("\tProcessed", n, "sites")
        #if n > 1000: break  # For debugging
        #if n < 35000: continue  # For debugging
        metadata = line.strip().split("\t")[:12]
        name = metadata[nameID]
        chr = int(metadata[lgID])
        pos = int(metadata[binID])
        ref_chrom, ref_pos = extract_chrom_pos(name, chrom_prefix)
        scaffold = find_scaffold(linker, ref_chrom, ref_pos)
        final = "\t".join([name,scaffold])
        map.append(final + "\n")
    IN.close()
    return map


def extract_chrom_pos(sitename, chrom_prefix):  # Turn a site name back into a chromosome location
    regex = re.search(string=sitename, pattern="S(\d+)_(\d+)")
    chrom = chrom_prefix + regex.group(1)
    pos = int(regex.group(2))
    #print("Extraction:", sitename, "becomes", chrom, ":", pos)  #for debugging
    return chrom, pos


def find_scaffold(linker, chrom, pos):
    target = np.where((linker["chrom"] == chrom) & (linker["start"] <= pos) & (linker["end"] >= pos))[0]
    if len(target) != 1:
        print("\t\tWARNING! Number of targets ==", len(target), "for chrom:", chrom, "pos:", pos)
    if len(target) == 0:
        return "unknown"
    #print("Finding:", chrom, ":", pos, "matches to scaffold", linker["scaffold"].iloc[target])  #for debugging
    #print(linker["scaffold"].iloc[target].iloc[0])
    return linker["scaffold"].iloc[target].iloc[0]

#
# def add_scaffold(map, lg, bin, scaffold):
#     #Create levels if needed
#     if lg not in map:
#         map[lg] = dict()
#     if bin not in map[lg]:
#         map[lg][bin] = dict()
#     #Increment scaffold if exists, otherwise set to 1
#     if scaffold in map[lg][bin]:
#         map[lg][bin][scaffold] += 1
#     else:
#         map[lg][bin][scaffold] = 1
#     return map


def output_map(map, outfile):
    out = open(outfile, "w")
    out.writelines("snp\tscaffold\n")
    out.writelines(map)
    out.close()


if __name__ == "__main__":
    main()