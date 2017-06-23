__author__ = 'jgw87'
"""
Take a hapmap output from my mapping pipeline and the linker file from the original concatenation and output the
original scaffolds in which bin they belong to
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys
import time

# sys.argv[1:] = ["-i", "SomData/4k_som_anchored_map_reordered_combined.hmp.txt",
#                 "-c", "chr",
#                 "-l", "/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/0_pseudochromosome_scaffold_key.txt",
#                 "-o", "999_test_scaffolds.txt"]

  # Columns in hapmap with needed data


def main():
    print("Assigning scaffolds to linkage positions")
    infile, linkfile, outfile, chrom_prefix = parse_args()
    linker = pd.read_csv(linkfile, sep="\t")
    map = map_scaffolds(infile, linker, chrom_prefix, outfile)
    #output_map(map, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--chrom_prefix")
    parser.add_argument("-i", "--infile")
    parser.add_argument("-l", "--linkfile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    return args.infile, args.linkfile, args.outfile, args.chrom_prefix


def map_scaffolds(infile, linker, chrom_prefix, outfile):
    """
    Take a hapmap and linker and return a nested dictionary of linkage groups -> bins -> scaffolds
    @param linker: pd.DataFrame
    """
    print("Parsing SNPs from", infile)
    mymap = dict()
    IN = open(infile, "r")
    OUT=open(outfile, "w")
    
    # Clear out header
    header=IN.readline()
    OUT.write(header)
    header = header.strip().split('\t')
    lgID, binID = header.index('chrom'), header.index('pos')
    nameID = header.index('rs#') if 'rs#' in header else header.index('tag')    # For hapmaps or tags

    # Process lines
    n = 0
    linker["chr_int"] = list(map( lambda x: int(x.lstrip("chr")), linker["chrom"]))
    linker = linker.sort(columns = ["chr_int", "start"])
    search_chr, search_start, search_end = np.array(linker["chr_int"]), np.array(linker["start"]), np.array(linker["end"])
    bigchrom = search_chr * 1e9 + search_start
    found, unmapped = 0, 0
    for line in IN:
        n += 1
        if n % 10000 == 0:
            print("\tProcessed", n, "sites")
        #if n > 10000: break  # For debugging
        #if n < 35000: continue  # For debugging
        data=line.strip().split("\t")
        metadata = data[:12]
        if metadata[lgID] == "*":   # Skip tags that couldn't be mapped
            continue
        name = metadata[nameID]
        lg = int(metadata[lgID])
        bin = float(metadata[binID])
        ref_chrom, ref_pos = extract_chrom_pos(name, chrom_prefix)

        target = np.searchsorted(bigchrom, int(ref_chrom.lstrip("chr")) * 1e9 + ref_pos, side="left") -1
        if search_chr[target] == int(ref_chrom.lstrip("chr")) and search_start[target] <= ref_pos and search_end[target] >= ref_pos:
            scaffold = linker["scaffold"].iloc[target]
            mypos = ref_pos - linker["start"].iloc[target]
            found+=1
        else:
            scaffold = "unknown"
            mypos = -1
            if ref_chrom=='chr0': unmapped+=1

        if scaffold != 'unknown':
            data[nameID] = scaffold + "_" + str(mypos)
            #data[lgID]=scaffold
            OUT.write("\t".join(data) + "\n")
        #mymap = add_scaffold(mymap, lg, bin, scaffold)
    IN.close()
    print("Matched",found,"of",n,"scaffolds (",unmapped,"of which were not mapped, so total of",found + unmapped,")")
    return mymap


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



def add_scaffold(map, lg, bin, scaffold):
    #Create levels if needed
    if lg not in map:
        map[lg] = dict()
    if bin not in map[lg]:
        map[lg][bin] = dict()
    #Increment scaffold if exists, otherwise set to 1
    if scaffold in map[lg][bin]:
        map[lg][bin][scaffold] += 1
    else:
        map[lg][bin][scaffold] = 1
    return map


def output_map(map, outfile):
    out = open(outfile, "w")
    out.writelines("linkage_group\tbin\tscaffold\tcount\n")
    for lg in sorted(map.keys()):
        for bin in sorted(map[lg].keys()):
            for scaffold in sorted(map[lg][bin]):
                out.writelines("\t".join((str(lg), str(bin), scaffold, str(map[lg][bin][scaffold]))) + "\n")
    out.close()


if __name__ == "__main__":
    main()