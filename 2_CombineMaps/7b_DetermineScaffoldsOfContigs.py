__author__ = 'jgw87'
"""
Take a SAM alignment of primer pairs (with _fwd and _rev in names) and determine which scaffold the pairing is on.
Only outputs those where both primers map to the same scaffold
    -s SAM file of aligned primers (assumes the -k option in Bowtie was used, giving multiple alignments)
    -l Linker file hooking pseudochromosomes to scaffolds
    -o Output file with primer name, scaffold number, and distance between primers (mostly for spot-checking)
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

# sys.argv[1:] = ["-s", "1a_primers_unpaired.sam",
#                 "-l", "../0_pseudochromosome_scaffold_key.txt",
#                 "-o", "1b_primer_scaffold_assignments.txt" ]

nameID, chrID, posID = 0, 2, 3  # SAM file columns with needed information


def main():
    linkfile, samfile, outfile = parse_args()
    linker = pd.read_csv(linkfile, sep="\t")
    contigs = parse_sam_file(samfile)
    assign_scaffolds(contigs, linker, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--linkfile")
    parser.add_argument("-s", "--samfile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    return args.linkfile, args.samfile, args.outfile


def parse_sam_file(samfile):
    print("Parsing SAM alignment", samfile)
    SAM = open(samfile, "r")
    contigs = dict()
    n = 0
    for line in SAM:
        if line.startswith("@"): continue  #Skip header lines
        n += 1
        data = line.strip().split("\t")
        name, chr, pos = data[nameID], data[chrID], data[posID]
        name = re.sub(string=name, pattern="_.+", repl="")
        if chr == "*":  #Skip unmapped reads
            continue
        #chr = int(chr.lstrip("chr"))
        pos = int(pos)
        contigs = add_contig(contigs, name, chr, pos)
    print("Parsed", n, "lines of input for", len(contigs), "unique contigs")
    return contigs


#Add a new primer to the dictionary, checking that various structures are in place
def add_contig(contigs, name, chr, pos):
    if name not in contigs:
        contigs[name] = dict()
        contigs[name]["chr"] = list()
        contigs[name]["pos"] = list()
    contigs[name]["chr"].append(chr)
    contigs[name]["pos"].append(pos)
    return contigs


def assign_scaffolds(contigs, linker, outfile):
    OUT = open(outfile, "w")
    OUT.writelines("contig\tscaffold\trank\n")
    for contig in contigs:
        # fwd = contigs[contig]["fwd"]  # extract these, which are dicts of "chr" and "pos", each a list
        # rev = contigs[contig]["rev"]
        # # i_fwd, i_rev = find_best_pairing(fwd, rev)
        # #print("Returned values",i_fwd, i_rev)
        # if i_fwd is None or i_rev is None:  #Skip if didn't find a match
        #     continue
        chr= contigs[contig]["chr"]
        pos = contigs[contig]["pos"]
        nmatches = len(chr)
        for i in range(nmatches):
            scaffold = find_scaffold(linker, chr[i], pos[i])
            rank = str(i+1)   #Make rank assignment
            if nmatches > 1:
                rank += " of " + str(nmatches)
            OUT.writelines("\t".join((contig, scaffold, rank)) + "\n")
        # rev_scaffold = find_scaffold(linker, rev["chr"][i_rev], rev["pos"][i_rev])
    OUT.close()


# #Take the fwd and reverse primer data nd find the best pairing, defined as the one with the lowest total position value
# def find_best_pairing(fwd, rev):
#     min_score = 1e6  # Way beyond anything that should come up
#     max_dist=1000   #Largest distance permitted
#     placeholder = None
#     #print("Fwd:", fwd_chr)
#     #print("Rev:", rev_chr)
#     for i in range(len(fwd["chr"])):
#         for j in range(len(rev["chr"])):
#             if fwd["chr"][i] == rev["chr"][j]:
#                 score = i + j
#                 distance = abs(fwd["pos"][i] - rev["pos"][j])
#                 if score < min_score and distance < max_dist:
#                     min_score = score
#                     placeholder = (i, j)
#     if placeholder is not None:
#         return placeholder
#     else:
#         return None, None


#Take a pseudochromosome and position and get the scaffold number from a Pandas dataframe with linker data
def find_scaffold(linker, chrom, pos):
    target = np.where((linker["chrom"] == chrom) & (linker["start"] <= pos) & (linker["end"] >= pos))[0]
    if len(target) != 1:
        print("\t\tWARNING! Number of targets ==", len(target), "for chrom:", chrom, "pos:", pos)
    if len(target) == 0:
        return "unknown"
    #print("Finding:", chrom, ":", pos, "matches to scaffold", linker["scaffold"].iloc[target])  #for debugging
    #print(linker["scaffold"].iloc[target].iloc[0])
    return linker["scaffold"].iloc[target].iloc[0]


if __name__ == "__main__":
    main()