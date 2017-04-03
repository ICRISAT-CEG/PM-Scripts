__author__ = 'jgw87'
"""
Take the mapped tags and determine their position in the original fasta files based off the SAM file
Output is set so can be processed by same script that takes the anchored SNPs
    -i Anchored tag file
    -m Map file (hapmap that tags were anchored to)
    -s SAM file of tag positions (if don't have TOPM for some reason)
    -t TOPM file (better than SAM)
    -o Output file
"""

import argparse
import Bio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

tagID, anchorID = 0, 3  # Columns in mapped tag file with needed info (tag sequence and anchor SNP name)
map_snpID, map_lgID, map_posID = 0, 2, 3    # Columns in anchor map (hapmap) with needed info
sam_tagID, sam_chrID, sam_posID = 9, 2, 3
topm_tagID, topm_chrID, topm_posID = 0, 3, 5

# sys.argv[1:] = ["-i", "SomData/5_som_anchored_tags_best.txt",
#                 "-s", "../03_SAM/refseq03_tags.sam",
#                 "-o", "999_temp_tags.txt", ]


def main():
    infile, mapfile, samfile, topmfile, outfile = parse_args()
    map = load_map(mapfile)  # Map = anchor map; used to extract LG and BIN
    key = make_binkey(infile, map)
    tagdata = None
    if topmfile:
        tagdata = extract_topm_position(key, topmfile)
    elif samfile:
        tagdata = extract_sam_position(key, samfile)    # Samfile provides the position data to create dummy SNP IDs to identify scaffolds
    else:
        print("Error: Need TOPM or SAM file to find tag positions")
        sys.exit(1)
    output_tagdata(tagdata, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-m", "--mapfile")
    parser.add_argument("-s", "--samfile", required=False)
    parser.add_argument("-t", "--topmfile", required=False)
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()

    if args.samfile is None and args.topmfile is None:
        print("Error: Must supply either a TOPM (preferred) or a SAM file")
        sys.exit(1)

    return args.infile, args.mapfile, args.samfile, args.topmfile, args.outfile


def load_map(mapfile):
    print("Loading anchor SNP map from",mapfile)
    map=dict()
    MAP = open(mapfile, "r")
    MAP.readline() #Clear header
    for line in MAP:
        metadata = line.split()[:12]
        snp, lg, pos = metadata[map_snpID], metadata[map_lgID], metadata[map_posID]
        if snp in map:
            print("WARNING!",snp,"already loaded into anchor map")
        map[snp] = lg + "\t" + pos
    print("\tLoaded",len(map),"anchor SNPs")
    return map


def make_binkey(infile, map):
    print("Reading in tags from", infile)
    binkey = dict()
    n = 0

    IN = open(infile, "r")
    IN.readline()
    for line in IN:
        n += 1
        if n % 100000 == 0:
            print("\tProcessed", n, "tags")
        #if n > 1000: break  #For debugging
        data = line.split()
        tag, anchor = data[tagID], data[anchorID]
        if tag in binkey:
            print("WARNING! Tag already loaded for", tag)
        if anchor not in map:
            print("WARNING!",anchor,"not found in reference map; skipping")
            continue
        binkey[tag] = map[anchor]   # Linkage group and bin from reference map

    IN.close()
    return (binkey)


def extract_topm_position(binkey, topmfile):
    print("Getting position data from TOPM file", topmfile)
    n = 0
    final = set()
    TOPM = open(topmfile, "r")
    header=TOPM.readline()
    unmapped=0
    for line in TOPM:
        n += 1
        if n % 1000000 == 0:
            print("\tProcessed", n, "lines (found",len(final),"tags)")
        #if n > 1e4: break   #For debugging
        data = line.split("\t")
        tag, chr, pos = data[topm_tagID], data[topm_chrID], data[topm_posID]
        #print(chr, ":",pos," - ", tag)
        if chr == "*":  # Skip unmapped tags
            #tag, inkey = find_tag_in_binkey(tag, binkey)
            if tag in binkey:
                unmapped+=1
            continue


        #if tag not in binkey:  # Check if tag is actually reverse-complemented
        #    tag = reverse_complement(tag)
        #if tag in binkey:
        if tag in binkey:
            chr = chr.lstrip("chr")
            dummy_snp = "dummy_S" + chr + "_" + pos
            tagdata = "\t".join((dummy_snp, tag, binkey[tag]))
            final.add(tagdata)
            #found+=1
    TOPM.close()
    print("Finished searching; final set has", len(final), "of", len(binkey), "input tags (",unmapped,"were not mapped )")
    return final


def extract_sam_position(binkey, samfile):
    print("Getting position data from SAM file", samfile)
    print("\tWarning: Tag matching for this is not as reliable as using a TOPM file")
    n = 0
    final = set()
    SAM = open(samfile, "r")
    unmapped=0
    for line in SAM:
        if line.startswith("@"):  # Skip header lines
            continue
        n += 1
        if n % 1000000 == 0:
            print("\tProcessed", n, "lines")
        #if n > 1e5: break   #For debugging
        data = line.split()
        tag, chr, pos = data[sam_tagID], data[sam_chrID], data[sam_posID]
        tag, tagfound = find_tag_in_binkey(tag, binkey)
        if chr == "*":  # Skip unmapped tags
            #tag, inkey = find_tag_in_binkey(tag, binkey)
            if tagfound:
                unmapped+=1
            continue
        
        
        #if tag not in binkey:  # Check if tag is actually reverse-complemented
        #    tag = reverse_complement(tag)
        #if tag in binkey:
        if tagfound:
            chr = chr.lstrip("chr")
            dummy_snp = "dummy_S" + chr + "_" + pos
            tagdata = "\t".join((dummy_snp, tag, binkey[tag]))
            final.add(tagdata)
            #found+=1
    SAM.close()
    print("Finished searching; final set has", len(final), "of", len(binkey), "input tags (",unmapped,"were not mapped )")
    return final

def find_tag_in_binkey(tag, binkey):
    revcomp = reverse_complement(tag)
    if tag in binkey:
        return tag, True
    if revcomp in binkey:
        return revcomp, True
    
    polyA="A"
    while len(tag) + len(polyA) <65:  # Add TASSEL's A padding and check again
        if tag + polyA in binkey:
            return tag + polyA, True
        if revcomp + polyA in binkey:
            return revcomp + polyA, True
        polyA+="A"

    return tag, False

def reverse_complement(seq):
    seq = seq[::-1]  # Obscure-looking code to reverse the string
    seq = seq.translate(str.maketrans("AGCT", "TCGA"))
    return seq


def output_tagdata(tagdata, outfile):
    """@type tagada: set"""
    print("Outputting results to", outfile)
    OUT = open(outfile, "w")
    OUT.writelines("anchor_snp\ttag\tchrom\tpos\n")
    for data in tagdata:
        if "\t" in data:    # Make sure data contains a tab, meaning the tag mapped
            OUT.writelines(data + "\n")
        else: print("NO TAG MAPPED FOR",data)
    OUT.close()


if __name__ == "__main__":
    main()