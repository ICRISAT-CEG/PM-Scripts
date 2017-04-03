__author__ = 'jgw87'
"""
Take the snp and tag scaffold assignments and consolidate into a unified framework. Only does those parts that are
supplied in the options
    -m Combined genetic map file
    -s snp assignment file (optional)
    -t tag assignment file (optional)
    -r raw output file (optional)
    -o polished output file (optional)
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

map_lgID, map_binID = 2, 3  # Columns in reference map hapmap file with the LG and bin data
lgID, binID, scaffoldID, countID = 0, 1, 2, 3

# sys.argv[1:] = ["-m", "SomData/4k_som_anchored_map_reordered_combined.hmp.txt",
#                 "-s", "SomData/6a_som_snp_scaffold_assignments.txt",
#                 "-t", "SomData/6c_som_tag_scaffold_assignments.txt",
#                 "-r", "SomData/6d_scaffold_raw_counts.txt",
#                 "-o", "SomData/6d_scaffold_bins.txt"]


def main():
    mapfile, snpfile, tagfile, rawfile, outfile = parse_args()
    mapkey = make_map_key(mapfile)
    scaffolds = dict()
    if snpfile:
        scaffolds = tabulate_scaffold_assignments(scaffolds, mapkey, snpfile, "snps")
    if tagfile:
        scaffolds = tabulate_scaffold_assignments(scaffolds, mapkey, tagfile, "tags")
    if rawfile:
        output_raw_data(scaffolds, rawfile)
    if outfile:
        output_clean_data(scaffolds, outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mapfile")
    parser.add_argument("-o", "--outfile", required=False)
    parser.add_argument("-r", "--rawfile", required=False)
    parser.add_argument("-s", "--snpfile", required=False)
    parser.add_argument("-t", "--tagfile", required=False)
    args = parser.parse_args()
    return args.mapfile, args.snpfile, args.tagfile, args.rawfile, args.outfile


def make_map_key(mapfile):
    print("Creating anchor map")
    anchor_bins = dict()
    target_bins = dict()
    #Read in reference map
    MAP = open(mapfile, "r")
    MAP.readline()  #Clear header
    for line in MAP:
        data = line.split("\t")[:4]
        lg, bin = int(data[map_lgID]), int(data[map_binID])
        if line.startswith("S"):  # Anchor snps
            anchor_bins = add_bin(anchor_bins, lg, bin)
        elif line.startswith("target"):  #Other SNPs
            target_bins = add_bin(target_bins, lg, bin)
        else:
            print("WARNING! Unrecognized SNP name for", line[:20], "...; Should begin with 'S' or 'target'")
    print("Found", count_unique_bins(anchor_bins), "unique anchor bins and", count_unique_bins(target_bins),
          "unique target bins")
    MAP.close()

    #Consolidate into a dict with LG -> given_bin -> anchor_bin
    print("Consolidating down to just anchor bins")
    final_bins = dict() # Dictionary to map things over
    binvals = dict()    # Stores positions of anchor SNPs so can find nearest
    multiples=0
    for lg in anchor_bins:  # Anchor bins; just copy over
        if lg not in final_bins:
            final_bins[lg] = dict()
            binvals[lg] = list()
        for bin in anchor_bins[lg]:
            final_bins[lg][bin] = bin
            binvals[lg].append(bin)
    for lg in binvals:
        binvals[lg] = np.array(binvals[lg]) # Convert to ndarray for easy math
    for lg in target_bins:  # Target bins; match to nearest anchor bin
        if lg not in final_bins:
            print("Warning! Target linkage group",lg,"not among anchor bins")
        for bin in target_bins[lg]:
            differences = np.abs( bin - binvals[lg])
            nearest = np.min(differences)
            nearest_index = np.where(differences == nearest)[0] # Returns tuple, so extract data out of tuple
            if len(nearest_index) > 1:
                #print("WARNING!",len(nearest_index),"nearest anchors found for",lg,":",bin," - ",binvals[lg][nearest_index])
                multiples+=1
                nearest_index = nearest_index[0]    # Arbitrarily take the one with the lower index
            #print("Nearest bin to",lg,":",bin, "is",lg,":", binvals[lg][nearest_index])
            final_bins[lg][bin] = int(binvals[lg][nearest_index])
    print("\tFound",multiples,"instances of equal distance to bins; first one chosen at random")
    return final_bins


def add_bin(group, lg, bin):
    if lg not in group:
        group[lg] = set()
    group[lg].add(bin)
    return group


def count_unique_bins(group):
    total=0
    for lg in group:
        total += len(group[lg])
    return total


def tabulate_scaffold_assignments(scaffolds, mapkey, infile, name):
    print("Tabulating scaffolds in",infile,"for",name)
    IN = open(infile, "r")
    IN.readline()   # Clear header
    for line in IN:
        data=line.strip().split("\t")
        lg, bin, scaffold, count = int(float(data[lgID])), int(float(data[binID])), data[scaffoldID], int(data[countID])
        anchorbin = mapkey[lg][bin]
        #print(scaffold,lg,":",bin,"->",lg,":",anchorbin,"; count=",count)
        if scaffold not in scaffolds:
            scaffolds[scaffold] = dict()
        if lg not in scaffolds[scaffold]:
            scaffolds[scaffold][lg] = dict()
        if anchorbin not in scaffolds[scaffold][lg]:
            scaffolds[scaffold][lg][anchorbin] = dict()
        if name in scaffolds[scaffold][lg][anchorbin]:
            #print("WARNING!",scaffold,"already has",name,"count data for",lg,":",anchorbin)
            scaffolds[scaffold][lg][anchorbin][name] += count
        else:
            scaffolds[scaffold][lg][anchorbin][name] = count
    return scaffolds



def output_raw_data(scaffolds, rawfile):
    scaffoldlist = list()
    lglist = list()
    binlist = list()
    snpcount = list()
    tagcount = list()
    for scaffold in scaffolds:
        for lg in scaffolds[scaffold]:
            for bin in scaffolds[scaffold][lg]:
                scaffoldlist.append(scaffold)
                lglist.append(lg)
                binlist.append(bin)
                if "snps" in scaffolds[scaffold][lg][bin]:
                    snpcount.append(scaffolds[scaffold][lg][bin]["snps"])
                else:
                    snpcount.append(0)
                if "tags" in scaffolds[scaffold][lg][bin]:
                    tagcount.append(scaffolds[scaffold][lg][bin]["tags"])
                else:
                    tagcount.append(0)
    crude=pd.DataFrame({"lg": lglist, "bin":binlist, "scaffold":scaffoldlist, "snp_count":snpcount, "tag_count":tagcount})
    crude = crude.sort(columns=("scaffold", "lg","bin"))
    crude.to_csv(rawfile, sep="\t", index=False)


def output_clean_data(scaffolds, outfile):
    clean = {"scaffold":list(), "lg":list(), "bin":list(), "binval":list()}
    ambiguous = {"scaffold":list(), "lg":list(), "bin":list(), "snps":list(), "tags":list(), "binval":list()}
    master_bins=dict()
    for scaffold in scaffolds:
        lg = list(scaffolds[scaffold].keys())
        if len(lg) == 1:
            lg = lg[0]
            b = find_best_bin(scaffolds[scaffold][lg])
            clean["scaffold"].append(scaffold)
            clean["lg"].append(lg)
            clean["bin"].append(b)
            master_bins = record_master_bin(master_bins, lg, b)
            #Record bin

        else:
            mylg = find_best_lg(scaffolds[scaffold])
            if len(mylg) == 1:  # Clean assignment
                mylg = mylg[0]
                b = find_best_bin(scaffolds[scaffold][mylg])
                clean["scaffold"].append(scaffold)
                clean["lg"].append(mylg)
                clean["bin"].append(b)
                master_bins = record_master_bin(master_bins, mylg, b)
            else:
                for lg in mylg:
                    b = find_best_bin(scaffolds[scaffold][lg])
                    ambiguous["scaffold"].append(scaffold)
                    ambiguous["lg"].append(lg)
                    ambiguous["bin"].append(b)
                    snptotal, tagtotal = 0, 0
                    for mybin in scaffolds[scaffold][lg]:
                        if "snps" in scaffolds[scaffold][lg][mybin]:
                            snptotal += scaffolds[scaffold][lg][mybin]["snps"]
                        if "tags" in scaffolds[scaffold][lg][mybin]:
                            tagtotal += scaffolds[scaffold][lg][mybin]["tags"]
                    ambiguous["snps"].append(snptotal)
                    ambiguous["tags"].append(tagtotal)
                    master_bins = record_master_bin(master_bins, lg, b)

    #Get mapping key of bin values to more "reasonable" ones like 1.001, 1.002, etc
    binmap = make_map_of_bins(master_bins)

    #Write out
    output_to_dataframe(clean, binmap, outfile)
    output_to_dataframe(ambiguous, binmap, re.sub(string=outfile, pattern=".txt", repl="_ambiguous.txt"))


def output_to_dataframe(data, binmap, outfile):
    for lg, bin in zip(data["lg"], data["bin"]):
        data["binval"].append(binmap[lg][bin])
    df = pd.DataFrame()
    for col in sorted(list(data.keys())):
        df[col] = data[col]
    df = df.sort(columns=["binval","scaffold"])
    df.to_csv(outfile, sep="\t", index=False)

def record_master_bin(master, lg, bin):
    if lg not in master:
        master[lg] = set()
    master[lg].add(bin)
    return master


# Find best bin in a single LG; returns a weighted median based on tag counts (or on SNP counts if no tags)
def find_best_bin(lg):
    #Assemble lists of bins, weighted by number of snps or number of tags
    tagbins = list()
    snpbins = list()
    for bin in lg:
        #print(lg,":",bin)
        if "snps" in lg[bin] and lg[bin]["snps"] > 0:
            snpbins += [bin] * lg[bin]["snps"]
        if "tags" in lg[bin] and lg[bin]["tags"] > 0:
            tagbins += [bin] * lg[bin]["tags"]

    #Determine median bin, preferring tags if possible
    median = None
    if len(tagbins) > 0:
        median = np.median(tagbins)
    elif len(snpbins) > 0:
        median = np.median(snpbins)
    else:
        print("WARNING! No tag or snp bins found!")

    #Return the median value. If an average, indicates which ones it should go between
    return median

def find_best_lg(scaffold):
    totals=dict()   # Keeps track of LG-specific SNPs and tags
    snptotal, tagtotal = 0, 0   # Grand totals
    # Count up SNPs and tags
    for lg in scaffold:
        totals[lg] = dict()
        totals[lg]["snps"] = 0
        totals[lg]["tags"] = 0
        for bin in scaffold[lg]:
            if "snps" in scaffold[lg][bin]:
                totals[lg]["snps"] += scaffold[lg][bin]["snps"]
            if "tags" in scaffold[lg][bin]:
                totals[lg]["tags"] += scaffold[lg][bin]["tags"]
        snptotal += totals[lg]["snps"]
        tagtotal += totals[lg]["tags"]
    #print("Totals: snps",snptotal,"tags:",tagtotal)
    #print("grandtotals:",totals)

    # Determine if any have more than half the total Tags (preferred) or SNPs
    #print("Determining best LG")
    for lg in totals:
        #print("\t",lg,totals[lg])
        if tagtotal > 0 and totals[lg]["tags"] > tagtotal/2:
            #print("\tLG",lg,"with",totals[lg]["tags"],"of",tagtotal,"tags has been chosen")
            return [lg]
        elif snptotal > 0 and totals[lg]["snps"] > snptotal/2:
            #print("\tLG",lg,"with",totals[lg]["snps"],"of",snptotal,"snps has been chosen")
            return [lg]
    #print("\tNo LG has > half the tags;")
    return list(totals.keys())  # If none have more than half, return list of all keys (= ambiguous assignation)


def make_map_of_bins(master):
    key = dict()
    for lg in master:
        key[lg] = dict()
        bins = np.array(list(master[lg]))
        bins.sort()
        #print(bins)
        binvals = lg + 0.001 * (np.array(range(len(bins))) + 1)  # Make it 1.001, 1.002, etc
        for bin, val in zip(bins, binvals):
            key[lg][bin] = val
    return key

if __name__ == "__main__":
    main()