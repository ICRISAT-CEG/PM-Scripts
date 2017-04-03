__author__ = 'jgw87'
"""
Take an existing scaffold map and add all the extended scaffolds in another map to it (that aren't already on it). Also
have to provide the core map to know which are the anchors.
closest core scaffold
    -e Som's extended map
    -c Core map (mostly used for anchor scaffold names)
    -o Output map
    -s The source map to be added to
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

sys.path.append('/home/jgw87/Software/Python')
import JasonUtils

# sys.argv[1:] = ["-c","5e_core_map_scaffolds.txt",
#                 "-s","5e_core_map_scaffolds.txt",
#                 "-e","../MakeNewMaps/SomData/6d_som_scaffold_bins.txt",
#                 "-o","6c_core_and_some_map_combined.txt" ]


def main():
    extendfile, corefile, sourcefile, outfile = parse_args()
    source, extended, anchors = load_maps(sourcefile, extendfile, corefile)
    finalmap = combine_maps(source, extended, anchors)
    finalmap.to_csv(outfile, sep="\t", index=True, index_label="scaffold")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--extendfile")
    parser.add_argument("-c", "--corefile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-s", "--sourcefile")
    args = parser.parse_args()
    return args.extendfile, args.corefile, args.sourcefile, args.outfile


def load_maps(sourcefile, extendfile, corefile):
    source = pd.read_csv(sourcefile, sep="\t", index_col=None)
    source = source.set_index("scaffold")

    extended = pd.read_csv(extendfile, sep="\t", index_col=None)
    extended = extended.set_index("scaffold")

    core = pd.read_csv(corefile, sep="\t", index_col=None)
    anchors = np.array(core["scaffold"], dtype=object)
    anchors = anchors[anchors != "unknown"]

    return source, extended, anchors



def combine_maps(source, extended, anchors):
    #Set up data frames
    tomap = find_scaffolds_to_map(source, extended)
    print("Tomap=",tomap)
    newmap = pd.DataFrame(data={"lg":0, "pos":0, "anchor":"NONE"}, index = tomap)
    if "anchor" not in source.columns:
        print("Note: 'anchor' not in source columns; adding with 'n/a'")
        source["anchor"] = "n/a"
    anchormap = source.loc[anchors,:]

    #Go through and match each new scaffold to its closest anchor
    print("Anchoring scaffolds")
    n=0
    for scaffold in newmap.index:
        #print("Anchoring scaffold",n,":",scaffold,"corresponding to index value",newmap.index[n])
        n+=1
        #if n <250: continue #For debugging
        if scaffold in source.index:
            print("Warning!",scaffold,"is already in the source map")
        if n % 500 == 0:
            print("\tProcessed",n,"scaffolds")
        #if n > 100: break   #For debugging
        anchor, lg, pos = determine_anchor(scaffold, extended, anchormap)
        #print("\tLocus:",newmap.loc[scaffold, "lg"],"to be",lg)
        newmap.loc[scaffold, "lg"] = lg
        newmap.loc[scaffold, "pos"] = pos
        newmap.loc[scaffold, "anchor"] = anchor

    # #Clean up
    finalmap = pd.concat([source, newmap])
    finalmap = finalmap.sort(columns=["lg","pos"])
    return finalmap

def find_scaffolds_to_map(source, extended):
    print("Determining which scaffolds to anchor")
    in_source = np.array(source.index, dtype=object)
    in_extended = np.array(extended.index, dtype=object)
    to_ignore = np.in1d(in_extended, in_source) #Find the ones in the extended map already been included in source
    print("\t Source has",len(in_source),"scaffolds")
    print("\t Extended has",len(in_extended),"scaffolds")
    print("\t",np.sum(to_ignore),"are in both and",np.sum(~to_ignore),"need to be added")
    return in_extended[~to_ignore]  #Return those not in the map


def determine_anchor(scaffold, extended, anchormap):

    #Check if is in the anchors (shouldn't be, but best to be careful)
    if scaffold in anchormap.index:
        print("Warning!",scaffold,"is an anchor SNP and shouldn't be mapped")
        return scaffold, int(anchormap.loc[scaffold, "lg"]), anchormap.loc[scaffold, "pos"]

    #Find nearest anchor on the same linkage group
    ##Subset out just the linkage group the scaffold mapped to
    oldlg = int(extended.loc[scaffold, "lg"])
    subset = extended.loc[extended["lg"] == oldlg,:]
    ##Determine which one is the closest and return its name, lg, and bin
    best_anchor = find_closest_anchor(scaffold, subset, anchormap)
    best_lg = anchormap.loc[best_anchor, "lg"]
    best_bin = anchormap.loc[best_anchor, "pos"]

    if isinstance(best_lg, pd.Series):
        best_lg = best_lg.iloc[0]
    if isinstance(best_bin, pd.Series):
        best_bin = best_bin.iloc[0]    
    return best_anchor, best_lg, best_bin

def find_closest_anchor(scaffold, subset, anchormap):
    mybin = subset.loc[scaffold, "bin"]
    best_dist=1e9
    best_anchor="NONE"
    #print(scaffold,"at",mybin)
    for anchor in anchormap.index:
        if anchor not in subset.index:  #Skip anchors that didn't make it into the source map (shouldn't be any)
            continue
        #print(subset.loc[anchor, :])
        #print(subset.loc[anchor, "bin"])
        dist = abs(float(subset.loc[anchor, "bin"]) - mybin)
        #print("\tAnchor",anchor,"has dist",dist,"and is at pos", subset.loc[anchor, "bin"])
        if dist < best_dist:
            #print("\t\tSet to best!")
            best_dist = dist
            best_anchor = anchor
    #print("\tBest anchor is",best_anchor,"with dist",best_dist)
    return best_anchor



if __name__ == "__main__":
    main()