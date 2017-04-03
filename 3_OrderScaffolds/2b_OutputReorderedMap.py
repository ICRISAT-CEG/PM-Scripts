__author__ = 'jgw87'
"""
Take output from step 2a_ and the original bin maps and create a new map
    -i Input ordering from Traveling Salesman algorithm
    -m Original map file
    -o Output file
"""

import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re
import sys

sys.path.append('/home/jgw87/Software/Python')
import JasonUtils

# sys.argv[1:] = ["-i", "2a_map_reordered_from_outer100_tsp.txt",
#                 "-m", "/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet/20140527_AlignToRefseqV1.1/ConsolidateMaps/7e_expanded_map2_som_841_renumbered.txt",
#                 "-o", "2b_new_map_from_outer100.txt"]


def main():
    infile, mapfile, outfile = parse_args()
    oldmap = load_old_map(mapfile)
    newmap = load_new_map(infile)
    newmap = extract_new_map(newmap)
    output_new_map(newmap, oldmap, outfile)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-m", "--mapfile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    return args.infile, args.mapfile, args.outfile


def load_old_map(mapfile):
    print("Loading previous map")
    MAP = open(mapfile, "r")
    header=MAP.readline().strip()
    scaffoldID, lgID ,posID = JasonUtils.get_headerIDs(header.strip().split(), ["scaffold","lg","pos"])

    oldmap=dict()
    n=0
    for line in MAP:
        data=line.strip().split("\t");
        #print(scaffoldID, lgID, posID, data)
        scaffold, lg, pos = data[scaffoldID], int(data[lgID]), float(data[posID])
        if lg not in oldmap:
            oldmap[lg] = dict()
        if pos not in oldmap[lg]:
            oldmap[lg][pos] = dict()
        oldmap[lg][pos][scaffold] = 1 # Flag to keep track of; used to check that all loaded
        n+=1
    MAP.close()
    print("\tLoaded",n,"scaffolds")
    return oldmap


def load_new_map(infile):
    print("Loading new map order")
    IN = open(infile, "r")
    header = IN.readline().strip().split("\t")
    scaffoldID, lgID, posID = JasonUtils.get_headerIDs(header, ["scaffolds", "lg", "pos"])

    newmap=dict()
    for line in IN:
        data=line.strip().split("\t")
        scaffolds, lg, pos = data[scaffoldID], int(data[lgID]), float(data[posID])
        if lg not in newmap:
            newmap[lg] = dict()
        newmap[lg][pos] = scaffolds.split(",")
    IN.close()
    return newmap


#Split each bin into scaffolds and get their orientation
def extract_new_map(newmap):
    print("Reformatting new map")
    for lg in newmap:
        for pos in newmap[lg]:
            mybin = newmap[lg][pos] #Extract list
            oldname, oldori = "null", "null"
            tempnames, tempdirs = [] ,[]
            for s in mybin:
                myname, myori = get_name_and_orientation(s)
                #print(s,"to",myname, myori)
                dir="?"
                #Determine orientation
                if myname == oldname:
                    if oldori == "left" and myori == "right":
                        tempdirs[-1] = "fwd"
                    elif oldori == "right" and myori == "left":
                        tempdirs[-1] ="rev"
                    else:
                        print("Error! Unable to determine orientation from",oldname, oldori,"to",myname, myori)
                if myname != oldname:
                    if myname not in tempnames: #Make sure not already added, in case scaffold got split (rare, but happens)
                        tempnames.append(myname)
                        tempdirs.append("?")
                oldname=myname
                oldori = myori
            newmap[lg][pos] = {'scaffold':tempnames, 'dir':tempdirs}

    #Check that everything got converted
    classes=dict()
    for lg in newmap:
        for pos in newmap[lg]:
            if newmap[lg][pos] is list():
                print("Error: ",lg,pos,"is still a list, not a dict()")
            else:
                mytype=type(newmap[lg][pos])
                if mytype not in classes:
                    classes[mytype]=0
                classes[mytype] +=1
    print("\tCheck for conversion (should only be dicts):\n\t",classes)
    return newmap


def output_new_map(newmap, oldmap, outfile):
    print("Writing new map to",outfile)
    OUT = open(outfile, "w")
    OUT.writelines("\t".join(["lg","bin","scaffold","orientation"]) + "\n")

    found, notfound = 0, 0
    for lg in sorted(newmap.keys()):
        for pos in sorted(newmap[lg].keys()):
            scaffolds = newmap[lg][pos]["scaffold"]
            dirs = newmap[lg][pos]["dir"]
            for i in range(len(scaffolds)):
                OUT.writelines("\t".join([str(lg),str(pos), scaffolds[i], dirs[i]]) + "\n")
                oldmap[lg][pos][scaffolds[i]] -= 1   #Subtract one (should set to 0, but lets see duplicates for <0)
                found+=1
            for oldscaff in oldmap[lg][pos]:
                mycount =oldmap[lg][pos][oldscaff]
                if mycount == 1:
                    notfound+=1
                    oldmap[lg][pos][oldscaff] -= 1
                    OUT.writelines("\t".join([str(lg),str(pos), oldscaff, "?"]) + "\n")
                if mycount <0:
                    print("\tWarning!",oldscaff,"found multiple times with a score of",mycount)
                    print("\t\t",scaffolds)
    print("Found",found,"scaffolds in refined order;",notfound,"were not found and added at the end of each bin")
    OUT.close()


def get_name_and_orientation(scaffold):
    #Add "scaffold" back into name if not already there (and not a contig)
    if scaffold.startswith("C") or scaffold.startswith("scaffold") or scaffold.startswith("unknown"):
        pass
    else:
        scaffold = "scaffold" + scaffold
    return scaffold.split("|")


if __name__ == "__main__":
    main()