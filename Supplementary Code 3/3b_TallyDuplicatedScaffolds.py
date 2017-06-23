__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
import sys

debug = False


def main():
    args = parse_args()
    tally = load_maps(args.infiles)
    duplicates = find_duplicates(tally)
    if len(duplicates) == 0:
        print("No duplciated scaffolds found; exiting")
        sys.exit(0)
    print_duplicates(duplicates, args.outfile)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="ALLMAPS input files")
    parser.add_argument("-o", "--outfile", help="Output file of duplicated scaffolds and their counts")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def load_maps(infiles):
    print("Loading maps from",len(infiles),"input files")
    tally=dict()
    for infile in infiles:
        print("\tLoading",infile)
        data = pd.read_csv(infile)
        for scaffold, lg in zip(np.array(data['Scaffold ID']), np.array(data['LG'])):
            # Set up nested dictionary to keep track of how often a scaffold is in a given input file
            if scaffold not in tally:
                tally[scaffold]=dict()
            if lg not in tally[scaffold]:
                tally[scaffold][lg] = dict()
            if infile not in tally[scaffold][lg]:
                tally[scaffold][lg][infile] = 0
            tally[scaffold][lg][infile] +=1
    # print(tally)
    return(tally)

def find_duplicates(tally):
    print("Finding duplicate scaffolds on separate linkage groups")
    dups = dict()
    for scaffold in tally:
        if len(tally[scaffold]) > 1:
            dups[scaffold] = tally[scaffold]
    print("\tFound",len(dups),"scaffolds on different linkage groups:",sorted(dups.keys()))
    return dups


def print_duplicates(dups, outfile):
    print("Outputting duplicate tallies to",outfile)
    OUT = open(outfile, 'w')
    OUT.write("scaffold\tlocations\n")
    for scaffold in sorted(dups.keys()):
        OUT.write(scaffold)
        for lg in dups[scaffold]:
            OUT.write("\tLG" + str(lg) + ":")
            for file in dups[scaffold][lg]:
                OUT.write(file + "(x" + str(dups[scaffold][lg][file]) + ")")
        OUT.write('\n')
    OUT.close()


if __name__ == '__main__': main()