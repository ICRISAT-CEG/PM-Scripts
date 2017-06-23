__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
from scipy.stats import mode

debug = False


def main():
    args = parse_args()
    ref_groups = load_scaffolds_majority(args.allfiles, lg_key='LG', scaffold_key='Scaffold ID', sep=',')
    # counts = match_tags(ref_groups, target_groups)
    # if args.tablefile:
    #     print("Outputting matched LG table to",args.tablefile)
    #     counts.to_csv(args.tablefile, sep='\t')
    # rename_key = get_equivalencies(counts)
    update_file(args.infile, args.outfile, ref_groups)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-a", "--allfiles", nargs="*", help="All ALLMAPS-formatted maps to be included to make scaffolds"
                                                            "map to single LGs")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_scaffolds_majority(infiles, lg_key="lg", scaffold_key="scaffold", sep='\t'):
    print("Loading scaffolds from",len(infiles)," input files and collapsing to best linkage group per scaffold")
    scaffolds, groups =dict(), dict()
    for infile in infiles:
        print("\tLoeading from",infile)
        IN = open(infile, 'r')
        header=IN.readline().strip().split(sep)
        lgID, scaffoldID = header.index(lg_key), header.index(scaffold_key)
        n=0
        for line in IN:
            data=line.strip().split(sep)
            myLG, myscaffold = data[lgID], data[scaffoldID]
            if myscaffold not in scaffolds:
                scaffolds[myscaffold] = list()
            if myLG not in groups:  # Not needed now, but for later
                groups[myLG] = set()
            scaffolds[myscaffold].append(myLG)
            n+=1
        print("\t\tLoaded",n,"scaffolds")

    #Collapse by most abundant linkage group for each scaffold
    for myscaffold in scaffolds:
        most_LG, most_count = mode(scaffolds[myscaffold])
        most_LG = most_LG[0]
        # print(myscaffold,":",scaffolds[myscaffold])
        # print(most_LG, most_count, len(scaffolds[myscaffold]))
        if most_count >= len(scaffolds[myscaffold])/2:  # Only add scaffold if the best one is >= half the total
            groups[most_LG].add(myscaffold)
        # else: print(myscaffold)

    # for c in sorted(groups.keys()):
    #     print("\t\tLG",c,":",len(groups[c]))
    return groups


# def match_tags(ref, uneak):
#     print("Matching tags")
#     counts = np.zeros(shape=(len(ref), len(uneak)), dtype=int)
#     ref_chroms = sorted(ref.keys())
#     uneak_chroms = sorted(uneak.keys())
#     for r in ref_chroms:
#         row = ref_chroms.index(r)
#         for u in uneak_chroms:
#             intersect = ref[r] & uneak[u]
#             col = uneak_chroms.index(u)
#             counts[row, col] = len(intersect)
#     counts = pd.DataFrame(counts, index=ref_chroms, columns=uneak_chroms)
#     print(counts)
#     return(counts)
#
# def get_equivalencies(counts):
#     print("Determining equivalencies")
#     key=dict()
#     for old in counts.columns:
#         ismax = np.argmax(np.array(counts[old]))
#         new = counts.index[ismax]
#         key[old]=new
#     print("\tEquivalency key is:")
#     for k in sorted(key.keys()):
#         print("\t\t",k,"->",key[k])
#     return key

def update_file(infile, outfile, groups):
    print("Removing scaffolds from minor linkage groups")
    IN = open(infile, 'r')
    OUT = open(outfile, 'w')
    header=IN.readline().strip().split(',')
    lgID, scaffoldID = header.index("LG"), header.index("Scaffold ID")
    OUT.write(','.join(header) + '\n')
    updated, skipped = 0, 0
    for line in IN:
        data=line.strip().split(',')
        myLG, myscaffold = data[lgID], data[scaffoldID]
        if myscaffold in groups[myLG]:
            updated +=1
            OUT.write(line)
        else:
            skipped +=1
    print("\tUpdated",updated,"total lines and removed",skipped,"total lines")


if __name__ == '__main__': main()