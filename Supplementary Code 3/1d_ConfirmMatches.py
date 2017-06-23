__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import mode

debug = False


def main():
    args = parse_args()
    if args.majority_rule:
        target_groups = [load_scaffolds_majority(infile) for infile in args.infiles]
    else:
        target_groups = [load_scaffolds(infile) for infile in args.infiles]


    # Set up plot
    nrow=ncol=len(target_groups)
    fig =plt.figure(figsize=(3 * ncol, 3 * nrow))
    grid=gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.5, wspace=0.5)

    for i in range(len(target_groups)):
        # ax = fig.add_subplot(grid[i,i], title=args.infiles[i])
        for j in range(i,len(target_groups)):
            if j >= len(target_groups): continue    # Safety if goes over
            print("Testing",args.infiles[i],"and",args.infiles[j])

            # Check that things are right
            matches = match_tags(target_groups[i], target_groups[j])
            equivalencies = get_equivalencies(matches)
            for e in equivalencies:
                if e != equivalencies[e]:
                    print("\tWarning! Linkage groups do not match!",e,"->",equivalencies[e])
                    # print(matches)

            # Plot
            title = args.infiles[i] if i==j else ""
            ax = fig.add_subplot(grid[i,j], title=title)
            ax.pcolor(np.array(matches))



    # counts = match_tags(target_groups)
    fig.savefig(args.outfile, dpi=100)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-m", "--majority-rule", default=False, action="store_true",
                        help="Whether a scaffold will be assigned to only the linkage group with the most instances of it")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_scaffolds(infile, lg_key="LG", scaffold_key="Scaffold ID", sep=','):
    print("Loading scaffolds from",infile)
    IN = open(infile, 'r')
    header=IN.readline().strip().split(sep)
    lgID, scaffoldID = header.index(lg_key), header.index(scaffold_key)
    groups=dict()
    n=0
    for line in IN:
        data=line.strip().split(sep)
        myLG, myscaffold = data[lgID], data[scaffoldID]
        if myLG not in groups:
            groups[myLG] = set()
        groups[myLG].add(myscaffold)
        n+=1

    print("\tLoaded",n,"total tags")
    # for c in sorted(groups.keys()):
    #     print("\t\tLG",c,":",len(groups[c]))
    return groups


def load_scaffolds_majority(infile, lg_key="LG", scaffold_key="Scaffold ID", sep=','):
    print("Loading scaffolds from",infile,"and collapsing by majority rule")
    IN = open(infile, 'r')
    header=IN.readline().strip().split(sep)
    lgID, scaffoldID = header.index(lg_key), header.index(scaffold_key)
    scaffolds, groups =dict(), dict()
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

    #Collapse by most abundant linkage group for each scaffold
    for myscaffold in scaffolds:
        most_LG, most_count = mode(scaffolds[myscaffold])
        most_LG = most_LG[0]
        if most_count >= len(scaffolds[myscaffold])/2:  # Only add scaffold if the best one is >= half the total
            groups[most_LG].add(myscaffold)

    print("\tLoaded",n,"total tags")
    # for c in sorted(groups.keys()):
    #     print("\t\tLG",c,":",len(groups[c]))
    return groups


def match_tags(a, b):
    print("Matching tags")
    counts = np.zeros(shape=(len(a), len(b)), dtype=int)
    ref_chroms = sorted(a.keys())
    uneak_chroms = sorted(b.keys())
    for r in ref_chroms:
        row = ref_chroms.index(r)
        for u in uneak_chroms:
            intersect = a[r] & b[u]
            col = uneak_chroms.index(u)
            counts[row, col] = len(intersect)
    counts = pd.DataFrame(counts, index=ref_chroms, columns=uneak_chroms)
    return(counts)

def get_equivalencies(counts):
    print("Determining equivalencies")
    key=dict()
    for old in counts.columns:
        ismax = np.argmax(np.array(counts[old]))
        new = counts.index[ismax]
        key[old]=new
    # print("\tEquivalency key is:")
    # for k in sorted(key.keys()):
    #     print("\t\t",k,"->",key[k])
    return key




if __name__ == '__main__': main()