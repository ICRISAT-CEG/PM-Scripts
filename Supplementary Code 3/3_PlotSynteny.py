__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

debug = False


def main():
    args = parse_args()
    a, a_offsets = get_gene_locs(args.a)
    b, b_offsets = get_gene_locs(args.b)

    print("Plotting synteny to",args.outfile)
    fig = plt.figure(figsize=[10,10])
    ax=fig.add_subplot(111, title="Synteny", xlabel=args.a, ylabel=args.b)

    x, y = list(), list()
    for gene in a:
        if gene not in b: continue
        x.append(a[gene])
        y.append(b[gene])

    ax.plot(x, y, "b.")
    for o in a_offsets:
        ax.axvline(a_offsets[o])
    for o in b_offsets:
        ax.axhline(b_offsets[o])
    ax.invert_yaxis()

    fig.savefig(args.outfile, dpi=100)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--a", help="file of gene locations 1")
    parser.add_argument("-b", "--b", help="second file of gene locations")
    parser.add_argument("-o", "--outfile", help="Output graphic")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def get_gene_locs(infile):
    print("Loading gene locations from",infile)
    data = pd.read_csv(infile, sep='\t')
    genes, chroms, poses = [np.array(data[c]) for c in ['gene','chr','pos']]    # Split b/c numpy iterates much faster than pandas
    chroms = np.array([str(x) for x in chroms])

    # Get maximum for each chrom to make an offset key
    chromlist = np.unique(chroms)
    chromlist = sorted(chromlist[chromlist != 'nan'])
    chromstart = [0] * len(chromlist)
    for i in range(1, len(chromlist)):
        previous = chromlist[i-1]
        chromstart[i] = max(poses[chroms==previous]) + chromstart[i-1]
    offset = {c:s for c,s in zip(chromlist, chromstart)}

    # Now go through and get a position for each gene
    genekey = dict()
    for mygene, mychrom, mypos in zip(genes, chroms, poses):
        if mygene in genekey: print("\tWarning! Gene",mygene,"already loaded")
        if mychrom == 'nan': continue
        if debug and infile.find("2_existing_pm_genes.txt") >=0 and mychrom != 1: continue
        genekey[mygene] = mypos + offset[mychrom]
    print("Loaded positions for",len(genekey),"genes")

    return genekey, offset




if __name__ == '__main__': main()