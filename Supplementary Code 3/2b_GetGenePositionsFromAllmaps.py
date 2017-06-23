__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd

debug = False


def main():
    args = parse_args()
    genes = make_gene_key(args.syntenymap)
    make_gff(args.infile, genes, args.outfile, args.chromname)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="ALLMAPS .agp file output (specifically, the .chr.agp for just mapped scaffolds)")
    parser.add_argument("-s", "--syntenymap", help="Comma-separated synteny map (to get scaffold positions of genes")
    parser.add_argument("-o", "--outfile", help="Output file in GFF-like format (chrom, gene, position) ")
    parser.add_argument("-c", "--chromname", help="Chromosome name (if want to replace)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def make_gene_key(infile):
    print("Making gene key from",infile)
    data = pd.read_csv(infile)
    genes, scaffolds, positions = [np.array(data[c]) for c in ['PM gene name', 'Scaffold ID', 'scaffold position']]

    key = dict()
    for mygene, myscaffold, mypos in zip(genes, scaffolds, positions):
        if myscaffold not in key:
            key[myscaffold]=dict()
        if mygene in key[myscaffold]:
            print("Error! Gene",mygene,"already loaded on scaffold",myscaffold)
        key[myscaffold][mygene] = mypos

    return(key)

def make_gff(infile, genes, outfile, chromname=None):
    print("Making a GFF file from map construction in",infile)
    IN = open(infile)
    OUT = open(outfile, "w")
    OUT.write("\t".join(["gene","chr","pos"]) + "\n")

    chrID, startID, endID, typeID, nameID, oriID = [None] * 6
    for line in IN:
        # parse header
        if line.startswith("# FIELDS"):
            line=line.replace("# FIELDS:", "")
            line=line.replace(" ", "")
            header = line.strip().split(',')
            chrID, startID, endID, typeID, nameID, oriID = \
                [header.index(c) for c in ["object", "object_beg", "object_end", "component_type", "component_id/gap_length", "orientation/linkage_evidence"]]
            continue
        elif line.startswith("#"): continue   # Skip comments

        # Parse actual input lines
        data=line.strip().split()
        if data[typeID] == "U": continue # Skip inserted gaps

        # Chromosone name
        mychrom = data[chrID]
        if chromname is not None: mychrom = chromname

        # Parse genes on scaffolds
        myscaffold = data[nameID]
        if myscaffold not in genes: continue    # Only match scaffolds that have genes can match against
        genedata = get_gene_positions(genes[myscaffold], int(data[startID]), int(data[endID]), data[oriID], mychrom)
        OUT.writelines(genedata)
    OUT.close()


def get_gene_positions(genes, start, stop, orientation, chrom):
    result=list()
    for gene in genes:
        scaffold_pos = genes[gene]
        assembly_pos=-1
        if orientation=="+" or orientation=="?":    # Assume forward orientation if ambiguous
            assembly_pos=start + scaffold_pos
        elif orientation=="-":
            assembly_pos=stop - scaffold_pos
        result.append("\t".join([gene, chrom, str(assembly_pos)]) + "\n")
    return result


if __name__ == '__main__': main()