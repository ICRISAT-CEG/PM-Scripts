__author__ = 'jason'

import argparse

debug = False

sam_name, sam_chr, sam_pos = 0, 2, 3    # Columns in sam file

def main():
    args = parse_args()
    passed = filter_contigs(args.infile, args.maxdist, args.multialign)
    output_amplicons(passed, args.outfile)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-m", "--multialign", default=False, action="store_true")
    parser.add_argument("-d", "--maxdist", default=500, type=int, help="Maximum distance between primer pairs")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def filter_contigs(infile, maxdist, multialign=False,):
    print("Filtering",infile,"for good amplicon contigs")
    # First load all the data in the SAM file
    SAM = open(infile, "r")
    all = dict()
    for line in SAM:
        if line.startswith("@"): continue   # Skip header rows
        data = line.strip().split("\t")
        name, chr, pos = data[sam_name], data[sam_chr], int(data[sam_pos])
        amplicon, direction = name.split("_")
        if amplicon not in all:
            all[amplicon]=dict()
        if (not multialign) and (direction in all[amplicon]):
            print("\tWarning! Primer",amplicon,direction,"already loaded!")
        if direction not in all[amplicon]:
            all[amplicon][direction]=list()
        all[amplicon][direction].append(dict(chr=chr, pos=pos, raw=line))
    print("Loaded data on",len(all),"amplicons")

    good, bad = dict(), list()
    for amplicon in all:
        # If only one primer of the pair
        if len(all[amplicon])!=2:
            print("Amplicon",amplicon,"does not have 2 primers listed!")
            bad.append(amplicon)
            continue
        # Go through each pair and determine if mapping is good
        for fwd, rev in zip(all[amplicon]["fwd"], all[amplicon]["rev"]):
             # if scaffolds mapped to don't match
            if fwd["chr"] != rev["chr"]:
                bad.append(amplicon)
                continue
            # if either primer is unmapped
            if fwd["chr"] == "*" or rev["chr"] == "*":
                bad.append(amplicon)
                continue
            # If too much distance between them
            if abs(fwd["pos"] - rev["pos"] > maxdist):
                bad.append(amplicon)
                continue

            # If both primers map to multiple valid locations ("XS:i:" tag from bowtie2); only do for single alignments
            if (not multialign) and ("XS:i:" in fwd["raw"]) and ("XS:i:" in rev["raw"]):
                bad.append(amplicon)
                continue

            # If made it this far, it's a good amplicon
            if amplicon not in good:
                good[amplicon] = list()
            good[amplicon].append(fwd["chr"])
    SAM.close()

    print("Found",len(good),"good amplicons with matched scaffolds and",len(bad),"bad ones with unmatched scaffolds")
    print("Good:",good,"\n\n")
    print("Bad:",bad,"\n\n")
    return good


def output_amplicons(passed, outfile):
    OUT = open(outfile, "w")
    OUT.write("amplicon\tscaffold\n")
    for amplicon in sorted(passed.keys()):
        for scaffold in passed[amplicon]:
            OUT.write(amplicon + "\t" + scaffold + "\n")
    OUT.close()


if __name__ == '__main__': main()