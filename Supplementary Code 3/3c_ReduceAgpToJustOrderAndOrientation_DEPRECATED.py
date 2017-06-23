__author__ = 'jgwall'

import argparse

debug = False


def main():
    args = parse_args()
    manuals = get_manual_set(args.manual_locations)
    reduce_agp(args.infile, args.outfile, manuals)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-m", "--manual-locations", nargs="*", help="List of manual linkage groups to put scaffolds on "
                                                                    "in form of 'scaffold:chr' (for removing them from wrong chromosomes)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def get_manual_set(locations):
    locs = dict()
    for l in locations:
        scaffold, chrom = l.split(':')
        if scaffold in locs:
            print("WARNING! Scaffold".scaffold,"told to be put in multiple locations! Only latest (",chrom,") will be used")
        locs[scaffold] = chrom
    print("Loaded manual locations for",len(locs),"scaffolds:",locs)
    return locs


def reduce_agp(infile, outfile, locations):
    print("Reducing AGP file",infile)
    IN = open(infile)
    OUT = open(outfile, "w")
    OUT.write("\t".join(["chr", "scaffold","orientation"]) + "\n")

    chrID, startID, endID, typeID, nameID, oriID = [None] * 6
    n, written =0 ,0
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
        n+=1
        mychrom = data[chrID]
        myscaffold = data[nameID]
        myori = data[oriID]

        if (myscaffold in locations) and (locations[myscaffold]) != mychrom:
            print("\tSkipping",myscaffold,"because is on chrom",mychrom,"but manually specified to be on chrom",locations[myscaffold])
            continue
        OUT.write("\t".join([mychrom, myscaffold, myori]) + '\n')
        written+=1
    OUT.close()
    print("\tLoaded",n,"scaffolds and wrote",written,"to",outfile)


if __name__ == '__main__': main()
