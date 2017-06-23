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

    chrID, startID, endID, typeID, nameID, lengthID, oriID = [None] * 7
    n, written =0, 0
    offset=0
    for line in IN:
        # parse header
        if line.startswith("# FIELDS"):
            OUT.write(line)
            line=line.replace("# FIELDS:", "")
            line=line.replace(" ", "")
            header = line.strip().split(',')
            chrID, startID, endID, typeID, nameID, lengthID, oriID = \
                [header.index(c) for c in ["object", "object_beg", "object_end", "component_type", "component_id/gap_length", "component_end/linkage", "orientation/linkage_evidence"]]
            continue
        elif line.startswith("#"): continue   # Skip comments

        # Parse actual input lines
        n+=1
        data=line.strip().split()
        # Test if need to delete scaffold
        mychrom, myscaffold = data[chrID], data[nameID]
        if (myscaffold in locations) and (locations[myscaffold]) != mychrom:
            gapsize = int(IN.readline().strip().split()[nameID])
            offset += int(data[lengthID]) + gapsize     # Add the length of this scaffold to the offset
            print("\tSkipping",myscaffold,"because is on chrom",mychrom,"but manually specified to be on chrom",locations[myscaffold],"; offset is now",offset)
            continue

        data[startID] = int(data[startID]) - offset
        data[endID] = int(data[endID]) - offset
        # data.append(data[startID] + offset);data.append(data[endID] + offset) # For debugging, to check that it's working okay
        data = [str(d) for d in data]
        OUT.write("\t".join(data) + '\n')
        written+=1
    OUT.close()
    print("\tLoaded",n,"scaffolds and wrote",written,"to",outfile)


if __name__ == '__main__': main()
