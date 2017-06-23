__author__ = 'jgwall'

import argparse

debug = False


def main():
    args = parse_args()
    key = make_key(args.agpfiles)
    remake_hapmap(args.infile, args.outfile, key)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Hapmap with genotpyes per scaffolds")
    parser.add_argument("-o", "--outfile", help="Reordered hapmap")
    parser.add_argument("-a", "--agpfiles", nargs="*", help=".agp files with new scaffold positions")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def make_key(infiles):
    print("Loading reorder keys from",len(infiles),"input AGP files")
    key=dict()
    for infile in infiles:
        print("\tLoading",infile)
        IN = open(infile)
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
            myscaffold = data[nameID]
            mychrom = data[chrID]
            mypos = data[startID]
            if myscaffold in key:
                print("\t\tWARNING! Scaffold",myscaffold,"already loaded on",key[myscaffold],"but will be replaced by",[mychrom, mypos])
            key[myscaffold] = [mychrom, mypos]
    print("\tLoaded",len(key),"total scaffolds")
    return key

def remake_hapmap(infile, outfile, key):
    print("Remaking hapmap from",infile)
    IN=open(infile, 'r')
    OUT = open(outfile, 'w')

    # Parse header
    header=IN.readline().strip().split('\t')
    scaffoldID, chrID, posID = [header.index(c) for c in ['snp', 'chrom', 'pos']]
    OUT.write('\t'.join(header) + "\n")

    # Redo hapmap
    skipped, moved = 0, 0
    for line in IN:
        data=line.strip().split('\t')
        if data[scaffoldID] not in key: # Skip scaffolds with unknown positions
            skipped+=1
            continue
        data[chrID], data[posID] = key[data[scaffoldID]]
        OUT.write("\t".join(data) + "\n")
        moved +=1
    print("\tMoved",moved,"recognized scaffolds and skipped",skipped,"unrecognized ones")

    IN.close()
    OUT.close()


if __name__ == '__main__': main()