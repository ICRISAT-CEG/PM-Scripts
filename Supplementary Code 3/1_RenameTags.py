__author__ = 'jgwall'

import argparse

debug = False


def main():
    args = parse_args()

    print("Building key to rename tags in",args.infile)
    key = dict()
    for line in open(args.keyfile):
        tag, name = line.strip().split()
        key[standardize_tag(tag)]=name
    print("\tLoaded",len(key),"tags to rename")

    # Process input file
    print("Renaming tags")
    IN = open(args.infile, 'r')
    OUT = open(args.outfile, 'w')

    # Deal with header
    header=IN.readline().strip().split('\t')
    tagID = header.index("tag")
    header[tagID] = "rs#"
    OUT.write("\t".join(header) + "\n")

    # Change each line
    found, missed = 0, 0
    for line in IN:
        data = line.strip().split('\t')
        mytag = standardize_tag(data[tagID])
        if mytag in key:
            data[tagID] = key[mytag]
            OUT.write("\t".join(data)+'\n')
            found+=1
        else:
            #print("Warning! Tag",mytag,"not found in input key")
            missed+=1
    print("\tFound",found,"tags that were renamed and",missed,"tags that could not be found and were not passed along")
    IN.close()
    OUT.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-k", "--keyfile", help="Two-column keyfile of tag (from SAM file) and new name")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


# Just a way to make sure tags are standardized; set as a separte function to make extensible
def standardize_tag(tag):
    return tag.rstrip("A")

if __name__ == '__main__': main()