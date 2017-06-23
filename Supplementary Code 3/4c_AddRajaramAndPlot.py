__author__ = 'jason'

import argparse
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

debug = False

# Z orders for the various plot components
goodline_z=11
badline_z=4
chrom_z=10
tick_z=20

# Color key for chromosome linkage
goodcolors= {1:"red", 2:"orange", 3:"yellow", 4:"green", 5:"darkturquoise", 6:"blue",7:"indigo"}
badcolors= {1:"lightcoral", 2:"moccasin", 3:"khaki", 4:"palegreen", 5:"aquamarine", 6:"lightblue",7:"mediumpurple"}

def main():
    args = parse_args()
    amplicons = pd.read_csv(args.infile, sep='\t')

    # Add Rajaram map positions
    amplicons = add_rajaram(amplicons, args.rajaramfile)

    if args.rescale_chroms:
        amplicons=rescale_chroms(amplicons, "chr", "pos")
        amplicons=rescale_chroms(amplicons, "raj_chr", "raj_pos")

    # Strip "chr" part of linkage group name
    amplicons["chr"] = [float(str(c).replace("chr", "")) for c in amplicons["chr"]]

    # Output text file
    if args.outfile:
        print("Outputting text results to",args.outfile)
        amplicons.to_csv(args.outfile, sep='\t', index=False)

    # Make graphic
    if args.graphfile:
        print("Outputting graphical results to",args.graphfile)
        output_graphic(amplicons, args.graphfile, args.flip_chr)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-r", "--rajaramfile")
    # parser.add_argument("--maplengths")
    # parser.add_argument("--rajlengths", help="Lengths of linakge groups for Rajaram")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-g", "--graphfile")
    parser.add_argument("--rescale-chroms", default=False, action="store_true",help="Whether to normalize chromosome lengths")
    parser.add_argument("-f", "--flip-chr", type=int, nargs="*", help="Assembly map chromosomes to flip")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def add_rajaram(amplicons, rajaramfile):
    print("Adding Rajaram et al's map locations")

    # Parse header
    IN = open(rajaramfile, "r")
    header = IN.readline().strip().split('\t')
    lgID, binID, markerID = header.index("lg"), header.index("bin"), header.index("marker")

    # Set up data structures for finding and storing results
    primerkey = [p.lower().lstrip("x") for p in amplicons["amplicon"]]
    primerlookup = set(primerkey)   # To hopefully speed things up
    amplicons["raj_chr"], amplicons["raj_pos"] = "unknown", -1

    # Go through and match up targets
    for line in IN:
        data = line.strip().split('\t')
        marker = data[markerID].lower().lstrip("x")
        if marker in primerlookup:
            index = primerkey.index(marker)
            amplicons["raj_chr"].iloc[index] = data[lgID]
            amplicons["raj_pos"].iloc[index] = data[binID]
    IN.close()
    # print(amplicons)

    # Quality check that all amplicons accounted for
    matched = np.sum(amplicons["raj_pos"] != -1)
    print("Matched",matched,"out of",len(amplicons),"amplicons to the Rajaram map")
    not_matched = sorted(amplicons["amplicon"].loc[amplicons["raj_pos"] == -1])
    print("\tWarning! The following",len(not_matched),"amplicons were not matched to Rajaram linkage map:", not_matched)
    return amplicons

def output_graphic(amplicons, graphfile, flip_chr):
    # Set up plotting parameters
    amplicons = amplicons.loc[amplicons["raj_chr"]!="unknown",:]
    amplicons = amplicons.loc[pd.notnull(amplicons["chr"]),:]
    amplicons["raj_chr"] = [int(r) for r in amplicons["raj_chr"]]
    amplicons["raj_pos"] = [float(r) for r in amplicons["raj_pos"]]
    #print(amplicons)


    nchr = len(np.unique(amplicons["chr"]))
    map_x = np.arange(nchr)  # Positions to plot each linkage group
    raj_x = map_x + 0.33
    map_key = {chrom:pos for chrom, pos in zip(sorted(np.unique(amplicons["chr"])), map_x)} # Lookup for chrom positions
    raj_key = {chrom:pos for chrom, pos in zip(sorted(np.unique(amplicons["raj_chr"])), raj_x)}

    fig = plt.figure(figsize=(15,5))
    ax = fig.add_subplot(111, title="Comparison of Linkage Maps", xlabel="Linkage Groups", ylabel="Position")

    # Add chromosome bars
    ax, map_lengths = add_linkage_bars(ax, amplicons, "", map_x, "dimgray", flip_chr)
    ax, raj_lengths = add_linkage_bars(ax, amplicons, "raj", raj_x, "darkgray")

    # TODO: Rescale maps so can be plotted on the same graph

    # Draw lines between the linkage groups
    som_offset, raj_offset = 0.07, 0.04
    for i in range(len(amplicons)):
        map_lg = amplicons["chr"].iloc[i]
        raj_lg = amplicons["raj_chr"].iloc[i]
        start_x = map_key[map_lg] + som_offset
        start_y = amplicons["pos"].iloc[i]
        # To flip to see if some linkage groups work better in opposite orientation
        if (flip_chr is not None) and (map_lg in flip_chr):
            start_y = map_lengths[map_lg] - start_y
        end_x = raj_key[raj_lg] - raj_offset
        end_y = amplicons["raj_pos"].iloc[i]
        color, linestyle, linewidth, zorder = goodcolors[map_lg], "-", 1.5, goodline_z  # Set values to "good", matching linkage groups by default
        # If linkage groups don't match, set to "bad" values
        if str(int(map_lg)) != str(int(raj_lg)):    # Sort of hackish way to make sure LGs match
            color=badcolors[map_lg]
            linestyle="--"
            zorder=badline_z
            linewidth=0.5
        myline = lines.Line2D([start_x, end_x], [start_y, end_y], linewidth=linewidth, linestyle=linestyle, color=color, zorder=zorder)
        ax.add_line(myline)

    # Set x/y limits
    ymin = np.min([np.min(amplicons["pos"]), np.min(amplicons["raj_pos"])])
    ymax = np.max([np.max(amplicons["pos"]), np.max(amplicons["raj_pos"])])
    yspan = ymax - ymin
    ymin-=yspan/10
    ymax+=yspan/10
    ax.set_xlim(left = np.min(map_x)-0.5, right=np.max(raj_x) + 0.5)
    ax.set_ylim(bottom = ymin, top=ymax)

    fig.savefig(graphfile, dpi=200)

def add_linkage_bars(ax, amplicons, stem, xvals, color, flip_map=None):
    print("Adding bars for",stem)
    chrID = (stem + "_chr").lstrip("_") # left-strip underscore in case the stem is ""
    posID = (stem + "_pos").lstrip("_")
    bar_width=0.05
    line_width = bar_width * 2

    # Get linkage group lengths
    groups, lengths = list(),list()
    for lg in sorted(np.unique(amplicons[chrID])):
        length = np.max(amplicons[posID].loc[amplicons[chrID] == lg])
        groups.append(lg)
        lengths.append(length)

    # Plot bars and ticks for matched markers
    for group, length, x in zip(groups, lengths, xvals):
        mygroup = amplicons.loc[amplicons[chrID]==group,:]
        # Add rectangle showing entire linkage group
        bar = patches.Rectangle(xy = (x, 0), width=bar_width, height=length, color=color, zorder=chrom_z)
        ax.add_patch(bar)
        ax.text(x, -15, stem + str(group), color=color)
        # Add lines for each marker position
        xleft = x - line_width/3
        for cm in mygroup[posID]:
            if (flip_map is not None) and group in (flip_map):
                cm = length - cm
            ax.add_line(lines.Line2D([xleft, xleft + line_width], [cm, cm], linewidth=1, color='black', zorder=tick_z))
    length_key = {group:length for group,length in zip(groups, lengths)}
    return ax,length_key


def rescale_chroms(amplicons, chrID, posID):
    print("Rescaling chromosome positions to 0-1 scale for columns",chrID, "&",posID)
    chroms = np.array(amplicons[chrID])
    poses = np.array([float(p) for p in amplicons[posID]])   # Numpy much faster and easier to slice
    chromset = set(chroms)
    for c in chromset:
        targets = chroms == c
        if np.sum(targets) >=1:
            mymax = max(poses[targets])
            poses[targets] = poses[targets] / mymax
    amplicons[posID] = poses
    return amplicons


if __name__ == '__main__': main()