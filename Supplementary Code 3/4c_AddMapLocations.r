#! /usr/bin/env Rscript

# Match the linkage group in the final map to the scaffolds 

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile")
parser$add_argument("-a", "--agpfile")
parser$add_argument("-o", "--outfile")
args=parser$parse_args()
# setwd("/home/jgwall/Projects/PM_maps/1_ReorderReferenceScaffolds/Allmaps/")
# args=parser$parse_args(c("-i","4b_good_amplicon_matches.txt","-a","3e_combined_assembly_key.agp","-o","99_test.txt"))

# Load up data
cat("Matching amplicons to linkage groups and positions\n")
amplicons=read.delim(args$infile)

# Match to the original LG and bin
snp_map=
contig_map=read.delim(args$agpfile, skip=1, header=F)
  names(contig_map) = c('chr', 'start','stop','n','code','scaffold','local_start','length','orientation')
contig_map$location=paste(contig_map$chr, contig_map$start, sep="|")
contigmatch = match(amplicons$scaffold, contig_map$scaffold)
amplicons$chr = contig_map$chr[contigmatch]
amplicons$pos = contig_map$start[contigmatch]

# Output
amplicons = amplicons[order(amplicons$chr, amplicons$pos),]
write.table(amplicons, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)

