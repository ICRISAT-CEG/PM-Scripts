#! /usr/bin/env Rscript

# Give an AGP file and tell which chromosomes to be flipped

library(argparse)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i","--infile")
parser$add_argument("-o","--outfile")
parser$add_argument("-f","--flip", nargs="*", help="Chromsomes to flip")
args=parser$parse_args()
#setwd('/home/jgwall/Projects/PM_maps/1_ReorderReferenceScaffolds/Allmaps/')
#args=parser$parse_args(c("-i","3e_combined_assembly_key.agp","-o","99_flipped.agp","-f","chr1", "chr2"))


# Load data
cat("Flipping chromosomes",args$flip,"in agp file",args$infile,"\n")
header=scan(args$infile, nlines=1, sep="\n", what=character())
agp = read.delim(args$infile, header=F, comment.char="#")
names(agp) = c('chrom','chr_start','chr_stop','num','code','name','x','length','orientation')

# Split by chromosome and flip the indicated ones
chroms = split(agp, agp$chrom)
chroms[args$flip] = lapply(chroms[args$flip], function(x){
  # Flip order
  x = x[nrow(x):1,]	
  
  # Change orientation
  new_ori = x$orientation
  new_ori[x$ori=="+"] = "-"	
  new_ori[x$ori=="-"] = "+"
  x$orientation = new_ori
  
  # Change start/end coords
  length = max(x$chr_stop)
  newstart = length - x$chr_stop + 1
  newstop = length - x$chr_start + 1
  x$chr_stop = newstop
  x$chr_start = newstart
  
  # Change number
  x$num = 1:nrow(x)
  return(x)
})

# Recombine
output = do.call(what=rbind, args=chroms)
write(header, file=args$outfile)
write.table(output, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=F, append=T)
