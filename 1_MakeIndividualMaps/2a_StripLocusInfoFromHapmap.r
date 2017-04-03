#! /usr/bin/Rscript

#Strip the scaffold and position numbers from a hapmap and combine into a single locus
#Arguments: (1) Input hapmap, (2) Output hapmap

args=commandArgs(T)
infile=args[1]
outfile=args[2]

#Read in data
header=scan(infile, nlines=1, what=character())
data=read.delim(infile, header=T)
names(data)=header	#Put right column names back

#Redo locus numbers
data$chrom=1
data$pos = (1:nrow(data) - 1) * 1000 + 1

#write back out
write.table(data, file=outfile, row.names=F, col.names=T, quote=F, sep="\t")