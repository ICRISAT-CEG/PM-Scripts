#! /usr/bin/Rscript

#Sort a hapmap based on chromosome, position, and then name

args=commandArgs(T)
infile=args[1]
outfile=args[2]

hmp = read.delim(infile, header=T)
header = scan(infile, nlines=1, what=character())
new_order = order(hmp$chrom, hmp$pos, hmp$rs.)
hmp = hmp[new_order,]
names(hmp) = header

write.table(hmp, file=outfile, sep="\t", quote=F, row.names=F, col.names=T)
