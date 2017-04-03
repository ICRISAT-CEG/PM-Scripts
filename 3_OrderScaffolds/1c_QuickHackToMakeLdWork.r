#! /usr/bin/Rscript

#Hack a hapmap file to make it so TASSEL will do LD (by fixing chromosome numbering, etc)
 
args=commandArgs(T)
options(stringsAsFactors=F)
x=read.delim(args[1])
x$chrom=1
x$pos=1:nrow(x)
write.table(x, file=args[2], quote=F, sep="\t", row.names=F, col.names=T)