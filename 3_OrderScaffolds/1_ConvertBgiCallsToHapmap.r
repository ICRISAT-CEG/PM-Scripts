#! /usr/bin/Rscript

#Convert BGI's SNP calls to hapmap format for easier slotting in with my tools
#Args: (1) BGI SNP matrix, (2) BGI sample name list, (3) Output file name (gzipped)

#setwd("/home/jgw87/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/BGI_PMiGAP_calls/")
options(stringsAsFactors=F)
args=commandArgs(T)

cat("Converting BGI calls\n")
infile=args[1]
samplenames=args[2]
outfilename=args[3]

#Read in data
bgi=read.delim(infile, sep="", header=F)
samples=scan(samplenames, what=character())

#Add names and swap missing data
names(bgi)[1:3] = c("scaffold", "pos", "ref")
names(bgi)[-(1:3)] = samples
bgi[bgi=="-"] = "N"

#Set up hapmap data frame
hmp=data.frame(rs=paste(bgi$scaffold,bgi$pos, sep="_"), alleles=bgi$ref, chrom=sub(bgi$scaffold, pattern="scaffold", repl=""), pos=bgi$pos)
hmp$strand=hmp$assembly=hmp$center=hmp$protLSID=hmp$assayLSID=hmp$panelLSID=hmp$QCcode=NA
hmp = cbind(hmp, bgi[,-(1:3)])

cat("Writing",nrow(hmp),"SNPs across",ncol(bgi)-3,"taxa to",outfilename,"in gzip format\n")
outfile=gzfile(outfilename, "w")
write.table(hmp, file=outfile, sep="\t", quote=F, row.names=F, col.names=T)
close(outfile)
