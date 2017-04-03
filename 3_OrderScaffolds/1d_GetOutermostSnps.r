#! /usr/bin/Rscript

#Take a hapmap and extract the list of SNPs in the outermost X of each chromosome/scaffolds
#Output file is 3 columns, with SNP name, scaffold, and side (left vs right)
#Args: (1) INput file, (2) Distance tolerance, (3) Output file prefix

args=commandArgs(T)
options(stringsAsFactors=F)
#setwd("/home/jgw87/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/BGI_PMiGAP_calls/")
#args=c("1b_bgi_pmigap_filt.hmp.txt.gz","1000","1d_outer1000")

infile=args[1]
dist=as.numeric(args[2])
outprefix=args[3]

#Read in data
cat("Finding tips within",dist,"bp of scaffold ends in", infile,"\n")
header=scan(infile, nlines=1, what=character())
colClasses=rep("NULL", length(header))
colClasses[c(1,3)] = "character"
colClasses[4] = "numeric"
hmp=read.delim(infile, colClasses=colClasses)

#Helper function to identify left-right SNPs on each chrom
find.tips=function(x, dist=1000){
	x$left = x$pos <= (min(x$pos) + dist)
	x$right = x$pos >= (max(x$pos) - dist)
	x=subset(x, x$left | x$right)
	return(x)
}

#Split by chromosome
chroms = split(hmp, f=hmp$chrom)
cat("\tSplit into",length(chroms),"unique scaffolds\n")
tips = lapply(chroms, FUN=find.tips, dist=dist)

#Get distributions of distances
before.dist = unlist(lapply(chroms, FUN=nrow))
after.dist=unlist(lapply(tips, FUN=nrow))

#Graph distributions
png(paste(outprefix,"_distributions.png", sep=""), width=600, height=800)
	par(mfrow=c(2,1), cex=1.1)
	nbreaks=40
	hist(before.dist, breaks=nbreaks, main="Raw scaffold size distribution")
	hist(after.dist, breaks=nbreaks, main="Scaffold (tips only) size distribution")
dev.off()

#Convert to single data frame and output
output=do.call(what=rbind, args=tips)
cat("Found",sum(output$left),"left tips and",sum(output$right),"right tips (",sum(output$left & output$right),"are both)\n")


lefts = subset(output, output$left)
rights = subset(output, output$right)
lefts=data.frame(snp=lefts$rs., scaffold=lefts$chrom, side="left")
rights = data.frame(snp=rights$rs., scaffold=rights$chrom, side="right")
output = rbind(lefts, rights)
output=output[order(output$scaffold, output$side, output$snp),]
write.table(output, file=paste(outprefix,"_tips.txt", sep=""), quote=F, row.names=F, col.names=T, sep="\t")
