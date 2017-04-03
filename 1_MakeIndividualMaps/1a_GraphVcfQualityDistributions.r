#Graph distribution of SNP calls by quality and % heterozygosity

#setwd("/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet_SomSekhar/Jason_Analysis/VCF_filter")
#infile="3_ref_snps_vcf_qualities.txt"
#outfile="3_ref_snps_vcf_dist"
args=commandArgs(T)
infile=args[1]
outfile=args[2]
binsize=1
#binsize=as.numeric(args[3])

cat("Graphing VCF quality scores and distributions from",infile,"\n")

data=read.delim(infile, header=T, colClasses=c("integer","character","NULL", "NULL"))
bins=seq(from=0, to=max(data$GC)+binsize, by=binsize)
data$bin = cut(x=data$GC, breaks=bins, labels=FALSE, right=FALSE)
data$n=1

#Calculate summary stats
counts=rep(0, length(bins))
hets=rep(0, length(bins))
for(mybin in unique(data$bin)){
	this.bin = data$bin==mybin
	counts[mybin] = sum(this.bin)	#Total calls
	hets[mybin] = sum(data$allele == "Het" & this.bin) / counts[mybin]	#Fraction heterozygous
}

#Plot data
png(paste(outfile,".png",sep=""), width=700, height=400)	
	par(mar=c(5,4,4,4), cex=1.3)
	plot(x=bins, y=counts, main="Distribution of VCF calls", xlab="GC quality score", ylab="Count", type="l", col="black")
	lines(x=bins, y=hets * max(counts), col="blue")
	axis(side=4, at=seq(from=0, to=1, by=0.25) * max(counts), labels=seq(from=0, to=1, by=0.25), col="blue")
	mtext(line=3, side=4, text="Fraction heterozygous", col="blue")
dev.off()

final=data.frame(GC=bins, nsites=counts, fract_het=hets)
write.table(final, file=paste(outfile,".txt",sep=""), quote=F, sep="\t", row.names=F, col.names=T)