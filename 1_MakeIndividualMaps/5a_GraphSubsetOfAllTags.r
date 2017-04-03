#! /usr/bin/rscript

#Graph a subset of mapped tags to get an idea of how good the mapping is
options(stringsAsFactors=F)
args=commandArgs(T)
infile=args[1]
outfile=args[2]

#Read in data
cat("Graphing pvalue distribution for tags\n")
data=read.delim(infile, skip=1, header=F)
tags=data[,1]
cat("\t",length(tags),"tags read. Processing...")

#Do -log10 conversion
pvals=-log10(as.matrix(data[,-1]))
for(i in 1:nrow(pvals)){
	infinite=pvals[i,]==Inf
	pvals[i,infinite] = max(pvals[i,!infinite]) + 1
}

#Helper function for rescaling data to show twice
rescale=function(x, in.range=c(0,1), out.range=c(0,1)){
	y = (x - min(in.range) ) / (max(in.range) - min(in.range))
	return (y * (max(out.range) - min(out.range)) + min(out.range))
}

#Format things for graphic
dim=ceiling(sqrt(length(tags)))
max.y=max(pvals, na.rm=T)
png(outfile, width=400*dim, height=300*dim)
	par(mfrow=c(dim,dim),cex=1.2)
	for(i in 1:length(tags)){
		plot(x=1:ncol(pvals), y=pvals[i,], main=paste("Pval dist. for tag",i), ylim=c(0, max.y), xlab=NA, ylab="-log10pval", type="l", lwd=2, lty="solid")
		y.rescale=rescale(pvals[i,], in.range=range(pvals[i,]), out.range=c(0,max.y))
		lines(x=1:ncol(pvals), y=y.rescale, lwd=1, lty="dotted")
		legend(x="topright", lty=c("solid","dotted"), legend=c("Raw","Rescaled"), bty="n", lwd=c(2,1))
	}
dev.off()
