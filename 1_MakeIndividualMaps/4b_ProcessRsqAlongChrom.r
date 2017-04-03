#! /usr/bin/Rscript

#Manuall check the assignment along chromosomes for the monsanto SNPs to see how my Java script works

options(stringsAsFactors=F)
args=commandArgs(T)
rsqfile=args[1]
n.subset=as.numeric(args[2])
min.rsq=as.numeric(args[3])
outprefix=args[4]

cat("Converting", rsqfile,"to visual output for",n.subset,"random snps\n")

#Read in data
header=scan(rsqfile, nlines=1, what=character())
colClasses=rep("numeric", length(header))
colClasses[1]="character"
rsq=as.matrix(read.delim(rsqfile, header=T, colClasses=colClasses, row.names=1))

#Find best position
rsq[is.na(rsq)]=0
best.rsq=apply(X=rsq, FUN=max, MARGIN=1, na.rm=T)
best.pos=apply(X=rsq, FUN=which.max, MARGIN=1)

#Pick random subset
mysub=sort(sample(1:nrow(rsq), size=n.subset))

#plot
dim=ceiling(sqrt(n.subset))
png(paste(outprefix,".png",sep=""), width=500*dim, height=200*dim)
	par(mfrow=c(dim, dim))
	for(i in 1:n.subset){
		toplot=mysub[i]
		plot(rsq[toplot,], type="l", main=paste("LD for",rownames(rsq)[toplot]), ylim=c(0,1))
		my.x=best.pos[toplot]
		my.y=best.rsq[toplot]
		points(x=my.x, y=my.y, pch=20, col="red", cex=2)
		abline(h=min.rsq, col="blue", lty="dashed")
	}
dev.off()

#Output filtered results
cat("Outputting parsed results for well-anchored SNP (Rsq >=",min.rsq,")\n")
tokeep = which(best.rsq >= min.rsq)
cat("   Total", length(tokeep), "sites pass filtering\n")
good.sites=rownames(rsq)[tokeep]
good.rsq=best.rsq[tokeep]
neighbors=colnames(rsq)[best.pos]
good.neighbors=neighbors[tokeep]
output=data.frame(target=good.sites, neighbor=good.neighbors, rsq=good.rsq)
write.table(output, file=paste(outprefix, "_goodsites.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)

#Write out bad sites
write(rownames(rsq)[-tokeep], file=paste(outprefix, "_badsites.txt", sep=""), ncolumns=1)