#Calculate and show # of rare alleles by taxon


cat("Testing for bad taxa via too many rare alleles\n")
args=commandArgs(TRUE)	
infile=args[1]
outfile=args[2]

nrows=1000	#Overestimate of how many rows need

#Read in data
header=scan(infile, what=character(), nlines=1)
data=scan(infile, what=character(), skip=3)
data=matrix(data, ncol=length(header), byrow=T)
taxa=data[,1]
data=matrix(as.numeric(data[,-1]), nrow=nrow(data))
rownames(data)=taxa
data=ceiling(data)	#Round hets up to 1 to just count "instances of a rare allele"

#Do significance test
rarecount = rowSums(data, na.rm=T)	#Count of # rare allele calls in each line
totalsnps = rowSums(!is.na(data))	#Count of total called SNPs in each line
rarefreq=rarecount/totalsnps
sig.test = pnorm(q=rarefreq, mean=mean(rarefreq, na.rm=T), sd=sd(rarefreq, na.rm=T), lower.tail=F)

#TODO: Switch this so take anything more than 2 SD away from the mean?


fdr = p.adjust(sig.test, method="fdr")	#False discovery rate 
col = ifelse(fdr < 0.05, yes="blue", no="black")
col = ifelse(fdr < 0.01, yes="red", no=col)

png(file=paste(outfile, ".png", sep=""), width=2000, height=1000)
	par(mfrow=c(1,2), cex=1.2)
	plot(x=1:length(rarefreq), y=rarefreq, pch=20, col=col, main="Rare alleles by line", xlab="line", ylab="Rare allele counts")
	problems = fdr < 0.05
	if(sum(problems, na.rm=T)>0){
		text(x=(1:length(rarefreq))[problems], y=rarefreq[problems], labels=rownames(data)[problems], pos=4, col=col[col!="black"], cex=0.7)	
		output=rownames(data)[problems]
	}else{
		text(x=1, y=min(rarefreq, na.rm=T), labels="No problem taxa identified", pos=4, col="darkgreen", cex=2)	
		output=""
	}
	plot(x=1:length(sort(rarefreq)), y=sort(rarefreq), pch=20, col=col[order(rarefreq)], main="Rare alleles by line (ordered)", xlab="line", ylab="Rare allele freq")
dev.off()

write.table(file=paste(outfile,".txt",sep=""), x=output, row.names=F, col.names=F, sep="\t", quote=F)