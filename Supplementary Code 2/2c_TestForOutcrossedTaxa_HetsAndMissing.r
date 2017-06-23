#Find highly heterozygous sites and output them as a file

#setwd("/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet/Old_pipelines/DorcusDiversityPop")
#Read in arguments
cat("Testing for potential outcrosses via excess heterozygosity\n")
args=commandArgs(TRUE)
infile=args[1]
outfile=args[2]
het.cutoff=as.numeric(args[3]) #Cutoff for inclusion
miss.cutoff=as.numeric(args[4])

#Load data
colClasses=rep("NULL", 9)
colClasses[c(2,5,7)] = c("character","numeric","numeric")
summary=read.delim(infile, colClasses=colClasses)
names(summary)=c("taxon","missing","het")

#Test for excess heterozygosity
too.het=summary$het>het.cutoff
too.het[is.na(too.het)]=F
#Test for excess missingnless
too.missing=summary$missing>miss.cutoff
too.missing[is.na(too.missing)]=F

badtaxa=summary$taxon[too.het | too.missing]
if(length(badtaxa)==0){	#Just in case nothing clears the cutoff
	badtaxa=""
}
write(badtaxa, file=outfile, ncol=1)

#Output graphic
png.name=sub(outfile, pattern=".txt", replacement=".png", fixed=T)
png(png.name, width=2000, height=1000)
	par(mfrow=c(1,2))
	#Het plot
	summary=summary[order(summary$het, decreasing=T),]
	plot(x=1:NROW(summary), y=summary$het, main="Filtering highly het taxa", xlab="", ylab="Freq. heterozygous", col=ifelse(summary$het>=het.cutoff, yes="red", no="black"), pch=20)
	if(any(too.het)){
		toremove=summary$het>het.cutoff
		text(x=(1:nrow(summary))[toremove], y=summary$het[toremove], labels=summary$taxon[toremove], pos=4, col="red", cex=1.5)	
	}
	abline(h=het.cutoff, col="red", lwd=2, lty=2)
	text(x=NROW(summary)/2, y=het.cutoff+.1, labels="Excluded")
	text(x=NROW(summary)/2, y=het.cutoff-.1, labels="Included")
	
	#Missing plot
	summary=summary[order(summary$missing),]
	plot(x=1:NROW(summary), y=summary$missing, main="Filtering highly missing taxa", xlab="", ylab="Freq. missing", col=ifelse(summary$missing>=miss.cutoff, yes="red", no="black"), pch=20)
	if(any(too.missing)){
		toremove=summary$missing>miss.cutoff
		text(x=(1:nrow(summary))[toremove], y=summary$missing[toremove], labels=summary$taxon[toremove], pos=2, col="red", cex=1.5)	
	}
	abline(h=miss.cutoff, col="red", lwd=2, lty=2)
	text(x=NROW(summary)/2, y=miss.cutoff+.1, labels="Excluded")
	text(x=NROW(summary)/2, y=miss.cutoff-.1, labels="Included")
	
	
dev.off()

