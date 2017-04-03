#Find highly heterozygous sites and output them as a file

#setwd("/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet/Old_pipelines/DorcusDiversityPop")
#Read in arguments
args=commandArgs(TRUE)
infile=args[1]
outfile=args[2]
cutoff=as.numeric(args[3]) #Cutoff for inclusion

#Load data
colClasses=rep("NULL", 18)
colClasses[c(2,14,16)] = c("character","numeric","numeric")
summary=read.delim(infile, colClasses=colClasses)
names(summary)=c("site","missing","het")
summary=summary[order(summary$het, decreasing=T),]

#Find cutoff and output
badsites=summary$site[summary$het>cutoff]
write(badsites, file=outfile, ncol=1)

#Output graphic
png.name=sub(outfile, pattern=".txt", replacement=".png", fixed=T)
png(png.name, width=1000, height=1000)
	plot(x=1:NROW(summary), y=summary$het, main="Filtering highly het sites", xlab="", ylab="Freq. heterozygous", col=ifelse(summary$het>=cutoff, yes="red", no="black"), pch=20)
	abline(h=cutoff, col="red", lwd=2, lty=2)
	text(x=NROW(summary)/2, y=cutoff+.1, labels="Excluded")
	text(x=NROW(summary)/2, y=cutoff-.1, labels="Included")
dev.off()
