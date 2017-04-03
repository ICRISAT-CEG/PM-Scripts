#! /usr/bin/Rscript
#Format my scaffold LD measures for use in the traveling salesman problem
#Args: (1) Input "distance" file, (2) Output file prefix

options(stringsAsFactors=F)
args=commandArgs(T)
#setwd("/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/BGI_PMiGAP_calls/")
#args=c("1f_bgi_pmigap_outer100.mean_ld.txt.gz", "2_outer100_formatted")

infile=args[1]
outprefix=args[2]

#Load data
firstrow=scan(infile, what=character(), nlines=1)
colClasses=rep("numeric", length(firstrow))
colClasses[1]="character"
distance = read.delim(infile, header=F, colClasses=colClasses, na.strings=c("NA","NaN",'�'), row.names=1)
distance = as.matrix(distance)
colnames(distance) = rownames(distance)

#Trim down to only those that have >90% relations present
nasums = rowSums(is.na(distance))
goodgroups = (nasums < ncol(distance) * 0.9)
trimmed = distance[goodgroups, goodgroups]

#Set NAs to -1 (= impossibly low LD) and then invert distances (so high LD = low distance)
trimmed[is.na(trimmed)] = -1
trimmed = 1 - trimmed 

#Modify so that left-right pairs are distance 0 from each other, and all sets are distance 0 from themselves
fixed=trimmed
groups=unique(sub(rownames(fixed), pattern="\\|.+", repl=""))
t1=proc.time()[3]
coords=list()
for(g in groups){
	mygroup = grep(rownames(fixed), pattern=paste("^",g,"\\|", sep=""))
	if(length(mygroup)==1){
		coords[[g]] = c(mygroup, mygroup)
	}else if(length(mygroup)==2){
		coords[[g]] = rbind(mygroup, rev(mygroup), c(min(mygroup), min(mygroup)), c(max(mygroup), max(mygroup)))
	}else{
		cat("Mygroup error for",g,":",rownames(fixed)[mygroup],"\n")
	}
}
cat("Required",proc.time()[3]-t1,"seconds to fix everything\t")
coords=do.call(what=rbind, args=coords)
fixed[coords]=0

write.table(fixed, file=paste(outprefix, ".txt",sep=""), sep="\t", quote=F, row.names=T, col.names=F)
save(fixed, file=paste(outprefix, ".rdata",sep=""))
