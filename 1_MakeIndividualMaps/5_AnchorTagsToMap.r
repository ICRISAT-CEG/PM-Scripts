#! /usr/bin/Rscript

#Map sequencing tags to a genetic map using binomial test
options(stringsAsFactors=F)
args=commandArgs(T)
library("parallel", lib.loc="/usr/lib/R/library")
#args=c("4a_uneak_core_sites_ordered.hmp.txt", "6_merged_tbt.txt", "6a_anchored_tags", "0.001")
mapfile=args[1]
tbtfile=args[2]
outfile=args[3]
pval=as.numeric(args[4])
ncores=as.numeric(args[5])
min.taxa=as.numeric(args[6])
nlines=1e5	#Number of lines of the TBT file to read in at once
metadata=1:11

cat("Anchoring tags to genetic map (only tags in at least",min.taxa,"taxa will be mapped\n")

#####
#Load up and process genetic map
#####
cat("Loading genetic map from",mapfile,"\n")
header=scan(mapfile, nlines=1, what=character(0))
map=read.delim(mapfile, header=T)
names(map)=header
meta=map[,metadata]
map=as.matrix(map[,-metadata])
rownames(map)=meta[,1]

##Process to numeric for binomial testing
alleles=do.call(what=rbind, args=strsplit(meta$alleles, split="/", fixed=T))
calls=matrix(NA, nrow=NROW(map), ncol=NCOL(map), dimnames=dimnames(map))
for(i in 1:NCOL(calls)){
	calls[,i] = ifelse(map[,i]==alleles[,1], yes=0, no=calls[,i])	#Major and minor alleles turned to 0/1 respectively; all others missing
	calls[,i] = ifelse(map[,i]==alleles[,2], yes=1, no=calls[,i])
}

#####
#Anchor tags
#####
cat("Loading TBT file from",tbtfile,"\n")

#Grab data from tbtfile and use to set parameters
ntags=scan(tbtfile, what=numeric(), nlines=1)[1]
if(ntags>nlines){
	cat("  Note: Too many tags to do at once (",ntags,") ; processing",nlines,"lines at a time\n")
}

##Helper function to process probs for boolean map site and TBT files; meant to be parallelized
#Returns a vector of probabilities for all tags against this one marker
get.probs=function(x, tbthit){
	#X should be a list of major and minor calls in boolean form
	min.tests=10
	nMajor = sum(x$major, na.rm=T)	#Precompute probabilities to save a little repetative computation
	nMinor = sum(x$minor, na.rm=T)
	nTotal= nMajor + nMinor	
	
	#Use R's fill-by-column to quickly compute intersections
	transtbt = t(tbthit)
	hitMajors = colSums(transtbt & x$major, na.rm=T)
	hitMinors = colSums(transtbt & x$minor, na.rm=T)
	
	hitTotals=hitMajors + hitMinors
	
	#Filter for minimum hits
	tooFew = which(hitTotals < min.tests)
	hitMajors[tooFew] = NA
	hitMinors[tooFew] = NA

	#Run binomial test	
	pMajors = 0.5 - abs(pbinom(q=hitMajors, size=hitTotals, prob=nMajor/nTotal) - 0.5) 	#Roundabout way to get standardized probability no matter if is at upper or lower end of distribution
	pMinors = 0.5 - abs(pbinom(q=hitMinors, size=hitTotals, prob=nMinor/nTotal) - 0.5)
	
	#Take larger value and turn ones with too few hits to p-values of 1
	probs=ifelse(pMajors > pMinors, yes=pMajors, no=pMinors)
	probs[is.na(probs)]=1
	return(probs)
}


#Previous functions; now trying to rework
###Helper function to process probs for boolean map site and TBT files; meant to be parallelized
##Returns a vector of probabilities for all tags against this one marker
#get.probs=function(x, tbthit){
#	#X should be a list of major and minor calls in boolean form
#	nMajor = sum(x$major, na.rm=T)	#Precompute probabilities to save a little repetative computation
#	nMinor = sum(x$minor, na.rm=T)
#	nTotal= nMajor + nMinor	
#	as.numeric(apply(X=tbthit, FUN=get.probs.one.tag, MARGIN=1, mycalls.list=x, nMajor=nMajor, nMinor=nMinor, nTotal=nTotal))
#}
#
###Another helper function, this to be applied by previous one over matrix of boolean tbt hits
##Returns a single probability of one tag against one marker
#get.probs.one.tag=function(x, mycalls.list, min.tests=10, nMajor=0, nMinor=0, nTotal=0){
#	hitMajor = sum(mycalls.list$major & x, na.rm=T)
#	hitMinor = sum(mycalls.list$minor & x, na.rm=T)
#	if(hitMajor + hitMinor < min.tests){return(1)}
#	pMajor = 0.5 - abs(pbinom(q=hitMajor, size=hitMajor + hitMinor, prob=nMajor/nTotal) - 0.5) 	#Roundabout way to get standardized probability no matter if is at upper or lower end of distribution
#	pMinor = 0.5 - abs(pbinom(q=hitMinor, size=hitMajor + hitMinor, prob=nMinor/nTotal) - 0.5)
#	return(max(pMajor, pMinor, na.rm=T))	#TODO: Depending on how this works, I may want to switch to the mean instead of the max (which is more conservative)
#}



##Helper function to process batches of tags
map.tags=function(tbtfile, meta, calls, nlines=1e5, pval=1e-3, colClasses=NA, outfile, append=F, ncores=4, header=NULL, min.taxa=10){
	#cat("Processing lines",skip,"to",skip+nlines,"\n")
	#Get header and taxa information
	#header=scan(tbtfile, what=character(), nlines=1, skip=1)
	#header=c("tag","length", header)
	colClasses=rep("integer", length(header))
	colClasses[which(! header %in% colnames(calls))]="NULL"	#Skip taxa that aren't in the call file
	colClasses[1] = "character"
	
	#Read in TBT file and convert to matrix
	cat("Reading in next batch of",nlines,"tags\n")
	pretime=Sys.time()
	#tbt=read.delim(tbtfile, skip=skip, nrows=nlines, colClasses=colClasses, header=F, sep="")
	tbt=read.delim(tbtfile, nrows=nlines, colClasses=colClasses, header=F, sep="")
	posttime=Sys.time()
	cat("\tTime to read in tbtfile:", (posttime-pretime),"\n")
	gc()
	names(tbt) = header[colClasses != "NULL"]
	tags=tbt[,1]
	tbt=as.matrix(tbt[,-1])
	rownames(tbt)=tags
	
	#Reorder columns in calls so match
	match=match(x=colnames(tbt), table=colnames(calls))
	calls=calls[,match]
	
	#Change to boolean and remove those with too few markers
	major= calls==0
	minor= calls==1
	tbthit = tbt >=1
	tbthit = tbthit[rowSums(tbthit)>=min.taxa,]	#Keep only those with 
	
	#Make into list for processing
	calls.list=list()
	for(i in 1:nrow(calls)){
		calls.list[[i]]=list()
		calls.list[[i]]$major=major[i,]
		calls.list[[i]]$minor=minor[i,]
	}
	
	#Actually process
	cat("\tCalculating probabilities\n")
	probs=mclapply(X=calls.list, FUN=get.probs, tbthit=tbthit, mc.cores=ncores)
	probs=do.call(what=cbind, args=probs)	#turn into matrix; tags are in rows and markers are in columns
	rownames(probs)=rownames(tbthit)
	colnames(probs)=rownames(calls)
	
	#find lowest probability
	best.probs = apply(X=probs, FUN=min, MARGIN=1, na.rm=T)
	best.probs.index = apply(X=probs, FUN=which.min, MARGIN=1)
	bests=data.frame(tag=rownames(tbthit), best_site=best.probs.index, pval=best.probs, meta[best.probs.index,c(1,3,4)])

	##Debugging code to test peaks, etc.
	#size=8
	#top=probs[best.probs<1e-3,]
	#top=top[sample(1:nrow(top), size=size),]
	#par(mfrow=c(size/2, 2))
	#for(i in 1:nrow(top)){
	#	plot(-log10(top[i,]), type="l")
	#}
	
	#Filter for cutoff 
	cat("\tFiltering and outputting\n")
	probs=probs[best.probs <=pval,]
	bests=bests[best.probs <=pval,]
	

	write.table(bests, file=paste(outfile,"_best.txt", sep=""), row.names=F, col.names=!append, quote=F, sep="\t", append=append)	#Use "append" to determine if adding column names, etc.
	write.table(probs, file=paste(outfile,"_all.txt", sep=""), row.names=T, col.names=!append, quote=F, sep="\t", append=append)	#Use "append" to determine if adding column names, etc.

	return(nlines)
}

#Actually process tags
done=0;
#Pre-open file to ease reading
tbtfilehandle=file(tbtfile, open="r")
header=scan(tbtfilehandle, what=character(), nlines=1, skip=1)
header=c("tag","length", header)
done=map.tags(tbtfilehandle, calls=calls, meta=meta, nlines=nlines, pval=pval, colClasses=colClasses, outfile=outfile, header=header, append=F, ncores=ncores, min.taxa=min.taxa)
while(done < ntags){
	cat("Processed",done,"lines\n")
	done=done + map.tags(tbtfilehandle, calls=calls, meta=meta, nlines=nlines, pval=pval, colClasses=colClasses, outfile=outfile, header=header, append=T, ncores=ncores, min.taxa=min.taxa)
	gc()
	#cat("Processed",done,"tags\n")
}
close(tbtfilehandle)