#Take a numeric output from TASSEL and cluster for LD analysis

#Read in data
args=commandArgs(T)
raw=as.matrix(read.delim(args[1], skip=2, header=T, row.names=1))
hmp.name=args[2]
header=scan(hmp.name, nlines=1, what=character())
hmp=read.delim(file=hmp.name, header=T)
names(hmp)=header
clusters=args[3]	#Comma-separated clusters to cut at
outprefix=args[4]

cutoff=0.8	#missing data cutoff, any markers with more than this fraction of missing or heterozygous data will be excluded
clusters=as.numeric(do.call(what=c, args=strsplit(clusters, split=",")))

#Remove het calls, and transpose so are looking at sites instaed of taxa
nohet=ifelse(raw==0.5, yes=NaN, no=raw)
nohet=t(nohet)

#Filter out sites with too high missing data
fract.missing = rowSums(is.na(nohet)) / NCOL(nohet)
filt=nohet[fract.missing < cutoff,]
goodsites=as.integer(sub(rownames(filt), pattern="^S", replacement="", perl=T)) + 1	#+1 is b/c TASSEL starts at 0 while R starts at 1

#Calc distance and use hierarchical clustering to find 
d=dist(filt)
clust.modes=c("ward")
clusts=list()
groups=matrix(NA, nrow=NROW(filt), ncol=length(clusters)*length(clust.modes))
group.names=rep("", NCOL(groups))
group.i=1
for(i in 1:length(clust.modes)){
	mode=clust.modes[i]
	clusts[[i]] = hclust(d, method=mode)
	for(size in clusters){
		groups[,group.i]=cutree(clusts[[i]], k=size)
		group.names[group.i] = paste(mode,"_k",size,sep="")
		group.i=group.i+1
	}
}
dimnames(groups)=list(rownames(filt), group.names)

#Loop over and export relevant data
hmp.filt=hmp[goodsites,]
for(i in 1:NCOL(groups)){
	hmp.sub = hmp.filt
	hmp.sub$chrom = groups[,i]
	hmp.sub=hmp.sub[order(hmp.sub$chrom, hmp.sub$pos),]
	filename=paste(outprefix, "_",colnames(groups)[i],".hmp.txt", sep="") 
	write.table(hmp.sub, file=filename, row.names=F, col.names=T, quote=F, sep="\t")
}