#Separate a hapmap into specified linkage groupings

cat("Separating hapmaps\n")
options(stringsAsFactors=F)
args=commandArgs(T)
infile=args[1]
outprefix=args[2]
groups=args[3]
chromID=3	#Column with chromosome ID

#Read in data
header=scan(infile, what=character(), nlines=1)
hmp=as.matrix(read.delim(infile, header=T, colClasses="character"))
colnames(hmp)=header
calls=12:NCOL(hmp)

##Define linkage groups
lg=strsplit(groups, split=",")[[1]]
lg=strsplit(lg, split="_")
	
#Go through and subset
for(i in 1:length(lg)){
	cat("\tParsing linkage group",i,"containing clusters",lg[[i]],"\n")
	sub=subset(hmp, hmp[,chromID] %in% lg[[i]])
	
	#filter out monomorphic sites
	include=rep(T, NROW(sub))
	for(j in 1:NROW(sub)){
		alleles = unique(sub[j,calls])
		if(sum(alleles %in% c("A","C","G","T")) < 2){	#If less than 2 non-missing, non-het alleles (note: not genotype calls, just alleles in general), count as false and filter out
			include[j]=F
		}
	}
	sub=sub[include,]
	
	outfile=paste(outprefix,"_LG",i,".hmp.txt",sep="")
	write.table(sub, file=outfile, append=F, row.names=F, col.names=T, sep="\t", quote=F)
}