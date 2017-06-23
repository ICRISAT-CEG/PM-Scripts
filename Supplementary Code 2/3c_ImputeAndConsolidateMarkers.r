#! /usr/bin/Rscript

#Identify redundant markers in a hapmap, use them to impute each other, and then reduce to just a single informative marker
options(stringsAsFactors=F)
library("parallel", lib.loc="/usr/lib/R/library")
args=commandArgs(T)
infile=args[1]	#Input hapmap file; should be able to use either normal Hapmap or ABH coded
outfile=args[2]
min.cor=as.numeric(args[3])	#Threshhold for how correlated two markers have to be to count as a single unit. Most conservative is 1
missing.val=args[4]	#Missing value
to.trim=as.logical(args[5])

cat("\nReading in",infile,"and consolidating markers.\n\n")

#Read in data
header=scan(infile, nlines=1, what=character())
hapmap=read.delim(infile, header=T)
names(hapmap)=header
datacols=12:ncol(hapmap)

#Subset just the genotype calls and get the alleles for each
##Helper function to trim out hets and missing values
trim.genos=function(x, to.remove=c("R","Y","S","W","K","M","H", "-", "N", NA)){
	x=x[! x %in% to.remove]
	return(sort(x))
}
genos=as.matrix(hapmap[,datacols])
rownames(genos)=hapmap[,1]
alleles=apply(X=genos, FUN=unique, MARGIN=1)
if(class(alleles)%in%c("matrix", "data.frame")){
	temp=list()
	for(i in 1:ncol(alleles)){
		temp[[i]]=alleles[,i]
	}
	alleles=temp
}
alleles=lapply(X=alleles, FUN=trim.genos)
alleles=do.call(what=rbind, args=alleles)

#Convert to a numeric 0/1 matrix
genos.num=matrix(NA, nrow=nrow(genos), ncol=ncol(genos), dimnames=dimnames(genos))
for(i in 1:ncol(genos)){
	genos.num[,i] = ifelse(genos[,i]==alleles[,1], yes=0, no=ifelse(genos[,i]==alleles[,2], yes=1, no=NA))	#By column, fill in numeric genotype based on which allele it matches
}

#Use hierarchical clustering to get the groups of perfect redundancy
##NOTE: This is a greedy algorithm, so different marker orders can yield different groups
d=dist(genos.num)
clust=hclust(d, method="ward")
assignments=cutree(clust, h=0)	#Cut at height 0 (=perfect association)
groups=list()
for(i in unique(assignments)){
	groups[[i]] = names(assignments)[assignments == i]
}
groupsize = sapply(X=groups, FUN=length)
groups[groupsize<2]=NULL	#Remove 1-snp groups

#Check that everything matches perfectly within each group
check.cor=function(marks, genos.num, min.cor=1){
	subgenos = genos.num[rownames(genos.num) %in% marks,]
	subcors=cor(t(subgenos), method="pearson", use="pairwise.complete.obs")
	#subcors = cors[rownames(cors) %in% marks, colnames(cors) %in% marks]
	if(any(abs(abs(subcors) - min.cor) > 1e-4)){
		return(FALSE)
	}else{
		return(T)
	}
}
good.groups=sapply(X=groups, FUN=check.cor, genos.num=genos.num, min.cor=min.cor)
if(any(!good.groups)){
	cat("Warning! Found some groups by clustering that don't pass quality criteria. Nullifying them.\n")
	groups[!good.groups]=NULL	#Just undo groups that don't fit
}else{
	cat("All",length(groups),"clustered groups pass quality-checking.\n")
}

#Impute within each group, using numeric genotypes
for(g in 1:length(groups)){
	mygroup=groups[[g]]
	mygenos = t(genos.num[rownames(genos.num) %in% mygroup,])
	d = as.matrix(dist(mygenos))
	haveMissing = rowSums(is.na(mygenos))
	toImpute = which(haveMissing>=1 & haveMissing < ncol(mygenos))
	for(i in toImpute){
		neighbors = which(d[i,]==0)
		for(col in 1:ncol(mygenos)){
			vals=unique(mygenos[neighbors, col])
			vals=vals[!is.na(vals)]
			if(length(vals==1)){	#Kludge to reduce processing time when only 1 value to compare
				mygenos[i, col] = vals[1]
			}else{
				values=sort(table(mygenos[neighbors, col]), decreasing=T)
				mygenos[i, col] = as.numeric(names(values)[1])
			}
		}
	}
	mygenos=t(mygenos)
	for(row in 1:nrow(mygenos)){
		genos.num[rownames(genos.num)==rownames(mygenos)[row],] = mygenos[row,]
	}
}

#Turn back into allele calls
genos.imputed=matrix(NA, nrow=nrow(genos), ncol=ncol(genos), dimnames=dimnames(genos))
for(i in 1:ncol(genos)){
	genos.imputed[,i] = ifelse(genos.num[,i]==0, yes=alleles[,1], no=ifelse(genos.num[,i]==1, yes=alleles[,2], no=missing.val))	#By column, fill in imputed genotype based on which allele it matches
}
genos.imputed[is.na(genos.imputed)]=missing.val

#Turn back into a hapmap
hmp.new = hapmap
hmp.new[,datacols] = genos.imputed
rownames(hmp.new)=hmp.new[,1]

#Remove redundant markers, keeping the one with the least amount of missing data
if(to.trim){
	for(g in 1:length(groups)){
		markers=hmp.new[,1]
		submap=hmp.new[hmp.new[,1]%in% groups[[g]], ]
		isMissing=rowSums(submap[,datacols]==missing.val)
		toKeep=which.min(isMissing)
		toExclude = submap[-toKeep,1]
		hmp.new=hmp.new[!rownames(hmp.new) %in% toExclude,]
	}
}

#Output
write.table(hmp.new, file=outfile, quote=F, sep="\t", col.names=T, row.names=F)