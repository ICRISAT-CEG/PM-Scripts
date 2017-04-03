#Parse the compiled, bootstrapped map into a (presumably) more stable configuration

options(stringsAsFactors=F)
args=commandArgs(T)
infile=args[1]
hapmap=args[2]
outfile=args[3]
mychr=args[4]

#Read in data
map=read.delim(infile, header=F, row.names=1)
hmp=read.delim(hapmap, header=T)
header=scan(hapmap, nlines=1, what=character())
names(hmp)=header

#Parse maps
##Split off SNP names and turn map into a matrix for easier handling
snps = row.names(map)
map=as.matrix(map)
##Swap values for those sets where the chromosome is arbitrarily ordered backwards relative to the first column
for(i in 2:NCOL(map)){
	if(cor(map[,1], map[,i], use="complete.obs") < -0.5){
		map[,i] = 1 - map[,i]
	}
}
##Get mean and standard deviation
means=rowMeans(map, na.rm=T)
vars=apply(X=map, FUN=var, MARGIN=1, na.rm=T)
stdev=sqrt(vars)
##Take markers in the inner 95% of variance
probs=pnorm(vars, mean=mean(vars), sd=sd(vars))
tokeep = probs >=0.025 & probs <=0.975

#Reduce map and output order
snps=snps[tokeep]
means=means[tokeep]
snps=snps[order(means)]
means=sort(means)
##Alter hapmap
hmp=hmp[hmp[,1] %in% snps,]
hmp=hmp[match(snps, table=hmp[,1]),]
hmp$pos = round(means * 1e4)
hmp$chrom=mychr
##Output
write.table(hmp, file=outfile, row.names=F, col.names=T, quote=F, sep="\t")