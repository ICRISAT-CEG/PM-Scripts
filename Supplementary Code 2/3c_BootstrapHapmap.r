#Bootstrap a hapmap for MSTmap ordering

#Arguments: (1) input file, (2) Output file, (3+) Names of lines to exclude from bootstrapping and add as-is (eg, parental lines)
options(stringsAsFactors=F)

args=commandArgs(T)
infile=args[1]
outfile=args[2]
n=args[3]
toKeep=args[-(1:3)]

cat("Bootstrapping",infile,n,"times\n")

#Load in data and identify lines to bootstrap
orig=read.delim(infile, header=T)
header=scan(infile, nlines=1, what=character())
names(orig)=header
datacols=12:NCOL(orig)
toExclude = which(names(orig) %in% toKeep)
datacols = datacols[!datacols %in% toExclude]

#Perform bootstraps
for(i in 1:n){
	out=sub(outfile, pattern="+", replacement=i, fixed=T)
	new.hmp=orig
	new.lines = sample(datacols, size=length(datacols), replace=T)
	new.hmp[,datacols] = orig[,new.lines]
	names(new.hmp)[datacols] = names(orig)[new.lines]
	write.table(new.hmp, file=out, sep="\t", row.names=F, col.names=T, quote=F)
}