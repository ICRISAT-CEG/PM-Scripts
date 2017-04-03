#! /usr/bin/Rscript

#Calculate total map length and output to screen and to a file
#Args: (1) Input file, (2) TSP-formatted distance file, (3) Output file

options(stringsAsFactors=F)
args=commandArgs(T)
#setwd("/home/jgw87/Desktop/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/BGI_PMiGAP_calls/")
#args=c("2b_new_map_from_outer100.txt","2_outer100_formatted_for_tsp.txt","2c_new_map_from_outer100_lengths.txt")

infile=args[1]
distfile=args[2]
outfile=args[3]

cat("Loading data...\n")
map=read.delim(infile)
	map=map[order(map$lg, map$bin),]
distance=NULL
	firstrow = scan(distfile, what=character(), nlines=1)
	colClasses=rep("numeric", length(firstrow))
	colClasses[1]="character"
	distance=as.matrix(read.delim(distfile, colClasses=colClasses, header=F, row.names=1))
	colnames(distance)=rownames(distance)

cat("Determining path\n")
lgs=split(map, f=map$lg)
maplengths=list()
for(i in 1:length(lgs)){
	#Reformat so matches distance matrix
	scaffnum = sub(lgs[[i]]$scaffold, pattern="scaffold", repl="")
	sides=data.frame(row.names=scaffnum, one=rep("left", length(scaffnum)), two="right")
		toflip = lgs[[i]]$orientation =="rev"
		sides[toflip,c("one","two")] = sides[toflip,c("two","one")]
	
	#Find the path through the map and reduce to ones where have distance measures	
	path1=paste(scaffnum, sides$one, sep="|")
	path2=paste(scaffnum, sides$two, sep="|")
	path=paste(path1, path2, sep=",", collapse=",")
	path=unlist(strsplit(path, split=","))
	path = path[path %in% rownames(distance)]
	
	#Get coordinates
	rowcoord = match(path, rownames(distance))
	towalk = cbind(rowcoord[-length(rowcoord)], rowcoord[-1])
	
	#Calc distances
	mydist = distance[towalk]
	maplengths[[i]] = data.frame(lg=names(lgs)[i], mapped=length(path), length=sum(mydist))
}

maplengths= do.call(what=rbind, args=maplengths)
maplengths=maplengths[order(maplengths$lg),]
maplengths=rbind(maplengths, data.frame(lg="total", mapped=sum(maplengths$mapped), length=sum(maplengths$length)))

cat("Map lengths from",infile,"are:\n")
print(maplengths)
cat("\n(smaller is better)\n")

write.table(maplengths, file=outfile, sep="\t", quote=F, row.names=F, col.names=T)