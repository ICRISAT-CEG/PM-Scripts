#! /usr/bin/Rscript
#Reorder the linkage map, first ordering by Traveling Salesman within each linkage bin and then ordering across linkage bins
#Args: (1) Linkage map file, (2) Input "distance" file, (3) Output map file with new order and orientation, (4) Directory in which concorde executable is found

options(stringsAsFactors=F)
args=commandArgs(T)
#setwd("/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/BGI_PMiGAP_calls/")
#args=c("/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet/20140527_AlignToRefseqV1.1/ConsolidateMaps/7e_expanded_map2_som_841_renumbered.txt", "1g_outer100_formatted.txt", "2a_map_reordered_from_outer100_tsp.txt",
#	"/media/STORAGE/Software/linkage_mapping/concorde/bin")

library("TSP", lib.loc="/media/STORAGE/Software/R/x86_64-pc-linux-gnu-library/2.13")
library("parallel", lib.loc="/media/STORAGE/Software/R/x86_64-pc-linux-gnu-library/2.13")
mapfile=args[1]
distfile=args[2]
outfile=args[3]
concorde.dir=args[4]
ncores=8

#Load data intocurrent enviroment and clean up for analysis
cat("Loading data for traveling salesman reordering\t")
map=read.delim(mapfile)
	map$all_anchors = map$best_anchor = map$anchor = NULL
firstrow = scan(distfile, what=character(), nlines=1)
	colClasses=rep("numeric", length(firstrow))
	colClasses[1]="character"
	distance=as.matrix(read.delim(distfile, colClasses=colClasses, header=F, row.names=1))
	colnames(distance)=rownames(distance)
concorde_path(concorde.dir)

#Split into individual bins and do TSP in parallel
cat("Performing traveling salesman reordering\t")
map$scaffold=sub(map$scaffold, pattern="scaffold", repl="")
maplist = split(map, f=list(map$lg, map$pos))
scaffold.names=sub(rownames(distance), pattern="\\|.+", repl="")

#Function wrapped to actually do TSP
do.tsp=function(mymap, mydist, scaffnames){
	if(nrow(mymap)==0){
		return(data.frame(scaffold1=character(), scaffold2=character(), distance=numeric(), lg=numeric(), pos=numeric()))
	}
	if(nrow(mymap)==1){	#Skip if map only has 1 scaffold
		newmap = data.frame(scaffold1=paste(mymap$scaffold,"left", sep="|"), scaffold2=paste(mymap$scaffold,"right", sep="|"), distance=0, lg=mymap$lg, pos=mymap$pos)
		return(newmap)
	}
	#Determine which scaffolds to keep
	tokeep = which(scaffnames %in% mymap$scaffold)
	subdist = mydist[tokeep, tokeep]
	if(class(subdist)!="matrix" || nrow(subdist)< 3){	#Need at least three nodes to do a TSP problem. 0 nodes makes a "numeric" object
		#Handle zero-row exceptions
		newmap = data.frame(scaffold1=paste(mymap$scaffold,"left", sep="|"), scaffold2=paste(mymap$scaffold,"right", sep="|"), distance=0, lg=mymap$lg, pos=mymap$pos)
		return(newmap)
	}
	mytsp = TSP(subdist)
	tour=solve_TSP(mytsp, method="concorde", control=list(precision=4) )
	
	#Find distances between each and break at the largest distance (= least LD)
	dist.coords = cbind(tour, c(tour[-1], tour[1]))	#Get the coordinates between the distances at each step
	between = subdist[dist.coords]
	tobreak = which.max(between)
	if(tobreak != length(tour)){	#If the break isn't the last point, remake. Otherwise just pass as is.
		new.order = c(tour[(tobreak+1):length(tour)], tour[1:tobreak])	#Break the loop at 'tobreak'
	}else{
		new.order=tour
	}
	
	#Create the new map & return it
	scaff1.coords = new.order[-length(new.order)] #only interested in the linear order, so leave off last-to-first loop
	scaff2.coords = new.order[-1]
	new.dist.coords = cbind(scaff1.coords, scaff2.coords)	
	new.between = subdist[new.dist.coords]
	newmap = data.frame(scaffold1 = rownames(subdist)[scaff1.coords], scaffold2 = rownames(subdist)[scaff2.coords], distance=new.between)
	newmap$lg = mymap$lg[1]	#Take first lg and position, which should all match anyway
	newmap$pos = mymap$pos[1]
	return(newmap)
}
newmaps = mclapply(maplist, FUN=do.tsp, mydist=distance, scaffnames=scaffold.names, mc.cores=ncores)


#Condense map to just scaffold order and position
cat("Condensing map")
condense.map=function(submap){
	if(nrow(submap) == 0){
		return(data.frame(lg=numeric(), pos=numeric(), scaffolds=character()))
	}
	path=c(submap$scaffold1, submap$scaffold2[nrow(submap)])
	path=paste(path, collapse=",")
	return(data.frame(lg=submap$lg[1], pos=submap$pos[1], scaffolds=path))
}
#Get the tips of each bin to determine relative orientation
get.bin.tips = function(binmap){
	if(nrow(binmap) == 0){
		return(data.frame(top=character(), bottom=character(), lg=numeric(), pos=numeric()))
	}
	return(data.frame(top=binmap$scaffold1[1], bottom=binmap$scaffold2[nrow(binmap)], lg=binmap$lg[1], pos=binmap$pos[1]))
}

#Make improved, condensed map
condensed.maps = lapply(newmaps, FUN=condense.map)
goodmap = as.data.frame(do.call(what=rbind, args=condensed.maps))
goodmap=goodmap[order(goodmap$lg, goodmap$pos),]
goodmap = goodmap[rowSums(is.na(goodmap))< ncol(goodmap),]	#Remove NA rows


#Find the tips of each bin to order relative to each other
bintips = lapply(newmaps, FUN=get.bin.tips)
bintips = do.call(what=rbind, args=bintips)
badrows=rowSums(is.na(bintips))
bintips = bintips[badrows!=ncol(bintips),]
bintips = bintips[order(bintips$lg, bintips$pos),]
lgs = split(bintips, f=bintips$lg)

get.dist=function(x, distance){
	i1=which(rownames(distance)==x[1])
	i2=which(colnames(distance)==x[2])
	if(length(i1)==0 || length(i2)==0){
		return(Inf)
	}else{
		return(distance[i1, i2])
	}
}

order.lg = function(lg, distance){
	#Order the first two bins relative to each other by testing which two have the minimum distance between them and putting them together
	fwdfwd = get.dist(c(lg$bottom[1], lg$top[2]), distance)
	fwdrev = get.dist(c(lg$bottom[1], lg$bottom[2]), distance)
	revfwd = get.dist(c(lg$top[1], lg$top[2]), distance)
	revrev = get.dist(c(lg$top[1], lg$bottom[2]), distance)
	best = which.min(c(fwdfwd, fwdrev, revfwd, revrev))
	if(best==1){
		#do nothing; current order is best
	}else if(best==2){
		lg[2,c("top", "bottom")] = lg[2,c("bottom","top")]
	}else if(best==3){
		lg[1,c("top", "bottom")] = lg[1,c("bottom","top")]
	}else if(best==4){
		lg[1,c("top", "bottom")] = lg[1,c("bottom","top")]
		lg[2,c("top", "bottom")] = lg[2,c("bottom","top")]
	}else{
		stop("Error: Best fit",best,"is outside permitted range of 1-4")
	}
	
	for(n in 3:nrow(lg)){	#Go through and order subsequent ones
		nextfwd = get.dist(c(lg$bottom[n-1], lg$top[n]), distance)
		nextrev = get.dist(c(lg$bottom[n-1], lg$bottom[n]), distance)
		if(nextrev < nextfwd){	#Flip if reversing the next bin gives a better distance measure
			lg[n,c("top", "bottom")] = lg[n,c("bottom","top")]
		}
	}
	return(lg)
}
newlg = lapply(lgs, FUN=order.lg, distance=distance)
newlg = do.call(what=rbind, args=newlg)
newlg = newlg[order(newlg$lg, newlg$pos),]


if(any(goodmap$lg != newlg$lg) || any(goodmap$pos != newlg$pos)){
	cat("\tWARNING! Some LG and positions do not match!!")
}


#Flip entire bins if necessary
scaffolds = strsplit(goodmap$scaffolds, split=",")
for(i in 1:length(scaffolds)){
	#If already in right order, skip to next
	if(scaffolds[[i]][1] == newlg$top[i] && scaffolds[[i]][length(scaffolds[[i]])] == newlg$bottom[i]){
		next	
	}
	#Flip if backwards
	if(scaffolds[[i]][1] == newlg$bottom[i] && scaffolds[[i]][length(scaffolds[[i]])] == newlg$top[i]){
		scaffolds[[i]] = rev(scaffolds[[i]])
		next
	}
	
	cat("Error at line",i,": Didn't match forward or reverse!\n")
}
newscaff = lapply(scaffolds, FUN=paste, collapse=",")
goodmap$scaffolds = unlist(newscaff)
write.table(goodmap, file=outfile, sep="\t", quote=F, row.names=F, col.names=T)