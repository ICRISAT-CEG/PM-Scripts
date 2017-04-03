#! /usr/bin/Rscript
#Use an R package for solving the Traveling Salesman Problem to try to find optimal grouping and ordering of scaffolds
#Args: (1) Input "distance" file, (2) Output file prefix

options(stringsAsFactors=F)
args=commandArgs(T)
#setwd("/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/BGI_PMiGAP_calls/")
#args=c("1g_outer100_formatted.rdata", "2_outer100_tsp.txt")

library("TSP", lib.loc="/media/STORAGE/Software/R/x86_64-pc-linux-gnu-library/2.13")
infile=args[1]
outprefix=args[2]

#Load data into current enviroment
load(infile)
distance=fixed
rm(fixed)
colnames(distance)=rownames(distance)

TSPdata = TSP(distance)

concorde_path("/home/jgw87/Software/linkage_mapping/concorde/bin/")
t1=proc.time()[3]
tour = solve_TSP(TSPdata, method="concorde", control=list(precision=4))
t2=proc.time()[3]
cat("Concorde requried", t2-t1,"seconds to solve\n")

#Turn into a map
dist.coords = cbind(tour, c(tour[-1], tour[1]))	#Get the coordinates between the distances at each step
between = distance[dist.coords]
map = data.frame(s1=names(tour), s2=c(names(tour)[-1], names(tour)[1]), dist=between)

#Find connections between scaffolds
s1 = sub(map$s1, pattern="\\|.+", repl="")
s2 = sub(map$s2, pattern="\\|.+", repl="")
same.scaff = s1==s2
joints=map[!same.scaff,]
goodjoints = joints[joints$dist<0.1,]	#Rsq >0.9
medjoints = joints[joints$dist>=0.1 & joints$dist <0.2,] #Rsq from 0.8 to 0.9

#Check how compares to my current map, namely, how often are linked scaffolds on same chr and at same pos
linkmap = read.delim("/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet/20140527_AlignToRefseqV1.1/ConsolidateMaps/7e_expanded_map2_som_841_renumbered.txt")
linkmap$snum = sub(linkmap$scaffold, pattern="scaffold", repl="")

check.match=function(x, linkmap){
	s1match=match(sub(x$s1, pattern="\\|.+", repl=""), linkmap$snum)
	s2match=match(sub(x$s2, pattern="\\|.+", repl=""), linkmap$snum)
	#Check chromosome
	chrmatch = linkmap$lg[s1match] == linkmap$lg[s2match]
	cat(sum(chrmatch, na.rm=T) / sum(is.finite(chrmatch)),"of linked hits are on the same chromosome\n")
	binmatch = abs(linkmap$pos[s1match] - linkmap$pos[s2match])
	cat(sum(binmatch<10, na.rm=T)/ sum(is.finite(binmatch)), "of hits on same chr are within 10 map units\n")
	hist(binmatch[chrmatch])
}


########Tests########
##Results
#nn = 0.1 sec; final length 771
#repetative_nn = 87.4 sec; final length 768
#nearest_insertion = 6.5 sec; final length 768
#cheapest_insertion = 6.7 sec; final length 763
#2-opt = 17.9 sec; final length 757
#concorde = 5.1 sec; final length 752	- Seems the best, and time is comprable to linkern (though that may change as nodes scales up)
#linkern = 4.8 sec; final length 754
concorde_path("/home/jgw87/Software/linkage_mapping/concorde/bin/")
subd = distance[1:1000, 1:1000]
TSPdata=TSP(subd)
results=list()
mymethods=c("nn","repetitive_nn", "nearest_insertion", "cheapest_insertion", "2-opt","concorde", "linkern")
for(method in mymethods){
	t1=proc.time()[3]
	tour = solve_TSP(TSPdata, method=method, control=list(precision=4))
	results[[method]] = tour
	t2=proc.time()[3]
	cat(method,"requried", t2-t1,"seconds to solve\n")
}

##Try different methods, store distances?