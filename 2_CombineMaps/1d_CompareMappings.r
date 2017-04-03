#! /usr/bin/rscript

#Compare the arrangements of up to 3 linkage maps with scaffolds
#Arguments: (1) output file, (2) comma-separated names of populations, (3) alpha transparency for lines, (4+) input files


options(stringsAsFactors=F)
args=commandArgs(T)
#setwd("/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet/20140324_AlignToRefseqV0.3/ConsolidateGeneticMaps/")
#args=c("99_map_comparison.png", "som,841", "4b_som_scaffold_bins.txt", "4b_841_scaffold_bins.txt")
outfile=args[1]
names=unlist(strsplit(args[2], split=","))	#Separate names
alpha=as.numeric(args[3])
infiles = args[-(1:3)]
cat("Comparing maps from",infiles,"\n")

#Load files
maps = list()
scaffolds = character()
for(i in 1:length(infiles)){
	maps[[i]] = read.delim(infiles[i])
	scaffolds = unique(c(scaffolds, maps[[i]]$scaffold))
}

#Create master
master = data.frame(scaffolds=scaffolds)
ylim=0
lg = 0
for(i in 1:length(maps)){
	lg_name = paste("lg", i,sep="")
	binname = paste("binval", i,sep="")
	mymatch = match(scaffolds, maps[[i]]$scaffold)
	master[lg_name] = as.numeric(maps[[i]]$lg[mymatch])
	master[binname] = as.numeric(maps[[i]]$binval[mymatch])
	ylim=range(c(ylim, master[binname]), na.rm=T)
	lg = unique(c(lg, master[,lg_name]))
}

##Calculate colors
ncolors = length(na.omit(lg))
colorvals = rainbow(ncolors, alpha=alpha)
colors = rep(NA, nrow(master))
for(i in 1:length(maps)){
	tofix = is.na(colors)
	lg_name = paste("lg", i,sep="")
	colors[tofix] = colorvals[master[tofix,lg_name]]
}


#plot graphic
xlim=c(1,length(maps))
png(outfile, height=500, width=300)
	plot(NA, NA, main="LG comparison", xlim=xlim, ylim=ylim)
	for(x in 2:max(xlim)){
		y0_name = paste("binval", x-1,sep="")
		y1_name = paste("binval", x,sep="")
		x0 = rep((x-1), nrow(master))
		x1 = rep(x, nrow(master))
		segments(x0=x0, x1=x1, y0=master[,y0_name], y1=master[,y1_name], col=colors)
	}
dev.off()