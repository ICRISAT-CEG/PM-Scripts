#Analyze map by R/qtl

#setwd("/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet_SomSekhar/Jason_Analysis/VCF_Filter2/")
library("snow", lib.loc="/media/STORAGE/Software/R/x86_64-pc-linux-gnu-library/2.13")
library("qtl", lib.loc="/media/STORAGE/Software/R/x86_64-pc-linux-gnu-library/2.13")
args=commandArgs(T)
infile=args[1]
outprefix=args[2]
winsize=as.numeric(args[3])

cat("Reading in genotypes to ripple from ",infile,"\n")

cross=read.cross(file=infile, format="csvr", sep="\t")
cross=est.rf(cross)
newmap=est.map(cross, map.function="kosambi", verbose=T, error.prob=0.01)
cross=replace.map(cross, newmap)

##########################
# Reordering linkage groups using the permutations from the ripples with the crossover count objectives - Script from Eli
##########################

findSwitches<-function(ripMat) {
        initial.xO <- ripMat[1,'obligXO']
        init.vec <- 1:(ncol(ripMat)-1)
        # Get the permutations with number of crossovers less than initial
        better.rows <- which(as.vector(ripMat[,'obligXO']) < initial.xO)
        wind <- attr(ripMat,"window")
        perm.mat <- matrix(0,nrow=length(better.rows),ncol=(wind+1))
        orig.mat <- matrix(0,nrow=length(better.rows),ncol=wind)
        count <- 1
        # Keep track of what's been switched - only keep the first instance of a permutation that includes a new marker
        switched <- c()
        for (i in better.rows)
        {
                test.vec <- as.vector(ripMat[i,1:(ncol(ripMat)-1)])
                perm <- test.vec[which(test.vec != init.vec)]
                if (!any(switched %in% perm)) {
                        perm.mat[count,1:length(perm)] <- perm
                        orig.mat[count,1:length(perm)] <- which(init.vec != test.vec)
                        perm.mat[count,(wind+1)] <- ripMat[i,ncol(ripMat)]
                        switched <- c(switched,perm[which(perm != 0)])
                        count <- count + 1
                }
        }
        keep.rows <- which(as.vector(perm.mat[,1]) != 0)
        perm.mat <- perm.mat[keep.rows,]
        orig.mat <- orig.mat[keep.rows,]
		if(class(perm.mat)!= "matrix"){
			perm.mat=matrix(perm.mat, nrow=1)
		}
		if(class(orig.mat)!= "matrix"){
			orig.mat=matrix(orig.mat, nrow=1)
		}
        return(list(perm=perm.mat,pos=orig.mat))
}
getSwitchVec<-function(switch.list,ripMat)
{
        init.vec <- 1:(ncol(ripMat)-1)
         if(nrow(switch.list$pos) == 0) return(init.vec)
        for (i in 1:nrow(switch.list['perm']$perm))
        {
                pos <- as.vector(switch.list['pos']$pos[i,which(as.vector(switch.list['pos']$pos[i,]) != 0)])
                init.vec[pos] <- as.vector(switch.list['perm']$perm[i,which(as.vector(switch.list['pos']$pos[i,]) != 0)])
        }
        return(init.vec)
}

######
#End Eli's functions
######

#Wrapper to make things easy
get.new.map=function(cross, window){
	ripples=list()
	for (chr in chrnames(cross)) {
		ripples[[chr]] <- ripple(cross,window=window,chr=chr,method='countxo',error.prob=0.01)
		switch.list <- findSwitches(ripples[[chr]])
		switch.vec <- getSwitchVec(switch.list,ripples[[chr]])
		cross <- switch.order(cross,chr=chr,order=switch.vec,error.prob=0.01,maxit=10000)
	}
	cross=est.rf(cross)
	newmap=est.map(cross, map.function="kosambi", verbose=T, error.prob=0.01)
	cross=replace.map(cross, newmap)
	return(cross)
}

cat("Rippling orders in window size",winsize,"(note: this may take a long time and a lot of memory)\n")
cross.new=get.new.map(cross=cross, window=winsize)


cat("Finished rippling. Outputting to",outprefix,"\n")

write.cross(cross.new,format="csvr",filestem=outprefix)
png(paste(outprefix, "_map.png", sep=""), width=200, height=400)
	plot.map(cross)
dev.off()

