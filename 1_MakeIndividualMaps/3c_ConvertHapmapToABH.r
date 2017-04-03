# Convert hapmap to Rqtl

#Load data and set up
options(stringsAsFactors=F)
args=commandArgs(T)
infile=args[1]
outfile=args[2]
parent1=args[3]
parent2=args[4]
hmp=read.delim(infile, header=T)
header=scan(infile, nlines=1, what=character())
names(hmp)=header
datacols=12:NCOL(hmp)	#columns with data
hets=c("S","W","R","Y","M","K","0")
hets.perl=paste(c("[",hets,"]"), collapse="")	#Expression for finding hets
cat("Converting",infile,"to A-B-H format\n")

#Determine which columns have the parents
p1=which(names(hmp)==parent1)
p2=which(names(hmp)==parent2)

#Impute parent lines (on the logic that if one is known and the other isn't, the unknown is the allele the other is not
new=hmp
problems=new[,p1]==new[,p2] | ((grepl(new[,p1], pattern=hets.perl, perl=T) & new[,p2]=="N") | ( new[,p1]=="N") & grepl(new[,p2], pattern=hets.perl, perl=T))	#Vector to keep track of problem sites where can't assign a parent to the genotypes; checks for parents with same calls or one with het call and other with missing
alleles=do.call(what=rbind, args=strsplit(hmp$alleles, split="/"))
for(r in 1:NROW(new)){
	if(problems[r]){	#If already flagged as a problem site, skip
		next;
	}
	if(new[r,p1] %in% alleles[r,] && new[r,p2] %in% alleles[r,]){	#If both are in the allele set, no need to impute
		next;
	}
	if(new[r,p1]=="N"){		#Impute Ns
		new[r,p1]=alleles[r,which(alleles[r,] != new[r,p2]) ]
	}else if(new[r,p2]=="N"){	
		new[r,p2]=alleles[r,which(alleles[r,] != new[r,p1]) ]
	}else if(new[r,p1] %in% hets){		#Fix erran het calls
		new[r,p1]=alleles[r,which(alleles[r,] != new[r,p2]) ]
	}else if(new[r,p2] %in% hets){
		new[r,p2]=alleles[r,which(alleles[r,] != new[r,p1]) ]
	}
}


#Go through data and alter for Rqtl
for(mycol in datacols){
	if(mycol==p1 || mycol==p2){next;}	#skip the parental columns
	new[,mycol] = ifelse(new[,mycol]=="N", yes="-", no=new[,mycol])	#Change missing
	new[,mycol] = ifelse(grepl(new[,mycol], pattern=hets.perl, perl=T), yes="H", no=new[,mycol])	#Change hets
	new[,mycol] = ifelse(new[,mycol]==new[,p1], yes="p1", no=new[,mycol])	#Match parent A (temp coding to prevent errors when parent B has Adenine
	new[,mycol] = ifelse(new[,mycol]==new[,p2], yes="p2", no=new[,mycol])	#Match parent B
}

#Replace temp codings with A and B
for(mycol in datacols){
	new[,mycol] = ifelse(new[,mycol]=="p1", yes="A", no=new[,mycol])
	new[,mycol] = ifelse(new[,mycol]=="p2", yes="B", no=new[,mycol])
}

new[,p1]="A"
new[,p2]="B"

write.table(new[!problems,], file=outfile, sep="\t", quote=F, row.names=F, col.names=T)