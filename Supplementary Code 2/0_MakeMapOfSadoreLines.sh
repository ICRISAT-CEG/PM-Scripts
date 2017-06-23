#! /bin/bash

#Turn called SNPs into a genetic map
#Note: this script is set up so that every step that requires user-mediated input from a previous step (that is, have to manually examine one step to proceed to the next) is separate, so
#      don't run this entire script at once unless you're really sure of what you're doing

pop="sadore"
popdir="Sadore"
parent1=3301
parent2=3304
maxprocs=8
keyfile=../0a_keyfile_for_refseq1.1.txt

##Setup directory and taxa list
if [ ! -d $popdir ]; then mkdir $popdir; fi
##Get list of taxa for this set and filter out of master file
grep "^C2C4LACXX" $keyfile | cut --fields=4 | tail -n +2 | tr '[:upper:]' '[:lower:]' | sed s/blank// > $popdir/0_${pop}_taxa.txt


##Initial VCF Filtering on taxa and depth
mindepth=2	#Minimum numberof reads for a snp call to be conisdered quality
./1_FilterVcfFilesComponent.sh $pop $popdir $mindepth

#####Filter on call qualitylsls
goodGCs="67,79,80,88,89,94,96,97,98,99"	#Acceptable GCs to include. MUST CONSULT RESULTS OF SCRIPT 1A TO DETERMINE	-> Note, not sure of should keep 99 for this pop. ~40% hets
perl 1b_ApplyAdHocVcfFilters.pl $popdir/1_${pop}.recode.vcf $popdir/1b_${pop}.recode.filtered.vcf $goodGCs


###Filtering on MAF, coverage, and heterozygosity
minMAF=0.25
maxMAF=0.75
minCov=50/100	#Coverage as fraction b/c bash doesn't do floating point arithmetic
./2_FilterOnMafAndGetLdComponent.sh $pop $popdir $minMAF $maxMAF $minCov


###Filter out highly heterozygous states
cutoff=0.4	#Must be found iteratively to get acceptable cutoff
./2b_FilterOutHighlyHetSitesComponent.sh $pop $popdir $cutoff

####Filter out potential outcrosses and highly missing taxa; last step is to strip hets
hetcutoff=0.55
missingcutoff=0.7
./2c_IdentifyAndRemoveOutcrossesComponent.sh $pop $popdir $hetcutoff $missingcutoff $parent1 $parent2 


#Cluster SNPs into linkage groups
clusters=$(seq --separator=, 6 16)	#Set range of clusters to investigate
echo "Clustering into clusters of $clusters"
./3_GroupMarkersByClusteringComponent.sh $pop $popdir $clusters

####################
#NOTE: This is the point you really need to stop and look at the clusterings to figure out what goes with what
####################

####Impute and bootstrap chromosome maps with MSTmap and create merged composites of each LG
best_cluster=8
groups="1,2,3,4,5,7,8"	#Clusters to join. Commas separate linkage groups, underscores separate clusters within the groups
nbootstraps=100
nomap_dist=30
nomap_size=0
./3b_BootstrapMapsComponent.sh $pop $popdir $best_cluster $groups ${parent1}d ${parent2}d $nbootstraps $nomap_dist $nomap_size	#Arbitrarily choose parent d to represent that parent (each was replicated 5 times)

#Run through R/qtl to ripple order and fine-tune
winsize=6	#warning: setting these too high can quickly run out of memory
maxprocs=3	
./3g_RippleMapsComponent.sh $pop $popdir $winsize $maxprocs

##Anchor lower-quality SNPs 
minMAF=0.25	#minimum and maximum minor allele freq to a SNP to be taken to anchor
maxMAF=0.75
ncheck=100	#Number of sites to spot-check
minRsq=0.6	#Bins in visual output of SNP LD
hetcutoff=0.2
min_non_het_count=100
vcfsource=../06_HapMap/refseq11_unfiltered.vcf
./4_AnchorSnpsOnCoreMapComponent.sh $pop $popdir $minMAF $maxMAF $ncheck $minRsq $hetcutoff $min_non_het_count $vcfsource


##Bootstrap, impute, and remake map
startchrom=1
endchrom=7
nbootstraps=100
nomap_dist=30
nomap_size=0
./4i_BootstrapImputeAndOrderFullMapComponent.sh $pop $popdir ${parent1}b ${parent2}b $startchrom $endchrom $nbootstraps $nomap_dist $nomap_size $maxprocs

#Anchor sequencing tags
tbtfile=../05_TBT/05b_MergedTBT/refseq03.tbt.txt
pval=0.0001
nToGraph=100
mintaxa=10
./5_MapSequencingTagsComponent.sh $pop $popdir $tbtfile $pval $maxprocs $nToGraph $mintaxa

#Extract groupings and give confidence to each (based on % coverage)
chromprefix="chr"
scaffold_key="/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v0.1/0_pseudochromosome_scaffold_key.txt"
samfile="../03_SAM/refseq11_tags.sam"
./6_AssignScaffoldsToLinkageGroupsComponent.sh $pop $popdir $chromprefix $scaffold_key $samfile