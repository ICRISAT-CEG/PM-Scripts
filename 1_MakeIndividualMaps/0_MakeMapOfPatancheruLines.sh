#! /bin/bash

#Turn called SNPs into a genetic map
#Note: this script is set up so that every step that requires user-mediated input from a previous step (that is, have to manually examine one step to proceed to the next) is separate, so
#      don't run this entire script at once unless you're really sure of what you're doing

pop="patancheru"
popdir="Patancheru"
parent1=841b
parent2=863b
maxprocs=7
keyfile=../0a_keyfile_for_refseq0.2_841_pop_only.txt

##Setup directory and taxa list
if [ ! -d $popdir ]; then mkdir $popdir; fi
##Get list of taxa for this set and filter out of master file
cut --fields=4 $keyfile | tail -n +2 | tr '[:upper:]' '[:lower:]' | sed s/blank// > $popdir/0_${pop}_taxa.txt


###Initial VCF Filtering on taxa and depth
mindepth=1	#Minimum numberof reads for a snp call to be conisdered quality; with deeped coverage, can set this higher
./1_FilterVcfFilesComponent.sh $pop $popdir $mindepth

#####Filter on call quality
goodGCs="66,67,79,80,88,89,94,96,97,98,99"	#Acceptable GCs to include. MUST CONSULT RESULTS OF SCRIPT 1A TO DETERMINE	-> Note, not sure of should keep 99 for this pop. ~40% hets
perl 1b_ApplyAdHocVcfFilters.pl $popdir/1_${pop}.recode.vcf $popdir/1b_${pop}.recode.filtered.vcf $goodGCs


###Filtering on MAF, coverage, and heterozygosity
minMAF=0.25
maxMAF=0.75
minCov=60/100	#Coverage as fraction b/c bash doesn't do floating point arithmetic
./2_FilterOnMafAndGetLdComponent.sh $pop $popdir $minMAF $maxMAF $minCov


###Filter out highly heterozygous states
cutoff=0.12	#Must be found iteratively to get acceptable cutoff
./2b_FilterOutHighlyHetSitesComponent.sh $pop $popdir $cutoff

####Filter out potential outcrosses and highly missing taxa; last step is to strip hets
hetcutoff=0.1
missingcutoff=0.5
./2c_IdentifyAndRemoveOutcrossesComponent.sh $pop $popdir $hetcutoff $missingcutoff $parent1 $parent2


#Cluster SNPs into linkage groups
clusters=$(seq --separator=, 6 16)	#Set range of clusters to investigate
echo "Clustering into clusters of $clusters"
./3_GroupMarkersByClusteringComponent.sh $pop $popdir $clusters

###################
#NOTE: This is the point you really need to stop and look at the clusterings to figure out what goes with what
###################

#####Impute and bootstrap chromosome maps with MSTmap and create merged composites of each LG
best_cluster=13
groups="1,2_5,3,6_11,7,8_13,10_12"	#Clusters to join. Commas separate linkage groups, underscores separate clusters within the groups
nbootstraps=100
nomap_dist=30
nomap_size=0
./3b_BootstrapMapsComponent.sh $pop $popdir $best_cluster $groups $parent1 $parent2 $nbootstraps $nomap_dist $nomap_size

#Run through R/qtl to ripple order and fine-tune
winsize=6
maxprocs=3	#warning: setting this too high can quickly run out of memory
./3g_RippleMapsComponent.sh $pop $popdir $winsize $maxprocs

##Anchor lower-quality SNPs 
minMAF=0.25	#minimum and maximum minor allele freq to a SNP to be taken to anchor
maxMAF=0.75
ncheck=100	#Number of sites to spot-check
minRsq=0.6	#Bins in visual output of SNP LD
hetcutoff=0.1
min_non_het_count=50
vcfsource=../06_HapMap/refseq11_unfiltered.vcf
./4_AnchorSnpsOnCoreMapComponent.sh $pop $popdir $minMAF $maxMAF $ncheck $minRsq $hetcutoff $min_non_het_count $vcfsource


##Bootstrap, impute, and remake map
startchrom=1
endchrom=7
nbootstraps=100
nomap_dist=30
nomap_size=0
./4i_BootstrapImputeAndOrderFullMapComponent.sh $pop $popdir $parent1 $parent2 $startchrom $endchrom $nbootstraps $nomap_dist $nomap_size $maxprocs

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