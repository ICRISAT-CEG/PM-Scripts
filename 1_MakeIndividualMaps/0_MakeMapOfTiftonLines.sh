#! /bin/bash

#Turn called SNPs into a genetic map
#Note: this script is set up so that every step that requires user-mediated input from a previous step (that is, have to manually examine one step to proceed to the next) is separate, so
#      don't run this entire script at once unless you're really sure of what you're doing

pop="tifton"
popdir="Tifton"
parent1=99_17_1
parent2=99b
maxprocs=7
keyfile=/media/STORAGE/Working_Files/GBS/Flowcells/Som/Som_keyfile.txt

#Setup directory and taxa list
if [ ! -d $popdir ]; then mkdir $popdir; fi
#Get list of taxa for this set and filter out of master file
cut --fields=4 $keyfile | tail -n +2 | tr '[:upper:]' '[:lower:]' | sed s/blank// > $popdir/0_${pop}_taxa.txt


##Initial VCF Filtering on taxa and depth
mindepth=5	#Minimum numberof reads for a snp call to be conisdered quality
./1_FilterVcfFilesComponent.sh $pop $popdir $mindepth

###Filter on call quality
goodGCs="67,80,89,94,96,97,98,99"	#Acceptable GCs to include. MUST CONSULT RESULTS OF SCRIPT 1A TO DETERMINE
perl 1b_ApplyAdHocVcfFilters.pl $popdir/1_${pop}.recode.vcf $popdir/1b_${pop}.recode.filtered.vcf $goodGCs


##Filtering on MAF, coverage, and heterozygosity
minMAF=0.25
maxMAF=0.75
minCov=60/100	#Coverage as fraction b/c bash doesn't do floating point arithmetic
./2_FilterOnMafAndGetLdComponent.sh $pop $popdir $minMAF $maxMAF $minCov


###Filter out highly heterozygous states
cutoff=0.15	#Must be found iteratively to get acceptable cutoff
./2b_FilterOutHighlyHetSitesComponent.sh $pop $popdir $cutoff

#####Filter out potential outcrosses and highly missing taxa; last step is to strip hets
hetcutoff=0.1
missingcutoff=0.5
./2c_IdentifyAndRemoveOutcrossesComponent.sh $pop $popdir $hetcutoff $missingcutoff $parent1 $parent2


##Cluster SNPs into linkage groups
#clusters=$(seq --separator=',' 6 15)
clusters="15,16,17"
echo $clusters
./3_GroupMarkersByClusteringComponent.sh $pop $popdir $clusters

####################
#NOTE: This is the point you really need to stop and look at the clusterings to figure out what goes with what
####################

#####Impute and bootstrap chromosome maps with MSTmap and create merged composites of each LG
best_cluster=15
groups="1_13,2_10,3_4,5_11,6_12,8_14,15"	#Clusters to join. Commas separate linkage groups, underscores separate clusters within the groups
nbootstraps=100
nomap_dist=30
nomap_size=0
./3b_BootstrapMapsComponent.sh $pop $popdir $best_cluster $groups $parent1 $parent2 $nbootstraps $nomap_dist $nomap_size


#Run through R/qtl to ripple order and fine-tune
winsize=7
maxprocs=4
./3g_RippleMapsComponent.sh $pop $popdir $winsize $maxprocs

###Anchor lower-quality SNPs 
minMAF=0.25	#minimum and maximum minor allele freq to a SNP to be taken to anchor
maxMAF=0.75
ncheck=100	#Number of sites to spot-check
minRsq=0.6	#Bins in visual output of SNP LD
hetcutoff=0.05
min_non_het_count=60
vcfsource=../06_HapMap/refseq11_unfiltered.vcf	#Which VCF file to pull the non-core SNPs from
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