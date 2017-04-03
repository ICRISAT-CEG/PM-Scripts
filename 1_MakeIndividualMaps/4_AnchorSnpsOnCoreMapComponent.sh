#! /bin/bash

#Anchor SNPs based their LD with the core markers

TASSEL="perl /home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms2g -Xmx12g"

pop=$1
popdir=$2
minmaf=$3
maxmaf=$4
ncheck=$5
minRsq=$6
hetcutoff=$7
minnonhetcount=$8
vcfsource=$9

##Get core sites and taxa
coremap=$popdir/3i_${pop}_rippled_combined.hmp.txt
tail -n +2 $coremap | cut --fields=1 > $popdir/4_${pop}_core_sites.txt
head -n 1 $coremap | sed -e "s/\t/\n/g" | tail -n +12 > $popdir/4_${pop}_good_taxa.txt
$TASSEL -fork1 -vcf $vcfsource -filterAlign -excludeSiteNamesInFile $popdir/4_${pop}_core_sites.txt -includeTaxaInFile $popdir/4_${pop}_good_taxa.txt \
	-fork2 -filterAlign -input1 -filterAlignMinFreq $minmaf -filterAlignMaxFreq $maxmaf -export $popdir/4_${pop}_noncore_sites.hmp.txt -runfork1
##Alter names so obvious which ones are anchored and which are core (also for script, which assumes anything starting with S is an anchor tag)
sed "s/^S/target_S/g" $popdir/4_${pop}_noncore_sites.hmp.txt > $popdir/4a_${pop}_noncore_sites_renamed.hmp.txt
#Recombine
$TASSEL -fork1 -h $coremap -fork2 -h $popdir/4a_${pop}_noncore_sites_renamed.hmp.txt -combine3 -input1 -input2 -intersect -export $popdir/4a_${pop}_all_sites_tomap.hmp.txt -runfork1 -runfork2 -runfork3

#Run custom anchoring java script
##Find best position for the target SNPs; again, note that the script expect to anchor all SNPs whose names don't start with "S" against those that do.
java -cp 4b_JavaScriptComponents/JasonScripts.jar Misc.FindRsqAlongChrom $popdir/4a_${pop}_all_sites_tomap.hmp.txt $popdir/4_${pop}_core_sites.txt $popdir/4b_${pop}_rsq_along_chrom.txt
#Visually spot-check assignment of some SNPs
Rscript 4b_ProcessRsqAlongChrom.r $popdir/4b_${pop}_rsq_along_chrom.txt $ncheck $minRsq $popdir/4c_${pop}_best_pos_for_targets

#Shuffle into map
anchored=$popdir/4d_${pop}_anchored_map.hmp.txt
perl 4d_ReorderTargetsIntoMap.pl $popdir/4a_${pop}_all_sites_tomap.hmp.txt $popdir/4c_${pop}_best_pos_for_targets_goodsites.txt $popdir/4c_${pop}_best_pos_for_targets_badsites.txt $anchored

#Filter out low-quality sites again (based on hets and % missing)
summary=$popdir/4e_${pop}_anchored_map.sitesummary.txt
badsites=$popdir/4e_${pop}_anchored_map.badsites.txt
filtered=$popdir/4f_${pop}_anchored_map_goodsites.hmp.txt
$TASSEL -fork1 -h $anchored -genotypeSummary site -export $summary -runfork1
Rscript 2b_FindHighlyHetSites.r $summary $badsites $hetcutoff
$TASSEL -fork1 -h $anchored -filterAlign -excludeSiteNamesInFile $badsites -export $filtered -runfork1

###Strip heterozygous calls
nohets=$popdir/4g_${pop}_anchored_map_goodsites_nohets.hmp.txt
newfilt=$popdir/4h_${pop}_anchored_map_goodsites_nohets_filt.hmp.txt
head -n 1 $filtered > $nohets
tail -n +2 $filtered | sed -e "s/\t[RYSWKM]/\tN/g" >> $nohets
$TASSEL -fork1 -h $nohets -filterAlign -filterAlignMinCount $minnonhetcount -export $newfilt -runfork1

