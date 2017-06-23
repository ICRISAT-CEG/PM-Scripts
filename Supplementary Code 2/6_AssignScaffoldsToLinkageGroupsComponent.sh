#! /bin/bash

#Take final orders and create documents 

TASSEL="/home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx8g"

pop=$1
popdir=$2
chromprefix=$3
scaffold_key=$4
samfile=$5

#core=$popdir/3i_${pop}_rippled_combined.hmp.txt
snps=$popdir/4k_${pop}_anchored_map_reordered_combined.hmp.txt
tags=$popdir/5_${pop}_anchored_tags_best.txt
refmap=$popdir/3i_${pop}_rippled_combined.hmp.txt

##SNPs
python3 6a_MapScaffoldsToLinkageGroups.py -i $snps -c $chromprefix -l $scaffold_key -o $popdir/6a_${pop}_snp_scaffold_assignments.txt

##Tags
python3 6b_GetPositionOfTags.py -i $tags -s $samfile -m $snps -o $popdir/6b_${pop}_tag_anchors.txt
python3 6a_MapScaffoldsToLinkageGroups.py -i $popdir/6b_${pop}_tag_anchors.txt -c $chromprefix -l $scaffold_key -o $popdir/6c_${pop}_tag_scaffold_assignments.txt

#Combine and output final assemblies
python3 6d_ConsolidateScaffoldsIntoBins.py -s $popdir/6a_${pop}_snp_scaffold_assignments.txt -t $popdir/6c_${pop}_tag_scaffold_assignments.txt -m $popdir/4k_${pop}_anchored_map_reordered_combined.hmp.txt \
	-r $popdir/6d_${pop}_scaffold_raw_counts.txt -o $popdir/6d_${pop}_scaffold_bins.txt


