#! /bin/bash

#Identify highly het sites and filter out

TASSEL="/home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx8g"

pop=$1
popdir=$2
cutoff=$3

#Get genotype summary
$TASSEL -fork1 -h $popdir/2a_${pop}_snps_goodmaf_renumber.hmp.txt -genotypeSummary site -export $popdir/2b_${pop}_snps_goodmaf_renumber.sitesummary.txt -runfork1

#Identify and export bad sites (guess-and-check to find appropriate cutoff)
Rscript 2b_FindHighlyHetSites.r $popdir/2b_${pop}_snps_goodmaf_renumber.sitesummary.txt $popdir/2b_${pop}_snps_goodmaf_renumber.badhet_sites.txt $cutoff

#Filter out bad sites
$TASSEL -fork1 -h $popdir/2a_${pop}_snps_goodmaf_renumber.hmp.txt -filterAlign -excludeSiteNamesInFile $popdir/2b_${pop}_snps_goodmaf_renumber.badhet_sites.txt -export $popdir/2c_${pop}_snps_goodmaf_goodsites.hmp.txt -runfork1
