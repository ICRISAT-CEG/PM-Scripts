#! /bin/bash

#Filter map files on MAF, coverage, and heterozygosit

TASSEL="/home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx8g"
TASSEL5="/home/jgw87/Software/tassel5-standalone/run_pipeline.pl -Xms4g -Xmx8g"


pop=$1
popdir=$2
minMAF=$3
maxMAF=$4
minCov=$5

#Filter for MAF and coverage
ntaxa=$(head -n 20 $popdir/1b_${pop}.recode.filtered.vcf | tail -n 1 | sed -r -e "s/\s+/\n/g" | wc -l)
ntaxa=$(($ntaxa-9))
mincount=$(($ntaxa * $minCov))
$TASSEL -fork1 -vcf $popdir/1b_${pop}.recode.filtered.vcf -filterAlign -filterAlignMinFreq $minMAF -filterAlignMaxFreq $maxMAF -filterAlignMinCount $mincount -export $popdir/2_${pop}_snps_goodmaf.hmp.txt -runfork1

#Renumber locus and position for easier work (and so can do LD plots)
Rscript 2a_StripLocusInfoFromHapmap.r $popdir/2_${pop}_snps_goodmaf.hmp.txt $popdir/2a_${pop}_snps_goodmaf_renumber.hmp.txt
$TASSEL5 -fork1 -h $popdir/2a_${pop}_snps_goodmaf_renumber.hmp.txt -subsetSites 5000 -ld -ldType All -ldHetTreatment Homozygous -ldd png -o $popdir/2a_${pop}_snps_goodmaf_renumber.ld.png -ldplotsize 2000 -runfork1

