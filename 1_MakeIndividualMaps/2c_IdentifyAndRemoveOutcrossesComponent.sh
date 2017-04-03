#! /bin/bash

#Check the number of rare alleles in each line and identify those with an excess of rare alleles, then filter out

TASSEL="perl /home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms2g -Xmx8g"

pop=$1
popdir=$2
hetcutoff=$3
missingcutoff=$4
parent1=$5
parent2=$6
keephets=$7	#optional argument

##Test for rare alleles
$TASSEL -fork1 -vcf $popdir/1_${pop}.recode.vcf -filterAlign -filterAlignMinFreq 0.001 -filterAlignMaxFreq 0.05 -export $popdir/2d_${pop}_rare_alleles.hmp.txt -runfork1
$TASSEL -fork1 -h $popdir/2d_${pop}_rare_alleles.hmp.txt -numericalGenoTransform collapse -export $popdir/2d_${pop}_rare_alleles.numeric.txt -runfork1
Rscript 2d_TestForOutcrossedTaxa_RareAlleles.r $popdir/2d_${pop}_rare_alleles.numeric.txt $popdir/2d_${pop}_potential_outcrosses_rare_alleles 

##Test for excess heterozygosity and missingness
$TASSEL -fork1 -h $popdir/2c_${pop}_snps_goodmaf_goodsites.hmp.txt -genotypeSummary taxa -export $popdir/2e_${pop}_snps_goodmaf_goodsites.taxasummary.txt -runfork1
Rscript 2c_TestForOutcrossedTaxa_HetsAndMissing.r $popdir/2e_${pop}_snps_goodmaf_goodsites.taxasummary.txt $popdir/2e_${pop}_potential_outcrosses_highly_het_or_missing.txt $hetcutoff $missingcutoff

##Remove bad taxa (ignoring parentals)
goodhap=$popdir/2g_${pop}_snps_goodmaf_goodsites_goodtaxa.hmp.txt
cat $popdir/2d_${pop}_potential_outcrosses_rare_alleles.txt $popdir/2e_${pop}_potential_outcrosses_highly_het_or_missing.txt | sed -r -e "s/$parent1//" -e "s/$parent2//" > $popdir/2g_${pop}_badtaxa_combined.txt
$TASSEL -fork1 -h $popdir/2c_${pop}_snps_goodmaf_goodsites.hmp.txt -excludeTaxaInFile $popdir/2g_${pop}_badtaxa_combined.txt -export $goodhap -runfork1

##Strip heterozygous calls
final=$popdir/2h_${pop}_snps_goodmaf_goodsites_goodtaxa_nohets.hmp.txt
if [ $keephets == "keephets" ]; then
	cp $goodhap $final
else
	head -n 1 $goodhap > $final
	tail -n +2 $goodhap | sed -e "s/\t[RYSWKM]/\tN/g" >> $final
fi
