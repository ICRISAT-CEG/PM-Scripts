#! /bin/bash

#First step of map-making: Filter VCF values for taxa and remove unreliable sites

VCFTOOLS="/home/jgw87/Software/Filetools/vcftools_0.1.11/bin/vcftools"

pop=$1
popdir=$2
mindepth=$3

#Filter for just the required taxa
$VCFTOOLS --vcf ../06_HapMap/refseq11_unfiltered.vcf --keep $popdir/0_${pop}_taxa.txt --out $popdir/1_${pop} --recode --recode-INFO-all --minDP $mindepth --maf 0.001 --max-maf 0.999	#Filter out any sites with MAF<0.001 b/c are basically monomorphic

#Graph distribution of quality scores and heterozygosity in order to make decision about what sort of quality filter to put on the calls
perl 1a_ExtractVcfSiteInfo.pl $popdir/1_${pop}.recode.vcf $popdir/1a_${pop}_vcf_qualities.txt
Rscript 1a_GraphVcfQualityDistributions.r $popdir/1a_${pop}_vcf_qualities.txt $popdir/1a_${pop}_vcf_qualities_parsed 
