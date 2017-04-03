#! /bin/bash

#Make tags-by-taxa file of raw data, convert to text, and then map in R

TASSEL="perl /home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms2g -Xmx12g"

pop=$1
popdir=$2
tbtfile=$3
pval=$4
ncores=$5
nToGraph=$6
mintaxa=$7

##Map tags
echo "Mapping tags with ${ncores} cores. WARNING: This will probably take a long time"
Rscript 5_AnchorTagsToMap.r $popdir/3i_${pop}_rippled_combined.hmp.txt $tbtfile $popdir/5_${pop}_anchored_tags $pval $ncores $mintaxa

#Output graphic of random sample
##Get random subset
perl 5a_RandomlySubsetFile.pl $popdir/5_${pop}_anchored_tags_all.txt $nToGraph header $popdir/5a_${pop}_anchored_tags_all_subset.txt
Rscript 5a_GraphSubsetOfAllTags.r $popdir/5a_${pop}_anchored_tags_all_subset.txt $popdir/5a_${pop}_anchored_tags_all_subset.png
