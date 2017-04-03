#! /bin/bash

#Convert final hapmap to numeric and then run clustering on it in R to identify clusters

TASSEL="perl /home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms2g -Xmx8g"
TASSEL5="perl /home/jgw87/Software/tassel5-standalone/run_pipeline.pl -Xms2g -Xmx12g"

pop=$1
popdir=$2
clusters=$3

hapmap=$popdir/2h_${pop}_snps_goodmaf_goodsites_goodtaxa_nohets.hmp.txt
numeric=$popdir/3_${pop}_snps_goodmaf_goodsites_goodtaxa_nohets.numeric.txt

#Filter and transform to numeric
$TASSEL -fork1 -h $hapmap -numericalGenoTransform collapse -export $numeric -runfork1

#Cluster in R
Rscript 3_GroupMarkersByClustering.r $numeric $hapmap $clusters $popdir/3a_${pop}_clustered

#Draw LD for each
maxprocs=6
procs=0
for hap in $popdir/3a_${pop}_clustered*.hmp.txt; do
	outpng=${hap/hmp.txt/ld.png}
	$TASSEL5 -fork1 -h $hap -subsetSites 5000 -ld -ldType All -ldHetTreatment Homozygous -ldd png -o $outpng -ldplotsize 2000 -runfork1 &
	procs=$(($procs + 1))
	if [ $procs -ge $maxprocs  ]; then
		wait
		procs=0
	fi
done

wait
echo -e  "\n\n\n######################"
echo "NOTE: Do NOT proceed with pipeline until you've really looked at the clustering and figured out which clusters go together"
echo "######################"
