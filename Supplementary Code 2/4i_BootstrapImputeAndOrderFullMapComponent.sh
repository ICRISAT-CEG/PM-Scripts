#! /bin/bash

#Take bootstrap samples of RILs and do MSTmap on them and see if that improves the order any

MSTMAP="/home/jgw87/Software/linkage_mapping/MSTmap/MSTMap.exe"
TASSEL="/home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx8g"
TASSEL_LOCAL="./4i_TASSEL_LOCAL/run_pipelinev4_local.pl -Xms4g -Xmx8g"

pop=$1
popdir=$2
parent1=$3
parent2=$4
startchrom=$5
endchrom=$6
nbootstraps=$7
nomap_dist=$8
nomap_size=$9
maxprocs=${10}

#Set up TASSEL commands to unify resulting hapmaps
forks=""
inputs=""
runs=""

##Do bootstrapping
hmp=$popdir/4h_${pop}_anchored_map_goodsites_nohets_filt.hmp.txt
bootdir=$popdir/4i_bootstraps
if [ ! -d $bootdir ]; then mkdir $bootdir; fi
iter=1
for chr in $(seq $startchrom $endchrom); do
	#Extract just this chromosome
	basename=4i_${pop}_anchored_map_LG${chr}
	$TASSEL_LOCAL -fork1 -h $hmp -filterAlign -filterAlignLocus $chr -export $popdir/$basename.hmp.txt -runfork1
	
	
	####Convert to ABH format
	abh=$bootdir/$basename.abh.txt
	Rscript 3c_ConvertHapmapToABH.r $popdir/$basename.hmp.txt $abh $parent1 $parent2 
	#Bootstrap. In this case bootstrap before imputing because want to destabilize the clustering used in imputing (since is greedy, different orders will impute different things together)
	Rscript 3c_BootstrapHapmap.r $abh ${abh/.abh.txt/_boot+.abh.txt} $nbootstraps $parent1 $parent2
	#Impute but do not consolidate into minimal marker set
	for boot in $bootdir/${basename}_boot*.abh.txt; do
		./4i2_ImputeAndOrderFullMapSubcomponent.sh $boot $parent1 $parent2 $nomap_dist $nomap_size &
		while [ $(pgrep -c 4i2_Imp) -ge $maxprocs ]; do	#Do a pgrep to count how many subprocesses have spawned
			sleep 1s
		done
	done
	wait
	
	##Compile results into one file
	outpos=$popdir/4j_${pop}_anchored_map_LG${chr}_bootstrapped_positions.txt
	perl 3e_CompileBootstrappedMaps.pl $bootdir/4i_${pop}_anchored_map_LG${chr}_boot*.mst_out.txt $outpos
	
	#Reconvert to hapmap (both abh and hmp format)
	outhap=$popdir/4j_${pop}_anchored_map_reordered_LG${chr}.hmp.txt
	outabh=${outhap/.hmp.txt/.abh.txt}
	Rscript 3f_ParseCompiledMap.r $outpos $abh $outabh $chr
	Rscript 3f_ParseCompiledMap.r $outpos $popdir/$basename.hmp.txt $outhap $chr
	$TASSEL -fork1 -h $outhap -ld -ldType All -ldHetTreatment Homozygous -ldd png -o ${outhap/.hmp.txt/.ld.png} -ldplotsize 2000 -runfork1
	
	#Compile commands for later unification with TASSEL
	forks="$forks -fork${iter} -h $outhap"
	inputs="$inputs -input${iter}"
	runs="$runs -runfork${iter}"
	iter=$(($iter + 1))
	
done

##Recombine (non-imputed) hapmaps into one and do LD on combined results
combinedhap=$popdir/4k_${pop}_anchored_map_reordered_combined.hmp.txt
$TASSEL $forks -combine${iter} $inputs -union -export $combinedhap $runs -runfork${iter}
##Take 5000 random lines to do LD on
subset=$popdir/4k_${pop}_anchored_map_reordered_combined_subset.hmp.txt
perl 5a_RandomlySubsetFile.pl $combinedhap 5000 header $subset
$TASSEL -fork1 -h $subset -ld -ldType All -ldHetTreatment Homozygous -ldd png -o ${subset/.hmp.txt/.ld.png} -ldplotsize 2000 -runfork1
