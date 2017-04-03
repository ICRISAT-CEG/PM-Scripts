#! /bin/bash

#Take bootstrap samples of RILs and do MSTmap on them and see if that improves the order any

MSTMAP="/home/jgw87/Software/linkage_mapping/MSTmap/MSTMap.exe"
TASSEL="/home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx8g"

pop=$1
popdir=$2
best_cluster=$3
groups=$4
parent1=$5
parent2=$6
nbootstraps=$7
nomap_dist=$8
nomap_size=$9


##Separate into individual hapmaps
Rscript 3b_SeparateClusteredHapmaps.r $popdir/3a_${pop}_clustered_ward_k${best_cluster}.hmp.txt $popdir/3b_${pop}_clustered_ward_k${best_cluster} $groups
combo=$popdir/3b_${pop}_clustered_ward_k${best_cluster}_combined.hmp.txt
head -n 1 $popdir/3b_${pop}_clustered_ward_k${best_cluster}_LG1.hmp.txt > $combo
tail -q -n +2 $popdir/3b_${pop}_clustered_ward_k${best_cluster}_LG*.hmp.txt >> $combo
$TASSEL -fork1 -h $combo -ld -ldType All -ldHetTreatment Homozygous -ldd png -o ${combo/.hmp.txt/.ld.png} -ldplotsize 2000 -runfork1

#Set up TASSEL commands to unify resulting hapmaps
forks=""
inputs=""
runs=""

##Do bootstrapping
bootdir=$popdir/3b_bootstraps
if [ ! -d $bootdir ]; then mkdir $bootdir; fi
iter=1
for chr_map in $popdir/3b_${pop}_clustered_ward_k*_LG*.hmp.txt; do
	basename=${chr_map/$popdir\//}
	basename=${basename/.hmp.txt/}
	
	##Convert to ABH format
	abh=$bootdir/$basename.abh.txt
	Rscript 3c_ConvertHapmapToABH.r $chr_map $abh $parent1 $parent2 
	##Impute and consolidate to minimal marker set
	abh_imp=$bootdir/${basename}_trimmed.abh.txt
	Rscript 3c_ImputeAndConsolidateMarkers.r $abh $abh_imp 1 "-" TRUE 
	##Bootstrap
	Rscript 3c_BootstrapHapmap.r $abh_imp ${abh_imp/.abh.txt/_boot+.abh.txt} $nbootstraps $parent1 $parent2
	
	##Setup and run MSTmap on all bootstraps
	pval=2	#Set to >1 so doesn't try to make linkage groups, just takes the ones I have
	for infile in $bootdir/${basename}_trimmed_boot*.abh.txt; do
		mst_in=${infile/.abh.txt/.mst_in.txt}
		perl 3d_FormatHapmapForMSTmap.pl $infile $mst_in $parent1 $parent2 DH $pval $nomap_dist $nomap_size
	
		#Run MSTmap
		mst_out=${mst_in/mst_in/mst_out}
		$MSTMAP $mst_in $mst_out & #Can put in background because MSTmap runs pretty quickly
	done
	wait
	
	##Compile results into one file
	outpos=${chr_map/.hmp.txt/_bootstrapped_positions.txt}
	outpos=${outpos/3b_/3e_}
	perl 3e_CompileBootstrappedMaps.pl $bootdir/${basename}_trimmed_boot*.mst_out.txt $outpos
	
	#Reconvert to hapmap (both abh and hmp format)
	outhap=$popdir/${basename/3b_/3f_}_reordered.hmp.txt
	outabh=${outhap/.hmp.txt/.abh.txt}
	mychr=${chr_map##*_LG}
	mychr=${mychr/.hmp.txt/}
	Rscript 3f_ParseCompiledMap.r $outpos $abh_imp $outabh $mychr
	Rscript 3f_ParseCompiledMap.r $outpos $chr_map $outhap $mychr
	$TASSEL -fork1 -h $outhap -ld -ldType All -ldHetTreatment Homozygous -ldd png -o ${outhap/.hmp.txt/.ld.png} -ldplotsize 2000 -runfork1
	
	#Compile commands for later unification with TASSEL
	forks="$forks -fork${iter} -h $outhap"
	inputs="$inputs -input${iter}"
	runs="$runs -runfork${iter}"
	iter=$(($iter + 1))

done

#Recombine hapmaps into one and do LD on combined results
combinedhap=$popdir/3f_${pop}_clustered_ward_k${best_cluster}_combined.hmp.txt
$TASSEL $forks -combine${iter} $inputs -union -export $combinedhap $runs -runfork${iter}
$TASSEL -fork1 -h $combinedhap -ld -ldType All -ldHetTreatment Homozygous -ldd png -o ${combinedhap/.hmp.txt/.ld.png} -ldplotsize 2000 -runfork1
