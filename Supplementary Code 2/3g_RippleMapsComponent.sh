#! /bin/bash

#Ripple hapmaps in R/qtl

TASSEL="perl /home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms2g -Xmx8g"

pop=$1
popdir=$2
winsize=$3
maxprocs=$4

##Do rippling in R/qtl
outfiles=""
for hap_in in $popdir/3f_${pop}_clustered_ward_k*_LG*_reordered.abh.txt; do
	
	#Convert to R/qtl format
	rqtl_in=${hap_in/3f_/3g_}
	rqtl_in=${rqtl_in/.abh.txt/.rqtl.csvr.txt}
	perl 3g_ConvertHapmapAbhToRqtl.pl $hap_in $rqtl_in
	rqtl_out=${rqtl_in/3g_/3h_}
	
	#Run rippling
	Rscript 3h_RippleBootstrapped.r $rqtl_in $rqtl_out $winsize &
	outfiles="$outfiles $rqtl_out.csv"
	
	#Sleep while waiting for processes to finish
	while [ $(pgrep -xc R) -ge $maxprocs ]; do
		sleep 5s
	done
	
done
wait

#Set up TASSEL commands to unify resulting hapmaps
forks=""
inputs=""
runs=""
iter=1
#Put back into hapmap
for outfile in $outfiles; do
	hapmap=${outfile/.rqtl.csvr.txt.csv/.hmp.txt}
	outpng=${hapmap/.hmp.txt/.ld.png}
	perl 3i_ConvertRqtlToHapmap.pl $outfile $hapmap
	$TASSEL -fork1 -h $hapmap -ld -ldType All -ldHetTreatment Homozygous -ldd png -o $outpng -ldplotsize 2000 -runfork1
	
	#Compile commands for later unification with TASSEL
	forks="$forks -fork${iter} -h $hapmap"
	inputs="$inputs -input${iter}"
	runs="$runs -runfork${iter}"
	iter=$(($iter + 1))
done

#Final TASSEL commands on combined data
combinedhap=$popdir/3i_${pop}_rippled_combined.hmp.txt
$TASSEL $forks -combine${iter} $inputs -union -export $combinedhap $runs -runfork${iter}
$TASSEL -fork1 -h $combinedhap -ld -ldType All -ldHetTreatment Homozygous -ldd png -o ${combinedhap/.hmp.txt/.ld.png} -ldplotsize 2000 -runfork1
