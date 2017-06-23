#! /bin/bash

#Subscript to impute and order a single bootstrapped map

MSTMAP="/home/jgw87/Software/linkage_mapping/MSTmap/MSTMap.exe"
TASSEL="/home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx8g"

bootfile=$1
parent1=$2
parent2=$3
nomap_dist=$4
nomap_size=$5

imputed=${bootfile/.abh.txt/.abh.imputed.txt}
echo "Imputing $bootfile to $imputed"
Rscript 3c_ImputeAndConsolidateMarkers.r $bootfile $imputed 1 "-" FALSE

##Convert to MSTmap format
pval=2	#Set to >1 so doesn't try to make linkage groups, just takes the ones I have
mst_in=${imputed/.abh.imputed.txt/.mst_in.txt}
perl 3d_FormatHapmapForMSTmap.pl $imputed $mst_in $parent1 $parent2 DH $pval $nomap_dist $nomap_size

#Run MSTmap
mst_out=${mst_in/mst_in/mst_out}
$MSTMAP $mst_in $mst_out  #Can put in background because MSTmap runs pretty quickly
