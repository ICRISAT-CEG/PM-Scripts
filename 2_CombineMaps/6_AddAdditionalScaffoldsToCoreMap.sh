#! /bin/bash

#Take the combined core map and add additional scaffolds by determining for each one which core scaffold they are in the most LD with
TASSEL4="perl /home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx20g"
TASSEL5="perl /home/jgw87/Software/tassel5-standalone/run_pipeline.pl -Xms4g -Xmx12g"
TASSEL5_LOCAL="perl /home/jgw87/Software/NetBeansProjects/run_pipelinev5_local.pl -Xms2g -Xmx10g"

hapmap="../06_HapMap/refseq11_unfiltered.hmp.txt"

###Filter down to just segregating sites in the appropriate populations
#$TASSEL4 -fork1 -h $hapmap \
#	-fork2 -includeTaxaInFile 4d_tifton_taxa.txt -input1 -filterAlign -filterAlignMinFreq 0.2 -filterAlignMinCount 40 -export 6_segregating_tifton.hmp.txt \
#	-fork3 -includeTaxaInFile 4d_patancheru_taxa.txt -input1 -filterAlign -filterAlignMinFreq 0.2 -filterAlignMinCount 40 -export 6_segregating_patancheru.hmp.txt \
#	-fork4 -includeTaxaInFile 4d_sadore_taxa.txt -input1 -filterAlign -filterAlignMinFreq 0.2 -filterAlignMinCount 80 -export 6_segregating_sadore.hmp.txt \
#	-runfork1 -runfork2 -runfork3 -runfork4

#minsites=3	#Minimum number of segregating sites a scaffold needs to have in a population to be considered
#for pop in patancheru tifton sadore; do
#	Rscript 6a_SortHapmap.r 6_segregating_${pop}.hmp.txt 6a_segregating_${pop}_sorted.hmp.txt
#	java -cp ./6b_GetMeanLdOfScaffolds/JasonScriptsTassel5.jar Misc.MakeDistanceMatrixOfScaffolds 6a_segregating_${pop}_sorted.hmp.txt 4_snp_scaffolds_all.txt 6b_segregating_${pop}_scaffold_distances.txt $minsites
#done

#infiles="6b_segregating_patancheru_scaffold_distances.txt,6b_segregating_sadore_scaffold_distances.txt,6b_segregating_tifton_scaffold_distances.txt"
#map=5e_core_map_scaffolds.txt
#cutoff=0.4	#Minimum LD between scaffolds to be counted
#python3 6c_FindBestAnchorForUnmappedScaffolds.py -i $infiles -m $map -o 6c_expanded_map_scaffolds.txt -c $cutoff

#cp 6c_expanded_map_scaffolds.txt 6d_expanded_map_scaffolds_all.txt
#tail -n +2 5e_core_map_scaffolds.txt >> 6d_expanded_map_scaffolds_all.txt
#lengths=/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/0_scaffold_lengths.txt
#python3 ../MakeNewMaps/6e_CountLengthOfAssembledScaffolds.py -i 6d_expanded_map_scaffolds_all.txt -a 6d_expanded_map_scaffolds_all.txt -l $lengths -o 6d_expanded_map_all_lengths.txt

#Iteratively add the extended maps to the core, working from best to worst populations
#python3 6e_AnchorExtendedMapToExistingMap.py -c 6d_expanded_map_scaffolds_all.txt -s 6d_expanded_map_scaffolds_all.txt -e ../MakeNewMaps/Tifton/6d_tifton_scaffold_bins.txt -o 6e_expanded_map1_tifton.txt
#python3 6e_AnchorExtendedMapToExistingMap.py -c 6d_expanded_map_scaffolds_all.txt -s 6e_expanded_map1_tifton.txt -e ../MakeNewMaps/Patancheru/6d_patancheru_scaffold_bins.txt -o 6e_expanded_map2_tifton_patancheru.txt
#python3 6e_AnchorExtendedMapToExistingMap.py -c 6d_expanded_map_scaffolds_all.txt -s 6e_expanded_map2_tifton_patancheru.txt -e ../MakeNewMaps/Sadore/6d_sadore_scaffold_bins.txt -o 6e_expanded_map3_tifton_patancheru_sadore.txt

##Calculate lengths
lengths=/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/0_scaffold_lengths.txt
for map in 6e*_map*.txt; do
	python3 ../MakeNewMaps/6e_CountLengthOfAssembledScaffolds.py -i $map -a $map -l $lengths -o ${map/_map/_map_lengths}
done
