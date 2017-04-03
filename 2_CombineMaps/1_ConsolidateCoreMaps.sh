#! /bin/bash

#Compare the linkage group assignments from the different maps and consolidate

MERGEMAP="/home/jgw87/Software/linkage_mapping/MergeMap/consensus_map.exe"
basedir="../MakeNewMaps"
scaffold_key="/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/0_pseudochromosome_scaffold_key.txt"

##Copy core maps over
#mapTifton=$basedir/TiftonData/3i_som_rippled_combined.hmp.txt
#mapPatancheru=$basedir/Patancherux863/3i_Patancheru_rippled_combined.hmp.txt
#mapSadore=$basedir/Sadore/3i_boubacar_rippled_combined.hmp.txt
#cp $mapTifton ./1_som_core_map.hmp.txt
#cp $mapPatancheru ./1_Patancheru_core_map.hmp.txt
#cp $mapSadore ./1_boubacar_core_map.hmp.txt

##Convert to scaffold assignments
#for pop in som Patancheru boubacar; do
#	#Set the correct population directory to get the final map
#	map=""
#	case $pop in
#	som)
#		map=$mapTifton ;;
#	Patancheru)
#		map=$mapPatancheru ;;
#	boubacar)
#		map=$mapSadore ;;
#	*)
#		echo "WARNING! Unknown population!" ;;
#	esac
#	
#	#Run analyses
#	python3 $basedir/6a_MapScaffoldsToLinkageGroups.py -i 1_${pop}_core_map.hmp.txt -c "chr" -l $scaffold_key -o 1a_${pop}_scaffold_assignments.txt
#	python3 $basedir/6d_ConsolidateScaffoldsIntoBins.py -s 1a_${pop}_scaffold_assignments.txt -m $map -o 1b_${pop}_scaffold_bins.txt
#done

###Make diagrams showing correspondences between maps
alpha=0.1
mysom=1b_som_scaffold_bins.txt
myPatancheru=1b_Patancheru_scaffold_bins.txt
myboubacar=1b_boubacar_scaffold_bins.txt
#Rscript 1d_CompareMappings.r 1c_map_comparison_som_Patancheru.png som,Patancheru $alpha $mysom $myPatancheru
#Rscript 1d_CompareMappings.r 1c_map_comparison_som_boubacar.png som,boubacar $alpha $mysom $myboubacar
#Rscript 1d_CompareMappings.r 1c_map_comparison_Patancheru_boubacar.png Patancheru,boubacar $alpha $myPatancheru $myboubacar


#####Make a manual file showing the equivalencies among the linkage groups
linker="1c_linkage_group_equivalencies_MANUAL.txt"
#python3 2_ConsolidateGeneticMapsByLinkageGroup.py -i $mysom,$myPatancheru,$myboubacar -l $linker -o 2_my_consolidated_linkage_groups_core.txt

###Filter input maps based on consensus
#python3 2a_FilterMergeMapInputByConsensus.py -i $mysom -l $linker -c 2_my_consolidated_linkage_groups_core.txt -o 2a_som_scaffolds_filtered.txt
#python3 2a_FilterMergeMapInputByConsensus.py -i $myPatancheru -l $linker -c 2_my_consolidated_linkage_groups_core.txt -o 2a_Patancheru_scaffolds_filtered.txt
#python3 2a_FilterMergeMapInputByConsensus.py -i $myboubacar -l $linker -c 2_my_consolidated_linkage_groups_core.txt -o 2a_boubacar_scaffolds_filtered.txt


##Format filtered scaffold assignments for mergemap; if pass option -g (group), only outputs markers on the specified group
##Also make configuration file for mergemap for each linkage group
#newlink="3_linkage_group_equivalencies_revised_names.txt"
#echo -e "lg\t2a_som_scaffolds_filtered.txt\t2a_Patancheru_scaffolds_filtered.txt\t2a_boubacar_scaffolds_filtered.txt" > $newlink
#tail -n +2 $linker >> $newlink
#for lg in $(seq 1 7); do
#	dir=2b_merge_chr${lg}
#	config=$dir/config.txt
#	if [ ! -d $dir ]; then mkdir $dir; fi
#	if [ -e $config ]; then rm $config; fi	#Clear out previous config files
#	for pop in som Patancheru boubacar; do
#		python3 3_FormatGeneticMapsForMergeMap.py -i 2a_${pop}_scaffolds_filtered.txt -l $newlink -o $dir/3_${pop}_mergemap_filtered_lg${lg}.txt -g $lg -b 100000
#		weight=1
#		if [ $pop == "som" ]; then weight=4; fi	#Weigh som's data better than the others
#		echo "$pop $weight 3_${pop}_mergemap_filtered_lg${lg}.txt" >> $config
#	done
#done

#Run Mergemap on each linkage group individually
for lg in $(seq 1 7); do
	cd 2b_merge_chr${lg}
	$MERGEMAP config.txt &
	cd ../
done
wait
echo "Finished!"

