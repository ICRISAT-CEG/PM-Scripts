#! /bin/bash

#Reorder the Pearl Millet linkage map using the Traveling Salesman problem as a framework

map=/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet/20140527_AlignToRefseqV1.1/ConsolidateMaps/7e_expanded_map2_som_841_renumbered.txt
concorde="/media/STORAGE/Software/linkage_mapping/concorde/bin"

for bin in 100 1000 10000; do
#	Rscript 2_FormatResultsForTravSalesman.r 1f_bgi_pmigap_outer${bin}.mean_ld.txt.gz 2_outer${bin}_formatted_for_tsp 2>&1 | tee -a $log
#	Rscript 2a_ReorderMapWithTravelingSalesman.r $map 2_outer${bin}_formatted_for_tsp.txt 2a_map_reordered_from_outer${bin}_tsp.txt $concorde 2>&1 | tee -a $log
#	python3 2b_OutputReorderedMap.py -i 2a_map_reordered_from_outer${bin}_tsp.txt -m $map -o 2b_new_map_from_outer${bin}.txt
	Rscript 2c_CalculateTotalMapLength.r 2b_new_map_from_outer${bin}.txt 2_outer${bin}_formatted_for_tsp.txt 2c_new_map_from_outer${bin}_lengths.txt
done

