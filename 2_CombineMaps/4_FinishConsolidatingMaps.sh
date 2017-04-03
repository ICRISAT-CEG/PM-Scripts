#! /bin/bash

#Finish consolidating the combined maps
TASSEL4="perl /home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx12g"
TASSEL5="perl /home/jgw87/Software/tassel5-standalone/run_pipeline.pl -Xms4g -Xmx12g"
TASSEL5_LOCAL="perl /home/jgw87/Software/NetBeansProjects/run_pipelinev5_local.pl -Xms2g -Xmx10g"


linker=/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/0_pseudochromosome_scaffold_key.txt
hapmap="../06_HapMap/refseq11_unfiltered.hmp.txt"

##Consolidate the final, linear core maps into a single file
#combined=4_core_map_combined.txt
#echo -e "lg\tscaffolds\tpos" > $combined
#for chr in $(seq 1 7); do
#	grep -e "^scaffold" -e "^C" 2b_merge_chr${chr}/linear_map_chart.txt | sed -r -e "s|^|$chr\t|"  >> $combined
#done

##Subset out the SNPs that fall on these scaffolds and put them in the new order
###Get list of all SNPs called and their source scaffold
#python3 4_MatchSnpsToScaffolds.py -i $hapmap -c chr -l $linker -o 4_snp_scaffolds_all.txt
##Get the SNPs that fall in the anchor scaffolds and order them according to the new map
#python3 4a_GetListOfSnpsInScaffolds.py -m 4_core_map_combined.txt -s 4_snp_scaffolds_all.txt -o 4a_consolidated_map_snps.txt
#tail -n +2 4a_consolidated_map_snps.txt | cut --fields=1 > 4a_snp_names.txt
#$TASSEL4 -fork1 -h $hapmap -includeSiteNamesInFile 4a_snp_names.txt -export 4b_core_scaffold_snps.hmp.txt -runfork1
#python3 4c_ReorderHapmapFromConsolidatedMap.py -i 4b_core_scaffold_snps.hmp.txt -s 4a_consolidated_map_snps.txt -o 4c_consolidated_map_reordered.hmp.txt -m 100

##Extract individual maps and check LD across them
basedir="../MakeNewMaps"
sort $basedir/SomData/0_som_taxa.txt | uniq > 4d_som_taxa.txt
sort $basedir/841x863/0_841_taxa.txt | uniq > 4d_841_taxa.txt
sort $basedir/Boubacar/0_boubacar_taxa.txt | uniq > 4d_boubacar_taxa.txt
$TASSEL5 -fork1 -h 4c_consolidated_map_reordered.hmp.txt \
	-fork2 -includeTaxaInFile 4d_som_taxa.txt -input1 -filterAlign -filterAlignMinFreq 0.2 -filterAlignMinCount 50 -export 4d_new_map_som.hmp.txt.gz \
	-fork3 -includeTaxaInFile 4d_841_taxa.txt -input1 -filterAlign -filterAlignMinFreq 0.2 -filterAlignMinCount 80 -export 4d_new_map_841.hmp.txt.gz \
	-fork4 -includeTaxaInFile 4d_boubacar_taxa.txt -input1 -filterAlign -filterAlignMinFreq 0.2 -filterAlignMinCount 100 -export 4d_new_map_boubacar.hmp.txt.gz \
	-runfork1 -runfork2 -runfork3 -runfork4

for pop in som 841 boubacar; do
	$TASSEL5 -fork1 -h 4d_new_map_${pop}.hmp.txt.gz -subsetSites 5000 -export 4e_new_map_${pop}_subset.hmp.txt.gz -runfork1
	$TASSEL5 -fork1 -h 4e_new_map_${pop}_subset.hmp.txt.gz -ld -ldType All -ldHetTreatment Homozygous -ldd png -ldplotsize 2000 -o 4e_new_map_${pop}_subset.ld.png -runfork1
done

##Determine amount of sequence in consolidated map
#lengths=/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/0_scaffold_lengths.txt
#cut --fields=2-4 4a_consolidated_map_snps.txt | uniq > 5e_core_map_scaffolds.txt
#python3 ../MakeNewMaps/6e_CountLengthOfAssembledScaffolds.py -i 5e_core_map_scaffolds.txt -a 5e_core_map_scaffolds.txt -l $lengths -o 5e_core_map_lengths_ignore_ambiguous.txt
