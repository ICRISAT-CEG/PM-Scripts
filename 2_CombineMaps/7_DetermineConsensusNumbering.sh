#! /bin/bash

#Align the primers from the consensus genetic map from Rajaram et al 2013 to the genome to identify linkage groups

BOWDIR=/home/jgw87/Software/Aligners/bowtie2-2.0.2
GENOME=/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/penn_glauc_v1.1
SCAFFKEY=/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/0_pseudochromosome_scaffold_key.txt
sourcedir=0_RajaramConsensusMapContigs

###Align to genome
#contigs=7_rajaram_contigs.fasta
#cat $sourcedir/ICMP_full_Length_17sequences_1.txt $sourcedir/IPES_full_Length_99sequences_1.txt > $contigs
#max_alignments=5
#perl $BOWDIR/bowtie2 -q -p 8 --very-sensitive-local -k $max_alignments -x $GENOME -U $contigs -f -S 7a_rajaram_contig_alignments.sam 	# Can do paired, but at this point benefit seems minimal


###Connect to scaffolds and my linkage groups
#python3 7b_DetermineScaffoldsOfContigs.py -s 7a_rajaram_contig_alignments.sam -l $SCAFFKEY -o 7b_rajaram_contig_scaffold_assignments.txt
#python3 7c_AnchorPrimersToLinkageGroups.py -i 7b_rajaram_contig_scaffold_assignments.txt -l 6e_expanded_map2_som_841.txt \
#                 -o 7c_contig_lg_assignments_map2_som_841.txt -c $sourcedir/0_Rajaman_consensus_map_ssr_assignments.txt

##Make heatmap to manually curate
#python3 7d_MakeHeatmapOfLgAssignments.py -i 7c_contig_lg_assignments_map2_som_841.txt -o 7d_map2_som_841_matrix

#Renumber map and try again
python3 7e_RenumberLinkageGroups.py -i 6e_expanded_map2_som_841.txt -l 7e_Rajaram_LG_equivalencies_MANUAL.txt -o 7e_expanded_map2_som_841_renumbered.txt
python3 7c_AnchorPrimersToLinkageGroups.py -i 7b_rajaram_contig_scaffold_assignments.txt -l 7e_expanded_map2_som_841_renumbered.txt \
                 -o 7f_contig_lg_assignments_map2_som_841_renumbered.txt -c $sourcedir/0_Rajaman_consensus_map_ssr_assignments.txt
python3 7d_MakeHeatmapOfLgAssignments.py -i 7f_contig_lg_assignments_map2_som_841_renumbered.txt -o 7g_map2_som_841_renumbered_matrix
