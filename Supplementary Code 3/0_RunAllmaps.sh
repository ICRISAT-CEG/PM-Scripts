#! /bin/bash

# Run ALLMAPS to join the different Pearl Millet maps together

TASSEL5="perl /home/jgwall/Software/TASSEL/tassel-5-standalone/run_pipeline.pl -Xms10g -Xmx40g"

workdir=Allmaps
# workdir=Allmaps_bak
if [ ! -e $workdir ]; then mkdir $workdir; fi

maxprocs=7
foxtail_genome=S_italica_genome
pm_version1_genome=/home/jgwall/0_PostdocFiles/Working_Files/GBS/Genomes/pennisetum_glaucum_reference_v1.0/pm_assembly_v1.0_20140816.mainchroms.fa.gz


# # copy synteny over
# cut -f1-4 -d, ../DevosLabOrdering/1e_map_with_gene_positions.txt > $workdir/1e_synteny_map.allmaps.txt

##############
# Formatting
##############

# # Copy original files over
# src=RemakeGbsMaps
# cp $src/Som/3i_som_rippled_combined.hmp.txt $workdir/1_som_map_core.hmp.txt
# cat $src/Som/4k_som_anchored_map_reordered_combined.hmp.txt | sed -r "s|^target_||" > $workdir/1_som_map_extended.hmp.txt   # strip target_ prefix
# cp $src/Som/5_som_anchored_tags_best.txt $workdir/1_som_tags.txt
# 
# cp $src/841/3i_841_rippled_combined.hmp.txt $workdir/1_841_map_core.hmp.txt
# cat $src/841/4k_841_anchored_map_reordered_combined.hmp.txt | sed -r "s|^target_||" > $workdir/1_841_map_extended.hmp.txt   # strip target_ prefix
# cp $src/841/5_841_anchored_tags_best.txt $workdir/1_841_tags.txt

# # Get scaffold for snps
# for infile in $workdir/1_*.hmp.txt; do
# # # #   # Map SNPs back to scaffolds
# # # #   outfile=${infile/1_/1a_}
# # # #   outfile=${outfile/.txt/.scaffolds.txt}
# # # #   echo python3 1_MapSnpsBackToScaffolds.py -i $infile -l 0_pseudochromosome_scaffold_key.txt -o $outfile -c "chr" # DEPRECATED; new map pipeline did this part already
#   
#   # Put in format for ALLMAPS
#   mapfile=${infile/1_/1b_}
#   mapfile=${mapfile/.txt/.scaffold_map.txt}
#   echo -e "Scaffold ID,scaffold position,LG,genetic position" > $mapfile
# #   cut -f1,3,4 $outfile | tail -n+2 | sed -r -e "s|_|,|g" -e "s|\t|,|g" >> $mapfile
#   cut -f1,3,4 $infile | tail -n+2 | sed -r -e "s|_|,|g" -e "s|\t|,|g" >> $mapfile
# done

# # Map tags to scaffolds
# samfile=../GBSv2/1b_tags.sam
# samtools view -S $samfile | cut -f1,3,4 | sed -r -e "s|^tagSeq=||" -e "s|chr(.+)\t(.+)|S\1_\2|" -e "s|\*\t0|S0_0|"  > $workdir/0_tag_scaffold_key.txt # '-99' is for 
# for tags in $workdir/1_*tags.txt; do
#   # Rename so can use same SNP-scaffolding script
#   renamed=${tags/.txt/.renamed.txt}
#   python3 1_RenameTags.py -i $tags -o $renamed -k $workdir/0_tag_scaffold_key.txt
#   
#   # Add scaffold information
#   outfile=${tags/1_/1a_}
#   outfile=${outfile/.txt/.scaffolds.txt}
#   python3 1_MapSnpsBackToScaffolds.py -i $renamed -l 0_pseudochromosome_scaffold_key.txt -o $outfile -c "chr"
#   
#   # Convert into ALLMAPS format
#   mapfile=${tags/1_/1b_}
#   mapfile=${mapfile/.txt/.scaffold_map.txt}
#   echo -e "Scaffold ID,scaffold position,LG,genetic position" > $mapfile
#   cut -f1,5,6 $outfile | tail -n+2 | sed -r -e "s|_|,|g" -e "s|\t|,|g" >> $mapfile
# done
# 
# # Renumber linkage groups to match consensus
# ref_key=0_reference_v1_scaffold_key.txt
# for map in $workdir/1b_*scaffold_map.txt; do
#   # Make equivalency key
#   outfile=${map/1b_/1c_}
#   outfile=${outfile/.scaffold_map.txt/.allmaps.txt}
#   table=${outfile/.allmaps.txt/.lg_equivalencies.txt}
#   majority=""
#   if [[ $map == *"tags"* ]] ; then majority="--majority-rule"; fi # Handle tags differently because have small amount of noise making their scaffolds map uniquely poorly
#   python3 1c_MatchLinkageGroupsToReference.py -i $map -r $ref_key -t  $table -o $outfile $majority
# done
# 
# # Confirm That everything's ordered the same
# python3 1d_ConfirmMatches.py -i $workdir/1c*.allmaps.txt -o $workdir/1d_confirm_lg_matches.any.png
# python3 1d_ConfirmMatches.py -i $workdir/1c*.allmaps.txt -o $workdir/1d_confirm_lg_matches.majority_rule.png --majority-rule
# 
# # Reduce scaffolds to only be on the most-supported linkage group
# for infile in $workdir/1c*.allmaps.txt; do
#   outfile=${infile/1c_/1e_}
#   python3 1e_RemoveScaffoldsFromMinorLGs.py -i $infile -a $workdir/1c*.allmaps.txt -o $outfile
# #   break
# done
# 
# # Split by chromosome to improve computation
# infiles=$workdir/1e*.allmaps.txt
# for infile in $infiles; do
#   outstem=${infile/1e_/1f_}
#   
#   chroms=$(tail -q -n +2 $infiles | cut -f3 -d, | sort --unique)
#   for c in $chroms; do
#     outfile=${outstem/.allmaps.txt/.chr$c.allmaps.txt}
#     Rscript -e "data=read.csv('$infile'); colnames(data)=scan('$infile', what=character(), nlines=1, sep=','); data=subset(data, data\$LG==$c); write.csv(data, file='$outfile', quote=F, row.names=F)"
#   done
# done


########################
# Run ALLMAPS on the assembled maps; previously installed the 'jcvi' module for python with pip2
###################


fasta=PM.scaffold.v1.1.fa
synteny=../DevosLabOrdering/1e_map_with_gene_positions.txt

# # Combine maps by chromosome
# for chrom in $(seq 1 7); do  # Reverse order to complement what doing on server
#   weights=$workdir/2_map_weights.chr$chrom.txt
#   infiles=$workdir/1f*.chr$chrom.allmaps.txt
#   combined=$workdir/2_combined_maps.chr$chrom.bed
#   
# #   python -m jcvi.assembly.allmaps merge -w $weights -o $combined $infiles
# #   
# #   # # Overwrite default weight file
# #   rm $weights
# #   for map in $infiles; do
# #     name=${map%%.*}
# #     weight=1
# #     if [[ $name == *"synteny"* ]]; then weight=10; fi
# #     if [[ $name == *"map_core"* ]]; then weight=5; fi
# #     if [[ $name == *"map_extended"* ]]; then weight=3; fi
# #     if [[ $name == *"tags"* ]]; then weight=1; fi
# #     echo -e "$name $weight" >> $weights
# #   done
# #   
# #   # Run ALLMAPS
# #   echo "Running ALLMAPS on chromosome $chrom"
# #   time python -m jcvi.assembly.allmaps path -w $weights --cpus $maxprocs --seed 1 --noplot $combined $fasta | tee $combined.log
# 
# #   Check things by synteny
#   python3 2b_GetGenePositionsFromAllmaps.py -i $workdir/2_combined_maps.chr$chrom.chr.agp --syntenymap $synteny -o $workdir/2a_gene_positions.chr$chrom.txt 
# done
# gzip -f $workdir/2_*.fasta $workdir/2_*.tour

# foxtail=../DevosLabOrdering/2_s_italica_gene_synteny.txt
# old_pm=../DevosLabOrdering/2_existing_pm_genes.txt
# combined_new=$workdir/3_assembly_gene_positions.txt
# head -q -n1 $workdir/2a_gene_positions.chr*.txt | head -n1 > $combined_new
# tail -q -n+2 $workdir/2a_gene_positions.chr*.txt >> $combined_new
# python3 3_PlotSynteny.py -a $combined_new -b $foxtail -o $workdir/3_synteny_new_and_foxtail.png
# python3 3_PlotSynteny.py -a $combined_new -b $old_pm -o $workdir/3_synteny_new_and_old.png
# python3 3_PlotSynteny.py -a $old_pm -b $foxtail -o $workdir/3_synteny_old_and_foxtail.png #--debug

# # Check LD
# for pop in som 841; do
#   rev ../MakeNewMapsByScaffold_MyCalls/$pop/7_${pop}_map_remade.hmp.txt| cut -f2- | rev > $workdir/3a_${pop}_scaffolds.hmp.txt # Removing an extra, non-genotype column at the end
#   python3 3a_ReorderScaffoldHapmap.py -i $workdir/3a_${pop}_scaffolds.hmp.txt -a $workdir/2_combined_maps.chr*.chr.agp -o $workdir/3a_${pop}_scaffolds.neworder.hmp.txt
#   $TASSEL5 -SortGenotypeFilePlugin -inputFile $workdir/3a_${pop}_scaffolds.neworder.hmp.txt -outputFile $workdir/3a_${pop}_scaffolds.neworder.sorted.hmp.txt
#   $TASSEL5 -h $workdir/3a_${pop}_scaffolds.neworder.sorted.hmp.txt -ld -ldType All -ldHetTreatment Homozygous -ldd png -o $workdir/3a_${pop}_scaffolds.neworder.sorted.ld.png -ldplotsize 2000
# done



################
# Some scaffolds got duplicated; I think I forgot to remove the scaffolds on the less-supported linakge groups from the synteny map (which may be the best one anyway), so time to do it after the fact
################

#python3 3b_TallyDuplicatedScaffolds.py -i $workdir/1e*.allmaps.txt -o $workdir/3b_duplciated_scaffolds.txt

# # # Based on manual inspection of the synteny stuff, I'm favoring putting all the duplicate ones where the genetic maps say. In some cases it's because the map used for synteny has no genotype calls
# # #    in others it's b/c my two independent mapping pops put then elsewhere
# scaffold1400    LG4
# scaffold1764    LG1
# scaffold4916    LG3
# scaffold8671    LG2
# scaffold9526    LG1

# manual_set="scaffold1400:chr4 scaffold1764:chr1 scaffold4916:chr3 scaffold8671:chr2 scaffold9526:chr1"
# for infile in $workdir/2_combined_maps.chr*.chr.agp; do
#   outfile=${infile/2_/3c_}
#   outfile=${outfile/.chr.agp/.scaffold_key.agp}
#   python3 3c_RemoveScaffoldsFromAgp.py -i $infile -o $outfile -m $manual_set
# done

# # Check LD when those scaffolds are removed (to see if it's any better)
# for pop in som 841; do
#   python3 3a_ReorderScaffoldHapmap.py -i $workdir/3a_${pop}_scaffolds.hmp.txt -a $workdir/3c_combined_maps.chr*.scaffold_key.agp -o $workdir/3d_${pop}_scaffolds.neworder.nodups.hmp.txt
#   $TASSEL5 -SortGenotypeFilePlugin -inputFile $workdir/3d_${pop}_scaffolds.neworder.nodups.hmp.txt -outputFile $workdir/3d_${pop}_scaffolds.neworder.nodups.sorted.hmp.txt
#   $TASSEL5 -h $workdir/3d_${pop}_scaffolds.neworder.nodups.sorted.hmp.txt -ld -ldType All -ldHetTreatment Homozygous -ldd png -o $workdir/3d_${pop}_scaffolds.neworder.nodups.sorted.ld.png -ldplotsize 2000
# done
# head -q -n1 $workdir/3c_*.scaffold_key.agp | head -n 1 > $workdir/3e_combined_assembly_key.agp
# tail -q -n +2 $workdir/3c_*.scaffold_key.agp >> $workdir/3e_combined_assembly_key.agp


################
# Orientation relative to Rajaram et al
################

# # Based off what I did for Som Punnuri to map his specific linkage map up
# BOWTIE="/home/jgwall/Software/Aligners/bowtie2-2.2.6/"
# INDEX=PM_contigs_unassembled/PM.scaffold.v0.4
# $BOWTIE/bowtie2 --very-sensitive-local -a -f --fr -p $maxprocs -x $INDEX -1 0_Rajaram_primers_fwd.fasta -2 0_Rajaram_primers_rev.fasta -S $workdir/4_Rajaram_primer_alignments.sam | tee $workdir/4_Rajaram_primer_alignments.all.log
# samtools view -hf 0x2 -S $workdir/4_Rajaram_primer_alignments.sam >  $workdir/4a_Rajaram_primer_alignments.filtered.sam
# python3 4b_GetGoodAmpliconScaffolds.py -i $workdir/4a_Rajaram_primer_alignments.filtered.sam -o $workdir/4b_good_amplicon_matches.txt --multialign
# Rscript 4c_AddMapLocations.r -i $workdir/4b_good_amplicon_matches.txt -a $workdir/3e_combined_assembly_key.agp -o $workdir/4c_good_amplicon_matches.all.map_locations.txt
# python3 4c_AddRajaramAndPlot.py -i $workdir/4c_good_amplicon_matches.all.map_locations.txt -r 0_Rajaram_map.txt -o $workdir/4d_matched_amplicon_locations.all.txt \
#   -g $workdir/4d_matched_amplicon_locations.all.png --rescale-chroms 
# python3 4c_AddRajaramAndPlot.py -i $workdir/4c_good_amplicon_matches.all.map_locations.txt -r 0_Rajaram_map.txt \
#   -g $workdir/4d_matched_amplicon_locations.all.flipped.png --rescale-chroms --flip-chr 1 2 5 6 7

# # Hard to tell about chromosome 3, but chroms 1, 2, 5, 6, and 7 all need to be flipped. Do so, then check synteny again
# foxtail=../DevosLabOrdering/2_s_italica_gene_synteny.txt
# synteny=../DevosLabOrdering/1e_map_with_gene_positions.txt
# Rscript 4e_FlipChroms.r -i $workdir/3e_combined_assembly_key.agp -o $workdir/4e_combined_assembly_key.flipped.agp -f chr1 chr2 chr5 chr6 chr7
# python3 2b_GetGenePositionsFromAllmaps.py -i $workdir/4e_combined_assembly_key.flipped.agp --syntenymap $synteny -o $workdir/4e_combined_assembly_key.flipped.gene_positions.txt
# python3 3_PlotSynteny.py -a $workdir/4e_combined_assembly_key.flipped.gene_positions.txt -b $foxtail -o $workdir/4e_combined_assembly_key.flipped.synteny.png

# Make a new FASTA assembly
fasta=PM.scaffold.v1.1.fa
python -m jcvi.formats.agp build $workdir/4e_combined_assembly_key.flipped.agp $fasta $workdir/4f_assembly1.1.fa

# perl 4f_agp2fasta.pl $workdir/4e_combined_assembly_key.flipped.agp $fasta $workdir/4f_assembly1.1.fa

samtools faidx $workdir/4f_assembly1.1.fa
