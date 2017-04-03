#! /bin/bash


#Global variables
TASSEL4="/home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx12g"
BOWDIR=/home/jgw87/Software/Aligners/bowtie2-2.0.2
GENOME=/media/STORAGE/Working_Files/GBS/Genomes/pennisetum_glaucum_v0.1/penn_glauc_v0.1	#Reference genome

keyfile=0_gbs_keyfile.txt	#Keyfile of samples
flowcells=/media/STORAGE/Working_Files/GBS/Flowcells	#Directory where fastq files are kept

##Step 1: Make directories (numbering is according to the step they are output from)
if [ ! -e 01_RawSequence ]; then ln -s $flowcells 01_RawSequence; fi	#Link to flowcells if haven't already. 
mkdir 02_TagCounts
mkdir 02_TagCounts/02a_IndividualCounts
mkdir 02_TagCounts/02b_MergedCounts
mkdir 02_TagCounts/02c_TagsToFastq
mkdir 03_SAM
mkdir 04_TOPM
mkdir 05_TBT
mkdir 05_TBT/05a_IndividualTBT
mkdir 05_TBT/05b_MergedTBT
mkdir 06_HapMap
mkdir 06_HapMap/06a_UnfilteredSNPs
mkdir 06_HapMap/06b_MergeDupSNPs
mkdir 06_HapMap/06c_HapMapFilteredSNPs

#Step 2: Do tag counts
$TASSEL4 -fork1 -FastqToTagCountPlugin -i 01_RawSequence -k $keyfile -e ApeKI -o 02_TagCounts/02a_IndividualCounts -endPlugin -runfork1
$TASSEL4 -fork1 -MergeMultipleTagCountPlugin -i 02_TagCounts/02a_IndividualCounts -o 02_TagCounts/02b_MergedCounts/refseq02_merged_tagcounts.cnt -c 5 -endPlugin -runfork1
$TASSEL4 -fork1 -TagCountToFastqPlugin -i 02_TagCounts/02b_MergedCounts/refseq03_merged_tagcounts.cnt -o 02_TagCounts/02c_TagsToFastq/refseq03_merged_tagcounts.fastq -c 5 -endPlugin -runfork1

##Step 3: Align tags to the genome scaffolds
perl $BOWDIR/bowtie2 -q -p 8 --very-sensitive-local -x $GENOME -U 02_TagCounts/02c_TagsToFastq/refseq03_merged_tagcounts.fastq -S 03_SAM/refseq11_tags.sam

#Step 4: Create TOPM (Tags On Physical Map) file
$TASSEL4 -fork1 -SAMConverterPlugin -i 03_SAM/refseq11_tags.sam -o 04_TOPM/refseq11_tags.topm.bin -endPlugin -runfork1

#Step 5: Make TBT (Tags By Taxa) files
$TASSEL4 -fork1 -FastqToTBTPlugin -i 01_RawSequence -k $keyfile -e ApeKI -o 05_TBT/05a_IndividualTBT -t 02_TagCounts/02b_MergedCounts/refseq02_merged_tagcounts.cnt -y -endPlugin -runfork1
$TASSEL4 -fork1 -MergeTagsByTaxaFilesPlugin -i 05_TBT/05a_IndividualTBT -o 05_TBT/05b_MergedTBT/refseq02.tbt.byte -x -endPlugin -runfork1
#Convert to text for later use
$TASSEL4 -fork1 -BinaryToTextPlugin -i 05_TBT/05b_MergedTBT/refseq02.tbt.byte -o 05_TBT/05b_MergedTBT/refseq02.tbt.txt -t TBTByte -endPlugin -runfork1


##Step 6: Make Hapmap and VCF Files
##Call SNPs
$TASSEL4 -fork1 -DiscoverySNPCallerPlugin -i 05_TBT/05b_MergedTBT/refseq03.tbt.byte -y -o 06_HapMap/06a_UnfilteredSNPs/refseq11_c+.hmp.txt -vcf -m 04_TOPM/refseq11_tags.topm.bin -sC 1 -eC 20 -mnMAF 0.01 -mnMAC 10 -mnLCov 0.1 -endPlugin -runfork1
$TASSEL4 -fork1 -MergeDuplicateSNPsPlugin -hmp 06_HapMap/06a_UnfilteredSNPs/refseq11_c+.hmp.txt -o 06_HapMap/06b_MergeDupSNPs/refseq11_c+.hmp.txt -sC 1 -eC 20 -callHets -misMat 0.1 -endPlugin -runfork1
$TASSEL4 -fork1 -MergeDuplicateSNPsPlugin -vcf 06_HapMap/06a_UnfilteredSNPs/refseq11_c+.vcf -o 06_HapMap/06b_MergeDupSNPs/refseq11_c+.vcf -sC 1 -eC 20 -callHets -misMat 0.1 -endPlugin -runfork1


###Step 7: Combine hapmap files into one
head -n 1 06_HapMap/06b_MergeDupSNPs/refseq11_c1.hmp.txt > 06_HapMap/refseq11_unfiltered.hmp.txt
tail -q -n +2 06_HapMap/06b_MergeDupSNPs/refseq11_c*.hmp.txt >> 06_HapMap/refseq11_unfiltered.hmp.txt
##Combine vcf files into one
grep --no-filename "^#" 06_HapMap/06b_MergeDupSNPs/refseq11_c1.vcf > 06_HapMap/refseq11_unfiltered.vcf	#Take all commented lines (begin with #)
grep --no-filename -v "^#" 06_HapMap/06b_MergeDupSNPs/refseq11_c*.vcf >> 06_HapMap/refseq11_unfiltered.vcf	#Take all non-commented lines


