#! /bin/bash

#Refine the pearl millet scaffold order using BGI's SNPs (which have WGS data on them, so are more complete than just GBS)

TASSEL4="perl /home/jgw87/Software/tassel4-standalone/run_pipeline.pl -Xms4g -Xmx20g"
TASSEL5="perl /home/jgw87/Software/tassel5-standalone/run_pipeline.pl -Xms4g -Xmx20g"


##Convert BGI's matrix to a hapmap for easier handling
#Rscript 1_ConvertBgiCallsToHapmap.r GBS_Genotype.matrix sample.name.order.list.txt 1_bgi_pmigap.hmp.txt.gz

##Order for TASSEL and then filter
#minCount=20	#Need at least 20 present calls to be able to call LD
#$TASSEL5 -SortGenotypeFilePlugin -inputFile 1_bgi_pmigap.hmp.txt.gz -outputFile 1a_bgi_pmigap_sorted.hmp.txt.gz
#$TASSEL5 -h 1a_bgi_pmigap_sorted.hmp.txt.gz -filterAlign -filterAlignMinCount $minCount -export 1b_bgi_pmigap_filt.hmp.txt.gz
##Remove sites with a lot of heterozygosity
#$TASSEL5 -h 1b_bgi_pmigap_filt.hmp.txt.gz -genotypeSummary site -export 1b_bgi_pmigap_filt.sitesummary.txt
#Rscript 0_FindHighlyHetSites.r 1b_bgi_pmigap_filt.sitesummary.txt 1b_bgi_pmigap_filt.badsites.txt 0.1
#$TASSEL5 -h 1b_bgi_pmigap_filt.hmp.txt.gz -excludeSiteNamesInFile 1b_bgi_pmigap_filt.badsites.txt -export 1c_bgi_pmigap_filt_goodsites.hmp.txt.gz

##Do a quick check to see if LD is actually present
#mkdir 0_ldcheck
#$TASSEL5 -h 1c_bgi_pmigap_filt_goodsites.hmp.txt.gz -subsetSites 4000 -export 0_ldcheck/1c_bgi_pmigap_subset.hmp.txt.gz
#Rscript 1c_QuickHackToMakeLdWork.r 0_ldcheck/1c_bgi_pmigap_subset.hmp.txt.gz 0_ldcheck/1c_bgi_pmigap_subset.hacked.hmp.txt
#$TASSEL5 -h 0_ldcheck/1c_bgi_pmigap_subset.hacked.hmp.txt -ld -ldType All -ldHetTreatment Homozygous -ldd png -ldplotsize 2000 -ldplotlabels false -o 0_ldcheck/1c_bgi_pmigap_subset.ld.png
#Rscript 1c_QuickHackToMakeLdWork.r 1c_bgi_pmigap_filt_goodsites.hmp.txt.gz 0_ldcheck/1c_bgi_pmigap_filt.hacked.hmp.txt
####$TASSEL5 -h 0_ldcheck/1c_bgi_pmigap_filt.hacked.hmp.txt -ld -ldType SlidingWindow -ldHetTreatment Homozygous -ldd png -ldplotsize 2000 -ldplotlabels false -o 0_ldcheck/1c_bgi_pmigap_filt.hacked.ld.png # -Doesn't work; too much memory to write PNG


##Get just the scaffold tips, testing different sizes
CLASSPATH="-cp /home/jgw87/Software/IntelliJ_idea-IC-135.690/Projects/JasonSandbox/out/artifacts/JasonSandbox_jar/JasonSandbox.jar CalcLdAmongScaffoldTips"
for dist in 100 1000 10000; do
	#Rscript 1d_GetOutermostSnps.r 1c_bgi_pmigap_filt_goodsites.hmp.txt.gz $dist 1d_outer$dist
	#tail -n +2 1d_outer${dist}_tips.txt | cut -f1 > 1d_outer${dist}_snpnames.txt
	#$TASSEL5 -h 1b_bgi_pmigap_filt.hmp.txt.gz -includeSiteNamesInFile 1d_outer${dist}_snpnames.txt -export 1e_bgi_pmigap_outer${dist}.hmp.txt.gz
	java $CLASSPATH 1e_bgi_pmigap_outer${dist}.hmp.txt.gz 1d_outer${dist}_tips.txt 1f_bgi_pmigap_outer${dist}.ld.txt.gz 1f_bgi_pmigap_outer${dist}.mean_ld.txt.gz
	#break
done
