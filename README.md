0_RunGbsPipeline.sh                                //GBS Tassel Commands
1_MakeIndividualMaps                               //Construction of Individual maps
├── 0_MakeMapOfPatancheruLines.sh                  //SNPs into a genetic map
├── 0_MakeMapOfSadoreLines.sh
├── 0_MakeMapOfTiftonLines.sh
├── 1a_ExtractVcfSiteInfo.pl                       //Get stats on the VCF SNPs
├── 1a_GraphVcfQualityDistributions.r              //Graph distribution of SNP calls
├── 1b_ApplyAdHocVcfFilters.pl                     //Apply ad-hoc filters to VCF file
├── 1_FilterVcfFilesComponent.sh                   //Filter VCF values
├── 2a_StripLocusInfoFromHapmap.r                  //Strip the scaffold and position numbers from a hapmap
├── 2b_FilterOutHighlyHetSitesComponent.sh         //Identify highly het sites and filter out
├── 2b_FindHighlyHetSites.r                        //Find highly heterozygous sites and output them
├── 2c_IdentifyAndRemoveOutcrossesComponent.sh     //Check the number of rare alleles in each line
├── 2c_TestForOutcrossedTaxa_HetsAndMissing.r      //Find highly heterozygous sites and output them
├── 2d_TestForOutcrossedTaxa_RareAlleles.r         //Calculate and show # of rare alleles by taxon
├── 2_FilterOnMafAndGetLdComponent.sh              //Filter map files on MAF, coverage, and heterozygosit
├── 3b_BootstrapMapsComponent.sh                   //Take bootstrap samples of RILs and do MSTmap on them
├── 3b_SeparateClusteredHapmaps.r                  //Separate a hapmap into specified linkage groupings
├── 3c_BootstrapHapmap.r                           //Bootstrap a hapmap for MSTmap ordering
├── 3c_ConvertHapmapToABH.r                        //Convert hapmap to Rqtl
├── 3c_ImputeAndConsolidateMarkers.r               //Identify redundant markers in a hapmap, use them to impute each other
├── 3d_FormatHapmapForMSTmap.pl                    //Reformat TASSEL's hapmap output for processing with MSTmap
├── 3e_CompileBootstrappedMaps.pl                  //Compile results of bootstrapped maps
├── 3f_ParseCompiledMap.r                          //Parse the compiled, bootstrapped map
├── 3g_ConvertHapmapAbhToRqtl.pl                   //Convert an MSTmap input and output file pair into an Rqtl file
├── 3g_RippleMapsComponent.sh                      //Ripple hapmaps in R/qtl
├── 3_GroupMarkersByClusteringComponent.sh         //Convert final hapmap to numeric and then run clustering
├── 3_GroupMarkersByClustering.r                   //Take a numeric output from TASSEL and cluster for LD analysis
├── 3h_RippleBootstrapped.r                        //Analyze map by R/qtl
├── 3i_ConvertRqtlToHapmap.pl                      //Convert an Rqtl "csvr" (rotated CSV) file back into a Hapmap
├── 4_AnchorSnpsOnCoreMapComponent.sh              //Anchor SNPs based their LD with the core markers
├── 4b_JavaScriptComponents_sourcecode
│   └── FindRsqAlongChrom.java
├── 4b_ProcessRsqAlongChrom.r                      //Manual check the assignment along chromosome
├── 4d_ReorderTargetsIntoMap.pl                    //Take the Target SNPs and reorder them back into the reference hapmap
├── 4i2_ImputeAndOrderFullMapSubcomponent.sh       //Subscript to impute and order a single bootstrapped map
├── 4i_BootstrapImputeAndOrderFullMapComponent.sh  //Take bootstrap samples of RILs and do MSTmap on them
├── 5a_GraphSubsetOfAllTags.r                      //Graph a subset of mapped tags
├── 5_AnchorTagsToMap.r                            //Map sequencing tags to a genetic map
├── 5a_RandomlySubsetFile.pl                       //Select a random subset of snps to generate an MDS plot
├── 5_MapSequencingTagsComponent.sh                //Make tags-by-taxa file of raw data, convert to text, and then map in R
├── 6a_MapScaffoldsToLinkageGroups.py   
├── 6_AssignScaffoldsToLinkageGroupsComponent.sh   //Take final orders and create documents
├── 6b_GetPositionOfTags.py                        //Take the mapped tags and determine their position in the original fasta files
├── 6d_ConsolidateScaffoldsIntoBins.py             //Take the snp and tag scaffold assignments and consolidate into a unified framework
└── 6e_CountLengthOfAssembledScaffolds.py          //Take inputs of binned scaffolds and calculate the total lengths
2_CombineMaps
├── 1_ConsolidateCoreMaps.sh                    //Compare the linkage group assignments from the different maps and consolidate
├── 1d_CompareMappings.r                        //Compare the arrangements of linkage maps with scaffolds
├── 2a_FilterMergeMapInputByConsensus.py        //Take scaffold assignments and filter them by the consensus
├── 2_ConsolidateGeneticMapsByLinkageGroup.py   //Take genetic maps generated and hook them together, calling consensus locations for scaffolds
├── 3_FormatGeneticMapsForMergeMap.py           //Take a genetic map and reformat it for MergeMap
├── 4a_GetListOfSnpsInScaffolds.py              //consolidated map and list of snp-scaffold assignments and get out the list of SNPs
├── 4c_ReorderHapmapFromConsolidatedMap.py      //hapmap and the list of SNP-scaffold assignments
├── 4_FinishConsolidatingMaps.sh                //Finish consolidating the combined maps
├── 4_MatchSnpsToScaffolds.py                   //determine which scaffold each SNP belongs on
├── 6_AddAdditionalScaffoldsToCoreMap.sh        //Take the combined core map and add additional scaffolds
├── 6a_SortHapmap.r                             //Sort a hapmap based on chromosome, position, and then name
├── 6c_FindBestAnchorForUnmappedScaffolds.py    //distance matrices for each scaffold combination and determine the best anchor scaffold 
├── 6e_AnchorExtendedMapToExistingMap.py        //Take an existing scaffold map and add all the extended scaffolds
├── 7b_DetermineScaffoldsOfContigs.py           //Take a SAM alignment of primer pairs (with _fwd and _rev in names) and determine which scaffold the pairing is on.
├── 7c_AnchorPrimersToLinkageGroups.py          //scaffold assignments from 1b_ and put with linkage groups
├── 7_DetermineConsensusNumbering.sh            //align the primers from the consensus genetic map to the genome to identify linkage groups
├── 7d_MakeHeatmapOfLgAssignments.py            //Take the list of LG assignments and turn it into a heatmap
├── 7e_RenumberLinkageGroups.py                 //Take a genetic map and renumber its linkage groups as needed
└── java_sourcecode
    └── MakeDistanceMatrixOfScaffolds.java
3_OrderScaffolds
├── 0_FindHighlyHetSites.r                        //Find highly heterozygous sites and output
├── 1_ConvertBgiCallsToHapmap.r                   //Convert SNP calls to hapmap format
├── 1c_QuickHackToMakeLdWork.r                    //hapmap file to make it so TASSEL will do LD
├── 1d_GetOutermostSnps.r                         //Take a hapmap and extract the list of SNPs in the outermost X of each chromosome/scaffolds
├── 1_RefineScaffoldOrder.sh                      //Refine the pearl millet scaffold order using SNPs
├── 2_ApplyTravelingSalesmanProblem.r             //Use an R package for solving the Traveling Salesman Problem to try to find optimal grouping and ordering of scaffolds
├── 2a_ReorderMapWithTravelingSalesman.r          //Reorder the linkage map, first ordering by Traveling Salesman within each linkage bin and then ordering across linkage bins
├── 2b_OutputReorderedMap.py                      //Take output from step 2a_ and the original bin maps and create a new map
├── 2c_CalculateTotalMapLength.r                  //Calculate total map length and output to screen and to a file
├── 2_FormatResultsForTravSalesman.r              //Format my scaffold LD measures for use in the traveling salesman problem
├── 2_ReorderMapWithTravelingSalesman.sh          //Reorder the Pearl Millet linkage map using the Traveling Salesman problem as a framework
└── java_sourcecode
    └── CalcLdAmongScaffoldTips.java

3 directories, 74 files
