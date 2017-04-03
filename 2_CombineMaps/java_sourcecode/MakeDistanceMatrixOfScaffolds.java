/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.logging.Level;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import static net.maizegenetics.analysis.popgen.LinkageDisequilibrium.getLDForSitePair;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.pipeline.TasselPipeline;
import net.maizegenetics.stats.statistics.FisherExact;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Utils;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.log4j.Logger;

/**
 *
 * @author jgw87
 */
public class MakeDistanceMatrixOfScaffolds {

    private static final Logger myLogger = Logger.getLogger(MakeDistanceMatrixOfScaffolds.class);
    private static final int snpID = 0, scaffoldID = 1;   //Columns with needed information in scaffold file
    private static int minsites;  //Minimum number of sits for a scaffold to be tested
    private static final int myMinTaxaForEstimate = 20;

    public static void main(String[] args) {
        /*//Initial stuff for debugging
        String workdir = "/media/STORAGE/Working_Files/GBS/Analysis/PearlMillet/20140324_AlignToRefseqV0.3/ConsolidateGeneticMaps/";
        String infile = workdir + "6a_segregating_841_sorted.hmp.txt";
        String scaffoldFile = workdir + "5_snp_scaffolds_all.txt";
        String outfile = workdir + "6b_segregating_841_scaffold_distances.txt";*/

        //Arguments; going positionally b/c lazy
        String infile = args[0];
        String scaffoldFile = args[1];
        String outfile = args[2];
        minsites = Integer.parseInt(args[3]);
        GenotypeTable hmp = ImportUtils.readFromHapmap(infile);
        HashMap<String, Integer> sitekey = getSitesOfSnpNames(hmp);
        HashMap<String, ArrayList> scaffolds = readInScaffolds(scaffoldFile);
        HashMap<String, Double> distances = calculateDistances(hmp, scaffolds, sitekey);
        outputDistanceMatrix(scaffolds, distances, outfile);
    }

    private static HashMap<String, ArrayList> readInScaffolds(String scaffoldFile) {
        System.out.println("Load SNP-scaffold assignments");
        HashMap<String, ArrayList> scaffolds = new HashMap();
        int n = 0;

        //Read in file
        try {
            BufferedReader reader = Utils.getBufferedReader(scaffoldFile);
            reader.readLine();   //Clear header
            String line = reader.readLine();
            //Go through each one and add
            while (line != null) {
                String[] data = line.split("\t");
                String snp = data[snpID];
                String scaff = data[scaffoldID];

                //Craete new array list if needed
                if (!scaffolds.containsKey(scaff)) {
                    scaffolds.put(scaff, new ArrayList<String>());
                }
                //Add SNP to ArrayList and move on
                scaffolds.get(scaff).add(snp);
                n++;
                line = reader.readLine();
            }

        } catch (IOException ex) {
            java.util.logging.Logger.getLogger(MakeDistanceMatrixOfScaffolds.class.getName()).log(Level.SEVERE, null, ex);
        }

        System.out.println("\tTotal " + n + " snps loaded in " + scaffolds.size() + " scaffolds");

        return scaffolds;
    }

    private static HashMap<String, Integer> getSitesOfSnpNames(GenotypeTable hmp) {
        HashMap<String, Integer> sitekey = new HashMap();
        for (int i = 0; i < hmp.numberOfSites(); i++) {
            sitekey.put(hmp.siteName(i), i);
        }
        return sitekey;
    }

    private static HashMap<String, Double> calculateDistances(GenotypeTable hmp, HashMap<String, ArrayList> scaffolds, HashMap<String, Integer> sitekey) {
        System.out.println("Calculating pairwise mean LD for each scaffold pair");
        String[] myscaffolds = scaffolds.keySet().toArray(new String[0]);
        HashMap<String, Double> distances = new HashMap(myscaffolds.length ^ 2 / 2);  //Appx how many comparisons there are max
        GenotypeTable hmp_nohets = GenotypeTableBuilder.getHomozygousInstance(hmp);

        //Loop over every combination of scaffolds and get all pairwise LD between them, then record as a double
        int n = 0;
        for (int i = 0; i < myscaffolds.length; i++) {
            String s0 = myscaffolds[i];
            Integer[] s0_sites = findSitesForScaffold(sitekey, scaffolds.get(s0));
            if (s0_sites.length < minsites) {
                //System.out.println("Skipping " + s0);
                continue;
            }
            for (int j = i + 1; j < myscaffolds.length; j++) {
                String s1 = myscaffolds[j];
                Integer[] s1_sites = findSitesForScaffold(sitekey, scaffolds.get(s1));
                if (s1_sites.length < minsites) {
                    //System.out.println("Skipping " + s1);
                    continue;
                }

                n++;
                if (n % 10000 == 0) {
                    System.out.println("\tProcessed " + n + " scaffold combinations");
                }
                /*if (n > 1000) { //For debugging
                 break;
                 }*/

                String key = s0 + "|" + s1;
                //System.out.println("\tCalculating for " + key);
                double mean_ld = calcMeanLdForScaffolds(hmp_nohets, s0_sites, s1_sites, myMinTaxaForEstimate);
                distances.put(key, mean_ld);
                //System.out.println("\t\tMean = " +  mean_ld);
            }
        }
        return distances;
    }

    //Get site indices, taking into account that not all loaded sites will be in each hapmap
    private static Integer[] findSitesForScaffold(HashMap<String, Integer> sitekey, ArrayList<String> snps) {
        ArrayList<Integer> indices = new ArrayList();
        for (int i = 0; i < snps.size(); i++) {
            String mysnp = snps.get(i);
            //System.out.println(mysnp);
            if (sitekey.containsKey(mysnp)) {
                indices.add(sitekey.get(mysnp));
                //System.out.println("\tis at index " + indices.get(indices.size()-1));
            }
        }
        Integer[] indices2 = indices.toArray(new Integer[0]);
        return indices2;

    }

    private static double calcMeanLdForScaffolds(GenotypeTable hmp, Integer[] s0_sites, Integer[] s1_sites, int minTaxa) {
        DescriptiveStatistics values = new DescriptiveStatistics();

        /*for (int site0 : s0_sites) {
            for (int site1 : s1_sites) {
                BitSet rMj = hmp.allelePresenceForAllTaxa(site0, GenotypeTable.WHICH_ALLELE.Major);
                BitSet rMn = hmp.allelePresenceForAllTaxa(site0, GenotypeTable.WHICH_ALLELE.Minor);
                BitSet cMj = hmp.allelePresenceForAllTaxa(site1, GenotypeTable.WHICH_ALLELE.Major);
                BitSet cMn = hmp.allelePresenceForAllTaxa(site1, GenotypeTable.WHICH_ALLELE.Minor);
                double r2 = getRsqForSitePair(rMj, rMn, cMj, cMn, 2, minTaxa, -1.0f);
                //Only take real number values
                //System.out.println("\t\tR2 for " + site0 + " by " + site1 + " is " + r2);
                if (Double.isNaN(r2)) {  //Skip NaN estimates
                    continue;
                }
                values.addValue(r2);
                //System.out.println("\t\t\tTotal: " + total + " among " + n);
            }
        }*/
        return values.getPercentile(50);    // 50th percentile = median; why they don't have a wrapper, I don't know
    }

    private static void outputDistanceMatrix(HashMap<String, ArrayList> scaffolds, HashMap<String, Double> distances, String outfile) {
        String[] scaffoldNames = scaffolds.keySet().toArray(new String[0]);
        try {
            BufferedWriter writer = Utils.getBufferedWriter(outfile);
            for (String s0 : scaffoldNames) {
                writer.append(s0);
                for (String s1 : scaffoldNames) {
                    //If is the same scaffold, LD is perfect so output 1.0 and skip to next
                    if (s1 == s0) {
                        writer.append("\t1.0");
                        continue;
                    }

                    //Make key and check that it (or reverse) is contained. If not, set to NA
                    String key = s0 + "|" + s1;
                    if (!distances.containsKey(key)) {    //If key isn't there, try reverse
                        key = s1 + "|" + s0;
                        if (!distances.containsKey(key)) {    //If still doesn't contain, write as missing
                            writer.append("\tNA");
                            continue;
                        }
                    }

                    //Output distance
                    writer.append("\t" + distances.get(key));
                }
                writer.append("\n");
            }
            writer.close();
        } catch (IOException ex) {
            java.util.logging.Logger.getLogger(MakeDistanceMatrixOfScaffolds.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    //Taken from linkage disequilibrium plugin and modified for what I need
    public static double getRsqForSitePair(BitSet rMj, BitSet rMn, BitSet cMj, BitSet cMn,
            int minMinorCnt, int minCnt, float minR2) {
        // float[] results = {Float.NaN, Float.NaN, Float.NaN, Float.NaN};
        double results = Double.NaN;
        int n = 0;
        int[][] contig = new int[2][2];
        n += contig[1][1] = (int) OpenBitSet.intersectionCount(rMn, cMn);
        n += contig[1][0] = (int) OpenBitSet.intersectionCount(rMn, cMj);
        if (contig[1][0] + contig[1][1] < minMinorCnt) {
            return results;
        }
        n += contig[0][1] = (int) OpenBitSet.intersectionCount(rMj, cMn);
        if (contig[0][1] + contig[1][1] < minMinorCnt) {
            return results;
        }
        n += contig[0][0] = (int) OpenBitSet.intersectionCount(rMj, cMj);
        //results.n = n;
        if (n < minCnt) {
            return results;
        }
        double rValue = LinkageDisequilibrium.calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minCnt);
        return rValue;
        /*
         results.r2 = (float) rValue;
         if (Double.isNaN(rValue)) {
         return results;
         }
         results.dprime = (float) LinkageDisequilibrium.calculateDPrime(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minCnt);;
         if (rValue < minR2) {
         return results;
         }
         double pValue = myFisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
         results.p = (float) pValue;
         return results;*/
    }

}

