/**
 * Created by jgw87 on 6/17/14.
 */

import net.maizegenetics.analysis.popgen.LDResult;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.stats.statistics.FisherExact;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.Utils;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

public class CalcLdAmongScaffoldTips {

    //private static final Logger myLogger = Logger.getLogger(CalcLdAmongScaffoldTips.class);
    private static final int nameID=0, scaffoldID=1, sideID=2;  //Columns in tip file

    public static void main(String[] args) {
        //args = new String[] {"/home/jgw87/Working_Files/GBS/Genomes/pennisetum_glaucum_v1.1/BGI_PMiGAP_calls/1e_bgi_pmigap_outer100.hmp.txt.gz", "/home/jgw87/1f_bgi_pmigap_outer100.ld2.txt"};

        //Deal with arguments
        printUsage();
        String infile = args[0];
        String tipfile = args[1];
        String outfile_raw = args[2];
        String outfile_mean = args[3];

        //Calculate raw Rsquared
        GenotypeTable genos = ImportUtils.readFromHapmap(infile);
        //genos = FilterGenotypeTable.getInstance(genos, 0, 1000);    //For debugging - limit to first 1000 sites
        String[] sitenames = getSiteNames(genos);
        float[][] rsq = calcRsq(genos);
        //outputRsq(outfile_raw, sitenames, rsq);

        //Calculate mean LD between each group
        HashMap<String, ArrayList> groups = loadGroups(tipfile);
        String[] groupnames = groups.keySet().toArray(new String[0]);
        float[][] mean_rsq = calcMeanRsq(groups, groupnames, rsq, sitenames);
        outputRsq(outfile_mean, groupnames, mean_rsq);
    }

    private static void printUsage() {
        System.out.println("Usage of this class with positional arguments:\n" +
                "\t(0) Name of input hapmap\n" +
                "\t(1) Output file of LD matrix measurements\n");
    }

    private static String[] getSiteNames(GenotypeTable genos) {
        String[] names = new String[genos.numberOfSites()];
        for (int i = 0; i < genos.numberOfSites(); i++) {
            names[i] = genos.siteName(i);
        }
        return names;
    }

    //Much of this taken from TASSEL's calculateBitLDForHaplotype class
    private static float[][] calcRsq(GenotypeTable genos) {
        //Set up homozygous-only sites and various support variables
        GenotypeTable workingAlignment = GenotypeTableBuilder.getHomozygousInstance(genos);
        int nsites = workingAlignment.numberOfSites();
        float[][] rsq = new float[nsites][nsites];
        int myMinTaxaForEstimate = 20;
        FisherExact myFisherExact = new FisherExact((2 * workingAlignment.numberOfTaxa()) + 10);

        //Actually perform LD calculation
        long total = nsites * nsites / 2;
        long n = 0;
        long marker = Math.round((double) total / 20);
        System.out.println("Performing " + total + " calculations");
        for (int row = 0; row < nsites; row++) {
            for (int col = row; col < nsites; col++) {
                n++;
                if (n % marker == 0) {
                    float percent = (float) (n * 100 / total);
                    System.out.println("\tCompleted " + Math.round(percent) + "%");
                }
                //If on the diagonal, set to 1 and skip to next
                if (row == col) {
                    rsq[row][col] = 1;
                    continue;
                }

                //Actual LD calculations
                BitSet rMj = workingAlignment.allelePresenceForAllTaxa(row, WHICH_ALLELE.Major);
                BitSet rMn = workingAlignment.allelePresenceForAllTaxa(row, WHICH_ALLELE.Minor);
                BitSet cMj = workingAlignment.allelePresenceForAllTaxa(col, WHICH_ALLELE.Major);
                BitSet cMn = workingAlignment.allelePresenceForAllTaxa(col, WHICH_ALLELE.Minor);
                LDResult ldr = LinkageDisequilibrium.getLDForSitePair(rMj, rMn, cMj, cMn, 2, myMinTaxaForEstimate, -1.0f, myFisherExact, row, col);
                rsq[row][col] = ldr.r2();
                rsq[col][row] = ldr.r2();

            }
        }

        return rsq;
    }

    private static void outputRsq(String outfile, String[] labels, float[][] rsq) {
        System.out.println("Writing Rsquared to " + outfile);
        try {
            GZIPOutputStream zipstream = new GZIPOutputStream(new FileOutputStream(new File(outfile)));
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(zipstream, "UTF-8"));
            DecimalFormat formatter = new DecimalFormat("0.0000");
            for (int row = 0; row < rsq.length; row++) {
                writer.append(labels[row]);
                for (int col = 0; col < rsq[0].length; col++) {
                    writer.append("\t" + formatter.format(rsq[row][col]));
                }
                writer.append("\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static HashMap<String, ArrayList> loadGroups(String tipfile) {
        BufferedReader reader = Utils.getBufferedReader(tipfile);
        HashMap <String, ArrayList> groups = new HashMap<>();
        int total=0;
        try {
            reader.readLine();  //Clear header
            String line = reader.readLine();
            while (line != null) {
                total++;
                String[] tokens = line.split("\t");
                String snp = tokens[nameID];
                String tip = tokens[scaffoldID] + "|" + tokens[sideID];

                if(!groups.containsKey(tip)){
                    groups.put(tip, new ArrayList());
                }
                groups.get(tip).add(snp);
                line = reader.readLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("\tLoaded " + groups.size() + " groups containing " + total + " sites\n");
        return groups;
    }

    private static float[][] calcMeanRsq(HashMap<String,ArrayList> groups, String[] groupnames, float[][] rsq, String[] sitenames){
        System.out.println("Calculating mean Rsq among groups");

        //Make a key for what name is at which index
        HashMap<String, Integer> sitekey = new HashMap<>();
        for(int i=0; i<sitenames.length; i++){
            //System.out.println("Putting " + sitenames[i] + " at key " + i);
            sitekey.put(sitenames[i], i);
        }

        //Calculate Rsquared
        float[][] means = new float[groupnames.length][groupnames.length];
        for(int a=0; a<groupnames.length;a++){
            ArrayList groupA = groups.get(groupnames[a]);
            for(int b=a; b<groupnames.length; b++){
                ArrayList groupB = groups.get(groupnames[b]);
                /*if(!(groupnames[a].equals("14|left") && groupnames[b].equals("14|left") )){ //For debugging
                    continue;
                }*/
                //System.out.println("Calculating mean Rsq for " + groupnames[a] + " and " + groupnames[b]);
                float meanrsq = getMeanRsqForGroups(groupA, groupB, rsq, sitekey);
                means[a][b] = meanrsq;
                means[b][a] = meanrsq;
            }
        }

        return means;
    }

    private static float getMeanRsqForGroups(ArrayList<String> groupA, ArrayList<String> groupB, float[][] rsq, HashMap<String, Integer> sitekey){
        float total=0;
        int n=0;
        for(String a: groupA){
            if(!sitekey.containsKey(a)){
                continue;
            }
            //System.out.println("GroupA = "+a);
            int indexA = sitekey.get(a);
            for(String b: groupB){
                //System.out.println("\tGroupB = "+b);
                if(!sitekey.containsKey(b)){
                    continue;
                }
                if(a.equals(b)){    //Skip comparing a SNP to itself
                    //System.out.println("\t\t\tSkipping because " +  a + " is the same site as " + b);
                    continue;
                }
                int indexB = sitekey.get(b);
                float myval = rsq[indexA][indexB];;
                if(! Float.isNaN(myval)) {
                    n++;
                    total += myval;
                }
               // System.out.println("\t\tN="+n+" ; rsq=" + rsq[indexA][indexB] + " ; total = " + total + " from " + a + " and " + b);
            }
        }
        if(n>0) {
            //System.out.println("\tReturnning average Rsq of " + (total/n));
            return total / n;
        }else{
            return Float.NaN;   //Return NaN if no SNPs counted
        }
    }
}
