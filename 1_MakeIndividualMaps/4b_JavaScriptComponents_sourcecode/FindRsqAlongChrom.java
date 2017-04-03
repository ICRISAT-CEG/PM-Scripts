package Misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.pal.statistics.FisherExact;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Utils;

/**
 *
 * @author jgw
 */
public class FindRsqAlongChrom {

	public static void main(String[] args) {
		String infile=args[0], anchorfile=args[1], outfile=args[2];
		FindRsqAlongChrom(infile, anchorfile, outfile);
	}

	public static void FindRsqAlongChrom(String infile, String anchorfile, String outfile) {
		//Import alignment
		System.out.println("Loading file\n");
		Alignment a = ImportUtils.readFromHapmap(infile, null);
		a = BitAlignment.getHomozygousNucleotideInstance(a, true);
		System.out.println("\tLoaded " + a.getSequenceCount() + " taxa with " + a.getSiteCount() + " sites");

		//Import list of anchor sites
		System.out.println("Loading list of anchor sites");
		BufferedReader reader = Utils.getBufferedReader(anchorfile);
		String site;
		HashMap <String, Boolean> anchors = new HashMap<>();
		try {
			site = reader.readLine();
			while(site != null){
				site=site.trim();
				anchors.put(site, true);
				site = reader.readLine();
			}
		} catch (IOException ex) {
			Logger.getLogger(FindRsqAlongChrom.class.getName()).log(Level.SEVERE, null, ex);
		}
		System.out.println("\tTotal " + anchors.size() + " anchor sites loaded");

		//Go through and identify target sites (those not beginning with S)
		System.out.println("Identifying target sites");
		String[] sites = a.getSNPIDs();
		ArrayList<Integer> siteFlags = new ArrayList<Integer>();
		ArrayList<Integer> nonFlags = new ArrayList<Integer>();
		for (int i = 0; i < sites.length; i++) {
			if(anchors.containsKey(sites[i])){	//If site is an anchor site, add to list of (non-flagged) anchors. Otherwise flag as something to be tested
				nonFlags.add(i);
			}else{
				siteFlags.add(new Integer(i));
			}
		}
		Integer[] targets = siteFlags.toArray(new Integer[0]);
		Integer[] nonTargets = nonFlags.toArray(new Integer[0]);
		//nonTargets = Arrays.copyOfRange(nonTargets, 0, 10000);	//For debugging, limit to first 10000 sites
		System.out.println("Found " + targets.length + " targets and " + nonTargets.length + " nontargets");


		//For each one of these targets, find its rsq values against all non-target sites and save
		int myMinTaxaForEstimate = 20;
		System.out.println("Calculating R-squared for combinations (NOTE: Must have at least " + myMinTaxaForEstimate + " joint comparisons to get a value; otherwise get NA)");
		FisherExact myFisherExact = new FisherExact((2 * a.getSequenceCount()) + 10);
		double[][] rsq = new double[targets.length][nonTargets.length];	//to store Rsquared values
		for (int i = 0; i < targets.length; i++) {
			if (i % 10000 == 0) {
				System.out.println("\tProcessing target # " + i);
			}
			int sitenum = targets[i];
			for (int n = 0; n < nonTargets.length; n++) {
				int compnum = nonTargets[n];
				BitSet rMj = a.getAllelePresenceForAllTaxa(sitenum, 0);
				BitSet rMn = a.getAllelePresenceForAllTaxa(sitenum, 1);
				BitSet cMj = a.getAllelePresenceForAllTaxa(compnum, 0);
				BitSet cMn = a.getAllelePresenceForAllTaxa(compnum, 1);
				double[] results = getRsqForSitePair(rMj, rMn, cMj, cMn, 2, myMinTaxaForEstimate, -1.0f, myFisherExact);
				rsq[i][n] = results[0];
			}
		}

		//Print out as a large matrix to check
		System.out.println("Outputting results to " + outfile);
		BufferedWriter out = Utils.getBufferedWriter(outfile);
		DecimalFormat df = new DecimalFormat("#.###");
		try {
			//Header
			out.write("site");
			for (int i = 0; i < nonTargets.length; i++) {
				out.write("\t" + sites[nonTargets[i]]);	//Output the site name
			}
			out.newLine();
			//Body
			for (int i = 0; i < targets.length; i++) {
				out.write(sites[targets[i]]);
				for (int n = 0; n < rsq[i].length; n++) {
					String myString = df.format(rsq[i][n]);
					if (myString.equals("�")) { //Odd character that gets returned sometimes for unknown or other error
						myString = "NA";
					}
					out.write("\t" + myString);
				}
				out.newLine();
			}
			out.close();
		} catch (IOException ex) {
			System.out.println("Error writing output file");
		}

		/*
		//Go through and calculate probabilities, going from each side
		System.out.println("Finding best positions");
		double min = 0.01, max = Double.MAX_VALUE;	//Minimum let "probability" fall to and maximum it can climb to
		double[][] minProbs = new double[rsq.length][rsq[0].length];
		int[] bestPos = new int[rsq.length];
		Arrays.fill(bestPos, -1);
		for (int i = 0; i < rsq.length; i++) {
			if (i % 10000 == 0) {
				System.out.println("\tProcessing target # " + i);
			}
			double[] leftProb = new double[rsq[i].length];
			double[] rightProb = new double[rsq[i].length];

			double left = 1, right = 1;	//Current "probabilities"
			for (int n = 0; n < rsq[i].length; n++) {
				if(n>=1 && !a.getLocusName(nonTargets[n]).equals(a.getLocusName(nonTargets[n-1]))){	//If n is on a new chromosome, reset probabilities
					left=1;
				}
				if(n<rsq[i].length -1 && !a.getLocusName(nonTargets[n]).equals(a.getLocusName(nonTargets[n+1]))){	//If n is on a new chromosome, reset probabilities
					right=1;
				}
				left = calcProbFromRsq(left, Math.sqrt(0.5 + rsq[i][n]), min, max);	//Ad-hoc adjustments to Rsquared to make it behave like I want. Centers around 1 and brings closer
				right = calcProbFromRsq(right, Math.sqrt(0.5 + rsq[i][rsq[i].length - n - 1]), min, max);
				//System.out.println("Left is " +  left + " and right is " + right + " from rsq " + rsq[i][n] + " and " + rsq[i][rsq[i].length - n - 1]);
				leftProb[n] = left;
				rightProb[rsq[i].length - n - 1] = right;
			}
			double maxProb = 0;
			for (int n = 0; n < rsq[i].length; n++) {	//Determine best position by finding the maximum point when take the minimum of the two sides
				minProbs[i][n] = Math.min(leftProb[n], rightProb[n]);
				//System.out.println("Minprob at " + n + " is " + minProbs[i][n] + " from left " + leftProb[n] + " and right " +  rightProb[n] + " compared to maxProb " + maxProb);
				if (minProbs[i][n] > maxProb) {
					maxProb = minProbs[i][n];
					//System.out.println("i = "  + i + " and n = " + n);
					bestPos[i] = nonTargets[n];
				}
			}
		}


		out = Utils.getBufferedWriter(bestfile);
		try {
			//Header
			out.write("site\tbest_position\tneighbor_name\tchr\tpos");
			out.newLine();
			//Body
			for (int i = 0; i < bestPos.length; i++) {
				out.write(sites[targets[i]] + "\t" + bestPos[i] + "\t" + a.getSNPID(bestPos[i])+ "\t" + a.getLocusName(bestPos[i])+ "\t" + a.getPositionInLocus(bestPos[i]));
				out.newLine();
			}
			out.close();
		} catch (IOException ex) {
			Logger.getLogger(JasonPipeline.class.getName()).log(Level.SEVERE, null, ex);
		}*/

	}

	public static double calcProbFromRsq(double prevProb, double currentProb, double min, double max) {
		if (Double.isNaN(currentProb)) {	//If get an Na result, make it so it doesn't change anything
			currentProb = 1;
		}
		double result = currentProb * prevProb;
		if (result < min) {
			return min;
		}
		if (result > max) {
			return max;
		}
		return result;
	}



	public static double[] getRsqForSitePair(BitSet rMj, BitSet rMn, BitSet cMj, BitSet cMn,
			int minMinorCnt, int minCnt, float minR2, FisherExact myFisherExact) {
		//Taken from LinkageDisequilibrium class and modified some
		//Array returns {Rsquare, Pvalue};
		double[] values = {Double.NaN, Double.NaN};
		int n = 0;
		int[][] contig = new int[2][2];
		n += contig[1][1] = (int) OpenBitSet.intersectionCount(rMn, cMn);
		n += contig[1][0] = (int) OpenBitSet.intersectionCount(rMn, cMj);
		if (contig[1][0] + contig[1][1] < minMinorCnt) {
			return values;
		}
		n += contig[0][1] = (int) OpenBitSet.intersectionCount(rMj, cMn);
		if (contig[0][1] + contig[1][1] < minMinorCnt) {
			return values;
		}
		n += contig[0][0] = (int) OpenBitSet.intersectionCount(rMj, cMj);
		if (n < minCnt) {
			return values;
		}
		double rValue = LinkageDisequilibrium.calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minCnt);
		if (Double.isNaN(rValue)) {
			return values;
		}
		values[0] = rValue;
		if (rValue < minR2) {
			return values;
		}
		values[1] = myFisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
		return values;
	}


}