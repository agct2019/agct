package test;

import gene.Gene;

import javax.swing.*;
import java.awt.*;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

public class Statistics implements Debuggable {

    public static String TOKEN_ALL_TESTS = "ALL_CHI2";
    public static String allTests;

    public static Vector annotationsP, annotationsF, annotationsC;
    // all annotations that pass the tests

    public static Vector annotationsPChi2, annotationsFChi2, annotationsCChi2;
    // all annotations with their Chi2 P-value

    public static Vector annotationsBH;
    // each element = Vetor with 1- Flag (P:, F:, C:), 2- Annotation, 3-
    // Clustering Number, 4- Cluster number, 5- P-value, 6- P*value (BH)
    // Vector is sorted according to increasing P-value (P* is monotonous)

    public static boolean allChi2Tests;
    public static int refNumClustering;
    public static Color CHI2_COLOR = Color.cyan;

    public static int ngau = 18;
    public static double y[] = {0.0021695375159141994, 0.011413521097787704,
            0.027972308950302116, 0.051727015600492421, 0.082502225484340941,
            0.12007019910960293, 0.16415283300752470, 0.21442376986779355,
            0.27051082840644336, 0.33199876341447887, 0.39843234186401943,
            0.46931971407375483, 0.54413605556657973, 0.62232745288031077,
            0.70331500465597174, 0.78649910768313447, 0.87126389619061517,
            0.95698180152629142};

    public static double w[] = {0.0055657196642445571, 0.012915947284065419,
            0.020181515297735382, 0.027298621498568734, 0.034213810770299537,
            0.040875750923643261, 0.047235083490265582, 0.053244713977759692,
            0.058860144245324798, 0.064039797355015485, 0.068745323835736408,
            0.072941885005653087, 0.076598410645870640, 0.079687828912071670,
            0.082187266704339706, 0.084078218979661945, 0.085346685739338721,
            0.085983275670394821};

    static final int SWITCH = 3000;

    static final double EPS = 6E-8;
    static final double FPMIN = 1E-30;
    static int ASWITCH = 100;
    static double gln;

    public static double LIMIT_P_MAGNITUDE = 0.20;
    public static double LIMIT_P_CHI2 = 0.01;
    private static double LIMIT_P_DELAUNAY = 0.20;

    public static double getLIMIT_P_DELAUNAY() {
        return LIMIT_P_DELAUNAY;
    }

    public static void setLIMIT_P_DELAUNAY(double lIMITPDELAUNAY) {
        LIMIT_P_DELAUNAY = lIMITPDELAUNAY;
    }

    public static boolean isPositiveDelauney(Gene gene, int neighborId) {
        double nangle = ((Double) gene.neighbor_angles.get(neighborId))
                .doubleValue();
        return nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY / 2.0;
    }


    public static boolean isNegativeDelauney(Gene gene, int neighborId) {
        double nangle = ((Double) gene.neighbor_angles.get(neighborId))
                .doubleValue();
        return nangle > Math.PI * (1.0 - Statistics.LIMIT_P_DELAUNAY / 2.0);
    }

//	private static double THR = Double.NEGATIVE_INFINITY;

//	public static boolean okPlot_Magnitude(Gene gene) {
//		if (THR == Double.NEGATIVE_INFINITY) {
//			ArrayList<Double> tmp = new ArrayList<Double>();
//			for (int j = 0; j < AGCT.getInstance().data.getMyDomain().numberSelectedGenes; j++) {
//				Gene _gene = (Gene) AGCT.getInstance().data
//						.getMyDomain()
//						.getGenes()
//						.get(
//								AGCT.getInstance().data.getMyDomain().selectedGeneNumberToGeneNumber[j]);
//				tmp.add(_gene.totalDistance);
//			}
//			THR = tmp.get((int) Math.round((tmp.size() - 1)
//					* (1 - Statistics.LIMIT_P_MAGNITUDE)));
//		}
//		if (gene.totalDistance >= THR) {
//			for (int k = 0; k < gene.neighbors.size(); k++) {
//				if (Statistics.isNegativeDelauney(gene, k)
//						|| Statistics.isPositiveDelauney(gene, k))
//					return true;
//			}
//			return false;
//		}
//		return false;
//	}

    static final double[] yyy = {0.0021695375159141994, 0.011413521097787704,
            0.027972308950302116, 0.051727015600492421, 0.082502225484340941,
            0.12007019910960293, 0.16415283300752470, 0.21442376986779355,
            0.27051082840644336, 0.33199876341447887, 0.39843234186401943,
            0.46931971407375483, 0.54413605556657973, 0.62232745288031077,
            0.70331500465597174, 0.78649910768313447, 0.87126389619061517,
            0.95698180152629142};

    static final double[] www = {0.0055657196642445571, 0.012915947284065419,
            0.020181515297735382, 0.027298621498568734, 0.034213810770299537,
            0.040875750923643261, 0.047235083490265582, 0.053244713977759692,
            0.058860144245324798, 0.064039797355015485, 0.068745323835736408,
            0.072941885005653087, 0.076598410645870640, 0.079687828912071670,
            0.082187266704339706, 0.084078218979661945, 0.085346685739338721,
            0.085983275670394821};

    private static double gser(double a, double x) {
        double sum, del, ap;
        gln = gammln(a);
        ap = a;
        del = sum = 1.0 / a;
        for (; ; ) {
            ++ap;
            del *= x / ap;
            sum += del;
            if (Math.abs(del) < Math.abs(sum) * EPS) {
                return sum * Math.exp(-x + a * Math.log(x) - gln);
            }
        }
    }

    private static void initAll() {
        annotationsP = annotationsF = annotationsC = null;
        annotationsPChi2 = annotationsFChi2 = annotationsCChi2 = null;
        annotationsBH = null;
        allChi2Tests = false;
        refNumClustering = -1;
    }

    /**
     * 不完全ガンマ関数
     *
     * @param a
     * @param x
     * @return 1-P(a,x)
     */
    private static double gammq(double a, double x) {
        if (x < 0.0 || a <= 0.0)
            Matrix.perror("bad args in incomplete gamma");
        if (x == 0.0)
            return 1.0;
        else if ((int) a >= ASWITCH)
            return gammpapprox(a, x, 0);
        else if (x < a + 1.0)
            return 1.0 - gser(a, x);
        else
            return gcf(a, x);
    }

    private static double gcf(double a, double x) {
        int i;
        double an, b, c, d, del, h;
        gln = gammln(a);
        b = x + 1.0 - a;
        c = 1.0 / FPMIN;
        d = 1.0 / b;
        h = d;
        for (i = 1; ; i++) {
            an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (Math.abs(d) < FPMIN)
                d = FPMIN;
            c = b + an / c;
            if (Math.abs(c) < FPMIN)
                c = FPMIN;
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) <= EPS)
                break;
        }
        return Math.exp(-x + a * Math.log(x) - gln) * h;
    }

    private static double gammpapprox(double a, double x, int psig) {
        int j;
        double xu, t, sum, ans, ret;
        double a1 = a - 1.0, lna1 = Math.log(a1), sqrta1 = Math.sqrt(a1);
        gln = gammln(a);
        if (x > a1)
            xu = Math.max(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
        else
            xu = Math.max(0., Math.min(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));
        sum = 0;
        for (j = 0; j < ngau; j++) {
            t = x + (xu - x) * y[j];
            sum += w[j] * Math.exp(-(t - a1) + a1 * (Math.log(t) - lna1));
        }
        ans = sum * (xu - x) * Math.exp(a1 * (lna1 - 1.) - gln);
        ret = (psig == 1) ? (ans > 0.0 ? 1.0 - ans : -ans) : (ans >= 0.0 ? ans
                : 1.0 + ans);

        return ret;
    }

    private static boolean contains(String large, String small) {
        if (small.length() > large.length())
            return false;
        int i;
        for (i = 0; i <= large.length() - small.length(); i++)
            if (small.equals(large.substring(i, i + small.length())))
                return true;
        return false;
    }

    public static void flushAll() {
        allTests = "";
        annotationsP = annotationsF = annotationsC = null;
        annotationsPChi2 = annotationsFChi2 = annotationsCChi2 = null;
        annotationsBH = null;
        allChi2Tests = false;
    }

    public static void flushAll(AGCT ap) {
        flushAll();
        ap.saveButtonStats.setEnabled(false);
    }

    public static void getChi2(JCheckBox sp, JCheckBox sf, JCheckBox sc,
                               JTextField tf, AGCT ag, int nc) {
        String ref = tf.getText();
        if (ref.equals(TOKEN_ALL_TESTS))
            getAllSignificantChi2(sp, sf, sc, ag, nc);
        else {
            //			getSinglePatternChi2(sp, sf, sc, tf, ag, nc);//なにもしない
        }
    }

    /*
        private static void getSinglePatternChi2(JCheckBox sp, JCheckBox sf,
                JCheckBox sc, JTextField tf, AGCT ag, int nc) {
            AGCTClustering_Algorithm algo = ag.getClustering(nc).myClusteringAlgorithm;
            Domain dom = ag.data.getMyDomain();
            Vector gano;
            int i, j, k;
            String ref = tf.getText(), label, str, curstr, dum, sdep, s0, s1, s2, s3, s4, s5, s6, p1, p2, p3, p4, summary = "Summary: ", insidec = "+ (tag inside cluster) ", outsidec = "- (tag outside cluster) ", indifferentc = "/ (indifferent tag) ", sitc = "";
            int cluti = -1, patti = -1, len, l1, l2, l3, l4;

            double[][] table = { { 0.0, 0.0 }, { 0.0, 0.0 } };
            double num, den, dchi, prob, ad, bc;
    //
    //		  Organisation of table: Cluster -> | YES | NO
    //
    //		  Pattern | v YES NO
    //

            sdep = "Clustering #" + nc + "\n";
            sdep += algo + "\n";
            sdep += "Contingency tables for pattern: " + ref;
            allTests += sdep + "\n";
            System.out.println(sdep);

            for (j = 0; j < algo.getNumberOfClusters(); j++) {
                table[0][0] = table[0][1] = table[1][0] = table[1][1] = 0.0;
                for (i = 0; i < dom.numberSelectedGenes; i++)
                    fillTable(algo, sp, sf, sc, ref, dom, table, i, j);

                dchi = chi2Yates(table);

                if (dchi >= 0.0) {

                    prob = gammq(0.5, 0.5 * dchi);

                    curstr = "";
                    l1 = (new String("" + table[0][0])).length();
                    l2 = (new String("" + table[0][1])).length();
                    l3 = (new String("" + table[1][0])).length();
                    l4 = (new String("" + table[1][1])).length();
                    len = (int) Math.max(l1, l2);
                    len = (int) Math.max(len, l3);
                    len = (int) Math.max(len, l4);

                    ad = table[0][0] * table[1][1];
                    bc = table[0][1] * table[1][0];

                    p1 = p2 = p3 = p4 = "";
                    for (k = 0; k < len - l1; k++)
                        p1 += " ";
                    for (k = 0; k < len - l2; k++)
                        p2 += " ";
                    for (k = 0; k < len - l3; k++)
                        p3 += " ";
                    for (k = 0; k < len - l4; k++)
                        p4 += " ";

                    len += 2;
                    if (len < 5)
                        len = 5;
                    for (k = 0; k < len; k++)
                        curstr += "-";

                    s0 = "2x2 contingency table for cluster " + j + ":";
                    s1 = "+" + curstr + "+" + curstr + "+";
                    s2 = "| " + p1 + table[0][0] + " | " + p2 + table[0][1] + " |";
                    s3 = "+" + curstr + "+" + curstr + "+";
                    s4 = "| " + p3 + table[1][0] + " | " + p4 + table[1][1] + " |";
                    s5 = "+" + curstr + "+" + curstr + "+";
                    s6 = " ** Prob[Chi2 > " + dchi + " ]= " + prob + " --> ";

                    sitc = "";
                    if (prob < LIMIT_P_CHI2) {
                        s6 += "I reject the independence assumption, tendency ";
                        if (ad > bc)
                            sitc = insidec;
                        else if (ad < bc)
                            sitc = outsidec;
                        else
                            sitc = indifferentc;
                        s6 += sitc;
                    } else {
                        s6 += "I keep the independence assumption ";
                        sitc = "none";
                    }
                    s6 += "(LIMIT_P_CHI2 = " + LIMIT_P_CHI2 + ")";

                    System.out.println(s0);
                    System.out.println(s1);
                    System.out.println(s2);
                    System.out.println(s3);
                    System.out.println(s4);
                    System.out.println(s5);
                    System.out.println(s6 + "\n");

                    allTests += (s0 + "\n" + s1 + "\n" + s2 + "\n" + s3 + "\n" + s4 + "\n"
                            + s5 + "\n" + s6 + "\n" + "\n");

                    summary += "[ " + j + " : " + sitc + "]";
                }
            }
            allTests += summary;
            System.out.println(summary);
        }
    */
    public static Hashtable genesToCluster(AGCTClustering_Algorithm algo,
                                           Domain dom) {
        Hashtable ret = new Hashtable();
        Vector l;
        int i, cl, k;
        for (i = 0; i < dom.numberSelectedGenes; i++) {
            cl = algo.majorityCluster(i);
            ret.put(new Integer(dom.selectedGeneNumberToGeneNumber[i]), new Integer(
                    cl));
            if (!Prototype.No_Reduction) {
                l = (Vector) Prototype.Closest_Center_To_Cluster_Points
                        .get(new Integer(dom.selectedGeneNumberToGeneNumber[i]));
                for (k = 0; k < l.size(); k++)
                    ret.put((Integer) l.get(k), new Integer(cl));
            }
        }
        return ret;
    }

    private static void fastFillTableExact(AGCTClustering_Algorithm algo,
                                           boolean sp, boolean sf, boolean sc, String ref, Domain dom,
                                           double table[][], int j, Hashtable gtc) {
        // makes exact pattern matching between annotations (faster)

        int i;
        int cluti = -1, patti = -1, k;
        Integer uu;
        Hashtable geneTags = null;
        if ((sf && sc) || (sf && sp) || (sc && sp))
            Matrix.perror("Statistics.class :: two out of three sp, sf, sc are true");

        Hashtable rig = null;
        Vector rigv;
        if (sf)
            rig = dom.annotationFToGenes;
        else if (sp)
            rig = dom.annotationPToGenes;
        else if (sc)
            rig = dom.annotationCToGenes;

        if (!rig.containsKey(ref))
            Matrix.perror("Statistics.class :: I do not find annotations for String "
                    + ref);
        rigv = (Vector) rig.get(ref);
        geneTags = new Hashtable();

        patti = 0;
        for (i = 0; i < rigv.size(); i++) {
            //			DoDebug.debug(rigv.get(i));
            if (!gtc.containsKey((Integer) rigv.get(i)))
                continue;
            k = ((Integer) gtc.get((Integer) rigv.get(i))).intValue();
            if (k == j)
                cluti = 0;
            else
                cluti = 1;
            table[patti][cluti] += 1.0;
            geneTags.put((Integer) rigv.get(i), new Boolean(true));
        }
        patti = 1;
        Enumeration extensions = gtc.keys();
        while (extensions.hasMoreElements()) {
            uu = (Integer) extensions.nextElement();
            if (!geneTags.containsKey(uu)) {
                k = ((Integer) gtc.get(uu)).intValue();
                if (k == j)
                    cluti = 0;
                else
                    cluti = 1;
                table[patti][cluti] += 1.0;
            }
        }

        geneTags = null;
    }

    public static void fillTable(AGCTClustering_Algorithm algo, JCheckBox sp,
                                 JCheckBox sf, JCheckBox sc, String ref, Domain dom, double table[][],
                                 int i, int j) {
        Gene gano;
        int cluti = -1, patti = -1, k;
        Vector l;
        if (algo.majorityClusterFor(i, j))
            cluti = 0;
        else
            cluti = 1;

        if (Prototype.No_Reduction) {
            gano = (Gene) dom.getGenes().get(dom.selectedGeneNumberToGeneNumber[i]);
            // if ( (gano.annotationsP != null) || (gano.annotationsF != null)
            // || (gano.annotationsC != null) )
            fillTableForGene(sp, sf, sc, ref, gano, table, cluti);
        } else {
            l = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(
                    dom.selectedGeneNumberToGeneNumber[i]));
            for (k = 0; k < l.size(); k++) {
                gano = (Gene) dom.getGenes().get(((Integer) l.get(k)).intValue());
                // if ( (gano.annotationsP != null) || (gano.annotationsF !=
                // null) || (gano.annotationsC != null) )
                fillTableForGene(sp, sf, sc, ref, gano, table, cluti);
            }
        }
    }

    public static void fillTableForGene(JCheckBox sp, JCheckBox sf, JCheckBox sc,
                                        String ref, Gene gg, double[][] table, int cluti) {
        int k, patti = -1;
        String label, str;
        Vector gano;
        boolean found;
        if ((sf != null) && (sf.isSelected())) {
            if (gg.annotationsF != null) {
                gano = gg.annotationsF;
                k = 0;
                found = false;
                do {
                    label = ((String[]) gano.get(k))[0];
                    if (!label.equals(Domain.ANNOTATION_F))
                        Matrix
                                .perror("Statistics.class :: annotation != F in F annotations");
                    str = ((String[]) gano.get(k))[1];
                    if (contains(str, ref))
                        found = true;
                    else
                        k++;
                } while ((!found) && (k < gano.size()));
                if (found)
                    patti = 0;
                else
                    patti = 1;
            } else
                patti = 1;
            table[patti][cluti] += 1.0;

			/*
             * for (k=0;k<gano.size();k++){ label = ((String[])gano.get(k))[0];
			 * str = ((String[])gano.get(k))[1]; if
			 * (!label.equals(Domain.ANNOTATION_F))
			 * Matrix.perror("Statistics.class :: annotation != F in F annotations"
			 * ); if (contains(str,ref)) patti = 0; else patti = 1;
			 * table[patti][cluti] += 1.0; }
			 */
        }
        if ((sp != null) && (sp.isSelected())) {
            if (gg.annotationsP != null) {
                gano = gg.annotationsP;
                k = 0;
                found = false;
                do {
                    label = ((String[]) gano.get(k))[0];
                    if (!label.equals(Domain.ANNOTATION_P))
                        Matrix
                                .perror("Statistics.class :: annotation != P in P annotations");
                    str = ((String[]) gano.get(k))[1];
                    if (contains(str, ref))
                        found = true;
                    else
                        k++;
                } while ((!found) && (k < gano.size()));
                if (found)
                    patti = 0;
                else
                    patti = 1;
            } else
                patti = 1;
            table[patti][cluti] += 1.0;

			/*
             * for (k=0;k<gano.size();k++){ label = ((String[])gano.get(k))[0];
			 * str = ((String[])gano.get(k))[1]; if
			 * (!label.equals(Domain.ANNOTATION_P))
			 * Matrix.perror("Statistics.class :: annotation != P in P annotations"
			 * ); if (contains(str,ref)) patti = 0; else patti = 1;
			 * table[patti][cluti] += 1.0; }
			 */
        }
        if ((sc != null) && (sc.isSelected())) {
            if (gg.annotationsC != null) {
                gano = gg.annotationsC;
                k = 0;
                found = false;
                do {
                    label = ((String[]) gano.get(k))[0];
                    if (!label.equals(Domain.ANNOTATION_C))
                        Matrix
                                .perror("Statistics.class :: annotation != C in C annotations");
                    str = ((String[]) gano.get(k))[1];
                    if (contains(str, ref))
                        found = true;
                    else
                        k++;
                } while ((!found) && (k < gano.size()));
                if (found)
                    patti = 0;
                else
                    patti = 1;
            } else
                patti = 1;
            table[patti][cluti] += 1.0;

			/*
			 * for (k=0;k<gano.size();k++){ label = ((String[])gano.get(k))[0];
			 * str = ((String[])gano.get(k))[1]; if
			 * (!label.equals(Domain.ANNOTATION_C))
			 * Matrix.perror("Statistics.class :: annotation != C in C annotations"
			 * ); if (contains(str,ref)){ patti = 0; }else patti = 1;
			 * table[patti][cluti] += 1.0;
			 *
			 * System.out.println(gg.name + ", " + k }
			 */
        }

    }

    public static double chi2Yates(double[][] table) {
        if ((table[0][0] < 0.0) || (table[0][1] < 0.0) || (table[1][0] < 0.0)
                || (table[1][1] < 0.0))
            Matrix.perror("Chi square with < 0 values in contingency table");
        double val = 0.0, tot, di1, di2, sum1, sum2, num, den, d1, d2, d3, d4, dchi;

        tot = table[0][0] + table[0][1] + table[1][0] + table[1][1];

        // Yates correction for ChiSquare

        di1 = table[0][0] * table[1][1];
        di2 = table[0][1] * table[1][0];

        sum1 = table[0][0] + table[1][1];
        sum2 = table[0][1] + table[1][0];

        num = tot * (Math.abs(di1 - di2) - (tot / 2.0))
                * (Math.abs(di1 - di2) - (tot / 2.0));

        d1 = table[0][0] + table[0][1];
        d2 = table[0][0] + table[1][0];
        d3 = table[1][0] + table[1][1];
        d4 = table[0][1] + table[1][1];

        if ((d1 == 0.0) || (d2 == 0.0) || (d3 == 0.0) || (d4 == 0.0)) {
            System.out
                    .println("Statistics.class warning :: a Chi2 table is degenerated");
            dchi = 0.0;
        } else {
            den = d1 * d2 * d3 * d4;

            if ((d1 != 0.0) && (d2 != 0.0) && (d3 != 0.0) && (d4 != 0.0))
                dchi = num / den;
            else
                dchi = -1.0;
        }

        if (dchi < 0.0)
            Matrix.perror("Statistics.class :: negative Chi2 distance");

        return dchi;
    }

    public static void getAllSignificantChi2(final JCheckBox sp,
                                             final JCheckBox sf, final JCheckBox sc, final AGCT ag, final int nc) {
        Thread t = new Thread() {
            public void run() {
                refNumClustering = nc;

                AGCTClustering_Algorithm algo = ag.getClustering(nc).myClusteringAlgorithm;
                Domain dom = ag.data.getMyDomain();
                Vector gano, allStrings, stringsByCluster = new Vector(), curStrings, vTable;
                int i, j, k, tp, uu, numclu, ntot = 0, cardtot = 0;
                String ref, label, str, curstr, dum, sdep, resString, saff, allBH = "";
                int cluti = -1, patti = -1;

                double[][] table = {{0.0, 0.0}, {0.0, 0.0}};
                double num, den, dchi, prob;
                boolean goForIt = false, bsp = false, bsf = false, bsc = false;

                annotationsP = annotationsF = annotationsC = null;
                annotationsPChi2 = annotationsFChi2 = annotationsCChi2 = null;
                Vector dumAnnotations = null, dumAnnotationsChi2 = null, littleVector;
                allChi2Tests = false;

                Hashtable gtc = genesToCluster(algo, dom);

				/*
				 * Organisation of table: Cluster -> | YES | NO
				 *
				 * Pattern | v YES NO
				 */

                sdep = "Clustering #" + nc + "\n";
                sdep += algo + "\n";
                sdep += "All significant Chi-square tests (p_limit = " + LIMIT_P_CHI2
                        + ")";
                allTests += sdep + "\n";
//				System.out.println(sdep);

                for (j = 0; j < algo.getNumberOfClusters(); j++)
                    stringsByCluster.addElement(new Vector());

                if ((sf != null) && (sf.isSelected()))
                    cardtot += ag.data.getMyDomain().selectedGenesAllAnnotationsF.size();
                if ((sp != null) && (sp.isSelected()))
                    cardtot += ag.data.getMyDomain().selectedGenesAllAnnotationsP.size();
                if ((sc != null) && (sc.isSelected()))
                    cardtot += ag.data.getMyDomain().selectedGenesAllAnnotationsC.size();

                AGCTCounter cc = new AGCTCounter("Computing all Chi2", cardtot);
                for (tp = 0; tp < 3; tp++) {
                    allStrings = null;
                    resString = null;
                    if ((tp == 0) && (sf != null) && (sf.isSelected())) {
                        goForIt = true;
                        allStrings = ag.data.getMyDomain().selectedGenesAllAnnotationsF;
                        resString = Domain.ANNOTATION_F;
                        annotationsF = new Vector();
                        annotationsFChi2 = new Vector();
                        for (j = 0; j < algo.getNumberOfClusters(); j++)
                            annotationsF.addElement(new Vector());
                        dumAnnotations = annotationsF;
                        dumAnnotationsChi2 = annotationsFChi2;
                        bsf = true;
                        bsp = bsc = false;
                    } else if ((tp == 1) && (sp != null) && (sp.isSelected())) {
                        goForIt = true;
                        allStrings = ag.data.getMyDomain().selectedGenesAllAnnotationsP;
                        resString = Domain.ANNOTATION_P;
                        annotationsP = new Vector();
                        annotationsPChi2 = new Vector();
                        for (j = 0; j < algo.getNumberOfClusters(); j++)
                            annotationsP.addElement(new Vector());
                        dumAnnotations = annotationsP;
                        dumAnnotationsChi2 = annotationsPChi2;
                        bsp = true;
                        bsf = bsc = false;
                    } else if ((tp == 2) && (sc != null) && (sc.isSelected())) {
                        goForIt = true;
                        allStrings = ag.data.getMyDomain().selectedGenesAllAnnotationsC;
                        resString = Domain.ANNOTATION_C;
                        annotationsC = new Vector();
                        annotationsCChi2 = new Vector();
                        for (j = 0; j < algo.getNumberOfClusters(); j++)
                            annotationsC.addElement(new Vector());
                        dumAnnotations = annotationsC;
                        dumAnnotationsChi2 = annotationsCChi2;
                        bsc = true;
                        bsp = bsf = false;
                    } else
                        goForIt = false;

                    if (goForIt) {
                        for (uu = 0; uu < allStrings.size(); uu++) {
                            cc.increment();

                            ref = (String) allStrings.get(uu);
                            numclu = 0;
                            for (j = 0; j < algo.getNumberOfClusters(); j++) {
                                curStrings = (Vector) stringsByCluster.get(j);

                                table[0][0] = table[0][1] = table[1][0] = table[1][1] = 0.0;
                                // for (i=0;i<dom.numberSelectedGenes;i++)
                                // fillTable(algo, sp, sf, sc, ref, dom, table,
                                // i, j);

                                // table2[0][0] = table2[0][1] = table2[1][0] =
                                // table2[1][1] = 0.0;
                                fastFillTableExact(algo, bsp, bsf, bsc, ref, dom, table, j, gtc);

                                // System.out.println("\n-----\nString " + ref +
                                // "\nTable : \n" + table[0][0] + " " +
                                // table[0][1] + "\n" + table[1][0] + " " +
                                // table[1][1] + "\n");
                                // System.out.println("Table2 : \n" +
                                // table2[0][0] + " " + table2[0][1] + "\n" +
                                // table2[1][0] + " " + table2[1][1] + "\n\n");

                                dchi = chi2Yates(table);

                                if (dchi >= 0.0) {
                                    prob = gammq(0.5, 0.5 * dchi);

                                    littleVector = new Vector();
                                    littleVector.addElement(new String(resString));
                                    littleVector.addElement(ref);
                                    littleVector.addElement(new Integer(nc));
                                    littleVector.addElement(new Integer(j));
                                    littleVector.addElement(new Double(prob));

                                    vTable = new Vector();
                                    vTable.addElement(new Double(table[0][0]));
                                    vTable.addElement(new Double(table[0][1]));
                                    vTable.addElement(new Double(table[1][0]));
                                    vTable.addElement(new Double(table[1][1]));

                                    littleVector.addElement(vTable);

                                    dumAnnotationsChi2.addElement(littleVector);

//									if (prob < LIMIT_P_CHI2) {
                                    //										ntot++;
                                    //										if (numclu == 0) {
                                    //											saff = "\nString " + ref + " (" + resString + ") : \n";
                                    //											allTests += saff;
                                    //											// System.out.print(saff);
                                    //										}
                                    //										saff = " cluster " + j + " : chi = " + dchi + " ; p-value = " + prob + "\n";
                                    //										allTests += saff;
                                    //										// System.out.print(saff);
                                    //										numclu++;
                                    //
                                    //										((Vector) dumAnnotations.get(j)).addElement(ref);
                                    //										curStrings.addElement(ref);
//									}
                                }
                            }
                        }
                    }
                }
                cc.end();

                if (ntot == 0) {
                } else {
                }
                allTests += "\n\n";

                fillAnnotationsBH(nc);

                for (i = annotationsBH.size() - 1; i >= 0; i--) {
                    littleVector = (Vector) Statistics.annotationsBH.get(i);
//					System.out.println(littleVector.get(0) + "\t" + littleVector.get(1)
//							+ "\t" + littleVector.get(2) + "\t" + littleVector.get(3) + "\t"
//							+ littleVector.get(4) + "\t" + littleVector.get(5));
                }

                allChi2Tests = true;
                ag.saveButtonStats.setEnabled(true);
            }
        };
        t.start();
    }

    public static void fillAnnotationsBH(int nc) {
        int i, t, siz;
        annotationsBH = new Vector();
        Vector dumVector = null, littleVector, lv;
        double cps = -1.0, ps;
        for (t = 0; t < 3; t++) {
            if (t == 0)
                dumVector = annotationsFChi2;
            else if (t == 1)
                dumVector = annotationsPChi2;
            else if (t == 2)
                dumVector = annotationsCChi2;
            if (dumVector != null) {
                for (i = 0; i < dumVector.size(); i++) {
                    lv = (Vector) dumVector.get(i);
                    littleVector = new Vector();
                    littleVector.addElement((String) lv.get(0));
                    littleVector.addElement((String) lv.get(1));
                    littleVector.addElement((Integer) lv.get(2));
                    littleVector.addElement((Integer) lv.get(3));
                    littleVector.addElement((Double) lv.get(4));
                    littleVector.addElement((Vector) lv.get(5));
                    annotationsBH.addElement(littleVector);
                }
            }
        }
        QuickSort.quicksortVectorDouble(annotationsBH, 4);
        siz = annotationsBH.size();
        for (i = 0; i < siz; i++) {
            littleVector = (Vector) annotationsBH.get(i);
            ps = ((Double) littleVector.get(4)).doubleValue();
            ps *= (double) siz;
            ps /= (double) (i + 1);
            if (ps >= 1.0)
                ps = 1.0;
            if ((i == 0) || (ps > cps))
                cps = ps;
            littleVector.addElement(new Double(cps));
        }
    }

    public static boolean containsInBH(String d, String ref) {
        int i = 0;
        String dumRef, dumString;
        if ((!allChi2Tests) || (annotationsBH == null)
                || (annotationsBH.size() == 0))
            return true;
        else {
            while ((i < annotationsBH.size())
                    && (((Double) ((Vector) annotationsBH.get(i)).get(5)).doubleValue() < LIMIT_P_CHI2)) {
                dumRef = (String) ((Vector) annotationsBH.get(i)).get(0);
                dumString = (String) ((Vector) annotationsBH.get(i)).get(1);
                if ((dumRef.equals(ref)) && (dumString.equals(d)))
                    return true;
                i++;
            }
            return false;
        }
    }

    public static double SQR(double x) {
        return x * x;
    }

    private static double gammln(double xx) {
        int j;
        double x, tmp, y, ser;
        double[] cof = {57.1562356658629235, -59.5979603554754912,
                14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
                .465236289270485756e-4, -.983744753048795646e-4,
                .158088703224912494e-3, -.210264441724104883e-3,
                .217439618115212643e-3, -.164318106536763890e-3,
                .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5};
        if (xx <= 0)
            Matrix.perror("bad arg in gammln");
        y = x = xx;
        tmp = x + 5.24218750000000000;
        tmp = (x + 0.5) * Math.log(tmp) - tmp;
        ser = 0.999999999999997092;
        for (j = 0; j < 14; j++)
            ser += cof[j] / ++y;
        return tmp + Math.log(2.5066282746310005 * ser / x);
    }

    public static double betaiapprox(double a, double b, double x) {
        int j;
        double xu, t, sum, ans;
        double a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);
        double lnmu = Math.log(mu), lnmuc = Math.log(1. - mu);
        t = Math.sqrt(a * b / ((a + b) * (a + b) * (a + b + 1.0)));
        if (x > a / (a + b)) {
            if (x >= 1.0)
                return 1.0;
            xu = Math.min(1., Math.max(mu + 10. * t, x + 5.0 * t));
        } else {
            if (x <= 0.0)
                return 0.0;
            xu = Math.max(0., Math.min(mu - 10. * t, x - 5.0 * t));
        }
        sum = 0;
        for (j = 0; j < 18; j++) {
            t = x + (xu - x) * yyy[j];
            sum += www[j]
                    * Math
                    .exp(a1 * (Math.log(t) - lnmu) + b1 * (Math.log(1 - t) - lnmuc));
        }
        ans = sum
                * (xu - x)
                * Math.exp(a1 * lnmu - gammln(a) + b1 * lnmuc - gammln(b)
                + gammln(a + b));
        return ans > 0.0 ? 1.0 - ans : -ans;
    }

    public static double betacf(double a, double b, double x) {
        int m, m2;
        double aa, c, d, del, h, qab, qam, qap;
        qab = a + b;
        qap = a + 1.0;
        qam = a - 1.0;
        c = 1.0;
        d = 1.0 - qab * x / qap;
        if (Math.abs(d) < FPMIN)
            d = FPMIN;
        d = 1.0 / d;
        h = d;
        for (m = 1; m < 10000; m++) {
            m2 = 2 * m;
            aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < FPMIN)
                d = FPMIN;
            c = 1.0 + aa / c;
            if (Math.abs(c) < FPMIN)
                c = FPMIN;
            d = 1.0 / d;
            h *= d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < FPMIN)
                d = FPMIN;
            c = 1.0 + aa / c;
            if (Math.abs(c) < FPMIN)
                c = FPMIN;
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) <= EPS)
                break;
        }
        return h;
    }

    public static double betai(double a, double b, double x) {
        double bt;
        if (a <= 0.0 || b <= 0.0)
            Matrix.perror("Bad a or b in routine betai");
        if (x < 0.0 || x > 1.0)
            Matrix.perror("Bad x in routine betai");
        if (x == 0.0 || x == 1.0)
            return x;
        if (a > SWITCH && b > SWITCH)
            return betaiapprox(a, b, x);
        bt = Math.exp(gammln(a + b) - gammln(a) - gammln(b) + a * Math.log(x) + b
                * Math.log(1.0 - x));
        if (x < (a + 1.0) / (a + b + 2.0))
            return bt * betacf(a, b, x) / a;
        else
            return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
    }

    /**
     * dataの算術平均と不偏分散が、それぞれ、aveETvar[0],aveETvar[1]に入る
     *
     * @param data
     * @param aveETvar
     */
    public static void avevar(double[] data, double[] aveETvar) {
        double s, ep, ave, var;
        int j, n = data.length;
        ave = 0.0;
        for (j = 0; j < n; j++)
            ave += data[j];
        ave /= n; // dataの算術平均
        var = ep = 0.0;
        for (j = 0; j < n; j++) {
            s = data[j] - ave;
            ep += s;
            var += s * s;
        }
        var = (var - ep * ep / n) / (n - 1);// two pass algorithm
        // http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Two-pass_algorithm
        aveETvar[0] = ave;
        aveETvar[1] = var;
    }

    public static double tutest(double[] data1, double[] data2) {
        double df, t, prob, ave1, ave2, var1, var2;
        double avar1[] = new double[2];
        double avar2[] = new double[2];
        int n1 = data1.length, n2 = data2.length;
        avevar(data1, avar1);
        avevar(data2, avar2);
        ave1 = avar1[0];
        var1 = avar1[1];
        ave2 = avar2[0];
        var2 = avar2[1];
        t = (ave1 - ave2) / Math.sqrt(var1 / n1 + var2 / n2);
        df = SQR(var1 / n1 + var2 / n2)
                / (SQR(var1 / n1) / (n1 - 1) + SQR(var2 / n2) / (n2 - 1));
        prob = betai(0.5 * df, 0.5, df / (df + SQR(t)));

        return prob;
    }
}