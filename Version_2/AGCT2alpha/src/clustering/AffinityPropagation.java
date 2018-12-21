package clustering;

import forDebug.Debug;
import matrix.AbstMatrix;
import matrix.MatrixFile;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;

public class AffinityPropagation {
    class Pair {
        int id;
        double val;

        Pair(int id, double val) {
            this.id = id;
            this.val = val;
        }
    }

    private final double dampfact;
    private final int convits;
    private final int maxits;
    private boolean endCalc = false;
    private final int n;
    private final AbstMatrix W;// cloneしているので、自由に変更してよい
    private int iteration;
    private double netSimilarity = Double.NEGATIVE_INFINITY;

    public Cluster[] clustering() {
        int maxK = 0;// 最大いくつの要素を非 -∞ としてもつか
        assert !endCalc;
        double[][] ss = new double[n][n];
        for (int i = 0; i < n; i++) {
            int tmp = 0;
            for (int j = 0; j < n; j++) {
                if (W.get(i, j) != 0)
                    tmp++;
                if (i != j) {
                    ss[i][j] = W.get(i, j);
                }
            }
            maxK = Math.max(maxK, tmp);
        }

        Debug.debug("n", n);
        Debug.debug("maxK", maxK);
        double[] preference = new double[n];
        for (int i = 0; i < n; i++) {
            {
                double[] ds = new double[n - 1];
                for (int j = 0, k = 0; j < n; j++) {
                    if (i != j) {
                        ds[k++] = ss[i][j];
                    }
                }
                Arrays.sort(ds);
                ss[i][i] = ds[(n - 1) / 2];
                preference[i] = ds[(n - 1) / 2];
            }
        }
        double[][] as = new double[n][n];
        double[][] rs = new double[n][n];
        int[] exemplar = new int[n];
        int noChange = 0;
        for (iteration = 0; iteration < maxits; iteration++) {
            for (int i = 0; i < n; i++) {
                Pair first = null;
                Pair second = null;
                for (int k2 = 0; k2 < n; k2++) {
                    double now = as[i][k2] + ss[i][k2];
                    if (first == null || first.val < now) {
                        second = first;
                        first = new Pair(k2, now);
                    } else if (second == null || second.val < now) {
                        second = new Pair(k2, now);
                    }
                }
                double[] list = new double[n];
                for (int k = 0; k < n; k++) {
                    double max = first.id != k ? first.val : second.val;
                    rs[i][k] = rs[i][k] * dampfact + (ss[i][k] - max) * (1 - dampfact);
                    list[k] = rs[i][k];
                }
                Arrays.sort(list);
                double thr = list[n - maxK];
                for (int k = 0; k < n; k++) {//大きい方からmaxK個のみを残す
                    if (rs[i][k] < thr)
                        rs[i][k] = Double.NEGATIVE_INFINITY;
                }
            }
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < n; k++) {
                    if (i != k) {
                        double sum = 0;
                        for (int i2 = 0; i2 < n; i2++) {
                            if (i2 != i && i2 != k) {
                                sum += Math.max(0, rs[i2][k]);
                            }
                        }
                        as[i][k] = as[i][k] * dampfact + Math.min(0, rs[k][k] + sum)
                                * (1 - dampfact);
                    } else {
                        double sum = 0;
                        for (int i2 = 0; i2 < n; i2++) {
                            if (i2 != k) {
                                sum += Math.max(0, rs[i2][k]);
                            }
                        }
                        as[k][k] = as[k][k] * dampfact + sum * (1 - dampfact);
                    }
                }
            }

            boolean changed = false;
            for (int i = 0; i < n; i++) {
                int id = 0;
                double max = Double.NEGATIVE_INFINITY;
                for (int j = 0; j < n; j++) {
                    double val = as[i][j] + rs[i][j];
                    if (val > max) {
                        id = j;
                        max = val;
                    }
                }
                if (max > netSimilarity) {
                    netSimilarity = max;
                }

                if (exemplar[i] != id) {
                    changed = true;
                    exemplar[i] = id;
                }
            }
            if (!changed) {
                if (++noChange >= convits)
                    break;
            } else {
                noChange = 0;
            }
        }
        TreeSet<Integer> set = new TreeSet<Integer>();
        for (int i = 0; i < n; i++) {
            set.add(exemplar[i]);
        }
        Integer[] captains = set.toArray(new Integer[0]);
        HashMap<Integer, Integer> rev = new HashMap<Integer, Integer>();
        int k = set.size();
        for (int i = 0; i < k; i++) {
            rev.put(captains[i], i);
        }
        Cluster[] res = new Cluster[k];
        for (int i = 0; i < k; i++) {
            ArrayList<Integer> list = new ArrayList<Integer>();
            for (int j = 0; j < n; j++) {
                if (exemplar[j] == captains[i] && exemplar[j] != j) {
                    list.add(j);
                }
            }
            int[] member = new int[list.size()];
            for (int j = 0; j < list.size(); j++) {
                member[j] = list.get(j);
            }
            res[i] = new Cluster(captains[i], member);
        }
        endCalc = true;
        return res;
    }

    /*
     * format of Affinity Propagation Web Application
     */
    public String summary() {
        assert endCalc;
        return String.format("\n" + "maxits=%d\n" + "convits=%d\n"
                + "dampfact=%f\n" + "number of data points: %d\n"
                + "Preferences: from the input file\n" + "\n"
                + "Number of identified clusters: %d\n" + "Net similarity:\n"
                + "	 %f\n" + "<grey>Sum of data point-to-exemplar similarities:\n"
                + "	 %f\n" + "</font><grey>Sum of exemplar preferences:\n" + "	 %f\n"
                + "</font>Number of iterations: %d\n" + "\n"
                + "Start:	August 03, 2010, Tuesday, 02:03:36 EDT\n"
                + "End:	August 03, 2010, Tuesday, 02:03:36 EDT\n", maxits, convits,
                dampfact, n, -1, -1.0, -1.0, -1.0, iteration);
    }

    private AffinityPropagation(double dampfact, int convits, int maxits,
//			boolean usePreferenceAsIs, 
                                AbstMatrix W) {
        this.dampfact = dampfact;
        this.convits = convits;
        this.maxits = maxits;
//		this.usePreferenceAsIs = usePreferenceAsIs;
        this.W = W.myClone();
        this.n = W.rowCount();
        assert W.rowCount() == W.colCount();
    }

    public AffinityPropagation(AbstMatrix W) {
        this(0.9, 200, 2000,
//		    true, 
                W);
    }

    public static void main(String[] args) {
        AbstMatrix mat = MatrixFile.getFrom("M.txt");
        Debug.debug(mat.colCount());
        Debug.debug(mat.density());
        AffinityPropagation ap = new AffinityPropagation(mat);
        ap.clustering();
    }
}
