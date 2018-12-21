package test;

import forDebug.Debug;

import java.awt.*;
import java.util.ArrayList;

public class AGCTClustering_AP extends AGCTClustering_Algorithm {

    private static double LAMBDA = 0.5;

    /**
     * -P で指定される. preferences.
     */
    final int valueP;

    private double vB;
    private double sKK; // s(k,k)

    //	private Matrix P;
    private double[][] P;
    // computes s(i,j)

    //	private AbstMatrix exists;
    //	private AbstMatrix respons;
    //	private AbstMatrix avail;
    //	private AbstMatrix randa;

    private double[] maxAS1, maxAS2;
    private double[] maxR;
    private boolean[] secondFound;

    /**
     * -infでないindexをおさめる
     */
    //	private ArrayList<Integer>[] indices;
    private int[][] indices;
    private double[][] respons;
    private double[][] avail;
    private double[][] randa;
    /**
     * trace[i] = i列における、i行に対応するid.
     */
    private int[] trace;

    private int[] indexAS1;
    /**
     * exemplar[i] : gene i's parent
     */
    private int[] exemplar;
    private int[] clusterNumber;

    private boolean readyToSave;

    AGCTClustering_AP(Clustering myc, AGCT mya, int valp, double valb, boolean bef, int nMI) {
        super(myc, mya, -1, bef, nMI);
        setMyReferenceName("AP");
        setSeedingType(-1);
        // seedingType = -1;
        setFinalClusteringObjectiveFunction(-1);
        // finalClusteringObjectiveFunction = -1;
        valueP = valp;
        vB = valb;
        P = null;
        respons = null;
        avail = null;
        randa = null;
        exemplar = null;
        clusterNumber = null;
        readyToSave = false;

        maxAS1 = maxAS2 = maxR = null;
        indexAS1 = null;
        secondFound = null;
    }

    public String getAdditionalSavingData() {
        if (!readyToSave)
            return super.getAdditionalSavingData();

        String val = "Exemplars for clustering:\n";

        int i, j;
        boolean found = false;
        for (i = 0; i < getNumberOfClusters(); i++) {
            j = 0;
            while (!found) {
                if (clusterNumber[j] == i)
                    found = true;
                else
                    j++;
            }
            val += " - for cluster " + i + " : " + getMyClustering().getRightGeneSelectedGeneNumberToGeneNumber(exemplar[j], getBefore()) + "\n";
        }
        val += "\n";

        return val;
    }

    private static double damp(double newval, double lastval) {
        return ((LAMBDA * lastval) + ((1 - LAMBDA) * newval));
    }

    private void initAvail() {
        System.out.println("AGCTClustering_AP.initAvail()");
        int iN = getMyClustering().getRightNumberSelectedGenes(getBefore());

        defaultMaxAS();

        //		AbstMatrix avail=new GeneralMatrix(iN,iN);

        for (int i = 0; i < iN; i++)
            for (int k = 0; k < indices[i].length; k++) {
                double maj = (P[i][k]);
                updateMaxAS(maj, i, indices[i][k]);
            }

        //		randa=new GeneralMatrix(iN,iN);
        //		randa = new double[iN][];
        //
        //		for(int j=0;j<iN;j++){
        //			randa[j] = new double[indices[j].length];
        //			//			for (int j2 = 0; j2 < iN; j2++) {
        //			for(int j2:indices[j]){
        //				randa.set(j,j2,(0));
        //			}
        //		}
    }

    private void initRespons() {
        //		System.out.println("AGCTClustering_AP.initRespons()");
        //		int iN=getMyClustering().getRightNumberSelectedGenes(getBefore());
        //		respons=new GeneralMatrix(iN,iN);
        //		for(int k=0;k<iN;k++){
        //			for(int k2:indices[k]){
        //				respons.set(k,k2,(0));
        //			}
        //		}
    }

    private void initExemplars() {
        //		System.out.println("AGCTClustering_AP.initExemplars()");
        exemplar = new int[getMyClustering().getRightNumberSelectedGenes(getBefore())];
    }

    private void initMax() {
//		System.out.println("AGCTClustering_AP.initMax()");
        int iN = getMyClustering().getRightNumberSelectedGenes(getBefore());
        maxAS1 = new double[iN];
        maxAS2 = new double[iN];
        maxR = new double[iN];

        indexAS1 = new int[iN];
        secondFound = new boolean[iN];
        defaultMaxAS();
    }

    private void defaultMaxAS() {
        int i, iN = getMyClustering().getRightNumberSelectedGenes(getBefore());
        for (i = 0; i < iN; i++) {
            maxAS1[i] = maxAS2[i] = 0.0;
            indexAS1[i] = -1;
            secondFound[i] = false;
        }
    }

    private void initEverything() {
        System.out.println("AGCTClustering_AP.initEverything()");
        computeSKK();
        toPref();
        initMax();
        initAvail();
        initRespons();
        initExemplars();
    }

    //	public void plot(JAGCTVisualizationPane vv) {
    //		plotClusterCenters(vv);
    //		plotAffinities(vv);
    //	}

    //	public void plotWithAnnotations(JAGCTVisualizationPane vv) {
    //		plotClusterCentersWithAnnotations(vv);
    //		plotAffinities(vv);
    //	}

    public void plotAffinities(JAGCTVisualizationPane vv) {
        Color c;
        int i, iN = getMyClustering().getRightNumberSelectedGenes(getBefore()), ncl;
        View3D v = vv.visualizationPanel.myView3D;
        Pnt3D pe, p;
        if (getPlotCentersAvailable())
            if (vv.myReferenceName.equals(JAGCTVisualizationPane.M_P) || (vv.myReferenceName.equals(JAGCTVisualizationPane.P_P))) {
                for (i = 0; i < iN; i++) {
                    ncl = clusterNumber[i];
                    if (vv.structPlottedCluster(ncl)) {
                        pe = vv.getRightPoint(exemplar[i]);
                        p = vv.getRightPoint(i);
                        v.drawShadowLineVsRef(pe, p, JAGCTGraphicsPanel.minColor, getMyClustering().colorCluster[ncl], false);
                    }
                }
            }
    }

    public String toString() {
        String val = "Affinity Propagation ";
        if (getNumberOfClusters() <= 0)
            val += "(clusters not computed)\n";
        else
            val += "with " + getNumberOfClusters() + " clusters\n";
        val += "Initial values for Preference: " + valueP + ", for Bregman p-norm divergence with p=" + vB + "\n";
        val += "Points are " + Clustering.referenceStringCenters[getMyClustering().indexReferenceStringCenters] + ".";
        return val;
    }

    private double dist(Pnt a, Pnt b) {
        return -Distortion.Bregman(a, b, 4, vB);
    }

    /**
     * dist(vals[i],vals[j]) を昇順に並べたときの、id番目の要素を返す。
     * O( n^2 d log たくさん)
     *
     * @param vals
     * @param id
     * @return
     */
    double computeSKK_sub(Pnt[] vals, int id) {
        int n = vals.length;
        if (id < 0 || n * (n - 1) / 2 <= id)
            throw new RuntimeException();

        double left = Double.MAX_VALUE, right = -Double.MAX_VALUE;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double d = dist(vals[i], vals[j]);
                left = Math.min(left, d);
                right = Math.max(right, d);
            }
        }
        // left <= right
        for (int i = 0; i < 1000; i++) {
            double mid = (right + left) / 2;
            int cnt = 0;
            double max = Double.NEGATIVE_INFINITY;
            for (int j = 0; j < n; j++) {
                for (int k = j + 1; k < n; k++) {
                    double d = dist(vals[j], vals[k]);
                    if (d < mid) {
                        cnt++;
                        max = Math.max(max, d);
                    }
                }
            }
            if (cnt == id + 1)
                return max;
            if (cnt > id + 1) {
                right = mid;
            } else {
                left = mid;
            }
        }
        return right;
    }

    /**
     * 全点間距離を計算し、小さい方からpp % をとって,それを sKK としている。
     * vp == -1 : pp = 0
     * vp == -2 : pp = 50
     * otherwise : pp = vp
     */
    private void computeSKK() {
        System.out.println("AGCTClustering_AP.computeSKK()");
        int i, iN = getMyClustering().getRightNumberSelectedGenes(getBefore()), pp, lim;

        Pnt[] vals = new Pnt[iN];
        for (i = 0; i < iN; i++) {
            vals[i] = getMyClustering().getRightPoint(i, getBefore());
        }

        //		double[] simil=new double[iN*(iN-1)/2];
        //		DoDebug.debug("simil.length",simil.length);
        //		for(i=0;i<iN*(iN-1)/2;i++)
        //			simil[i]=0.0;

        //		k=0;
        //		for(i=0;i<iN-1;i++)
        //			for(j=i+1;j<iN;j++){
        //				simil[k]=-Distortion.Bregman(getMyClustering().getRightPoint(i,getBefore()),getMyClustering().getRightPoint(j,getBefore()),4,vB);
        //				k++;
        //			}

        //		QuickSort.quicksort(simil);

        if (valueP == -1)
            //			sKK=simil[0];
            sKK = computeSKK_sub(vals, 0);
        else {
            if (valueP == -2)
                pp = 50;
            else {
                if (valueP < 1)
                    Matrix.perror("invalid value for P");
                pp = valueP;
            }
            lim = (int) ((((long) iN * (iN - 1) / 2) * pp) / 100);
            //			DoDebug.debug("lim",lim);
            //			sKK=simil[lim];
            sKK = computeSKK_sub(vals, lim);
        }
        //		simil=null;

        // System.out.println("sKK = " + sKK + " (vP = " + vP + ")");
    }

    /**
     * Matrix Pをつくる。
     */
    private void toPref() {
        System.out.println("AGCTClustering_AP.toPref()");
        Matrix P = new Matrix("Preferences", getMyClustering().getRightNumberSelectedGenes(getBefore()), getMyClustering().getRightNumberSelectedGenes(getBefore()));
        P.toP(getMyAGCT().data.getMyDomain(), getMyClustering(), vB, sKK);
        //		exists = new GeneralMatrix(P.rowCount(),P.colCount());

        //		indices=new ArrayList[P.rowCount()];
        indices = new int[P.rowCount()][];
        this.P = new double[P.rowCount()][];
        respons = new double[P.rowCount()][];
        avail = new double[P.rowCount()][];
        randa = new double[P.rowCount()][];
        trace = new int[P.rowCount()];

        for (int i = 0; i < P.rowCount(); i++) {
            //			indices[i]=new ArrayList<Integer>();
            ArrayList<Integer> list = new ArrayList<Integer>();
            for (int j = 0; j < P.colCount(); j++) {
                if (P.get(i, j) != 0 || i == j) {
                    if (i == j) {
                        trace[i] = list.size();
                    }
                    list.add(j);
                }
            }

            indices[i] = new int[list.size()];
            this.P[i] = new double[list.size()];
            respons[i] = new double[list.size()];
            avail[i] = new double[list.size()];
            randa[i] = new double[list.size()];
            for (int j = 0; j < list.size(); j++) {
                indices[i][j] = list.get(j);
                this.P[i][j] = P.get(i, indices[i][j]);
                //必要ないかも
                if (this.P[i][j] == 0)
                    this.P[i][j] = Double.NEGATIVE_INFINITY;
            }
        }

        for (int i = 0; i < this.P.length; i++) {
            for (int j = 0; j < this.P[i].length; j++) {
                assert this.P[i][j] != 0;
            }
            //			DoDebug.debug("this.P["+i+"].length",this.P[i].length);
        }

        P = null;
    }
    /**
     * matrixに入っている値を実際の値に変換する
     * @param d
     * @return
     */
    //	private double tod(double d){
    //		if(d==0)
    //			return Double.NEGATIVE_INFINITY;
    //		if(d==Double.NEGATIVE_INFINITY)
    //			return 0;
    //		return d;
    //	}
    /**
     * 実際の値をmatrixに入れる値に変換する
     * @param d
     * @return
     */
    //	private double tom(double d){
    //		if(d==Double.NEGATIVE_INFINITY)
    //			return 0;
    //		if(d==0)
    //			return Double.NEGATIVE_INFINITY;
    //		return d;
    //	}

    /**
     * maxRを更新していく。exists(i,k)==0なる部分はみない。exists(i,k)==1なるインデックスだけを持っておけば、速い。
     *
     * @param it
     */
    private void updateRespons(int it) {
        System.out.println("AGCTClustering_AP.updateRespons()");
        int iN = getMyClustering().getRightNumberSelectedGenes(getBefore());
        double dum;

        for (int i = 0; i < iN; i++)
            maxR[i] = 0.0;

        for (int i = 0; i < iN; i++) {
            for (int k = 0; k < indices[i].length; k++) {

                if (indexAS1[i] == indices[i][k])
                    dum = maxAS2[i];
                else
                    dum = maxAS1[i];

                if (it > 0)
                    respons[i][k] = AGCTClustering_AP.damp(respons[i][k], P[i][k] - dum);
                else
                    respons[i][k] = P[i][k] - dum;

                maxR[indices[i][k]] += Math.max(0.0, respons[i][k]);
            }
        }
    }

    /**
     * １、２番目に大きいrの値を更新する　O(1)
     *
     * @param maj
     * @param i
     * @param j
     */
    private void updateMaxAS(double maj, int i, int j) {
        if ((indexAS1[i] == -1) || (maj > maxAS1[i]) || ((maj == maxAS1[i]) && (indexAS1[i] != j))) {
            if ((indexAS1[i] != -1) && ((maj > maxAS2[i]) || (!secondFound[i]))) {
                maxAS2[i] = maxAS1[i];
                secondFound[i] = true;
            }
            maxAS1[i] = maj;
            indexAS1[i] = j;
        } else if ((indexAS1[i] != -1) && ((maj > maxAS2[i]) || (!secondFound[i]))) {
            maxAS2[i] = maj;
            secondFound[i] = true;
        }
    }

    private void updateAvail(int it) {
        System.out.println("AGCTClustering_AP.updateAvail()");
        int i, iN = getMyClustering().getRightNumberSelectedGenes(getBefore());
        double dum, maj;
        defaultMaxAS();

        for (i = 0; i < iN; i++) {
            for (int k = 0; k < indices[i].length; k++) {

                if (indices[i][k] != i)
                    dum = (respons[indices[i][k]][trace[indices[i][k]]]) + maxR[indices[i][k]] - Math.max(0.0, (respons[i][k]))
                            - Math.max(0.0, (respons[indices[i][k]][trace[indices[i][k]]]));
                else
                    dum = maxR[indices[i][k]] - Math.max(0.0, (respons[indices[i][k]][trace[indices[i][k]]]));

                if (indices[i][k] != i)
                    dum = Math.min(0.0, dum);

                if (it > 0)
                    avail[i][k] = ((AGCTClustering_AP.damp((avail[i][k]), dum)));
                else
                    avail[i][k] = ((dum));

                maj = (avail[i][k]) + (P[i][k]);
                updateMaxAS(maj, i, indices[i][k]);
            }
        }

		/*
         * for (i=0;i<iN;i++) System.out.print("  MaxAS1[" + i + "] = "
		 * +maxAS1[i]); System.out.println(""); for (i=0;i<iN;i++)
		 * System.out.print("  indexAS1[" + i + "] = " + indexAS1[i]);
		 * System.out.println(""); for (i=0;i<iN;i++)
		 * System.out.print("  MaxAS2[" + i + "] = " +maxAS2[i]);
		 * System.out.println(""); for (i=0;i<iN;i++)
		 * System.out.print("  secondFound[" + i + "] = " + secondFound[i]);
		 * System.out.println("");
		 */

    }

    /**
     * O(-infでない要素数)
     *
     * @param it
     */
    private void updateRanda(int it) {
        System.out.println("AGCTClustering_AP.updateRanda()");
        int i, iN = getMyClustering().getRightNumberSelectedGenes(getBefore());
        for (i = 0; i < iN; i++)
            for (int k = 0; k < indices[i].length; k++)
                randa[i][k] = (((avail[i][k]) + (respons[i][k])));
    }

    private void updateExemplars() {
        System.out.println("AGCTClustering_AP.updateExemplars()");
        int i, iN = getMyClustering().getRightNumberSelectedGenes(getBefore());
        double maxar;

        for (i = 0; i < iN; i++) {
            maxar = Double.NEGATIVE_INFINITY;
            for (int k = 0; k < indices[i].length; k++)
                if (maxar < randa[i][k]) {
                    maxar = randa[i][k];
                    exemplar[i] = indices[i][k];
                }
        }

		/*
         * for (i=0;i<iN;i++) System.out.print(" [" + i + " --> " + exemplar[i]
		 * + "] "); System.out.println("");
		 */
    }

    private void computeClusterNumber() {
        int i, ref, replace, ni = -1, iN = getMyClustering().getRightNumberSelectedGenes(getBefore());
        clusterNumber = new int[iN];
        boolean larger = true, foundcluster = false;

        for (i = 0; i < iN; i++)
            clusterNumber[i] = exemplar[i];

        ref = replace = 0;
        setNumberOfClusters(0);
        // numberOfClusters = 0;

        for (i = 0; i < iN; i++)
            if ((i == 0) || (clusterNumber[i] < replace))
                replace = clusterNumber[i];

        do {
            foundcluster = false;
            ni = iN;
            for (i = 0; i < iN; i++) {
                if ((!foundcluster) && (clusterNumber[i] == replace)) {
                    foundcluster = true;
                    setNumberOfClusters(getNumberOfClusters() + 1);
                    // numberOfClusters++;
                }
                if ((clusterNumber[i] > replace) && (clusterNumber[i] < ni))
                    ni = clusterNumber[i];

                if (clusterNumber[i] == replace)
                    clusterNumber[i] = ref;
            }
            if ((ni > replace) && (ni < iN)) {
                replace = ni;
                ref++;
            } else
                larger = false;
        } while (larger == true);
        Debug.debug("AGCTClustering_AP : computeClusterNumber(),clusterNumber", clusterNumber);


        getMyClustering().nclusters = getNumberOfClusters();


        captain = new int[getMyClustering().nclusters];
        for (int j = 0; j < iN; j++) {
            if (exemplar[j] == j) {
                assert captain[clusterNumber[j]] == 0;
                captain[clusterNumber[j]] = j;
            }
        }
//		JOptionPane.showConfirmDialog(null, Arrays.toString(captain));
    }

    /**
     * captain[i] = captain of cluster i
     */
    int[] captain;

    private void toMemberships() {
        initSoftMemberships();
        int i, j;
        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++) {
            for (j = 0; j < getNumberOfClusters(); j++) {
                getSoft_memberships().set(j, i, 0.0);
                if (clusterNumber[i] == j)
                    getSoft_memberships().set(j, i, 1.0);
            }
        }
        setSoft_memberships_available(true);
        //		soft_memberships_available = true;
    }

    public void saveSpace() {
        P = null;
        respons = avail = randa = null;
        maxAS1 = maxAS2 = maxR = null;
        secondFound = null;
        indexAS1 = null;
    }

    private void AP() {
        System.out.println("AGCTClustering_AP.AP()");
        boolean stop = false;
        int i, niter = 0, iN = getMyClustering().getRightNumberSelectedGenes(getBefore()), nid = 0, ham;
        int[] lastexemplar = new int[iN];

        AGCTCounter ccc;

        if (getnMaxIterations() > 0)
            ccc = new AGCTCounter("Clustering AP", getnMaxIterations());
        else
            ccc = new AGCTCounter("Clustering AP", Affinity_Propagation_Id_Max);

        Debug.debug("getnMaxIterations()", getnMaxIterations());
        Debug.debug("Affinity_Propagation_Id_Max", Affinity_Propagation_Id_Max);
        System.gc();
        Debug.debug("finish gc");
        do {

            updateRespons(niter);
            updateAvail(niter);
            updateRanda(niter);
            updateExemplars();
            // ccc.increment();

            if (niter > 0) {
                ham = Matrix.Hamming(lastexemplar, exemplar);
                if (ham == 0)
                    nid++;
                else
                    nid = 0;

                Debug.debug("Hamming", ham);
            }

            if (getnMaxIterations() > 0)
                ccc.increment();
            else
                ccc.setBound(nid);

            for (i = 0; i < iN; i++)
                lastexemplar[i] = exemplar[i];

            niter++;
            if (getnMaxIterations() > 0)
                if (niter >= getnMaxIterations())
                    stop = true;

            if (nid >= Affinity_Propagation_Id_Max)
                stop = true;

            Debug.debug("niter", niter);
            Debug.debug("nid", nid);

        } while (!stop);
        ccc.end();

        setPlotCentersAvailable(true);
        //		plotCentersAvailable = true;
    }

    public void toClustering() {
        System.out.println("AGCTClustering_AP.toClustering()");
        //		dimRef = getDimension(getBefore());

        initEverything();
        AP();


        computeClusterNumber();
        getMyClustering().toClusterColors();

        toMemberships();
        toClusterChoices();
        setSoft_clustering_computed(true);
        // soft_clustering_computed = true;
        readyToSave = true;

        toHardMemberships();

        toSelectedGeneToHardCluster();
        setHard_clustering_computed(true);
        setOk_for_statistics(true);
        ControlProcess.put("hardClusteringProcessed", true);
        ControlProcess.put("ClusteringProcessed_AP", true);
    }
}