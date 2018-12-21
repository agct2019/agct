package test;

import java.util.Vector;


public class AGCTClustering_NM extends AGCTClustering_Algorithm {

    int n_Variables;
    // number of columns of the data matrix (typically = number of description
    // variables)

    int n_Method;
    // method to cast matrix V = WH in positive form
    // = 0 : method of Kim and Tidor, NAR 2008 (see below)
    // = 1 : shift my the min

    int t_Update;
    // type of update
    // = 0 : Euclidean update
    // = 1 : KL update

    static int DEFAULT_N_METHOD = 0;
    static int DEFAULT_T_UPDATE = 1;
    static int DEFAULT_N_CLUSTERS = 2;

    Matrix V_NM, W_NM, H_NM;

    public static String name_Method(AGCTClustering_NM a) {
        int i = a.n_Method;
        String ret = "Casting entries into positive values following ";
        if (i == 0)
            ret += "Kim & Tidor's duplicate method (NAR 2008)";
        else if (i == 1)
            ret += "Min shift for all values";
        ret += ".";
        return ret;
    }

    public static String name_Update(AGCTClustering_NM a) {
        int i = a.t_Update;
        String ret = "Update rule minimizes distortion " + Distortion.DIVERGENCES_NAME[i] + ".";
        return ret;
    }

    AGCTClustering_NM(Clustering myc, AGCT mya, int nc, boolean bef, int nMI, int nMet, int tUpd) {
        super(myc, mya, nc, bef, nMI);
        // myReferenceName = "NM";
        setMyReferenceName("NM");
        setFinalClusteringObjectiveFunction(-1);
        // finalClusteringObjectiveFunction = -1;
        V_NM = W_NM = H_NM = null;

        if (nMet == -1)
            n_Method = DEFAULT_N_METHOD;
        else if (nMet < 0)
            Matrix.perror("AGCTClustering_NM :: Bad value for the positive defining method");
        else
            n_Method = nMet;

        if (tUpd == -1)
            t_Update = DEFAULT_T_UPDATE;
        else if ((tUpd != 0) && (tUpd != 1))
            Matrix.perror("AGCTClustering_NM :: Bad value for the update rule method");
        else
            t_Update = tUpd;

        n_Variables = getMyClustering().getRightPoint(0, getBefore()).coordinates.length;
        // System.out.println("n_Variables = " + n_Variables + ", nMethod = " +
        // nMet);
    }

    public void safeCheck() {
        // Rule taken from
        // "Subsystem Identification Through Dimensionality Reduction of Large-Scale Gene Expression Data"
        // Philip M. Kim and Bruce Tidor, Nucleics Acid Research 2008

        int mm = getMyClustering().getRightNumberSelectedGenes(getBefore());
        int nn = n_Variables;
        int kk = getNumberOfClusters();
        boolean passTest = false;

        if ((n_Method == 0) && (kk * (mm + 2 * nn) < 2 * nn * mm))
            passTest = true;
        else if ((n_Method != 0) && (kk * (mm + nn) < nn * mm))
            passTest = true;

        if (!passTest) {
            if (AGCT.DoDebug)
                System.out.println("Too many clusters for the dimension :: switching to full Manifold coordinates");
            getMyClustering().indexReferenceStringCenters = 0;
        }
    }

    public void initNMF() {
        safeCheck();
        init_V_NM();
        init_W_NM();
        init_H_NM();
        finalInit();
    }

    public void init_V_NM() {
        int i, j, dimX, dimY;
        double dum, min = -1.0;
        if (n_Method == 0)
            dimX = 2 * n_Variables;
        else
            dimX = n_Variables;
        dimY = getMyClustering().getRightNumberSelectedGenes(getBefore());

        V_NM = new Matrix("Data Matrix (NMF)", dimX, dimY);

        if (n_Method == 0) {
            for (i = 0; i < n_Variables; i++) {
                for (j = 0; j < getMyClustering().getRightNumberSelectedGenes(getBefore()); j++) {
                    dum = getMyClustering().getRightPoint(j, getBefore()).coordinates[i];
                    if (dum < 0.0) {
                        V_NM.set(i, j, 0.0);
                        V_NM.set(i + n_Variables, j, -dum);
                    } else {
                        V_NM.set(i, j, dum);
                        V_NM.set(i + n_Variables, j, 0.0);
                    }
                }
            }
        } else if (n_Method == 1) {
            for (i = 0; i < n_Variables; i++) {
                for (j = 0; j < getMyClustering().getRightNumberSelectedGenes(getBefore()); j++) {
                    dum = getMyClustering().getRightPoint(j, getBefore()).coordinates[i];
                    if (((i == 0) && (j == 0)) || (dum < min))
                        min = dum;
                }
            }
            for (i = 0; i < n_Variables; i++) {
                for (j = 0; j < getMyClustering().getRightNumberSelectedGenes(getBefore()); j++) {
                    dum = getMyClustering().getRightPoint(j, getBefore()).coordinates[i];
                    V_NM.set(i, j, dum - min);
                }
            }
        }
    }

    public void init_W_NM() {
        int i, j, dimX, dimY;
        double dum;
        dimX = V_NM.rowCount();
        dimY = getNumberOfClusters();

        W_NM = new Matrix("Meta Variables (NMF)", dimX, dimY);

        for (i = 0; i < dimX; i++) {
            for (j = 0; j < dimY; j++) {
                W_NM.set(i, j, AGCT.RandomGenerator.nextDouble());
                if (W_NM.get(i, j) == 0.0)
                    W_NM.set(i, j, 0.5);

                if (W_NM.get(i, j) < 0.0)
                    W_NM.set(i, j, -W_NM.get(i, j));
            }
        }
    }

    public void init_H_NM() {
        int i, j, dimX, dimY;
        double dum;
        dimX = getNumberOfClusters();
        dimY = getMyClustering().getRightNumberSelectedGenes(getBefore());

        H_NM = new Matrix("Meta Clustering (NMF)", dimX, dimY);

        for (i = 0; i < dimX; i++) {
            for (j = 0; j < dimY; j++) {
                H_NM.set(i, j, AGCT.RandomGenerator.nextDouble());
                if (H_NM.get(i, j) == 0.0)
                    H_NM.set(i, j, 0.5);

                if (H_NM.get(i, j) < 0.0)
                    H_NM.set(i, j, -H_NM.get(i, j));
            }
        }
    }

    public void finalInit() {
        if (n_Method == 0)
            n_Variables = 2 * n_Variables;
    }

//	public void plot(JAGCTVisualizationPane vv) {
//		plotClusterCenters(vv);
//	}

//	public void plotWithAnnotations(JAGCTVisualizationPane vv) {
//		plotClusterCentersWithAnnotations(vv);
//	}

    public String toString() {
        String val = "Non Negative Matrix Factorization with " + getNumberOfClusters() + " clusters, over " + n_Variables + " description variables.\n";
        val += "Points are " + Clustering.referenceStringCenters[getMyClustering().indexReferenceStringCenters] + ".\n";
        val += AGCTClustering_NM.name_Method(this) + ".\n" + AGCTClustering_NM.name_Update(this) + ".\n";
        return val;
    }

    public double update_L22(Matrix nextW, Matrix Wt, Matrix WtW, Matrix nextH, Matrix Ht, Matrix HHt, Matrix WtWH, Matrix WtV, Matrix VHt, Matrix WHHt) {
        int i, j;
        double dum;

        for (i = 0; i < H_NM.rowCount(); i++)
            for (j = 0; j < H_NM.colCount(); j++) {
                dum = H_NM.get(i, j) * WtV.get(i, j) / WtWH.get(i, j);
                nextH.set(i, j, dum);
            }

        H_NM.copyOnThis(nextH);
        Ht.transpose(H_NM);
        HHt.dot(H_NM, Ht);
        VHt.dot(V_NM, Ht);
        WHHt.dot(W_NM, HHt);

        for (i = 0; i < W_NM.rowCount(); i++)
            for (j = 0; j < W_NM.colCount(); j++) {
                dum = W_NM.get(i, j) * VHt.get(i, j) / WHHt.get(i, j);
                nextW.set(i, j, dum);
            }

        W_NM.copyOnThis(nextW);
        Wt.transpose(W_NM);
        WtW.dot(Wt, W_NM);
        WtWH.dot(WtW, H_NM);
        WtV.dot(Wt, V_NM);

        return Matrix.BregmanError(V_NM, W_NM, H_NM, 0, -1);
    }

    public double update_KL(Matrix nextW, Matrix nextH, Matrix WH) {
        int i, j, k;
        double dum, numer, denom;
        for (i = 0; i < H_NM.rowCount(); i++)
            for (j = 0; j < H_NM.colCount(); j++) {
                numer = 0.0;
                for (k = 0; k < V_NM.rowCount(); k++)
                    numer += (W_NM.get(k, i) * V_NM.get(k, j) / WH.get(k, j));
                denom = 0.0;
                for (k = 0; k < W_NM.rowCount(); k++)
                    denom += W_NM.get(k, i);
                nextH.set(i, j, H_NM.get(i, j) * numer / denom);
            }

        H_NM.copyOnThis(nextH);
        WH.dot(W_NM, H_NM);

        for (i = 0; i < W_NM.rowCount(); i++)
            for (j = 0; j < W_NM.colCount(); j++) {
                numer = 0.0;
                for (k = 0; k < V_NM.colCount(); k++)
                    numer += (H_NM.get(j, k) * V_NM.get(i, k) / WH.get(i, k));
                denom = 0.0;
                for (k = 0; k < H_NM.colCount(); k++)
                    denom += H_NM.get(j, k);
                nextW.set(i, j, W_NM.get(i, j) * numer / denom);
            }

        W_NM.copyOnThis(nextW);
        WH.dot(W_NM, H_NM);

        return Matrix.BregmanError(V_NM, W_NM, H_NM, 1, -1);
    }

    public void NMF() {
        Matrix Wt, Ht, WtW, WtWH, WtV, VHt, HHt, WHHt, WH;

        Wt = new Matrix("Wt", W_NM.colCount(), W_NM.rowCount());
        Wt.transpose(W_NM);
        Ht = new Matrix("Ht", H_NM.colCount(), H_NM.rowCount());
        Ht.transpose(H_NM);
        HHt = new Matrix("HHt", H_NM.rowCount(), Ht.colCount());
        HHt.dot(H_NM, Ht);
        WtW = new Matrix("WtW", Wt.rowCount(), W_NM.colCount());
        WtW.dot(Wt, W_NM);
        WtWH = new Matrix("WtWH", WtW.rowCount(), H_NM.colCount());
        WtWH.dot(WtW, H_NM);
        WtV = new Matrix("WtV", Wt.rowCount(), V_NM.colCount());
        WtV.dot(Wt, V_NM);
        VHt = new Matrix("VHt", V_NM.rowCount(), Ht.colCount());
        VHt.dot(V_NM, Ht);
        WHHt = new Matrix("WHHt", W_NM.rowCount(), HHt.colCount());
        WHHt.dot(W_NM, HHt);
        WH = new Matrix("WH", W_NM.rowCount(), H_NM.colCount());
        WH.dot(W_NM, H_NM);

        Matrix nextW, nextH;

        int i, niter = 0;
        double delta;
        double dprev = -1, dcur, ratio;
        boolean stop = false;
        AGCTCounter ccc;

        if (getnMaxIterations() > 0)
            ccc = new AGCTCounter("Clustering NM", getnMaxIterations());
        else
            ccc = new AGCTCounter("Clustering NM", 0);

        nextW = new Matrix("Copy W", W_NM.rowCount(), W_NM.colCount());
        nextH = new Matrix("Copy H", H_NM.rowCount(), H_NM.colCount());

        do {
            nextW.copyOnThis(W_NM);
            nextH.copyOnThis(H_NM);

            if (t_Update == 0)
                dcur = update_L22(nextW, Wt, WtW, nextH, Ht, HHt, WtWH, WtV, VHt, WHHt);
            else
                dcur = update_KL(nextW, nextH, WH);

            setFinalClusteringObjectiveFunction(dcur);
            // finalClusteringObjectiveFunction = dcur;

            // System.out.println("NMF distortions :: " + dprev + " --> " +
            // dcur);

            if (niter > 0) {
                if (getnMaxIterations() > 0)
                    ccc.increment();
                else
                    ccc.setPercent((int) (100.0 * AGCTClustering_Algorithm.EPSILON_KM / (dprev - dcur)));

                if (dprev < dcur)
                    Matrix.perror("Distortion increase :: " + dprev + " --> " + dcur);

                ratio = (dprev - dcur) / dprev;
                // getMyAGCT().myInformationFrame.appendText(niter + " : " +
                // ratio + "\n");
                if (ratio < AGCTClustering_Algorithm.EPSILON_KM)
                    stop = true;
            }
            dprev = dcur;
            niter++;

            if (getnMaxIterations() > 0)
                if (niter >= getnMaxIterations())
                    stop = true;

        } while (!stop);
        JInformationFrame.getInstance().setText("ok.");
        ccc.end();

        toCentersAndMemberships();
        nextW = nextH = Wt = Ht = HHt = WtW = WtWH = WtV = VHt = WHHt = WH = null;
    }

    public void toCentersAndMemberships() {
        setHard_memberships_centers(new Vector<Pnt>());
        // hard_memberships_centers = new Vector();

        int[] affect = new int[getMyClustering().getRightNumberSelectedGenes(getBefore())];
        int i, j, k, imax;
        double tot, vmax = -1.0, vcur;
        Pnt p;

        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++) {
            imax = -1;
            for (j = 0; j < getNumberOfClusters(); j++) {
                vcur = H_NM.get(j, i);
                if (vcur < 0.0)
                    Matrix.perror("AGCTClustering_NM.class :: negative value");
                if ((j == 0) || (vcur > vmax)) {
                    imax = j;
                    vmax = vcur;
                }
            }
            affect[i] = imax;
        }

        for (i = 0; i < getNumberOfClusters(); i++)
            if (containsIndex(affect, i)) {
                tot = 0.0;
                p = new Pnt(getDimension(getBefore()));
                for (j = 0; j < affect.length; j++)
                    if (affect[j] == i) {
                        tot += 1.0;
                        for (k = 0; k < p.coordinates.length; k++)
                            p.coordinates[k] += getMyClustering().getRightPoint(j, getBefore()).coordinates[k];
                    }
                for (k = 0; k < p.coordinates.length; k++)
                    p.coordinates[k] /= tot;
                getHard_memberships_centers().add(p);
            }

        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++)
            for (j = 0; j < getNumberOfClusters(); j++)
                if (affect[i] == j)
                    getSoft_memberships().set(j, i, 1.0);
                else
                    getSoft_memberships().set(j, i, 0.0);
        setSoft_memberships_available(true);
//		soft_memberships_available = true;
        affect = null;
    }

    public void toClustering() {
        initNMF();
        NMF();
        toClusterChoices();
        toSelectedGeneToHardCluster();
        setPlotCentersAvailable(true);
//		plotCentersAvailable = true;
        setHard_clustering_computed(true);
//		hard_clustering_computed = true;
        setSoft_clustering_computed(false);
        // soft_clustering_computed = false;
        setOk_for_statistics(true);
//		ok_for_statistics = true;
        ControlProcess.put("hardClusteringProcessed", true);
        ControlProcess.put("ClusteringProcessed_NM", true);
    }
}