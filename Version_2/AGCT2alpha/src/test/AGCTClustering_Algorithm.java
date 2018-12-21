package test;

import filtering.JLeftIndicator;
import forDebug.Debug;
import gene.Gene;
import matrix.MatrixConfig;

import java.awt.*;
import java.util.Vector;

public class AGCTClustering_Algorithm implements Debuggable {

    protected static final double EPSILON_KM = 10E-10;//KM,NM only.
    protected static final double EPSILON_EM = 10E-20;//EM only.
    public static final String REFERENCE_NAME_KM = "KM";
    public static final String REFERENCE_NAME_EM = "EM";
    public static final String REFERENCE_NAME_CP = "CP";
    public static final String REFERENCE_NAME_HC = "HC";
    public static final String REFERENCE_NAME_AP = "AP";
    public static final String REFERENCE_NAME_NM = "NM";
    public static final String ALL_CLUSTERS = "All";
    static final String CLUSTER_NAME = "C";
    public String myReferenceName;
    // One of KM, CP, EM, HC, AP, NM (see above)

    public String getMyReferenceName() {
        return myReferenceName;
    }

    protected void setMyReferenceName(String myReferenceName) {
        this.myReferenceName = myReferenceName;
    }

    public Clustering myClustering;//child only

    protected Clustering getMyClustering() {
        return myClustering;
    }

    private AGCT myAGCT;//child only

    protected AGCT getMyAGCT() {
        return myAGCT;
    }

    public int numberOfClusters;

    public int getNumberOfClusters() {
        return numberOfClusters;
    }

    protected void setNumberOfClusters(int numberOfClusters) {//child only.
        this.numberOfClusters = numberOfClusters;
    }

    private final int nMaxIterations;//AP,EM,KM,NM only.

    protected int getnMaxIterations() {
        return nMaxIterations;
    }

    private int seedingType;//AP,CP,EM,HC,KM,NM only.
    // -1 = N/A
    // 0 = Forgy
    // 1 = Arthur - Vassilvitskii
    // 2 = Bregman

    protected int getSeedingType() {//EM,KM only.
        return seedingType;
    }

    protected void setSeedingType(int seedingType) {
        this.seedingType = seedingType;
    }

    protected static final String SEEDING_NAME[] = {"Forgy", "Arthur-Vassilvitskii", "Bregman (general AV)"};//EM,KM only.

    private double finalClusteringObjectiveFunction;

    public double getFinalClusteringObjectiveFunction() {
        Debug.debug("#### AGCTClustering_Algorithm.getFinalClusteringObjectiveFunction() ####");
        return finalClusteringObjectiveFunction;
    }

    protected void setFinalClusteringObjectiveFunction(double finalClusteringObjectiveFunction) {//AP,EM,KM,NM only.
        this.finalClusteringObjectiveFunction = finalClusteringObjectiveFunction;
    }

    private double[][] referencesByCluster;//AGCTFileWriter only. other classes don't rewrite.

    public double[][] getReferencesByCluster() {
        return referencesByCluster;
    }

    /**
     * soft_memberships.get(clusterId,geneId) = probability of the gene is belong to the cluster.
     */
    public Matrix soft_memberships;//child , JClusterFrame,Matrix

    /**
     * soft_memberships.get(clusterId,geneId) = probability of the gene is belong to the cluster.
     */
    public Matrix getSoft_memberships() {
        return soft_memberships;
    }

    void setSoft_memberships(Matrix softMemberships) {//CP only.
        soft_memberships = softMemberships;
    }

    public Vector<Pnt> hard_memberships_centers;//EM,KM,NM only.
    // hard_memberships_centers = used in clustering algorithms

    protected Vector<Pnt> getHard_memberships_centers() {
        return hard_memberships_centers;
    }

    protected void setHard_memberships_centers(Vector<Pnt> hardMembershipsCenters) {//KM,NM only.
        hard_memberships_centers = hardMembershipsCenters;
    }

    private Vector<Pnt> hard_memberships_centers_manifold_total;
    public Vector<Pnt3D> hard_memberships_centers_manifold_Pnt3D;
    // hard_memberships_centers_manifold_total = computed from soft_memberships
    // for display
    // hard_memberships_centers_manifold_Pnt3D = keeps the current 3D view of
    // cluster centers (faster than making computation each time !)
    private Vector<Pnt> hard_memberships_centers_pca_total;
    public Vector<Pnt3D> hard_memberships_centers_pca_Pnt3D;
    // hard_memberships_centers_pca_total = computed from soft_memberships for
    // display
    // hard_memberships_centers_pca_Pnt3D = keeps the current 3D view of cluster
    // centers (faster than making computation each time !)
    private boolean soft_clustering_computed;//child only.

    protected void setSoft_clustering_computed(boolean soft_clustering_computed) {
        this.soft_clustering_computed = soft_clustering_computed;
    }

    private boolean hard_clustering_computed;//child only.

    protected void setHard_clustering_computed(boolean hard_clustering_computed) {
        this.hard_clustering_computed = hard_clustering_computed;
    }

    private boolean soft_memberships_available;//child only.

    protected void setSoft_memberships_available(boolean soft_memberships_available) {
        this.soft_memberships_available = soft_memberships_available;
    }

    private boolean ok_for_statistics;

    boolean getOk_for_statistics() {// JVisualizationPane only.
        return ok_for_statistics;
    }

    protected void setOk_for_statistics(boolean ok_for_statistics) {//child only.
        this.ok_for_statistics = ok_for_statistics;
    }

    public boolean plotCentersAvailable;//child only.

    boolean getPlotCentersAvailable() {//AP only.
        return plotCentersAvailable;
    }

    protected void setPlotCentersAvailable(boolean plotCentersAvailable) {
        this.plotCentersAvailable = plotCentersAvailable;
    }

    final private boolean before;//child , Matrix.
    // true iff we make clutering on data before prototype selection

    protected boolean getBefore() {
        return before;
    }

    private String[] clusterChoices;// AGCT,JAGCTVisualizationPane.

    String[] getClusterChoices() {
        return clusterChoices;
    }

    private int[] selectedGeneToHardCluster;//JClusteringProfileFrame.
    // gives for each gene the most probable cluster to which it belongs

    public int[] getSelectedGeneToHardCluster() {
        return selectedGeneToHardCluster;
    }

    //	public int[] clusterSizesPrototypes;// my only.
    public int[] clusterSizesWhole;
    // Sizes of clusters as prototypes or for all points

    public AGCTClustering_Algorithm(Clustering myc, AGCT mya, int _numberOfClusters, boolean bef, int nMI) {
        myReferenceName = "";

        numberOfClusters = _numberOfClusters;
        nMaxIterations = nMI;
        myClustering = myc;
        myAGCT = mya;

        soft_memberships = null;
        hard_memberships_centers = null;
        hard_memberships_centers_manifold_total = null;
        hard_memberships_centers_manifold_Pnt3D = null;

        hard_memberships_centers_pca_total = null;
        hard_memberships_centers_pca_Pnt3D = null;

        referencesByCluster = null;

        soft_clustering_computed = false;
        hard_clustering_computed = false;
        plotCentersAvailable = false;
        soft_memberships_available = false;
        ok_for_statistics = false;

        before = bef;
    }

    public void toClusteringLite(Vector allPoints) {
        Matrix.perror("AGCTClustering_Algorithm.class :: no such algorithm");
    }

    public void saveSpace() {
    }

    public String getAdditionalSavingData() {
        return "None";
    }

    final public void plotClusterCentersWithAnnotations(JAGCTVisualizationPane vv) {
        int i;
        View3D v = vv.visualizationPanel.myView3D;
        String sz;
        if (plotCentersAvailable) {

            if (vv.myReferenceName.equals(JAGCTVisualizationPane.M_P)) {
                softMembershipsToHardCenters_Manifold_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
                for (i = 0; i < numberOfClusters; i++) {
                    if (vv.structPlottedCluster(i)) {
                        v.g.setColor(myClustering.colorCluster[i]);
                        if (Prototype.No_Reduction) {
                            sz = "" + JLeftIndicator.getInstance().clusterMemberCount(i);
                        } else {
                            sz = "" + JLeftIndicator.getInstance().clusterMemberCount(i) + "|" + clusterSizesWhole[i];
                        }
                        v.drawCenterVsRefWithAnnotations((Pnt3D) hard_memberships_centers_manifold_Pnt3D.elementAt(i), i, sz);
                    }
                }
            } else if (vv.myReferenceName.equals(JAGCTVisualizationPane.P_P)) {
                softMembershipsToHardCenters_Pca_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
                for (i = 0; i < numberOfClusters; i++) {
                    if (vv.structPlottedCluster(i)) {
                        v.g.setColor(myClustering.colorCluster[i]);
                        if (Prototype.No_Reduction) {
                            sz = "" + JLeftIndicator.getInstance().clusterMemberCount(i);
                        } else {
                            sz = "" + JLeftIndicator.getInstance().clusterMemberCount(i) + "|" + clusterSizesWhole[i];
                        }
                        v.drawCenterVsRefWithAnnotations((Pnt3D) hard_memberships_centers_pca_Pnt3D.elementAt(i), i, sz);
                    }
                }
            }
        }
    }

    final public int getDimension(boolean before) {
        return myClustering.getRightPoint(0, before).coordinates.length;
    }

    final public void initSoftMemberships() {
        soft_memberships = new Matrix("Soft_Memberships_" + myReferenceName, numberOfClusters, myClustering.getRightNumberSelectedGenes(before));
        // Nrows
        // =
        // number
        // of
        // clusters
        int i, j;

        for (i = 0; i < myClustering.getRightNumberSelectedGenes(before); i++) {
            for (j = 0; j < numberOfClusters; j++) {
                soft_memberships.set(j, i, 0.0);
            }
        }
    }

    /**
     * ランダムに（重複は無いように）numberOfClusters個の遺伝子をseedとして選ぶ。
     */
    public void initHardMembershipForgy() { // Same
        int i, j, k;
        boolean deja;
        int[] indexes_to_pick = new int[numberOfClusters];

        for (i = 0; i < numberOfClusters; i++) {
            indexes_to_pick[i] = 0;
        }

        for (i = 0; i < numberOfClusters; i++) {
            do {
                j = AGCT.RandomGenerator.nextInt(myClustering.getRightNumberSelectedGenes(before));
                k = 0;
                deja = false;
                do {
                    if (indexes_to_pick[k] == j) {
                        deja = true;
                    }
                    k++;
                } while ((k < i) && (deja == false));
            } while (deja == true);
            indexes_to_pick[i] = j;
        }

        for (i = 0; i < numberOfClusters; i++) {
            hard_memberships_centers.add(new Pnt(myClustering.getRightPoint(indexes_to_pick[i], before))); // Correction
        }        // !
    }

    public void softMembershipsToHardCenters_Manifold_Total() {
        int i, j, k, dim = AGCT.Number_Of_Manifold_Components;
        hard_memberships_centers_manifold_total = new Vector();
        Pnt p;
        Gene gg;
        double w;
        for (i = 0; i < numberOfClusters; i++) {
            p = new Pnt(dim);
            for (k = 0; k < dim; k++) {
                p.coordinates[k] = 0.0;
            }
            w = 0.0;
            for (j = 0; j < myClustering.getRightNumberSelectedGenes(before); j++) {
                gg = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(j, before);
                w += soft_memberships.get(i, j);
                for (k = 0; k < dim; k++) {
                    p.coordinates[k] += (soft_memberships.get(i, j) * gg.manifold_components_total.coordinates[k]);
                }
            }

            if (w != 0.0) {
                for (k = 0; k < dim; k++) {
                    p.coordinates[k] /= w;
                }
            } else {
                for (k = 0; k < dim; k++) {
                    p.coordinates[k] = 0.0;
                }
            }

            hard_memberships_centers_manifold_total.add(p);
        }
    }

    public void softMembershipsToHardCenters_Pca_Total() {
        int i, j, k, dim = myAGCT.data.getDimFeatures();
        hard_memberships_centers_pca_total = new Vector();
        Pnt p;
        Gene gg;
        double w;
        for (i = 0; i < numberOfClusters; i++) {
            p = new Pnt(dim);
            for (k = 0; k < dim; k++) {
                p.coordinates[k] = 0.0;
            }
            w = 0.0;
            for (j = 0; j < myClustering.getRightNumberSelectedGenes(before); j++) {
                gg = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(j, before);
                w += soft_memberships.get(i, j);
                for (k = 0; k < dim; k++) {
                    p.coordinates[k] += (soft_memberships.get(i, j) * gg.pca_components.coordinates[k]);
                }
            }

            if (w != 0.0) {
                for (k = 0; k < dim; k++) {
                    p.coordinates[k] /= w;
                }
            } else {
                for (k = 0; k < dim; k++) {
                    p.coordinates[k] = 0.0;
                }
            }

            hard_memberships_centers_pca_total.add(p);
        }
    }

    public void softMembershipsToHardCenters_Manifold_Pnt3D(int xAxis, int yAxis, int zAxis) {
        int i;
        if (hard_memberships_centers_manifold_Pnt3D == null) {
            hard_memberships_centers_manifold_Pnt3D = new Vector();
            for (i = 0; i < numberOfClusters; i++) {
                hard_memberships_centers_manifold_Pnt3D.add(new Pnt3D());
            }
        }
        if (hard_memberships_centers_manifold_total == null) {
            softMembershipsToHardCenters_Manifold_Total();
        }

        Pnt p;
        for (i = 0; i < numberOfClusters; i++) {
            p = (Pnt3D) hard_memberships_centers_manifold_Pnt3D.elementAt(i);
            p.coordinates[0] = ((Pnt) hard_memberships_centers_manifold_total.elementAt(i)).coordinates[xAxis];
            p.coordinates[1] = ((Pnt) hard_memberships_centers_manifold_total.elementAt(i)).coordinates[yAxis];
            p.coordinates[2] = ((Pnt) hard_memberships_centers_manifold_total.elementAt(i)).coordinates[zAxis];
        }
    }

    public void softMembershipsToHardCenters_Pca_Pnt3D(int xAxis, int yAxis, int zAxis) {
        int i;
        if (hard_memberships_centers_pca_Pnt3D == null) {
            hard_memberships_centers_pca_Pnt3D = new Vector<Pnt3D>();
            for (i = 0; i < numberOfClusters; i++) {
                hard_memberships_centers_pca_Pnt3D.add(new Pnt3D());
            }
        }
        if (hard_memberships_centers_pca_total == null) {
            softMembershipsToHardCenters_Pca_Total();
        }

        Pnt p;
        for (i = 0; i < numberOfClusters; i++) {
            p = (Pnt3D) hard_memberships_centers_pca_Pnt3D.elementAt(i);
            p.coordinates[0] = ((Pnt) hard_memberships_centers_pca_total.elementAt(i)).coordinates[xAxis];
            p.coordinates[1] = ((Pnt) hard_memberships_centers_pca_total.elementAt(i)).coordinates[yAxis];
            p.coordinates[2] = ((Pnt) hard_memberships_centers_pca_total.elementAt(i)).coordinates[zAxis];
        }
    }

    public void toClustering() {
        // does Nothing !
        Matrix.perror("Instance of AGCTClustering_Algorithm badly defined");
    }

    public void toClusterChoices() {
        clusterChoices = new String[numberOfClusters + 1];
        clusterChoices[0] = AGCTClustering_Algorithm.ALL_CLUSTERS;
        int i;
        for (i = 0; i < numberOfClusters; i++) {
            clusterChoices[i + 1] = new String(AGCTClustering_Algorithm.CLUSTER_NAME + i);
        }
    }

    public void fillGeneMemberships() {
        Gene gg;
        int i, j;
        Double[] memb;
        if (soft_memberships != null) {
            for (i = 0; i < myClustering.getRightNumberSelectedGenes(before); i++) {
                gg = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(i, before);
                memb = new Double[numberOfClusters];
                for (j = 0; j < numberOfClusters; j++) {
                    memb[j] = new Double(soft_memberships.get(j, i));
                }
                gg.addClusterMemberships(memb);
            }
        }
        memb = null;
    }

    /**
     * クラスタ中心と,geneの距離の総和.
     * soft_clusteringの場合,重みつき平均をとる.
     *
     * @return
     */
    public double totalDistortion() {
        if ((hard_memberships_centers == null) || (!soft_memberships_available)) {
            Matrix.perror("No hard membership distortion available");
        }

        double val = 0.0;
        for (int geneId = 0; geneId < myClustering.getRightNumberSelectedGenes(before); geneId++) {
            for (int clusterId = 0; clusterId < numberOfClusters; clusterId++) {
                val += (soft_memberships.get(clusterId, geneId) * Distortion.distortion_l22(myClustering.getRightPoint(geneId, before), (Pnt) hard_memberships_centers.elementAt(clusterId)));
            }
        }
        return val;
    }

    private double averageIntraClusterSimilarity() {
        if (!soft_memberships_available) {
            Matrix.perror("soft_memberships unavailable");
        }
        Gene gi, gj;

        int i, j, k, sizz;
        double simK, val = 0.0;
        for (k = 0; k < numberOfClusters; k++) {
            simK = 0.0;
            sizz = 0;
            for (i = 0; i < myClustering.getRightNumberSelectedGenes(before) - 1; i++) {
                for (j = i + 1; j < myClustering.getRightNumberSelectedGenes(before); j++) {
                    if ((majorityClusterFor(i, k)) && (majorityClusterFor(j, k))) {
                        sizz++;
                        gi = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(i, before);
                        gj = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(j, before);
                        simK += Gene.getSimilarity(gi, gj);
                    }
                }
            }
            if (sizz == 0) {
                Matrix.perror("empty cluster");
            }
            simK /= (double) sizz;
            val += simK;
        }
        val /= (double) numberOfClusters;
        return val;
    }

    public boolean containsIndex(int[] affect, int id) {
        int i;
        for (i = 0; i < affect.length; i++) {
            if (affect[i] == id) {
                return true;
            }
        }
        return false;
    }

    /**
     * 各クラスタ毎に,そのクラスタに属する確率が最大のgeneに関して,以下の値(そのようなgeneがなければ0)
     * を求めて,その和を返す.
     * 2 * (sum{i!=j} gene_i.dot(gene_j)) / (number of genes)
     *
     * @return
     */
    public double kernelSim() {

        Debug.debug("#### AGCTClustering_Algorithm.kernelSim() ####");
        // Cf Zass + Shashua ICCV 06, (2)
        if (!soft_memberships_available) {
            Matrix.perror("soft_memberships unavailable");
        }

        double res = 0.0;// 結果を格納する場所
        for (int clusterId = 0; clusterId < numberOfClusters; clusterId++) {// 全くらすたについて以下を計算し、和をとる。
            double simK = 0.0;// 
            int sizz = 0;
            for (int geneId = 0; geneId < myClustering.getRightNumberSelectedGenes(before); geneId++) {
                // getRightNumberSe...を知らないのでチェックしてください。
                if (majorityClusterFor(geneId, clusterId)) {
                    sizz++;
                    Gene gi = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(geneId, before);
                    for (int geneId2 = 0; geneId2 < myClustering.getRightNumberSelectedGenes(before); geneId2++) {
                        if (majorityClusterFor(geneId2, clusterId)) {
                            Gene gj = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(geneId2, before);
                            // gi,gj ともに,clusterIdに属している確率が最大.
                            simK += gi.dot(gj);// 内積を加算
                            // study the dot func.

                        }
                    }
                }
            }

            // 重複を除去(case i = j)


            // want to implement with Scala
            for (int geneId = 0; geneId < myClustering.getRightNumberSelectedGenes(before); geneId++) {
                if (majorityClusterFor(geneId, clusterId)) {
                    Gene gi = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(geneId, before);
                    simK -= gi.dot(gi);
                }
            }

            // normalization of simK
            if (sizz != 0) {
                simK /= ((double) sizz);
                res += simK;
            }
        }
        return res;
    }

    public int majorityCluster(int ex) {
        int j, index = 0;
        for (j = 1; j < numberOfClusters; j++) {
            if (soft_memberships.get(j, ex) > soft_memberships.get(index, ex)) {
                index = j;
            }
        }
        return index;
    }

    public void toSelectedGeneToHardCluster() {
        selectedGeneToHardCluster = new int[myClustering.getRightNumberSelectedGenes(before)];
        int i, j, nnc, adden, geneId, nid;
        Vector vid;
        Gene gene, ggn;

        if (AGCT.Referenced_Available) {
            if (AGCT.Max_Reference == -1) {
                Matrix.perror("AGCTClustering_Algorithm.class :: bad max reference");
            }
            referencesByCluster = new double[AGCT.Max_Reference + 1][];
            for (i = 0; i < AGCT.Max_Reference + 1; i++) {
                referencesByCluster[i] = new double[numberOfClusters];
                for (j = 0; j < numberOfClusters; j++) {
                    referencesByCluster[i][j] = 0.0;
                }
            }
        }

        for (i = 0; i < myClustering.getRightNumberSelectedGenes(before); i++) {
            selectedGeneToHardCluster[i] = majorityCluster(i);
        }

//		clusterSizesPrototypes = new int[numberOfClusters];
        clusterSizesWhole = new int[numberOfClusters];
        for (i = 0; i < numberOfClusters; i++) {
//			clusterSizesPrototypes[i] = 0;
            clusterSizesWhole[i] = 0;
        }
        for (i = 0; i < myClustering.getRightNumberSelectedGenes(before); i++) {
            nnc = selectedGeneToHardCluster[i];
            adden = 0;
            if ((before) || (Prototype.No_Reduction)) {
                adden = 1;
            } else {
                geneId = myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i];
                vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(geneId));
                adden = vid.size();
            }

            if (!before) {
                geneId = myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i];
                gene = myAGCT.data.getMyDomain().getGenes().get(geneId);

                if (gene.isReferenced()) {
                    referencesByCluster[gene.getTypeReferenced()][nnc] += 1.0;
                }

                if (!Prototype.No_Reduction) {
                    vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(geneId));
                    for (j = 1; j < vid.size(); j++) {
                        nid = ((Integer) vid.elementAt(j)).intValue();
                        ggn = (Gene) myAGCT.data.getMyDomain().getGenes().get(nid);

                        if (ggn.isReferenced()) {
                            referencesByCluster[ggn.getTypeReferenced()][nnc] += 1.0;
                        }
                    }
                }
            }

//			clusterSizesPrototypes[nnc] += 1;
            clusterSizesWhole[nnc] += adden;
        }

        myClustering.nclusters = numberOfClusters;
    }

    /**
     * geneがclusterに属している確率が最大の時, trueを返す.
     *
     * @param geneId
     * @param clusterId
     * @return
     */
    public boolean majorityClusterFor(int geneId, int clusterId) {
        // ex MUST be between 0 and numberSelectedGenes

        int tmpClusterId;
        for (tmpClusterId = 0; tmpClusterId < numberOfClusters; tmpClusterId++) {
            if ((tmpClusterId != clusterId) && (soft_memberships.get(tmpClusterId, geneId) > soft_memberships.get(clusterId, geneId))) {
                return false;
            }
        }
        return true;
    }

    public boolean majorityGeneCluster(int ng, int nc) {
        if ((soft_clustering_computed) || (hard_clustering_computed)) {
            int i;
            for (i = 0; i < numberOfClusters; i++) {
                if ((i != nc) && (soft_memberships.get(i, ng) > soft_memberships.get(nc, ng))) // System.out.println("Ng = " + ng + ", Nc = " + nc +
                // ", i = " + i);
                {
                    return false;
                }
            }
        }
        return true;
    }

    public void toHardMemberships() {
        Debug.debug("AGCTClustering_Algortithm # toHardMemberships");
        hard_memberships_centers = new Vector();
        int i, j, k, dd = getDimension(before);
        double tot;
        Pnt p;
        for (i = 0; i < numberOfClusters; i++) {
            p = new Pnt(dd);
            for (k = 0; k < dd; k++) {
                p.coordinates[k] = 0.0;
            }

            tot = 0.0;
            for (j = 0; j < myClustering.getRightNumberSelectedGenes(before); j++) {
                tot += soft_memberships.get(i, j);
                for (k = 0; k < dd; k++) {
                    p.coordinates[k] += (soft_memberships.get(i, j) * myClustering.getRightPoint(j, before).coordinates[k]);
                }
            }

            if (tot == 0.0) {
                Debug.debug("AGCTClusterig_Algortihm # toHardMemberships, tot == 0.0");
                Matrix.perror("No members for cluster " + i);
            }

            for (k = 0; k < dd; k++) {
                p.coordinates[k] /= tot;
            }
            hard_memberships_centers.add(p);
        }

        finalClusteringObjectiveFunction = totalDistortion();
    }

    // Bregman Arthur Vassilvitskii init.
    public void update_for_chosen(int newCenterGeneId, boolean[] chosen, double[] min_distortion, int dist_used, double parametre) {
        Pnt p, c = myClustering.getRightPoint(newCenterGeneId, before);
        int j, dim = myClustering.getRightNumberSelectedGenes(before);
        double dist;
        min_distortion[newCenterGeneId] = 0.0;

        for (j = 0; j < dim; j++) {
            if (!chosen[j]) {
                p = myClustering.getRightPoint(j, before);

                // System.out.println(dist_used);

                if (dist_used == -1) {
                    dist = Distortion.distortion_l22(p, c);
                } else {
                    dist = Distortion.Bregman(p, c, dist_used, parametre);
                }

                if ((dist < 0.0) && (dist > -MatrixConfig.Precision_For_Eigensystems)) {
                    dist = 0.0;
                }

                if (dist < 0.0) {
                    Matrix.perror("Negative distortion (" + dist + ") for j = " + j + " in AV seeding");
                }

                if ((min_distortion[j] == -1.0) || (dist < min_distortion[j])) {
                    min_distortion[j] = dist;
                }
            }
        }
    }

    /**
     * http://www.stanford.edu/~darthur/kMeansPlusPlus.pdfを参照
     * 次のように選ぶ
     * 1.ランダムに一つcenterを選ぶ
     * 2.pに対し、既に選ばれた中心のうち、最も近いものとの距離をD(p)とする。
     * この時、xが選ばれる確率が、D(x)に比例するように、ランダムに次の中心を選ぶ
     * (* 論文ではD(x)^2に比例と書いてあるが、距離として既に^2したもの(L2)を使っているのでこれでよい*)
     * 3.2をK-1回繰り返す(K is number of clusters)
     * (*つまり、距離関数として、L2でのBregmanを使用するということ*)
     *
     * @param dist_used Distortion.javaに対応.
     *                  -1 -> distortion_l22
     *                  i >= 0 -> {"L22","Kullback-Leibler","Itakura-Saito","Amari (alpha)","p-norm (p)"}[i].
     * @param parametre
     */
    public void initHardMembershipArthurVassilvitskii(int dist_used, double parametre) {
        int i, j, numberOfGenes = myClustering.getRightNumberSelectedGenes(before);
        int[] centerGeneIds = new int[numberOfClusters];
        // 既に選ばれた中心のうち、最も近いものとの距離
        double[] D = new double[numberOfGenes];
        boolean[] chosen = new boolean[numberOfGenes];
        Gene c, p;
        double dist, totdist, ra, binf, bsup;

        for (i = 0; i < numberOfClusters; i++) {
            centerGeneIds[i] = -1;
        }

        for (i = 0; i < numberOfGenes; i++) {
            D[i] = -1.0;
        }

        //	1.ランダムに一つcenterを選ぶ
        centerGeneIds[0] = AGCT.RandomGenerator.nextInt(numberOfGenes);

        AGCTCounter cc = new AGCTCounter("AV Initialization", numberOfClusters - 1);
        for (i = 1; i < numberOfClusters; i++) {
            chosen[centerGeneIds[i - 1]] = true;
            update_for_chosen(centerGeneIds[i - 1], chosen, D, dist_used, parametre);

            totdist = 0.0;
            for (j = 0; j < numberOfGenes; j++) {
                if (!chosen[j]) {
                    totdist += D[j];
                }
            }
            if (totdist == 0.0) {
                Matrix.perror("Zero total distortion in AV seeding");
            }

            ra = (AGCT.RandomGenerator.nextDouble()) * totdist;

            if (ra == 0.0) {
                Matrix.perror("No more seed to select");
            }

            binf = 0.0;
            j = 0;
            while ((j < numberOfGenes) && ((chosen[j]) || (ra < binf) || (ra > binf + D[j]))) {
                if (!chosen[j]) {
                    binf += D[j];
                }
                j++;
            }
            // xが選ばれる確率が、D(x)に比例するように、ランダムに次の中心を選ぶ　
            centerGeneIds[i] = j;
            cc.increment();
        }
        cc.end();

        for (i = 0; i < numberOfClusters - 1; i++) {
            for (j = i + 1; j < numberOfClusters; j++) {
                if (centerGeneIds[i] == centerGeneIds[j]) {
                    Matrix.perror("AV seeded twice the same point");
                }
            }
        }

        for (i = 0; i < numberOfClusters; i++) {
            hard_memberships_centers.add(new Pnt(myClustering.getRightPoint(centerGeneIds[i], before))); // Correction
        }        // !
    }
}