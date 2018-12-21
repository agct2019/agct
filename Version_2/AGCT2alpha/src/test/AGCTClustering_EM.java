package test;


import forDebug.Debug;

public class AGCTClustering_EM extends AGCTClustering_Algorithm {

    private int dimensionOfPoint;
    /* numberOfClusters要素で, 要素は,dimensionOfPoint * dimensionOfPoint のマトリクス.
     * 対角要素のみ 1e-10 で初期化.
     */
    private Matrix[] sig;
    /*
     * numberOfClusters要素.  1.0/numberOfClustersで初期化.
     */
    private double[] frac;
    /*
     *
     */
    private double[] lndets;
    private double loglike;
    private boolean stop_choleski;

    AGCTClustering_EM(Clustering myClustering, AGCT myAGCT, int numClusters, int seeding, boolean bef, int maxIterations) {
        super(myClustering, myAGCT, numClusters, bef, maxIterations);
        setSeedingType(seeding);
        setMyReferenceName("EM");
        setFinalClusteringObjectiveFunction(-1);

        sig = null;
        lndets = null;
        frac = null;
        stop_choleski = false;
    }

//	 public void plot(JAGCTVisualizationPane vv) {
//		plotClusterCenters(vv);
//	}

//	 public void plotWithAnnotations(JAGCTVisualizationPane vv) {
//		plotClusterCentersWithAnnotations(vv);
//	}

    public String toString() {
        String val = "Expectation Maximization with " + getNumberOfClusters() + " Gaussians.\n";
        val += "Initialization : " + AGCTClustering_Algorithm.SEEDING_NAME[getSeedingType()] + ".\n";
        // val += "Distortion : " + Distortion.DISTORTION_NAME[distortionType] +
        // ".\n";
        val += "Points are " + Clustering.referenceStringCenters[getMyClustering().indexReferenceStringCenters] + ".";
        return val;
    }

    /*
     * choose initial centers
     * initialize sig,lndets,and frac.
     */
    private void initSoftMembershipsEM() {
        Debug.debug("AGCTClustering_EM # initSoftMembershipsEM()");
        // At first, choose initial centers.
        if (getSeedingType() == 0)    // Forgy
            initHardMembershipForgy();
        else if (getSeedingType() == 1)
            initHardMembershipArthurVassilvitskii(-1, -1.0);
        else
            Matrix.perror("Invalid seeding type");

        final int numberOfClusters = getNumberOfClusters();
        final int dimensionOfPoint = this.dimensionOfPoint;

        sig = new Matrix[numberOfClusters];
        for (int k = 0; k < numberOfClusters; k++)
            sig[k] = new Matrix("Covariance of cluster " + k, dimensionOfPoint, dimensionOfPoint);
        lndets = new double[numberOfClusters];
        frac = new double[numberOfClusters];

        for (int k = 0; k < numberOfClusters; k++) {
            // 1.0/numberOfClustersで初期化
            frac[k] = 1.0 / (double) numberOfClusters;
            for (int i = 0; i < dimensionOfPoint; i++) {
                for (int j = 0; j < dimensionOfPoint; j++)
                    sig[k].set(i, j, 0.0);
                sig[k].set(i, i, 1.0E-10);
            }
        }
    }

    private double EStep() {
        final int numberOfClusters = getNumberOfClusters();
        final int dimensionOfPoint = this.dimensionOfPoint;
        final int numberOfGenes = getMyClustering().getRightNumberSelectedGenes(getBefore());
        double tmp, max, oldloglike = loglike, valret;
        double[] u = new double[dimensionOfPoint], v = new double[dimensionOfPoint];

        for (int k = 0; k < numberOfClusters; k++) {
            boolean pass_choleski = Matrix.test_choltmp(sig[k]);
            if (!pass_choleski) {
                stop_choleski = true;
                break;
            }
        }

        if (!stop_choleski) {
            for (int k = 0; k < numberOfClusters; k++) {
                Matrix.choltmp(sig[k]);
                lndets[k] = sig[k].logdet();
                for (int n = 0; n < numberOfGenes; n++) {
                    for (int m = 0; m < dimensionOfPoint; m++)
                        u[m] = getMyClustering().getRightPoint(n, getBefore()).coordinates[m] - getHard_memberships_centers().elementAt(k).coordinates[m];
                    sig[k].elsolve(u, v);
                    double sum = 0.0;
                    for (int m = 0; m < dimensionOfPoint; m++)
                        sum += (v[m] * v[m]);
                    getSoft_memberships().set(k, n,
                            -0.5 * (sum + lndets[k]) + Math.log(frac[k])
                    );
                }
            }
            loglike = 0.0;
            for (int geneId = 0; geneId < numberOfGenes; geneId++) {
                max = -99.9E99;
                for (int clusterId = 0; clusterId < numberOfClusters; clusterId++)
                    if (getSoft_memberships().get(clusterId, geneId) > max)
                        max = getSoft_memberships().get(clusterId, geneId);
                double sum = 0.0;
                for (int clusterId = 0; clusterId < numberOfClusters; clusterId++)
                    sum += Math.exp(getSoft_memberships().get(clusterId, geneId) - max);
                tmp = max + Math.log(sum);
                for (int k = 0; k < numberOfClusters; k++)
                    getSoft_memberships().set(k, geneId,
                            Math.exp(getSoft_memberships().get(k, geneId) - tmp));
                loglike += tmp;
            }
            setFinalClusteringObjectiveFunction(loglike);
            valret = loglike - oldloglike;
        } else
            valret = 0.0;

        return valret;
    }

    private void MStep() {
        int numberOfClusters = getNumberOfClusters();
        int dimensionOfPoint = this.dimensionOfPoint;
        int numberOfGenes = getMyClustering().getRightNumberSelectedGenes(getBefore());
        for (int clusterId = 0; clusterId < numberOfClusters; clusterId++) {
            double wgt = 0.0;
            for (int geneId = 0; geneId < numberOfGenes; geneId++)
                wgt += getSoft_memberships().get(clusterId, geneId);
            frac[clusterId] = wgt / (double) numberOfGenes;
            for (int pointId = 0; pointId < dimensionOfPoint; pointId++) {
                double sum = 0;
                for (int geneId2 = 0; geneId2 < numberOfGenes; geneId2++) {
                    sum += getSoft_memberships().get(clusterId, geneId2) * getMyClustering().getRightPoint(geneId2, getBefore()).coordinates[pointId];
                }
                (getHard_memberships_centers().elementAt(clusterId)).coordinates[pointId] = sum / wgt;
                for (int pointId2 = 0; pointId2 < dimensionOfPoint; pointId2++) {
                    double sum2 = 0;
                    for (int geneId2 = 0; geneId2 < numberOfGenes; geneId2++) {
                        sum2 +=
                                getSoft_memberships().get(clusterId, geneId2)

                                        * (getMyClustering().getRightPoint(geneId2, getBefore()).coordinates[pointId] -
                                        getHard_memberships_centers().elementAt(clusterId).coordinates[pointId])

                                        * (getMyClustering().getRightPoint(geneId2, getBefore()).coordinates[pointId2] -
                                        getHard_memberships_centers().elementAt(clusterId).coordinates[pointId2]);
                    }
                    sig[clusterId].set(pointId, pointId2, sum2 / wgt);
                }
            }
        }
    }

    private void EM() {
        Debug.debug("AGCTClustering_EM # EM");
        Debug.debug("nMaxIterations = " + getnMaxIterations());
        int niter = 0;
        double delta;
        boolean stop = false;
        AGCTCounter ccc;

        if (getnMaxIterations() > 0)
            ccc = new AGCTCounter("Clustering EM", getnMaxIterations());
        else
            ccc = new AGCTCounter("Clustering EM", 0);

        do {
            delta = Math.abs(EStep());

            if (getnMaxIterations() > 0)
                ccc.increment();
            else
                ccc.setPercent((int) (100.0 * AGCTClustering_Algorithm.EPSILON_EM / delta));

            if ((stop_choleski) || (delta < AGCTClustering_Algorithm.EPSILON_EM)) {
                if (stop_choleski) {
                    Debug.debug("stop_choleski == true");
                }
                if (delta < AGCTClustering_Algorithm.EPSILON_EM) {
                    Debug.debug("delta = " + delta);
                    Debug.debug("EPSILON_EM = " + AGCTClustering_Algorithm.EPSILON_EM);
                }
                stop = true;
            } else
                MStep();

            niter++;

            if (getnMaxIterations() > 0)
                if (niter >= getnMaxIterations())
                    stop = true;

        } while (!stop);
        ccc.end();

        setPlotCentersAvailable(true);
    }

    public void toClustering() {
        // クラスタリングに使用する点の次数
        dimensionOfPoint = getDimension(getBefore());

        Debug.debug("BEGIN initSoftMembershipsEM()");
        initSoftMembershipsEM();
        Debug.debug("END initSoftMembershipsEM()");
        Debug.debug("BEGIN EM()");
        EM();
        Debug.debug("END EM()");
        Debug.debug("BEGIN toClusterChoices()");
        toClusterChoices();
        Debug.debug("END toClusterChoices()");
        setSoft_memberships_available(true);
        setSoft_clustering_computed(true);
        setHard_clustering_computed(false);
        Debug.debug("BEGIN toHardMemberships");
        toHardMemberships();
        Debug.debug("END toHardMemberships");

        Debug.debug("BEGIN toSelectedGeneToHardCluster");
        toSelectedGeneToHardCluster();
        Debug.debug("END toSelectedGeneToHardCluster");

        setOk_for_statistics(true);
        ControlProcess.put("softClusteringProcessed", true);
        ControlProcess.put("ClusteringProcessed_EM", true);
    }
}