package test;

import gene.Gene;

import java.util.Random;

public class AGCTClustering_CP extends AGCTClustering_Algorithm {

    public static Matrix myCW;
    public static boolean CopyOfSimilarityMatrixIsAlreadyComputed;

    public static void init() {
        if (AGCT.MYDEBUG)
            AGCT.debug("AGCTClustering_CP.init()");
        myCW = null;
        CopyOfSimilarityMatrixIsAlreadyComputed = false;
    }

    /**
     * myCWをCと同じサイズに初期化している
     *
     * @param C
     */
    public static void toCW(Matrix C) {
        myCW = new Matrix("Copy of" + C.name, C.rowCount(), C.colCount());
        CopyOfSimilarityMatrixIsAlreadyComputed = true;
    }

    AGCTClustering_CP(Clustering myc, AGCT mya, int nc, boolean bef) {
        super(myc, mya, nc, bef, 0);
        setSeedingType(-1);
        //		seedingType = -1;
        //		myReferenceName = "CP";
        setMyReferenceName("CP");
    }

    //	public void plot(JAGCTVisualizationPane vv) {
    //		plotClusterCenters(vv);
    //	}

    //	public void plotWithAnnotations(JAGCTVisualizationPane vv) {
    //		plotClusterCentersWithAnnotations(vv);
    //	}

    public String toString() {
        String val = "CP-factorization with " + getNumberOfClusters() + " clusters.";
        return val;
    }

    /**
     * まずランダムにK個(geneの中から)中心点を選ぶ。
     * それから、各geneに対し、それに最も近い中心と、二番目に近い中心を選び、その二つに対する所属度をともに 1/sqrt(2) にする。
     * ただし、clusterが二つの場合は、二つに対する所属度をランダムに割り振る。(* 実際は、
     * v = (1 + nextInt(9))/10.0 とし、(sqrt(v),sqrt(1-v)) としている*)
     * どの遺伝子も、∑所属度^2 = 1 を満たす。(どの中心点も、ではない)
     */
    public void initSoftMembershipForgy() {
        int i, j, k, ib1, ib2;
        double dbcur = 0.0, db1 = 0.0, db2 = 0.0, vd1, vd2;
        boolean deja = false;
        Gene gg1, gg2, ggj, ggC;
        Random r = new Random(AGCT.DEFAULT_RANDOMSEED);

        int[] indexes_to_pick = new int[getNumberOfClusters()];
        for (i = 0; i < getNumberOfClusters(); i++)
            indexes_to_pick[i] = 0;

        for (i = 0; i < getNumberOfClusters(); i++) {
            do {
                //ランダムにindexを選ぶ
                j = r.nextInt(getMyClustering().getRightNumberSelectedGenes(getBefore()));
                k = 0;
                deja = false;
                do {
                    if (indexes_to_pick[k] == j) {
                        deja = true;
                    }
                    k++;
                } while ((k < i) && (deja == false));
            } while (deja == true);
            //ランダムにK個indexをえらぶ
            indexes_to_pick[i] = j;
        }

        for (i = 0; i < getNumberOfClusters(); i++)
            for (j = 0; j < getMyClustering().getRightNumberSelectedGenes(getBefore()); j++)
                // まずは、全て所属度０にする
                getSoft_memberships().set(i, j, 0.0);

        for (j = 0; j < getMyClustering().getRightNumberSelectedGenes(getBefore()); j++) {
            ib1 = 0;
            ib2 = 1;

            gg1 = getMyClustering().getRightGeneSelectedGeneNumberToGeneNumber(indexes_to_pick[ib1], getBefore());
            gg2 = getMyClustering().getRightGeneSelectedGeneNumberToGeneNumber(indexes_to_pick[ib2], getBefore());
            ggj = getMyClustering().getRightGeneSelectedGeneNumberToGeneNumber(j, getBefore());

            db1 = Gene.squaredMagnitudeDifference(gg1, ggj);// 中心０からの距離^2
            db2 = Gene.squaredMagnitudeDifference(gg2, ggj);// 中心１からの距離^2

            if (db1 > db2) {
                int it = ib1;
                ib1 = ib2;
                ib2 = it;

                double dt = db1;
                db1 = db2;
                db2 = dt;
            }// db1 <= db2

            if (getNumberOfClusters() > 2) {
                for (i = 2; i < getNumberOfClusters(); i++) {
                    ggC = getMyClustering().getRightGeneSelectedGeneNumberToGeneNumber(indexes_to_pick[i], getBefore());
                    dbcur = db2 = Gene.squaredMagnitudeDifference(ggC, ggj);// 中心iからの距離^2

                    // dbcur < db1 <= db2 なら
                    if (dbcur < db1) {
                        ib2 = ib1;
                        db2 = db1;

                        ib1 = i;
                        db1 = dbcur;//db1 = dbcur,db2 = db1 に置き換え
                    }
                    // db1 <= dbcur < db2 なら
                    else if (dbcur < db2) {
                        ib2 = i;
                        db2 = dbcur;// db2 = dbcur に置き換え
                    }
                }// for文終わり

                // つまり、db1が最小距離、db2が二番目の最小距離、ib1,ib2はそのindexとなっている

                getSoft_memberships().set(ib1, j, 1.0 / Math.sqrt(2.0));
                getSoft_memberships().set(ib2, j, 1.0 / Math.sqrt(2.0));
            } else {
                i = 1 + r.nextInt(9);
                vd1 = ((double) i) / 10.0;
                vd2 = 1.0 - vd1;

                // System.out.println("vd1 = " + vd1 + ", vd2 = " + vd2);

                if ((vd1 <= 0.0) || (vd2 <= 0.0))
                    Matrix.perror("Not all strictly positive");

                getSoft_memberships().set(0, j, Math.sqrt(vd1));
                getSoft_memberships().set(1, j, Math.sqrt(vd2));
            }
        }
    }

    /**
     * 1.受け取ったCopyOfSimilarityMatrix にfastDoublyStochasticApproximationを施し,
     * 二重確率行列(行,列の和が1)にする.
     * 2.その行列を使って,soft_membershipsをcompletePositiveFactorization
     *
     * @param CopyOfSimilarityMatrix
     */
    public void clustering_Zass_Shashua(Matrix CopyOfSimilarityMatrix) {
        if (!AGCTClustering_CP.CopyOfSimilarityMatrixIsAlreadyComputed) {
            AGCTClustering_CP.toCW(CopyOfSimilarityMatrix);
            myCW.show(10, 10);
            AGCTClustering_CP.myCW.fastDoublyStochasticApproximation();
            myCW.show(10, 10);

            Matrix.assertDoublyStochastic(AGCTClustering_CP.myCW);
            Matrix.assertSymmetric(AGCTClustering_CP.myCW);
        }
        getSoft_memberships().completePositiveFactorization(getMyAGCT().data.getMyDomain(), AGCTClustering_CP.myCW);

        getSoft_memberships().normalize();
        Matrix.assertColumnStochastic(getSoft_memberships());
        reduce();
        Matrix.assertColumnStochastic(getSoft_memberships());
        setSoft_memberships_available(true);
        //		soft_memberships_available = true;
        setPlotCentersAvailable(true);
        //		plotCentersAvailable = true;
    }

    public void reduce() {
        boolean notEmpty[] = new boolean[getNumberOfClusters()];
        int i, j, fc = 0;
        double tot;

        for (i = 0; i < getNumberOfClusters(); i++) {
            tot = 0.0;
            for (j = 0; j < getMyClustering().getRightNumberSelectedGenes(getBefore()); j++) {
                tot += getSoft_memberships().get(i, j);
            }
            if (tot == 0.0) {
                notEmpty[i] = false;
                if (AGCT.DoDebug)
                    System.out.println("No members for cluster " + i);
            } else {
                notEmpty[i] = true;
                fc++;
            }
        }

        if (fc < getNumberOfClusters()) {
            Matrix nsm = new Matrix(getSoft_memberships().name, fc, getSoft_memberships().colCount());

            fc = 0;
            for (i = 0; i < getNumberOfClusters(); i++)
                if (notEmpty[i] == true) {
                    for (j = 0; j < getMyClustering().getRightNumberSelectedGenes(getBefore()); j++)
                        nsm.set(fc, j, getSoft_memberships().get(i, j));
                    fc++;
                }

            setSoft_memberships(nsm);
            //			soft_memberships = nsm;
            nsm = null;
            setNumberOfClusters(getSoft_memberships().rowCount());
            //			numberOfClusters = getSoft_memberships().rowCount();
            getMyClustering().nclusters = getNumberOfClusters();
        }

        notEmpty = null;
    }

    public void toClustering() {
        initSoftMembershipForgy();
        clustering_Zass_Shashua(getMyAGCT().data.getCopyOfSimilarityMatrix()/*Wのコピー*/);
        toClusterChoices();
        setSoft_clustering_computed(true);
        //		soft_clustering_computed = true;
        setHard_clustering_computed(false);
        //		hard_clustering_computed = false;

        toHardMemberships();

        toSelectedGeneToHardCluster();

        setOk_for_statistics(true);
        //		ok_for_statistics = true;

        ControlProcess.put("softClusteringProcessed", true);
        ControlProcess.put("ClusteringProcessed_CP", true);
    }
}
