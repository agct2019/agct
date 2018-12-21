package clustering;

import test.*;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

public class AGCTClustering_KM extends AGCTClustering_Algorithm {

    private int distortionType;
    // values are given in Distortion.java

    private double dist_param;
    // Amari = alpha
    // p-norm = p

    private int n_IterativeKMeans;
    // plateau size for n_IterativeKMeans (or -1)

    private boolean sameDistortion;

    private int[] affect;

    // gives for each gene its cluster

    public AGCTClustering_KM(Clustering clustering, AGCT agct, int _numberOfClusters, int seedingType, boolean keepDist, int di, double va, boolean bef, int nMI, int nI) {
        super(clustering, agct, _numberOfClusters, bef, nMI);
        setSeedingType(seedingType);
        // seedingType = ti;
        distortionType = di;
        sameDistortion = keepDist;
        n_IterativeKMeans = 0;

        if ((AGCT.DoDebug) && (nI >= 1000))
            System.out.println("Value for Bregman Iterative K-Means = " + nI + " is too large...");

        if ((nI > 0) && (nI < 1000))
            n_IterativeKMeans = nI;

        if ((distortionType == 3) && (va != -1.0))
            dist_param = va;
        else if (distortionType == 3)
            dist_param = 0.0;

        if ((distortionType == 4) && (va != -1.0))
            dist_param = va;
        else if (distortionType == 4)
            dist_param = 2.0;
        setMyReferenceName("KM");
        setFinalClusteringObjectiveFunction(-1);
    }

    public String toString() {
        int a = 0;
        if (sameDistortion)
            a = 1;

        String val = "K-Means with " + getNumberOfClusters() + " centers.\n";
        val += "Initialization : " + AGCTClustering_Algorithm.SEEDING_NAME[getSeedingType() + a] + ".\n";
        val += "Distortion : " + Distortion.DIVERGENCES_NAME[distortionType];
        if (distortionType == 3)
            val += " (with alpha = " + dist_param + ")\n";
        else if (distortionType == 4)
            val += " (with p = " + dist_param + ")\n";
        else
            val += "\n";
        val += "Points are " + Clustering.referenceStringCenters[getMyClustering().indexReferenceStringCenters] + ".";
        return val;
    }

    public int closestCenterIndex(int id) {
        Pnt p = getMyClustering().getRightPoint(id, getBefore());
        int i, v = -1;
        double dcur, dmin = -1;
        for (i = 0; i < getHard_memberships_centers().size(); i++) {
            dcur = Distortion.Bregman((Pnt) p, (Pnt) getHard_memberships_centers().elementAt(i), distortionType, dist_param);
            if ((i == 0) || (dcur < dmin)) {
                dmin = dcur;
                v = i;
            }
        }
        return v;
    }

    /**
     * **
     * public int closestPointIndexExclusive(Hashtable t, int id){ //The main
     * difference with closestPointIndex is that we prevent a point to be chosen
     * as the representative of 2+ clusters Pnt p = (Pnt)
     * getHard_memberships_centers().elementAt(id); int i, v = -1, trueInd,
     * testv = -1; double dcur, dmin = -1.0, testdmin = -1.0; boolean firstTest
     * = true;
     * <p/>
     * for (i=0;i<myClustering.getRightNumberSelectedGenes(getBefore());i++){
     * dcur = Distortion.Bregman(myClustering.getRightPoint(i, getBefore()), p,
     * distortionType, dist_param); trueInd =
     * myClustering.getRightIndexSelectedGeneNumberToGeneNumber(i, getBefore());
     * <p/>
     * if (trueInd == -1)
     * Matrix.perror("AGCTClustering_KM.class :: wrong trueInd = -1");
     * <p/>
     * if ( (i==0) || (dcur < testdmin) ){ testdmin = dcur; testv = trueInd; }
     * <p/>
     * if ( (!Prototype.currentCentersContain(trueInd)) && ( firstTest || (dcur
     * < dmin) ) ){ firstTest = false; dmin = dcur; v = trueInd; } }
     * <p/>
     * if (v != testv) System.out.println("Cluster " + id + " :: " + testv +
     * " --> " + v);
     * <p/>
     * return v; }
     * ***
     */

    public int closestPointIndexExclusive(Hashtable t, int id, Hashtable bum) {
        // The main difference with closestPointIndex is that we prevent a point
        // to be chosen as the representative of 2+ clusters
        // also, faster than the old method.
        Pnt p = (Pnt) getHard_memberships_centers().elementAt(id);
        int i, v = -1, trueInd, testv = -1;
        double dcur, dmin = -1.0, testdmin = -1.0;
        boolean firstTest = true;

        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++) {
            dcur = Distortion.Bregman(getMyClustering().getRightPoint(i, getBefore()), p, distortionType, dist_param);
            trueInd = getMyClustering().getRightIndexSelectedGeneNumberToGeneNumber(i, getBefore());

            if (trueInd == -1)
                Matrix.perror("AGCTClustering_KM.class :: wrong trueInd = -1");

            if ((i == 0) || (dcur < testdmin)) {
                testdmin = dcur;
                testv = trueInd;
            }

            if ((!bum.containsKey(new Integer(trueInd))) && (firstTest || (dcur < dmin))) {
                firstTest = false;
                dmin = dcur;
                v = trueInd;
            }
        }

        if (v != testv)
            System.out.println("Cluster " + id + " :: " + testv + " --> " + v);

        return v;
    }

    public int closestPointIndex(int id) {
        Pnt p = (Pnt) getHard_memberships_centers().elementAt(id);
        int i, v = -1;
        double dcur, dmin = -1;
        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++) {
            dcur = Distortion.Bregman(getMyClustering().getRightPoint(i, getBefore()), p, distortionType, dist_param);
            if ((i == 0) || (dcur < dmin)) {
                dmin = dcur;
                v = getMyClustering().getRightIndexSelectedGeneNumberToGeneNumber(i, getBefore());
            }
        }
        return v;
    }

    public boolean Iterative_KMeans(int[] affect) {
        Pnt p;
        int i;
        int[] cursizes = new int[getNumberOfClusters()];
        int[] ret = new int[2];
        int[] step_chosen = new int[getMyClustering().getRightNumberSelectedGenes(getBefore())];
        boolean stop = false;
        boolean atleastone = false;

        reduce_affect();
        i = 0;
        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++)
            step_chosen[i] = -1;
        for (i = 0; i < getNumberOfClusters(); i++)
            cursizes[i] = 0;
        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++)
            cursizes[affect[i]]++;

        i = 0;
        do {
            IterativeKMeans_next_point(affect, step_chosen, cursizes, getHard_memberships_centers(), ret);

            if (ret[0] != -1) {
                /*
                 * System.out.println(" Confirmation --- my best bet :: " +
				 * ret[0] + " in cluster " + affect[ret[0]] +
				 * " goes to cluster " + ret[1]);
				 * System.out.println(" Point :: " +
				 * myClustering.getRightPoint(ret[0], getBefore()));
				 * System.out.print(" Cluster sizes before :: "); for
				 * (i=0;i<cursizes.length;i++) System.out.print(" (" + i + ", "
				 * + cursizes[i] + ") "); System.out.println();
				 * System.out.println(" Center before  :: " + (Pnt)
				 * getHard_memberships_centers().elementAt(ret[1]));
				 */

                // creation de cluster ?
                if (cursizes[ret[1]] == 0) {
                    if (getHard_memberships_centers().size() != ret[1])
                        Matrix.perror("JAGCTClustering_KM :: cluster size mismatch");
                    p = new Pnt(getDimension(getBefore()));
                    getHard_memberships_centers().add(p);
                }
                center_add(getMyClustering().getRightPoint(ret[0], getBefore()), (Pnt) getHard_memberships_centers().elementAt(ret[1]), cursizes[ret[1]]);
                center_remove(getMyClustering().getRightPoint(ret[0], getBefore()), (Pnt) getHard_memberships_centers().elementAt(affect[ret[0]]), cursizes[affect[ret[0]]]);
                cursizes[affect[ret[0]]]--;
                cursizes[ret[1]]++;
                affect[ret[0]] = ret[1];

				/*
				 * System.out.print("\n Cluster sizes after :: "); for
				 * (i=0;i<cursizes.length;i++) System.out.print(" (" + i + ", "
				 * + cursizes[i] + ") "); System.out.println();
				 * System.out.println(" Center after  :: " + (Pnt)
				 * getHard_memberships_centers().elementAt(ret[1]));
				 */
            }
            i++;
            if ((ret[0] != -1) && (ret[1] == -1))
                Matrix.perror("JAGCTClustering_KM :: bad correspondence 1");
            if ((ret[0] == -1) && (ret[1] != -1))
                Matrix.perror("JAGCTClustering_KM :: bad correspondence 2");

            if (ret[0] != -1)
                atleastone = true;

            if ((ret[0] == -1) || (i >= n_IterativeKMeans))
                stop = true;
            ret[0] = ret[1] = -1;
        } while (!stop);

        return atleastone;
    }

    public void IterativeKMeans_next_point(int[] affect, int[] step_chosen, int[] cursizes, Vector curcent, int[] ret) {
        // step_chosen = -1 (jamais bouge), ou bien contient le num�ro de
        // cluster ou il bascule si >= 0
        // cursizes = tailles courantes des clusters
        // curcent = vector de Pnt(s)
        // BEWARE : do not modify affect

        int i, j, index_best_bet_i = -1, index_best_bet_cluster = -1;
        double cur_bet, best_bet = 0.0;
        boolean stop;
        if (step_chosen.length != getMyClustering().getRightNumberSelectedGenes(getBefore()))
            Matrix.perror("AGCTClustering_KM :: Length of step_chosen does not match");
        if (cursizes.length != getNumberOfClusters())
            Matrix.perror("AGCTClustering_KM :: Length of cursizes does not match");

        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++)
            if (step_chosen[i] == -1) {
                j = 0;
                stop = false;
                do {
                    if (affect[i] != j) {
                        cur_bet = delta_div(getMyClustering().getRightPoint(i, getBefore()), affect[i], j, cursizes, curcent);

                        // System.out.println(" " + i + " :: " + affect[i] +
                        // " -> " + j + " == " + cur_bet);

                        if (((index_best_bet_i == -1) && (cur_bet > 0.0)) || ((index_best_bet_i != -1) && (cur_bet > best_bet))) {
                            index_best_bet_i = i;
                            index_best_bet_cluster = j;
                            best_bet = cur_bet;
                        }
                    }
                    if ((j >= getNumberOfClusters() - 1) || (cursizes[j] == 0))
                        stop = true;
                    j++;
                } while (!stop);
            }

        ret[0] = index_best_bet_i;
        ret[1] = index_best_bet_cluster;

		/*
		 * if (ret[0] != -1) System.out.println(" My best bet :: " + ret[0] +
		 * " in cluster " + affect[ret[0]] + " goes to cluster " + ret[1]);
		 */
    }

    public void reduce_affect() {
        // before running Iterative k means, we need to modify affect to remove
        // the centers that do not exist;
        int i = 0, j, max_value;
        do {
            if (!containsIndex(affect, i)) {
                for (j = 0; j < getMyClustering().getRightNumberSelectedGenes(getBefore()); j++)
                    if (affect[j] > i)
                        affect[j]--;
            } else
                i++;
        } while (i < getNumberOfClusters());

        // Check that the max value of affect is == size of hard_membership
        // centers
        max_value = -1;
        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++)
            if ((max_value == -1) || (affect[i] > max_value))
                max_value = affect[i];
        if (max_value != getHard_memberships_centers().size() - 1)
            Matrix.perror("AGCTClustering_Algorithm :: bad sizes (max_value = " + max_value + ", should be = " + (getHard_memberships_centers().size() - 1) + ")");
    }

    public double delta_div(Pnt pt, int clusterori, int clusterfin, int[] cursizes, Vector curcent) {
        double A, B, C, D, val;
        A = B = C = D = 0.0;
        int i, tailleori = cursizes[clusterori];
        if (tailleori == 0)
            Matrix.perror("AGCTClustering_KM :: Cluster empty");
        int taillefin = cursizes[clusterfin];
        Pnt oriap = new Pnt(getDimension(getBefore()));
        for (i = 0; i < oriap.coordinates.length; i++)
            oriap.coordinates[i] = ((Pnt) curcent.elementAt(clusterori)).coordinates[i];
        Pnt finap = new Pnt(getDimension(getBefore()));
        ;
        for (i = 0; i < finap.coordinates.length; i++)
            if (taillefin > 0.0)
                finap.coordinates[i] = ((Pnt) curcent.elementAt(clusterfin)).coordinates[i];
            else
                finap.coordinates[i] = 0.0;

        center_remove(pt, oriap, tailleori);
        center_add(pt, finap, taillefin);

        A = ((double) (taillefin + 1)) * Distortion.phi(finap, distortionType, dist_param);
        if (tailleori == 1)
            B = 0.0;
        else
            B = ((double) (tailleori - 1.0)) * Distortion.phi(oriap, distortionType, dist_param);
        C = ((double) tailleori) * Distortion.phi((Pnt) curcent.elementAt(clusterori), distortionType, dist_param);
        if (taillefin > 0.0)
            D = ((double) taillefin) * Distortion.phi((Pnt) curcent.elementAt(clusterfin), distortionType, dist_param);
        else
            D = 0.0;

        oriap = finap = null;

        val = A + B - C - D;
        return val;
    }

    public void center_add(Pnt p1, Pnt curcent, int curtaille) {
        // recomputes the center if we add point p1;
        int i;
        double tmp;
        if (curcent == null) {
            if (curtaille > 0)
                Matrix.perror("AGCTClustering_KM :: The cluster should be empty");
            curcent = new Pnt(getDimension(getBefore()));
        }
        if (curcent.coordinates.length != p1.coordinates.length)
            Matrix.perror("AGCTClustering_KM :: Dimension mismatch");
        for (i = 0; i < p1.coordinates.length; i++) {
            tmp = curcent.coordinates[i];
            tmp *= ((double) curtaille);
            tmp += p1.coordinates[i];
            tmp /= ((double) (curtaille + 1));
            curcent.coordinates[i] = tmp;
        }
    }

    public void center_remove(Pnt p1, Pnt curcent, int curtaille) {
        // recomputes the center if we remove point p1;
        int i;
        double tmp;
        if ((curtaille == 0) || (curcent == null))
            Matrix.perror("AGCTClustering_KM :: The cluster should not be empty");

        if (curcent.coordinates.length != p1.coordinates.length)
            Matrix.perror("AGCTClustering_KM :: Dimension mismatch");

        if (curtaille > 1)
            for (i = 0; i < p1.coordinates.length; i++) {
                tmp = curcent.coordinates[i];
                tmp *= ((double) curtaille);
                tmp -= p1.coordinates[i];
                tmp /= ((double) (curtaille - 1));
                curcent.coordinates[i] = tmp;
            }
        else
            curcent = null;
    }

    public void KMeans_points_to_centers(int[] affect) {
        int i;
        if (affect.length != getMyClustering().getRightNumberSelectedGenes(getBefore()))
            Matrix.perror("Dimension mismatch");
        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++)
            affect[i] = closestCenterIndex(i);
    }

    public void KMeans_new_centers(int[] affect) {
        setHard_memberships_centers(new Vector<Pnt>());
        // hard_memberships_centers = new Vector();
        int i, j, k;
        double tot;
        Pnt p;
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
                for (k = 0; k < p.coordinates.length; k++) {
                    if (tot == 0.0) {
                        System.err.println("devided by zero!");
                        System.exit(0);
                    }

                    p.coordinates[k] /= tot;
                }
                getHard_memberships_centers().add(p);
            }
    }

    public double totalDistortion() {
        double val = 0.0;
        int i;
        for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++)
            val += Distortion.Bregman(getMyClustering().getRightPoint(i, getBefore()), (Pnt) getHard_memberships_centers().elementAt(closestCenterIndex(i)), distortionType, dist_param);

        return val;
    }

    public void KMeans() {
        int i, j, niter = 0;
        double dprev = -1, dcur, dcurI, ratio;
        boolean stop = false;
        boolean iterative = false;
        affect = new int[getMyClustering().getRightNumberSelectedGenes(getBefore())];

        AGCTCounter ccc;
        if (getnMaxIterations() > 0)
            ccc = new AGCTCounter("Clustering KM", getnMaxIterations());
        else
            ccc = new AGCTCounter("Clustering KM", 0);

        JInformationFrame.getInstance().setText("Iteration : decrease (%)\n");
        do {
            KMeans_points_to_centers(affect);
            KMeans_new_centers(affect);
            dcur = totalDistortion();
            setFinalClusteringObjectiveFunction(dcur);
            // finalClusteringObjectiveFunction = dcur;

            // System.out.println("KMeans distortion :: " + dcur);

            if (niter > 0) {

                if (getnMaxIterations() > 0)
                    ccc.increment();
                else
                    ccc.setPercent((int) (100.0 * AGCTClustering_Algorithm.EPSILON_KM / (dprev - dcur)));

                if (dprev < dcur)
                    Matrix.perror("Distortion increase :: " + dprev + " --> " + dcur);

                ratio = (dprev - dcur) / dprev;
                JInformationFrame.getInstance().appendText(niter + " : " + (ratio * 100.0) + "\n");
                if (ratio < AGCTClustering_Algorithm.EPSILON_KM)
                    stop = true;

				/*
				 * if (stop == false)
				 * System.out.println("Improvement on KMeans :: " + dprev +
				 * " -> " + dcur);
				 */

				/*
				 * for (i=0;i<getHard_memberships_centers().size();i++)
				 * System.out.println("Center before (" + i + ") :: " + (Pnt)
				 * getHard_memberships_centers().elementAt(i));
				 */

                // Bregman Iterative K-Means !
                if ((stop == true) && (n_IterativeKMeans > 0)) {
                    iterative = Iterative_KMeans(affect);
                    dcurI = totalDistortion();
                    setFinalClusteringObjectiveFunction(dcurI);
                    // finalClusteringObjectiveFunction = dcurI;

                    if (dcurI > dcur)
                        Matrix.perror("AGCTClustering_KM :: the cost increases from " + dcur + " to " + dcurI + " in incremental BKM");
                    // System.out.println("Distortion before :: " + dcur +
                    // ", after " + dcurI);

                    dcur = dcurI;
                }
                if (iterative)
                    stop = false; // Improvement !

				/*
				 * for (i=0;i<getHard_memberships_centers().size();i++)
				 * System.out.println("Center after (" + i + ") :: " + (Pnt)
				 * getHard_memberships_centers().elementAt(i));
				 */

            }
            dprev = dcur;
            niter++;

            if (getnMaxIterations() > 0)
                if (niter >= getnMaxIterations())
                    stop = true;

        } while (!stop);
        JInformationFrame.getInstance().setText("ok.");
        ccc.end();

        if (getHard_memberships_centers().size() < getNumberOfClusters())
            JInformationFrame.getInstance().setText("Bregman KMeans reduced the number of clusters from " + getNumberOfClusters() + " to " + getHard_memberships_centers().size());
        setNumberOfClusters(getHard_memberships_centers().size());
        // numberOfClusters = getHard_memberships_centers().size();

        if (!getBefore()) {
            for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++)
                for (j = 0; j < getNumberOfClusters(); j++)
                    if (closestCenterIndex(i) == j)
                        getSoft_memberships().set(j, i, 1.0);
                    else
                        getSoft_memberships().set(j, i, 0.0);
            setSoft_memberships_available(true);
            // soft_memberships_available = true;
            setPlotCentersAvailable(true);
            // plotCentersAvailable = true;
        }
    }

    public void toClustering() {
        if (getSeedingType() == 0)
            initHardMembershipForgy();
        else if (getSeedingType() == 1)
            if (!sameDistortion)
                initHardMembershipArthurVassilvitskii(-1, -1.0);
            else
                initHardMembershipArthurVassilvitskii(distortionType, dist_param);
        else
            Matrix.perror("Invalid seeding type");

        KMeans();
        toClusterChoices();

        if (!getBefore()) {
            toSelectedGeneToHardCluster();
            setHard_clustering_computed(true);
            // hard_clustering_computed = true;
            setSoft_clustering_computed(false);
            // soft_clustering_computed = false;
            setOk_for_statistics(true);
            // ok_for_statistics = true;
            ControlProcess.put("hardClusteringProcessed", true);
            ControlProcess.put("ClusteringProcessed_KM", true);
        }
    }

    public void toClusteringLite(Vector allPoints) {
        affect = new int[getMyAGCT().data.getMyDomain().numberSelectedGenes];
        int i, j, pind, pclu;
        Vector element;
        double dcur;

        for (i = 0; i < getMyAGCT().data.getMyDomain().numberSelectedGenes; i++) {
            element = (Vector) allPoints.elementAt(i);
            pind = ((Integer) element.elementAt(0)).intValue();
            pclu = ((Integer) element.elementAt(2)).intValue();
            affect[pind] = pclu;

            for (j = 0; j < getNumberOfClusters(); j++)
                getSoft_memberships().set(j, pind, 0.0);
            getSoft_memberships().set(pclu, pind, 1.0);
        }
        KMeans_new_centers(affect);
        dcur = totalDistortion();
        setFinalClusteringObjectiveFunction(dcur);
        // finalClusteringObjectiveFunction = dcur;
        setNumberOfClusters(getHard_memberships_centers().size());
        // numberOfClusters = getHard_memberships_centers().size();

        setSoft_memberships_available(true);
        // soft_memberships_available = true;
        setPlotCentersAvailable(true);
        // plotCentersAvailable = true;

        toClusterChoices();
        toSelectedGeneToHardCluster();
        setHard_clustering_computed(true);
        // hard_clustering_computed = true;
        setSoft_clustering_computed(false);
        // soft_clustering_computed = false;
        setOk_for_statistics(true);
        // ok_for_statistics = true;
        ControlProcess.put("hardClusteringProcessed", true);
        ControlProcess.put("ClusteringProcessed_KM", true);
    }

    public void toClusterNumberToClosestCenter(Hashtable t, Vector v) {
        int i, id;
        Hashtable bum = new Hashtable();

        AGCTCounter cc = new AGCTCounter("Prototype selection, data structure #1", getNumberOfClusters());
        for (i = 0; i < getNumberOfClusters(); i++) {
            id = closestPointIndexExclusive(t, i, bum);
            t.put(new Integer(i), new Integer(id));
            v.addElement(new Integer(id));

            bum.put(new Integer(id), new Boolean(true));
            // REMOVE
            // System.out.println(" toClusterNumberToClosestCenter: center of cluster #"
            // + i + " is gene " + id);

            cc.increment();
        }
        bum = null;
        cc.end();
    }

    public void toClosestCenterToClusterNumber(Hashtable ref, Hashtable t) {
        Enumeration extensions = ref.keys();
        Integer R, S;
        if (extensions != null) {
            AGCTCounter cc = new AGCTCounter("Prototype selection, data structure #2", ref.size());
            while (extensions.hasMoreElements()) {
                R = (Integer) extensions.nextElement();
                S = (Integer) ref.get(R);
                t.put(new Integer(S.intValue()), new Integer(R.intValue()));

                // REMOVE
                // System.out.println(" toClosestCenterToClusterNumber: gene " +
                // S + " is the center of cluster #" + R);

                cc.increment();
            }
            cc.end();
        }
    }

    public void toClosestCenterToClusterPoints(Hashtable cccn, Hashtable cncc, Hashtable t, Hashtable distortions) {
        Enumeration extensions = cccn.keys();
        Vector cSet, cDist;
        Integer I, J;
        int i;
        Pnt p, q;
        double dist;

        if (extensions != null) {
            AGCTCounter cc = new AGCTCounter("Prototype selection, data structure #3", cccn.size());
            while (extensions.hasMoreElements()) {
                I = (Integer) extensions.nextElement();

                cSet = new Vector();
                cDist = new Vector();
                t.put(I, cSet);
                distortions.put(I, cDist);
                cc.increment();
            }
            cc.end();

            cc = new AGCTCounter("Prototype selection, data structure #4", getMyClustering().getRightNumberSelectedGenes(getBefore()));
            for (i = 0; i < getMyClustering().getRightNumberSelectedGenes(getBefore()); i++) {
                I = new Integer(affect[i]);
                J = (Integer) cncc.get(I);
                cSet = (Vector) t.get(J);
                cDist = (Vector) distortions.get(J);

                cSet.addElement(new Integer(getMyClustering().getRightIndexSelectedGeneNumberToGeneNumber(i, getBefore())));

                p = getMyClustering().getRightPoint(i, getBefore());
                q = getMyClustering().getRightClosestCenterPoint(J.intValue(), getBefore());
                dist = Distortion.Bregman(p, q, distortionType, dist_param);

                cDist.addElement(new Double(dist));
                cc.increment();
            }
            cc.end();
        }
    }
}

/****
 * if(extensions != null) { AGCTCounter cc = new
 * AGCTCounter(getMyAGCT().myInformationFrame,
 * "Prototype selection, data structure #3", cccn.size()); while
 * (extensions.hasMoreElements()) { I = (Integer) extensions.nextElement(); id =
 * I.intValue(); clid = ( (Integer) cccn.get(I) ).intValue();
 *
 * cSet = new Vector(); cDist = new Vector(); t.put(I, cSet); distortions.put(I,
 * cDist); cc.increment(); } cc.end();
 *
 * cc = new AGCTCounter(getMyAGCT().myInformationFrame,
 * "Prototype selection, data structure #4",
 * getMyClustering().getRightNumberSelectedGenes(getBefore())); for
 * (i=0;i<getMyClustering().getRightNumberSelectedGenes(getBefore());i++){ I =
 * new Integer(affect[i]); J = (Integer) cncc.get(I); cSet = (Vector) t.get(J);
 * cDist = (Vector) distortions.get(J);
 *
 * cSet.addElement(new
 * Integer(getMyClustering().getRightIndexSelectedGeneNumberToGeneNumber(i,
 * getBefore())));
 *
 * p = getMyClustering().getRightPoint(i, getBefore()); q =
 * getMyClustering().getRightClosestCenterPoint(J.intValue(), getBefore()); dist
 * = Distortion.Bregman(p, q, distortionType, dist_param);
 *
 * cDist.addElement(new Double(dist)); cc.increment(); } cc.end(); }
 *****/
