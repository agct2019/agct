package clusteringProfile;

import forDebug.Debug;
import gene.Gene;
import gene.GeneList;
import test.*;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Vector;

//import javax.swing.JScrollPane;

/**
 * @author Takashi
 */
public class JClusteringProfileGridPlot extends JPanel/* JScrollPane */ implements
        Debuggable {

    /**
     *
     */
    private static final long serialVersionUID = 5517626235527757513L;
    JClusteringProfileFrame myFrame;
    AGCT myAGCT;
    private int nclusters;
    public byte percentile;

    /**
     * @return
     */
    public int getNumberOfClusters() {
        return nclusters;
    }

    /**
     *
     */
    public JClusteringProfileStampByChart[][] all_stamps;
    boolean clusterVSgenes;

    // true -> we plot cluster profiles, else gene profiles
    public JClusteringProfileGridPlot(JClusteringProfileFrame jf, AGCT ma) {
        // super(VERTICAL_SCROLLBAR_ALWAYS, HORIZONTAL_SCROLLBAR_ALWAYS);
        myFrame = jf;
        myAGCT = ma;
        all_stamps = null;
        clusterVSgenes = true;
        // setMinimumSize(new Dimension(500, 300));
    }

    @SuppressWarnings("unchecked")
    public void initAverageProfiles() {

        System.out.println("JClusteringProfileGridPlot.initAverageProfiles()");

        int np, k, l, index, cId, l25, l50, l75, l100;
        Clustering clustering;
        Domain d = myAGCT.data.getMyDomain();

        double[] ddum;

        Vector clusterElements, element;

        all_stamps = null;
        // percentile = myFrame.percentilePlots.getSelectedIndex();
        switch (myFrame.percentilePlots.getSelectedIndex()) {
            case 0:
                percentile = JClusteringProfileStamp.AVERAGE;
                break;
            case 1:
                percentile = JClusteringProfileStamp.Q25;
                break;
            case 2:
                percentile = JClusteringProfileStamp.Q50;
                break;
            case 3:
                percentile = JClusteringProfileStamp.Q75;
                break;
            case 4:
                percentile = JClusteringProfileStamp.ALL;
                break;
            default:
                break;
        }

        if ((myFrame.listPlots.getSelectedIndex() >= 0) && (myFrame.onePlot)) {
            System.out
                    .println("JClusteringProfileGridPlot.initAverageProfiles() : onePlot");

            np = (myFrame.keysForCombo.get(myFrame.listPlots
                    .getSelectedIndex())).intValue();
            clustering = (Clustering) myAGCT.data.getAllClusterings().get(np);
            int experimentId = clustering.getExperimentId();

            nclusters = clustering.getNumberOfClusters();

            int nligand = d.getLigands().numberOfSelectedLigands();
            int ntime = d.getTimes().length;
            double[][][] all_curves_average_profile = new double[nclusters][nligand][ntime];
            double[][][] all_curves_standard_deviations = new double[nclusters][nligand][ntime];
            boolean[][][] all_curves_undetermined_profile = new boolean[nclusters][nligand][ntime];

            double[][][] all_curves_Q25 = new double[nclusters][nligand][ntime];
            double[][][] all_curves_Q50 = new double[nclusters][nligand][ntime];
            double[][][] all_curves_Q75 = new double[nclusters][nligand][ntime];
            double[][][] all_curves_Q100 = new double[nclusters][nligand][ntime];
            double[][][] ex2 = new double[nclusters][nligand][ntime];
            double[][][] nb_values = new double[nclusters][nligand][ntime];

            // 変数の初期化終わり

            GeneList genes = new GeneList();
            int[] clusterIdToNumberOfGenes = new int[nclusters];
            for (int geneId = 0; geneId < d.numberSelectedGenes; geneId++) {
                Gene g = d.getGenes().get(
                        d.selectedGeneNumberToGeneNumber[geneId]);

                genes.add(g);

                index = 0;
                cId = clustering.myClusteringAlgorithm.getSelectedGeneToHardCluster()[geneId];// gの属するclusterId
                clusterIdToNumberOfGenes[cId]++;
                for (int doseId = 0; doseId < d.getLigands().size(); doseId++) {
                    if (d.getLigands().get(doseId).isChecked()) {
                        for (k = 0; k < nb_values[cId][index].length; k++) {
                            if (!g.getUndeterminedCoordinate(index)[k]) {// 値が定まっている時間
                                all_curves_average_profile[cId][index][k] += g
                                        .getRawCoordinate(index)[k];// ここでデータを入れている
                                ex2[cId][index][k] += (g
                                        .getRawCoordinate(index)[k] * g
                                        .getRawCoordinate(index)[k]);
                                nb_values[cId][index][k] += 1.0;
                            }
                        }
                        index++;
                    }
                }
            }

            for (int clusterId = 0; clusterId < nclusters; clusterId++) {
                for (int doseId = 0; doseId < all_curves_average_profile[clusterId].length; doseId++) {
                    for (k = 0; k < all_curves_average_profile[clusterId][doseId].length; k++) {
                        if (nb_values[clusterId][doseId][k] > 0.0) {
                            all_curves_average_profile[clusterId][doseId][k] /= nb_values[clusterId][doseId][k];
                            ex2[clusterId][doseId][k] /= nb_values[clusterId][doseId][k];
                            all_curves_undetermined_profile[clusterId][doseId][k] = false;
                        } else {
                            all_curves_undetermined_profile[clusterId][doseId][k] = true;
                        }
                    }
                }
            }
            for (int clusterId = 0; clusterId < nclusters; clusterId++) {
                for (int doseId = 0; doseId < all_curves_average_profile[clusterId].length; doseId++) {
                    for (k = 0; k < all_curves_average_profile[clusterId][doseId].length; k++) {
                        if (nb_values[clusterId][doseId][k] > 0.0) {
                            double average = all_curves_average_profile[clusterId][doseId][k];
                            all_curves_standard_deviations[clusterId][doseId][k] = Math
                                    .sqrt(ex2[clusterId][doseId][k]
                                            - (average * average));
                        }
                    }
                }
            }
            ArrayList<ArrayList<ArrayList<double[]>>> allplot = new ArrayList<ArrayList<ArrayList<double[]>>>();
            // Quartiles
            clusterElements = new Vector();
            for (int clusterId = 0; clusterId < nclusters; clusterId++) {
                clusterElements.addElement(new Vector());
            }
            for (int clusterId = 0; clusterId < d.numberSelectedGenes; clusterId++) {
                element = (Vector) clusterElements
                        .get(clustering.myClusteringAlgorithm
                                .getSelectedGeneToHardCluster()[clusterId]);
                element.add(new Integer(clusterId));
            }

            for (int clusterId = 0; clusterId < nclusters; clusterId++) {
                index = 0;
                element = (Vector) clusterElements.get(clusterId);
                allplot.add(new ArrayList<ArrayList<double[]>>());
                int doseId = -1;
                for (int realDoseId = 0; realDoseId < d.getLigands().size(); realDoseId++) {
                    allplot.get(clusterId).add(new ArrayList<double[]>());
                    if (d.getLigands().get(realDoseId).isChecked()) {
                        doseId++;
                        for (k = 0; k < nb_values[clusterId][index].length; k++) {

                            Vector clusterValues = new Vector();// cluster
                            // iに属するgeneのligand
                            // jに対するtime
                            // kでのactivationを入れたもの
                            for (l = 0; l < element.size(); l++) {
                                // cluster i に属するgeneがでてくる
                                Gene g = d
                                        .getGenes()
                                        .get(d.selectedGeneNumberToGeneNumber[((Integer) element
                                                .get(l)).intValue()]);
                                allplot.get(clusterId).get(realDoseId)
                                        .add(g.getRawCoordinate(index));
                                if (!g.getUndeterminedCoordinate(index)[k]) {
                                    // add activation
                                    clusterValues.addElement(new Double(g
                                            .getRawCoordinate(index)[k]));
                                }
                            }
                            ddum = new double[clusterValues.size()];

                            for (l = 0; l < clusterValues.size(); l++) {
                                ddum[l] = ((Double) clusterValues.get(l))
                                        .doubleValue();
                            }
                            QuickSort.quicksort(ddum);


                            l25 = ddum.length / 4;
                            if (l25 > 0) {
                                l25--;
                            }
                            Debug.debug("clusterId", clusterId);
                            Debug.debug("doseId", doseId);
                            Debug.debug("k", k);
                            Debug.debug("l25", l25);
                            all_curves_Q25[clusterId][doseId][k] = ddum[l25];// 小さい方から1/4

                            l50 = ddum.length / 2;
                            if (l50 > 0) {
                                l50--;
                            }
                            all_curves_Q50[clusterId][doseId][k] = ddum[l50];

                            l75 = ddum.length - (ddum.length / 4);
                            if (l75 > 0) {
                                l75--;
                            }
                            all_curves_Q75[clusterId][doseId][k] = ddum[l75];

                            l100 = ddum.length;
                            if (l100 > 0) {
                                l100--;
                            }
                            all_curves_Q100[clusterId][doseId][k] = ddum[l100];
                        }
                        index++;
                    }
                }
            }

            // 枠の大きさを決める。
//			平均値のグラフの最大・最小をとっている
            double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
            for (int clusterId = 0; clusterId < nclusters; clusterId++) {
                for (int doseId = 0; doseId < all_curves_average_profile[clusterId].length; doseId++) {
                    for (int i = 0; i < all_curves_average_profile[clusterId][doseId].length; i++) {
                        if ((all_curves_average_profile[clusterId][doseId][i] < min)) {
                            min = all_curves_average_profile[clusterId][doseId][i];
                        }
                        if ((all_curves_average_profile[clusterId][doseId][i] > max)) {
                            max = all_curves_average_profile[clusterId][doseId][i];
                        }
                    }
                }
            }


            double min75 = Double.POSITIVE_INFINITY, max75 = Double.NEGATIVE_INFINITY;
            for (int clusterId = 0; clusterId < nclusters; clusterId++) {
                for (int doseId = 0; doseId < all_curves_Q75[clusterId].length; doseId++) {
                    for (int i = 0; i < all_curves_Q75[clusterId][doseId].length; i++) {
                        if ((all_curves_Q75[clusterId][doseId][i] < min)) {
                            min75 = all_curves_Q75[clusterId][doseId][i];
                        }
                        if ((all_curves_Q75[clusterId][doseId][i] > max)) {
                            max75 = all_curves_Q75[clusterId][doseId][i];
                        }
                    }
                }
            }

            double min25 = Double.POSITIVE_INFINITY, max25 = Double.NEGATIVE_INFINITY;
            for (int clusterId = 0; clusterId < nclusters; clusterId++) {
                for (int doseId = 0; doseId < all_curves_Q25[clusterId].length; doseId++) {
                    for (int i = 0; i < all_curves_Q25[clusterId][doseId].length; i++) {
                        if ((all_curves_Q25[clusterId][doseId][i] < min)) {
                            min25 = all_curves_Q25[clusterId][doseId][i];
                        }
                        if ((all_curves_Q25[clusterId][doseId][i] > max)) {
                            max25 = all_curves_Q25[clusterId][doseId][i];
                        }
                    }
                }
            }

            double min100 = Double.POSITIVE_INFINITY, max100 = Double.NEGATIVE_INFINITY;
            for (int clusterId = 0; clusterId < nclusters; clusterId++) {
                for (int doseId = 0; doseId < all_curves_Q100[clusterId].length; doseId++) {
                    for (int i = 0; i < all_curves_Q100[clusterId][doseId].length; i++) {
                        if ((all_curves_Q100[clusterId][doseId][i] < min)) {
                            min100 = all_curves_Q100[clusterId][doseId][i];
                        }
                        if ((all_curves_Q100[clusterId][doseId][i] > max)) {
                            max100 = all_curves_Q100[clusterId][doseId][i];
                        }
                    }
                }
            }


            double ymin = 0, ymax = 0;
            switch (myFrame.percentilePlots.getSelectedIndex()) {
                case 0:
                    ymin = min;
                    ymax = max;
                    break;
                case 1:
                    ymin = 0;
                    ymax = max100;
                    break;
                case 2:
                    ymin = 0;
                    ymax = max100;

                    break;
                case 3:
                    ymin = 0;
                    ymax = max100;

                    break;
                case 4:
                    ymin = 0;
                    ymax = max100;
                    break;
                default:
                    break;
            }


            all_stamps = new JClusteringProfileStampByChart[nclusters][];// cluster
            // id,lignad
            // id ->
            // JClusteringProfileStamp
            for (int clusterId = 0; clusterId < nclusters; clusterId++) {
                all_stamps[clusterId] = new JClusteringProfileStampByChart[all_curves_average_profile[clusterId].length];
                for (int doseId = 0; doseId < all_curves_average_profile[clusterId].length; doseId++) {
                    all_stamps[clusterId][doseId] = new JClusteringProfileStampByChart(
                            clusterId,
                            doseId,
                            experimentId,

                            all_curves_average_profile[clusterId][doseId],
                            all_curves_standard_deviations[clusterId][doseId],
                            all_curves_undetermined_profile[clusterId][doseId],
                            all_curves_Q25[clusterId][doseId],
                            all_curves_Q50[clusterId][doseId],
                            all_curves_Q75[clusterId][doseId], allplot,
                            genes, ymin, ymax, percentile, clusterIdToNumberOfGenes[clusterId]);
                }
            }
        }
    }

    /**
     * 　ここと、下のメソッドが再描画のためのメソッド
     */
    public void makeUpdateLayout() {
        System.out.println("makeUpdateLayout");
        int i, j;

        removeAll();
        // JPanel jp = new JPanel();
        // FIXME 下、jpがついてなかったかもしれない
        /* jp. */
        setLayout(new GridLayout(nclusters, myAGCT.data.getMyDomain()
                .getLigands().numberOfSelectedLigands()));

        for (i = 0; i < nclusters; i++) {
            for (j = 0; j < myAGCT.data.getMyDomain().getLigands()
                    .numberOfSelectedLigands(); j++) {
                /* jp. */
                add(all_stamps[i][j]);//描画リストに追加
            }
        }

        // jp.setPreferredSize(new Dimension(200,200));

		/* set the size with the dimension of clusters */
        setPreferredSize(new Dimension(72 * myAGCT.data.getMyDomain()
                .getLigands().numberOfSelectedLigands(), 128 * nclusters));

        // this.setViewportView(jp);
        // jp.setVisible(true);
        // jpの内容が伝わらないのでは？

        repaintEverything();
        // repaint();
    }

    /**
     *
     */
    public void repaintEverything() {
        System.out.println("repaintEverything");
        int i, j;
        for (i = 0; i < nclusters; i++) {
            for (j = 0; j < myAGCT.data.getMyDomain().getLigands()
                    .numberOfSelectedLigands(); j++) {
                all_stamps[i][j].repaint();
                revalidate();
            }
        }
//		repaint();
    }

    /**
     * @param g
     */
    public void paintComponent(Graphics g) {
        if ((!clusterVSgenes)
                && ((myFrame.getListGenes() == null) || (myFrame.getListGenes()
                .size() < 1))) {
            removeAll();
            myFrame.updateBorder("");
            g.drawString("No visualisation possible: no gene buffered", 40, 40);
        } else if ((!clusterVSgenes)
                || ((myFrame.onePlot) && (myFrame.plotAuthorized))) {
//			removeAll();

            System.out.println("JClusteringProfileGridPlot # paintComponent");
            setBackground(Color.white);

            if (myFrame.moving) {
                makeUpdateLayout();
                myFrame.moving = false;
            }

            repaintEverything();
        }
    }

    /**
     * @return
     */
    public Rectangle getCaptureRectangle() {
        System.out.println("getCaptureRectangle");
        Rectangle bounds = getBounds();
        bounds.setLocation(getLocationOnScreen());
        return bounds;
    }
}
