package test;

import forDebug.Debug;
import gene.Gene;
import gui.JConfig;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.Vector;

class KeyChange extends KeyAdapter {
    public static int DELTA = 50;
    public static double FACTOR = 2.0;

    JClusteringFrame myg;

    KeyChange(JClusteringFrame g) {
        myg = g;
    }

    public void keyPressed(KeyEvent e) {

        if (myg.graphPanel.zoomAuthorized) {

            switch (e.getKeyCode()) {
                case KeyEvent.VK_ENTER:
                    myg.graphPanel.magnification = 1.0;
                    myg.graphPanel.Xoffset = 0;
                    myg.graphPanel.Yoffset = 0;
                    break;
                case KeyEvent.VK_LEFT:
                    myg.graphPanel.Xoffset += KeyChange.DELTA;
                    break;
                case KeyEvent.VK_RIGHT:
                    myg.graphPanel.Xoffset -= KeyChange.DELTA;
                    break;
                case KeyEvent.VK_DOWN:
                    myg.graphPanel.Yoffset -= KeyChange.DELTA;
                    break;
                case KeyEvent.VK_UP:
                    myg.graphPanel.Yoffset += KeyChange.DELTA;
                    break;
                case KeyEvent.VK_P:
                    myg.graphPanel.magnification *= KeyChange.FACTOR;
                    break;
                case KeyEvent.VK_M:
                    myg.graphPanel.magnification /= KeyChange.FACTOR;
                    break;
            }


            if (myg.listPlots.getSelectedIndex() >= 0) {
                myg.onePlot = true;
                myg.graphPanel.update = true;
                myg.graphPanel.repaint();
            }
        }
    }
}


class JClusteringPlot extends JPanel implements MouseMotionListener, MouseListener, Debuggable {

    public static Color[] colors = {Color.blue, Color.red, Color.pink, Color.cyan,
            Color.black, Color.yellow};

    public static Color bgColor = Color.white;
    public static Color lineColor = Color.red;
    public static Color structColor = Color.green;
    public static Color pointColor = Color.blue;
    public static Color textColor = Color.pink;

    public static Color highlightColor = Color.magenta;

    public static int pointRadius = 3;
    public static int highlightRadius = 10;

    public static double leftMarginPercent = 10.0,
            rightMarginPercent = 10.0,
            upMarginPercent = 10.0,
            downMarginPercent = 10.0;

    JClusteringFrame myFrame;
    AGCT myAGCT;
    boolean update;

    double magnification;
    boolean zoomAuthorized;
    int Xoffset, Yoffset, pMouseX, pMouseY, cMouseX, cMouseY;

    JClusteringPlot(JClusteringFrame jf, AGCT ma) {
        myFrame = jf;
        myAGCT = ma;
        update = zoomAuthorized = false;
        magnification = 1.0;
        Xoffset = Yoffset = 0;
    }

    public void mousePressed(MouseEvent e) {
        pMouseX = e.getX();
        pMouseY = e.getY();
        requestFocus();
    }

    public void mouseEntered(MouseEvent e) {
        pMouseX = e.getX();
        pMouseY = e.getY();
    }

    public void mouseExited(MouseEvent e) {
    }

    public void mouseClicked(MouseEvent e) {
    }

    public void mouseReleased(MouseEvent e) {
    }


    public void mouseDragged(MouseEvent e) {
        cMouseX = e.getX();
        cMouseY = e.getY();

        int deltaX = (cMouseX - pMouseX);
        int deltaY = (cMouseY - pMouseY);

        Xoffset -= deltaX;
        Yoffset -= deltaY;

        pMouseX = cMouseX;
        pMouseY = cMouseY;

        if (myFrame.listPlots.getSelectedIndex() >= 0) {
            myFrame.onePlot = true;
            update = true;
            repaint();
        }
    }

    public void mouseMoved(MouseEvent e) {
    }

    public Rectangle getCaptureRectangle() {
        Rectangle bounds = getBounds();
        bounds.setLocation(getLocationOnScreen());
        return bounds;
    }

    /**
     * 要チェック
     * プロットのトリガ
     *
     * @param g
     */
    public void paintComponent(Graphics g) {
        if ((myFrame.onePlot) && (myFrame.plotAuthorized) && (update)) {
            plotClustering(g);
        }
        revalidate();
    }

    public int refX(AGCTClustering_Algorithm aa, int kORp) {
        int ktarg;
        if ((kORp == 1) && (aa.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP))) {
            ktarg = ((AGCTClustering_AP) aa).valueP;
        } else {
            ktarg = aa.getNumberOfClusters();
        }
        return ktarg;
    }

    /**
     * kORpORn = 0 -> -K numberOfClusters
     * 1 -> -P preferences of AP.
     * <p/>
     * dORs = 0 -> All_Distortions
     *
     * @param aa
     * @param kORpORn
     * @param dORs
     * @return
     */
    public double refY(AGCTClustering_Algorithm aa, int kORpORn, int dORs) {
        // kORpORn = -K(number of clusters) or -P(preferences) or Clustering Number
        // dORs = Distortion or Similarity
        if ((kORpORn == 1) && (aa.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP))) {
            return (double) aa.getNumberOfClusters();
        } else if ((kORpORn == KM_NUMBER) || (dORs == ALL_DISTORTIONS)) {
            return aa.getFinalClusteringObjectiveFunction();
        } else {
            return aa.kernelSim();
        }
    }

    /**
     * @param ss
     * @param kORpORn
     * @param dORs
     * @return
     */
    public Vector getData(String ss, int kORpORn, int dORs) {
        /*
         * Returns a vector over all clusterings ss containing in this order:
         * listx = new Integer[indices.size()]; : all X data values (k or p or number) 
         * listy = new Double[indices.size()]; : all Y data (averages)
         * lists = new Double[indices.size()]; : all Sigmas over Y
         * listnb = new Integer [indices.size()]; : the # Y values for each X
         * listminmax = new Vector(4); : xmin, xmax, ymin, ymax;
         * listValues = new Vector(); : all values as double for each X
         * 
         * Tables ordered following increasing list?x
         */

        // kORpORn = 0 for -K = X & distort = Y, 
        //           1 for -P = X & -K = Y, 
        //           2 for Clustering Number = X, distort = Y

        Vector ret = new Vector();
        Vector<Integer> indices = new Vector<Integer>();
        Vector minmax = new Vector();
        Vector allValues = null;
        int j, kcur, ktarg = -1, ni, refVal;
        boolean add;
        Integer idum;
        double dni, dni2, miny = -1.0, maxy = -1.0;

        for (int clusterId = 0; clusterId < myAGCT.data.getAllClusterings().size(); clusterId++) {
            AGCTClustering_Algorithm clusteringAlgorithm = ((Clustering) myAGCT.data.getAllClusterings().get(clusterId)).myClusteringAlgorithm;
            if (clusteringAlgorithm.getMyReferenceName().equals(ss)) {
                if (kORpORn < 2) {
                    ktarg = refX(clusteringAlgorithm, kORpORn);
                } else if (kORpORn == 2) {
                    ktarg = clusterId;
                } else {
                    Matrix.perror("JClusteringFrame.class :: bad value for kORpORn");
                }
                add = true;
                j = 0;
                if (indices.size() > 0) {
                    while ((add == true) && (j < indices.size())) {
                        kcur = ((Integer) indices.get(j)).intValue();
                        if (kcur == ktarg) {
                            add = false;
                        } else {
                            j++;
                        }
                    }
                }
                if (add) {
                    indices.add(new Integer(ktarg));

                }
            }
        }
        Integer listx[] = new Integer[indices.size()];
        Double listy[] = new Double[indices.size()];
        Double listDelta[] = new Double[indices.size()];
        Integer listnb[] = new Integer[indices.size()];
        Double sumSquared[] = new Double[indices.size()];

        for (int i = 0; i < indices.size() - 1; i++) {
            listx[i] = null;
            listy[i] = null;
            listDelta[i] = null;
            listnb[i] = null;
        }

        if (kORpORn == KM_NUMBER) {
            allValues = new Vector();
            for (int i = 0; i < indices.size(); i++) {
                AGCTClustering_Algorithm clusteringAlgorithm = ((Clustering) myAGCT.data.getAllClusterings().get(((Integer) indices.get(i)).intValue())).myClusteringAlgorithm;
                dni = refY(clusteringAlgorithm, kORpORn, dORs);
                listx[i] = new Integer(((Integer) indices.get(i)).intValue());
                listy[i] = new Double(dni);
                listDelta[i] = new Double(0.0);
                listnb[i] = new Integer(1);
                allValues.addElement(new Integer(clusteringAlgorithm.getNumberOfClusters()));

                if ((i == 0) || (dni < miny)) {
                    miny = dni;
                }

                if ((i == 0) || (dni > maxy)) {
                    maxy = dni;
                }
            }
            minmax.add(new Integer(((Integer) indices.get(0)).intValue()));
            minmax.add(new Integer(((Integer) indices.get(indices.size() - 1)).intValue()));
            minmax.add(new Double(miny));
            minmax.add(new Double(maxy));
        } else if (kORpORn < 2) {
            for (int clusterId = 0; clusterId < indices.size() - 1; clusterId++) {
                for (j = clusterId + 1; j < indices.size(); j++) {
                    if (((Integer) indices.get(clusterId)).intValue() > ((Integer) indices.get(j)).intValue()) {
                        idum = (Integer) indices.get(clusterId);
                        indices.setElementAt(indices.get(j), clusterId);
                        indices.setElementAt(idum, j);
                    }
                }
            }

            for (int clusterId = 0; clusterId < indices.size(); clusterId++) {
                listx[clusterId] = (Integer) indices.get(clusterId);
            }
            minmax.add(new Integer(listx[0].intValue()));
            minmax.add(new Integer(listx[indices.size() - 1].intValue()));

            for (j = 0; j < indices.size(); j++) {
                ktarg = ((Integer) indices.get(j)).intValue();
                for (int clusterId = 0; clusterId < myAGCT.data.getAllClusterings().size(); clusterId++) {
                    AGCTClustering_Algorithm clusteringAlgorithm = ((Clustering) myAGCT.data.getAllClusterings().get(clusterId)).myClusteringAlgorithm;
                    refVal = refX(clusteringAlgorithm, kORpORn);

                    if ((clusteringAlgorithm.getMyReferenceName().equals(ss)) && (refVal == ktarg)) {

                        if (allValues == null) {
                            allValues = new Vector();
                            allValues.addElement(new Vector());
                        } else if (allValues.size() == j) {
                            allValues.addElement(new Vector());
                        }

                        ((Vector) allValues.get(j)).addElement(new Double(refY(clusteringAlgorithm, kORpORn, dORs)));

                        if (listnb[j] == null) {
                            listnb[j] = new Integer(1);
                        } else {
                            ni = 1 + listnb[j].intValue();
                            listnb[j] = new Integer(ni);
                        }

                        if (listy[j] == null) {
                            listy[j] = new Double(refY(clusteringAlgorithm, kORpORn, dORs));
                        } else {
                            dni = refY(clusteringAlgorithm, kORpORn, dORs) + listy[j].doubleValue();
                            listy[j] = new Double(dni);
                        }

                        if (sumSquared[j] == null) {
                            sumSquared[j] = new Double(refY(clusteringAlgorithm, kORpORn, dORs) * refY(clusteringAlgorithm, kORpORn, dORs));
                        } else {
                            dni = (refY(clusteringAlgorithm, kORpORn, dORs) * refY(clusteringAlgorithm, kORpORn, dORs)) + sumSquared[j].doubleValue();
                            sumSquared[j] = new Double(dni);
                        }
                    }
                }
            }

            for (j = 0; j < indices.size(); j++) {
                dni = (listy[j].doubleValue() / (double) listnb[j].intValue());
                listy[j] = new Double(dni);

                if ((j == 0) || (dni < miny)) {
                    miny = dni;
                }

                if ((j == 0) || (dni > maxy)) {
                    maxy = dni;
                }

                dni2 = (sumSquared[j].doubleValue() / (double) listnb[j].intValue()) - (dni * dni);
                dni = Math.sqrt(dni2);
                listDelta[j] = new Double(dni);
            }
            minmax.add(new Double(miny));
            minmax.add(new Double(maxy));
        }

        ret.add(listx);
        ret.add(listy);
        ret.add(listDelta);
        ret.add(listnb);
        ret.add(minmax);
        ret.add(allValues);// ここがsetOfValuesへと流れていく。

        indices = null;
        sumSquared = null;

        return ret;
    }

    private static final int ALL_DISTORTIONS = 0;// 下に使うための値

    /**
     * dORs = 0 -> All_Distortions
     *
     * @param graphics
     * @param dORs
     */
    public void plotAllClusterings(Graphics graphics, int dORs) {
        JInformationFrame.getInstance().setText("Processing All Clusterings\nThis may take some time...");

        int[] rX = new int[2];
        double[] rY = new double[2];
        double rat;
        boolean tested;
        for (int clusterId = 0; clusterId < myAGCT.data.getAllClusterings().size(); clusterId++) {
            AGCTClustering_Algorithm aa = ((Clustering) myAGCT.data.getAllClusterings().get(clusterId)).myClusteringAlgorithm;
            if ((clusterId == 0) || (aa.getNumberOfClusters() < rX[0])) {
                rX[0] = aa.getNumberOfClusters();
            }
            if ((clusterId == 0) || (aa.getNumberOfClusters() > rX[1])) {
                rX[1] = aa.getNumberOfClusters();
            }

            // rX[0] : min{numberOfClusters}
            // rX[1] : max{numberOfClusters]

            rat = refY(aa, 0, dORs);
            if ((clusterId == 0) || (rat < rY[0])) {
                rY[0] = rat;
            }
            if ((clusterId == 0) || (rat > rY[1])) {
                rY[1] = rat;
            }

            // rY[0] : min{refY(aa,0,dORs)}
            // rY[1] : max{refY(aa,0,dORs)}
        }

        for (int clusteringAlgorithmId = 0; clusteringAlgorithmId < Clustering.clusteringAlgorithms.length; clusteringAlgorithmId++) {
            tested = false;
            for (int clusterId = 0; clusterId < myAGCT.data.getAllClusterings().size(); clusterId++) {
                AGCTClustering_Algorithm aa = ((Clustering) myAGCT.data.getAllClusterings().get(clusterId)).myClusteringAlgorithm;
                if ((!tested) && (aa.getMyReferenceName().equals(Clustering.clusteringAlgorithms[clusteringAlgorithmId]))) {
                    tested = true;
                    plotClustering(graphics, aa.getMyReferenceName(),
                            JClusteringPlot.colors[clusteringAlgorithmId],
                            false, // verticalLines 
                            rX, // referencesX
                            rY, // referencesY
                            true, // useReferences
                            false, // makeStates
                            dORs);
                }
            }
        }
        if (dORs == ALL_DISTORTIONS) {
            myFrame.updateBorder(JClusteringFrame.String_All_Distortions);
        } else {
            myFrame.updateBorder(JClusteringFrame.String_All_Similarities);
        }
        JInformationFrame.getInstance().appendText(" done.");
    }

    public void initPanel(Graphics g) {
        g.setColor(bgColor);
        super.paintComponent(g);
        g.fillRect(0, 0, getWidth(), getHeight());
    }

    private static final int AP_KP = 1;
    private static final int KM_NUMBER = 2;

    public void plotClustering(Graphics g, String kfc, Color color, boolean verticalLines, int[] referencesX, double[] referencesY, boolean useReferences, boolean makeStats, int dORs) {
        int i, j, k, xx, yy, xx2, yy2, diff, si, sj, xmin, xmax, numberOfClusterings, kORpORn;
        boolean keepIt;
        double dymin = -1, dymax = -1, delta, prob, dii, djj;
        String sX, sY;
        double listy[], listDelta[], datai[], dataj[];
        int listx[], listnb[];
        Vector data, setOfValues;

        if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM)) {
            diff = 3 * pointRadius;
            myFrame.updateBorder((String) myFrame.listPlots.getSelectedItem());
            // myFrame.updateBorder(JClusteringFrame.String_KM);
        } else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM)) {
            diff = 3 * pointRadius;
            myFrame.updateBorder(JClusteringFrame.String_EM);
        } else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM)) {
            diff = 3 * pointRadius;
            myFrame.updateBorder(JClusteringFrame.String_NM);
        } else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP)) {
            diff = 3 * pointRadius;
            myFrame.updateBorder(JClusteringFrame.String_CP);
        } else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)) {
            diff = 3 * pointRadius;
            myFrame.updateBorder(JClusteringFrame.String_HC_HC);
        } else if (((String) myFrame.listPlots.getSelectedItem()).equals(JClusteringFrame.String_AP_KP)) {
            diff = 3 * pointRadius;
            myFrame.updateBorder(JClusteringFrame.String_AP_KP);
        } else {
            diff = 3 * pointRadius;
            myFrame.updateBorder(JClusteringFrame.String_AP);
        }

        kORpORn = 0;
        if ((kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP)) && (((String) myFrame.listPlots.getSelectedItem()).equals(JClusteringFrame.String_AP_KP))) {
            kORpORn = AP_KP;
        } else if (((String) myFrame.listPlots.getSelectedItem()).equals(JClusteringFrame.String_KM_Number)) {
            kORpORn = KM_NUMBER;
        }

        if (kORpORn == 0) {
            sX = "k = ";
            sY = "m = ";
        } else if (kORpORn == AP_KP) {
            sX = "p = ";
            sY = "k = ";
        } else {
            sX = "#";
            sY = "d = ";
        }

        data = getData(kfc, kORpORn, dORs);
        numberOfClusterings = ((Integer[]) data.get(0)).length;
        listx = new int[numberOfClusterings];// ソートされている
        listy = new double[numberOfClusterings];
        listDelta = new double[numberOfClusterings];// 誤差棒の範囲
        listnb = new int[numberOfClusterings];// xに対応するデータの個数

        for (int clusteringId = 0; clusteringId < numberOfClusterings; clusteringId++) {
            // System.out.print(" " + ((Integer []) data.get(0))[i]);

            listx[clusteringId] = ((Integer[]) data.get(0))[clusteringId].intValue();
            listy[clusteringId] = ((Double[]) data.get(1))[clusteringId].doubleValue();
            listDelta[clusteringId] = ((Double[]) data.get(2))[clusteringId].doubleValue();
            listnb[clusteringId] = ((Integer[]) data.get(3))[clusteringId].intValue();
        }


        if (useReferences) {
            xmin = referencesX[0];
            xmax = referencesX[referencesX.length - 1];

            dymin = referencesY[0];
            dymax = referencesY[referencesY.length - 1];
        } else {
            xmin = ((Integer) ((Vector) data.get(4)).get(0)).intValue();
            xmax = ((Integer) ((Vector) data.get(4)).get(1)).intValue();

            dymin = ((Double) ((Vector) data.get(4)).get(2)).doubleValue();
            dymax = ((Double) ((Vector) data.get(4)).get(3)).doubleValue();
        }
        //このうえ、したのすべての値がおかしい。（1通りしかない）

        setOfValues = (Vector) data.get(5);
        data = null;
        Debug.debug("##### debug Of JClusteringFrame.getData() #####");
        Debug.debug("listx", listx);
        Debug.debug("listy", listy);// 平均値
        Debug.debug("lists", listDelta);
        Debug.debug("listnb", listnb);
        Debug.debug("xmin", "xmax", "dymin", "dymax", xmin, xmax, dymin, dymax);
        Debug.debug("setOfValues", setOfValues);// 重複を含めての値

        //[[0.6127513176977166], [0.5232526906030248], [0.46689223673069963], [0.41167072305622776, 0.4236546730371733], [0.36890672080470116], [0.3879329659800131]];
                
                
        
    /* 出力例
         * [##### debug Of JClusteringFrame.getData()]
        [listx, [5, 6, 7, 8, 9, 10]]
        [listy, [0.6127513176977166, 0.5232526906030248, 0.46689223673069963, 0.41766269804670053, 0.36890672080470116, 0.3879329659800131]]
        [lists, [0.0, 0.0, 0.0, 0.00599197499047213, 0.0, 0.0]]
        [listnb, [1, 1, 1, 2, 1, 1]]
        [xmin, xmax, dymin, dymax, 5, 21, 0.18482774254360607, 0.6127513176977166]
        [setOfValues, [[0.6127513176977166], [0.5232526906030248], [0.46689223673069963], [0.41167072305622776, 0.4236546730371733], [0.36890672080470116], [0.3879329659800131]]]
         */

        // Structure

        // 緑色でx,y軸を書く -> 黒に変更tk

        g.setColor(structColor);
        g.setFont(JConfig.SMALL_FONT);

        xx = (int) X(((double) xmin), ((double) xmin), ((double) xmax));
        yy = (int) Y(dymin, dymin, dymax);
        xx2 = (int) X(((double) xmax), ((double) xmin), ((double) xmax));
        yy2 = (int) Y(dymin, dymin, dymax);
        g.drawLine(xx, yy, xx2, yy2);


        xx2 = (int) X(((double) xmin), ((double) xmin), ((double) xmax));
        yy2 = (int) Y(dymax, dymin, dymax);
        g.drawLine(xx, yy, xx2, yy2);
// xy軸終わり

        for (int clusteringId = 0; clusteringId < numberOfClusterings; clusteringId++) {
            xx = (int) X(((double) listx[clusteringId]), ((double) xmin), ((double) xmax));
            yy = (int) Y(dymin, dymin, dymax);

            g.fillOval(xx - pointRadius,
                    yy - pointRadius, pointRadius + pointRadius, pointRadius + pointRadius);
            if (listx[clusteringId] == xmin) {
                g.drawString(sX + listx[clusteringId], xx - diff, yy + 2 * diff);
            } else {
                g.drawString("" + listx[clusteringId], xx - diff, yy + 2 * diff);
            }

            if (verticalLines) {
                xx2 = (int) X(((double) listx[clusteringId]), ((double) xmin), ((double) xmax));
                yy2 = (int) Y(listy[clusteringId], dymin, dymax);
                g.drawLine(xx, yy, xx2, yy2);
            }
        }

        // Data

        for (i = 0; i < numberOfClusterings; i++) {

            g.setColor(color);
            if (listDelta[i] < 0.0) {
                Matrix.perror("Negative std dev.");
            }

            if (listDelta[i] > 0.0) {
                delta = 2 * listDelta[i];

                xx = (int) X(((double) listx[i]), ((double) xmin), ((double) xmax));
                yy = (int) Y(listy[i] - delta, dymin, dymax);

                xx2 = (int) X(((double) listx[i]), ((double) xmin), ((double) xmax));
                yy2 = (int) Y(listy[i] + delta, dymin, dymax);

                g.drawLine(xx, yy, xx2, yy2);
                g.drawLine(xx - pointRadius, yy, xx + pointRadius, yy);
                g.drawLine(xx2 - pointRadius, yy2, xx2 + pointRadius, yy2);

            }
            xx = (int) X(((double) listx[i]), ((double) xmin), ((double) xmax));
            yy = (int) Y(listy[i], dymin, dymax);

            g.fillOval(xx - pointRadius, yy - pointRadius, pointRadius + pointRadius, pointRadius + pointRadius);

            if ((i == 0) && (useReferences)) {
                g.setFont(JConfig.HIGHLIGHT_FONT);
                g.drawString(kfc, xx, yy - diff);
            }


//            FIXME このへんでデザインを変更する。
            if (!useReferences) {
                g.setColor(textColor);
            }
            g.setFont(JConfig.SMALL_FONT);

            if ((kORpORn == 0) || (kORpORn == 2)) {
                if (listnb[i] > 1) {
                    g.drawString(sY + DF.format(listy[i]) + "¥n ; s = " + DF.format(listDelta[i]) + "¥n (n = " + listnb[i] + ")", xx, yy + 1);
                } else if (kORpORn == 0) {
                    g.drawString(sY + DF.format(listy[i]) + " (n = " + listnb[i] + ")", xx, yy + 1);
                } else {
                    g.drawString(sY + DF.format(listy[i]) + " (k = " + setOfValues.get(i) + ")", xx, yy + 1);
                }
            } else {
                g.drawString(sY + DF.format(listy[i]), (int) X((double) xmin, (double) xmin, (double) xmax) + diff, yy + 1);
            }

            // g.drawString("" + DF.format(dyy) + "(n = " + (int) ((Double)
            // nbValues.get(i)).doubleValue() + ")", xx, yy);
        }

        if (!useReferences) {
            g.setColor(lineColor);
        }
        for (i = 0; i < numberOfClusterings - 1; i++) {
            xx = (int) X(((double) listx[i]), ((double) xmin), ((double) xmax));
            yy = (int) Y(listy[i], dymin, dymax);

            xx2 = (int) X(((double) listx[i + 1]), ((double) xmin), ((double) xmax));
            yy2 = (int) Y(listy[i + 1], dymin, dymax);

            g.drawLine(xx, yy, xx2, yy2);
        }

        // statistical tests

        if ((kORpORn == 0) && (makeStats)) {
            for (i = 0; i < numberOfClusterings; i++) {
                dii = listy[i];
                keepIt = true;
                j = 0;

                while ((keepIt) && (j < i)) {
                    djj = listy[j];
                    if (j != i) {
                        si = listnb[i];
                        sj = listnb[j];

                        datai = new double[si];
                        dataj = new double[sj];

                        for (k = 0; k < si; k++) {
                            datai[k] = ((Double) ((Vector) setOfValues.get(i)).get(k)).doubleValue();
                        }
                        for (k = 0; k < sj; k++) {
                            dataj[k] = ((Double) ((Vector) setOfValues.get(j)).get(k)).doubleValue();
                        }
                        prob = Statistics.tutest(datai, dataj);

                        // System.out.println("i = " + i + ", dii = " + dii +
                        // ", j = " + j + ", djj = " + djj + ", prob = " +
                        // prob);

                        if ((kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM)) || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP))
                                || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP)) || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC))
                                || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM))) {
                            if ((djj < dii) || (prob > Statistics.LIMIT_P_CHI2)) {
                                keepIt = false;
                            }
                        }

                        if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM)) {
                            if ((djj > dii) || (prob > Statistics.LIMIT_P_CHI2)) {
                                keepIt = false;
                            }
                        }
                    }
                    j++;
                }

                if (listnb[i] == 1.0) {
                    keepIt = false;
                }

                if (keepIt) {
                    g.setColor(highlightColor);
                    xx = (int) X(((double) listx[i]), ((double) xmin), ((double) xmax));
                    yy = (int) Y(listy[i], dymin, dymax);

                    g.drawArc(xx - highlightRadius,
                            yy - highlightRadius,
                            highlightRadius + highlightRadius,
                            highlightRadius + highlightRadius, 0, 360);
                }
            }
        }

        data = setOfValues = null;
    }


    /**
     * gが何か調べています。
     *
     * @param g
     */
    public void plotClustering(Graphics g) {
        initPanel(g);

        int ind;
        String kfc = "", hcs = "";
        boolean ok = true;

        try {
            ind = myFrame.listPlots.getSelectedIndex();
            kfc = (String) myFrame.keysForCombo.get(ind);
            hcs = (String) myFrame.listPlots.getSelectedItem();
        } catch (NullPointerException e) {
            ok = false;
        }

        if (ok) {
            if ((zoomAuthorized == false) && (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)))
                zoomAuthorized = true;
            if ((zoomAuthorized == true) && (!kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)))
                zoomAuthorized = false;

            if ((kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM))
                    || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM))
                    || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP))
                    || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP))
                    || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM))) {
                plotClustering(g, kfc, pointColor, true, null, null, false, true, 0);
            } else if (hcs.equals(JClusteringFrame.String_HC_HC)) {
                plotClustering(g, kfc, pointColor, true, null, null, false, false, 0);
            } else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)) {
                plotClusteringHC(g);
            } else if ((kfc.equals(AGCTClustering_Algorithm.ALL_CLUSTERS)) && (hcs.equals(JClusteringFrame.String_All_Distortions))) {
                plotAllClusterings(g, ALL_DISTORTIONS);
            } else if ((kfc.equals(AGCTClustering_Algorithm.ALL_CLUSTERS)) && (hcs.equals(JClusteringFrame.String_All_Similarities))) {
                plotAllClusterings(g, 1);
            }
            update = false;
        }
    }

    /**
     * @param g
     */
    public void plotClusteringHC(Graphics g) {
        myFrame.updateBorder(JClusteringFrame.String_HC);

        int ind = myFrame.listPlots.getSelectedIndex(),
                indhc = ((Integer) myFrame.indexHC.elementAt(ind)).intValue(),
                i, j, k, g1, g2, tper, cli, clj, minp = -1;

        double xx, yy, oxx, oyy, xx1, oxx1, yy1, oyy1, xx2, oxx2, yy2, oyy2, nxx, onxx,
                minDisplayX = 0.0, maxDisplayX = (double) getWidth(), minDisplayY = 0.0, maxDisplayY = (double) getHeight();

        double dxmax = ((AGCTClustering_HC) ((Clustering) myAGCT.data.getAllClusterings().elementAt(indhc)).myClusteringAlgorithm).totdistort;
        double delta = dxmax * 0.10;
        dxmax += delta;

        boolean first = true;

        int nclust = ((Clustering) myAGCT.data.getAllClusterings().elementAt(indhc)).getNumberOfClusters();

        double dxmin = 0.0;
        double dymax = (double) myAGCT.data.getMyDomain().numberSelectedGenes;
        double dymin = 0.0;

        double dxx, dyy, dxx1, dyy1, dxx2, dyy2;

        Vector allX = new Vector(), allY = new Vector();
        Vector dummyHierarchy = ((AGCTClustering_HC) ((Clustering) myAGCT.data.getAllClusterings().elementAt(indhc)).myClusteringAlgorithm).dummyHierarchy;
        Vector nonDummyHierarchy = ((AGCTClustering_HC) ((Clustering) myAGCT.data.getAllClusterings().elementAt(indhc)).myClusteringAlgorithm).nonDummyHierarchy;
        Cluster ccur, ci, cj;
        Vector v1, v2;

        int[] place = new int[myAGCT.data.getMyDomain().numberSelectedGenes];

        Matrix sf = ((AGCTClustering_HC) ((Clustering) myAGCT.data.getAllClusterings().elementAt(indhc)).myClusteringAlgorithm).soft_memberships;
        Gene gg;
        boolean isTop;

        for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++)
            place[i] = i;

        for (i = dummyHierarchy.size() - 1; i >= 0; i--) {
            ccur = (Cluster) dummyHierarchy.elementAt(i);
            if ((ccur.indexPrevious1 != -1) || (ccur.indexPrevious2 != -1)) {
                v1 = ((Cluster) dummyHierarchy.elementAt(ccur.indexPrevious1)).allGenes;
                v2 = ((Cluster) dummyHierarchy.elementAt(ccur.indexPrevious2)).allGenes;

                for (j = 0; j < v1.size(); j++) {
                    g1 = ((Integer) v1.elementAt(j)).intValue();

                    for (k = 0; k < v2.size(); k++) {
                        g2 = ((Integer) v2.elementAt(k)).intValue();
                        if ((k == 0) || (place[minp] > place[g2]))
                            minp = g2;
                    }

                    if (place[g1] > place[minp]) {
                        tper = place[g1];
                        place[g1] = place[minp];
                        place[minp] = tper;
                    }
                }
            }
        }

        for (i = 0; i < nonDummyHierarchy.size(); i++) {
            ccur = (Cluster) nonDummyHierarchy.elementAt(i);
            isTop = ccur.top;
            dxx = dxmax - ccur.curdist;

            if ((ccur.indexPrevious1 == -1) || (ccur.indexPrevious2 == -1)) {
                dyy = (double) place[i];
            } else {
                dyy = ((Double) allY.elementAt(ccur.indexPrevious1)).doubleValue()
                        + ((Double) allY.elementAt(ccur.indexPrevious2)).doubleValue();
                dyy /= 2.0;

                g.setColor(structColor);

                dxx1 = ((Double) allX.elementAt(ccur.indexPrevious1)).doubleValue();
                dyy1 = ((Double) allY.elementAt(ccur.indexPrevious1)).doubleValue();
                xx = X(dxx, dxmin, dxmax);
                xx1 = X(dxx1, dxmin, dxmax);
                yy1 = Y(dyy1, dymin, dymax);

                oxx = scaledX(xx, magnification, Xoffset);
                oxx1 = scaledX(xx1, magnification, Xoffset);
                oyy1 = scaledY(yy1, magnification, Yoffset);

                g.setColor(((Clustering) myAGCT.data.getAllClusterings().elementAt(indhc)).colorCluster[ccur.clusterIndex]);

                // horizontal line down
                g.drawLine((int) oxx, (int) oyy1, (int) oxx1, (int) oyy1);

                dxx2 = ((Double) allX.elementAt(ccur.indexPrevious2)).doubleValue();
                dyy2 = ((Double) allY.elementAt(ccur.indexPrevious2)).doubleValue();
                xx2 = X(dxx2, dxmin, dxmax);
                yy2 = Y(dyy2, dymin, dymax);

                oxx2 = scaledX(xx2, magnification, Xoffset);
                oyy2 = scaledY(yy2, magnification, Yoffset);

                // horizontal line up
                g.drawLine((int) oxx, (int) oyy2, (int) oxx2, (int) oyy2);

                // vertical line
                g.drawLine((int) oxx, (int) oyy1, (int) oxx, (int) oyy2);

                if ((first) || (oxx < minDisplayX))
                    minDisplayX = oxx;
                if ((first) || (oxx > maxDisplayX))
                    maxDisplayX = oxx;

                if ((first) || (oxx1 < minDisplayX))
                    minDisplayX = oxx1;
                if ((first) || (oxx1 > maxDisplayX))
                    maxDisplayX = oxx1;

                if ((first) || (oxx2 < minDisplayX))
                    minDisplayX = oxx2;
                if ((first) || (oxx2 > maxDisplayX))
                    maxDisplayX = oxx2;

                if ((first) || (oyy1 < minDisplayY))
                    minDisplayY = oyy1;
                if ((first) || (oyy1 > maxDisplayY))
                    maxDisplayY = oyy1;

                if ((first) || (oyy2 < minDisplayY))
                    minDisplayY = oyy2;
                if ((first) || (oyy2 > maxDisplayY))
                    maxDisplayY = oyy2;
            }

            allX.add(new Double(dxx));
            allY.add(new Double(dyy));

            xx = X(dxx, dxmin, dxmax);
            yy = Y(dyy, dymin, dymax);

            oxx = scaledX(xx, magnification, Xoffset);
            oyy = scaledY(yy, magnification, Yoffset);

            if ((first) || (oxx < minDisplayX))
                minDisplayX = oxx;
            if ((first) || (oxx > maxDisplayX))
                maxDisplayX = oxx;

            if ((first) || (oyy < minDisplayY))
                minDisplayY = oyy;
            if ((first) || (oyy > maxDisplayY))
                maxDisplayY = oyy;

            if (isTop) {
                g.setColor(((Clustering) myAGCT.data.getAllClusterings().elementAt(indhc)).colorCluster[ccur.clusterIndex]);
                g.setFont(JConfig.SMALL_FONT);
                nxx = X(dxmin, dxmin, dxmax);

                onxx = scaledX(nxx, magnification, Xoffset);

                if ((first) || (onxx < minDisplayX))
                    minDisplayX = onxx;
                if ((first) || (onxx > maxDisplayX))
                    maxDisplayX = onxx;

                g.drawLine((int) onxx, (int) oyy, (int) oxx, (int) oyy);
                g.drawString(Clustering.Center_Name + ccur.clusterIndex, (int) onxx, (int) oyy - (pointRadius));
            }

            g.setColor(((Clustering) myAGCT.data.getAllClusterings().elementAt(indhc)).colorCluster[ccur.clusterIndex]);

            //g.setColor(pointColor);
            //g.fillOval(xx - pointRadius,
            //	       yy - pointRadius,
            //       pointRadius + pointRadius,
            //       pointRadius + pointRadius);

            if ((ccur.indexPrevious1 == -1) || (ccur.indexPrevious2 == -1)) {
                gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);

                g.setFont(JConfig.SMALL_FONT);
//                Debug.debug("Drawing the gene name " + gg.myID + "(" + gg.asciiName + ")");
//                g.drawString(gg.myID + "(" + gg.asciiName + ")", (int) Math.round(xx), (int) Math.round(yy + pointRadius));
                g.drawString(gg.myID + "(" + gg.asciiName + ")", (int)scaledX(xx, magnification, Xoffset), (int)scaledY(yy + pointRadius,magnification, Yoffset));
            }

            first = false;
        }

        int dX = getWidth() / 20;
        int dY = getHeight() / 20;
        int lim = 5;

        double ratX = (double) getWidth() / ((double) maxDisplayX - (double) minDisplayX);
        if (ratX > 1.0)
            ratX = 1.0;
        double ratY = (double) getHeight() / ((double) maxDisplayY - (double) minDisplayY);
        if (ratY > 1.0)
            ratY = 1.0;

        ratX *= (double) dX;
        ratY *= (double) dY;

        if (ratX < lim)
            ratX = lim;

        if (ratY < lim)
            ratY = lim;

        double off = 15.0;

        int pX = (int) (((-minDisplayX + ((double) getWidth() / 2.0)) / ((double) maxDisplayX - (double) minDisplayX)) * (double) dX);

        int pY = (int) (((-minDisplayY + ((double) getHeight() / 2.0)) / ((double) maxDisplayY - (double) minDisplayY)) * (double) dY);

        double mX, mY, MX, MY;
        mX = (double) pX - (ratX / 2.0);
        MX = (double) pX + (ratX / 2.0);
        mY = (double) pY - (ratY / 2.0);
        MY = (double) pY + (ratY / 2.0);

        if (mX < 0.0)
            mX = 0.0;
        if (mX > dX)
            mX = dX;
        if (MX < 0.0)
            MX = 0.0;
        if (MX > dX)
            MX = dX;

        if (mY < 0.0)
            mY = 0.0;
        if (mY > dY)
            mY = dY;
        if (MY < 0.0)
            MY = 0.0;
        if (MY > dY)
            MY = dY;

        Graphics2D gg2 = (Graphics2D) g;

        g.setColor(Color.black);

        gg2.fill(new java.awt.geom.Rectangle2D.Double(off + 0.0, off + 0.0,
                (double) dX,
                (double) dY));

        g.setColor(Color.white);

        gg2.fill(new java.awt.geom.Rectangle2D.Double(off + mX, off + mY, MX - mX, MY - mY));


        g.setColor(Color.black);

        gg2.draw(new java.awt.geom.Rectangle2D.Double(off + 0.0, off + 0.0,
                (double) dX,
                (double) dY));

        if (pX < 0)
            pX = 0;

        if (pY < 0)
            pY = 0;

        if (pX > dX)
            pX = dX;

        if (pY > dY)
            pY = dY;

	/*g.fillOval((int) off + pX - pointRadius,
           (int) off + pY - pointRadius,
	           pointRadius + pointRadius,
	           pointRadius + pointRadius);*/
    }

    public double scaledX(double xx, double mag, int Xoffset) {
        return (((((double) xx - ((double) getWidth() / 2.0)) * mag) + (double) getWidth() / 2.0) - (double) Xoffset);
    }

    public double scaledY(double yy, double mag, int Yoffset) {
        return (((((double) yy - ((double) getHeight() / 2.0)) * mag) + (double) getHeight() / 2.0) - (double) Yoffset);
    }

    public double X(double xx, double xmin, double xmax) {
        double w = ((double) getWidth() * (100.0 - leftMarginPercent - rightMarginPercent) / 100.0);
        double delta = ((double) getWidth() * leftMarginPercent) / 100.0;
        double v = delta + (((xx - xmin) * w) / (xmax - xmin));
        return v;
    }

    public double Y(double yy, double ymin, double ymax) {
        double h = ((double) getHeight() * (100.0 - upMarginPercent - downMarginPercent) / 100);
        double delta = ((double) getHeight() * upMarginPercent) / 100.0;
        double v = h + delta - (((yy - ymin) * h) / (ymax - ymin));
        return v;
    }
}

class JClusteringFrame extends JFrame {

    public static String Default_Border = "(no plot selected)";
    public static String Default_Combo[] = {"(no clustering)"};


    // TODO Here what you looking for.
    public static String String_KM = "All K-Means Distortions = f(K) (regardless of -Point)";
    public static String String_KM_Number = "All K-Means Distortions = f(Clustering Number) (regardless of -Point)";
    public static String String_EM = "All EM Log-likelihoods = f(K) (regardless of -Point)";
    public static String String_CP = "All CP Distortions = f(K) (regardless of -Point)";
    public static String String_AP = "All AP Distortions = f(K) (regardless of -Point)";
    public static String String_NM = "All NM Distortions = f(K) (regardless of -Point)";
    public static String String_HC_HC = "All HC Distortions = f(K) (regardless of -Point)";
    public static String String_HC = "Dendrogram for Hierarchical Clustering ";
    public static String String_AP_KP = "K = f(P) for All AP";
    /**
     * plotAllClusterings(graphics,ALL_DISTORTIONS) がよばれる.
     */
    public static String String_All_Distortions = "All Distortions = f(K) (regardless of -Point, for all algorithms)";
    public static String String_All_Similarities = "Kernel similarities = f(K) (initial space, for all algorithms)";
    String borderString;
    JComboBox listPlots;
    JClusteringPlot graphPanel;
    JButton goButton, cameraButton;
    Box selectionBox;
    Vector keysForCombo;
    // contains the vector of reference clustering class names for the combo
    // indices
    Vector indexHC;
    // contains Integers = index of the HC in the Clusterings in AGCT
    AGCT myAGCT;
    boolean onePlot, plotAuthorized, oneAlgo;

    public void flushAll() {
        listPlots.removeAllItems();
        listPlots.addItem(JClusteringFrame.Default_Combo);
        keysForCombo = null;
        indexHC = null;
        listPlots.setEnabled(false);
        onePlot = oneAlgo = false;
        plotAuthorized = false;
    }

    public void captureAndSave() {
        Rectangle rect = graphPanel.getCaptureRectangle();
        BufferedImage screencapture = null, framecapture = null;
        try {
            framecapture = new Robot().createScreenCapture(rect);
        } catch (AWTException e) {
        }

        JFileChooser chooser = new JFileChooser();
        ExampleFileFilter filter = new ExampleFileFilter();
        filter.addExtension("png");
        filter.setDescription(".png Files");
        chooser.setFileFilter(filter);
        chooser.setApproveButtonText("Save");
        int returnVal = chooser.showSaveDialog(this);
        if ((chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION)) {
            try {
                ImageIO.write(framecapture, "png", chooser.getSelectedFile());
            } catch (IOException e) {
            }
            JInformationFrame.getInstance().setText("Saving clustering pane to file: " + chooser.getSelectedFile().getName());
            ControlProcess.put("frameCapturedClustering", true);
        }
    }

    public void displayInfo() {
        cameraButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/camera.png"))));
        cameraButton.setToolTipText("capture the visible clustering plot");
        cameraButton.setActionCommand("capture");
        cameraButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                captureAndSave();
            }
        });
        cameraButton.setEnabled(false);

        goButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/go.png"))));
        goButton.setToolTipText("run clustering");
        goButton.setActionCommand("clustering");
        goButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (listPlots.getSelectedIndex() >= 0) {
                    onePlot = true;
                    graphPanel.update = true;
                    graphPanel.repaint();
                }
            }
        });
        goButton.setEnabled(false);

        graphPanel = new JClusteringPlot(this, myAGCT);
        graphPanel.setBorder(BorderFactory.createTitledBorder(JClusteringFrame.Default_Border));

        graphPanel.addMouseListener(graphPanel);
        graphPanel.addMouseMotionListener(graphPanel);


        listPlots = new JComboBox(JClusteringFrame.Default_Combo);
        listPlots.setToolTipText("select which clustering data to plot");
        listPlots.setSelectedIndex(0);
        listPlots.setEnabled(false);

        KeyChange listener = new KeyChange(this);
        cameraButton.addKeyListener(listener);
        goButton.addKeyListener(listener);
        listPlots.addKeyListener(listener);
        graphPanel.addKeyListener(listener);

        selectionBox = Box.createHorizontalBox();
        selectionBox.add(cameraButton);
        selectionBox.add(listPlots);
        selectionBox.add(goButton);

        Container pane = getContentPane();
        pane.setLayout(new BorderLayout());
        pane.add(selectionBox, BorderLayout.NORTH);
        pane.add(graphPanel, BorderLayout.CENTER);

//		addWindowListener(new FermetureListener("Closing AGCT's ClusteringFrame\n"));
        setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        Image img = Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/information.png"));
        setIconImage(img);
    }

    public void updateBorder(String t) {
        graphPanel.setBorder(BorderFactory.createTitledBorder(t));
    }

    public void updateClusteringLF() {
        int i, nhc = 0;
        AGCTClustering_Algorithm aa;
        flushAll();

        if (myAGCT.data.getAllClusterings() != null) {
            listPlots.removeAllItems();
            if (!isVisible()) {
                setVisible(true);
                myAGCT.myFrame.toFront();
            }

            listPlots.setEnabled(true);
            keysForCombo = new Vector();
            indexHC = new Vector();
            cameraButton.setEnabled(true);
            goButton.setEnabled(true);

            for (i = 0; i < myAGCT.data.getAllClusterings().size(); i++) {
                aa = ((Clustering) myAGCT.data.getAllClusterings().get(i)).myClusteringAlgorithm;
                if (oneAlgo == false) {
                    keysForCombo.add(AGCTClustering_Algorithm.ALL_CLUSTERS);
                    listPlots.addItem(JClusteringFrame.String_All_Distortions);
                    keysForCombo.add(AGCTClustering_Algorithm.ALL_CLUSTERS);
                    listPlots.addItem(JClusteringFrame.String_All_Similarities);
                    oneAlgo = true;
                    indexHC.add(new Integer(-1));
                    indexHC.add(new Integer(-1));
                }
                if ((!aa.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)) && (!keysForCombo.contains(aa.getMyReferenceName()))) {
                    keysForCombo.add(aa.getMyReferenceName());
                    indexHC.add(new Integer(-1));
                    if (aa.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM)) {
                        listPlots.addItem(JClusteringFrame.String_KM);
                        keysForCombo.add(aa.getMyReferenceName());
                        indexHC.add(new Integer(-1));
                        listPlots.addItem(JClusteringFrame.String_KM_Number);
                    } else if (aa.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM)) {
                        listPlots.addItem(JClusteringFrame.String_EM);
                    } else if (aa.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP)) {
                        listPlots.addItem(JClusteringFrame.String_CP);
                    } else if (aa.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM)) {
                        listPlots.addItem(JClusteringFrame.String_NM);
                    } else if (aa.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP)) {
                        listPlots.addItem(JClusteringFrame.String_AP);

                        indexHC.add(new Integer(-1));
                        keysForCombo.add(aa.getMyReferenceName());
                        listPlots.addItem(JClusteringFrame.String_AP_KP);
                    }
                }
            }

            oneAlgo = false;
            for (i = 0; i < myAGCT.data.getAllClusterings().size(); i++) {
                aa = ((Clustering) myAGCT.data.getAllClusterings().get(i)).myClusteringAlgorithm;
                if (aa.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)) {
                    if (!oneAlgo) {
                        keysForCombo.add(aa.getMyReferenceName());
                        listPlots.addItem(JClusteringFrame.String_HC_HC);
                        indexHC.add(new Integer(-1));
                        oneAlgo = true;
                    }

                    keysForCombo.add(aa.getMyReferenceName());
                    listPlots.addItem(JClusteringFrame.String_HC + nhc);
                    indexHC.add(new Integer(i));
                    nhc++;
                }
            }
            plotAuthorized = true;
        }

    }

    /**
     * @param name
     * @param agct
     */
    JClusteringFrame(String name, AGCT agct) {
        super(name);
        if (AGCT.MYDEBUG) {
            AGCT.debug("JCluteringFrame(name,agct)");
        }
        setSize(AGCT.WindowWidth, 600);
        myAGCT = agct;
        displayInfo();
        setVisible(false);
        setResizable(true);
        onePlot = false;

        addComponentListener(new ComponentAdapter() {

            public void componentResized(ComponentEvent evt) {
                graphPanel.update = true;
                if (onePlot) {
                    graphPanel.repaint();
                }
            }

            public void componentMoved(ComponentEvent evt) {
                graphPanel.update = true;
                if (onePlot) {
                    graphPanel.repaint();
                }
            }
        });
    }
}