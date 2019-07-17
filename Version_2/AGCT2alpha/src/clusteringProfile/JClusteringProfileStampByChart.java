/**
 *
 */
package clusteringProfile;

import forDebug.Debug;
import gene.GeneList;
import gene.JGeneSaveFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;

/**
 * @author takashi
 *
 */

/**
 * @author takashi
 */
public class JClusteringProfileStampByChart extends JPanel {

    /**
     *
     */
    private static final long serialVersionUID = -1199420031885940667L;
    private GeneList geneList;
    final public static byte Q25 = 25;
    final public static byte Q50 = 50;
    final public static byte Q75 = 75;
    final public static byte AVERAGE = 100;
    public static final byte ALL = 125;
    private byte percentile;
    String myClusterName, myLigandName;
    /**
     * 時間の最小、最大
     */
    double xmin, xmax;
    /**
     * 平均のactivationの最大、最小
     */
    double ymin, ymax;
    boolean isNull;
    private ClusteringProfile clusteringProfile;

    static boolean errorBarEnabled;

    public static void setErrorBarEnabled(boolean enabled) {
        errorBarEnabled = enabled;
        if (renderer == null) return;
        renderer.setDrawXError(errorBarEnabled);
        renderer.setDrawYError(errorBarEnabled);
    }

    public ClusteringProfile getClusteringProfile() {
        return this.clusteringProfile;
    }

    /**
     * @param clusterId
     * @param doseId
     * @param averageProfiles
     * @param standardDeviations
     * @param undeterminedProfiles
     * @param q25
     * @param q50
     * @param q75
     * @param allplot
     * @param allGenes
     * @param ymin
     * @param ymax
     * @param percentile
     * @param numGenes
     */
    public JClusteringProfileStampByChart(int clusterId, int doseId, int experimentId, double[] averageProfiles,
                                          double[] standardDeviations, boolean[] undeterminedProfiles, double[] q25, double[] q50,
                                          double[] q75, ArrayList<ArrayList<ArrayList<double[]>>> allplot,
                                          GeneList allGenes, double ymin, double ymax, byte percentile, int numGenes) {
        Debug.debug("ByChart!! " + clusterId + " " + doseId);
        setLayout(new BorderLayout());
        this.ymin = ymin;
        this.ymax = ymax;
        this.geneList = allGenes;
        myClusterName = "C" + clusterId;
        myLigandName = JClusteringProfileFrame.Ligand_Name_Summaries[doseId];

        this.percentile = percentile;

        // 新しく追加されている命令とクラス（dataクラスとかんがえられる）
        clusteringProfile = new ClusteringProfile(clusterId, doseId, experimentId, averageProfiles, standardDeviations, undeterminedProfiles, q25, q50,
                q75, allplot, allGenes, numGenes);
        isNull = false;

        xmin = ((Double) JClusteringProfileFrame.Time_Stamps_Summaries[doseId]
                .get(0)).doubleValue();// 時間の最小
        xmax = ((Double) JClusteringProfileFrame.Time_Stamps_Summaries[doseId]
                .get(JClusteringProfileFrame.Time_Stamps_Summaries[doseId].size() - 1))
                .doubleValue();// 時間の最大
        // for (int i = 0; i < q75.length; i++) {
        // if ((i == 0) || (q25[i] < ymin))
        // ymin = genes_profile[i];
        // if ((i == 0) || (q75[i] > ymax))
        // ymax = q75[i];
        // }

        setMinimumSize(new Dimension(300, 500));
        setPreferredSize(new Dimension(300, 500));

        addMouseListener(new MyMouseListener2());

        initialize();

    }


    public void initialize() {
        int i1, i2;
        boolean stop = false;
        double[] subject_drawing = null;
        String label = null;

        switch (percentile) {// nullのままになってしまった
            case Q25:
                subject_drawing = clusteringProfile.q25_profiles;
                label = "q25";
                break;
            case Q50:
                subject_drawing = clusteringProfile.q50_profiles;
                label = "q50";
                break;
            case Q75:
                subject_drawing = clusteringProfile.q75_profiles;
                label = "Q75";
                break;
            case AVERAGE:
                subject_drawing = clusteringProfile.average_profiles;
                label = "Average";
                break;
            case ALL:
                label = "All";
                break;// switch in another place
            default:
                System.err.println("percentile is not set");
                break;
        }

        if (!isNull) {
            if (percentile == ALL) {
                allGenePlot();
            } else {
                partGenePlot(stop, subject_drawing, label);
            }
        }
    }

    static XYErrorRenderer renderer;

    private void partGenePlot(boolean stop, double[] subject_drawing, String label) {
        XYIntervalSeriesCollection dataset = createDataset(stop, subject_drawing, label);

        String chartName = createChartName(clusteringProfile);
        ///
        NumberAxis axisX = new NumberAxis("Time");
        NumberAxis axisY = new NumberAxis("Expression");
        if (renderer == null)
            renderer = new XYErrorRenderer();
        renderer.setErrorPaint(Color.RED);
        renderer.setSeriesPaint(0, Color.BLUE);
        renderer.setDrawXError(errorBarEnabled);
        renderer.setDrawYError(errorBarEnabled);
        XYLineAndShapeRenderer xyLineAndShapeRenderer
                = new XYLineAndShapeRenderer();
        xyLineAndShapeRenderer.setSeriesPaint(0, Color.BLUE);
        XYPlot plot = new XYPlot(dataset, axisX, axisY, renderer);
        plot.clearAnnotations();
        plot.setDataset(1, dataset);
        plot.setRenderer(1, xyLineAndShapeRenderer);
        plot.setBackgroundPaint(Color.lightGray);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        JFreeChart chart = new JFreeChart(chartName, plot);
        chart.removeLegend();
        chart.setBackgroundPaint(Color.WHITE);
        chart.getTitle().setFont(new Font(Font.DIALOG, Font.PLAIN, 20));
        ChartPanel panel = new ChartPanel(chart);
        add(panel, BorderLayout.CENTER);
    }

    private XYIntervalSeriesCollection createDataset(boolean stop, double[] subject_drawing, String label) {
        int i1;
        int i2;
        XYIntervalSeriesCollection dataset = new XYIntervalSeriesCollection();
        XYIntervalSeries series = new XYIntervalSeries(label);

        i1 = 0;
        // これより下、subject_drawingを導入する前はそこに全てclusteringProfile.average_profilesが入っていた。

        while ((i1 < subject_drawing.length)
                && (clusteringProfile.undetermined_profiles[i1])) {
            i1++;
        }
        do {

            i2 = i1 + 1;
            while ((i2 < subject_drawing.length)
                    && (clusteringProfile.undetermined_profiles[i2])) {
                i2++;
            }
            if (i2 < subject_drawing.length) {// 値が定まっているところだけプロットする
                double time1 = clusteringProfile.time(i1);
                double sd1 = clusteringProfile.standard_deviations[i1];
                double value1 = subject_drawing[i1];
                series.add(time1, time1, time1,
                        value1, value1 - sd1, value1 + sd1);
                double time2 = clusteringProfile.time(i2);
                double value2 = subject_drawing[i2];
                double sd2 = clusteringProfile.standard_deviations[i2];
                series.add(time2, time2, time2,
                        value2, value2 - sd2, value2 + sd2);
                i1 = i2;
            } else {
                stop = true;
            }
        } while (!stop);

        dataset.addSeries(series);
        return dataset;
    }

    private void allGenePlot() {
        boolean stop;
        int i1;
        int i2;
        XYSeriesCollection dataset = new XYSeriesCollection();

        ArrayList<XYSeries> serieses = new ArrayList<XYSeries>();
        ArrayList<XYSeriesCollection> datasetlist = new ArrayList<XYSeriesCollection>();
//				XYSeries oneseries = new XYSeries("");

        JFreeChart chart = ChartFactory.createXYLineChart("Cluster --#"
                + clusteringProfile.myClusterId + "("
                + clusteringProfile.myLigandId + ") [" + clusteringProfile.all_genes_plot().size() + "genes]", "time", "expression",
                dataset, PlotOrientation.VERTICAL, false, false, false);


        int i_cp = 0;
        for (double[] drawingTarget : clusteringProfile
                .all_genes_plot()) {
            XYSeries oneseries = new XYSeries("");
            dataset = new XYSeriesCollection();

            stop = false;
            i1 = 0;
            while ((i1 < drawingTarget.length)
                    && (clusteringProfile.undetermined_profiles[i1])) {
                i1++;
            }
            do {
                i2 = i1 + 1;
                while ((i2 < drawingTarget.length)
                        && (clusteringProfile.undetermined_profiles[i2])) {
                    i2++;
                }
                if (i2 < drawingTarget.length) {// 値が定まっているところだけプロットする
                    oneseries.add(clusteringProfile.time(i1),
                            drawingTarget[i1]);
                    oneseries.add(clusteringProfile.time(i2),
                            drawingTarget[i2]);


                    i1 = i2;
                } else {
                    stop = true;
                }
            } while (!stop);

            dataset.addSeries(oneseries);

            XYItemRenderer renderer = new StandardXYItemRenderer();
            renderer.setBaseOutlinePaint(Color.red);
            XYPlot plot = chart.getXYPlot();
            plot.setRenderer(i_cp, renderer);
            plot.setDataset(i_cp, dataset);

            i_cp++;
        }
        chart.getTitle().setFont(new Font(Font.DIALOG, Font.PLAIN, 20));
        ChartPanel cpanel = new ChartPanel(chart);
        add(cpanel, BorderLayout.CENTER);
    }

    private String createChartName(ClusteringProfile clusteringProfile) {
        int numGenes = clusteringProfile.getNumberOfGenes();
        return "Cluster "
                + clusteringProfile.myClusterId
                + "[" + numGenes + "genes]";
    }

    // 上から下に移動したtk
    // //////////// output function part /////////////
    private class MyMouseListener2 extends MouseAdapter {
        private JPopupMenu popupMenu = new JPopupMenu();

        private void init() {
            JMenuItem geneSaveMenuItem = new JMenuItem("save geneList");
            geneSaveMenuItem.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    new JGeneSaveFrame(geneList);
                }
            });
            popupMenu.add(geneSaveMenuItem);
        }

        private MyMouseListener2() {
            init();
        }

        private void showPopup(MouseEvent e) {
            if (e.isPopupTrigger()) {
                popupMenu.show(e.getComponent(), e.getX(), e.getY());
            }
        }

        @Override
        public void mousePressed(MouseEvent e) {
            showPopup(e);
        }

        @Override
        public void mouseReleased(MouseEvent e) {
            showPopup(e);
        }
    }

    double[] X = {0, 1, 3};
    double[] Y = {2, 1, 4};
    double[] minusY = {1, 1, 0.5};
    double[] plusY = {0.5, 2, 2};


    private IntervalXYDataset createIntervalXYDataset() {
        XYIntervalSeriesCollection collection = new XYIntervalSeriesCollection();
        XYIntervalSeries series = new XYIntervalSeries("");
        for (int i = 0; i < X.length; i++) {
            series.add(X[i], X[i], X[i], Y[i], Y[i] - minusY[i], Y[i] + plusY[i]);
        }
        collection.addSeries(series);
        return collection;
    }

    private JFreeChart createIntervalXYChart(IntervalXYDataset intervalXYDataset) {
        NumberAxis axisX = new NumberAxis("X");
        NumberAxis axisY = new NumberAxis("Y");
        XYErrorRenderer renderer = new XYErrorRenderer();
        renderer.setErrorPaint(Color.RED);
        renderer.setSeriesPaint(0, Color.BLUE);
        XYLineAndShapeRenderer xyLineAndShapeRenderer
                = new XYLineAndShapeRenderer();
        xyLineAndShapeRenderer.setSeriesPaint(0, Color.BLUE);
        XYPlot plot = new XYPlot(intervalXYDataset, axisX, axisY, renderer);
        plot.setDataset(1, intervalXYDataset);
        plot.setRenderer(1, xyLineAndShapeRenderer);
        plot.setBackgroundPaint(Color.lightGray);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);


        JFreeChart chart = new JFreeChart("XYErrorRenderer Demo 1", plot);
        chart.removeLegend();
        chart.setBackgroundPaint(Color.white);
        return chart;
    }

    private JPanel createPanel() {
        return new ChartPanel(createIntervalXYChart(createIntervalXYDataset()));
    }

    void run() {
        JPanel panel = createPanel();
        JFrame frame = new JFrame("TestFrame");
        frame.setSize(400, 600);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(panel, BorderLayout.CENTER);
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        new JClusteringProfileStampByChart().run();
    }

    private JClusteringProfileStampByChart() {
    }
}
