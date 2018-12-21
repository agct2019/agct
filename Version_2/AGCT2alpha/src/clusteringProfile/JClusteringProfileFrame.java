package clusteringProfile;

/*****
 * Classes JClusteringProfile* :
 * Display the average expression profile for each cluster
 * in the current selected clutering
 *****/

import forDebug.Debug;
import gene.Gene;
import test.*;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.Vector;

public class JClusteringProfileFrame extends JFrame {

    public static String Default_Border = "(no clustering plotted)";
    public static String Default_Combo[] = {"(no clustering)"};
    public static String Default_Percentile[] = {"(default: average)"};
    /**
     * domain.getTimes() と同じものが入っている
     * timeSeries
     */
    public static Vector[] Time_Stamps_Summaries;
    // We only store time stamps for selected ligands (saves time)
    public static String[] Ligand_Name_Summaries;
    public JComboBox percentilePlots;
    public JComboBox listPlots;
    public JClusteringProfileGridPlot graphPanel;// 
    JButton goButton, cameraButton, saveProfileButton;
    Box selectionBox;
    public Vector<Integer> keysForCombo;
    // contains the vector of Integers = index in allClusterings of the cluster
    // to be plotted
    AGCT myAGCT;
    boolean onePlot, plotAuthorized, moving;
    private Vector<Integer> listGenes;

    public Vector<Integer> getListGenes() {
        return listGenes;
    }

    // list of genes whose profile can be plotted
    int maxNumberGenes;

    public static void fillTime_Stamps_Summaries(Domain d) {
        if (AGCT.MYDEBUG) {
            AGCT.debug("JClusteringProfileFrame.fillTime_Stamps_Summaries(domain)");
        }
        int i, j, index = 0;
        double cd;
        Vector<Double> current;
        Time_Stamps_Summaries = new Vector[d.getLigands().numberOfSelectedLigands()];
        Ligand_Name_Summaries = new String[d.getLigands().numberOfSelectedLigands()];
        for (i = 0; i < d.getLigands().size(); i++) {
            if (d.getLigands().get(i).isChecked()) {
                current = new Vector<Double>();// d.getTimes()をVectorに入れただけ
                for (j = 0; j < d.getTimes().length; j++) {
                    cd = d.getTimes()[j];
                    current.addElement(cd);
                }
                Time_Stamps_Summaries[index] = current;
                Ligand_Name_Summaries[index] = d.getLigands().get(i).getName();
                index++;
            }
        }
    }

    public boolean atLeastTwoGenes() {
        return listGenes != null && listGenes.size() >= 2;
    }

    public void addGene(int num) {
        if (listGenes == null) {
            listGenes = new Vector<Integer>();
        }
        if ((maxNumberGenes < AGCT.DEFAULT_Max_Number_Genes) || (listGenes.size() <= AGCT.DEFAULT_Max_Number_Genes)) {
            maxNumberGenes = AGCT.DEFAULT_Max_Number_Genes;
        } else {
            Vector<Integer> trim = new Vector<Integer>();
            int i;
            for (i = 0; i < AGCT.DEFAULT_Max_Number_Genes; i++) {
                trim.addElement(listGenes.get(listGenes.size() - AGCT.DEFAULT_Max_Number_Genes + i));
            }
            listGenes = trim;
            maxNumberGenes = AGCT.DEFAULT_Max_Number_Genes;
        }
        if (listGenes.size() >= maxNumberGenes) {
            listGenes.removeElementAt(0);
        }
        if (!listGenes.contains(new Integer(num))) {
            listGenes.addElement(num);
        }
    }

    public void flushGenes() {
        listGenes = new Vector<Integer>();
    }

    public void flushAll() {
        listPlots.removeAllItems();
        listPlots.addItem(JClusteringProfileFrame.Default_Combo);
        keysForCombo = null;
        listPlots.setEnabled(false);
        onePlot = plotAuthorized = false;
        flushGenes();
    }

    public void captureAndSave() {
        Rectangle rect = graphPanel.getCaptureRectangle();
        BufferedImage framecapture = null;
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
            JInformationFrame.getInstance().setText("Saving visible clustering profile pane to file: " + chooser.getSelectedFile().getName());
            ControlProcess.put("frameCapturedClusteringProfile", true);
        }
    }

    public void displayInfo() {
        cameraButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/camera.png"))));
        cameraButton.setToolTipText("capture the visible cluster profiles plot");
        cameraButton.setActionCommand("capture");
        cameraButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                captureAndSave();
            }
        });
        cameraButton.setEnabled(false);
        final JCheckBox errorBarButton = new JCheckBox("SD Bar");
        errorBarButton.setEnabled(false);
        errorBarButton.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent itemEvent) {
                boolean selected = itemEvent.getStateChange() == ItemEvent.SELECTED;
                if (selected) {
                    if (JClusteringProfileStampByChart.renderer != null) {
                        JClusteringProfileStampByChart.setErrorBarEnabled(true);
                    }
                } else {
                    JClusteringProfileStampByChart.setErrorBarEnabled(false);
                }
            }
        });
        goButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/go.png"))));
        goButton.setToolTipText("run / update cluster profiles");
        goButton.setActionCommand("clustering");
        goButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                if (listPlots.getSelectedIndex() >= 0) {
                    errorBarButton.setEnabled(true);
                    goPlotClusters();// ここからクラスタリングの評価に入る。
                }
            }
        });
        goButton.setEnabled(false);


        saveProfileButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/diskProfile.png"))));
        saveProfileButton.setToolTipText("save current profiles");
        saveProfileButton.setActionCommand("save_profiles");
        saveProfileButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                if (listPlots.getSelectedIndex() >= 0) {
                    graphPanel.initAverageProfiles();
                    myAGCT.requestSave_ClusterProfile();
                }
            }
        });
        saveProfileButton.setEnabled(false);

        graphPanel = /*new JScrollPane*/new JClusteringProfileGridPlot(this, myAGCT);// TODO revalidateがうまくいっていません。
        graphPanel.setBorder(BorderFactory.createTitledBorder(JClusteringProfileFrame.Default_Border));

        percentilePlots = new JComboBox(JClusteringProfileFrame.Default_Percentile);
        percentilePlots.setToolTipText("select which percentile to plot");
        percentilePlots.setSelectedIndex(0);

        percentilePlots.addItem("q25");
        percentilePlots.addItem("q50");
        percentilePlots.addItem("q75");
        percentilePlots.setEnabled(false);

        listPlots = new JComboBox(JClusteringProfileFrame.Default_Combo);
        listPlots.setToolTipText("select which clustering data to plot");
        listPlots.setSelectedIndex(0);
        listPlots.setEnabled(false);

        selectionBox = Box.createHorizontalBox();
        // selectionBox.add(cameraButton);
        selectionBox.add(listPlots);
        selectionBox.add(percentilePlots);
//        JToggleButton errorBarButton = new JToggleButton()
        selectionBox.add(errorBarButton);
        selectionBox.add(goButton);
        selectionBox.add(saveProfileButton);

        Container pane = getContentPane();
        pane.setLayout(new BorderLayout());
        pane.add(selectionBox, BorderLayout.NORTH);
        pane.add(new JScrollPane(graphPanel, ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS, ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS)/*graphPanel*/, BorderLayout.CENTER);

//		addWindowListener(new FermetureListener("Closing AGCT's ClusteringProfileFrame\n"));
        setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);

        Image img = Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/information.png"));
        setIconImage(img);
    }


    // クラスタ描画の水源コードです。
    public void goPlotClusters() {
        Debug.debug("BEGIN goPlotClusters");

        graphPanel.clusterVSgenes = true;
        graphPanel.initAverageProfiles();

        onePlot = moving = true;
        updateBorder(currentClusterString() + "   Time-ExpressionLevel Plot");
        Debug.debug("END goPlotClusters");
        repaint();
    }


    /**
     * @param title
     */
    public void updateBorder(String title) {
        if (graphPanel.clusterVSgenes) {
            graphPanel.setBorder(BorderFactory.createTitledBorder(title));
        } else {
            String tt = "";
            int i, numgen;
            Gene gg;
            Domain d = myAGCT.data.getMyDomain();
            if ((listGenes != null) && (listGenes.size() > 0)) {
                tt = "Profiles for genes: ";
                for (i = 0; i < listGenes.size(); i++) {
                    numgen = (listGenes.get(i)).intValue();
                    gg = (Gene) d.getGenes().get(d.selectedGeneNumberToGeneNumber[numgen]);

                    tt += gg.asciiName + " ";
                }
                graphPanel.setBorder(BorderFactory.createTitledBorder(tt));
            } else {
                graphPanel.setBorder(BorderFactory.createTitledBorder(tt));
            }
        }
    }

    /**
     * provide L&F of JClusteringProfileFrame
     *
     * @param cc
     * @param indexClustering
     */
    public void addClusteringLF(Clustering cc, int indexClustering) {

        if (myAGCT.data.getAllClusterings() != null) {
            if (onePlot == false) {
                listPlots.removeAllItems();
                onePlot = plotAuthorized = true;
                keysForCombo = new Vector<Integer>();

                listPlots.setEnabled(true);
                percentilePlots.setEnabled(true);
                cameraButton.setEnabled(true);
                goButton.setEnabled(true);
                saveProfileButton.setEnabled(true);
            }

            listPlots.addItem("C#" + indexClustering + " (" + cc.myClusteringAlgorithm.getMyReferenceName() + ")");
            keysForCombo.addElement(new Integer(indexClustering));

            if (!isVisible()) {
                setVisible(true);
                myAGCT.myFrame.toFront();
            }
        }
    }

    public String currentClusterString() {
        int id = listPlots.getSelectedIndex();
        String rn = ((Clustering) myAGCT.data.getAllClusterings().get(id)).myClusteringAlgorithm.getMyReferenceName();

        return "C#" + id + " (" + rn + ")";
    }

    public JClusteringProfileFrame(String name, AGCT c) {
        super(name);
        setSize(AGCT.WindowWidth, 600);
        myAGCT = c;
        displayInfo();
        setVisible(false);
        setResizable(true);
        onePlot = plotAuthorized = false;
        moving = false;
        listGenes = new Vector<Integer>();
        maxNumberGenes = AGCT.DEFAULT_Max_Number_Genes;
        setMinimumSize(new Dimension(300, 100));
    }
}