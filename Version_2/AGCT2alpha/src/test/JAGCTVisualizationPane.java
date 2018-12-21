package test;

import forDebug.Debug;
import gene.Gene;
import middleMan.Clustering_painter;
import timer.Timer.TimerKey;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

public class JAGCTVisualizationPane extends JPanel implements ActionListener,
        Debuggable {
    static List<JAGCTVisualizationPane> instances = new ArrayList<JAGCTVisualizationPane>();

    static interface ClusterSelectedAction {
        void update(int clusterId);
    }

    ClusterSelectedAction clusteringSelectedAction;
    public static final String M_P = "manifoldPanel";

    @Override
    public String getName() {
        // TODO Auto-generated method stub
        return myReferenceName;
    }

    public static final String P_P = "pcaPanel";

    final static String C_C_P = "correlationPanel";

    final static String allManifoldButtonString = "All";
    final static String onlyPointedButtonString = "Pointed";
    final static String noPointedButtonString = "None";
    static String softClusteringString = "Show";

    static String geneIDString = "ID";
    static String geneTrueNameString = "Name";
    static String geneFString = "F:";
    static String genePString = "P:";
    static String geneCString = "C:";
    static String geneStatString = "S:";
    static String statPString = "P:";
    static String statFString = "F:";
    static String statCString = "C:";
    static String starString = "Star";
    static String statTextString = "(Pattern)";
    static int DEFAULT_HIERARCHY_LEVEL = -1;

    static String[] DEFAULT_SEARCH_TEXT = {"Gene to find", "Feature to find"};
    static String[] DEFAULT_STRING_CLUSTERING = {"C#0"};
    static String[] DEFAULT_CLUSTERS = {"---"};

    static int maxRadiusMemberships = 65;
    static int minRadiusMemberships = 5;
    static int numberTicksAlphaStar = 500;

    static double factorRadiusPrototypeCluster = 1.5;

    static int Max_Shadow_Color = Clustering.RGB_Max - 30;

    static final String[] allMembershipViews = {"Disks", "Circles"};

    static final String[] allViews = {"Up XYZ", "Up Z", "Down Z", "Up Y",
            "Down Y", "Up X", "Down X"};
    static final Pnt3D[][] allViewPoints = {{Pnt3D.ijk, null},
            {Pnt3D.k.scale(-1), null}, {Pnt3D.k, null},
            {Pnt3D.j.scale(-1), Pnt3D.i}, {Pnt3D.j, Pnt3D.i},
            {Pnt3D.i.scale(-1), null}, {Pnt3D.i, null}};

    JTextField xText, yText, zText, searchText, statText;
    JLabel xLabel, yLabel, zLabel, searchLabel, hiLabel;
    JButton manifoldEigenButton, delauConsistencyButton, pcaEigenButton,
            saveClusteringButton, loadClusteringButton, delauButton,
            statButton, hiPlusButton, hiMinusButton;
    boolean noStatButton, noHiButtons;

    public JAGCTGraphicsPanel visualizationPanel;
    JComboBox viewSelect, membershipViewSelect, clusteringSelect,
            clusterChoice;

    String informationManifold;

    JRadioButton allManifoldButton, onlyPointedButton, noPointedButton;
    JCheckBox clusterStructure, memberships, geneID, geneTrueName, geneF,
            geneP, geneC, geneStat, statF, statP, statC, stars, keepManifold;
    JSlider sliderRadiusMemberships, sliderAlphaStars;

    AGCT myAGCT;
    Domain myDomain;
    String myTabName;

    public String myReferenceName;

    String searchString;
    ImageIcon myImageIcon, triangulationIdleIcon, triangulationRunningIcon;

    public int xAxis;

    public int yAxis;

    public int zAxis;

    int radiusMemberships, searchIndex, currentClustering;

    boolean manifoldAvailable, pcaAvailable, correlationAvailable;
    // false when no dataset is loaded (no selection of ligands or genes is
    // possible)

    int[] depthOrderSelectedGenes;

    // for Graphical outputs: used to order genes from the deepmost

    public static String stringClustering(int n) {
        String val;
        val = "C#" + n;
        return val;
    }

    public static String stringClustering(String s) {
        return s;
    }

    static Pnt3D[] getViewPoint(String s) {
        boolean trouve = false;
        int i = 0;
        Pnt3D pnt[] = null;
        do {
            if (allViews[i].equals(s))
                trouve = true;
            else
                i++;
        } while ((trouve == false) && (i < allViews.length));
        if (trouve == false)
            Matrix.perror("String not found");
        else {
            pnt = new Pnt3D[2];
            pnt = allViewPoints[i];
        }
        return pnt;
    }

    public JAGCTVisualizationPane(AGCT a, Domain d, String name, ImageIcon ic,
                                  String referenceName) {

        statP = statF = statC = null;

        informationManifold = JAGCTVisualizationPane.allManifoldButtonString;
        sliderRadiusMemberships = null;

        searchString = "-";

        myReferenceName = referenceName;
        myTabName = name;
        manifoldAvailable = false;
        noStatButton = true;
        noHiButtons = true;

        myAGCT = a;
        myDomain = d;
        myImageIcon = ic;

        xAxis = 0;
        yAxis = 1;
        zAxis = 2;
        radiusMemberships = 10;
        searchIndex = -1;
        currentClustering = -1;

        triangulationIdleIcon = new ImageIcon(
                Toolkit.getDefaultToolkit()
                        .getImage(
                                java.net.URLClassLoader
                                        .getSystemResource("Images/set4_triangulate.gif")));
        triangulationIdleIcon.setDescription("Idle");

        triangulationRunningIcon = new ImageIcon(Toolkit.getDefaultToolkit()
                .getImage(
                        java.net.URLClassLoader
                                .getSystemResource("Images/cog.png")));
        triangulationRunningIcon.setDescription("Running");

        visualizationPanel = new JAGCTGraphicsPanel(this, myAGCT, myDomain,
                "visualizationPanel");
        visualizationPanel.setBackground(Color.WHITE);

        visualizationPanel.addMouseListener(visualizationPanel);
        visualizationPanel.addMouseMotionListener(visualizationPanel);
        visualizationPanel.addKeyListener(visualizationPanel);
        visualizationPanel.addComponentListener(visualizationPanel);
        visualizationPanel.addMouseWheelListener(visualizationPanel);

        setLayout(new BorderLayout());
        add(visualizationPanel, BorderLayout.CENTER);

        instances.add(this);
    }

    public Rectangle getCaptureRectangle() {
        return visualizationPanel.getCaptureRectangle();
    }

    public void disableButtons() {
        if ((myReferenceName.equals(JAGCTVisualizationPane.M_P))
                || (myReferenceName.equals(JAGCTVisualizationPane.P_P))) {
            setManifoldButtons(false);
            setGenesButtons(false);
            setStatButton(false);
            setHiButtons(false);
        }
        setMembershipsButtons(false);
    }

    public void toTrash() {
        if (ControlProcess.hasTrue("filtersSelected")) {
            disableButtons();
            if (myReferenceName.equals(JAGCTVisualizationPane.M_P))
                keepManifold.setEnabled(true);
        }

        depthOrderSelectedGenes = null;
        informationManifold = JAGCTVisualizationPane.allManifoldButtonString;
        sliderRadiusMemberships = null;
        searchString = "-";
        statP = statF = statC = null;

        xAxis = 0;
        yAxis = 1;
        zAxis = 2;
        radiusMemberships = 10;
        searchIndex = -1;
        currentClustering = -1;

        manifoldAvailable = pcaAvailable = correlationAvailable = false;
        removeAll();
        visualizationPanel = new JAGCTGraphicsPanel(this, myAGCT, myDomain,
                "visualizationPanel");
        visualizationPanel.setBackground(Color.WHITE);

        visualizationPanel.addMouseListener(visualizationPanel);
        visualizationPanel.addMouseMotionListener(visualizationPanel);
        visualizationPanel.addKeyListener(visualizationPanel);
        visualizationPanel.addComponentListener(visualizationPanel);
        visualizationPanel.addMouseWheelListener(visualizationPanel);

        noStatButton = true;
        noHiButtons = true;

        setLayout(new BorderLayout());
        add(visualizationPanel, BorderLayout.CENTER);
    }

    /***********************************************************
     * Methods to give the right point and vars names, etc
     *****/

    /**
     * 計算された点(座標変換されない)を返す.
     *
     * @param iSel
     * @return
     */
    public Pnt3D getRightPoint(int iSel) {
        Pnt3D val = null;
        if (myReferenceName.equals(JAGCTVisualizationPane.M_P)) {
            val = ((Gene) myDomain.getGenes().get(
                    myDomain.selectedGeneNumberToGeneNumber[iSel])).manifold_Pnt3D;

            if (myAGCT.data.isReflectionZ())
                val.coordinates[2] = -val.coordinates[2];

            if (myAGCT.data.isReflectionX())
                val.coordinates[0] = -val.coordinates[0];

            if (myAGCT.data.isReflectionY())
                val.coordinates[1] = -val.coordinates[1];

            // System.out.println("!");
        } else if (myReferenceName.equals(JAGCTVisualizationPane.P_P))
            val = ((Gene) myDomain.getGenes().get(
                    myDomain.selectedGeneNumberToGeneNumber[iSel])).pca_Pnt3D;
        else if (myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
            val = (Pnt3D) myAGCT.data.getFeature_pca_Pnt3D().get(iSel);
        else
            Matrix.perror("No Right point to give");
        return val;
    }

    public int getRightIndex(int iSel) {
        int val = 0;
        if (myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
            val = iSel;
        else
            val = depthOrderSelectedGenes[iSel];
        return val;
    }

    public void drawLabelVsRef(View3D aView3D, int iSel,
                               JAGCTGraphicsPanel pane, boolean highlightRef) {
        Gene gg;
        String val;
        Pnt3D rp = getRightPoint(iSel);
        if ((myReferenceName.equals(JAGCTVisualizationPane.M_P))
                || (myReferenceName.equals(JAGCTVisualizationPane.P_P))) {
            gg = (Gene) myDomain.getGenes().get(
                    myDomain.selectedGeneNumberToGeneNumber[iSel]);
            aView3D.drawGeneVsRef(gg, rp, pane, false, highlightRef);
        } else {
            val = myAGCT.getFeatureName(iSel);
            aView3D.drawStringVsRef(val, rp);
        }
    }

    public Gene getDirectGene(int iGene) {
        return (Gene) myDomain.getGenes().get(iGene);
    }

    public String getFeatureName(int iSel) {
        String val = "";
        int i;
        Gene gg;
        if ((myReferenceName.equals(JAGCTVisualizationPane.M_P))
                || (myReferenceName.equals(JAGCTVisualizationPane.P_P))) {
            gg = (Gene) myDomain.getGenes().get(
                    myDomain.selectedGeneNumberToGeneNumber[iSel]);
            if ((geneID.isSelected()) && (geneID.isEnabled()))
                val += gg.getName();
        } else if (myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
            val = myAGCT.getFeatureName(iSel);
        else
            Matrix.perror("No Right point to use");
        return val;
    }

    public boolean testOnFeature(int iSel, String s) {
        int i, gid, nid;
        Vector vid;

        Gene gg, gveci;
        if ((myReferenceName.equals(JAGCTVisualizationPane.M_P))
                || (myReferenceName.equals(JAGCTVisualizationPane.P_P))) {
            gg = (Gene) myDomain.getGenes().get(
                    myDomain.selectedGeneNumberToGeneNumber[iSel]);
            if ((geneID.isSelected()) && (geneID.isEnabled())) {
                if (s.equals(gg.getName()))
                    return true;
            }
            if ((geneTrueName.isSelected()) && (geneTrueName.isEnabled())) {

                if ((gg.asciiName != null)
                        && (!gg.asciiName.equals(new String("")))
                        && ((gg.asciiName.lastIndexOf(s)) != -1))
                    return true;
            }

            if (!Prototype.No_Reduction) {
                gid = myDomain.selectedGeneNumberToGeneNumber[iSel];
                vid = (Vector) Prototype.Closest_Center_To_Cluster_Points
                        .get(new Integer(gid));
                if (vid.size() >= 1) {
                    for (i = 0; i < vid.size(); i++) {
                        nid = ((Integer) vid.get(i)).intValue();
                        gveci = getDirectGene(nid);

                        if ((geneID.isSelected()) && (geneID.isEnabled())) {
                            if (s.equals(gveci.getName()))
                                return true;
                        }
                        if ((geneTrueName.isSelected())
                                && (geneTrueName.isEnabled())) {

                            if ((gveci.asciiName != null)
                                    && (!gveci.asciiName.equals(new String("")))
                                    && ((gveci.asciiName.lastIndexOf(s)) != -1))
                                return true;
                        }
                    }
                }
            }
        } else if (myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) {
            if (s.equals(myAGCT.getFeatureName(iSel)))
                return true;
        } else
            Matrix.perror("No Right point to use");
        return false;
    }

    public int getTotalFeatures() {
        int val = -1;
        if ((myReferenceName.equals(JAGCTVisualizationPane.M_P))
                || (myReferenceName.equals(JAGCTVisualizationPane.P_P)))
            val = myDomain.numberSelectedGenes;
        else if (myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
            val = myAGCT.data.getDimFeatures();
        return val;
    }

    /**
     * **
     * End of these methods
     * **********************************************************
     */

    public void modificationSearchText(String s) {
        if ((s != null) && (!s.equals("")) && (!s.equals("-"))) {
            int i = 0;
            boolean found = false;
            do {
                if (testOnFeature(i, s)) // (s.equals(getFeatureName(i)))
                    found = true;
                else
                    i++;
            } while ((found == false) && (i < getTotalFeatures()));
            if (found) {
                searchIndex = i;
                searchString = s;
                visualizationPanel.reference_Pnt3D = getRightPoint(i);
                visualizationPanel.repaint();
            } else {
                searchString = "-";
            }
            searchText.setText(searchString);
        }
    }

    public void initCorrelation3D() {
        myAGCT.initCorrelation3D(xAxis, yAxis, zAxis);
        myAGCT.computeMinMax_CorrelationPnt3D();
    }

    public void updateManifold3D() {
        int i, j;
        Gene g;
        for (i = 0; i < myDomain.numberSelectedGenes; i++) {
            g = ((Gene) myDomain.getGenes().get(
                    myDomain.selectedGeneNumberToGeneNumber[i]));
            g.updateManifold_Pnt3D(xAxis, yAxis, zAxis);
        }
        myDomain.computeMinMax_ManifoldPnt3D();
        repaint();
    }

    public void updatePCA3D() {
        int i, j;
        Gene g;
        for (i = 0; i < myDomain.numberSelectedGenes; i++) {
            g = ((Gene) myDomain.getGenes().get(
                    myDomain.selectedGeneNumberToGeneNumber[i]));
            g.updatePCA_Pnt3D(xAxis, yAxis, zAxis);
        }
        myDomain.computeMinMax_PCAPnt3D();
        repaint();
    }

    public void initManifold3D() {
        System.out.println("JAGCTVisualizationPane.initManifold3D()");

        int i;
        Gene g;
        for (i = 0; i < myDomain.numberSelectedGenes; i++) {
            g = ((Gene) myDomain.getGenes().get(
                    myDomain.selectedGeneNumberToGeneNumber[i]));
            g.initManifold3D(xAxis, yAxis, zAxis);
        }
        myDomain.computeMinMax_ManifoldPnt3D();

        /*****
         * No more rescaling for (i=0;i<myDomain.numberSelectedGenes;i++){ g = (
         * (Gene) myDomain.domainGenes.get(myDomain.
         * selectedGeneNumberToGeneNumber[i])); g.initManifold_Pnt3D(); }
         *****/
        repaint();

        System.out.println("END_JAGCTVisualizationPane.initManifold3D()");
    }

    public void initPCA3D() {
        int i, j;
        Gene g;
        for (i = 0; i < myDomain.numberSelectedGenes; i++) {
            g = ((Gene) myDomain.getGenes().get(
                    myDomain.selectedGeneNumberToGeneNumber[i]));
            g.initPCA3D(xAxis, yAxis, zAxis);
        }
        myDomain.computeMinMax_PCAPnt3D();

        /*****
         * No more rescaling for (i=0;i<myDomain.numberSelectedGenes;i++){ g = (
         * (Gene) myDomain.domainGenes.get(myDomain.
         * selectedGeneNumberToGeneNumber[i])); g.initManifold_Pnt3D(); }
         *****/

        repaint();
    }

    public void hiPlus() {
        visualizationPanel.hiPlus();
        setHiLabel();
        repaint();
    }

    public void hiMinus() {
        visualizationPanel.hiMinus();
        setHiLabel();
        repaint();
    }

    public void modificationXText(String s) {
        int v = Integer.parseInt(s);
        int i;
        Gene g;

        if (myTabName.equals(JMainFrameTabbedPane.THE_MANIFOLD_PANE_NAME))
            if ((v >= 0) && (v < myAGCT.data.getDimElements())
                    && (v < AGCT.Number_Of_Manifold_Components) && (v != yAxis)
                    && (v != zAxis)) {
                xAxis = v;
                updateManifold3D();
            } else {
                xText.setText(new Integer(xAxis).toString());
            }
        else if ((myTabName.equals(JMainFrameTabbedPane.THE_PCA_PANE_NAME))
                || (myTabName
                .equals(JMainFrameTabbedPane.THE_CORRELATION_PANE_NAME)))
            if ((v >= 0) && (v < myAGCT.data.getDimFeatures()) && (v != yAxis)
                    && (v != zAxis)) {
                xAxis = v;
                if (myTabName.equals(JMainFrameTabbedPane.THE_PCA_PANE_NAME))
                    updatePCA3D();
                else
                    myAGCT.updateCorrelationComponents(xAxis, yAxis, zAxis);
            } else {
                xText.setText(new Integer(xAxis).toString());
            }
        repaint();
    }

    public void modificationYText(String s) {
        int v = Integer.parseInt(s);
        int i;
        Gene g;
        if (myTabName.equals(JMainFrameTabbedPane.THE_MANIFOLD_PANE_NAME))
            if ((v >= 0) && (v < myAGCT.data.getDimElements())
                    && (v < AGCT.Number_Of_Manifold_Components) && (v != xAxis)
                    && (v != zAxis)) {
                yAxis = v;
                updateManifold3D();
            } else {
                yText.setText(new Integer(yAxis).toString());
            }
        else if ((myTabName.equals(JMainFrameTabbedPane.THE_PCA_PANE_NAME))
                || (myTabName
                .equals(JMainFrameTabbedPane.THE_CORRELATION_PANE_NAME)))
            if ((v >= 0) && (v < myAGCT.data.getDimFeatures()) && (v != xAxis)
                    && (v != zAxis)) {
                yAxis = v;
                if (myTabName.equals(JMainFrameTabbedPane.THE_PCA_PANE_NAME))
                    updatePCA3D();
                else
                    myAGCT.updateCorrelationComponents(xAxis, yAxis, zAxis);
            } else {
                yText.setText(new Integer(yAxis).toString());
            }
        repaint();
    }

    public void modificationZText(String s) {
        int v = Integer.parseInt(s);
        int i;
        Gene g;
        if (myTabName.equals(JMainFrameTabbedPane.THE_MANIFOLD_PANE_NAME))
            if ((v >= 0) && (v < myAGCT.data.getDimElements())
                    && (v < AGCT.Number_Of_Manifold_Components) && (v != yAxis)
                    && (v != xAxis)) {
                zAxis = v;
                updateManifold3D();
            } else {
                zText.setText(new Integer(zAxis).toString());
            }
        else if ((myTabName.equals(JMainFrameTabbedPane.THE_PCA_PANE_NAME))
                || (myTabName
                .equals(JMainFrameTabbedPane.THE_CORRELATION_PANE_NAME)))
            if ((v >= 0) && (v < myAGCT.data.getDimFeatures()) && (v != yAxis)
                    && (v != xAxis)) {
                zAxis = v;
                if (myTabName.equals(JMainFrameTabbedPane.THE_PCA_PANE_NAME))
                    updatePCA3D();
                else
                    myAGCT.updateCorrelationComponents(xAxis, yAxis, zAxis);
            } else {
                zText.setText(new Integer(zAxis).toString());
            }
        repaint();
    }

    public void modificationInfoManifold(String s) {
        informationManifold = s;
        repaint();
    }

    public void modificationRadius(int i) {
        radiusMemberships = i;
        repaint();
    }

    public void modificationAlpha(int i) {
        repaint();
    }

    public String getXLabel() {
        String s = "X-axis: ";
        return s;
    }

    public void setHiLabel() {
        hiLabel.setText(visualizationPanel.hierarchyF + "");
    }

    public String getYLabel() {
        String s = "Y-axis: ";
        return s;
    }

    public String getZLabel() {
        String s = "Z-axis: ";
        return s;
    }

    public void setDomain(Domain d) {
        myDomain = d;
        visualizationPanel.setDomain(d);
    }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        visualizationPanel.paintComponent(g);
    }

    public void displayComponents() {
        if (myReferenceName.equals(JAGCTVisualizationPane.M_P))
            displayManifoldComponents();
        else if (myReferenceName.equals(JAGCTVisualizationPane.P_P))
            displayPCAComponents();
        else if (myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
            displayCorrelationComponents();

        if ((myReferenceName.equals(JAGCTVisualizationPane.M_P))
                || (myReferenceName.equals(JAGCTVisualizationPane.P_P)))
            displayInfoFilter();
        else
            displayInfoFilterCorrelation();
    }

    public void displayInfoFilter() {
        while (!ControlProcess.hasTrue("selectedGenesAllAnnotationsOK"))
            ;

        JToolBar toolbarManifold = new JToolBar();
        JToolBar toolbarGenes = new JToolBar();
        JToolBar toolbarMemberships = new JToolBar();
        JToolBar toolbarManifoldInside = new JToolBar();

        toolbarManifold.setFloatable(false);
        toolbarGenes.setFloatable(false);
        toolbarMemberships.setFloatable(false);
        toolbarManifoldInside.setFloatable(false);

        ButtonGroup manifoldGroup = new ButtonGroup();

        allManifoldButton = new JRadioButton(allManifoldButtonString);
        onlyPointedButton = new JRadioButton(onlyPointedButtonString);
        noPointedButton = new JRadioButton(noPointedButtonString);
//		allManifoldButton.setSelected(true);
        noPointedButton.setSelected(true);
        setManifoldButtons(false);

        manifoldGroup.add(allManifoldButton);
        manifoldGroup.add(onlyPointedButton);
        manifoldGroup.add(noPointedButton);

        toolbarManifoldInside.add(allManifoldButton);
        toolbarManifoldInside.add(onlyPointedButton);
        toolbarManifoldInside.add(noPointedButton);
        toolbarManifoldInside.setBorder(BorderFactory.createTitledBorder(""));

        allManifoldButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                modificationInfoManifold(allManifoldButton.getText());
            }
        });

        onlyPointedButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                modificationInfoManifold(onlyPointedButton.getText());
            }
        });

        noPointedButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                modificationInfoManifold(noPointedButton.getText());
            }
        });

        // toolbarManifold.add(allManifoldButton);
        // toolbarManifold.add(onlyPointedButton);
        // toolbarManifold.add(noPointedButton);
        toolbarManifold.add(toolbarManifoldInside);
        toolbarManifold.setBorder(BorderFactory
                .createTitledBorder("Manifold edges"));

        geneID = new JCheckBox(JAGCTVisualizationPane.geneIDString);
        geneID.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });

        geneID.setSelected(false);
        geneID.setEnabled(false);

        geneTrueName = new JCheckBox(JAGCTVisualizationPane.geneTrueNameString);
        geneTrueName.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });

        geneTrueName.setSelected(false);
        geneTrueName.setEnabled(false);

        geneF = new JCheckBox(JAGCTVisualizationPane.geneFString);
        geneF.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });
        geneF.setSelected(false);
        geneF.setEnabled(false);

        geneP = new JCheckBox(JAGCTVisualizationPane.genePString);
        geneP.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });
        geneP.setSelected(false);
        geneP.setEnabled(false);

        geneC = new JCheckBox(JAGCTVisualizationPane.geneCString);
        geneC.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });
        geneC.setSelected(false);
        geneC.setEnabled(false);

        geneStat = new JCheckBox(JAGCTVisualizationPane.geneStatString);
        geneStat.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });
        geneStat.setSelected(false);
        geneStat.setEnabled(false);

        stars = new JCheckBox(JAGCTVisualizationPane.starString);
        stars.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });
        stars.setSelected(false);
        stars.setEnabled(false);

        sliderAlphaStars = new JSlider(JSlider.HORIZONTAL, 0,
                JAGCTVisualizationPane.numberTicksAlphaStar, 0);
        sliderAlphaStars.setPreferredSize(new Dimension(50, 10));
        sliderAlphaStars.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                modificationAlpha(sliderAlphaStars.getValue());
            }
        });
        // sliderAlphaStars.setMajorTickSpacing(30);
        // sliderAlphaStars.setMinorTickSpacing(15);
        sliderAlphaStars.setPaintTicks(false);
        sliderAlphaStars.setPaintLabels(false);
        sliderAlphaStars.setEnabled(true);

        String tit = "Gene annotations";

        if (!Prototype.No_Reduction) {
            JToolBar jta = new JToolBar();
            jta.setBorder(BorderFactory.createTitledBorder(""));

            jta.add(stars);
            jta.add(sliderAlphaStars);

            toolbarGenes.add(jta);
            tit += " & stars";
        }

        JToolBar jtbnames = new JToolBar();
        jtbnames.setFloatable(false);
        JToolBar jtbanno = null;

        if (myAGCT.data.isAnnotationsExist())
            jtbnames.setBorder(BorderFactory.createTitledBorder(""));

        hiPlusButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit()
                .getImage(
                        java.net.URLClassLoader
                                .getSystemResource("Images/arrow_up.png"))));
        hiPlusButton.setBorderPainted(false);
        hiPlusButton.setToolTipText("hierarchy level +1");
        hiPlusButton.setActionCommand("hiPlus");
        hiPlusButton.addActionListener(this);

        hiMinusButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit()
                .getImage(
                        java.net.URLClassLoader
                                .getSystemResource("Images/arrow_down.png"))));
        hiMinusButton.setBorderPainted(false);
        hiMinusButton.setToolTipText("hierarchy level -1");
        hiMinusButton.setActionCommand("hiMinus");
        hiMinusButton.addActionListener(this);

        hiLabel = new JLabel("");
        hiLabel.setToolTipText("actual hierarchy value");
        setHiLabel();

        jtbnames.add(geneID);
        if (myAGCT.data.isAnnotationsExist()) {
            jtbnames.add(geneTrueName);

            if ((myAGCT.data.isAnnotationsExistF())
                    || (myAGCT.data.isAnnotationsExistP())
                    || (myAGCT.data.isAnnotationsExistC())) {
                jtbanno = new JToolBar();
                jtbanno.setFloatable(false);
                jtbanno.setBorder(BorderFactory.createTitledBorder(""));

                if (myDomain.maxHierarchy > 1) {
                    jtbanno.add(hiPlusButton);
                    jtbanno.add(hiLabel);
                    jtbanno.add(hiMinusButton);
                }
            }

            if (myAGCT.data.isAnnotationsExistF())
                jtbanno.add(geneF);
            if (myAGCT.data.isAnnotationsExistP())
                jtbanno.add(geneP);
            if (myAGCT.data.isAnnotationsExistC())
                jtbanno.add(geneC);

            // if ( (myAGCT.annotationsExistF) || (myAGCT.annotationsExistP) ||
            // (myAGCT.annotationsExistC) )
            // jtbanno.add(geneStat);
        }

        toolbarGenes.add(jtbnames);
        if ((myAGCT.data.isAnnotationsExistF())
                || (myAGCT.data.isAnnotationsExistP())
                || (myAGCT.data.isAnnotationsExistC())) {
            toolbarGenes.add(jtbanno);
        }

        if ((myAGCT.data.isAnnotationsExistF())
                || (myAGCT.data.isAnnotationsExistP())
                || (myAGCT.data.isAnnotationsExistC()))
            toolbarGenes.setBorder(BorderFactory
                    .createTitledBorder("Gene names & annotations"));
        else
            toolbarGenes.setBorder(BorderFactory
                    .createTitledBorder("Gene names"));

        memberships = new JCheckBox(JAGCTVisualizationPane.softClusteringString);
        memberships.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });

        memberships.setSelected(false);
        memberships.setEnabled(false);

        clusterStructure = new JCheckBox("Show");
        clusterStructure.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });

        clusterStructure.setSelected(false);
        clusterStructure.setEnabled(false);

        clusteringSelectedAction = new ClusterSelectedAction() {
            @Override
            public void update(int clusterId) {
                setCurrentClustering(clusterId);
                updateStatButtons();
                repaint();
            }
        };
        makeClusteringSelect();

        clusterChoice = new JComboBox(JAGCTVisualizationPane.DEFAULT_CLUSTERS);
        clusterChoice.setToolTipText("select the cluster(s)");
        clusterChoice.setSelectedIndex(0);
        clusterChoice.setEnabled(false);
        clusterChoice.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                repaint();
            }
        });

        membershipViewSelect = new JComboBox(
                JAGCTVisualizationPane.allMembershipViews);
        membershipViewSelect.setToolTipText("select the membership L&F");
        membershipViewSelect.setSelectedIndex(0);
        membershipViewSelect.setEnabled(false);
        membershipViewSelect.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                repaint();
            }
        });

        sliderRadiusMemberships = new JSlider(JSlider.HORIZONTAL,
                JAGCTVisualizationPane.minRadiusMemberships,
                JAGCTVisualizationPane.maxRadiusMemberships, radiusMemberships);
        sliderRadiusMemberships.setPreferredSize(new Dimension(100, 10));
        sliderRadiusMemberships.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                modificationRadius(sliderRadiusMemberships.getValue());
            }
        });
        // sliderRadiusMemberships.setMajorTickSpacing(30);
        // sliderRadiusMemberships.setMinorTickSpacing(15);
        sliderRadiusMemberships.setPaintTicks(false);
        sliderRadiusMemberships.setPaintLabels(false);
        sliderRadiusMemberships.setEnabled(false);

        JToolBar jtb1 = new JToolBar();
        JToolBar jtb2 = new JToolBar();
        JToolBar jtb3 = new JToolBar();

        jtb1.add(clusterStructure);
        jtb1.add(clusteringSelect);

        jtb3.add(memberships);
        jtb2.add(clusterChoice);
        jtb3.add(membershipViewSelect);
        jtb3.add(sliderRadiusMemberships);

        jtb1.setBorder(BorderFactory.createTitledBorder(""));
        jtb2.setBorder(BorderFactory.createTitledBorder(""));
        jtb3.setBorder(BorderFactory.createTitledBorder(""));

        toolbarMemberships.add(jtb1);
        toolbarMemberships.add(jtb2);
        toolbarMemberships.add(jtb3);

        JToolBar superBar = new JToolBar();
        superBar.setFloatable(false);

        superBar.add(toolbarMemberships);

        if ((myDomain.someAnnotationsF) || (myDomain.someAnnotationsP)
                || (myDomain.someAnnotationsC)) {
            JToolBar statBar = new JToolBar();
            createIncludeStatButtons(statBar);
            statBar.setBorder(BorderFactory
                    .createTitledBorder("Chi2 tests (clusters)"));
            superBar.add(statBar);
            noStatButton = false;
            noHiButtons = false;
        } else {
            noStatButton = true;
            noHiButtons = true;
        }

		/*
         * toolbarMemberships.add(clusterStructure);
		 * toolbarMemberships.add(clusteringSelect);
		 * toolbarMemberships.add(memberships);
		 * toolbarMemberships.add(membershipViewSelect);
		 * toolbarMemberships.add(sliderRadiusMemberships);
		 */

        toolbarMemberships.setBorder(BorderFactory
                .createTitledBorder("Clusterings, Clusters and Memberships"));

        // JPanel filterPanel = new JPanel();
        // filterPanel.setLayout(new BorderLayout());
        // filterPanel.add(toolbarManifold, BorderLayout.WEST);
        // filterPanel.add(toolbarGenes, BorderLayout.CENTER);
        // filterPanel.add(superBar, BorderLayout.EAST);

        JToolBar megaTool = new JToolBar();
        megaTool.setFloatable(true);
        megaTool.add(toolbarManifold);
        megaTool.add(Box.createHorizontalGlue());
        megaTool.add(toolbarGenes);
        megaTool.add(Box.createHorizontalGlue());
        megaTool.add(superBar);
        add(megaTool, BorderLayout.SOUTH);
    }

    private void makeClusteringSelect() {
        clusteringSelect = new JComboBox(
                new String[]{JAGCTVisualizationPane.stringClustering(0)});
        clusteringSelect.setToolTipText("select the clustering");
        clusteringSelect.setSelectedIndex(0);
        clusteringSelect.setEnabled(false);
        clusteringSelect.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                int id = clusteringSelect.getSelectedIndex();
                Debug.debug("JAGCTVisualizationPane # clusteringSelectActionPerformed", id);
                for (JAGCTVisualizationPane instance : instances) {
                    if (instance != JAGCTVisualizationPane.this) {
                        if (instance.clusteringSelect.getItemCount() > id) {
                            instance.clusteringSelect.setSelectedIndex(id);
                            instance.clusteringSelectedAction.update(id);
                        }
                    }
                }
                clusteringSelectedAction.update(id);
            }
        });
    }


    // TODO refactor
    public void doit() {
        if (statText == null)
            return;
        // while(statText==null)System.out.println("statText = null");
        statText.setText(Statistics.TOKEN_ALL_TESTS);
        statButton.doClick();
    }

    public void updateStatButtons() {

        if (!noStatButton)
            if ((myAGCT.data.getNbClustering() < 1)
                    || ((clusteringSelect.getSelectedIndex() >= 0)
                    && (((Clustering) myAGCT.data.getAllClusterings()
                    .get(clusteringSelect.getSelectedIndex())).myClusteringAlgorithm
                    .getOk_for_statistics() == false) && (statButton
                    .isEnabled())))
                setStatButton(false);
            else if ((clusteringSelect.getSelectedIndex() >= 0)
                    && (((Clustering) myAGCT.data.getAllClusterings().get(
                    clusteringSelect.getSelectedIndex())).myClusteringAlgorithm
                    .getOk_for_statistics() == true)
                    && (((Clustering) myAGCT.data.getAllClusterings().get(
                    clusteringSelect.getSelectedIndex())).myClusteringAlgorithm
                    .getNumberOfClusters() > 1)
                    && (!statButton.isEnabled()))
                setStatButton(true);
    }

    public void updateHiButtons() {
        if (!noHiButtons) {
            if ((myTabName.equals(JMainFrameTabbedPane.THE_MANIFOLD_PANE_NAME))
                    && (ControlProcess.hasTrue("manifoldProcessed")))
                setHiButtons(true);
            else if ((myTabName.equals(JMainFrameTabbedPane.THE_PCA_PANE_NAME))
                    && (ControlProcess.hasTrue("pcaProcessed")))
                setHiButtons(true);
            else
                setHiButtons(false);
        }
    }

    public void displayInfoFilterCorrelation() {
        JToolBar toolbarMemberships = new JToolBar();

        memberships = new JCheckBox(JAGCTVisualizationPane.softClusteringString);
        memberships.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                repaint();
            }
        });
        memberships.setSelected(true);
        memberships.setEnabled(false);

        clusteringSelectedAction = new ClusterSelectedAction() {
            @Override
            public void update(int clusteringId) {
                setCurrentClustering(clusteringId);
                repaint();
            }
        };
        makeClusteringSelect();

        toolbarMemberships.add(clusteringSelect);
        toolbarMemberships.add(memberships);
        toolbarMemberships.setBorder(BorderFactory
                .createTitledBorder("Clusterings"));

        JPanel filterPanel = new JPanel();
        filterPanel.setLayout(new BorderLayout());
        filterPanel.add(toolbarMemberships, BorderLayout.EAST);

        add(filterPanel, BorderLayout.SOUTH);
    }

    public void setManifoldButtons(boolean b) {
        allManifoldButton.setEnabled(b);
        onlyPointedButton.setEnabled(b);
        noPointedButton.setEnabled(b);
    }

    public void setGenesButtons(boolean b) {
        Debug.debug("JAGCTVisualizationPane # setGenesButtons", b);
        geneID.setEnabled(b);
        geneTrueName.setEnabled(b);
        if (myAGCT.data.getAnnotationsExist()) {
            geneF.setEnabled(b);
            geneP.setEnabled(b);
            geneC.setEnabled(b);
        }
        geneStat.setEnabled(b);
        stars.setEnabled(b);
    }

    public void setHiButtons(boolean b) {
        hiPlusButton.setEnabled(b);
        hiMinusButton.setEnabled(b);
    }

    public void setStatButton(boolean b) {
        if (myAGCT.data.getAnnotationsExist()) {
            statButton.setEnabled(b);
            statF.setEnabled(b);
            statP.setEnabled(b);
            statC.setEnabled(b);
            statText.setEnabled(b);
        }
    }

    public void setMembershipsButtons(boolean b) {
        memberships.setEnabled(b);
        if (sliderRadiusMemberships != null) {
            clusterStructure.setEnabled(b);
            clusterChoice.setEnabled(b);
            membershipViewSelect.setEnabled(b);
            sliderRadiusMemberships.setEnabled(b);
        }
    }

    public void putBasicButtons(JToolBar toolbar) {
        viewSelect = new JComboBox(JAGCTVisualizationPane.allViews);
        viewSelect.setToolTipText("select the view");
        viewSelect.setSelectedIndex(0);
        viewSelect.setEnabled(false);
        viewSelect.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                visualizationPanel.changeView(viewSelect.getSelectedIndex());
            }
        });

        toolbar.add(viewSelect);

        saveClusteringButton = new JButton(
                new ImageIcon(
                        Toolkit.getDefaultToolkit()
                                .getImage(
                                        java.net.URLClassLoader
                                                .getSystemResource("Images/diskClustering.png"))));
        saveClusteringButton.setToolTipText("save visible clustering");
        saveClusteringButton.setActionCommand("save_clustering");
        saveClusteringButton.addActionListener(this);
        saveClusteringButton.setEnabled(true);

        loadClusteringButton = new JButton(new ImageIcon(Toolkit
                .getDefaultToolkit().getImage(
                        java.net.URLClassLoader
                                .getSystemResource("Images/folder_go.png"))));
        loadClusteringButton.setToolTipText("load clustering");
        loadClusteringButton.setActionCommand("load_clustering");
        loadClusteringButton.addActionListener(this);
        loadClusteringButton.setEnabled(false);

        xText = new JTextField(new Integer(xAxis).toString());
        xText.setEditable(true);
        xText.setEnabled(false);
        xText.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                modificationXText(xText.getText());
            }
        });
        xLabel = new JLabel(getXLabel());

        yText = new JTextField(new Integer(yAxis).toString());
        yText.setEditable(true);
        yText.setEnabled(false);
        yText.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                modificationYText(yText.getText());
            }
        });
        yLabel = new JLabel(getYLabel());

        zText = new JTextField(new Integer(zAxis).toString());
        zText.setEditable(true);
        zText.setEnabled(false);
        zText.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                modificationZText(zText.getText());
            }
        });
        zLabel = new JLabel(getZLabel());

        if ((myReferenceName.equals(JAGCTVisualizationPane.M_P))
                || (myReferenceName.equals(JAGCTVisualizationPane.P_P)))
            searchText = new JTextField(
                    JAGCTVisualizationPane.DEFAULT_SEARCH_TEXT[0]);
        else if (myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
            searchText = new JTextField(
                    JAGCTVisualizationPane.DEFAULT_SEARCH_TEXT[1]);

        searchText.setEditable(true);
        searchText.setEnabled(false);
        searchText.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                modificationSearchText(searchText.getText());
            }
        });
        searchLabel = new JLabel(new ImageIcon(Toolkit.getDefaultToolkit()
                .getImage(
                        java.net.URLClassLoader
                                .getSystemResource("Images/eye.png"))));

        toolbar.add(Box.createHorizontalGlue());

        toolbar.addSeparator();
        toolbar.add(loadClusteringButton);
        toolbar.add(saveClusteringButton);
        toolbar.addSeparator();
        toolbar.add(xLabel);
        toolbar.add(xText);
        toolbar.addSeparator();
        toolbar.add(yLabel);
        toolbar.add(yText);
        toolbar.addSeparator();
        toolbar.add(zLabel);
        toolbar.add(zText);
        toolbar.addSeparator();
        toolbar.add(searchLabel);
        toolbar.add(searchText);
    }

    public void displayCorrelationComponents() {
        correlationAvailable = true;
        JToolBar buttonBar = new JToolBar();

        putBasicButtons(buttonBar);

        add(buttonBar, BorderLayout.NORTH);
    }

    public void displayPCAComponents() {
        pcaAvailable = true;

        pcaEigenButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit()
                .getImage(
                        java.net.URLClassLoader
                                .getSystemResource("Images/pca.gif"))));
        pcaEigenButton.setToolTipText("proceed to pca eigensystem");
        pcaEigenButton.setActionCommand("compute_pca_eigensystem");

        pcaEigenButton.addActionListener(this);
        pcaEigenButton.setEnabled(true);

        JToolBar buttonBar = new JToolBar();
        buttonBar.add(pcaEigenButton);
        buttonBar.addSeparator();

        putBasicButtons(buttonBar);

        add(buttonBar, BorderLayout.NORTH);
    }

    public void displayManifoldComponents() {
        manifoldAvailable = true;

        manifoldEigenButton = new JButton(new ImageIcon(Toolkit
                .getDefaultToolkit().getImage(
                        java.net.URLClassLoader
                                .getSystemResource("Images/manifold.gif"))));
        manifoldEigenButton.setToolTipText("proceed to manifold eigensystem");
        manifoldEigenButton.setActionCommand("compute_manifold_eigensystem");
        // manifoldEigenButton.setEnabled(false);

        keepManifold = new JCheckBox("Keep");
        keepManifold.setToolTipText("Keep the manifold in scenario");
        keepManifold.setEnabled(true);
        keepManifold.setSelected(true);

        delauButton = new JButton(triangulationIdleIcon);
        delauButton.setToolTipText("proceed to triangulation");
        delauButton.setActionCommand("compute_triangulation");
        // delauButton.setEnabled(false);

        delauConsistencyButton = new JButton(new ImageIcon(Toolkit
                .getDefaultToolkit().getImage(
                        java.net.URLClassLoader
                                .getSystemResource("Images/analyze.png"))));
        delauConsistencyButton
                .setToolTipText("make some tests to the current triangulation");
        delauConsistencyButton.setActionCommand("consistency_triangulation");
        // delauConsistencyButton.setEnabled(false);

        manifoldEigenButton.addActionListener(this);
        delauButton.addActionListener(this);
        delauConsistencyButton.addActionListener(this);
        manifoldEigenButton.setEnabled(true);
        delauButton.setEnabled(false);
        delauConsistencyButton.setEnabled(false);

        JToolBar buttonBar = new JToolBar();
        buttonBar.add(manifoldEigenButton);
        buttonBar.add(keepManifold);
        buttonBar.addSeparator();
        buttonBar.add(delauButton);
        buttonBar.add(delauConsistencyButton);
        buttonBar.addSeparator();

        putBasicButtons(buttonBar);

        add(buttonBar, BorderLayout.NORTH);
        // add(JResponsiveGeneSetter.getInstance(),BorderLayout.WEST);
    }

    public void createIncludeStatButtons(JToolBar j) {
        statButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit()
                .getImage(
                        java.net.URLClassLoader
                                .getSystemResource("Images/go.png"))));
        statButton.setToolTipText("statistical tests about clusters");
        statButton.setActionCommand("stats");
        statButton.addActionListener(this);

        statText = new JTextField(JAGCTVisualizationPane.statTextString);
        statText.setEditable(true);

        j.add(statText);

        if (myDomain.someAnnotationsF) {
            statF = new JCheckBox(JAGCTVisualizationPane.statFString);
            statF.setSelected(true);
            j.add(statF);
        }

        if (myDomain.someAnnotationsP) {
            statP = new JCheckBox(JAGCTVisualizationPane.statPString);
            statP.setSelected(true);
            j.add(statP);
        }

        if (myDomain.someAnnotationsC) {
            statC = new JCheckBox(JAGCTVisualizationPane.statCString);
            statC.setSelected(true);
            j.add(statC);
        }

        j.add(statButton);
        updateStatButtons();
    }

    public void initDepthOrderSelectedGenes() {
//		AGCT.debug("numberSelectedGenes : "
//				+ myAGCT.data.getMyDomain().numberSelectedGenes);
        while (myAGCT.data.getMyDomain().numberSelectedGenes < 0)
            ;
        depthOrderSelectedGenes = new int[myAGCT.data.getMyDomain().numberSelectedGenes];

        for (int i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++)
            depthOrderSelectedGenes[i] = i;
    }

    public void updateDepth(View3D view) {
        System.err.println("JAGCTVisualizationPane.updateDepth()");
        Debug.debug(depthOrderSelectedGenes);
        double[] depths = new double[myAGCT.data.getMyDomain().numberSelectedGenes];
        int i;
        for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
            depths[i] = view.depth(getRightPoint(depthOrderSelectedGenes[i]));
        }
        updateDepth(depths);

		/*
         * for (i=0;i<(myAGCT.myDomain.numberSelectedGenes-1)/2;i++){ tmp =
		 * depthOrderSelectedGenes[i]; depthOrderSelectedGenes[i] =
		 * depthOrderSelectedGenes[myAGCT.myDomain.numberSelectedGenes-1-i];
		 * depthOrderSelectedGenes[myAGCT.myDomain.numberSelectedGenes-1-i] =
		 * tmp; }
		 */
        depths = null;
    }

    public void updateDepth(double[] depths) {
        QuickSort.quicksort(depths, depthOrderSelectedGenes);
    }

    /**
     * Manifoldの始まり。ここでtimerのstart finishを行っている
     */
    public void goManifold() {
        System.out.println("JAGCTVisualizationPane.goManifold()");

        AGCT.getInstance().getTimer().start(TimerKey.Manifold);

        if (AGCT.Loading_From_Scenario_Manifold_Begin) {
            while (!AGCT.Loading_From_Scenario_Manifold_End)
                ;
        } else
            myAGCT.manifoldEigensystem();

        AGCT.getInstance().getTimer().finish(TimerKey.Manifold);

        // この時点で、PCApaneが使えるようになっている.
        if (AGCT.AUTO) {
            AGCT.getInstance().myTabbedPane.myManifoldPane.delauButton
                    .doClick();
        }

        Debug.debug("END_goManifold.");
    }

    /**
     * PCAの始まり timerのstart finishを行っている
     */
    public void goPCA() {
        System.out.println("JAGCTVisualizationPane.goPCA()");
        AGCT.getInstance().getTimer().start(TimerKey.PCA);

        initDepthOrderSelectedGenes();
        myAGCT.pcaEigensystem();

        AGCT.getInstance().getTimer().finish(TimerKey.PCA);
    }

    public void goTriangulation() {
        myAGCT.triangulation();
    }

    public void actionPerformed(ActionEvent e) {
        boolean bdum;
        if (e.getActionCommand().equals("compute_manifold_eigensystem")) {
            if (!ControlProcess.hasTrue("manifoldProcessed")) {
                ControlProcess.put("manifoldProcessed", true);
                bdum = keepManifold.isSelected();
                keepManifold.setEnabled(false);
                goManifold();
                if (!bdum)
                    Scenario.add("AGCT_Manifold_Processed");
                Debug.debug("END_goManifold.");
                JOptionPane.showMessageDialog(this, "Validation is finished.");
            }
        } else if ((e.getActionCommand().equals("compute_triangulation"))
                && (ControlProcess.hasTrue("neighborsComputed"))) {
            if (!ControlProcess.hasTrue("manifoldTriangulated")) {
                goTriangulation();
                Scenario.add("AGCT_Manifold_Triangulated");
            }
        } else if (e.getActionCommand().equals("compute_pca_eigensystem")) {
            if (!ControlProcess.hasTrue("pcaProcessed")) {
                goPCA();
                Scenario.add("AGCT_PCA_Processed");

                // この時点で、myClusteringPaneが使える
                if (AGCT.AUTO) {
                    int i;
                    for (i = 0; i < JAGCTClusteringPane.allChoices.length; i++) {
                        if (JAGCTClusteringPane.allChoices[i]
                                .equals("Affinity Propagation"))
                            break;
                    }
                    AGCT.getInstance().myTabbedPane.myClusteringPane.algoSelect
                            .setSelectedIndex(i);// APにする
//					AGCT.getInstance().myTabbedPane.myClusteringPane.goButton
//							.doClick();
                }
            }
        } else if (e.getActionCommand().equals("stats")) {
            if ((((statP != null) && (statP.isSelected()))
                    || ((statF != null) && (statF.isSelected())) || ((statC != null) && (statC
                    .isSelected()))) && (statText.getText() != ""))
                Statistics.getChi2(statP, statF, statC, statText, myAGCT,
                        currentClustering);
        } else if (e.getActionCommand().equals("hiPlus"))
            hiPlus();
        else if (e.getActionCommand().equals("hiMinus"))
            hiMinus();
        else if (e.getActionCommand().equals("consistency_triangulation"))
            myAGCT.triangleConsistency();
        else if (e.getActionCommand().equals("save_clustering"))
            myAGCT.requestSave_ThisClustering(currentClustering);
        else if (e.getActionCommand().equals("load_clustering"))
            myAGCT.requestLoad_Clustering();
    }

    public void setCurrentClustering(int n) {
        currentClustering = n;

        if ((!myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) && (n >= 0)) {
            int i;
            Clustering cc = myAGCT.getClustering(n);
            String[] col = cc.myClusteringAlgorithm.getClusterChoices();

            if (cc.myClusteringAlgorithm.getMyReferenceName().equals(
                    AGCTClustering_Algorithm.REFERENCE_NAME_KM))
                saveClusteringButton.setEnabled(true);

            if (myReferenceName.equals(JAGCTVisualizationPane.M_P)) {
                myAGCT.myTabbedPane.myManifoldPane.clusterChoice
                        .removeAllItems();
                for (i = 0; i < col.length; i++)
                    myAGCT.myTabbedPane.myManifoldPane.clusterChoice
                            .addItem(col[i]);
            }

            if (myReferenceName.equals(JAGCTVisualizationPane.P_P)) {
                myAGCT.myTabbedPane.myPCAPane.clusterChoice.removeAllItems();
                for (i = 0; i < col.length; i++)
                    myAGCT.myTabbedPane.myPCAPane.clusterChoice.addItem(col[i]);
            }
        }
    }

    public void plotClustering(Graphics g) {
        if ((myAGCT.data.getAllClusterings() != null)
                && (myAGCT.data.getAllClusterings().size() > 0))
            if ((currentClustering == Statistics.refNumClustering)
                    && (Statistics.allChi2Tests))
                Clustering_painter.plotWithAnnotations((Clustering) myAGCT.data
                        .getAllClusterings().get(currentClustering), this);
            else
                Clustering_painter.plot((Clustering) myAGCT.data
                        .getAllClusterings().get(currentClustering), this);
    }

    public boolean genePlottedCluster(int ng) {
        if ((myAGCT.data.getAllClusterings() != null)
                && (myAGCT.data.getAllClusterings().size() > 0)
                && (clusterChoice != null)
                && (clusterChoice.getSelectedIndex() != 0)
                && (currentClustering > -1)) {
            int nc = clusterChoice.getSelectedIndex() - 1;
            return myAGCT.getClustering(currentClustering).myClusteringAlgorithm
                    .majorityGeneCluster(ng, nc);
        }
        return true;
    }

    public boolean structPlottedCluster(int nc) {
        if ((clusterChoice.getSelectedIndex() != 0)
                && (clusterChoice.getSelectedIndex() - 1 != nc))
            return false;
        return true;
    }
}