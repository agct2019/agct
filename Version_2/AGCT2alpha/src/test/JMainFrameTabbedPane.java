package test;

import filtering.JLeftIndicator;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;


public class JMainFrameTabbedPane extends JTabbedPane implements ActionListener, Debuggable, ChangeListener {
    AGCT myAGCT;
    Domain myDomain;

    public final static String THE_SELECTION_PANE_NAME = "Selection", THE_MANIFOLD_PANE_NAME = "Manifold", THE_PCA_PANE_NAME = "PCA", THE_CORRELATION_PANE_NAME = "Correlation",
            THE_CLUSTERING_PANE_NAME = "Clustering";

    JAGCTSelectionPane mySelectionPane;

    JAGCTVisualizationPane myManifoldPane;
    JAGCTVisualizationPane myPCAPane, myCorrelationPane;
    JAGCTClusteringPane myClusteringPane;

    private String mySelectionPaneName, myManifoldPaneName, myPCAPaneName, myCorrelationPaneName, myClusteringPaneName;

    JMainFrameTabbedPane(AGCT agct, Domain domain) {
        super();
        if (AGCT.MYDEBUG)
            AGCT.debug("JMainFrameTabbedPane.JMainFrameTabbedPane(agct,domain)");
        myAGCT = agct;
        myDomain = domain;

        mySelectionPaneName = JMainFrameTabbedPane.THE_SELECTION_PANE_NAME;
        mySelectionPane = new JAGCTSelectionPane(myAGCT, myDomain, mySelectionPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(
                java.net.URLClassLoader.getSystemResource("Images/book_open.png"))));
        mySelectionPane.setBackground(Color.WHITE);
        mySelectionPane.addComponentListener(new ComponentAdapter() {
            public void componentShown(ComponentEvent evt) {
                lookUpComponentShown(evt);
                myAGCT.saveButtonWeb.setEnabled(false);
            }
        });

        myManifoldPaneName = JMainFrameTabbedPane.THE_MANIFOLD_PANE_NAME;
        myManifoldPane = new JAGCTVisualizationPane(myAGCT, myDomain, myManifoldPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(
                java.net.URLClassLoader.getSystemResource("Images/mapManifold.png"))), JAGCTVisualizationPane.M_P);
        myManifoldPane.setBackground(Color.WHITE);
        myManifoldPane.addComponentListener(new ComponentAdapter() {
            public void componentShown(ComponentEvent evt) {
                lookUpComponentShown(evt);
                if (ControlProcess.hasTrue("manifoldProcessed"))
                    myAGCT.saveButtonWeb.setEnabled(true);
                else
                    myAGCT.saveButtonWeb.setEnabled(false);
            }
        });
        // myManifoldPane.add(JResponsiveGeneSetter.getInstance());

        myPCAPaneName = "PCA";
        myPCAPane = new JAGCTVisualizationPane(myAGCT, myDomain, myPCAPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/mapPCA.png"))),
                JAGCTVisualizationPane.P_P);
        myPCAPane.setBackground(Color.WHITE);
        myPCAPane.addComponentListener(new ComponentAdapter() {
            public void componentShown(ComponentEvent evt) {
                lookUpComponentShown(evt);
                if (ControlProcess.hasTrue("pcaProcessed"))
                    myAGCT.saveButtonWeb.setEnabled(true);
                else
                    myAGCT.saveButtonWeb.setEnabled(false);
            }
        });

        myCorrelationPaneName = "Correlation";
        myCorrelationPane = new JAGCTVisualizationPane(myAGCT, myDomain, myCorrelationPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(
                java.net.URLClassLoader.getSystemResource("Images/mapCorrelation.gif"))), JAGCTVisualizationPane.C_C_P);
        myCorrelationPane.setBackground(Color.WHITE);
        myCorrelationPane.addComponentListener(new ComponentAdapter() {
            public void componentShown(ComponentEvent evt) {
                lookUpComponentShown(evt);
                myAGCT.saveButtonWeb.setEnabled(false);
            }
        });


        myClusteringPaneName = JMainFrameTabbedPane.THE_CLUSTERING_PANE_NAME;
        myClusteringPane = new JAGCTClusteringPane(myAGCT, myDomain, myClusteringPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(
                java.net.URLClassLoader.getSystemResource("Images/clustering.png"))));
        myClusteringPane.setBackground(Color.WHITE);
        myClusteringPane.addComponentListener(new ComponentAdapter() {
            public void componentShown(ComponentEvent evt) {
                lookUpComponentShown(evt);
                myAGCT.saveButtonWeb.setEnabled(false);
            }
        });

        // plug all tabs

        addTab(mySelectionPane.myName, mySelectionPane.myImageIcon, mySelectionPane);
        addTab(myManifoldPane.myTabName, myManifoldPane.myImageIcon, myManifoldPane);
        addTab(myPCAPane.myTabName, myPCAPane.myImageIcon, myPCAPane);
        addTab(myCorrelationPane.myTabName, myCorrelationPane.myImageIcon, myCorrelationPane);
        addTab(myClusteringPane.myName, myClusteringPane.myImageIcon, myClusteringPane);

        addChangeListener(new ChangeListener() {

            @Override
            public void stateChanged(ChangeEvent e) {
                // TODO Auto-generated method stub
                String name = getSelectedComponent().getName();
                if (JLeftIndicator.getInstance().getLoaded())
                    if (name != null && (name.equals(JAGCTVisualizationPane.P_P) || name.equals(JAGCTVisualizationPane.M_P))) {
                        ((JAGCTVisualizationPane) getSelectedComponent()).add(JLeftIndicator.getInstance(), BorderLayout.WEST);
                    }
            }
        });
        Scenario.tickSwitch();
//		AGCT.getInstance().scenarioRecordButton.doClick();
    }

    public void setAllEnabled(boolean b) {
        for (int i = 0; i < getTabCount(); i++) {
            setEnabledAt(i, b);
        }
    }

    public Rectangle getCaptureRectangle() {
        // return ;
        //
        int i = getSelectedIndex();
        Rectangle rect = null;
        // Rectangle rect =
        // ((JAGCTAbstractPane)getSelectedComponent()).getCaptureRectangle();
        if (i == 0)
            rect = mySelectionPane.getCaptureRectangle();
        else if (i == 1)
            rect = myManifoldPane.getCaptureRectangle();
        else if (i == 2)
            rect = myPCAPane.getCaptureRectangle();
        else if (i == 3)
            rect = myCorrelationPane.getCaptureRectangle();
        else if (i == 4)
            rect = myClusteringPane.getCaptureRectangle();

        if (rect == null)
            Matrix.perror("No Rectangle");

        return rect;
    }

    public void setDomain(Domain d) {
        myDomain = d;
        mySelectionPane.setDomain(d);
        myManifoldPane.setDomain(d);
        myPCAPane.setDomain(d);
        myCorrelationPane.setDomain(d);
        myClusteringPane.setDomain(d);
    }

    /**
     * 使用可能になる。
     */
    public void filterData() {
        System.out.println("JMainFrameTabbedPane.filterData()");
        myDomain.initializationFilters();
        mySelectionPane.selectionAvailable = true;
        mySelectionPane.displayComponents();
        validate();
        // setSelectedIndex(0);
        // setVisible(true);
    }

    public void toManifold() {
        myManifoldPane.displayComponents();
    }

    public void toPCA() {
        myPCAPane.displayComponents();
        System.err.println("JMainFrameTabbedPane.toPCA()");
//		for(;;);
        myCorrelationPane.displayComponents();
    }

    public void toClustering() {
        myClusteringPane.displayComponents();
    }

    private void lookUpComponentShown(ComponentEvent evt) {
        // getSelectedComponent().requestFocus();

		/*
         * String cn = evt.getComponent().getClass().getName(); if
		 * (cn.equals(mySelectionClassName)) mySelectionPane.requestFocus();
		 * else if (cn.equals(myManifoldClassName))
		 * myManifoldPane.requestFocus();
		 */
    }

    public void actionPerformed(ActionEvent e) {
    }

    @Override
    public void stateChanged(ChangeEvent e) {
        // TODO Auto-generated method stub
        String name = getSelectedComponent().getName();
        if (JLeftIndicator.getInstance().getLoaded())
            if (name != null && (name.equals(JAGCTVisualizationPane.P_P) || name.equals(JAGCTVisualizationPane.M_P))) {
                ((JAGCTVisualizationPane) getSelectedComponent()).add(JLeftIndicator.getInstance(), BorderLayout.WEST);
            }
    }


}