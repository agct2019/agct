import java.awt.event.*;
import java.awt.Color;
import java.awt.KeyboardFocusManager;
import java.awt.Component;
import java.awt.Rectangle;
import java.awt.Toolkit;
import javax.swing.*;

class JMainFrameTabbedPane extends JTabbedPane implements ActionListener, Debuggable{
    AGCT myAGCT;
    Domain myDomain;

    public static String THE_SELECTION_PANE_NAME = "Selection",
	THE_MANIFOLD_PANE_NAME = "Manifold",
	THE_PCA_PANE_NAME = "PCA",
	THE_CORRELATION_PANE_NAME = "Correlation",
	THE_CLUSTERING_PANE_NAME = "Clustering";

    JAGCTSelectionPane mySelectionPane;
    JAGCTVisualizationPane myManifoldPane, myPCAPane, myCorrelationPane;
    JAGCTClusteringPane myClusteringPane;

    String mySelectionPaneName, myManifoldPaneName, myPCAPaneName, myCorrelationPaneName, myClusteringPaneName;
    String mySelectionClassName, myManifoldClassName, myPCAClassName, myCorrelationClassName, myClusteringClassName;

    JMainFrameTabbedPane(AGCT a, Domain d){
	super();

	myAGCT = a;
	myDomain = d;

	mySelectionPaneName = JMainFrameTabbedPane.THE_SELECTION_PANE_NAME;
	mySelectionClassName = "JAGCTSelectionPane";
	mySelectionPane = new JAGCTSelectionPane(myAGCT, myDomain, mySelectionPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/book_open.png"))));
	mySelectionPane.setBackground(Color.WHITE);
	mySelectionPane.addComponentListener(new ComponentAdapter() {
		    public void componentShown(ComponentEvent evt) {
			lookUpComponentShown(evt);
			myAGCT.saveButtonWeb.setEnabled(false);
		    }
		});

	myManifoldPaneName = JMainFrameTabbedPane.THE_MANIFOLD_PANE_NAME;
	myManifoldClassName = "JAGCTVisualizationPane";
	myManifoldPane = new JAGCTVisualizationPane(myAGCT, myDomain, myManifoldPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/mapManifold.png"))), JAGCTVisualizationPane.M_P);
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

	myPCAPaneName = "PCA";
	myPCAClassName = JMainFrameTabbedPane.THE_PCA_PANE_NAME;
	myPCAPane = new JAGCTVisualizationPane(myAGCT, myDomain, myPCAPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/mapPCA.png"))), JAGCTVisualizationPane.P_P);
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
	myCorrelationClassName = JMainFrameTabbedPane.THE_CORRELATION_PANE_NAME;
	myCorrelationPane = new JAGCTVisualizationPane(myAGCT, myDomain, myCorrelationPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/mapCorrelation.gif"))), JAGCTVisualizationPane.C_C_P);
	myCorrelationPane.setBackground(Color.WHITE);
	myCorrelationPane.addComponentListener(new ComponentAdapter() {
		    public void componentShown(ComponentEvent evt) {
			lookUpComponentShown(evt);
			myAGCT.saveButtonWeb.setEnabled(false);
		    }
		});

	myClusteringPaneName = JMainFrameTabbedPane.THE_CLUSTERING_PANE_NAME;
	myClusteringClassName = "JAGCTClusteringPane";
	myClusteringPane = new JAGCTClusteringPane(myAGCT, myDomain, myClusteringPaneName, new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/clustering.png"))));
	myClusteringPane.setBackground(Color.WHITE);
	myClusteringPane.addComponentListener(new ComponentAdapter() {
		    public void componentShown(ComponentEvent evt) {
			lookUpComponentShown(evt);
			myAGCT.saveButtonWeb.setEnabled(false);
		    }
		});
	

	//plug all tabs

	addTab(mySelectionPane.myName, mySelectionPane.myImageIcon, mySelectionPane);
	addTab(myManifoldPane.myTabName, myManifoldPane.myImageIcon, myManifoldPane);
	addTab(myPCAPane.myTabName, myPCAPane.myImageIcon, myPCAPane);
	addTab(myCorrelationPane.myTabName, myCorrelationPane.myImageIcon, myCorrelationPane);
	addTab(myClusteringPane.myName, myClusteringPane.myImageIcon, myClusteringPane);
    }

    public void setAllEnabled(boolean b){
	int i;
	for (i=0;i<getTabCount();i++){
	    setEnabledAt(i, b);
	}
    }

    public Rectangle getCaptureRectangle(){
	int i = getSelectedIndex();
	Rectangle rect = null;
	if (i==0)
	    rect = mySelectionPane.getCaptureRectangle();
	else if (i==1)
	    rect = myManifoldPane.getCaptureRectangle();
	else if (i==2)
	    rect = myPCAPane.getCaptureRectangle();
	else if (i==3)
	    rect = myCorrelationPane.getCaptureRectangle();
	else if (i==4)
	    rect = myClusteringPane.getCaptureRectangle();

	if (rect == null)
	    Matrix.perror("No Rectangle");

	return rect;
    }


    public void setDomain(Domain d){
	myDomain = d;
	mySelectionPane.setDomain(d);
	myManifoldPane.setDomain(d);
	myPCAPane.setDomain(d);
	myCorrelationPane.setDomain(d);
	myClusteringPane.setDomain(d);
    }

    public void filterData(){
	mySelectionPane.filterData();
	mySelectionPane.displayComponents();
	validate();
	//setSelectedIndex(0);
	//setVisible(true);
    }

    public void toManifold(){
	myManifoldPane.displayComponents();
    }

    public void toPCA(){
	myPCAPane.displayComponents();
	myCorrelationPane.displayComponents();
    }

    public void toClustering(){
	myClusteringPane.displayComponents();
    }

    private void lookUpComponentShown(ComponentEvent evt) {
	//getSelectedComponent().requestFocus();

	/*String cn = evt.getComponent().getClass().getName();
	if (cn.equals(mySelectionClassName))
	mySelectionPane.requestFocus();
	else if (cn.equals(myManifoldClassName))
	myManifoldPane.requestFocus();*/
    }

    public void actionPerformed (ActionEvent e) {
    }
}