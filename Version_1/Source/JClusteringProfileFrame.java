/*****
 * Classes JClusteringProfile* :
 * Display the average expression profile for each cluster
 * in the current selected clutering
 *****/

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;
import javax.imageio.ImageIO;

class JClusteringProfileStamp extends JPanel implements Debuggable{
    
    public static Color LINE_COLOR = Color.red, 
	STRUCT_COLOR = Color.green,
	STDDEV_COLOR = Color.pink;

    public static double leftMarginPercent = 5.0,
	rightMarginPercent = 5.0,
	upMarginPercent = 5.0,
	downMarginPercent = 5.0;

    double [] average_profiles;
    double [] sigma_profiles;
    boolean [] undetermined_profiles;
    //says where there are no values at all;
    double [] q25_profiles;
    double [] q50_profiles;
    double [] q75_profiles;

    int myCluster, myLigand;
    //myLigand = in 0 .. SELECTED ligands - 1
    //myCluster = in 0 .. current nclusters

    String myClusterName, myLigandName;

    double xmin, xmax, ymin, ymax;

    boolean isNull;

    JClusteringProfileStamp(int nc, int nl, double [] av, double [] si, boolean [] un, double [] q25, double [] q50, double [] q75){
	myCluster = nc;
	myLigand = nl;

	myClusterName = "C" + nc;
	myLigandName = JClusteringProfileFrame.Ligand_Name_Summaries[nl];

	average_profiles = av;
	sigma_profiles = si;
	undetermined_profiles = un;
	q25_profiles = q25;
	q50_profiles = q50;
	q75_profiles = q75;
	isNull = false;

	xmin = ((Double) JClusteringProfileFrame.Time_Stamps_Summaries[nl].elementAt(0)).doubleValue();
	xmax = ((Double) JClusteringProfileFrame.Time_Stamps_Summaries[nl].elementAt(JClusteringProfileFrame.Time_Stamps_Summaries[nl].size() - 1)).doubleValue();

	int i;
	for (i=0;i<av.length;i++){
	    if ( (i==0) || (av[i] < ymin) )
		ymin = av[i];
	    if ( (i==0) || (av[i] > ymax) )
		ymax = av[i];
	}
    }

    JClusteringProfileStamp(int ng, int nl, double [] genes_profile, boolean [] genes_undetermined_profile, String genes_names){
	myCluster = ng;
	myLigand = nl;
	
	myClusterName = genes_names;
	myLigandName = JClusteringProfileFrame.Ligand_Name_Summaries[nl];

	average_profiles = genes_profile;
	undetermined_profiles = genes_undetermined_profile;
	q25_profiles = q50_profiles = q75_profiles = null;
	isNull = false;

	xmin = ((Double) JClusteringProfileFrame.Time_Stamps_Summaries[nl].elementAt(0)).doubleValue();
	xmax = ((Double) JClusteringProfileFrame.Time_Stamps_Summaries[nl].elementAt(JClusteringProfileFrame.Time_Stamps_Summaries[nl].size() - 1)).doubleValue();

	int i;
	for (i=0;i<genes_profile.length;i++){
	    if ( (i==0) || (genes_profile[i] < ymin) )
		ymin = genes_profile[i];
	    if ( (i==0) || (genes_profile[i] > ymax) )
		ymax = genes_profile[i];
	}
    }

    JClusteringProfileStamp(){
	isNull = true;
    }

    public double time(int u){
	//System.out.println("myLigand = " + myLigand + ", u = " + u);
	return ((Double) JClusteringProfileFrame.Time_Stamps_Summaries[myLigand].elementAt(u)).doubleValue();
    }

    public void structure(Graphics g){
	g.setColor(Color.white);
	g.fillRect(X(xmin, xmin, xmax),
		   Y(ymax, ymin, ymax),
		   Math.abs(X(xmax, xmin, xmax) - X(xmin, xmin, xmax)),
		   Math.abs(Y(ymin, ymin, ymax) - Y(ymax, ymin, ymax)));

	g.setColor(JClusteringProfileStamp.STRUCT_COLOR);
	g.drawLine(X(xmin, xmin, xmax),
		   Y(ymin, ymin, ymax),
		   X(xmax, xmin, xmax),
		   Y(ymin, ymin, ymax));
	g.drawLine(X(xmax, xmin, xmax),
		   Y(ymin, ymin, ymax),
		   X(xmax, xmin, xmax),
		   Y(ymax, ymin, ymax));
	g.drawLine(X(xmin, xmin, xmax),
		   Y(ymax, ymin, ymax),
		   X(xmax, xmin, xmax),
		   Y(ymax, ymin, ymax));
	g.drawLine(X(xmin, xmin, xmax),
		   Y(ymin, ymin, ymax),
		   X(xmin, xmin, xmax),
		   Y(ymax, ymin, ymax));

	g.setFont(JAGCTGraphicsPanel.SMALL_FONT);
	g.drawString(myLigandName + " " + myClusterName, X(xmin, xmin, xmax), Y(ymin, ymin, ymax));
    }

    public void paintComponent (Graphics g){
	double mix, max;
	int i1, i2;
	boolean stop = false;
	if (!isNull){
	    
	    setBackground(Color.white);
	    structure(g);

	    i1 = 0;
	    while( (i1<average_profiles.length) && (undetermined_profiles[i1]) )
		i1++;	    
	    do{
		i2 = i1 + 1;
		while( (i2<average_profiles.length) && (undetermined_profiles[i2]) )
		    i2++;
		if (i2 < average_profiles.length) {
		    g.setColor(JClusteringProfileStamp.LINE_COLOR);
		    g.drawLine(X(time(i1), xmin, xmax),
			       Y(average_profiles[i1], ymin, ymax),
			       X(time(i2), xmin, xmax),
			       Y(average_profiles[i2], ymin, ymax));

		    i1 = i2;
		}else
		    stop = true;
	    }while(!stop);
	}
    }


    public int X(double xx, double xmin, double xmax){
	double w = ( (double) getWidth() * (100.0 - leftMarginPercent - rightMarginPercent) / 100.0 );
	double delta = ( (double) getWidth() * leftMarginPercent) / 100.0;
	double v = delta + ( ( (xx - xmin) * w ) / (xmax - xmin) );
	return (int) v;
    }

    public int Y(double yy, double ymin, double ymax){
	double h = ( (double) getHeight() * (100.0 - upMarginPercent - downMarginPercent) / 100 );
	double delta = ( (double) getHeight() * upMarginPercent) / 100.0;
	double v = h + delta - ( ( (yy - ymin) * h ) / (ymax - ymin) );
	return (int) v;
    }
}

class JClusteringProfileGridPlot extends JPanel implements Debuggable{

    public static Color[] colors = {Color.red,Color.pink};
    public static int H_GAP = 1, V_GAP = 1;

    JClusteringProfileFrame myFrame;
    AGCT myAGCT;

    int nclusters;

    JClusteringProfileStamp [][] all_stamps;
    boolean clusterVSgenes;
    //true -> we plot cluster profiles, else gene profiles

    JClusteringProfileGridPlot(JClusteringProfileFrame jf, AGCT ma){
	int i;

	myFrame = jf;
	myAGCT = ma;
	all_stamps = null;
	clusterVSgenes = true;
    }
    
    public void switchToClustering(){
	clusterVSgenes = true;
    }

    public void initGeneProfiles(Vector allGenes){
	//components = Integers to map using selectedGeneToGeneIndex
	if ( (allGenes == null) || (allGenes.size()<1) )
	    Matrix.perror("JClusteringProfileFrame.class :: no gene to plot");
	int ng = allGenes.size(), numgen, index;
	nclusters = ng;

	Gene gg;
	int i, j, k;
	double [][][] genes_profile;
	// X = #gene, Y = #ligand, Z = #time stamp, val = profile
	boolean [][][] genes_undetermined_profile;
	//idem but says if there does not exist some value for the (X, Y, Z) of the profile at hand;
	String [] genes_name;
	double nb_values [][][];
	Domain d = myAGCT.myDomain;

	genes_profile = new double [ng][][];
	nb_values = new double [ng][][];
	genes_undetermined_profile = new boolean [ng][][];
	genes_name = new String [ng];
	for (i=0;i<ng;i++){
	    genes_profile[i] = new double [d.numberSelectedLigands][];
	    nb_values[i] = new double [d.numberSelectedLigands][];
	    genes_undetermined_profile[i] = new boolean [d.numberSelectedLigands][];
	    numgen = ((Integer) allGenes.elementAt(i)).intValue();
	    gg = (Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[numgen]);
	    genes_name[i] = gg.name + "|" + gg.asciiName;
	}

	index = 0;
	for (i=0;i<d.numberLigands;i++){
	    if (d.selectedLigand[i]){
		for (j=0;j<ng;j++){
		    genes_profile[j][index] = new double [( (Integer) d.numberTimes.elementAt(i)).intValue()];
		    nb_values[j][index] = new double [( (Integer) d.numberTimes.elementAt(i)).intValue()];
		    genes_undetermined_profile[j][index] = new boolean [( (Integer) d.numberTimes.elementAt(i)).intValue()];
		}
		index++;
	    }
	}
	for (i=0;i<ng;i++){
	    for (j=0;j<genes_profile[i].length;j++){
		for (k=0;k<genes_profile[i][j].length;k++){
			genes_profile[i][j][k] = 0.0;
			nb_values[i][j][k] = 0.0;
		}
	    }
	}

	for (i=0;i<ng;i++){
	    numgen = ((Integer) allGenes.elementAt(i)).intValue();
	    gg = (Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[numgen]);
	    
	    //System.out.println(g.rawCoordinates.length + " " + g.undeterminedCoordinates.length);
	    
	    index = 0;
	    for (j=0;j<d.numberLigands;j++){
		if (d.selectedLigand[j]){
		    for (k=0;k<nb_values[i][index].length;k++){
			if (!gg.undeterminedCoordinates[index][k]){
			    genes_profile[i][index][k] = gg.rawCoordinates[index][k];
			    nb_values[i][index][k] += 1.0;
			}
		    }
		    index++;
		}
	    }
	}
	
	for (i=0;i<ng;i++){
	    for (j=0;j<genes_profile[i].length;j++){
		for (k=0;k<genes_profile[i][j].length;k++){
		    if (nb_values[i][j][k] > 0.0)
			genes_undetermined_profile[i][j][k] = false;
		    else
			genes_undetermined_profile[i][j][k] = true;
		}
	    }
	}
	
	all_stamps = new JClusteringProfileStamp[ng][];
	for (i=0;i<ng;i++){
	    all_stamps[i] = new JClusteringProfileStamp[genes_profile[i].length];
	    for (j=0;j<genes_profile[i].length;j++){
		all_stamps[i][j] = new JClusteringProfileStamp(i, j, genes_profile[i][j], genes_undetermined_profile[i][j], genes_name[i]);
	    }
	}
    }

    public void initAverageProfiles(boolean quartiles){
	int np, i, j, k, l, index, nclu, l25, l50, l75;
	Clustering cc;
	Domain d = myAGCT.myDomain;
	Gene g;

	double nb_values [][][];
	double ex2 [][][];

	double [][][] all_curves_average_profile;
	// X = #cluster, Y = #ligand, Z = #time stamp, val = average profile
	double [][][] all_curves_sigma_profile;
	//idem for sigma
	boolean [][][] all_curves_undetermined_profile;
	//idem but says if there does not exist some value for the (X, Y, Z) of the profile at hand;

	double [][][] all_curves_Q25;
	double [][][] all_curves_Q50;
	double [][][] all_curves_Q75;
	//idem but for the 25% and 75% quartile
	double [] ddum;

	Vector clusterElements, element, clusterValues;

	all_stamps = null;

	if ( (myFrame.listPlots.getSelectedIndex()>=0) && (myFrame.onePlot) ){
	    np = ((Integer) myFrame.keysForCombo.elementAt(myFrame.listPlots.getSelectedIndex())).intValue();
	    cc = (Clustering) myAGCT.allClusterings.elementAt(np);
	    nclusters = cc.nclusters;

	    all_curves_average_profile = new double [nclusters][][];
	    all_curves_sigma_profile = new double [nclusters][][];
	    all_curves_undetermined_profile = new boolean [nclusters][][];
	    all_curves_Q25 = new double [nclusters][][];
	    all_curves_Q50 = new double [nclusters][][];
	    all_curves_Q75 = new double [nclusters][][];
	    ex2 = new double [nclusters][][];
	    nb_values = new double [nclusters][][];
	    for (i=0;i<nclusters;i++){
		all_curves_average_profile[i] = new double [d.numberSelectedLigands][];
		all_curves_sigma_profile[i] = new double [d.numberSelectedLigands][];
		all_curves_undetermined_profile[i] = new boolean [d.numberSelectedLigands][];
		all_curves_Q25[i] = new double [d.numberSelectedLigands][];
		all_curves_Q50[i] = new double [d.numberSelectedLigands][];
		all_curves_Q75[i] = new double [d.numberSelectedLigands][];
		ex2[i] = new double [d.numberSelectedLigands][];
		nb_values[i] = new double [d.numberSelectedLigands][];		
	    }
	    index = 0;
	    for (i=0;i<d.numberLigands;i++){
		if (d.selectedLigand[i]){
		    for (j=0;j<nclusters;j++){
			all_curves_average_profile[j][index] = new double [( (Integer) d.numberTimes.elementAt(i)).intValue()];
			all_curves_sigma_profile[j][index] = new double [( (Integer) d.numberTimes.elementAt(i)).intValue()];
			all_curves_undetermined_profile[j][index] = new boolean [( (Integer) d.numberTimes.elementAt(i)).intValue()];
			all_curves_Q25[j][index] = new double [( (Integer) d.numberTimes.elementAt(i)).intValue()];
			all_curves_Q50[j][index] = new double [( (Integer) d.numberTimes.elementAt(i)).intValue()];
			all_curves_Q75[j][index] = new double [( (Integer) d.numberTimes.elementAt(i)).intValue()];
			ex2[j][index] = new double [( (Integer) d.numberTimes.elementAt(i)).intValue()];
			nb_values[j][index] = new double [( (Integer) d.numberTimes.elementAt(i)).intValue()];
		    }
		    index++;
		}
	    }
	    for (i=0;i<nclusters;i++){
		for (j=0;j<all_curves_average_profile[i].length;j++){
		    for (k=0;k<all_curves_average_profile[i][j].length;k++){
			all_curves_average_profile[i][j][k] = 0.0;
			all_curves_sigma_profile[i][j][k] = 0.0;
			all_curves_Q25[i][j][k] = 0.0;
			all_curves_Q50[i][j][k] = 0.0;
			all_curves_Q75[i][j][k] = 0.0;
			ex2[i][j][k] = 0.0;
			nb_values[i][j][k] = 0.0;
		    }
		}
	    }

	    for (i=0;i<d.numberSelectedGenes;i++){
		g = (Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[i]);

		//System.out.println(g.rawCoordinates.length + " " + g.undeterminedCoordinates.length);

		index = 0;
		nclu = cc.myClusteringAlgorithm.selectedGeneToHardCluster[i];
		for (j=0;j<d.numberLigands;j++){
		    if (d.selectedLigand[j]){
			for (k=0;k<nb_values[nclu][index].length;k++){
			    if (!g.undeterminedCoordinates[index][k]){
				//System.out.println(" ex: " + i + " :: " + g.rawCoordinates[j][k]);
				all_curves_average_profile[nclu][index][k] += g.rawCoordinates[index][k];
				ex2[nclu][index][k] += (g.rawCoordinates[index][k]*g.rawCoordinates[index][k]);
				nb_values[nclu][index][k] += 1.0;
			    }
			}
			index++;
		    }
		}
	    }

	    for (i=0;i<nclusters;i++){
		for (j=0;j<all_curves_average_profile[i].length;j++){
		    for (k=0;k<all_curves_average_profile[i][j].length;k++){
			if (nb_values[i][j][k] > 0.0){
			    all_curves_average_profile[i][j][k] /= nb_values[i][j][k];
			    ex2[i][j][k] /= nb_values[i][j][k];
			    all_curves_undetermined_profile[i][j][k] = false;
			}else
			    all_curves_undetermined_profile[i][j][k] = true;
		    }
		}
	    }
	    for (i=0;i<nclusters;i++){
		for (j=0;j<all_curves_average_profile[i].length;j++){
		    for (k=0;k<all_curves_average_profile[i][j].length;k++){
			if (nb_values[i][j][k] > 0.0)
			    all_curves_sigma_profile[i][j][k] = Math.sqrt(ex2[i][j][k] - (all_curves_average_profile[i][j][k] * all_curves_average_profile[i][j][k]));
		    }
		}
	    }

	    // Quartiles
	    if (quartiles){
		clusterElements = new Vector();
		for (i=0;i<nclusters;i++)
		    clusterElements.addElement(new Vector());
		for (i=0;i<d.numberSelectedGenes;i++){
		    element = (Vector) clusterElements.elementAt(cc.myClusteringAlgorithm.selectedGeneToHardCluster[i]);
		    element.add(new Integer(i));
		}
		for (i=0;i<nclusters;i++){
		    index = 0;
		    element = (Vector) clusterElements.elementAt(i);
		    for (j=0;j<d.numberLigands;j++){
			if (d.selectedLigand[j]){
			    for (k=0;k<nb_values[i][index].length;k++){
				clusterValues = new Vector();
				for (l=0;l<element.size();l++){
				    g = (Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[((Integer) element.elementAt(l)).intValue()]);
				    if (!g.undeterminedCoordinates[index][k]){
					clusterValues.addElement(new Double(g.rawCoordinates[index][k]));
				    }
				}
				ddum = new double [clusterValues.size()];
				for (l=0;l<clusterValues.size();l++)
				    ddum[l] = ( (Double) clusterValues.elementAt(l) ).doubleValue();
				QuickSort.quicksort(ddum);

				l25 = ddum.length / 4;
				if (l25>0)
				    l25--;
				all_curves_Q25[i][j][k] = ddum[l25];

				l50 = ddum.length / 2;
				if (l50>0)
				    l50--;
				all_curves_Q50[i][j][k] = ddum[l50];

				l75 = ddum.length - (ddum.length / 4);
				if (l75>0)
				    l75--;
				all_curves_Q75[i][j][k] = ddum[l75];
			    }
			    index++;
			}
		    }
		}
	    }
	    // And finally... the stamps !
	    
	    all_stamps = new JClusteringProfileStamp[nclusters][];
	    for (i=0;i<nclusters;i++){
		all_stamps[i] = new JClusteringProfileStamp[all_curves_average_profile[i].length];
		for (j=0;j<all_curves_average_profile[i].length;j++){
		    if (quartiles)
			all_stamps[i][j] = new JClusteringProfileStamp(i, j, all_curves_average_profile[i][j], all_curves_sigma_profile[i][j], all_curves_undetermined_profile[i][j], all_curves_Q25[i][j], all_curves_Q50[i][j], all_curves_Q75[i][j]);
		    else
			all_stamps[i][j] = new JClusteringProfileStamp(i, j, all_curves_average_profile[i][j], all_curves_sigma_profile[i][j], all_curves_undetermined_profile[i][j], null, null, null);
		}
	    }
	}

	nb_values = ex2 = null;
	all_curves_average_profile = all_curves_sigma_profile = null;
	all_curves_undetermined_profile = null;
	clusterElements = null;
    }

    public void makeUpdateLayout(){
	int i, j;
	JLabel jl;
	JClusteringProfileStamp js;

	removeAll();
	setLayout(new GridLayout(nclusters, myAGCT.myDomain.numberSelectedLigands));

	for (i=0;i<nclusters;i++){
	    for (j=0;j<myAGCT.myDomain.numberSelectedLigands;j++){
		add(all_stamps[i][j]);
	    }
	}
	repaint();
    }

    public void repaintEverything(){
	int i, j;
	for (i=0;i<nclusters;i++){
	    for (j=0;j<myAGCT.myDomain.numberSelectedLigands;j++){
		all_stamps[i][j].repaint();
		revalidate();
	    }
	}
    }

    public void paintComponent (Graphics g) {
	if ( (!clusterVSgenes) && ( (myFrame.listGenes == null) || (myFrame.listGenes.size() < 1) ) ){
	    removeAll();
	    myFrame.updateBorder("");
	    g.drawString("No visualisation possible: no gene buffered", 40, 40);
	}else if ( (!clusterVSgenes) || ( (myFrame.onePlot) && (myFrame.plotAuthorized) ) ){
	    setBackground(Color.white);

	    if (myFrame.moving){
		makeUpdateLayout();
		myFrame.moving = false;
	    }
	    repaintEverything();
	}
    }

    public Rectangle getCaptureRectangle(){
	Rectangle bounds = getBounds();
        bounds.setLocation(getLocationOnScreen());
	return bounds;
    }
}

class JClusteringProfileFrame extends JFrame{
    
    public static String Default_Border = "(no clustering plotted)";
    public static String Default_Combo[] = {"(no clustering)"};

    public static Vector [] Time_Stamps_Summaries;
    //We only store time stamps for selected ligands (saves time)

    public static String [] Ligand_Name_Summaries;
    //idem for ligand names

    String borderString;
    JComboBox listPlots;
    JClusteringProfileGridPlot graphPanel;
    JButton goButton, cameraButton, saveProfileButton;
    Box mainBox, selectionBox;

    Vector keysForCombo;
    //contains the vector of Integers = index in allClusterings of the cluster to be plotted
    
    AGCT myAGCT;
    boolean onePlot, plotAuthorized, moving;

    Vector listGenes;
    //list of genes whose profile can be plotted
    int maxNumberGenes;

    public static void fillTime_Stamps_Summaries(Domain d){
	int i,j,index = 0;
	double cd;
	Vector current;
	Time_Stamps_Summaries = new Vector [d.numberSelectedLigands];
	Ligand_Name_Summaries = new String [d.numberSelectedLigands];
	for (i=0;i<d.numberLigands;i++){
	    if (d.selectedLigand[i]){
		current = new Vector();
		for (j=0;j<( (Integer) d.numberTimes.elementAt(i)).intValue();j++){
		    cd = ((Double) ((Vector) d.domainTimes.elementAt(i)).elementAt(j+1)).doubleValue();
		    current.addElement(new Double(cd));
		}
		Time_Stamps_Summaries [index] = current;
		Ligand_Name_Summaries [index] = ((String) ((Vector) d.domainTimes.elementAt(i)).elementAt(0));
		index++;
	    }
	}
    }

    public boolean atLeastTwoGenes(){
	if (listGenes == null)
	    return false;
	if (listGenes.size() < 2)
	    return false;
	return true;
    }

    public void addGene(int num){
	if (listGenes == null)
	    listGenes = new Vector();
	if ( (maxNumberGenes < AGCT.DEFAULT_Max_Number_Genes) || (listGenes.size() <= AGCT.DEFAULT_Max_Number_Genes) ){
	    maxNumberGenes = AGCT.DEFAULT_Max_Number_Genes;
	}else{
	    Vector trim = new Vector();
	    int i;
	    for (i=0;i<AGCT.DEFAULT_Max_Number_Genes;i++)
		trim.addElement(listGenes.elementAt(listGenes.size()-AGCT.DEFAULT_Max_Number_Genes+i));
	    listGenes = trim;
	    trim = null;
	    maxNumberGenes = AGCT.DEFAULT_Max_Number_Genes;
	}
	if (listGenes.size() >= maxNumberGenes)
	    listGenes.removeElementAt(0);
	if (!listGenes.contains(new Integer(num)))
	    listGenes.addElement(new Integer(num));
    }

    public void flushGenes(){
	listGenes = new Vector();
    }

    public void flushAll(){
	listPlots.removeAllItems();
	listPlots.addItem(JClusteringProfileFrame.Default_Combo);
	keysForCombo = null;
	listPlots.setEnabled(false);
	onePlot = plotAuthorized  = false;
	flushGenes();
    }

    public void captureAndSave(){
	Rectangle rect = graphPanel.getCaptureRectangle();
	BufferedImage screencapture = null, framecapture = null;
	try{
	    framecapture = new Robot().createScreenCapture(rect);
	}catch(AWTException e){}
	
	JFileChooser chooser = new JFileChooser();
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("png");
	filter.setDescription(".png Files");
	chooser.setFileFilter(filter);
	chooser.setApproveButtonText("Save");
	int returnVal = chooser.showSaveDialog(this);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    try{
		ImageIO.write(framecapture, "png", chooser.getSelectedFile());
	    }catch(IOException e){}
	    myAGCT.myInformationFrame.setText("Saving visible clustering profile pane to file: " +
				       chooser.getSelectedFile().getName());
	    ControlProcess.put("frameCapturedClusteringProfile",true);
	}
    }

    public void displayInfo(){
	cameraButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/camera.png"))));
	cameraButton.setToolTipText("capture the visible cluster profiles plot");
        cameraButton.setActionCommand("capture");
	cameraButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){
		    captureAndSave();
		}
	    });
	cameraButton.setEnabled(false);

	goButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/go.png"))));
	goButton.setToolTipText("run / update cluster profiles");
        goButton.setActionCommand("clustering");
	goButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){
		    if (listPlots.getSelectedIndex()>=0){
			goPlotClusters();
		    }
		}
	    });
	goButton.setEnabled(false);

	saveProfileButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/diskProfile.png"))));
	saveProfileButton.setToolTipText("save current profiles");
        saveProfileButton.setActionCommand("save_profiles");
	saveProfileButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){
		    if (listPlots.getSelectedIndex()>=0){
			graphPanel.initAverageProfiles(true);
			myAGCT.requestSave_ClusterProfile();
		    }
		}
	    });
	saveProfileButton.setEnabled(false);


	graphPanel = new JClusteringProfileGridPlot(this, myAGCT);
	graphPanel.setBorder(BorderFactory.createTitledBorder(JClusteringProfileFrame.Default_Border));

	listPlots = new JComboBox(JClusteringProfileFrame.Default_Combo);
	listPlots.setToolTipText("select which clustering data to plot");
	listPlots.setSelectedIndex(0);
	listPlots.setEnabled(false);

	selectionBox = Box.createHorizontalBox();
	//selectionBox.add(cameraButton);
	selectionBox.add(listPlots);
	selectionBox.add(goButton);
	selectionBox.add(saveProfileButton);

        Container pane = getContentPane();
        pane.setLayout(new BorderLayout());
	pane.add(selectionBox, BorderLayout.NORTH);
	pane.add(graphPanel, BorderLayout.CENTER);

	addWindowListener(new FermetureListener("Closing AGCT's ClusteringProfileFrame\n"));

	Image img = Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/information.png"));
	setIconImage( img );
    }

    public void goPlotClusters(){
	graphPanel.clusterVSgenes = true;
	graphPanel.initAverageProfiles(false);
	onePlot = moving = true;
	updateBorder(currentClusterString());
	repaint();
    }

    public void goPlotGenes(){
	graphPanel.clusterVSgenes = false;
	moving = true;
	graphPanel.initGeneProfiles(listGenes);
	updateBorder("");
	repaint();
    }

    public void updateBorder(String t){
	if (graphPanel.clusterVSgenes)
	    graphPanel.setBorder(BorderFactory.createTitledBorder(t));
	else{
	    String tt = "";
	    int i, numgen;
	    Gene gg;
	    Domain d = myAGCT.myDomain;
	    if ( (listGenes != null) && (listGenes.size()>0) ){
		tt = "Profiles for genes: ";
		for (i=0;i<listGenes.size();i++){
		    numgen = ((Integer) listGenes.elementAt(i)).intValue();
		    gg = (Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[numgen]);

		    tt += gg.asciiName + " ";
		}
		graphPanel.setBorder(BorderFactory.createTitledBorder(tt));
	    }else
		graphPanel.setBorder(BorderFactory.createTitledBorder(tt));
	}
    }

    public void addClusteringLF(Clustering cc, int indexClustering){

	if (myAGCT.allClusterings != null){
	    if (onePlot == false){
		listPlots.removeAllItems();
		onePlot = plotAuthorized = true;
		keysForCombo = new Vector();

		listPlots.setEnabled(true);
		cameraButton.setEnabled(true);
		goButton.setEnabled(true);
		saveProfileButton.setEnabled(true);
	    }
	
	    listPlots.addItem("C#" + indexClustering + " (" + cc.myClusteringAlgorithm.myReferenceName + ")");
	    keysForCombo.addElement(new Integer(indexClustering));

	    if (!isVisible()){
		setVisible(true);
		myAGCT.myFrame.toFront();
	    }
	}
    }

    public String currentClusterString(){
	int id = listPlots.getSelectedIndex();
	String rn = ((Clustering) myAGCT.allClusterings.elementAt(id)).myClusteringAlgorithm.myReferenceName;

	return "C#" + id + " (" + rn + ")";
    }
    
    JClusteringProfileFrame(String name, AGCT c){
	super(name);
	setSize(AGCT.WindowWidth,600);
	myAGCT = c;
	displayInfo();
	setVisible(false);
	setResizable(true);
	onePlot = plotAuthorized = false;
	moving = false;
	listGenes = new Vector();
	maxNumberGenes = AGCT.DEFAULT_Max_Number_Genes;
    }
}