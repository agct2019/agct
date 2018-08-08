import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import java.util.*;
import java.text.DecimalFormat;
import javax.swing.Box.Filler;

class JAGCTClusteringPane extends JPanel implements Debuggable {

    static String [] allChoices = { "None",
				    "Complete Positive Factorization",
				    "Bregman (Iterative) K-Means",
				    "Expectation Maximization",
				    "Hierarchical Clustering",
				    "Affinity Propagation",
				    "Non-Negative Matrix Factorization"};

    static String [] defaultStrings = { "",
					"CP -Point Manifold3D -K 3",
					"KM -Point Manifold3D -Iterative 0 -Generator L2 -Seed Forgy -I 0 -K 3",
					"EM -Point Manifold3D -Seed Forgy -I 0 -K 3",
					"HC -Point Manifold3D -Distance Ward -K 3",
					"AP -Point Manifold3D -P -2 -B 2 -I 100",
					"NM -Point Manifold3D -M MS -U L2 -I 0 -K 2"};

    static String rs = "\n\n--------------------\nRemarks:\n\n"
	+ "*1: instead of flag -K X, you can use flag\n-L X Y\nto indicate loop from X to Y clusters\n(ex: -L 4 8 makes runs for K=4, 5, 6, 7 and 8 clusters)\n\n"
	+ "*2: you can also use flag -R X\nwhere X >= 1 integer, to indicate that\nthe previous specifications\nhave to be repeated X times\n--------------------";

    static String rs2 = "\n\n--------------------\nRemarks:\n\n"
	+ "*1: instead of flag -P X, you can use flag\n-L X Y (for X and Y > 0)\nto indicate loop from X% to Y% quantiles\n(ex: -L 1 4 makes runs for p=1, 2, 3, and 4%)\n\n"
	+ "*2: you can also use flag -R X\nwhere X >= 1 integer, to indicate that\nthe previous specifications\nhave to be repeated X times\n--------------------";
    
    static String rs3 = "\n\n(*) If the number of clusters is too large, switch to full Manifold coordinates.\n(Rule from Kim & Tidor, NAR 2008: k(m+n)<mn)";

    static String [] searchKeywords = {"-Help", "-Size", "-Type", "-Best", "-BestRatio"};
    static String [] helpSearchStrings = {"Prints this help\n\n * Search is multicriteria (the following can be combined)", 
					  "-Size XX (integer) restricts search to clusterings having exactly XX clusters",
					  "-Type YY (in {CP, KM, EM, HC, AP}) restricts search to clusterings of type YY",
					  "-Best XX (integer in 1, 2, ..., 100) restrict search to clustering belonging to the XX% of the best clustering objective functions",
					  "-BestRatio XX (integer in 1, 2, ..., 100) does the same as -Best *but* for normalized clustering function (/ by the number of clusters)"};
    static String exampleSearchString = " * Example:\nWriting:\n-Size 3 -Type KM\nreturns the indexes of clusterings of Bregman K-Means type, and having 3 clusters";

    static String [] helpStrings = { "",
				     "Complete Positive factorization heuristic\n\n"
				     + "-Point :\n Manifold -> use manifold points\n Manifold3D -> use the CURRENT 3D view of manifold points\n PCA -> use the CURRENT 3D view of principal components\n\n (Point type is useful ONLY to compute distortions in CP)\n\n"
				     + "-K : number of clusters\n\n" + rs,

				     "Bregman K-means\n\n"
				     + "-Point :\n Manifold -> use manifold points\n Manifold3D -> use the CURRENT 3D view of manifold points\n PCA -> use the CURRENT 3D view of principal components\n\n"
				     + "-Generator : type of generator used for the Bregman divergence\n L2 -> use squared Euclidean distance (regular K-means, no constraint)\n KL -> use Kullback-Leibler divergence (only with non negative entries)\n IS -> Itakura-Saito divergence (only with non negative entries)\n AM -> use Amari alpha divergence (only with non-negative entries)\n PN -> use p-norm (no constraint)\n\n Note :\n with AM you can specify flag -G [real value] to specify -1.0<alpha<1.0 (default = 0.0)\n with PN you can specify flag -G [real value] to specify p>1.0 (default = 2.0)\n\n"
				     + "-Seed : type of seeding of the initial centers\n Forgy -> use Forgy initialization\n Arthur -> use Arthur-Vassilvitskii initialization (L22 divergence)\n Bregman -> use Arthur-Vassilvitskii initialization with the *same* generator as in -Generator\n\n" 
				     + "-K : number of clusters\n\n" 
				     + "-I : max number of iterations\n  0 -> do not use max number\n >0 -> if the number of iterations exceeds this bound, stop\n\n" 
				     + "-Iterative : Bregman Iterative K-Means alternating with Bregman K-Means\n  0 -> Does not perform Bregman Iterative K-Means\n >0 -> Alternates K-Means and (this number of iterations of Iterative Bregman K-Means) until no significant change" + rs,

				     "Expectation Maximization\n\n"
				     + "-Point :\n Manifold -> use manifold points\n Manifold3D -> use the CURRENT 3D view of manifold points\n PCA -> use the CURRENT 3D view of principal components\n\n"
				     + "-Seed : type of seeding of the initial centers\n Forgy -> use Forgy initialization\n Arthur -> use Arthur-Vassilvitskii initialization (L22 divergence)\n\n" 
				     + "-K : number of clusters\n\n" 
				     + "-I : max number of iterations\n  0 -> do not use max number\n >0 -> if the number of iterations exceeds this bound, stop" + rs,

				     "Hierarchical Clustering\n\n"
				     + "-Point :\n Manifold -> use manifold points\n Manifold3D -> use the CURRENT 3D view of manifold points\n PCA -> use the CURRENT 3D view of principal components\n\n"
				     + "-K : number of clusters\n\n"
				     + "-Distance :\n Ward -> use Ward criterion for the aggregation\n Single_Linkage -> use single linkage criterion for the aggregation\n\n" + rs,

				     "Affinity Propagation\n\n"
				     + "-Point :\n Manifold -> use manifold points\n Manifold3D -> use the CURRENT 3D view of manifold points\n PCA -> use the CURRENT 3D view of principal components\n\n"
				     + "-P : preferences\n -2 -> use the median of the similarities\n -1 -> use the minimum of the similarities\n 1 - 99 -> use the p% quantile, which splits similarities in p% - (100-p)%\n\n"
				     + "-B : Bregman divergence with generator p-norm divergence, p=B \n (real > 1.0)\n\n" 
				     + "-I : max number of iterations\n  0 -> do not use max number\n >0 -> if the number of iterations exceeds this bound, stop" + rs2,

				     "Non-Negative Matrix Factorization\n\n"
				     + "-Point :\n Manifold -> use manifold points\n Manifold3D -> use the CURRENT 3D view of manifold points\n PCA -> use the CURRENT 3D view of principal components\n\n"
				     + "-M : method to initialize the positive matrix\n KT -> method of Kim and Tidor, Nucleic Acids Research 2008\n MS -> Shifts data by the Min\n\n"
				     + "-U : type of error to be minimized during NMF\n L2 -> squared Euclidean distance\n KL -> Kullback-Leibler divergence\n\n"
				     + "-K : number of clusters, must be < " + AGCT.Number_Of_Manifold_Components + " (*)\n\n" 
				     + "-I : max number of iterations\n  0 -> do not use max number\n >0 -> if the number of iterations exceeds this bound, stop" + rs + rs3};
 
    AGCT myAGCT;
    Domain myDomain;
    String myName;
    ImageIcon myImageIcon;

    JTextField stringChoice, stringSearch;
    JButton goButton, goSearchButton;
    JComboBox algoSelect;
    JTextArea helpPane;
    JTextArea listClusteringPane;

    boolean clusteringAvailable;
    
    JAGCTClusteringPane(AGCT a, Domain d, String name, ImageIcon ic){
	myName = name;
	myAGCT = a;
	myDomain = d;
	myImageIcon = ic;

	clusteringAvailable = false;
    }


    public Rectangle getCaptureRectangle(){
	Rectangle bounds = getBounds();
        bounds.setLocation(getLocationOnScreen());
	return bounds;
    }


    public void setDomain(Domain d){
	myDomain = d;
    }

    public void toTrash(){
	clusteringAvailable = false;
	removeAll();
    }


    public void paintComponent (Graphics g) {
	if (clusteringAvailable == false) {
	    super.paintComponent(g);
	    g.drawString("No clustering algorithm available: no data processed", myAGCT.myTabbedPane.mySelectionPane.xMargin, myAGCT.myTabbedPane.mySelectionPane.yMargin);
	}
    }

    public void displayComponents(){
	Box selectionBox, stringBox, informationBox, clusterBox, superBox, searchBox, generalBox;
	
	JPanel upPanel;

	algoSelect = new JComboBox(JAGCTClusteringPane.allChoices);
	algoSelect.setToolTipText("select the clustering algorithm");
	algoSelect.setSelectedIndex(0);
	algoSelect.setEnabled(false);
	algoSelect.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){
		    changeAlgo(algoSelect.getSelectedIndex());
		}
	    });

	goButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/go.png"))));
	goButton.setToolTipText("run clustering");
        goButton.setActionCommand("clustering");
	goButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){
		    myAGCT.clustering(stringChoice.getText());
		    Scenario.add("JAGCTClusteringPane_Run_Clustering",stringChoice.getText());
		}
	    });
	goButton.setEnabled(false);

	goSearchButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/magnifier.png"))));
	goSearchButton.setToolTipText("search");
        goSearchButton.setActionCommand("search");
	goSearchButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){
		    displaySearch(stringSearch.getText());
		}
	    });
	goSearchButton.setEnabled(false);

	selectionBox = Box.createHorizontalBox();
	selectionBox.setBorder(BorderFactory.createTitledBorder("Clustering type"));
	selectionBox.add(algoSelect);
	selectionBox.add(goButton);

	stringSearch = new JTextField(JAGCTClusteringPane.searchKeywords[0], 30);
	stringSearch.setEditable(false);
	searchBox = Box.createHorizontalBox();
	searchBox.setBorder(BorderFactory.createTitledBorder("Search clusterings"));
	searchBox.add(stringSearch);
	searchBox.add(goSearchButton);

	stringChoice = new JTextField(JAGCTClusteringPane.defaultStrings[algoSelect.getSelectedIndex()], 30);
	stringChoice.setEditable(false);
	stringBox = Box.createHorizontalBox(); 
	stringBox.setBorder(BorderFactory.createTitledBorder("Command line"));
	stringBox.add(stringChoice);

	generalBox = Box.createVerticalBox();
	generalBox.setBorder(BorderFactory.createTitledBorder("Compute clusterings"));
	generalBox.add(selectionBox);
	generalBox.add(stringBox);

	helpPane = new JTextArea();
	helpPane.setRows(1);
	helpPane.setText(JAGCTClusteringPane.helpStrings[algoSelect.getSelectedIndex()]);
	helpPane.setEditable(false);
	informationBox = Box.createHorizontalBox(); 
	informationBox.setBorder(BorderFactory.createTitledBorder("Help"));
	informationBox.add(new JScrollPane(helpPane));

	listClusteringPane = new JTextArea();
	listClusteringPane.setRows(1);
	listClusteringPane.setText(clusteringList());
	listClusteringPane.setEditable(false);
	clusterBox = Box.createHorizontalBox();
	clusterBox.setBorder(BorderFactory.createTitledBorder("Ordered list of clusterings"));
	clusterBox.add(new JScrollPane(listClusteringPane));

	superBox = Box.createVerticalBox();
	superBox.add(informationBox);
	superBox.add(clusterBox);

	upPanel = new JPanel();
	upPanel.setLayout(new BorderLayout());
	upPanel.add(generalBox, BorderLayout.NORTH);
	upPanel.add(searchBox, BorderLayout.SOUTH);

        setLayout(new BorderLayout());
	add(upPanel,BorderLayout.NORTH);
	add(superBox, BorderLayout.CENTER);
    }

    public String clusteringList(){
	String val = "";
	int i;
	if (myAGCT.allClusterings == null){
	    val += "None";
	    goSearchButton.setEnabled(false);
	}else{
	    goSearchButton.setEnabled(true);
	    for (i=0;i<myAGCT.allClusterings.size();i++)
		val += "Clustering #" + i + ": " + ((Clustering) myAGCT.allClusterings.elementAt(i)).toString() + "\n\n";
	}
	return val;
    }

    public void activate(){
	algoSelect.setEnabled(true);
	goButton.setEnabled(true);
	stringChoice.setEnabled(true);
	stringChoice.setEditable(true);
    }

    public void activateSearch(){
	goSearchButton.setEnabled(true);
 	stringSearch.setEnabled(true);
	stringSearch.setEditable(true);
    }

    public void changeAlgo(int i){
	stringChoice.setText(JAGCTClusteringPane.defaultStrings[i]);
	helpPane.setText("(Help on Clustering type)\n\n" + JAGCTClusteringPane.helpStrings[i]);
	helpPane.setCaretPosition(0);
    }

    public void displaySearch(String ref){
	String s, st;
	StringTokenizer t;
	boolean foundOneKeyword = false, foundCurrentKeyword, notHelp = false;
	int i, dumSize = 0, dumBest = 0, dumBestRatio = 0;
	String dumType = "";

	boolean activateSize = false, activateType = false, activateBest = false, activateBestRatio = false;

	if ( (ref != null) && (ref.length() > 0) ){
	    t = new StringTokenizer(ref, " ");
	    while (t.hasMoreTokens()){
		foundCurrentKeyword = false;
		s = t.nextToken();
		
		i=0;
		while( (i<JAGCTClusteringPane.searchKeywords.length) && (foundCurrentKeyword == false) )
		    if (s.equals(JAGCTClusteringPane.searchKeywords[i])){
			foundOneKeyword = true;
			foundCurrentKeyword = true;
		    }else
			i++;
		
		if (foundCurrentKeyword){
		    if (i==0)
			searchHelp();
		    else{
			notHelp = true;
			if (i==1){
			    if (t.hasMoreTokens()){
				dumSize = Integer.parseInt(t.nextToken());
				activateSize = true;
			    }
			}else if (i==2){
			    if (t.hasMoreTokens()){
				dumType = t.nextToken();
				activateType = true;
			    }
			}else if (i==3){
			    if (t.hasMoreTokens()){
				dumBest = Integer.parseInt(t.nextToken());
				activateBest = true;
			    }
			}else if (i==4){
			    if (t.hasMoreTokens()){
				dumBestRatio = Integer.parseInt(t.nextToken());
				activateBestRatio = true;
			    }
			}
		    }
		}
	    }

	    if ( (foundOneKeyword) && (notHelp) )
		search(activateSize, dumSize, activateType, dumType, activateBest, dumBest, activateBestRatio, dumBestRatio);

	    if (!foundOneKeyword)
		stringSearch.setText(JAGCTClusteringPane.searchKeywords[0]);
	}
    }

    public void search(boolean bSize, int dSize, boolean bType, String sType, boolean bBest, int dBest, boolean bBestRatio, int dBestRatio){
	int i, numb = 0;
	if ( (myAGCT.allClusterings != null) && (myAGCT.allClusterings.size() > 0) && (bSize || bType || bBest || bBestRatio) ){
	    helpPane.setText(headerSearch(bSize, dSize, bType, sType, bBest, dBest, bBestRatio, dBestRatio));
	    for (i=0;i<myAGCT.allClusterings.size();i++)
		if (globalTest(i, bSize, dSize, bType, sType, bBest, dBest, bBestRatio, dBestRatio)){
		    if (numb == 0)
			helpPane.append("Indexes: ");
		    helpPane.append(i + " ");
		    numb ++;
		}
	}
	if (numb == 0)
	    helpPane.setText("None found");
    }

    public String headerSearch(boolean bSize, int dSize, boolean bType, String sType, boolean bBest, int dBest, boolean bBestRatio, int dBestRatio){
	String val = "Searching clusterings matching the following condition: ";
	int n = 0;
	if (bSize){
	    val += "(number of clusters = " + dSize + ")";
	    n++;
	}
	if (bType){
	    if (n>0)
		val += " AND ";
	    val += "(type of clusterings = " + sType + ")";
	    n++;
	}
	if (bBest){
	    if (n>0)
		val += " AND ";
	    val += "(belongs to the " + dBest + " % of best clustering wrt objective function)";
	    n++;
	}
	if (bBestRatio){
	    if (n>0)
		val += " AND ";
	    val += "(belongs to the " + dBestRatio + " % of best clustering wrt objective function / number of clusters)";
	    n++;
	}
	val += "\n";
	return val;
    }
    
    public boolean globalTest(int nCluster, boolean bSize, int dSize, boolean bType, String sType, boolean bBest, int dBest, boolean bBestRatio, int dBestRatio){
	return ( ( (!bSize) || (testSize(nCluster, dSize)) )
		 && ( (!bType) || (testType(nCluster, sType)) )
		 && ( (!bBest) || (testBest(nCluster, dBest)) )
		 && ( (!bBestRatio) || (testBestRatio(nCluster, dBestRatio)) ) );
    }

    public boolean testBest(int nCluster, int dumProp){
	double baseline = ((Clustering) myAGCT.allClusterings.elementAt(nCluster)).myClusteringAlgorithm.finalClusteringObjectiveFunction;
	double nBest = 0.0, nTot = (double) myAGCT.allClusterings.size();
	int i, prop;
	for (i=0;i<myAGCT.allClusterings.size();i++)
	    if ( ((Clustering) myAGCT.allClusterings.elementAt(i)).myClusteringAlgorithm.finalClusteringObjectiveFunction <= baseline)
		nBest += 1.0;
	prop = (int) (100.0 * nBest / nTot);

	if (prop < dumProp)
	    return true;
	else
	    return false;
    }

    public boolean testBestRatio(int nCluster, int dumProp){
	double baseline = ((Clustering) myAGCT.allClusterings.elementAt(nCluster)).myClusteringAlgorithm.finalClusteringObjectiveFunction
	    / (double) ((Clustering) myAGCT.allClusterings.elementAt(nCluster)).nclusters;
	double nBest = 0.0, nTot = (double) myAGCT.allClusterings.size();
	int i, prop;
	for (i=0;i<myAGCT.allClusterings.size();i++)
	    if ( (((Clustering) myAGCT.allClusterings.elementAt(i)).myClusteringAlgorithm.finalClusteringObjectiveFunction /
		  (double) ((Clustering) myAGCT.allClusterings.elementAt(i)).nclusters) <= baseline)
		nBest += 1.0;
	prop = (int) (100.0 * nBest / nTot);

	if (prop < dumProp)
	    return true;
	else
	    return false;
    }

    public boolean testSize(int nCluster, int dumSize){
	if ( ((Clustering) myAGCT.allClusterings.elementAt(nCluster)).nclusters == dumSize )
	    return true;
	else
	    return false;
    }

    public boolean testType(int nCluster, String dumType){
	if ( ((Clustering) myAGCT.allClusterings.elementAt(nCluster)).myClusteringAlgorithm.myReferenceName.equals(dumType) )
	    return true;
	else
	    return false;
    }

    public void searchHelp(){
	int i;
	helpPane.setText("(Help on Search Clustering)\n\n");
	for (i=0;i<JAGCTClusteringPane.searchKeywords.length;i++)
	    helpPane.append(JAGCTClusteringPane.searchKeywords[i] + ":\n " + JAGCTClusteringPane.helpSearchStrings[i] + "\n\n");
	helpPane.append(JAGCTClusteringPane.exampleSearchString);
	helpPane.setCaretPosition(0);
    }
}
