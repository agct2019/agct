import java.lang.*;
import java.util.*;
import java.io.*;

class Scenario {
    /*******************************************************************************************
     * Class that records and plays scenarios
     *****/

    public static int METHOD_SONY = 0, METHOD_UAG = 1;

    public static int MATRIX_STORING_METHOD = METHOD_UAG; 
    // METHOD_UAG = all coordinates are stored
    // METHOD_SONY = couples (coordinate, value) are stored, coordinate giving the y index where to put value (allows reduced storage for W in particular)

    public static String allKeywords [] = {"Load", 
					   "JAGCTSelectionPane_Method_F",                         
					   "JAGCTSelectionPane_Validate",                        
					   "JAGCTCheckBox_StateChanged",                         
					   "AGCT_Modification_Method_F",
					   "AGCT_Modification_Method_W",
					   "AGCT_Modification_Method_S",
					   "AGCT_Modification_Method_N",
					   "AGCT_Modification_Number_Of_Neighbors",
					   "AGCT_Modification_Number_Of_Manifold_Components",
					   "AGCT_Modification_Number_Of_Triangulation_Dimensions", // 10
					   "AGCT_Modification_T_Heat_Kernel",
					   "AGCT_Modification_T_Heat_Kernel_Similarity_Neighbors",
					   "AGCT_Modification_Sparsify_Statistic",
					   "JDSScaling_Modification_Bregman_Divergence",
					   "Statistics_Modification_Limit_P",
					   "AGCT_Modification_Use_Debug",
					   "AGCT_Modification_Warning",
					   "AGCT_Modification_Perspective",
					   "AGCT_Modification_Sort_Depth",
					   "AGCT_Modification_Use_Shadow",                         // 20
					   "AGCT_Manifold_Processed",
					   "AGCT_PCA_Processed",
					   "AGCT_Manifold_Triangulated",
					   "JAGCTClusteringPane_Run_Clustering",
					   "AGCT_Modification_Number_Of_Wavelet_Stamps",          
					   "JAGCTSelectionPane_Max_Number_Of_Features",           
					   "JAGCTSelectionPane_Feature_Selection_Method",         
					   "JAGCTSelectionPane_Max_Number_Of_Prototypes",         
					   "JAGCTSelectionPane_Prototype_Selection_Method",
					   "JAGCTSelectionPane_Keep_Prototypes",                    // 30
					   "JAGCTSelectionPane_Prototypes_Begin",
					   "JAGCTSelectionPane_Prototypes_End",
					   "JAGCTManifoldPane_Keep_Manifold",                  
					   "JAGCTManifoldPane_Manifold_Begin",
					   "JAGCTManifoldPane_Manifold_End",
					   "AGCT_Save_Constant_Wavelet"};

    //Manifold = ORDERED list of gene names, W (last step), M (ready for components)

    public static boolean runningForPrototypes, bClosest_Center_Ordered_List, bClosest_Center_To_Cluster_Number, bCluster_Number_To_Closest_Center, bClosest_Center_To_Cluster_Points, bClosest_Center_To_Normalized_Distortions;

    public static boolean runningForManifold, bOrdered_List_Names, bMatrix_W, bManifold_Eigenvalues, bMatrix_M;

    public static String myToken = "|";

    public static String commentString = "//";

    //Load | datafile = loads a datafile
    //JAGCTSelectionPane_Method_F | Integer = codes the Method_F variable
    //JAGCTSelectionPane_Validate | "" = validates on JAGCTSelectionPane
    //JAGCTCheckBox_StateChanged | Integer Type = gives a checkbox whose checked status has changed

    public static Vector allStrings = new Vector();
    //Each string = an action;

    public static boolean scenarioTick = false, isRunning = false;

    public static void block(int i){
	if ( ( (i>=1) && (i<=3) ) || ( (i>= 21) && (i<=29) ) ){
	    while(!ControlProcess.hasTrue("displaySelectionPane"));
	    //while(!ControlProcess.hasTrue("dataLoaded"));
	}
	if (i==21)
	    while(!ControlProcess.hasTrue("featuresSelected"));
	if (i==21)
	    while(!ControlProcess.hasTrue("prototypeSelected"));
	if (i==21)
	    while(!ControlProcess.hasTrue("featuresAndPrototypesSelected"));

	if ( (i==34) || (i==35) )
	    while(!ControlProcess.hasTrue("featuresSelected"));
	if ( (i==34) || (i==35) )
	    while(!ControlProcess.hasTrue("prototypeSelected"));
	if ( (i==34) || (i==35) )
	    while(!ControlProcess.hasTrue("featuresAndPrototypesSelected"));
	if ( (i==34) || (i==35) )
	    while(!ControlProcess.hasTrue("annotationsLoaded"));

	//if (i==31)
	//  while(!ControlProcess.hasTrue("selectedGenesAllAnnotations"));
	//  while(!ControlProcess.hasTrue("annotationsLoaded"));

	if ( (i==22) || (i==23) )
	    while(!ControlProcess.hasTrue("manifoldProcessed"));
	if (i==24)
	    while(!ControlProcess.hasTrue("clusteringAvailable"));
    }

    public static void flushAll(){
	runningForPrototypes = runningForManifold = false;
	allStrings = null;
	scenarioTick = false;
	isRunning = false;
    }

    public static void tickOn(){
	scenarioTick = true;
    }

    public static void tickOff(){
	scenarioTick = false;
    }

    public static void tickSwitch(){
	scenarioTick = !scenarioTick;
    }

    public static boolean isTicking(){
	return scenarioTick;
    }

    public static boolean isRunning(){
	return isRunning;
    }

    public static void newScenario(){
	allStrings = new Vector();
    }

    public static void addThis(String t){
	if (allStrings == null)
	    newScenario();
	allStrings.addElement(t);
    }

    public static void affiche(){
	if (allStrings != null){
	    int i;
	    for (i=0;i<allStrings.size();i++)
		System.out.println(allStrings.elementAt(i));
	}
    }

    public static boolean containsKwd(String t){
	if (t == null)
	    return false;

	int i;
	for (i=0;i<allKeywords.length;i++)
	    if (allKeywords[i].equals(t))
		return true;
	return false;
    }

    public static int kwdIndex(String t){
	if (t == null)
	    return -1;

	int i;
	for (i=0;i<allKeywords.length;i++)
	    if (allKeywords[i].equals(t))
		return i;
	return -1;
    }

    public static void add(String kwd){
	add(kwd, "");
    }

    public static void add(String kwd, int i){
	add(kwd, i + "");
    }

    public static void add(String kwd, double i){
	add(kwd, i + "");
    }

    public static void add(String kwd, boolean i){
	add(kwd, boolean2int(i) + "");
    }

    public static void add(String kwd, String arg){
	if (scenarioTick){
	    if (!containsKwd(kwd))
		Matrix.perror("Scenario.class :: Keyword " + kwd + " not registered");
	    addThis(kwd + myToken + arg);
	}
    }

    public static boolean isComment(String s){
	if (s.length() < Scenario.commentString.length())
	    return false;
	String t = s.substring(0,2);
	if (t.equals(Scenario.commentString))
	    return true;
	return false;
    }


    public static boolean isKeepPrototypes(String s){
	int ref = allKeywords[30].length();

	if (s.length() < ref)
	    return false;
	String t = s.substring(0,ref);
	if (t.equals(allKeywords[30]))
	    return true;
	return false;
    }

    public static boolean isKeepManifold(String s){
	int ref = allKeywords[33].length();

	if (s.length() < ref)
	    return false;
	String t = s.substring(0,ref);
	if (t.equals(allKeywords[33]))
	    return true;
	return false;
    }


    public static void loadExecute(File rf, AGCT ap){
	FileReader e;
	BufferedReader br;
	String s;

	try{
	    e = new FileReader(rf);
	}catch(FileNotFoundException ex){
	    ap.myInformationFrame.setText("The scenario you try to load does not exist");
	    return;
	}

	newScenario();
	br = new BufferedReader(e);
	try{
	    while ( (s=br.readLine()) != null){
		if ( (!isComment(s)) && (!s.equals("")) )
		    addThis(s);
	    }
	    e.close();
	}catch(IOException ex){
	    ap.myInformationFrame.setText("IOError when loading scenario");
	    ap.rawDomainExists = false;
	    return;
	}
	e = null;

	if ( (allStrings != null) && (allStrings.size() > 0) )
	    executeScenario(ap);
    }

    public static void executeScenario(final AGCT ap){
	Thread t = new Thread(){
		public void run(){
		    int i;
		    isRunning = true;
		    System.out.print("Executing scenario (" + allStrings.size() + " lines)... please wait... ");

		    ap.myTabbedPane.setAllEnabled(false);

		    AGCTCounter cc = new AGCTCounter(ap.myInformationFrame, "Scenario...", allStrings.size());
		    for (i=0;i<allStrings.size();i++){
			executeLine((String) allStrings.elementAt(i), ap);
			cc.setText("Scenario...");
			cc.increment();
		    }
		    cc.end();

		    flushAll(); //isRunning = false; (saves memory)

		    ap.myTabbedPane.setAllEnabled(true);

		    ap.scenarioLoadButton.setIcon(ap.scenarioLoadIdleIcon);
		}
	    };
	t.start();
    }

    public static int index(String s){
	int i=0;
	while (i < allKeywords.length){
	    if (s.equals(allKeywords[i]))
		return i;
	    else
		i++;
	}
	return -1;
    }

    public static int indexFirstToken(String ex){
	StringTokenizer t = new StringTokenizer(ex, Scenario.myToken);
	return index(t.nextToken());
    }

    public static boolean int2boolean(int dum){
	if (dum == 0)
	    return false;
	else if (dum == 1)
	    return true;
	else
	    Matrix.perror("Scenario :: not a boolean");
	return false;
    }

    public static int boolean2int(boolean dum){
	if (dum)
	    return 1;
	return 0;
    }

    public static void displayInformationFrame(AGCT ap){
	String s = "";

	if ( (allStrings == null) || (allStrings.size() == 0) )
	    s = "No scenario to display";
	else{
	    s = "Currently recorded scenario:\n";
	    int i;
	    for (i=0;i<allStrings.size();i++){
		s += ((String) allStrings.elementAt(i)) + "\n";
	    }
	}
	ap.myInformationFrame.setText(s);
    }

    public static void executeLine(String ex, AGCT ap){

	if (runningForPrototypes)
	    executeLinePrototypes(ex, ap);
	else if ( (runningForManifold) && ( (indexFirstToken(ex) == 34) || (indexFirstToken(ex) == -1) ) )
	    executeLineManifold(ex, ap);
	else
	    executeLineNormal(ex, ap);
    }

    

    public static boolean flagChange_Prototypes(String ex){
	int i;
	i = Prototype.indexInTokens_Prototypes(ex);
	if (i == -1)
	    return false;
	else{
	    if (i == 0)
		bClosest_Center_To_Cluster_Number = true;
	    else if (i==1)
		bClosest_Center_To_Cluster_Number = false;
	    else if (i==2)
		bClosest_Center_To_Cluster_Points = true;
	    else if (i==3)
		bClosest_Center_To_Cluster_Points = false;
	    else if (i==4){
		bClosest_Center_To_Normalized_Distortions = true;

		/*Enumeration extensions = Prototype.Closest_Center_To_Cluster_Number.keys();
		Integer R, vp;
		while (extensions.hasMoreElements()) {
		    R = (Integer) extensions.nextElement();
		    vp = (Integer) Prototype.Closest_Center_To_Cluster_Number.get(R);
		    System.out.print("(" + R + ":: " + vp + ")");
		    }*/
	    }else if (i==5)
		bClosest_Center_To_Normalized_Distortions = false;
	    else if (i==6)
		bCluster_Number_To_Closest_Center = true;
	    else if (i==7)
		bCluster_Number_To_Closest_Center = false;
	    else if (i==8){
		runningForPrototypes = false;
		Prototype.Loading_From_Scenario_End = true;
		Prototype.Prototypes_Selected = true;
		ControlProcess.put("prototypeSelected",true);

		/*Enumeration extensions = Prototype.Closest_Center_To_Cluster_Number.keys();
		Integer R, vp;
		while (extensions.hasMoreElements()) {
		    R = (Integer) extensions.nextElement();
		    vp = (Integer) Prototype.Closest_Center_To_Cluster_Number.get(R);
		    System.out.print("(" + R + ": " + vp + ")");
		}
		System.exit(0);*/
	    }else if (i==9)
		bClosest_Center_Ordered_List = true;
	    else if (i==10)
		bClosest_Center_Ordered_List = false;
	    

	    //System.out.println("Bluk proto = " + i);

	    return true;
	}
    }

    public static boolean flagChange_Manifold(String ex, AGCT ap){
	int i;
	i = AGCT.indexInTokens_Manifold(ex);
	if (i == -1)
	    return false;
	else{
	    if (i == 0)
		bOrdered_List_Names = true;
	    else if (i==1)
		bOrdered_List_Names = false;
	    else if (i==2)
		bMatrix_W = true;
	    else if (i==3)
		bMatrix_W = false;
	    else if (i==4)
		bManifold_Eigenvalues = true;
	    else if (i==5)
		bManifold_Eigenvalues = false;
	    else if (i==6)
		bMatrix_M = true;
	    else if (i==7)
		bMatrix_M = false;
	    else if (i==8){
		runningForManifold = false;
		AGCT.Loading_From_Scenario_Manifold_End = true;
		//ControlProcess.put("manifoldProcessed",true);
		ap.justAfterManifold();
	    }else if (i==9){
		System.out.println("W loaded. On to rest of manifold (computing)");

		runningForManifold = false;
		AGCT.Loading_From_Scenario_Manifold_End = true;
		ap.manifoldEigensystem(false);
	    }

	    return true;
	}
    }


    public static void executeLinePrototypes(String ex, AGCT ap){
	if (!flagChange_Prototypes(ex)){
	    Prototype.executeScenario(ex, ap);
	}
    }

    public static void executeLineManifold(String ex, AGCT ap){
	if (!flagChange_Manifold(ex, ap)){
	    AGCT.executeManifold(ex, ap);
	}
    }

    public static void executeLineNormal(String ex, AGCT ap){
	StringTokenizer t;
	String s, nam;
	int dum, mI, mJ, sel, index;
	double ddum;
	boolean bdum, isSel;
	t = new StringTokenizer(ex, Scenario.myToken);
	s = t.nextToken();
	index = index(s);
	//System.out.println("Index = " + index + "Bulk line = " + ex);

	//System.out.println(ControlProcess.hasTrue("featuresSelected") + " " + ControlProcess.hasTrue("prototypeSelected") + " " + ControlProcess.hasTrue("featuresAndPrototypesSelected"));

	//blocks until conditions are met to run the process
	block(index);
	//System.out.println("Bulk line = " + ex);

	if (index ==0){
	    ap.goLoadDomain(t.nextToken());
	}else if (index == 1){
	    dum = Integer.parseInt(t.nextToken());
	    ap.myTabbedPane.mySelectionPane.setMethod_F(dum);
	    ap.myTabbedPane.mySelectionPane.fm.setSelectedIndex(dum);
	}else if (index == 2){
	    ap.myTabbedPane.mySelectionPane.goValidate();
	}else if (index == 3){
	    nam = t.nextToken();
	    mI = Integer.parseInt(t.nextToken());
	    mJ = Integer.parseInt(t.nextToken());
	    sel = Integer.parseInt(t.nextToken());

	    isSel = false;
	    if (sel == 1)
		isSel = true;
	    else if (sel == 0)
		isSel = false;
	    else
		Matrix.perror("Scenario :: non boolean selection");

	    if (nam.equals(JAGCTCheckBox.refGene))
		((JAGCTCheckBox) ap.myTabbedPane.mySelectionPane.checkBoxGene.elementAt(mI)).setSelected(isSel);
	    else if (nam.equals(JAGCTCheckBox.refLigand))
		((JAGCTCheckBox) ap.myTabbedPane.mySelectionPane.checkBoxLigand.elementAt(mI)).setSelected(isSel);
	    else if (nam.equals(JAGCTCheckBox.refGroup))
		((JAGCTCheckBox) ((Vector) ap.myTabbedPane.mySelectionPane.checkBoxGroup.elementAt(mI)).elementAt(mJ)).setSelected(isSel);
	}else if (index == 4){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationMethod_F(dum);
	}else if (index == 5){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationMethod_W(dum);
	}else if (index == 6){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationMethod_S(dum);
	}else if (index == 7){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationMethod_N(dum);
	}else if (index == 8){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationNumber_Of_Neighbors(dum);
	}else if (index == 9){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationNumber_Of_Manifold_Components(dum);
	}else if (index == 10){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationNumber_Of_Triangulation_Dimensions(dum);
	}else if (index == 11){
	    ddum = Double.parseDouble(t.nextToken());
	    ap.requestModificationT_Heat_Kernel(ddum);
	}else if (index == 12){
	    ddum = Double.parseDouble(t.nextToken());
	    ap.requestModificationT_Heat_Kernel_Similarity_Neighbors(ddum);
	}else if (index == 13){
	    ddum = Double.parseDouble(t.nextToken());
	    ap.requestModificationSparsify_Statistic(ddum);
	}else if (index == 14){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationBregDiv(dum);
	}else if (index == 15){
	    ddum = Double.parseDouble(t.nextToken());
	    ap.requestModificationLimit_P_Chi2(ddum);
	}else if (index == 16){
	    dum = Integer.parseInt(t.nextToken());
	    bdum = int2boolean(dum);
	    ap.requestModificationUseDebug(bdum);
	}else if (index == 17){
	    dum = Integer.parseInt(t.nextToken());
	    bdum = int2boolean(dum);
	    ap.requestModificationUseWarning(bdum);
	}else if (index == 18){
	    dum = Integer.parseInt(t.nextToken());
	    bdum = int2boolean(dum);
	    ap.requestModificationPerspective(bdum);
	}else if (index == 19){
	    dum = Integer.parseInt(t.nextToken());
	    bdum = int2boolean(dum);
	    ap.requestModificationSort_Depth(bdum);
	}else if (index == 20){
	    dum = Integer.parseInt(t.nextToken());
	    bdum = int2boolean(dum);
	    ap.requestModificationUse_Shadow(bdum);
	}else if (index == 21){
	    ap.myTabbedPane.myManifoldPane.goManifold();
	}else if (index == 22){
	    ap.myTabbedPane.myPCAPane.goPCA();
	}else if (index == 23){
	    ap.myTabbedPane.myManifoldPane.goTriangulation();
	}else if (index == 24){
	    ap.nonThreadClustering(t.nextToken());
	}else if (index == 25){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationNumber_Of_Wavelet_Stamps(dum);
	    ap.myTabbedPane.mySelectionPane.moreNumber_Of_Wavelet_Stamps();
	}else if (index == 26){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationMax_Number_Of_Features(dum);
	    ap.myTabbedPane.mySelectionPane.moreMax_Number_Of_Features();
	}else if (index == 27){
	    dum = Integer.parseInt(t.nextToken());
	    ap.myTabbedPane.mySelectionPane.playMethod_FS(dum);
	}else if (index == 28){
	    dum = Integer.parseInt(t.nextToken());
	    ap.requestModificationMax_Number_Of_Prototypes(dum);
	    ap.myTabbedPane.mySelectionPane.moreMax_Number_Of_Prototypes();
	}else if (index == 29){
	    dum = Integer.parseInt(t.nextToken());
	    ap.myTabbedPane.mySelectionPane.playMethod_PS(dum);
	}else if (index == 31){
	    runningForPrototypes = true;
	    Prototype.flushAll();
	    Prototype.Loading_From_Scenario_Begin = true;
	    bClosest_Center_Ordered_List = bClosest_Center_To_Cluster_Number = bCluster_Number_To_Closest_Center = bClosest_Center_To_Cluster_Points = bClosest_Center_To_Normalized_Distortions = false;
	}else if (index == 32){
	    runningForPrototypes = false;
	}else if (index == 34){
	    runningForManifold = true;
	    AGCT.Loading_From_Scenario_Manifold_Begin = true;
	    bOrdered_List_Names = bMatrix_W = bManifold_Eigenvalues = bMatrix_M = false;
	    ap.justBeforeManifold();
	}else if (index == 35){
	    ap.justAfterManifold();
	    runningForManifold = false;
	}else if (index == 36){
	    dum = Integer.parseInt(t.nextToken());
	    bdum = int2boolean(dum);
	    ap.requestModificationSaveConstantWavelet(bdum);
	}
    }
}
