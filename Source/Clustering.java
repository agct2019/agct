import java.util.*;
import java.awt.Color;

 class Clustering implements Debuggable{

    static String Center_Name = "C#_";
    static String Keyword_Loop = "-L";
    static String Keyword_Repeat = "-R";
    static String Keyword_Generator = "-Generator";
    static String Keyword_Parametre = "-G";

    AGCT myAGCT;

    Matrix VV;
    // (nClusters x nFeatures) matrix which gives the principal components of clusters

    public static final String clusteringAlgorithms[] = {"CP", "KM", "EM", "HC", "AP", "NM"};
    static final String referenceStringCenters[] = {"Manifold", "Manifold3D", "PCA", "FinalCoordinates_1D"};
    static final String referenceSeeding[] = {"Forgy", "Arthur", "Bregman"};
    static final String referenceGenerator[] = {"L2", "KL", "IS", "AM", "PN", "NM"};
    static final String referenceNMFPositiveGenerator[] = {"KT", "MS"};

    static final String HC_DISTANCE[] = {"Ward", "Single_Linkage"};

    double [] redPropClusters, greenPropClusters, bluePropClusters;
    int [] redCluster, greenCluster, blueCluster;
    Color [] colorCluster;

    int nclusters, indexReferenceStringCenters, clusteringAlgorithmId, myNumber, valueP, seeding, generator, nMaxIterations, nIterative, nVariables, nMethod, tUpdate;
    double parametre, valueBregman;
    int indexDistance;
    //myNumber = index of the Clustering in AGCT
    String bornString;

    public static int RGB_Max = 250;

    AGCTClustering_Algorithm myClusteringAlgorithm;

    Clustering(AGCT mya, int numb){
	myNumber = numb;
	myAGCT = mya;

	VV = null;
	redPropClusters = null;
	greenPropClusters = null;
	bluePropClusters = null;
	bornString = null;
	seeding = 0;
	generator = 0;
	parametre = -1.0;
	nMaxIterations = 0;
    }

    public void saveSpace(){
	myClusteringAlgorithm.saveSpace();
    }

    public void plot(JAGCTVisualizationPane v){
	if (myClusteringAlgorithm != null)
	    if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM))
		((AGCTClustering_KM) myClusteringAlgorithm).plot(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM))
		((AGCTClustering_EM) myClusteringAlgorithm).plot(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP))
		((AGCTClustering_CP) myClusteringAlgorithm).plot(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC))
		((AGCTClustering_HC) myClusteringAlgorithm).plot(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP))
		((AGCTClustering_AP) myClusteringAlgorithm).plot(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM))
		((AGCTClustering_NM) myClusteringAlgorithm).plot(v);
    }

    public void plotWithAnnotations(JAGCTVisualizationPane v){
	if (myClusteringAlgorithm != null)
	    if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM))
		((AGCTClustering_KM) myClusteringAlgorithm).plotWithAnnotations(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM))
		((AGCTClustering_EM) myClusteringAlgorithm).plotWithAnnotations(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP))
		((AGCTClustering_CP) myClusteringAlgorithm).plotWithAnnotations(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC))
		((AGCTClustering_HC) myClusteringAlgorithm).plotWithAnnotations(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP))
		((AGCTClustering_AP) myClusteringAlgorithm).plotWithAnnotations(v);
	    else if (myClusteringAlgorithm.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM))
		((AGCTClustering_NM) myClusteringAlgorithm).plotWithAnnotations(v);
    }

    public String toString(){
	return myClusteringAlgorithm.toString();
    }

    public static boolean isAlgo(String ref, String thisAlgo){
	String s;
	StringTokenizer t = new StringTokenizer(ref, " ");
	while (t.hasMoreTokens()){
	    s = t.nextToken();
	    if (s.equals(thisAlgo))
		return true;
	}
	return false;
    }

    public static int getIndex(String s){
	int val = -1;
	int i = 0;
	do{
	    if (s.equals(Clustering.clusteringAlgorithms[i]))
		val = i;
	    i++;
	}while ( (val == -1) && (i < Clustering.clusteringAlgorithms.length) );
	return val;
    }

    public static boolean containsLoop(String ref){
	StringTokenizer t = new StringTokenizer(ref, " ");
	String s;
	while(t.hasMoreTokens()){
	    s = t.nextToken();
	    if (s.equals(Clustering.Keyword_Loop))
		return true;
	}
	return false;
    }

    public static boolean containsRepeat(String ref){
	StringTokenizer t = new StringTokenizer(ref, " ");
	String s;
	while(t.hasMoreTokens()){
	    s = t.nextToken();
	    if (s.equals(Clustering.Keyword_Repeat))
		return true;
	}
	return false;
    }

    public static int [] getBoundsNclusters(String ref){
	int [] bds = new int[2];
	StringTokenizer t = new StringTokenizer(ref, " ");
	String s;

	bds[0] = bds[1] = -1;
	while(t.hasMoreTokens()){
	    s = t.nextToken();
	    if (s.equals(Clustering.Keyword_Loop)){
		s = t.nextToken();
		bds[0] = Integer.parseInt(s);
		s = t.nextToken();
		bds[1] = Integer.parseInt(s);
	    }	
	}
	return bds;
    }

    public static String getReplacementString(int ncl, String ref){
	StringTokenizer t = new StringTokenizer(ref, " ");
	String s, sret = "";

	while(t.hasMoreTokens()){
	    s = t.nextToken();
	    if (s.equals(Clustering.Keyword_Loop)){
		s = t.nextToken();
		s = t.nextToken();
		
		if (isAlgo(ref, Clustering.clusteringAlgorithms[4]))
		    sret += "-P " + ncl + " ";
		else
		    sret += "-K " + ncl + " ";
	    }else
		sret += s + " ";
	}
	return sret;
    }

    public static String flushRepeatString(String ref){
	StringTokenizer t = new StringTokenizer(ref, " ");
	String s, sret = "";

	while(t.hasMoreTokens()){
	    s = t.nextToken();
	    if (s.equals(Clustering.Keyword_Repeat)){
		s = t.nextToken();
	    }else
		sret += s + " ";
	}
	return sret;
    }

    public static int getRepeatTimes(String ref){
	StringTokenizer t = new StringTokenizer(ref, " ");
	String s;
	int val = -1;
	boolean found = false;

	while(t.hasMoreTokens()){
	    s = t.nextToken();
	    if (s.equals(Clustering.Keyword_Repeat)){
		s = t.nextToken();
		val = Integer.parseInt(s);
		found = true;
	    }
	}

	if ( (!found) || (val == -1) )
	    Matrix.perror("Negative or unspecified repeat times");

	return val;
    }

    public void getOptions(String ref){
	bornString = ref;
	//Tokens recenses: -K, -I, -Iterative, -Point, -P, Clustering.Keyword_Generator -G, -B, -Distance, -M, -U

	String s, st;
	boolean foundAlgo = false;
	int id = -1, i;
	indexReferenceStringCenters = -1;
	indexDistance = -1;
	nMaxIterations = -1;
	nIterative = 0;
	nMethod = tUpdate = -1;

	StringTokenizer t = new StringTokenizer(ref, " ");
	while (t.hasMoreTokens()){
	    s = t.nextToken();

	    if (!foundAlgo){
		id = Clustering.getIndex(s);
		if (id != -1)
		    foundAlgo = true;
	    }

	    if (s.equals(new String("-K"))){
		nclusters = Integer.parseInt(t.nextToken());
	    }

	    if (s.equals(new String("-I"))){
		nMaxIterations = Integer.parseInt(t.nextToken());
	    }

	    if (s.equals(new String("-Iterative"))){
		nIterative = Integer.parseInt(t.nextToken());
	    }

	    if (s.equals(new String("-Point"))){
		st = t.nextToken();
		for (i=0;i<Clustering.referenceStringCenters.length;i++)
		    if (st.equals(Clustering.referenceStringCenters[i]))
			indexReferenceStringCenters = i;
	    }

	    if (s.equals(new String("-P"))){
		valueP = Integer.parseInt(t.nextToken());
	    }

	    if (s.equals(new String("-Seed"))){
		st = t.nextToken();
		if (st.equals(Clustering.referenceSeeding[0]))
		    seeding = 0;
		else if (st.equals(Clustering.referenceSeeding[1]))
		    seeding = 1;
		else if (st.equals(Clustering.referenceSeeding[2]))
		    seeding = 2;
		else
		    seeding = 0;
	    }

	    if (s.equals(new String("-M"))){
		st = t.nextToken();
		if (st.equals(Clustering.referenceNMFPositiveGenerator[0]))
		    nMethod = 0;
		else if (st.equals(Clustering.referenceNMFPositiveGenerator[1]))
		    nMethod = 1;
		else
		    nMethod = 0;
	    }

	    if (s.equals(new String("-U"))){
		st = t.nextToken();
		if (st.equals(Clustering.referenceGenerator[0]))
		    tUpdate = 0;
		else if (st.equals(Clustering.referenceGenerator[1]))
		    tUpdate = 1;
		else
		    tUpdate = 0;
	    }

	    if (s.equals(Clustering.Keyword_Generator)){
		st = t.nextToken();
		if (st.equals(Clustering.referenceGenerator[0]))
		    generator = 0;
		else if (st.equals(Clustering.referenceGenerator[1]))
		    generator = 1;
		else if (st.equals(Clustering.referenceGenerator[2]))
		    generator = 2;
		else if (st.equals(Clustering.referenceGenerator[3]))
		    generator = 3;
		else if (st.equals(Clustering.referenceGenerator[4]))
		    generator = 4;
		else if (st.equals(Clustering.referenceGenerator[5]))
		    generator = 5;
		else
		    generator = 0;
	    }

	    if (s.equals(new String("-G"))){
	        parametre = Double.parseDouble(t.nextToken());
	    }

	    if (s.equals(new String("-B"))){
	        valueBregman = Double.parseDouble(t.nextToken());
	    }

	    if (s.equals(new String("-Distance"))){
		st = t.nextToken();
		for (i=0;i<Clustering.HC_DISTANCE.length;i++)
		    if (st.equals(Clustering.HC_DISTANCE[i]))
			indexDistance = i;
	    }
	}

	if (id == -1)
	    Matrix.perror("No such clustering algorithm");

	initOptions(id);
    }

    public void toClusterColors(){
	int i;

	redPropClusters = new double [nclusters];
	greenPropClusters = new double [nclusters];
	bluePropClusters = new double [nclusters];
	
	for (i=0;i<nclusters;i++){
	    redPropClusters[i] = (double) i / ( ( (double) nclusters ) - 1.0 );
	    greenPropClusters[i] = ( (redPropClusters[i]<=0.5) ?  (2.0*redPropClusters[i]) : (2.0*(1.0-redPropClusters[i])) );
	    bluePropClusters[i] = 1.0 - redPropClusters[i];
	}

	redCluster = new int [nclusters];
	greenCluster = new int [nclusters];
	blueCluster = new int [nclusters];
	colorCluster = new Color [nclusters];

	for (i=0;i<nclusters;i++){
	    redCluster[i] = (int) ( (double) Clustering.RGB_Max * (1.0 - redPropClusters[i]) );
	    redCluster[i] = Clustering.RGB_Max - redCluster[i];

	    greenCluster[i] = (int) ( (double) Clustering.RGB_Max * (1.0 - greenPropClusters[i]) );
	    greenCluster[i] = Clustering.RGB_Max - greenCluster[i];

	    blueCluster[i] = (int) ( (double) Clustering.RGB_Max * (1.0 - bluePropClusters[i]) );
	    blueCluster[i] = Clustering.RGB_Max - blueCluster[i];

	    colorCluster[i] = new Color(redCluster[i], greenCluster[i], blueCluster[i]);
	}
    }

    public void initOptions(int id){

	clusteringAlgorithmId = id;
	int i, j;
	boolean kd;

	if (nclusters == -1)
	    nclusters = AGCT.Number_Of_Clusters;

	if (id != 4)
	    toClusterColors();

	if (id == 0){
	    myClusteringAlgorithm = new AGCTClustering_CP(this, myAGCT, nclusters, false);
	    myClusteringAlgorithm.initSoftMemberships();
	}else if (id == 1){
	    if (seeding > 1){
		seeding = 1;
		kd = true;
	    }else
		kd = false;

	    myClusteringAlgorithm = new AGCTClustering_KM(this, myAGCT, nclusters, seeding, kd, generator, parametre, false, nMaxIterations, nIterative);
	    myClusteringAlgorithm.initSoftMemberships();
	    myClusteringAlgorithm.initHardMemberships();
	}else if (id == 2){
	    if (seeding > 1)
		seeding = 1;//No General Bregman allocation

	    myClusteringAlgorithm = new AGCTClustering_EM(this, myAGCT, nclusters, seeding, false, nMaxIterations);
	    myClusteringAlgorithm.initSoftMemberships();
	    myClusteringAlgorithm.initHardMemberships();
	}else if (id == 3){
	    myClusteringAlgorithm = new AGCTClustering_HC(this, myAGCT, nclusters, indexDistance, false);
	    myClusteringAlgorithm.initSoftMemberships();
	}else if (id == 4){
	    myClusteringAlgorithm = new AGCTClustering_AP(this, myAGCT, valueP, valueBregman, false, nMaxIterations);
	}else if (id == 5){
	    myClusteringAlgorithm = new AGCTClustering_NM(this, myAGCT, nclusters, false, nMaxIterations, nMethod, tUpdate);
	    myClusteringAlgorithm.initSoftMemberships();
	    myClusteringAlgorithm.initHardMemberships();
	}
    }

    public String getName(){
	String val = "";
	if (clusteringAlgorithmId == 0)
	    val = AGCTClustering_Algorithm.REFERENCE_NAME_CP;
	else if (clusteringAlgorithmId == 1)
	    val = AGCTClustering_Algorithm.REFERENCE_NAME_KM;
	else if (clusteringAlgorithmId == 2)
	    val = AGCTClustering_Algorithm.REFERENCE_NAME_EM;
	else if (clusteringAlgorithmId == 3)
	    val = AGCTClustering_Algorithm.REFERENCE_NAME_HC;
	else if (clusteringAlgorithmId == 4)
	    val = AGCTClustering_Algorithm.REFERENCE_NAME_AP;
	else if (clusteringAlgorithmId == 5)
	    val = AGCTClustering_Algorithm.REFERENCE_NAME_NM;
	return val;
    }

    public void toClustering(){
	myClusteringAlgorithm.toClustering();
	generate_VV();
	myAGCT.myTabbedPane.myManifoldPane.setMembershipsButtons(true);
	myAGCT.myTabbedPane.myPCAPane.setMembershipsButtons(true);
	myAGCT.myTabbedPane.myCorrelationPane.setMembershipsButtons(true);

	myClusteringAlgorithm.fillGeneMemberships();
    }

    public int getRightNumberSelectedGenes(boolean before){
	if (before)
	    return myAGCT.myDomain.numberSelectedGenesBefore;
	else
	    return myAGCT.myDomain.numberSelectedGenes;
    }

    public int getRightIndexSelectedGeneNumberToGeneNumber(int j, boolean before){
	if (before)
	    return myAGCT.myDomain.selectedGeneNumberToGeneNumberBefore[j];
	else
	    return myAGCT.myDomain.selectedGeneNumberToGeneNumber[j];
    }

    public Gene getRightGeneSelectedGeneNumberToGeneNumber(int j, boolean before){
	if (before)
	    return (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumberBefore[j]);
	else
	    return (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[j]);
    }

    public Pnt getRightClosestCenterPoint(int i, boolean before){
	//To be used only with Prototype
	Gene gg;

	if (before)
	    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(i);
	else
	    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(i);

	if ( (before) && (indexReferenceStringCenters < 3) )
	    Matrix.perror("Coordinates not defined");

	if (indexReferenceStringCenters == 0)
	    return gg.manifold_components_total;
	else if (indexReferenceStringCenters == 1)
	    return gg.manifold_Pnt3D;
	else if (indexReferenceStringCenters == 2)
	    return gg.pca_Pnt3D;
	else if (indexReferenceStringCenters == 3)
	    return gg.finalCoordinates_1D;

	Matrix.perror("No point to choose");

	return null;
    }

    public Pnt getRightPoint(int i, boolean before){
	Gene gg = getRightGeneSelectedGeneNumberToGeneNumber(i, before);

	if ( (before) && (indexReferenceStringCenters < 3) )
	    Matrix.perror("Coordinates not defined");

	if (indexReferenceStringCenters == 0)
	    return gg.manifold_components_total;
	else if (indexReferenceStringCenters == 1)
	    return gg.manifold_Pnt3D;
	else if (indexReferenceStringCenters == 2)
	    return gg.pca_Pnt3D;
	else if (indexReferenceStringCenters == 3)
	    return gg.finalCoordinates_1D;

	Matrix.perror("No point to choose :: indexReferenceStringCenters = " + indexReferenceStringCenters);

	return null;
    }

    public void generate_VV(){
	VV = new Matrix("VV", nclusters, myAGCT.dimFeatures);
	VV.toVV(myAGCT.myDomain, this);
	myAGCT.myInformationFrame.setText("Matrix N (10 x 10 upperleft block):\n" + VV.affiche(10));
    }
}
