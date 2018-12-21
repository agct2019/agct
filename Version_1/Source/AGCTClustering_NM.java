import java.util.*;

 class AGCTClustering_NM extends AGCTClustering_Algorithm{

    int n_Variables;
    //number of columns of the data matrix (typically = number of description variables)

    int n_Method;
    //method to cast matrix V = WH in positive form
    // = 0 : method of Kim and Tidor, NAR 2008 (see below)
    // = 1 : shift my the min

    int t_Update;
    //type of update
    // = 0 : Euclidean update
    // = 1 : KL update

    static int DEFAULT_N_METHOD = 0;
    static int DEFAULT_T_UPDATE = 1;
    static int DEFAULT_N_CLUSTERS = 2;

    Matrix V_NM, W_NM, H_NM;

    public static String name_Method(AGCTClustering_NM a){
	int i = a.n_Method;
	String ret = "Casting entries into positive values following ";
	if (i == 0)
	    ret += "Kim & Tidor's duplicate method (NAR 2008)";
	else if (i == 1)
	    ret += "Min shift for all values";
	ret += ".";
	return ret;
    }

    public static String name_Update(AGCTClustering_NM a){
	int i = a.t_Update;
	String ret = "Update rule minimizes distortion " + Distortion.DIVERGENCES_NAME[i] + ".";
	return ret;
    }

    AGCTClustering_NM(Clustering myc, AGCT mya, int nc, boolean bef, int nMI, int nMet, int tUpd){
	super(myc, mya, nc, bef, nMI);
	myReferenceName = "NM";
	finalClusteringObjectiveFunction = -1;
	V_NM = W_NM = H_NM = null;

	if (nMet == -1)
	    n_Method = DEFAULT_N_METHOD;
	else if (nMet < 0)
	    Matrix.perror("AGCTClustering_NM :: Bad value for the positive defining method");
	else
	    n_Method = nMet;	

	if (tUpd == -1)
	    t_Update = DEFAULT_T_UPDATE;
	else if ( (tUpd != 0) && (tUpd != 1) )
	    Matrix.perror("AGCTClustering_NM :: Bad value for the update rule method");
	else
	    t_Update = tUpd;	

	n_Variables = myClustering.getRightPoint(0, before).coordinates.length;
	//System.out.println("n_Variables = " + n_Variables + ", nMethod = " + nMet);
    }

    public void safeCheck(){
	//Rule taken from "Subsystem Identification Through Dimensionality Reduction of Large-Scale Gene Expression Data"
	//Philip M. Kim and Bruce Tidor, Nucleics Acid Research 2008

	int mm = myClustering.getRightNumberSelectedGenes(before);
	int nn = n_Variables;
	int kk = nclusters;
	boolean passTest = false;

	if ( (n_Method == 0) && ( kk * ( mm + 2 * nn ) < 2 * nn * mm ) )
	    passTest = true;
	else if ( (n_Method != 0) && ( kk * ( mm + nn ) < nn * mm ) )
	    passTest = true;

	if (!passTest){
	    if (AGCT.Debug)
		System.out.println("Too many clusters for the dimension :: switching to full Manifold coordinates");
	    myClustering.indexReferenceStringCenters = 0;
	}
    }

    public void initNMF(){
	safeCheck();
	init_V_NM();
	init_W_NM();
	init_H_NM();
	finalInit();
    }

    public void init_V_NM(){
	int i, j, dimX, dimY;
	double dum, min = -1.0;
	if (n_Method == 0)
	    dimX = 2 * n_Variables;
	else
	    dimX = n_Variables;
	dimY = myClustering.getRightNumberSelectedGenes(before);

	V_NM = new Matrix("Data Matrix (NMF)", dimX, dimY);

	if (n_Method == 0) {
	    for (i=0;i<n_Variables;i++){
		for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++){
		    dum = myClustering.getRightPoint(j, before).coordinates[i];
		    if (dum < 0.0){
			V_NM.coordinates[i][j] = 0.0;
			V_NM.coordinates[i + n_Variables][j] = -dum;
		    }else{
			V_NM.coordinates[i][j] = dum;
			V_NM.coordinates[i + n_Variables][j] = 0.0;
		    }
		}
	    }
	}else if (n_Method == 1) {
	    for (i=0;i<n_Variables;i++){
		for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++){
		    dum = myClustering.getRightPoint(j, before).coordinates[i];
		    if ( ( (i==0) && (j==0) ) || (dum < min) )
			min = dum;
		}
	    }
	    for (i=0;i<n_Variables;i++){
		for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++){
		    dum = myClustering.getRightPoint(j, before).coordinates[i];
		    V_NM.coordinates[i][j] = dum - min;
		}
	    }
	}
    }

    public void init_W_NM(){
	int i, j, dimX, dimY;
	double dum;
	dimX = V_NM.dimX;
	dimY = nclusters;

	W_NM = new Matrix("Meta Variables (NMF)", dimX, dimY);

	for (i=0;i<dimX;i++){
	    for (j=0;j<dimY;j++){
		W_NM.coordinates[i][j] = AGCT.RG.nextDouble();
		if (W_NM.coordinates[i][j] == 0.0)
		    W_NM.coordinates[i][j] = 0.5;

		if (W_NM.coordinates[i][j] < 0.0)
		    W_NM.coordinates[i][j] = -W_NM.coordinates[i][j];
	    }
	}
    }

    public void init_H_NM(){
	int i, j, dimX, dimY;
	double dum;
	dimX = nclusters;
	dimY = myClustering.getRightNumberSelectedGenes(before);

	H_NM = new Matrix("Meta Clustering (NMF)", dimX, dimY);

	for (i=0;i<dimX;i++){
	    for (j=0;j<dimY;j++){
		H_NM.coordinates[i][j] = AGCT.RG.nextDouble();
		if (H_NM.coordinates[i][j] == 0.0)
		    H_NM.coordinates[i][j] = 0.5;

		if (H_NM.coordinates[i][j] < 0.0)
		    H_NM.coordinates[i][j] = -H_NM.coordinates[i][j];
	    }
	}
    }

    public void finalInit(){
	if (n_Method == 0)
	    n_Variables = 2*n_Variables;
    }

    public void plot(JAGCTVisualizationPane vv){
	plotClusterCenters(vv);
    }

    public void plotWithAnnotations(JAGCTVisualizationPane vv){
	plotClusterCentersWithAnnotations(vv);
    }

    public String toString(){
	String val = "Non Negative Matrix Factorization with " + nclusters + " clusters, over " + n_Variables + " description variables.\n";
	val += "Points are " + Clustering.referenceStringCenters[myClustering.indexReferenceStringCenters] + ".\n";
	val += AGCTClustering_NM.name_Method(this) + ".\n" + AGCTClustering_NM.name_Update(this) + ".\n";
	return val;
    }

    public double update_L22(Matrix nextW, Matrix Wt, Matrix WtW, Matrix nextH, Matrix Ht, Matrix HHt, Matrix WtWH, Matrix WtV, Matrix VHt, Matrix WHHt){
	int i, j;
	double dum;

	for (i=0;i<H_NM.dimX;i++)
	    for (j=0;j<H_NM.dimY;j++){
		dum = H_NM.coordinates[i][j] * WtV.coordinates[i][j] / WtWH.coordinates[i][j];
		nextH.coordinates[i][j] = dum;
	    }

	H_NM.copyOnThis(nextH);
	Ht.transpose(H_NM);
	HHt.dot(H_NM,Ht);
	VHt.dot(V_NM,Ht);
	WHHt.dot(W_NM,HHt);

	for (i=0;i<W_NM.dimX;i++)
	    for (j=0;j<W_NM.dimY;j++){
		dum = W_NM.coordinates[i][j] * VHt.coordinates[i][j] / WHHt.coordinates[i][j];
		nextW.coordinates[i][j] = dum;
	    }

	W_NM.copyOnThis(nextW);
	Wt.transpose(W_NM);
	WtW.dot(Wt,W_NM);
	WtWH.dot(WtW,H_NM);
	WtV.dot(Wt,V_NM);

	return Matrix.BregmanError(V_NM, W_NM, H_NM, 0, -1);
    }

    public double update_KL(Matrix nextW, Matrix nextH, Matrix WH){
	int i, j, k;
	double dum, numer, denom;
	for (i=0;i<H_NM.dimX;i++)
	    for (j=0;j<H_NM.dimY;j++){
		numer = 0.0;
		for (k=0;k<V_NM.dimX;k++)
		    numer += (W_NM.coordinates[k][i] * V_NM.coordinates[k][j] / WH.coordinates[k][j]);
		denom = 0.0;
		for (k=0;k<W_NM.dimX;k++)
		    denom += W_NM.coordinates[k][i];
		nextH.coordinates[i][j] = H_NM.coordinates[i][j] * numer / denom;
	    }

	H_NM.copyOnThis(nextH);
	WH.dot(W_NM,H_NM);

	for (i=0;i<W_NM.dimX;i++)
	    for (j=0;j<W_NM.dimY;j++){
		numer = 0.0;
		for (k=0;k<V_NM.dimY;k++)
		    numer += (H_NM.coordinates[j][k] * V_NM.coordinates[i][k] / WH.coordinates[i][k]);
		denom = 0.0;
		for (k=0;k<H_NM.dimY;k++)
		    denom += H_NM.coordinates[j][k];
		nextW.coordinates[i][j] = W_NM.coordinates[i][j] * numer / denom;
	    }

	W_NM.copyOnThis(nextW);
	WH.dot(W_NM,H_NM);

	return Matrix.BregmanError(V_NM, W_NM, H_NM, 1, -1);
    }

    public void NMF(){
	Matrix Wt, Ht, WtW, WtWH, WtV, VHt, HHt, WHHt, WH;

	Wt = new Matrix("Wt", W_NM.dimY, W_NM.dimX);
	Wt.transpose(W_NM);
	Ht = new Matrix("Ht", H_NM.dimY, H_NM.dimX);
	Ht.transpose(H_NM);
	HHt = new Matrix("HHt", H_NM.dimX, Ht.dimY);
	HHt.dot(H_NM,Ht);
	WtW = new Matrix("WtW", Wt.dimX, W_NM.dimY);
	WtW.dot(Wt,W_NM);
	WtWH = new Matrix("WtWH", WtW.dimX, H_NM.dimY);
	WtWH.dot(WtW,H_NM);
	WtV = new Matrix("WtV", Wt.dimX, V_NM.dimY);
	WtV.dot(Wt,V_NM);
	VHt = new Matrix("VHt", V_NM.dimX, Ht.dimY);
	VHt.dot(V_NM,Ht);
	WHHt = new Matrix("WHHt", W_NM.dimX, HHt.dimY);
	WHHt.dot(W_NM,HHt);
	WH = new Matrix("WH", W_NM.dimX, H_NM.dimY);
	WH.dot(W_NM,H_NM);

	Matrix nextW, nextH;

	int i, niter = 0;
	double delta;
	double dprev = -1, dcur, ratio;
	boolean stop = false;
	AGCTCounter ccc;

	if (nMaxIterations>0)
	    ccc = new AGCTCounter(myAGCT.myInformationFrame, "Clustering NM", nMaxIterations);
	else
	    ccc = new AGCTCounter(myAGCT.myInformationFrame, "Clustering NM", 0);

	nextW = new Matrix("Copy W", W_NM.dimX, W_NM.dimY);
	nextH = new Matrix("Copy H", H_NM.dimX, H_NM.dimY);

	do{
	    nextW.copyOnThis(W_NM);
	    nextH.copyOnThis(H_NM);

	    if (t_Update == 0)
		dcur = update_L22(nextW, Wt, WtW, nextH, Ht, HHt, WtWH, WtV, VHt, WHHt);
	    else 
		dcur = update_KL(nextW, nextH, WH);

	    finalClusteringObjectiveFunction = dcur;

	    //System.out.println("NMF distortions :: " + dprev + " --> " + dcur);

	    if (niter > 0){
		if (nMaxIterations > 0)
		    ccc.increment();
		else
		    ccc.setPercent( (int) ( 100.0 * AGCTClustering_Algorithm.EPSILON_KM / (dprev - dcur) ) );
		
		if (dprev < dcur)
		    Matrix.perror("Distortion increase :: " + dprev + " --> " + dcur);
		
		ratio = (dprev - dcur) / dprev;
		//myAGCT.myInformationFrame.appendText(niter + " : " + ratio + "\n");
		if (ratio < AGCTClustering_Algorithm.EPSILON_KM)
		    stop = true;
	    }
	    dprev = dcur;
	    niter++;
	    

	    if (nMaxIterations>0)
		if ( niter >= nMaxIterations )
		    stop = true;

	}while(!stop);
	myAGCT.myInformationFrame.setText("ok.");
	ccc.end();

	toCentersAndMemberships();
	nextW = nextH = Wt = Ht = HHt = WtW = WtWH = WtV = VHt = WHHt = WH = null;
    }

    public void toCentersAndMemberships(){
	hard_memberships_centers = new Vector();

	int [] affect = new int[myClustering.getRightNumberSelectedGenes(before)];
	int i, j, k, imax;
	double tot, vmax = -1.0, vcur;
	Pnt p;

	for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++){
	    imax = -1;
	    for (j=0;j<nclusters;j++){
		vcur = H_NM.coordinates[j][i];
		if (vcur < 0.0)
		    Matrix.perror("AGCTClustering_NM.class :: negative value");
		if ( (j==0) || (vcur > vmax) ){
		    imax = j;
		    vmax = vcur;
		}
	    }
	    affect[i] = imax;
	}

	for (i=0;i<nclusters;i++)
	    if (containsIndex(affect,i)){
		tot = 0.0;
		p = new Pnt(getDimension(before)); 
		for (j=0;j<affect.length;j++)
		    if (affect[j] == i){
			tot += 1.0;
			for (k=0;k<p.coordinates.length;k++)
			    p.coordinates[k] += myClustering.getRightPoint(j, before).coordinates[k];
		    }
		for (k=0;k<p.coordinates.length;k++)
		    p.coordinates[k] /= tot;
		hard_memberships_centers.add(p);
	    }

	for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++)
	    for (j=0;j<nclusters;j++)
		if (affect[i] == j)
		    soft_memberships.coordinates[j][i] = 1.0;
		else
		    soft_memberships.coordinates[j][i] = 0.0;
	    
	soft_memberships_available = true;
	affect = null;
    }


    public void toClustering(){
	initNMF();
	NMF();
	toClusterChoices();
	toSelectedGeneToHardCluster();
	plotCentersAvailable = true;
	hard_clustering_computed = true;
	soft_clustering_computed = false;
	ok_for_statistics = true;
	ControlProcess.put("hardClusteringProcessed",true);
	ControlProcess.put("ClusteringProcessed_NM",true);
    }
}
