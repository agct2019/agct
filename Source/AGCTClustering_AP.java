import java.util.*;
import java.awt.Color;

 class AGCTClustering_AP extends AGCTClustering_Algorithm{

    public static double LAMBDA = 0.5;

    int dimRef;

    int vP;
    double vB;
    double sKK; //s(k,k)

    Matrix P; //computes s(i,j)
    double [][] respons;
    double [][] avail;
    double [][] randa;

    double [] maxAS1, maxAS2;
    double [] maxR;
    boolean [] secondFound;

    int [] indexAS1;

    int[] exemplar;
    int[] clusterNumber;

    boolean readyToSave;

    AGCTClustering_AP(Clustering myc, AGCT mya, int valp, double valb, boolean bef, int nMI){
	super(myc, mya, -1, bef, nMI);
	myReferenceName = "AP";
	seedingType = -1;
	finalClusteringObjectiveFunction = -1;
	vP = valp;
	vB = valb;
	P = null;
	respons = null;
	avail = null;
	randa = null;
	exemplar = null;
	clusterNumber = null;
	readyToSave = false;

	maxAS1 = maxAS2 = maxR = null;
	indexAS1 = null;
	secondFound = null;
    }

    public String getAdditionalSavingData(){
	if (!readyToSave)
	    return super.getAdditionalSavingData();

	String val = "Exemplars for clustering:\n";
	
	int i, j;
	boolean found = false;
	for (i=0;i<nclusters;i++){
	    j=0;
	    while(!found){
		if (clusterNumber[j] == i)
		    found = true;
		else
		    j++;
	    }
	    val += " - for cluster " + i + " : " + myClustering.getRightGeneSelectedGeneNumberToGeneNumber(exemplar[j], before) + "\n";   
	}	
	val += "\n";

	return val;
    }

    public static double damp(double newval, double lastval){
	return ( (LAMBDA * lastval) + ( (1 - LAMBDA) * newval) );
    }

    public void initAvail(){
	int i, k, iN = myClustering.getRightNumberSelectedGenes(before);
	double maj;

	defaultMaxAS();

	avail = new double[iN][];
	for (i=0;i<iN;i++)
	    avail[i] = new double[iN];

	for (i=0;i<iN;i++)
	    for (k=0;k<iN;k++){
		avail[i][k] = 0.0;

		maj = P.coordinates[i][k];
		updateMaxAS(maj, i, k);
	    }

	randa = new double[iN][];
	for (i=0;i<iN;i++)
	    randa[i] = new double[iN];

	for (i=0;i<iN;i++)
	    for (k=0;k<iN;k++)
		randa[i][k] = 0.0;

    }

    public void initRespons(){
	int i, j, iN = myClustering.getRightNumberSelectedGenes(before);
	respons = new double[iN][];
	for (i=0;i<iN;i++)
	    respons[i] = new double[iN];
    }

    public void initExemplars(){
	exemplar = new int[myClustering.getRightNumberSelectedGenes(before)];
    }

    public void initMax(){
	int i, iN = myClustering.getRightNumberSelectedGenes(before);
	maxAS1 = new double[iN];
	maxAS2 = new double[iN];
	maxR = new double[iN];

	indexAS1 = new int[iN];
	secondFound = new boolean[iN];
	defaultMaxAS();
    }

    public void defaultMaxAS(){
	int i, iN = myClustering.getRightNumberSelectedGenes(before);
	for (i=0;i<iN;i++){
	    maxAS1[i] = maxAS2[i] = 0.0;
	    indexAS1[i] = -1;
	    secondFound[i] = false;
	}
    }

    public void initEverything(){
	computeSKK();
	toPref();
	initMax();
	initAvail();
	initRespons();
	initExemplars();
    }

    public void plot(JAGCTVisualizationPane vv){
	plotClusterCenters(vv);
	plotAffinities(vv);
    }

    public void plotWithAnnotations(JAGCTVisualizationPane vv){
	plotClusterCentersWithAnnotations(vv);
	plotAffinities(vv);
    }

    public void plotAffinities(JAGCTVisualizationPane vv){
	Color c;
	int i, iN = myClustering.getRightNumberSelectedGenes(before), ncl;
	View3D v = vv.visualizationPanel.myView3D;
	Pnt3D pe, p;
	if ( plotCentersAvailable )
	    if (vv.myReferenceName.equals(JAGCTVisualizationPane.M_P) || (vv.myReferenceName.equals(JAGCTVisualizationPane.P_P)) ){
		for (i=0;i<iN;i++){
		    ncl = clusterNumber[i];
		    if (vv.structPlottedCluster(ncl)){
			pe = vv.getRightPoint(exemplar[i]);
			p = vv.getRightPoint(i);
			v.drawShadowLineVsRef(pe, p, JAGCTGraphicsPanel.minColor, myClustering.colorCluster[ncl], false);
		    }
		}
	    }
    }

    public String toString(){
	String val = "Affinity Propagation ";
	if (nclusters <= 0)
	    val += "(clusters not computed)\n";
	else
	    val += "with " + nclusters + " clusters\n";
	val += "Initial values for Preference: " + vP + ", for Bregman p-norm divergence with p=" + vB + "\n";
	val += "Points are " + Clustering.referenceStringCenters[myClustering.indexReferenceStringCenters] + ".";
	return val;
    }

    public void computeSKK(){
	int i, j, k, iN = myClustering.getRightNumberSelectedGenes(before), pp, lim;

    	double [] simil = new double [iN * (iN - 1) / 2];
	for (i=0;i<iN * (iN - 1) / 2;i++)
	    simil[i] = 0.0;

	k = 0;
	for (i=0;i<iN-1;i++)
	    for (j=i+1;j<iN;j++){
		simil[k] = - Distortion.Bregman(myClustering.getRightPoint(i, before),
						myClustering.getRightPoint(j, before),
					      4,
					      vB);
		k++;
	    }

	QuickSort.quicksort(simil);

	if (vP == -1)
	    sKK = simil[0];
	else{
	    if (vP == -2)
		pp = 50;
	    else{
		if (vP<1)
		    Matrix.perror("invalid value for P");
		pp = vP;
	    }

	    if (iN < 1000)
		lim = ( (iN * (iN - 1) / 2) * pp ) / 100;
	    else
		lim = (iN * (iN - 1) / 200) * pp;

	    try{
		sKK = simil[lim];

	    }catch(ArrayIndexOutOfBoundsException ee){
		System.out.println("Affinity propagation : problem with index " + lim + " with iN = " + iN + " and pp = " + pp);
		System.out.println(" iN * (iN - 1) / 2 = " + (iN * (iN - 1) / 2) );
		System.out.println(" (iN * (iN - 1) / 2) * pp = " + ( (iN * (iN - 1) / 2) * pp ) );

		lim = 0;
		sKK = simil[lim];
	    }
	}

	//System.out.println("sKK = " + sKK + " (vP = " + vP + ")");
    }


    public void toPref(){
	P = new Matrix("Preferences", myClustering.getRightNumberSelectedGenes(before), myClustering.getRightNumberSelectedGenes(before));
	P.toP(myAGCT.myDomain, myClustering, vB, sKK);
    }

    public void updateRespons(int it){
	int i, k, kp, km = -1, iN = myClustering.getRightNumberSelectedGenes(before);
	double dum;

	for (i=0;i<iN;i++)
	    maxR[i] = 0.0;

	for (i=0;i<iN;i++){
	    for (k=0;k<iN;k++){

		if (indexAS1[i] == k)
		    dum = maxAS2[i];
		else
		    dum = maxAS1[i];
		
		if (it > 0)
		    respons[i][k] = AGCTClustering_AP.damp(respons[i][k], P.coordinates[i][k] - dum);
		else
		    respons[i][k] = P.coordinates[i][k] - dum;

		maxR[k] += Math.max(0.0, respons[i][k]);
	    }
	}
    }

    public void updateMaxAS(double maj, int i, int k){
	if ( (indexAS1[i] == -1) || (maj > maxAS1[i]) || ( (maj == maxAS1[i]) && (indexAS1[i] != k) ) ){
	    if ( (indexAS1[i] != -1) && ( (maj > maxAS2[i]) || (!secondFound[i]) ) ){
		maxAS2[i] = maxAS1[i];
		secondFound[i] = true;
	    }
	    maxAS1[i] = maj;
	    indexAS1[i] = k;
	}else if ( (indexAS1[i] != -1) && ( (maj > maxAS2[i]) || (!secondFound[i]) ) ){
	    maxAS2[i] = maj;
	    secondFound[i] = true;
	}
    }

    public void updateAvail(int it){
	int i, k, ip, iN = myClustering.getRightNumberSelectedGenes(before);
	double dum, maj;
	defaultMaxAS();

	for (i=0;i<iN;i++){
	    for (k=0;k<iN;k++){

		if (k!=i)
		    dum = respons[k][k] + maxR[k] - Math.max(0.0, respons[i][k]) - Math.max(0.0, respons[k][k]);
		else
		    dum = maxR[k] - Math.max(0.0, respons[k][k]);

		if (k != i)
		    dum = Math.min(0.0, dum);

		if (it > 0)
			avail[i][k] = AGCTClustering_AP.damp(avail[i][k], dum);
		else
			avail[i][k] = dum;

		maj = avail[i][k] + P.coordinates[i][k];
		updateMaxAS(maj, i, k);
	    }
	}

	/*for (i=0;i<iN;i++)
	  System.out.print("  MaxAS1[" + i + "] = " +maxAS1[i]);
	  System.out.println("");
	  for (i=0;i<iN;i++)
	  System.out.print("  indexAS1[" + i + "] = " + indexAS1[i]);
	  System.out.println("");
	  for (i=0;i<iN;i++)
	  System.out.print("  MaxAS2[" + i + "] = " +maxAS2[i]);
	  System.out.println("");
	  for (i=0;i<iN;i++)
	  System.out.print("  secondFound[" + i + "] = " + secondFound[i]);
	  System.out.println("");*/

    }

    public void updateRanda(int it){
	int i, k, iN = myClustering.getRightNumberSelectedGenes(before);
	for (i=0;i<iN;i++)
	    for (k=0;k<iN;k++)
		randa[i][k] = avail[i][k] + respons[i][k];

    }

    public void updateExemplars(){
	int i, k, iN = myClustering.getRightNumberSelectedGenes(before);
	double maxar;

	for (i=0;i<iN;i++){
	    maxar = 0.0;
	    for (k=0;k<iN;k++)
		if ( (k==0) || (maxar < randa[i][k]) ){
		    maxar = randa[i][k];
		    exemplar[i] = k;
		}
	}

	/*for (i=0;i<iN;i++)
	  System.out.print(" [" + i + " --> " + exemplar[i] + "] ");
	  System.out.println("");*/
    }

    public void computeClusterNumber(){
	int i, ref, replace, ni = -1, iN = myClustering.getRightNumberSelectedGenes(before);
	clusterNumber = new int[iN];
	boolean larger = true, foundcluster = false;
	
	for (i=0;i<iN;i++)
	    clusterNumber[i] = exemplar[i];
	
	ref = replace = 0;
	nclusters = 0;

	for (i=0;i<iN;i++)
	    if ( (i==0) || (clusterNumber[i] < replace) )
		replace = clusterNumber[i];

	do{
	    foundcluster = false;
	    ni = iN;
	    for (i=0;i<iN;i++){
		if ( (!foundcluster) && (clusterNumber[i] == replace) ){
		    foundcluster = true;
		    nclusters++;
		}
		if ( (clusterNumber[i] > replace) && (clusterNumber[i] < ni) )
		    ni = clusterNumber[i];

		if (clusterNumber[i] == replace)
		    clusterNumber[i] = ref;
	    }
	    if ( (ni > replace) && (ni <iN) ){
		replace = ni;
		ref++;
	    }else
		larger = false;
	}while(larger == true);

	myClustering.nclusters = nclusters;
    }

    public void toMemberships(){
	initSoftMemberships();
	int i, j;
	for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++){
	    for (j=0;j<nclusters;j++){
		soft_memberships.coordinates[j][i] = 0.0;
		if (clusterNumber[i] == j)
		    soft_memberships.coordinates[j][i] = 1.0;
	    }    
	}
	soft_memberships_available = true;
    }

    public void saveSpace(){
	P = null;
	respons = avail = randa = null;
	maxAS1 = maxAS2 = maxR = null;
	secondFound = null;
	indexAS1 = null;
    }

    public void AP(){
	boolean stop = false;
	int i, niter = 0, iN = myClustering.getRightNumberSelectedGenes(before), nid = 0, ham;
	int [] lastexemplar = new int[iN];

	AGCTCounter ccc;

	if (nMaxIterations>0)
	    ccc = new AGCTCounter(myAGCT.myInformationFrame, "Clustering AP", nMaxIterations);
	else
	    ccc = new AGCTCounter(myAGCT.myInformationFrame, "Clustering AP", Affinity_Propagation_Id_Max);

	do{
	    updateRespons(niter);
	    updateAvail(niter);
	    updateRanda(niter);
	    updateExemplars();
	    //ccc.increment();

	    if (niter>0){
		ham = Matrix.Hamming(lastexemplar, exemplar);
		if (ham == 0)
		    nid ++;
		else
		    nid = 0;

		//System.out.println("Iter = " + iter + " ; Hamming = " + ham);
	    }

	    if (nMaxIterations > 0)
		ccc.increment();
	    else
		ccc.setBound(nid);

	    for (i=0;i<iN;i++)
		lastexemplar[i] = exemplar[i];

	    niter++;
	    if (nMaxIterations>0)
		if ( niter >= nMaxIterations )
		    stop = true;

	    if (nid >= Affinity_Propagation_Id_Max)
		stop = true;

	}while (!stop);
	ccc.end();

	plotCentersAvailable = true;
    }

    public void toClustering(){
	dimRef = getDimension(before);

	initEverything();
	AP();

	computeClusterNumber();
	myClustering.toClusterColors();

	toMemberships();
	toClusterChoices();

	soft_clustering_computed = true;
	readyToSave = true;

	toHardMemberships();

	toSelectedGeneToHardCluster();

	hard_clustering_computed = true;
	ok_for_statistics = true;
	ControlProcess.put("hardClusteringProcessed",true);
	ControlProcess.put("ClusteringProcessed_AP",true);
    }
}
