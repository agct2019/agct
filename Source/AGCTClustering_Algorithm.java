import java.util.*;
import java.awt.Color;

 class AGCTClustering_Algorithm implements Debuggable{

    static final double EPSILON_KM = 10E-10;
    static final double EPSILON_EM = 10E-20;

    static final String REFERENCE_NAME_KM = "KM", 
	REFERENCE_NAME_EM = "EM",  
	REFERENCE_NAME_CP = "CP",  
	REFERENCE_NAME_HC = "HC",  
	REFERENCE_NAME_AP = "AP",
	REFERENCE_NAME_NM = "NM";  

    static String ALL_CLUSTERS = "All";
    static String CLUSTER_NAME = "C";

    static int N_MAX_ITERATIONS_AP = 1000;
    // to fasten AP...

    String myReferenceName;
    // One of KM, CP, EM, HC, AP, NM (see above)

    Clustering myClustering;
    AGCT myAGCT;

    int nclusters, nMaxIterations;
    int seedingType;
    // -1 = N/A
    // 0 = Forgy
    // 1 = Arthur - Vassilvitskii
    // 2 = Bregman

    public static String SEEDING_NAME[] = {"Forgy", "Arthur-Vassilvitskii", "Bregman (general AV)"};

    double finalClusteringObjectiveFunction;
    double [][] referencesByCluster;


    Matrix soft_memberships;
    Vector hard_memberships_centers;
    //hard_memberships_centers = used in clustering algorithms

    Vector hard_memberships_centers_manifold_total, hard_memberships_centers_manifold_Pnt3D;
    //hard_memberships_centers_manifold_total = computed from soft_memberships for display
    //hard_memberships_centers_manifold_Pnt3D = keeps the current 3D view of cluster centers (faster than making computation each time !)

    Vector hard_memberships_centers_pca_total, hard_memberships_centers_pca_Pnt3D;
    //hard_memberships_centers_pca_total = computed from soft_memberships for display
    //hard_memberships_centers_pca_Pnt3D = keeps the current 3D view of cluster centers (faster than making computation each time !)

    boolean soft_clustering_computed, hard_clustering_computed, soft_memberships_available, ok_for_statistics;

    boolean plotCentersAvailable;

    boolean before;
    // true iff we make clutering on data before prototype selection

    String [] clusterChoices;

    int [] selectedGeneToHardCluster;
    //gives for each gene the most probable cluster to which it belongs
    int [] clusterSizesPrototypes;
    int [] clusterSizesWhole;
    //Sizes of clusters as prototypes or for all points

    double [] hard_average_distance;
    double [] hard_min_distance;
    double [] hard_max_distance;
    double [] hard_median_distance;
    // If a distance matrix is given, computes the average per-cluster distance matrix;
    boolean averageDistanceComputed;

    double [][] dimension;
    // computes the dimension of each cluster (X), for each K (Y)
    double [] average_dimension_clusters;
    // computes the average dimension of each cluster
    Vector average_dimension_histogram;
    // used to store histograms of dimension
    public static int NUMBER_BINS = 50;

     public static int DELTA_DIMENSION = 5; //XXX 10
    // computes the HARMONIC average dimension between DELTA_DIMENSION and AGCT.Max_K_Dimension_Estimation - DELTA_DIMENSION

    AGCTClustering_Algorithm(Clustering myc, AGCT mya, int nc, boolean bef, int nMI){
	myReferenceName = "";

	nclusters = nc;
	nMaxIterations = nMI;
	myClustering = myc;
	myAGCT = mya;

	soft_memberships = null;
	hard_memberships_centers = null;
	hard_memberships_centers_manifold_total = null;
	hard_memberships_centers_manifold_Pnt3D = null;

	hard_memberships_centers_pca_total = null;
	hard_memberships_centers_pca_Pnt3D = null;

	referencesByCluster = null;

	dimension = null;
	average_dimension_clusters = null;

	soft_clustering_computed = false;
	hard_clustering_computed = false;
	plotCentersAvailable = false;
	soft_memberships_available = false;
	ok_for_statistics = false;

	hard_average_distance = hard_min_distance = hard_median_distance = hard_max_distance = null;
	averageDistanceComputed = false;

	before = bef;
    }

    public void toClusteringLite(Vector allPoints){
	Matrix.perror("AGCTClustering_Algorithm.class :: no such algorithm");
    }

    public void saveSpace(){
    }

    public String getAdditionalSavingData(){
	return "None";
    }

    public String getDistanceInformation(){
	String ret = "";
	int i;

	ret += "(Cl #\tSize\tMin\tAvg\tMed\tMax)\n";
	for (i=0;i<nclusters;i++)
	    ret += "C# " + i + "\t" + clusterSizesPrototypes[i] + "\t" + Matrix.DF.format(hard_min_distance[i]) + "\t" + Matrix.DF.format(hard_average_distance[i]) + "\t" + Matrix.DF.format(hard_median_distance[i]) + "\t" + Matrix.DF.format(hard_max_distance[i]) + "\n";

	return ret;
    }

    public double returnFinalObjectiveFunction(){
	return finalClusteringObjectiveFunction;
    }

    public String distances(int i){
	return "" + Matrix.DF.format(hard_min_distance[i]) + " < " + Matrix.DF.format(hard_average_distance[i]) + " (" +  Matrix.DF.format(hard_median_distance[i]) + ") < " + Matrix.DF.format(hard_max_distance[i]);
    }

    public static String dimension(double [][] dim, int nc, int nk){
	if (dim == null)
	    return "?";
	return "" + ((int) dim[nc][nk]);
    }
    
    public static String average_dimension_clusters(double [] dim, int nc){
	if (dim == null)
	    return "?";
	return "" + DF.format(dim[nc]);
    }
    
    public void plotClusterCenters(JAGCTVisualizationPane vv){
	Color c;
	int i;
	View3D v = vv.visualizationPanel.myView3D;
	String sz;
	if ( plotCentersAvailable ){

	    if (vv.myReferenceName.equals(JAGCTVisualizationPane.M_P)){
		softMembershipsToHardCenters_Manifold_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
		for (i=0;i<nclusters;i++)
		    if (vv.structPlottedCluster(i)){
			v.g.setColor(myClustering.colorCluster[i]);

			sz = "" + AGCTClustering_Algorithm.average_dimension_clusters(average_dimension_clusters,i) + "Davg_" + AGCTClustering_Algorithm.dimension(dimension,i,AGCT.K_Dimension_Estimation) + "D, ";

			if (Prototype.No_Reduction)
			    sz += "" + clusterSizesPrototypes[i] + "G";
			else
			    sz += "" + clusterSizesPrototypes[i] + "P|" + clusterSizesWhole[i] + "G";

			if ( ( (AGCT.Filter_Similarities_With_Distances) || (AGCT.Filter_Triangulation_With_Distances) ) && (averageDistanceComputed) )
			    sz += "|" + distances(i);

			v.drawCenterVsRef((Pnt3D) hard_memberships_centers_manifold_Pnt3D.elementAt(i), i, sz);
		    }
	    }else if (vv.myReferenceName.equals(JAGCTVisualizationPane.P_P)){
		softMembershipsToHardCenters_Pca_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
		for (i=0;i<nclusters;i++)
		    if (vv.structPlottedCluster(i)){
			v.g.setColor(myClustering.colorCluster[i]);
			if (Prototype.No_Reduction)
			    sz = "" + clusterSizesPrototypes[i];
			else
			    sz = "" + distances(i);

			if ( ( (AGCT.Filter_Similarities_With_Distances) || (AGCT.Filter_Triangulation_With_Distances) ) && (averageDistanceComputed) )
			    sz += "|" + Matrix.DF.format(hard_average_distance[i]);

			v.drawCenterVsRef((Pnt3D) hard_memberships_centers_pca_Pnt3D.elementAt(i), i, sz);
		    }
	    }
	}
    }
    
    public void plotClusterCentersWithAnnotations(JAGCTVisualizationPane vv){
	Color c;
	int i;
	View3D v = vv.visualizationPanel.myView3D;
	String sz;
	if ( plotCentersAvailable ){

	    if (vv.myReferenceName.equals(JAGCTVisualizationPane.M_P)){
		softMembershipsToHardCenters_Manifold_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
		for (i=0;i<nclusters;i++)
		    if (vv.structPlottedCluster(i)){
			v.g.setColor(myClustering.colorCluster[i]);
			if (Prototype.No_Reduction)
			    sz = "" + clusterSizesPrototypes[i];
			else
			    sz = "" + clusterSizesPrototypes[i] + "|" + clusterSizesWhole[i];
			v.drawCenterVsRefWithAnnotations((Pnt3D) hard_memberships_centers_manifold_Pnt3D.elementAt(i), i, sz);
		    }
	    }else if (vv.myReferenceName.equals(JAGCTVisualizationPane.P_P)){
		softMembershipsToHardCenters_Pca_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
		for (i=0;i<nclusters;i++)
		    if (vv.structPlottedCluster(i)){
			v.g.setColor(myClustering.colorCluster[i]);
			if (Prototype.No_Reduction)
			    sz = "" + clusterSizesPrototypes[i];
			else
			    sz = "" + clusterSizesPrototypes[i] + "|" + clusterSizesWhole[i];
			v.drawCenterVsRefWithAnnotations((Pnt3D) hard_memberships_centers_pca_Pnt3D.elementAt(i), i, sz);
		    }
	    }
	}
    }
    

    public int getDimension(boolean before){
	return myClustering.getRightPoint(0, before).coordinates.length;
    }

    public void initSoftMemberships(){
	soft_memberships = new Matrix("Soft_Memberships_" + myReferenceName,nclusters,myClustering.getRightNumberSelectedGenes(before)); //Nrows = number of clusters
	int i, j;
	
	for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++){
	    for (j=0;j<nclusters;j++)
		soft_memberships.coordinates[j][i] = 0.0;
	}
    }

    public void initHardMemberships(){
	hard_memberships_centers = new Vector();
    }

    public void initHardMembershipForgy(){
	int i, j, k;
	boolean deja;
	int [] indexes_to_pick = new int [nclusters];

	for (i=0;i<nclusters;i++)
	    indexes_to_pick[i] = 0;

	for (i=0;i<nclusters;i++){
	    do{
		j = AGCT.RG.nextInt(myClustering.getRightNumberSelectedGenes(before));
		k=0;
		deja= false;
		do{
		    if (indexes_to_pick[k] == j){
			deja = true;
		    }
		    k++;
		}while ( (k<i) && (deja == false) );
	    }while (deja == true);
	    indexes_to_pick[i] = j;
	}

	for (i=0;i<nclusters;i++)
	    hard_memberships_centers.add(new Pnt(myClustering.getRightPoint(indexes_to_pick[i], before))); //Correction !
    }

    public void softMembershipsToHardCenters_Manifold_Total(){
	int i, j, k, dim = AGCT.Number_Of_Manifold_Components;
	hard_memberships_centers_manifold_total = new Vector();
	Pnt p;
	Gene gg;
	double w;
	for (i=0;i<nclusters;i++){
	    p = new Pnt(dim);
	    for (k=0;k<dim;k++)
		p.coordinates[k] = 0.0;
	    w = 0.0;
	    for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++){
		gg = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(j, before);
		w += soft_memberships.coordinates[i][j];
		for (k=0;k<dim;k++){
		    p.coordinates[k] += (soft_memberships.coordinates[i][j] * gg.manifold_components_total.coordinates[k]);
		}
	    }

	    if (w != 0.0)
		for (k=0;k<dim;k++)
		    p.coordinates[k] /= w;
	    else
		for (k=0;k<dim;k++)
		    p.coordinates[k] = 0.0;

	    hard_memberships_centers_manifold_total.add(p);
	}
    }

    public void softMembershipsToHardCenters_Pca_Total(){
	int i, j, k, dim = myAGCT.dimFeatures;
	hard_memberships_centers_pca_total = new Vector();
	Pnt p;
	Gene gg;
	double w;
	for (i=0;i<nclusters;i++){
	    p = new Pnt(dim);
	    for (k=0;k<dim;k++)
		p.coordinates[k] = 0.0;
	    w = 0.0;
	    for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++){
		gg = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(j, before);
		w += soft_memberships.coordinates[i][j];
		for (k=0;k<dim;k++){
		    p.coordinates[k] += (soft_memberships.coordinates[i][j] * gg.pca_components.coordinates[k]);
		}
	    }

	    if (w != 0.0)
		for (k=0;k<dim;k++)
		    p.coordinates[k] /= w;
	    else
		for (k=0;k<dim;k++)
		    p.coordinates[k] = 0.0;

	    hard_memberships_centers_pca_total.add(p);
	}
    }


    public void softMembershipsToHardCenters_Manifold_Pnt3D(int xAxis, int yAxis, int zAxis){
	int i;
	if (hard_memberships_centers_manifold_Pnt3D == null){
	    hard_memberships_centers_manifold_Pnt3D = new Vector();
	    for (i=0;i<nclusters;i++)
		hard_memberships_centers_manifold_Pnt3D.add(new Pnt3D());
	}
	if (hard_memberships_centers_manifold_total == null)
	    softMembershipsToHardCenters_Manifold_Total();

	Pnt p;
	for (i=0;i<nclusters;i++){
	    p = (Pnt3D) hard_memberships_centers_manifold_Pnt3D.elementAt(i);
	    p.coordinates[0] = ( (Pnt) hard_memberships_centers_manifold_total.elementAt(i) ).coordinates[xAxis];
	    p.coordinates[1] = ( (Pnt) hard_memberships_centers_manifold_total.elementAt(i) ).coordinates[yAxis];
	    p.coordinates[2] = ( (Pnt) hard_memberships_centers_manifold_total.elementAt(i) ).coordinates[zAxis];
	}
    }

    public void softMembershipsToHardCenters_Pca_Pnt3D(int xAxis, int yAxis, int zAxis){
	int i;
	if (hard_memberships_centers_pca_Pnt3D == null){
	    hard_memberships_centers_pca_Pnt3D = new Vector();
	    for (i=0;i<nclusters;i++)
		hard_memberships_centers_pca_Pnt3D.add(new Pnt3D());
	}
	if (hard_memberships_centers_pca_total == null)
	    softMembershipsToHardCenters_Pca_Total();

	Pnt p;
	for (i=0;i<nclusters;i++){
	    p = (Pnt3D) hard_memberships_centers_pca_Pnt3D.elementAt(i);
	    p.coordinates[0] = ( (Pnt) hard_memberships_centers_pca_total.elementAt(i) ).coordinates[xAxis];
	    p.coordinates[1] = ( (Pnt) hard_memberships_centers_pca_total.elementAt(i) ).coordinates[yAxis];
	    p.coordinates[2] = ( (Pnt) hard_memberships_centers_pca_total.elementAt(i) ).coordinates[zAxis];
	}
    }

    public void toClustering(){
	//does Nothing !
	Matrix.perror("Instance of AGCTClustering_Algorithm badly defined");
    }

    public void toClusterChoices(){
	clusterChoices = new String[nclusters + 1];
	clusterChoices[0] = AGCTClustering_Algorithm.ALL_CLUSTERS;
	int i;
	for (i=0;i<nclusters;i++)
	    clusterChoices[i+1] = new String(AGCTClustering_Algorithm.CLUSTER_NAME + i);
    }

    public void compute_dimension(){
	Gene gg;
	int i, j, k;
	boolean allpos, dimcomp = false;
	Vector dumv = null, dumh = null;
	double [] dumd;
	double mindim = -1.0, maxdim = -1.0, curdim;

	if (soft_memberships == null)
	    Matrix.perror("Cannot compute dimensions : soft-memberships = null");
	
	double [] average_dimension = new double[nclusters];
	double [] total_weight = new double[nclusters];

	if (hard_clustering_computed){
	    average_dimension_histogram = new Vector();
	    dumh = new Vector();
	    for (i=0;i<nclusters;i++){
		average_dimension_histogram.addElement(new Vector());
		dumh.addElement(new Vector());
	    }
	}
	
	dimension = new double [nclusters][];
	for (i=0;i<nclusters;i++)
	    dimension[i] = new double[AGCT.Max_K_Dimension_Estimation];

	double tot = 0.0;

	for (i=0;i<nclusters;i++){
	    average_dimension[i] = 0.0;
	    total_weight[i] = 0.0;
	    for (j=0;j<AGCT.Max_K_Dimension_Estimation;j++)
		dimension[i][j] = 0.0;
	}

	for (k=0;k<AGCT.Max_K_Dimension_Estimation;k++){
	    allpos = true;
	    i = 0;
	    do{
		gg = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(i, before);
		if (gg.local_dimension[k] == -1.0)
		    allpos = false;
		i++;
	    }while( (i<myClustering.getRightNumberSelectedGenes(before)) && (allpos) );

	    if (!allpos){
		for (j=0;j<nclusters;j++)
		    dimension[j][k] = -1.0;
	    }else{
		for (i=0;i<nclusters;i++){
		    average_dimension[i] = 0.0;
		    total_weight[i] = 0.0;
		}

		for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++){
		    tot = 0.0;
	     
		    gg = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(i, before);
		    if ( (!dimcomp) && (hard_clustering_computed) ){
			curdim = gg.averageDimension;

			dumv = (Vector) dumh.elementAt(selectedGeneToHardCluster[i]);
			dumv.addElement(new Double(curdim));

			if ((i==0) || (curdim > maxdim))
			    maxdim = curdim;
			if ((i==0) || (curdim < mindim))
			    mindim = curdim;
		    }

		    for (j=0;j<nclusters;j++){
			tot += soft_memberships.coordinates[j][i];
			average_dimension[j] += (soft_memberships.coordinates[j][i] * gg.local_dimension[k]);
			total_weight[j] += soft_memberships.coordinates[j][i];
		    }
		    
		    if (tot != 1.0)
			Matrix.perror("soft memberships do not sum to one");
		}
		dimcomp = true;

		for (j=0;j<nclusters;j++){
		    if (total_weight[j] == 0.0)
			Matrix.perror("Cluster " + total_weight[j] + " has zero total weight");
		    dimension[j][k] = average_dimension[j] / total_weight[j];
		}
	    }
	} 

	average_dimension_clusters = new double [nclusters]; 
	for (j=0;j<nclusters;j++){
	    average_dimension_clusters[j] = 0.0;
	}
	if (AGCTClustering_Algorithm.DELTA_DIMENSION > AGCT.Max_K_Dimension_Estimation - AGCTClustering_Algorithm.DELTA_DIMENSION - 2)
	    Matrix.perror("AGCTClustering_Algorithm.class :: computation of the dimension impossible --- K too small");
	    
	for (j=0;j<nclusters;j++){
	    tot = 0.0;
	    for (k=AGCTClustering_Algorithm.DELTA_DIMENSION;k<AGCT.Max_K_Dimension_Estimation - AGCTClustering_Algorithm.DELTA_DIMENSION;k++)
		tot += (1.0 / dimension[j][k]);
	    tot = 1.0 / tot;
	    tot *= ( (double) (AGCT.Max_K_Dimension_Estimation - 2 * AGCTClustering_Algorithm.DELTA_DIMENSION) );
	    average_dimension_clusters[j] = tot;
	}

	if (hard_clustering_computed){
	    double step, begd, endd, dd;
	    int card = 0;
	    Vector vv, vv2;

	    step = ( maxdim - mindim ) / AGCTClustering_Algorithm.NUMBER_BINS;

	    for (i=0;i<nclusters;i++){
		begd = mindim;
		endd = mindim + step;

		dumv = (Vector) dumh.elementAt(i);
		dumd = new double [dumv.size()];
		for (j=0;j<dumv.size();j++)
		    dumd[j] = ((Double) dumv.elementAt(j)).doubleValue();
		QuickSort.quicksort(dumd);
		for (k=0;k<AGCTClustering_Algorithm.NUMBER_BINS;k++){
		    card = 0;
		    for (j=0;j<dumd.length;j++){
			dd = dumd[j];
			if ( (dd >= begd) && (dd <= endd) )
			    card++;
		    }
		    vv = new Vector();
		    vv.addElement(new Double((begd + endd) / 2.0));
		    vv.addElement(new Integer(card));

		    vv2 = (Vector) average_dimension_histogram.elementAt(i);
		    vv2.addElement(vv);
		    
		    begd = endd;
		    endd += step;
		}
	    }
	}
    }

    public void fillGeneMemberships(){
	Gene gg;
	int i, j;
	Double [] memb;

	if (soft_memberships != null){
	    for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++){
		gg = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(i, before);
		memb = new Double [nclusters];
		for (j=0;j<nclusters;j++)
		    memb[j] = new Double(soft_memberships.coordinates[j][i]);

		gg.addClusterMemberships(memb);
	    }
	}
	memb = null;

	if (AGCT.LOCAL_DIMENSION_COMPUTED)
	    compute_dimension();
    }

    public double totalDistortion(){
	if ( (hard_memberships_centers == null) || (!soft_memberships_available) )
	    Matrix.perror("No hard membership distortion available");

	double val = 0.0;
	int i,j;
	for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++)
	    for (j=0;j<nclusters;j++)
		val += (soft_memberships.coordinates[j][i] * Distortion.distortion_l22(myClustering.getRightPoint(i, before), (Pnt) hard_memberships_centers.elementAt(j)));
	return val;
    }

    public double averageIntraClusterSimilarity(){
	if (!soft_memberships_available)
	    Matrix.perror("soft_memberships unavailable");
	Gene gi, gj;
	
	int i,j,k,sizz;
	double simK, val = 0.0;
	for (k=0;k<nclusters;k++){
	    simK = 0.0;
	    sizz = 0;
	    for (i=0;i<myClustering.getRightNumberSelectedGenes(before)-1;i++)
		for (j=i+1;j<myClustering.getRightNumberSelectedGenes(before);j++)
		    if ( (majorityClusterFor(i,k)) && (majorityClusterFor(j,k)) ){
			sizz++;
			gi = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(i, before);
			gj = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(j, before);
			simK += Gene.getSimilarity(gi, gj);
		    }
	    if (sizz == 0)
		Matrix.perror("empty cluster");
	    simK /= (double) sizz;
	    val += simK;
	}
	val /= (double) nclusters;
	return val;
    }

    public boolean containsIndex(int [] affect, int id){
	int i;
	for (i=0;i<affect.length;i++)
	    if (affect[i] == id)
		return true;
	return false;
    }

    public double kernelSim(){
	//Cf Zass + Shashua ICCV 06, (2)
	if (!soft_memberships_available)
	    Matrix.perror("soft_memberships unavailable");
	Gene gi, gj;
	
	int i,j,k,sizz,neff = 0;
	double simK, val = 0.0;
	for (k=0;k<nclusters;k++){
	    simK = 0.0;
	    sizz = 0;
	    for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++){
		if (majorityClusterFor(i,k)){
		    sizz++;
		    gi = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(i, before);
		    for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++){
			if (majorityClusterFor(j,k)){
			    gj = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(j, before);
			    simK += Gene.dot(gi, gj);
			}
		    }
		}
	    }
	    for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++)
		if (majorityClusterFor(i,k)){
		    gi = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(i, before);
		    simK -= Gene.dot(gi, gi);
		}
	    if (sizz != 0){
		neff ++;
		simK /= ((double) sizz);
		val += simK;
	    }
	}
	return val;
    }

    public int majorityCluster(int ex){
	int j, index = 0;
	for (j=1;j<nclusters;j++)
	    if (soft_memberships.coordinates[j][ex] > soft_memberships.coordinates[index][ex])
		index = j;
	return index;
    }
	 

    public void fillAverageDistance(){
	hard_average_distance = new double [nclusters];
	hard_min_distance = new double [nclusters];
	hard_max_distance = new double [nclusters];
	hard_median_distance = new double [nclusters];
	int i, j, k, nb, nb2;
	int [] indexcluster;
	double [] median;
	double dist, totdist, mindist = -1.0, maxdist = -1.0;
	for (j=0;j<nclusters;j++){
	    nb = nb2 = 0;
	    totdist = 0.0;
	    mindist = maxdist = -1.0;
	    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++)
		if (selectedGeneToHardCluster[i] == j)
		    nb++;
	    if (nb == 0)
		Matrix.perror("No point in cluster " + j);
	    indexcluster = new int [nb];
	    k = 0;
	    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++)
		if (selectedGeneToHardCluster[i] == j){
		    indexcluster[k] = i;
		    k++;
		}
	    if (nb > 1){
		median = new double [nb*(nb-1)/2];
		for (i=0;i<indexcluster.length-1;i++)
		    for (k=i+1;k<indexcluster.length;k++){
			dist = myAGCT.W.distance( ((Vector) ((Vector) myAGCT.allGenesCoordinates.elementAt(myAGCT.geneToCoord[indexcluster[i]])).elementAt(1)), ((Vector) ((Vector) myAGCT.allGenesCoordinates.elementAt(myAGCT.geneToCoord[indexcluster[k]])).elementAt(1)));
			totdist += dist;
			
			if ( (nb2 == 0) || (dist < mindist) )
			    mindist = dist;

			if ( (nb2 == 0) || (dist > maxdist) )
			    maxdist = dist;

			median[nb2] = dist;
			nb2 ++;
		    }
		totdist /= ((double) nb2);
		hard_average_distance[j] = totdist;
		hard_min_distance[j] = mindist;
		hard_max_distance[j] = maxdist;

		QuickSort.quicksort(median);
		hard_median_distance[j] = median[median.length / 2];
	    }else{
		hard_average_distance[j] = 0.0;
		hard_min_distance[j] = 0.0;
		hard_max_distance[j] = 0.0;
		hard_median_distance[j] = 0.0;
	    }
	}
	averageDistanceComputed = true;
    }

   
    public void toSelectedGeneToHardCluster(){
	selectedGeneToHardCluster = new int[myClustering.getRightNumberSelectedGenes(before)];
	int i, j, nnc, adden, gid, nid;
	Vector vid;
	Gene gg, ggn;

	if (AGCT.Referenced_Available){
	    if (AGCT.Max_Reference == -1)
		Matrix.perror("AGCTClustering_Algorithm.class :: bad max reference");
	    referencesByCluster = new double [AGCT.Max_Reference+1][];
	    for (i=0;i<AGCT.Max_Reference+1;i++){
		referencesByCluster[i] = new double [nclusters];
		for (j=0;j<nclusters;j++)
		    referencesByCluster[i][j] = 0.0;
	    }
	}

	for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++)
	    selectedGeneToHardCluster[i] = majorityCluster(i);

	if ( (AGCT.Filter_Similarities_With_Distances) || (AGCT.Filter_Triangulation_With_Distances) )
	    fillAverageDistance();

	clusterSizesPrototypes = new int[nclusters];
	clusterSizesWhole = new int[nclusters];
	for (i=0;i<nclusters;i++){
	    clusterSizesPrototypes[i] = 0;
	    clusterSizesWhole[i] = 0;
	}
	for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++){
	    nnc = selectedGeneToHardCluster[i];
	    adden = 0;
	    if ( (before) || (Prototype.No_Reduction) )
		adden = 1;
	    else{
		gid = myAGCT.myDomain.selectedGeneNumberToGeneNumber[i];
		vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
		adden = vid.size();
	    }

	    if (!before){
		gid = myAGCT.myDomain.selectedGeneNumberToGeneNumber[i];
		gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(gid);

		if (gg.referenced)
		    referencesByCluster[gg.typeReferenced][nnc] += 1.0;

		if (!Prototype.No_Reduction){
		    vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
		    for (j=1;j<vid.size();j++){
			nid = ( (Integer) vid.elementAt(j) ).intValue();
			ggn =  (Gene) myAGCT.myDomain.domainGenes.elementAt(nid);
			
			if (ggn.referenced)
			    referencesByCluster[ggn.typeReferenced][nnc] += 1.0;
		    }
		}
	    }

	    clusterSizesPrototypes[nnc]+=1;
	    clusterSizesWhole[nnc]+=adden;
	}

	myClustering.nclusters = nclusters;
    }

    public boolean majorityClusterFor(int ex, int cl){
	//ex MUST be between 0 and numberSelectedGenes

	int j;
	for (j=0;j<nclusters;j++)
	    if ( (j!=cl) && (soft_memberships.coordinates[j][ex] > soft_memberships.coordinates[cl][ex]) ){
		return false;
	    }
	return true;
    }

    public boolean majorityGeneCluster(int ng, int nc){
	if ( (soft_clustering_computed) || (hard_clustering_computed) ){
	    int i;
	    for (i=0;i<nclusters;i++){
		if ( (i!=nc) && (soft_memberships.coordinates[i][ng] > soft_memberships.coordinates[nc][ng]) )
		    //System.out.println("Ng = " + ng + ", Nc = " + nc + ", i = " + i);
		    return false;
	    }
	}
	return true;
    }


    public void toHardMemberships(){
	hard_memberships_centers = new Vector();
	int i, j, k, dd = getDimension(before);
	double tot;
	Pnt p;
	for (i=0;i<nclusters;i++){
	    p = new Pnt(dd);
	    for (k=0;k<dd;k++)
		p.coordinates[k] = 0.0;
 
	    tot = 0.0;
	    for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++){
		tot += soft_memberships.coordinates[i][j];
		for (k=0;k<dd;k++)
		    p.coordinates[k] += (soft_memberships.coordinates[i][j] * myClustering.getRightPoint(j, before).coordinates[k]);
	    }

	    if (tot == 0.0)
		Matrix.perror("No members for cluster " + i);

	    for (k=0;k<dd;k++)
		p.coordinates[k] /= tot;
	    hard_memberships_centers.add(p);
	}

	finalClusteringObjectiveFunction = totalDistortion();

	//System.out.println("Cluster centers : " + hard_memberships_centers + "; objective function = " + finalClusteringObjectiveFunction);

    }

    //Bregman Arthur Vassilvitskii init.

    public void update_for_chosen(int ch, boolean [] chosen, double [] min_distortion, int dist_used, double parametre){
	Pnt p, c = myClustering.getRightPoint(ch, before);
	int j, dim = myClustering.getRightNumberSelectedGenes(before);
	double dist;
	min_distortion[ch] = 0.0;

	for (j=0;j < dim;j++){
	    if (!chosen[j]){
		p = myClustering.getRightPoint(j, before);
		
		//System.out.println(dist_used);

		if (dist_used == -1)
		    dist = Distortion.distortion_l22(p, c);
		else
		    dist = Distortion.Bregman(p, c, dist_used, parametre);
		
		if ( (dist < 0.0) && (dist > -Precision_For_Eigensystems) )
		    dist = 0.0;
		
		if (dist < 0.0)
		    Matrix.perror("Negative distortion (" + dist + ") for j = " + j + " in AV seeding");
		
		if ( (min_distortion[j] == -1.0) || ( dist < min_distortion[j] ) )
		    min_distortion[j] = dist;
	    }
	}
    }

    public void initHardMembershipArthurVassilvitskii(int dist_used, double parametre){
	int i, j, k, dim = myClustering.getRightNumberSelectedGenes(before);
	boolean deja;
	int [] indexes_to_pick = new int [nclusters];
	double [] min_distortion = new double [dim];
	boolean [] chosen = new boolean [dim];
	Gene c, p;
	double dist, totdist, ra, binf, bsup;

	for (i=0;i<nclusters;i++)
	    indexes_to_pick[i] = -1;

	for (i=0;i<dim;i++){
	    chosen[i] = false;
	    min_distortion[i] = -1.0;
	}

	indexes_to_pick[0] = AGCT.RG.nextInt(dim);

	AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "AV Initialization", nclusters - 1);
	for (i=1;i<nclusters;i++){
	    chosen[indexes_to_pick[i-1]] = true;
	    update_for_chosen(indexes_to_pick[i-1], chosen, min_distortion, dist_used, parametre);

	    totdist = 0.0;
	    for (j=0;j<dim;j++)
		if (!chosen[j])
		    totdist += min_distortion[j];
	    if (totdist == 0.0)
		    Matrix.perror("Zero total distortion in AV seeding");

	    ra = (AGCT.RG.nextDouble())*totdist;

	    if (ra == 0.0)
		Matrix.perror("No more seed to select");

	    binf = 0.0;
	    j = 0;
	    while( (j < dim) && ( (chosen[j]) || (ra<binf) || (ra>binf + min_distortion[j]) ) ){
		if (!chosen[j])
		    binf += min_distortion[j];
		j++;
	    }
	    indexes_to_pick[i] = j;
	    cc.increment();
	}
	cc.end();

	for (i=0;i<nclusters-1;i++)
	    for (j=i+1;j<nclusters;j++)
		if (indexes_to_pick[i] == indexes_to_pick[j])
		    Matrix.perror("AV seeded twice the same point");

	for (i=0;i<nclusters;i++)
	    hard_memberships_centers.add(new Pnt(myClustering.getRightPoint(indexes_to_pick[i], before))); //Correction !
    }

}
