import java.util.*;

 class AGCTClustering_CP extends AGCTClustering_Algorithm{

    public static Matrix myCW;
    public static boolean cwComputed;
    
    public static void init(){
	myCW = null;
	cwComputed = false;
    }

    public static void toCW(Matrix C){
	myCW = new Matrix("Copy of" + C.name, C.dimX, C.dimY);
	cwComputed = true;
    }

    AGCTClustering_CP(Clustering myc, AGCT mya, int nc, boolean bef){
	super(myc, mya, nc, bef, 0);
	seedingType = -1;
	myReferenceName = "CP";
    }

    public void plot(JAGCTVisualizationPane vv){
	plotClusterCenters(vv);
    }

    public void plotWithAnnotations(JAGCTVisualizationPane vv){
	plotClusterCentersWithAnnotations(vv);
    }

    public String toString(){
	String val = "CP-factorization with " + nclusters + " clusters.";
	return val;
    }

    public void initSoftMembershipForgy(){
	int i, j, k, ib1, ib2, it, nvaria = myAGCT.dimFeatures;
	double dbcur=0.0, db1=0.0, db2=0.0, dt, vd1, vd2;
	boolean deja = false;
	Gene gg1, gg2, ggj, ggC;
	Random r = new Random();

	int [] indexes_to_pick = new int [nclusters];
	for (i=0;i<nclusters;i++)
	    indexes_to_pick[i] = 0;

	for (i=0;i<nclusters;i++){
	    do{
		j = r.nextInt(myClustering.getRightNumberSelectedGenes(before));
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
	    for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++)
		soft_memberships.coordinates[i][j] = 0.0;
  
	for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++){
	    ib1 = 0;
	    ib2 = 1;

	    gg1 = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(indexes_to_pick[ib1], before);
	    gg2 = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(indexes_to_pick[ib2], before);
	    ggj = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(j, before);

	    db1 = Gene.squaredMagnitudeDifference(gg1, ggj);
	    db2 = Gene.squaredMagnitudeDifference(gg2, ggj);

	    if (db1 > db2){
		it = ib1;
		ib1 = ib2;
		ib2 = it;
		
		dt = db1;
		db1 = db2;
		db2 = dt;
	    }

	    if (nclusters > 2){
		for (i=2;i<nclusters;i++){
		    ggC = myClustering.getRightGeneSelectedGeneNumberToGeneNumber(indexes_to_pick[i], before);
		    dbcur = db2 = Gene.squaredMagnitudeDifference(ggC, ggj);

		    if (dbcur < db1){
			ib2 = ib1;
			db2 = db1;
			
			ib1 = i;
			db1 = dbcur;
		    }else if (dbcur < db2){
			ib2 = i;
			db2 = dbcur;
		    }
		}

		soft_memberships.coordinates[ib1][j] = 1.0 / Math.sqrt(2.0);
		soft_memberships.coordinates[ib2][j] = 1.0 / Math.sqrt(2.0);
	    }else{
		i = 1 + r.nextInt(9);
		vd1 = ((double) i) / 10.0;
		vd2 = 1.0 - vd1;

		//System.out.println("vd1 = " + vd1 + ", vd2 = " + vd2);

		if ( (vd1 <= 0.0) || (vd2 <= 0.0) )
		    Matrix.perror("Not all strictly positive");

		soft_memberships.coordinates[0][j] = Math.sqrt(vd1);
		soft_memberships.coordinates[1][j] = Math.sqrt(vd2);
	    }
	}
    }

    public void clustering_Zass_Shashua(Matrix DS){
	if (!AGCTClustering_CP.cwComputed){
	    AGCTClustering_CP.toCW(DS);
	    AGCTClustering_CP.myCW.fastDoublyStochasticApproximation(myAGCT.myDomain);

	    Matrix.checkDoublyStochastic(AGCTClustering_CP.myCW);
	    Matrix.checkSymmetric(AGCTClustering_CP.myCW);
	}

	soft_memberships.completePositiveFactorization(myAGCT.myDomain, AGCTClustering_CP.myCW);

	soft_memberships.normalize();
	Matrix.checkColumnStochastic(soft_memberships);
	reduce();
	Matrix.checkColumnStochastic(soft_memberships);
	soft_memberships_available = true;
	plotCentersAvailable = true;
    }

    public void reduce(){
	boolean notEmpty [] = new boolean[nclusters];
	int i, j, fc = 0;
	double tot;

	for (i=0;i<nclusters;i++){
	    tot = 0.0;
	    for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++){
		tot += soft_memberships.coordinates[i][j];
	    }
	    if (tot == 0.0){
		notEmpty[i] = false;
		if (AGCT.Debug)
		    System.out.println("No members for cluster " + i);
	    }else{
		notEmpty[i] = true;
		fc++;
	    }
	}

	if (fc < nclusters){
	    Matrix nsm = new Matrix(soft_memberships.name, fc, soft_memberships.dimY);

	    fc = 0;
	    for (i=0;i<nclusters;i++)
		if (notEmpty[i] == true){
		    for (j=0;j<myClustering.getRightNumberSelectedGenes(before);j++)
			nsm.coordinates[fc][j] = soft_memberships.coordinates[i][j];
		    fc++;
		}

	    soft_memberships = nsm;
	    nsm = null;
	    nclusters = soft_memberships.dimX;
	    myClustering.nclusters = nclusters;
	}

	notEmpty = null;
    }

    public void toClustering(){
	initSoftMembershipForgy();
	clustering_Zass_Shashua(myAGCT.CW);
	toClusterChoices();
	soft_clustering_computed = true;
	hard_clustering_computed = false;

	toHardMemberships();

	toSelectedGeneToHardCluster();

	ok_for_statistics = true;

	ControlProcess.put("softClusteringProcessed",true);
	ControlProcess.put("ClusteringProcessed_CP",true);
    }
}
