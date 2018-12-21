import java.util.*;

 class AGCTClustering_EM extends AGCTClustering_Algorithm{

    int dimRef;

    Matrix [] sig;
    double [] frac, lndets;
    double loglike;
    boolean stop_choleski;

    AGCTClustering_EM(Clustering myc, AGCT mya, int nc, int ti, boolean bef, int nMI){
	super(myc, mya, nc, bef, nMI);
	seedingType = ti;
	myReferenceName = "EM";
	finalClusteringObjectiveFunction = -1;

	sig = null;
	lndets = null;
	frac = null;
	stop_choleski = false;
    }

    public void plot(JAGCTVisualizationPane vv){
	plotClusterCenters(vv);
    }

    public void plotWithAnnotations(JAGCTVisualizationPane vv){
	plotClusterCentersWithAnnotations(vv);
    }

    public String toString(){
	String val = "Expectation Maximization with " + nclusters + " Gaussians.\n";
	val += "Initialization : " + AGCTClustering_Algorithm.SEEDING_NAME[seedingType] + ".\n";
	//val += "Distortion : " + Distortion.DISTORTION_NAME[distortionType] + ".\n";
	val += "Points are " + Clustering.referenceStringCenters[myClustering.indexReferenceStringCenters] + ".";
	return val;
    }

    public void initSoftMembershipsEM(){
	if (seedingType == 0)
	    initHardMembershipForgy();
	else if (seedingType == 1)
	    initHardMembershipArthurVassilvitskii(-1,-1.0);
	else
	    Matrix.perror("Invalid seeding type");

	int i, j, k;
	int kk = nclusters, mm = dimRef;

	sig = new Matrix[kk];
	for (k=0;k<kk;k++)
	    sig[k] = new Matrix("Covariance of cluster " + k, mm, mm);
	lndets = new double[kk];
	frac = new double[kk];

	for (k=0;k<kk;k++){
	    frac[k] = 1.0/(double) kk;
	    for (i=0;i<mm;i++){
		for (j=0;j<mm;j++)
		    sig[k].coordinates[i][j] = 0.0;
		sig[k].coordinates[i][i] = 1.0E-10;
	    }
	}
    }

    public double EStep(){
	int k, m, n, i;
	int kk = nclusters, mm = dimRef, nn = myClustering.getRightNumberSelectedGenes(before);
	double tmp, sum, max, oldloglike = loglike, valret;
	double [] u = new double[mm], v = new double[mm];
	Gene gg;
	boolean pass_choleski;

	for (k=0;k<kk;k++){
	    pass_choleski = Matrix.test_choltmp(sig[k]);
	    if (!pass_choleski){
		stop_choleski = true;
		break;
	    }
	}

	if (!stop_choleski){
	    for (k=0;k<kk;k++){
		Matrix.choltmp(sig[k]);
		lndets[k] = sig[k].logdet();
		for (n=0;n<nn;n++){
		    for (m=0;m<mm;m++)
			u[m] = myClustering.getRightPoint(n, before).coordinates[m] - ((Pnt) hard_memberships_centers.elementAt(k)).coordinates[m];
		    sig[k].elsolve(u, v);
		    for (sum = 0.0, m = 0; m < mm; m ++)
			sum += (v[m] * v[m]);
		    soft_memberships.coordinates[k][n] = -0.5 * (sum + lndets[k]) + Math.log(frac[k]);
		}
	    }
	    loglike = 0.0;
	    for (n=0;n<nn;n++){
		max = -99.9E99;
		for (k=0;k<kk;k++)
		    if(soft_memberships.coordinates[k][n] > max)
			max = soft_memberships.coordinates[k][n];
		for (sum=0.0, k=0;k<kk;k++)
		    sum+=Math.exp(soft_memberships.coordinates[k][n]-max);
		tmp = max + Math.log(sum);
		for (k=0;k<kk;k++)
		    soft_memberships.coordinates[k][n] = Math.exp(soft_memberships.coordinates[k][n] - tmp);
		loglike += tmp;
	    }
	    finalClusteringObjectiveFunction = loglike;
	    
	    valret = loglike - oldloglike;
	}else
	    valret = 0.0;

	return valret;
    }

    public void MStep(){
	int kk = nclusters, mm = dimRef, nn = myClustering.getRightNumberSelectedGenes(before);
	int j, n, k, m;
	double wgt, sum;
	for (k=0;k<kk;k++){
	    wgt = 0.0;
	    for (n=0;n<nn;n++)
		wgt += soft_memberships.coordinates[k][n];
	    frac[k] = wgt/ (double) nn;
	    for (m=0;m<mm;m++){
		for (sum=0.0, n=0;n<nn;n++)
		    sum += soft_memberships.coordinates[k][n]*myClustering.getRightPoint(n, before).coordinates[m];
		((Pnt) hard_memberships_centers.elementAt(k)).coordinates[m] = sum / wgt;
		for (j=0;j<mm;j++){
		    for (sum=0.0, n=0;n<nn;n++){
			sum += soft_memberships.coordinates[k][n]*
			    (myClustering.getRightPoint(n, before).coordinates[m] - ((Pnt) hard_memberships_centers.elementAt(k)).coordinates[m]) 
			    * (myClustering.getRightPoint(n, before).coordinates[j] - ((Pnt) hard_memberships_centers.elementAt(k)).coordinates[j]);
		    }
		    sig[k].coordinates[m][j] = sum/wgt;
		}
	    }
	}
    }

    public void EM(){
	int i, niter = 0;
	double delta;
	boolean stop = false;
	AGCTCounter ccc;

	if (nMaxIterations>0)
	    ccc = new AGCTCounter(myAGCT.myInformationFrame, "Clustering EM", nMaxIterations);
	else
	    ccc = new AGCTCounter(myAGCT.myInformationFrame, "Clustering EM", 0);

	do{
	    delta = Math.abs(EStep());

	    if (nMaxIterations > 0)
		ccc.increment();
	    else
		ccc.setPercent( (int) ( 100.0 * AGCTClustering_Algorithm.EPSILON_EM / delta) );

	    if ( (stop_choleski) || (delta < AGCTClustering_Algorithm.EPSILON_EM) )
		stop = true;
	    else
		MStep();

	    niter++;

	    if (nMaxIterations>0)
		if ( niter >= nMaxIterations )
		    stop = true;

	}while (!stop);
	ccc.end();

	plotCentersAvailable = true;
    }

    public void toClustering(){
	dimRef = getDimension(before);

	initSoftMembershipsEM();
	EM();
	toClusterChoices();
	soft_memberships_available = true;
	soft_clustering_computed = true;
	hard_clustering_computed = false;
	toHardMemberships();

	toSelectedGeneToHardCluster();

	ok_for_statistics = true;
	ControlProcess.put("softClusteringProcessed",true);
	ControlProcess.put("ClusteringProcessed_EM",true);
    }
}
