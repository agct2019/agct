import java.util.*;
import java.awt.Color;

class Cluster{
    Vector allGenes;
    //ints à mapper par selectedGeneNumberToGeneNumber

    Pnt center;
    Pnt manifold_center;
    Pnt pca_center;

    Pnt3D manifold_Pnt3D;
    Pnt3D pca_Pnt3D;

    int index;
    int indexPrevious1, indexPrevious2;
    //points to the previous clusters that were merged;

    int clusterIndex;
    boolean top;

    double curdist;
    public Cluster(Clustering cc, int ng){
	Gene gg = cc.getRightGeneSelectedGeneNumberToGeneNumber(ng, cc.myClusteringAlgorithm.before);
	allGenes = new Vector();
	allGenes.add(new Integer(ng));
	center = new Pnt(cc.getRightPoint(ng, cc.myClusteringAlgorithm.before));
	manifold_center = new Pnt(gg.manifold_components_total);
	pca_center = new Pnt(gg.pca_components);
	indexPrevious1 = -1;
	indexPrevious2 = -1;
	clusterIndex = -1;
	index = ng;
	curdist = 0.0;

	manifold_Pnt3D = new Pnt3D();
	pca_Pnt3D = new Pnt3D();

	manifold_Pnt3D.coordinates[0] = manifold_center.coordinates[0];
	manifold_Pnt3D.coordinates[1] = manifold_center.coordinates[1];
	manifold_Pnt3D.coordinates[2] = manifold_center.coordinates[2];

	pca_Pnt3D.coordinates[0] = pca_center.coordinates[0];
	pca_Pnt3D.coordinates[1] = pca_center.coordinates[1];
	pca_Pnt3D.coordinates[2] = pca_center.coordinates[2];

	top = false;
    }

    public Cluster(Cluster cop){
	int i;

	if (cop.allGenes != null){
	    allGenes = new Vector();
	    for (i=0;i<cop.allGenes.size();i++)
		allGenes.add(new Integer( ( (Integer) cop.allGenes.elementAt(i) ).intValue() ));
	}

	if (cop.center != null)
	    center = new Pnt(cop.center);

	if (cop.manifold_center != null)
	    manifold_center = new Pnt(cop.manifold_center);

	if (cop.pca_center != null)
	    pca_center = new Pnt(cop.pca_center);

	manifold_Pnt3D = new Pnt3D(cop.manifold_Pnt3D.coordinates[0],
				   cop.manifold_Pnt3D.coordinates[1],
				   cop.manifold_Pnt3D.coordinates[2]);
	pca_Pnt3D = new Pnt3D(cop.pca_Pnt3D.coordinates[0],
			      cop.pca_Pnt3D.coordinates[1],
			      cop.pca_Pnt3D.coordinates[2]);

	indexPrevious1 = cop.indexPrevious1;
	indexPrevious2 = cop.indexPrevious2;

	index = cop.index;
	clusterIndex = cop.clusterIndex;
	curdist = cop.curdist;

	top = false;
    }

    public void putInside(Cluster cc, int p1, int p2, int np, double dcur){
	int i;
	double ms = (double) allGenes.size(), hs = (double) cc.allGenes.size();
	indexPrevious1 = p1;
	indexPrevious2 = p2;
	index = np;

	for (i=0;i<cc.allGenes.size();i++)
	    allGenes.add((Integer) cc.allGenes.elementAt(i));

	if (center.dimension() != cc.center.dimension())
	    Matrix.perror("Dimension mismatch");
	if (manifold_center.dimension() != cc.manifold_center.dimension())
	    Matrix.perror("Dimension mismatch Manifold");
	if (pca_center.dimension() != cc.pca_center.dimension())
	    Matrix.perror("Dimension mismatch PCA");

	for (i=0;i<center.dimension();i++)
	    center.coordinates[i] = ( (ms*center.coordinates[i]) + (hs*cc.center.coordinates[i]) ) / (ms + hs);
	for (i=0;i<manifold_center.dimension();i++)
	    manifold_center.coordinates[i] = ( (ms*manifold_center.coordinates[i]) + (hs*cc.manifold_center.coordinates[i]) ) / (ms + hs);
	for (i=0;i<pca_center.dimension();i++)
	    pca_center.coordinates[i] = ( (ms*pca_center.coordinates[i]) + (hs*cc.pca_center.coordinates[i]) ) / (ms + hs);
	
	manifold_Pnt3D.coordinates[0] = manifold_center.coordinates[0];
	manifold_Pnt3D.coordinates[1] = manifold_center.coordinates[1];
	manifold_Pnt3D.coordinates[2] = manifold_center.coordinates[2];

	pca_Pnt3D.coordinates[0] = pca_center.coordinates[0];
	pca_Pnt3D.coordinates[1] = pca_center.coordinates[1];
	pca_Pnt3D.coordinates[2] = pca_center.coordinates[2];

	curdist = dcur;
    }

    public void printToString(){
	System.out.println("Cluster " + index);
	System.out.println("Center = " + center);
	System.out.println("Manifold Center = " + manifold_center);
	System.out.println("PCA Center = " + pca_center);
	System.out.println("All genes = " + allGenes);
	System.out.println("index Previous 1 = " + indexPrevious1 + ", index Previous 2 = " + indexPrevious2);
    }

    public double distance_Ward(Cluster cc){
	double v;
	double ms = (double) allGenes.size(), hs = (double) cc.allGenes.size();
	v = (ms*hs/(ms+hs))*Distortion.distortion_l22(center, cc.center);
	return v;
    }

    public double distance_Single_Linkage(Clustering cl, Cluster cc){
	double dmin = -1, dcur;
	int i, j;
	Pnt mp, hp;
	for (i=0;i<allGenes.size();i++)
	    for (j=0;j<cc.allGenes.size();j++){
		mp = cl.getRightPoint( ( (Integer) allGenes.elementAt(i) ).intValue() , cl.myClusteringAlgorithm.before);
		hp = cl.getRightPoint( ( (Integer) cc.allGenes.elementAt(j) ).intValue() , cl.myClusteringAlgorithm.before);
		dcur = Distortion.distortion_l22(mp, hp);
		if ( ( (i==0) && (j==0) ) || (dcur < dmin) )
		    dmin = dcur;
	    }
	if (dmin == -1)
	    Matrix.perror("No distance computed !");
	return dmin;
    }
}

 class AGCTClustering_HC extends AGCTClustering_Algorithm{

    Vector nonDummyClusters;
    Vector nonDummyHierarchy;

    Vector dummyClusters;
    Vector dummyHierarchy;

    int dummyIndex;
    int nonDummyIndex;

    int distanceType;

    boolean treeAvailable;
    double totdistort;

    public static double getDistance(AGCTClustering_HC ref, Cluster c1, Cluster c2){
	double val = -1;

	if (ref.distanceType == 0)
	    val = c1.distance_Ward(c2);
	else if (ref.distanceType == 1)
	    val = c1.distance_Single_Linkage(ref.myClustering, c2);

	if (val == -1)
	    Matrix.perror("Distance type not valid!");

	return val;
    }


    AGCTClustering_HC(Clustering myc, AGCT mya, int nc, int dt, boolean bef){
	super(myc, mya, nc, bef, 0);
	seedingType = -1;
	distanceType = dt;
	if ( (distanceType < 0) || (distanceType > 1) )
	    Matrix.perror("Distance type not valid!");
	myReferenceName = "HC";
	dummyIndex = 0;
	nonDummyIndex = 0;

	treeAvailable = false;
    }

    public void plot(JAGCTVisualizationPane vv){
	plotClusterCenters(vv);
	plotTree(vv);
    }

    public void plotWithAnnotations(JAGCTVisualizationPane vv){
	plotClusterCentersWithAnnotations(vv);
	plotTree(vv);
    }

    public void plotTree(JAGCTVisualizationPane vv){
	Color c;
	int i, i1, i2;
	Pnt3D p, p1, p2;
	Cluster cc, c1, c2;
	View3D v = vv.visualizationPanel.myView3D;
	if ( treeAvailable ){
	    if (vv.myReferenceName.equals(JAGCTVisualizationPane.M_P)){
		treePlot_Manifold_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
		for (i=0;i<nonDummyHierarchy.size();i++){
			cc = (Cluster) nonDummyHierarchy.elementAt(i);
			if (vv.structPlottedCluster(cc.clusterIndex)){
			    i1 = cc.indexPrevious1;
			    i2 = cc.indexPrevious2;
			    if ( (i1 != -1) && (i2 != -1) ){
				c1 = (Cluster) nonDummyHierarchy.elementAt(i1);
				c2 = (Cluster) nonDummyHierarchy.elementAt(i2);
				p = cc.manifold_Pnt3D;
				p1 = c1.manifold_Pnt3D;
				p2 = c2.manifold_Pnt3D;
				
				v.drawShadowLineVsRef(p, p1, JAGCTGraphicsPanel.minColor, myClustering.colorCluster[cc.clusterIndex], false);
				v.drawShadowLineVsRef(p, p2, JAGCTGraphicsPanel.minColor, myClustering.colorCluster[cc.clusterIndex], false);
			    }
			}
		    }
	    }else if (vv.myReferenceName.equals(JAGCTVisualizationPane.P_P)){
		treePlot_Pca_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
		for (i=0;i<nonDummyHierarchy.size();i++){
			cc = (Cluster) nonDummyHierarchy.elementAt(i);
			if (vv.structPlottedCluster(cc.clusterIndex)){
			    i1 = cc.indexPrevious1;
			    i2 = cc.indexPrevious2;
			    if ( (i1 != -1) && (i2 != -1) ){
				c1 = (Cluster) nonDummyHierarchy.elementAt(i1);
				c2 = (Cluster) nonDummyHierarchy.elementAt(i2);
				p = cc.pca_Pnt3D;
				p1 = c1.pca_Pnt3D;
				p2 = c2.pca_Pnt3D;
				
				v.drawShadowLineVsRef(p, p1, JAGCTGraphicsPanel.minColor, myClustering.colorCluster[cc.clusterIndex], false);
				v.drawShadowLineVsRef(p, p2, JAGCTGraphicsPanel.minColor, myClustering.colorCluster[cc.clusterIndex], false);
			    }
			}
		}
	    }
	}
    }

    public void moveIndexes(){
	int i, i1, i2, in1, in2;
	Cluster cc, c1, c2;

	for (i=0;i<nonDummyHierarchy.size();i++){
	    cc = (Cluster) nonDummyHierarchy.elementAt(i);
	    i1 = cc.indexPrevious1;
	    i2 = cc.indexPrevious2;
	    if ( (i1 != -1) && (i2 != -1) ){
		c1 = (Cluster) nonDummyHierarchy.elementAt(i1);
		c2 = (Cluster) nonDummyHierarchy.elementAt(i2);
		in1 = c1.clusterIndex;
		in2 = c2.clusterIndex;
		
		if (in1 != in2)
		    Matrix.perror("Mismatch between clusters");
		else
		    cc.clusterIndex = in1;
	    }else
		if (cc.clusterIndex == -1)
		    Matrix.perror("Negative cluster index");
	}
    }

    public void treePlot_Manifold_Pnt3D(int xAxis, int yAxis, int zAxis){
	Pnt3D p;
	int i;
	for (i=0;i<nonDummyHierarchy.size();i++){
	    p = ( (Cluster) nonDummyHierarchy.elementAt(i) ).manifold_Pnt3D;
	    p.coordinates[0] = ( (Cluster) nonDummyHierarchy.elementAt(i) ).manifold_center.coordinates[xAxis];
	    p.coordinates[1] = ( (Cluster) nonDummyHierarchy.elementAt(i) ).manifold_center.coordinates[yAxis];
	    p.coordinates[2] = ( (Cluster) nonDummyHierarchy.elementAt(i) ).manifold_center.coordinates[zAxis];
	}
    }

    public void treePlot_Pca_Pnt3D(int xAxis, int yAxis, int zAxis){
	Pnt3D p;
	int i;
	for (i=0;i<nonDummyHierarchy.size();i++){
	    p = ( (Cluster) nonDummyHierarchy.elementAt(i) ).pca_Pnt3D;
	    p.coordinates[0] = ( (Cluster) nonDummyHierarchy.elementAt(i) ).pca_center.coordinates[xAxis];
	    p.coordinates[1] = ( (Cluster) nonDummyHierarchy.elementAt(i) ).pca_center.coordinates[yAxis];
	    p.coordinates[2] = ( (Cluster) nonDummyHierarchy.elementAt(i) ).pca_center.coordinates[zAxis];
	}
    }

    public String toString(){
	String val = "Hierarchical clustering with " + nclusters + " centers.\n";
	val += "Points are " + Clustering.referenceStringCenters[myClustering.indexReferenceStringCenters] + ".\n";
	val += "Distance is " + Clustering.HC_DISTANCE[distanceType] + ".";
	return val;
    }

    public void init_CAH(){
	nonDummyClusters = new Vector();
	nonDummyHierarchy = new Vector();
	dummyClusters = new Vector();
	dummyHierarchy = new Vector();
	Cluster cc;
	int i;
	for (i=0;i<myClustering.getRightNumberSelectedGenes(before);i++){
	    cc = new Cluster(myClustering, i);
	    nonDummyClusters.add(cc);
	    cc = new Cluster(myClustering, i);
	    dummyClusters.add(cc);
	    cc = new Cluster(myClustering, i);
	    nonDummyHierarchy.add(cc);
	    nonDummyIndex++;
	    cc = new Cluster(myClustering, i);
	    dummyHierarchy.add(cc);
	    dummyIndex++;
	}
    }

    public void one_step_CAH(boolean dummy){
	int i, j;
	Cluster dci, dcj;
	double ddist, ddopt = -1;
	int dri = -1, drj = -1;

	Vector curClusters;
	Vector curHierarchy;
	if (dummy){
	    curClusters = dummyClusters;
	    curHierarchy = dummyHierarchy;
	}else{
	    curClusters = nonDummyClusters;
	    curHierarchy = nonDummyHierarchy;
	}

	for (i=0;i<curClusters.size()-1;i++){
	    for (j=i+1;j<curClusters.size();j++){
		dci = (Cluster) curClusters.elementAt(i);
		dcj = (Cluster) curClusters.elementAt(j);
		ddist = AGCTClustering_HC.getDistance(this, dci, dcj);
		if ( ( (i==0) && (j==1) ) || (ddist < ddopt) ){
		    dri = i;
		    drj = j;
		    ddopt = ddist;
		}
	    }
	}

	if (ddopt == -1)
	    Matrix.perror("Negative optimal distance!");

	if ( (dri == -1) || (drj == -1) )
	    Matrix.perror("Negative optimal distance!");
	
	if (!dummy)
	    totdistort += ddopt;

	dci = (Cluster) curClusters.elementAt(dri);
	dcj = (Cluster) curClusters.elementAt(drj);

	if (!dummy){
	    dci.putInside(dcj, dci.index, dcj.index, nonDummyIndex, totdistort);
	    nonDummyHierarchy.add(new Cluster(dci));
	    nonDummyIndex++;
	}else{
	    dci.putInside(dcj, dci.index, dcj.index, dummyIndex, totdistort);
	    dummyHierarchy.add(new Cluster(dci));
	    dummyIndex++;
	}

	if (!dummy)
	    nonDummyClusters.removeElementAt(drj);
	else
	    dummyClusters.removeElementAt(drj);
    }

    public void trueTops(){
	int i;
	boolean [] plotted = new boolean[nclusters];
	Cluster cl;
	for (i=0;i<nclusters;i++)
	    plotted[i] = false;
	for (i=nonDummyHierarchy.size() - 1; i>=0;i--){
	    cl = ((Cluster) nonDummyHierarchy.elementAt(i));
	    if (!plotted[cl.clusterIndex]){
		cl.top = true;
		plotted[cl.clusterIndex] = true;
	    }else
		cl.top = false;
	}
	plotted = null;
    }

    public void CAH(){
	int i, j, k, id, ntot = nonDummyClusters.size();
	Cluster cc;
	totdistort = 0.0;
	int nbtc = ntot;

	AGCTCounter ccc = new AGCTCounter(myAGCT.myInformationFrame, "Clustering HC", ntot - nclusters);
	while(nbtc > 1){
	    one_step_CAH(true);
	    if (nbtc > nclusters)
		one_step_CAH(false);

	    nbtc--;
	    ccc.increment();
	}
	ccc.end();

	/*System.out.println("Dummy Hierarchy size = " + dummyHierarchy.size());
	for (i=0;i<dummyHierarchy.size();i++)
	((Cluster) dummyHierarchy.elementAt(i)).printToString();*/

	for (i=nonDummyClusters.size()-1;i>=0;i--){
	    cc = (Cluster) nonDummyClusters.elementAt(i);
	    for (j=0;j<cc.allGenes.size();j++){
		id = ( (Integer) cc.allGenes.elementAt(j) ).intValue();
		for (k=0;k<nclusters;k++){
		    if (k==i)
			soft_memberships.coordinates[k][id] = 1.0;
		    else
			soft_memberships.coordinates[k][id] = 0.0;
		}
		( (Cluster) nonDummyHierarchy.elementAt(id) ).clusterIndex = i;
	    }
	}

	moveIndexes();
	trueTops();

	plotCentersAvailable = true;
	treeAvailable = true;
	soft_memberships_available = true;
    }

    public void toClustering(){
	init_CAH();
	CAH();
	toClusterChoices();
	hard_clustering_computed = true;
	soft_clustering_computed = false;
	toHardMemberships();

	toSelectedGeneToHardCluster();

	ok_for_statistics = true;
	ControlProcess.put("hardClusteringProcessed",true);
	ControlProcess.put("ClusteringProcessed_HC",true);
    }
}
