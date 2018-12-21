import java.util.*;

class Gene implements Debuggable {

    public static Hashtable FinalIndex;
    // returns an Integer [][] with Integers X et Y
    // Keys are Integers = indexes

    public static boolean FinalIndex_computed = false;

    AGCT myAGCT;
    Domain myDomain;

    String name;
    String lightAnnotation;

    String asciiName;
    Vector annotationsF, annotationsP, annotationsC;
    boolean annotationsFilled;

    boolean referenced;
    int typeReferenced;

    double [][] rawCoordinates;
    //coordinates computed from the Domain features
    //X is #ligand, Y is time stamp
    boolean [][] undeterminedCoordinates;
    //yes iff raw coordinate is "?" for (X, Y)

    double [][] finalCoordinates;
    //coordinates used after processing with AGCT.Method_F

    Pnt finalCoordinates_1D;

    Vector cluster_memberships;
    //each element = Double[] with the cluster memberships

    /**********************************************************************************
     * Components before processing
     *****/

    int [] index_neighbors_dimension;
    double [] distance_neighbors_dimension;
    //index* = index among selected Genes of neighbors of nearest neighbors according to squaredMagnitudeDifference
    //distance* = distance from this to the gene in index*

    double [] local_dimension;
    //computes the local dimension of manifold for each value of k = 2, 3, ..., AGCT.Max_K_Dimension_Estimation
    //even when k=1 is not performed, the dimension is the same as for the other two vectors (faster processing, reduced risk of errors)

    double averageDimension;

    /**********************************************************************************
     * Components after processing
     *****/

    Vector neighbors;
    //Contains Integers = index among selected Genes of neighbors in the triangulation
    Vector neighbor_angles;
    //Contains the angle wrt neighbors

    Pnt manifold_components_total;
    //components on the manifold for the Gene
    Pnt manifold_components_triangulation;
    //components used for the triangulation
    Pnt manifold_point3D;
    //idem (in 4D), but with a Pnt3D, should be deprecated (no more scaling)
    Pnt3D manifold_Pnt3D;

    Pnt pca_components;
    //principal components for the Gene
    Pnt pca_point3D;
    //idem, but with a Pnt3D
    Pnt3D pca_Pnt3D;

    boolean selected, manifold_components_computed, pca_components_computed;

    public static void test_computeFinalIndex(double [][] finalC){
	if (!FinalIndex_computed)
	    computeFinalIndex(finalC);
    }

    public static void flushAll(){
	FinalIndex = null;
	FinalIndex_computed = false;
    }

    public static void computeFinalIndex(double [][] finalC){
	//on the basis of finalC, computes the Hashtables

	Integer [] I;
	int rx, ry, i = 0;
	boolean arret = false;

	FinalIndex_computed = false;
	FinalIndex = null;
	FinalIndex = new Hashtable();

	rx = ry = 0;
	do{
	    I = new Integer[2];
	    I[0] = new Integer(rx);
	    I[1] = new Integer(ry);
	    
	    FinalIndex.put(new Integer(i), I);

	    ry++;
	    if (ry == finalC[rx].length){
		rx++;
		ry = 0;

		if (rx == finalC.length)
		    arret = true;
	    }
	    i++;

	}while(!arret);

	FinalIndex_computed = true;
    }


    public static void printFinalIndex(){
	Enumeration extensions = FinalIndex.keys();
	Integer I;
	if(extensions != null) {
	    while (extensions.hasMoreElements()) {
		I = (Integer) extensions.nextElement();
		System.out.print(I + " -> " + ((Integer[]) FinalIndex.get(I))[0] + ", " + ((Integer[]) FinalIndex.get(I))[1] + "; ");
	    }
	}
    }

    public void computeAverageDimension(){
	int k;
	averageDimension = 0.0;
	for (k=AGCTClustering_Algorithm.DELTA_DIMENSION;k<AGCT.Max_K_Dimension_Estimation - AGCTClustering_Algorithm.DELTA_DIMENSION;k++)
		averageDimension += (1.0 / local_dimension[k]);
	averageDimension = 1.0 / averageDimension;
	averageDimension *= ( (double) (AGCT.Max_K_Dimension_Estimation - 2 * AGCTClustering_Algorithm.DELTA_DIMENSION) );
    }

    public void fillNeighborsDimension(double [] distances, int [] indexes){
	int i, j, k, l;

	if ( (finalCoordinates_1D == null) || (finalCoordinates_1D.coordinates == null) || (finalCoordinates_1D.coordinates.length == 0) )
	    Matrix.perror("finalCoordinates_1D of gene " + name + " / " + asciiName + " is null or of length 0");

	double sum, vallog, minlog = 1.0 / ((double) finalCoordinates_1D.coordinates.length);
	boolean compte = false, allequal = true;
	index_neighbors_dimension = new int [AGCT.Max_K_Dimension_Estimation];
	distance_neighbors_dimension = new double [AGCT.Max_K_Dimension_Estimation];
	local_dimension = new double [AGCT.Max_K_Dimension_Estimation];

	if (distances.length < AGCT.Max_K_Dimension_Estimation)
	    Matrix.perror("Not enough genes in data to fill neighbors for dimension estimation (" + distances.length + " < " + AGCT.Max_K_Dimension_Estimation);

	if ( ( (name != null) && (name.toLowerCase().contains(new String("YGR010W").toLowerCase())) ) || ( (asciiName != null) && (asciiName.toLowerCase().contains(new String("YGR010W").toLowerCase())) ) ){
	    System.out.print(name + " / " + asciiName + " (dim = " + finalCoordinates_1D.coordinates.length + ") --- DIST : ");
	    for (i=0;i<AGCT.Max_K_Dimension_Estimation;i++)
		System.out.print(distances[i] + " ");
	    System.out.println("");
	}

	for (i=0;i<AGCT.Max_K_Dimension_Estimation;i++){
	    index_neighbors_dimension[i] = indexes[i];
	    distance_neighbors_dimension[i] = distances[i];

	    if ( (i == 0) || (distance_neighbors_dimension[i] == 0.0) )
		local_dimension[i] = -1.0;
	    else{
		j = 0;
		while(local_dimension[j] == 0.0){
		    j++;
		};
		if (j > i-3)
		    local_dimension[i] = -1.0;
		else{
		    k = 0; //k = i - j - 1;
		    sum = 0.0;
		    for (l = j; l < i ; l++){
			vallog = Math.log(distance_neighbors_dimension[i] / distance_neighbors_dimension[l]);

			//sum += Math.log(distance_neighbors_dimension[i] / distance_neighbors_dimension[l]);
			if (vallog > minlog){//(Math.abs(distance_neighbors_dimension[i] - distance_neighbors_dimension[l]) > Statistics.EPS_DIM){
			    sum += vallog;
			    k += 1;
			    allequal = false;
			    if ( ( (name != null) && (name.toLowerCase().contains(new String("YGR010W").toLowerCase())) ) || ( (asciiName != null) && (asciiName.toLowerCase().contains(new String("YGR010W").toLowerCase())) ) )
			      System.out.println(name + " / " + asciiName + " DELTA = " + Math.abs(distance_neighbors_dimension[i] - distance_neighbors_dimension[l]));
			}
		    }
		    if (!allequal){
			sum /= (double) k;
			local_dimension[i] = 1.0 / sum;
		    }else
			local_dimension[i] = 1.0;
		    compte = true;
		}
	    }
	}

	if ( ( (name != null) && (name.toLowerCase().contains(new String("YGR010W").toLowerCase())) ) || ( (asciiName != null) && (asciiName.toLowerCase().contains(new String("YGR010W").toLowerCase())) ) ){ 
	    System.out.print(name + " / " + asciiName + " DIM : ");
	    for (i=0;i<AGCT.Max_K_Dimension_Estimation;i++)
		System.out.print(local_dimension[i] + " ");
	    System.out.println("");
	}

	if (!compte){
	    for (i=0;i<AGCT.Max_K_Dimension_Estimation;i++)
		local_dimension[i] = AGCT.MIN_DIMENSION;
	    System.out.println(" ** Warning :: Gene " + name + " / " + asciiName + " has zero dimension estimation");
	}else{
	    while ( (AGCT.Min_K_Dimension_Estimation < AGCT.Max_K_Dimension_Estimation) && (local_dimension[AGCT.Min_K_Dimension_Estimation-1] == -1.0) ){
		AGCT.Min_K_Dimension_Estimation ++;
	    }

	    if (AGCT.K_Dimension_Estimation < AGCT.Min_K_Dimension_Estimation)
		AGCT.K_Dimension_Estimation = AGCT.Min_K_Dimension_Estimation;

	    if (AGCT.Min_K_Dimension_Estimation == AGCT.Max_K_Dimension_Estimation)
		Matrix.perror("AGCT.Min_K_Dimension_Estimation is = AGCT.Max_K_Dimension_Estimation for Gene " + name + " / " + asciiName);
	}
    }


    public String toLightHTMLString(){
	String val = "<HTML>", sc, sn;

	if (name != null){
	    val += name;
	    if ( (typeReferenced == -1) || (myAGCT.myDomain.orderedReferences == null) ){
		if (typeReferenced != -1)
		    val += "(referenced as " + typeReferenced + ")";
	    }else{
		sc = Integer.toHexString(JAGCTGraphicsPanel.referencedColors[typeReferenced].getRGB());
		sn = sc.substring(2, sc.length());
		val += "(referenced as <FONT COLOR=#" + sn + ">";
		val += myAGCT.myDomain.orderedReferences.elementAt(typeReferenced);
		val += "</FONT>)";
	    }
	}
	if ( (asciiName != null) && (asciiName != "") )
	    val += "_" + asciiName + "";
	val += "</HTML>";
	return val;
    }

    public String toString(){
	String val = "";
	int i, j;
	Gene nn;
	Double [] cc;

	if (name != null)
	    val += name + "(" + referenced + "|" + typeReferenced + ")";
	if ( (asciiName != null) && (asciiName != "") )
	    val += "_" + asciiName + "";
	val += "\n";

	if (annotationsF != null){
	    val += " -- Annotations (F):\n";
	    for (i=0;i<annotationsF.size();i++)
		val += " " + ((String[])annotationsF.elementAt(i))[0] + "  " + ((String[])annotationsF.elementAt(i))[1] + "\n";
	    val += "\n";
	}
	if (annotationsP != null){
	    val += " -- Annotations (P):\n";
	    for (i=0;i<annotationsP.size();i++)
		val += " " + ((String[])annotationsP.elementAt(i))[0] + "  " + ((String[])annotationsP.elementAt(i))[1] + "\n";
	    val += "\n";
	}
	if (annotationsC != null){
	    val += " -- Annotations (C):\n";
	    for (i=0;i<annotationsC.size();i++)
		val += " " + ((String[])annotationsC.elementAt(i))[0] + "  " + ((String[])annotationsC.elementAt(i))[1] + "\n";
	    val += "\n";
	}

	if (neighbors != null){
	    val += " -- Natural Neighbors (triangulation):\n";
	    for (i=0;i<neighbors.size();i++){
		nn = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[( (Integer) neighbors.elementAt(i)).intValue()]);
		val += " " + nn.name;
		if ( (nn.asciiName != null) && (nn.asciiName != "") )
		    val += "_" + nn.asciiName + "";
		val += "\n";
	    }
	}

	if (!AGCTFileWriter.LIGHT_SAVING){
	    if (cluster_memberships != null){
		val += " -- Cluster memberships:\n";
		for (i=0;i<cluster_memberships.size();i++){
		    cc = getClusterMemberships(i);
		    val += " - Memberships [Clustering #" + i + "]: ";
		    for (j=0;j<cc.length;j++)
			val += " Cluster_" + j + " : " + DF.format(getClusterMemberships(i,j)) + " --";
		    val += "\n";
		}
	    }
	}

	if ( (lightAnnotation != null) && (lightAnnotation != "") )
	    val += "\n -- light Annotation (in gene file): " + lightAnnotation + "\n";

	val += "\n\n";
 	
	return val;
    }
    
    public Double[] getClusterMemberships(int i){
	return (Double []) cluster_memberships.elementAt(i);
    }

    public double getClusterMemberships(int i, int j){
	return getClusterMemberships(i)[j].doubleValue();
    }
	
    public void addClusterMemberships(Double [] cm){
	if (cluster_memberships == null)
	    cluster_memberships = new Vector();
	cluster_memberships.add(cm);
    }

    public void putReferenced(){
	referenced = true;
    }

    public void typeReferenced(int i){
	typeReferenced = i;
    }

    public boolean geneOrNeighborsIsReferencedDelaunay(){
	boolean ref = false;
	int j;
	Gene nn;
	if (neighbors != null){
	    for (j=0;j<neighbors.size();j++){
		nn = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[( (Integer) neighbors.elementAt(j)).intValue()]);
		if (nn.referenced)
		    return true;
	    }
	}
	return referenced;
    }

    Gene(String name, int nl){
	cluster_memberships = null;
	asciiName = null;
	annotationsF = annotationsP = annotationsC = null;
	annotationsFilled = false;
	referenced = false;
	typeReferenced = -1;
	neighbors = null;

	this.name = name;
	finalCoordinates = null;
	lightAnnotation = "";
	selected = false;
	manifold_components_computed = false;

	rawCoordinates = new double [nl][];
	undeterminedCoordinates = new boolean [nl][];
	manifold_components_triangulation = new Pnt(AGCT.Number_Of_Triangulation_Dimensions);
	manifold_components_total = new Pnt(AGCT.Number_Of_Manifold_Components);

	manifold_point3D = new Pnt(4);
	pca_point3D = new Pnt(4);

	manifold_Pnt3D = new Pnt3D();
	pca_Pnt3D = new Pnt3D();
    }

    public static boolean dimCoincide(Gene a, Gene b){
	return (a.finalCoordinates.length == b.finalCoordinates.length);
    }

    public static double squaredMagnitudeDifference(Gene a, Gene b){
	//returns || a - b ||^2
	double v = 0.0, fca, fcb;
	int i, dim = Feature.getNumberOfFeaturesAfterSelection();
	if (!dimCoincide(a, b))
	    Matrix.perror("Dimension mismatch between genes " + a.name + " and " + b.name);

	for (i=0;i<dim;i++){
	    fca = a.getFinalCoordinates(i,true);
	    fcb = b.getFinalCoordinates(i,true);
	    v += ( (fca - fcb) * (fca - fcb) );
	}
	return v;
    }

    public static double l22(Gene a){
	// returns || a ||^2
	double v = 0.0;
	int i, dim = Feature.getNumberOfFeaturesAfterSelection();
	for (i=0;i<dim;i++)
	    v += ( a.getFinalCoordinates(i,true) * a.getFinalCoordinates(i,true) );

	return v;
    }

    public static double dot(Gene a, Gene b){
	// returns <a, b>
	double v = 0.0;
	int i, dim = Feature.getNumberOfFeaturesAfterSelection();
	if (!dimCoincide(a, b))
	    Matrix.perror("Dimension mismatch between genes " + a.name + " and " + b.name);

	for (i=0;i<dim;i++)
	    v += ( a.getFinalCoordinates(i,true) * b.getFinalCoordinates(i,true) );
	return v;
    }

   public static double getAbsoluteCosineSimilarity(Gene a, Gene b){
	if (!dimCoincide(a, b))
	    Matrix.perror("Dimension mismatch between genes " + a.name + " and " + b.name);

	double ps = dot(a, b);
	double dt1 = Math.sqrt(l22(a));
	double dt2 = Math.sqrt(l22(b));
	double val=0.0;
	double linf = 0.0;
	double lsup = 1.0;

	if ( (dt1 == 0.0) || (dt2 == 0.0) )
	    val = 0.0;
	else{
	    val =  Math.abs( ( (ps) / (dt1 * dt2) ) );

	    if ( (dt1 < 0.0) || (dt2 < 0.0) )
		Matrix.perror("Negative norms !");
	}

	if ( (val < linf) && (val > - Precision_For_Eigensystems) )
	    val = linf;
	if ( (val > lsup) && ( val < lsup + Precision_For_Eigensystems ) )
	    val = lsup;
	if ( (val < linf) || (val > lsup) )
	    Matrix.perror("Absolute cosine similarity outside bounds");
	return val;
    }

   public static double cosine(Gene a, Gene b){
	if (!dimCoincide(a, b))
	    Matrix.perror("Dimension mismatch between genes " + a.name + " and " + b.name);

	double ps = dot(a, b);
	double dt1 = Math.sqrt(l22(a));
	double dt2 = Math.sqrt(l22(b));
	double val=0.0;
	double linf = -1.0;
	double lsup = 1.0;

	if ( (dt1 == 0.0) || (dt2 == 0.0) )
	    val = 0.0;
	else{
	    val = ( (ps) / (dt1 * dt2) );

	    if ( (dt1 < 0.0) || (dt2 < 0.0) )
		Matrix.perror("Negative norms !");
	}

	if ( (val < linf) && (val > linf - Precision_For_Cosine) )
	    val = linf;
	if ( (val > lsup) && ( val < lsup + Precision_For_Cosine ) )
	    val = lsup;
	if ( (val < linf) || (val > lsup) ){
	    Matrix.perror("cosine " + val + " outside bounds");
	}return val;
    }

    public static double getCosineSimilarity(Gene a, Gene b){
	if (!dimCoincide(a, b))
	    Matrix.perror("Dimension mismatch between genes " + a.name + " and " + b.name);

	double ps = dot(a, b);
	double dt1 = Math.sqrt(l22(a));
	double dt2 = Math.sqrt(l22(b));
	double val=0.0;
	double linf = 0.0;
	double lsup = 1.0;

	if ( (dt1 == 0.0) || (dt2 == 0.0) )
	    val = 0.0;
	else{
	    val = (1.0 + ( (ps) / (dt1 * dt2) ) ) / 2.0;

	    if ( (dt1 < 0.0) || (dt2 < 0.0) )
		Matrix.perror("Negative norms !");
	}

	if ( (val < linf) && (val > - Precision_For_Eigensystems) )
	    val = 0.0;
	if ( (val > lsup) && ( val < lsup + Precision_For_Eigensystems ) )
	    val = lsup;
	if ( (val < linf) || (val > lsup) ){
	    Matrix.perror("Cosine similarity outside bounds");
	}
	return val;
    }

    public static double maxSim(){
	if (AGCT.Method_S == 0)
	    return 1.0;
	if (AGCT.Method_S == 1)
	    return 1.0;
	if (AGCT.Method_S == 2)
	    return 1.0;
	Matrix.perror("Bad value for Method_S");
	return -1.0;
    }

    public static double getSimilarity(Gene a, Gene b){
	if (a.name.equals(b.name))
	    return maxSim();

	double v = 0.0;
	if (AGCT.Method_S == 0)
	    v = Math.exp(-squaredMagnitudeDifference(a, b)/AGCT.T_Heat_Kernel);
	else if (AGCT.Method_S == 1)
	    v = getCosineSimilarity(a, b);
	else if (AGCT.Method_S == 2)
	    v = getAbsoluteCosineSimilarity(a, b);
	return v;
    }

    public String coordinatesString(){
	//we assume rawCoordinates, undeterminedCoordinates AND finalCoordinates computed
	//we also assume ligands selected ARE computed
	String val = " Gene " + name + " ---\n";
	int i, j;
	val += "Raw Coordinates by " +  AGCTFileWriter.DATA_Ligand + " : \n";
	for (i=0;i<rawCoordinates.length;i++){
	    val += " [ " + i + " : ";
	    for (j=0;j<rawCoordinates[i].length;j++){
		if (undeterminedCoordinates[i][j] == false)
		    val += Gene.DF.format(rawCoordinates[i][j]);
		else
		    val += "(?)";
		if (j < rawCoordinates[i].length - 1)
		    val += ", ";
	    }
	    val += "]\n";
	}
	val += "\n\nFinal Coordinates by " +  AGCTFileWriter.DATA_Ligand + " : \n";
	for (i=0;i<finalCoordinates.length;i++){
	    val += " [ " + i + " : ";
	    for (j=0;j<finalCoordinates[i].length;j++){
		val += Gene.DF.format(finalCoordinates[i][j]);
		if (j < finalCoordinates[i].length - 1)
		    val += ", ";
	    }
	    val += "]\n";
	}

	return val;
    }

    public double getFinalCoordinates(int x, boolean after){
	if (after) 
	    return getFinalCoordinatesAfterSelection(x);
	else
	    return getFinalCoordinatesBeforeSelection(x);
    }

    public void computeFinalCoordinates_1D(){
	int i, dim = Feature.getNumberOfFeaturesAfterSelection();
	finalCoordinates_1D = new Pnt(dim);
	for (i=0;i<dim;i++)
	    finalCoordinates_1D.coordinates[i] = getFinalCoordinates(i, true);
    }


    public double getFinalCoordinatesBeforeSelection(int x){
	//new, uses Hashtable FinalIndex
	Integer i = new Integer(x);
	if (!Gene.FinalIndex_computed)
	    Matrix.perror("FinalIndex not computed in Gene (perhaps a Thread problem ?)");

	return finalCoordinates[((Integer[]) Gene.FinalIndex.get(i))[0]][((Integer[]) Gene.FinalIndex.get(i))[1]];
    }

    public double getFinalCoordinatesAfterSelection(int x){
	//new, uses Hashtable FinalIndex
	Integer i = (Integer) Feature.featureIndexToNonTrashedIndex.get(new Integer(x));
	if (!Gene.FinalIndex_computed)
	    Matrix.perror("FinalIndex not computed in Gene (perhaps a Thread problem ?)");

	if (x > Feature.featureIndexToNonTrashedIndex.size())
	    Matrix.perror("feature coordinate outside bounds");

	return finalCoordinates[((Integer[]) Gene.FinalIndex.get(i))[0]][((Integer[]) Gene.FinalIndex.get(i))[1]];
    }

    public void finishVariables(AGCT a, Domain d){
	myAGCT = a;
	myDomain = d;
    }

    public void putLightAnnotation(String s){
	lightAnnotation = s;
    }

    public void addLightAnnotation(String s){
	lightAnnotation += s;
    }

    public void addAnnotation(boolean f, boolean p, boolean c, String s){
	String [] an;

	if ( (f && p) || (f && c) || (p && c) )
	    Matrix.perror("Gene.class :: at least two of F, P and C annotations are true for " + s);
	if ((!f) && (!p) && (!c))
	    Matrix.perror("Gene.class :: F, P and C annotations are false for " + s);
	else{
	    if ( (f) && (annotationsF == null) )
		annotationsF = new Vector();
	    else if ( (p) && (annotationsP == null) )
		annotationsP = new Vector();
	    else if ( (c) && (annotationsC == null) )
		annotationsC = new Vector();

	    an = new String[2];
	    if (f)
		an[0] = Domain.ANNOTATION_F;
	    else if (p)
		an[0] = Domain.ANNOTATION_P;
	    else
		an[0] = Domain.ANNOTATION_C;
		
	    an[1] = s;
	    if (f)
		annotationsF.add(an);
	    else if (p)
		annotationsP.add(an);
	    if (c)
		annotationsC.add(an);

	    annotationsFilled = true;
	}   
    }

    public void recomputeSelectedFeatures(){
	selected = true;

	double [][] rc2 = new double [myDomain.numberSelectedLigands][];
	boolean [][] uc2 = new boolean [myDomain.numberSelectedLigands][];

	int i, j, index = 0;
	for (i=0;i<myDomain.numberLigands;i++){
	    if (myDomain.selectedLigand[i]){
		rc2[index] = new double[( (Integer) myDomain.numberTimes.elementAt(i)).intValue()];
		uc2[index] = new boolean[( (Integer) myDomain.numberTimes.elementAt(i)).intValue()];
		for (j=0;j<( (Integer) myDomain.numberTimes.elementAt(i)).intValue();j++){
		    rc2[index][j] = rawCoordinates[i][j];
		    uc2[index][j] = undeterminedCoordinates[i][j];
		}
		index++;
	    }
	}

	rawCoordinates = rc2;
	undeterminedCoordinates = uc2;
	finalCoordinates = new double[myDomain.numberSelectedLigands][];
    }

    public void fillManifoldComponents(Matrix m, int n){
	int i, start = 0;

	if ( (myAGCT.Method_N == 0) || (myAGCT.Method_N == 1) )
	    start = 1;

	for (i=0;i<manifold_components_total.coordinates.length;i++){
	    if (Statistics.isNaN(m.coordinates[i+start][n]))
		Matrix.perror("Coordinate " + i + " of gene " + name + " / " + asciiName + " is NaN");

	    manifold_components_total.coordinates[i] = m.coordinates[i+start][n];
	    if (i<3){
		manifold_point3D.coordinates[i] = manifold_components_total.coordinates[i];
		manifold_Pnt3D.coordinates[i] = manifold_components_total.coordinates[i];
	    }
	    if (i<AGCT.Number_Of_Triangulation_Dimensions)
		manifold_components_triangulation.coordinates[i] = manifold_components_total.coordinates[i];
	}
	manifold_point3D.coordinates[3] = 1.0;
	manifold_components_computed = true;
    }

    public void fillPCAComponents(Matrix E, double [] average, double [] sigma){
    	double acp, slo, ave, sig, curc;
	int y, z;

	pca_components = new Pnt(myAGCT.dimFeatures);

	for (y=0;y<myAGCT.dimFeatures;y++)
	    pca_components.coordinates[y] = 0.0;;

	for (y=0;y<myAGCT.dimFeatures;y++){
	    for (z=0;z<myAGCT.dimFeatures;z++){
		acp = E.coordinates[y][z];
		slo = getFinalCoordinates(z,true);
		ave = average[z];
		sig = sigma[z];

		curc = pca_components.coordinates[y];
		if (sig > 0.0){
		    curc += ( acp * ( (slo - ave ) / sig ) );
		    pca_components.coordinates[y] = curc;
		}
	    }
	}

	pca_components_computed = true;

	pca_point3D.coordinates[0] = pca_components.coordinates[0];
	pca_point3D.coordinates[1] = pca_components.coordinates[1];
	pca_point3D.coordinates[2] = pca_components.coordinates[2];

	pca_point3D.coordinates[3] = 1.0;

	pca_Pnt3D.coordinates[0] = pca_components.coordinates[0];
	pca_Pnt3D.coordinates[1] = pca_components.coordinates[1];
	pca_Pnt3D.coordinates[2] = pca_components.coordinates[2];
    }

    public void updateManifold_Pnt3D(int xAxis, int yAxis, int zAxis){
	manifold_Pnt3D.coordinates[0] = manifold_components_total.coordinates[xAxis];
	manifold_Pnt3D.coordinates[1] = manifold_components_total.coordinates[yAxis];
	manifold_Pnt3D.coordinates[2] = manifold_components_total.coordinates[zAxis];
    }

    public void updatePCA_Pnt3D(int xAxis, int yAxis, int zAxis){
	pca_Pnt3D.coordinates[0] = pca_components.coordinates[xAxis];
	pca_Pnt3D.coordinates[1] = pca_components.coordinates[yAxis];
	pca_Pnt3D.coordinates[2] = pca_components.coordinates[zAxis];
    }


    public void initManifold3D(int xAxis, int yAxis, int zAxis){
	manifold_point3D.coordinates[0] = manifold_components_total.coordinates[xAxis];
	manifold_point3D.coordinates[1] = manifold_components_total.coordinates[yAxis];
	manifold_point3D.coordinates[2] = manifold_components_total.coordinates[zAxis];
	manifold_point3D.coordinates[3] = 1.0;
    }

    public void initPCA3D(int xAxis, int yAxis, int zAxis){
	pca_point3D.coordinates[0] = pca_components.coordinates[xAxis];
	pca_point3D.coordinates[1] = pca_components.coordinates[yAxis];
	pca_point3D.coordinates[2] = pca_components.coordinates[zAxis];
	pca_point3D.coordinates[3] = 1.0;
    }

    public void initManifold_Pnt3D(){
	//Makes the assumption that manifold_point3D has been updated !

	int j;
	double dum;

	for (j=0;j<3;j++){
	    dum = manifold_point3D.coordinates[j];
	    manifold_Pnt3D.coordinates[j] = (2.0 * (dum - myDomain.min_ManifoldPnt3D.coordinates[j]) / (myDomain.max_ManifoldPnt3D.coordinates[j] - myDomain.min_ManifoldPnt3D.coordinates[j]) ) - 1.0;
	}
    }

    public void putNeighbor(Integer u){
	//if (name.equals(new String("658")))
	//System.out.println("Neighbor : " + ((Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[u.intValue()])).name);
	double ac;

	if (neighbors == null){
	    neighbors = new Vector();
	    neighbor_angles = new Vector();
	}
	if (!neighbors.contains(u)){
	    neighbors.add(u);
	    ac = Math.acos(Gene.cosine(this,((Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[u.intValue()]))));
	    neighbor_angles.add(new Double(ac));
	    //System.out.println(ac + " --- " + Math.PI * Statistics.LIMIT_P/2.0 + ", " + Math.PI * (1.0 - (Statistics.LIMIT_P/2.0)));
	}

	if (neighbors.size() != neighbor_angles.size())
	    Matrix.perror("Size mismatch");
    }
}
