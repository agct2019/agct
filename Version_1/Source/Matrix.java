/*
 * Class Matrix
 */
import java.io.*;
import java.util.*;
import java.util.Map.Entry;

 class Matrix implements Debuggable{
    
    String name;

    double[][] coordinates;

    public double get(int x, int y) {
        return coordinates[x][y];
    }

    public void set(int x, int y, double value) {
       coordinates[x][y] = value;
    }

    public void times(double v, int row){
	int j;
	for (j=0;j<dimY;j++)
	    coordinates[row][j] *= v;
    }

    public static double FrobAbs(Matrix a, Matrix b, int maxRow){
	if (a.dimX != b.dimX)
	    Matrix.perror(a + " and " + b + " have dimX mismatch");
	if (a.dimY != b.dimY)
	    Matrix.perror(a + " and " + b + " have dimY mismatch");

	double v = 0.0;
	int i, j;
	for (i=0;i<maxRow;i++)
	    for (j=0;j<a.dimY;j++)
		v += ( (Math.abs(a.coordinates[i][j]) - Math.abs(b.coordinates[i][j])) * (Math.abs(a.coordinates[i][j]) - Math.abs(b.coordinates[i][j])) );
	return v;
    }

    public static double LinftyAbs(Matrix a, Matrix b, int maxRow){
	if (a.dimX != b.dimX)
	    Matrix.perror(a + " and " + b + " have dimX mismatch");
	if (a.dimY != b.dimY)
	    Matrix.perror(a + " and " + b + " have dimY mismatch");

	double v = 0.0;
	int i, j;
	for (i=0;i<maxRow;i++)
	    for (j=0;j<a.dimY;j++)
		if ( ( (i==0) && (j==0) ) || (v < Math.abs(Math.abs(a.coordinates[i][j]) - Math.abs(b.coordinates[i][j]))) )
		    v = Math.abs(Math.abs(a.coordinates[i][j]) - Math.abs(b.coordinates[i][j]));
	return v;
    }


    Matrix choleski_decomposition;

    int dimX;//number of rows
    int dimY;//number of columns

    boolean isVector; //true iff dimX = 1 (vectors are ROW)
    boolean hasDefaultContent; //true iff the content has defaultValue
    boolean isSquared; //true iff dimX = dimY
    boolean isSymmetric; //true iff isSquared && symmetric

    static double DefaultValue = 0;
    static double Sqrarg;
    static double Tiny = 1.0E-30;

    /****************************************************************************************************
     * Class Methods
     *****/

    public static boolean checkEigensystem(Domain dd, Matrix a, Matrix v, double [] eigval, int k){
	if (AGCT.Debug) System.out.print("Checking the eigensystem: " + v.name + " (" + v.dimX + " x " + v.dimY + ") contains the (row) eigenvectors of " + a.name + "... ");
	dd.myAGCT.myInformationFrame.setTextProgressBar("Check Eigs (" + v.name + ", " + a.name + ")");
	
	checkSameDimensions(a, v);
	checkSquare(a);
	
	double val, m;
	int x,y,xxmax, n = a.dimX, nmade = 0, percent, vm;
	xxmax = (k==-1) ? n : k;
	vm = xxmax;

	double [] fv = new double [n];
	double [] copy_v = new double [n];
	double [] delta = new double [n];
	
	for (x=0;x<xxmax;x++){
	    percent = (int) (100.0 * ((double) nmade / (double) (vm)));
	    dd.myAGCT.myInformationFrame.setValueProgressBar(percent);


	    for (y=0;y<n;y++)
		copy_v[y] = v.coordinates[x][y];

	    Matrix.dot(a, copy_v, fv);

	    for (y=0;y<n;y++){
		val = copy_v[y];
		copy_v[y] = val * eigval[x];
	    }
	    
	    Matrix.substract(copy_v, fv, delta);
	    val = Matrix.l22(delta);

	    if ( val > Precision_For_Eigensystems ){
		if (AGCT.Debug) System.out.print("\n Checking Eigensystem problem for vector " + x + " :: ||Av - lambda v|| = " + val + " (limit precision = " + Precision_For_Eigensystems + ")");
		dd.myAGCT.myInformationFrame.setValueProgressBar(0);
		dd.myAGCT.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);

		return false;
	    }
	    nmade ++;
	}
	
	if (AGCT.Debug) System.out.print("ok.\n");
	
	dd.myAGCT.myInformationFrame.setValueProgressBar(0);
	dd.myAGCT.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);

	return true;
    }
    
    public static double magnitude (double [] v){
	return Math.sqrt(dot(v, v));
    }
    
    public static double l22(double [] v){
	return dot(v, v);
    }

    public static int Hamming(int[] x, int [] y){
	if (x.length != y.length)
	    perror("Dimension mismatch");

	int val = 0;
	int i;
	for (i=0;i<x.length;i++)
	    if (x[i] != y[i])
		val ++;

	return val;
    }

    public static Matrix identity(int d){
	Matrix m = new Matrix("Identity", d, d);
	int i,j;

	for (i=0;i<d;i++)
	    for (j=0;j<d;j++)
		m.coordinates[i][j] = 0.0;
	for (i=0;i<d;i++)
	    m.coordinates[i][i] = 1.0;

	return m;
    }

    public static Matrix dot(Matrix v, Matrix w, String name){
	//makes the product inside new matrix
	int i, j, k;
	double sum;
	if (v.dimY != w.dimX)
	    Matrix.perror("Matrix.class :: Dimension mismatch");

	Matrix m = new Matrix(name, v.dimX, w.dimY);

	for (i=0;i<v.dimX;i++){
	    for (j=0;j<w.dimY;j++){
		sum = 0.0;
		for (k=0;k<v.dimY;k++)
		    sum += (v.coordinates[i][k] * w.coordinates[k][j]);
		m.coordinates[i][j] = sum;
	    }
	}

	return m;
    }

    public void dot(Matrix w, Matrix h){
	//makes the product inside this
	if (dimX != w.dimX)
	    Matrix.perror("Matrix.class :: Dimension mismatch (this,W)");
	if (dimY != h.dimY)
	    Matrix.perror("Matrix.class :: Dimension mismatch (this,H)");
	if (w.dimY != h.dimX)
	    Matrix.perror("Matrix.class :: Dimension mismatch (W,H)");
	int i, j, k;
	double sum;

	for (i=0;i<w.dimX;i++){
	    for (j=0;j<h.dimY;j++){
		sum = 0.0;
		for (k=0;k<w.dimY;k++)
		    sum += (w.coordinates[i][k] * h.coordinates[k][j]);
		coordinates[i][j] = sum;
	    }
	}

	isVector = w.isVector;
	hasDefaultContent = false;
	if (dimX == dimY)
	    isSquared = true;
	else
	    isSquared = false;
	isSymmetric = false;
    }

    public static double dot(double [] a, double [] b){
	double v = 0.0;
	int i;
	if (a.length != b.length)
	    Matrix.perror("Dot product between different dimensions");
	for (i=0;i<a.length;i++)
	    v += (a[i]*b[i]);
	return v;
    }

    public static void dot(Matrix a, double [] v, double [] res){
	int x, y;
	if (v.length != res.length)
	    Matrix.perror("Dot product between different dimensions");	double val;
	for (x=0;x<a.dimX;x++){
	    val = 0.0;
	    for (y=0;y<a.dimY;y++)
		val += (a.coordinates[x][y]*v[y]);
	    res[x] = val;
	}
    }

    public static void substract(double [] a, double [] b, double [] res){
	int i;
	if (a.length != b.length)
	    Matrix.perror("Difference between different dimensions");
	for (i=0;i<a.length;i++)
	    res[i] = a[i] - b[i];

    }

    public static double sqr(double a){
	return ((Sqrarg=(a)) == 0.0 ? 0.0 : Sqrarg*Sqrarg);
    }

    public static double sign(double a, double b){
	return ((b) >= 0.0 ? Math.abs(a) : -Math.abs(a));
    }
    
    public static double pythag(double a, double b){
	double absa, absb;
	absa=Math.abs(a);
	absb=Math.abs(b);
	if (absa > absb) return absa*Math.sqrt(1.0+Matrix.sqr(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*Math.sqrt(1.0+Matrix.sqr(absa/absb)));
    }

    public static void checkSymmetric(Matrix m){
	checkSquare(m);

	if (m.isSymmetric() == false) 
	    Matrix.perror("Matrix " + m.name + " is not symmetric...");
    }

    public static void checkSquare(Matrix m){
	if (m.dimX != m.dimY)
	    Matrix.perror("Matrix " + m.name + " is not square...");
    }

    public static void checkSameDimensions(Matrix m, Matrix n){
	if ( (m.dimX != n.dimX) || (m.dimY != n.dimY) )
	    Matrix.perror("Matrices " + m.name + " and " + n.name + " do not have the same dimensions...");
    }

    public static void checkRowNormal(Matrix m){
	int col = m.dimY;
	int lig = m.dimX;
	int i,j;
	double norm, v;
	boolean ok = true;
	String out = "";

	for (i=0;i<lig;i++){
	    norm = 0.0;
	    for (j=0;j<col;j++){
		v = m.coordinates[i][j];
		norm += (v*v);
	    }

	    if ( (norm < 1.0 - Precision_For_Stochasticity) || (norm > 1.0 + Precision_For_Stochasticity) ){
		ok = false;
		out += "norm for col " + i + " = " + norm + " != 1. ";
	    }
	}

	if (ok == false)
	    Matrix.perror("Matrix " + m.name + " is not row normal: " + out + " (required precision: " + Precision_For_Stochasticity + " )");
    }

    public static void checkColumnStochastic(Matrix m){
	int col = m.dimY;
	int lig = m.dimX;
	int i,j;
	double sum;
	boolean ok = true;
	String out = "";

	for (i=0;i<col;i++){
	    sum = 0.0;
	    for (j=0;j<lig;j++){
		if (m.coordinates[j][i] < 0.0){
		    ok = false;
		    out += "coordinate[" + i + "][" + j + "] <0. ";
		}
		sum += m.coordinates[j][i];
	    }
	    if ( (sum < 1.0 - Precision_For_Stochasticity) || (sum > 1.0 + Precision_For_Stochasticity) ){
		ok = false;
		out += "sum for col " + i + " = " + sum + " != 1. ";
	    }
	}

	if (ok == false)
	    Matrix.perror("Matrix " + m.name + " is not column stochastic: " + out + " (required precision: " + Precision_For_Stochasticity + " )");
    }

    public static void checkRowStochastic(Matrix m){
	int col = m.dimY;
	int lig = m.dimX;
	int i,j;
	double sum;
	boolean ok = true;
	String out = "";

	for (i=0;i<lig;i++){
	    sum = 0.0;
	    for (j=0;j<col;j++){
		if (m.coordinates[i][j] < 0.0){
		    ok = false;
		    out += "coordinate[" + i + "][" + j + "] <0. ";
		}
		sum += m.coordinates[i][j];
	    }
	    if ( (sum < 1.0 - Precision_For_Stochasticity) || (sum > 1.0 + Precision_For_Stochasticity) ){
		ok = false;
		out += "sum for row " + i + " = " + sum + " != 1. ";
	    }
	}

	if (ok == false)
	    Matrix.perror("Matrix " + m.name + " is not row stochastic: " + out + " (required precision: " + Precision_For_Stochasticity + " )");
    }

    public static void checkDoublyStochastic(Matrix m){
	checkRowStochastic(m);
	checkColumnStochastic(m);
    }

    public static void perror(String error_text){
        System.out.println(error_text);
        System.out.println("...now exiting to system...\n");
        System.exit(1);
    }

    
    public static double multiplyComponent (Matrix W, Matrix H, int r, int c){
	//returns (WH)_{rc}
	double v=0.0;
	int i;
	for (i=0;i<W.dimY;i++)
	    v+=(W.coordinates[r][i]*H.coordinates[i][c]);
	return v;
    }

    public static double BregmanError (Matrix V, Matrix W, Matrix H, int param, int alpha){
	//returns Bregman(V||WH)
	int i,j;
	Pnt p = new Pnt(1);
	Pnt q = new Pnt(1);
	double val=0.0, dum;

	if (V.dimX != W.dimX)
	    perror("Matrix.class :: Dimension mismatch (V,W)");
	if (V.dimY != H.dimY)
	    perror("Matrix.class :: Dimension mismatch (V,H)");
	if (W.dimY != H.dimX)
	    perror("Matrix.class :: Dimension mismatch (W,H)");

	for (i=0;i<V.dimX;i++)
	    for (j=0;j<V.dimY;j++){
		p.coordinates[0] = V.coordinates[i][j];
		q.coordinates[0] = multiplyComponent(W,H,i,j);
		dum = Distortion.Bregman(p, q, param, alpha);
		val += dum;
	    }

	return val;
    }

    public static double [] multiplyVector(Matrix m, double [] p){
	//returns vector mp
	if (m.dimY != p.length)
	    perror("Matrix.class :: Dimension mismatch");
	double [] ret = new double [p.length];
	int i, j;
	for (i=0;i<m.dimX;i++){
	    ret[i] = 0.0;
	    for (j=0;j<m.dimY;j++)
		ret[i] += m.coordinates[i][j] * p[j];
	}
	return ret;
    }

    public static Matrix xRotateMatrix(double alpha){
	Matrix m = new Matrix("xRotate",4,4);
	m.coordinates = new double [][]{
	    {1.0, 0.0, 0.0, 0.0},
	    {0.0, Math.cos(alpha), Math.sin(alpha), 0.0},
	    {0.0, -(Math.sin(alpha)), Math.cos(alpha), 0.0},
	    {0.0, 0.0, 0.0, 1.0} };
	return m;
    }

    public static Matrix yRotateMatrix(double alpha){
	Matrix m = new Matrix("yRotate",4,4);
	m.coordinates = new double [][]{
	    {Math.cos(alpha), 0.0, -Math.sin(alpha), 0.0},
	    {0.0, 1.0, 0.0, 0.0},
	    {Math.sin(alpha), 0.0, Math.cos(alpha), 0.0},
	    {0.0, 0.0, 0.0, 1.0} };
	return m;
    }

    public static Matrix zRotateMatrix(double alpha){
	Matrix m = new Matrix("zRotate",4,4);
	m.coordinates = new double [][]{
	    {Math.cos(alpha), Math.sin(alpha), 0.0, 0.0},
	    {-Math.sin(alpha), Math.cos(alpha), 0.0, 0.0},
	    {0.0, 0.0, 1.0, 0.0},
	    {0.0, 0.0, 0.0, 1.0} };
	return m;
    }

    public Matrix (String nm, int dX, int dY){

	if ( (dX == 0) || (dY == 0) ){
	    Matrix.perror("Error(s) in Matrix Constructor");
	}
	
	name = nm;

	if (dX == 1)
	    isVector = true;
	else
	    isVector = false;

	hasDefaultContent = true;

	isSymmetric = true;

	if (dX == dY)
	    isSquared = true;
	else{
	    isSquared = false;
	    isSymmetric = false;
	}

	dimX = dX;
	dimY = dY;
	int i,j;
	
	coordinates = new double[dimX][];
	for (i=0;i<dimX;i++)
	    coordinates[i] = new double[dimY];
	
	for (i=0;i<dimX;i++)
	    for (j=0;j<dimY;j++)
		coordinates[i][j] = DefaultValue;
    }

    public static boolean test_choltmp(Matrix a){
	if (a.dimX != a.dimY)
	    return false;
	int i, j, k, n = a.dimX;
	double sum;
	Matrix test_choleski = new Matrix("Choleski decomposition of " + a.name, n, n);
	for (i=0;i<n;i++)
	    for (j=0;j<n;j++)
		test_choleski.coordinates[i][j] = a.coordinates[i][j]; 

	for (i=0;i<n;i++){
	    for (j=i;j<n;j++){
		for (sum=test_choleski.coordinates[i][j],k=i-1;k>=0;k--)
		    sum -= test_choleski.coordinates[i][k] * test_choleski.coordinates[j][k];
		if (i==j){
		    if (sum <= 0.0)
			return false;
		    test_choleski.coordinates[i][i] = Math.sqrt(sum);
		}else
		    test_choleski.coordinates[j][i] = sum / test_choleski.coordinates[i][i];
	    }
	}
	for (i=0;i<n;i++)
	    for (j=0;j<i;j++)
		test_choleski.coordinates[j][i] = 0.0;
	test_choleski = null;
	return true;
    }

    public static void choltmp(Matrix a){
	if (a.dimX != a.dimY)
	    perror("Dimension mismatch");
	int i, j, k, n = a.dimX;
	double sum;
	a.choleski_decomposition = new Matrix("Choleski decomposition of " + a.name, n, n);
	for (i=0;i<n;i++)
	    for (j=0;j<n;j++)
		a.choleski_decomposition.coordinates[i][j] = a.coordinates[i][j];

	for (i=0;i<n;i++){
	    for (j=i;j<n;j++){
		for (sum=a.choleski_decomposition.coordinates[i][j],k=i-1;k>=0;k--)
		    sum -= a.choleski_decomposition.coordinates[i][k] * a.choleski_decomposition.coordinates[j][k];
		if (i==j){
		    if (sum <= 0.0)
			perror(a.choleski_decomposition.name + " is not positive definite");
		    a.choleski_decomposition.coordinates[i][i] = Math.sqrt(sum);
		}else
		    a.choleski_decomposition.coordinates[j][i] = sum / a.choleski_decomposition.coordinates[i][i];
	    }
	}
	for (i=0;i<n;i++)
	    for (j=0;j<i;j++)
		a.choleski_decomposition.coordinates[j][i] = 0.0;
    }

    /*********************************************************************************************************
     * Instance methods
     *****/


    public String afficheComplet(int n){
	//prints a n x n submatrix sans arrondi
	int i, j, v;
	String val = "";

	if (n!=-1)
	    v = Math.min(n, Math.min(dimX, dimY));
	else
	    v = Math.min(dimX, dimY);

	for (i=0;i<v;i++){
	    for (j=0;j<v;j++)
		val += coordinates[i][j] + " ";
	    val += "\n";
	}
	return val;
    }

    public String affiche(int n){
	//prints a n x n submatrix
	int i, j, v;
	String val = "";

	if (n!=-1)
	    v = Math.min(n, Math.min(dimX, dimY));
	else
	    v = Math.min(dimX, dimY);

	for (i=0;i<v;i++){
	    for (j=0;j<v;j++)
		val += Matrix.DF.format(coordinates[i][j]) + " ";
	    val += "\n";
	}
	return val;
    }

    public String toString(){
	int i, j;
	String val = "Matrix " + name + "\n";

	for (i=0;i<dimX;i++){
	    for (j=0;j<dimY;j++)
		val += Matrix.DF.format(coordinates[i][j]) + " ";
	    val += "\n";
	}
	return val;
    }

    public void outPrint(int n){
	//prints a n x n submatrix
	int i, j, v;
	String val = "";

	System.out.println("Matrix " + name);

	if (n!=-1)
	    v = Math.min(n, Math.min(dimX, dimY));
	else
	    v = Math.min(dimX, dimY);

	for (i=0;i<v;i++){
	    for (j=0;j<v;j++)
		System.out.print(coordinates[i][j] + " ");
	    System.out.print("\n");
	}

	System.out.println("\n\n");
    }

    public void copyOnThis(Matrix m){
	if ( (dimX != m.dimX) || (dimY != m.dimY) )
	    Matrix.perror("Dimension mismatch to copy matrices");

	int i,j;
	for (i=0;i<dimX;i++)
	    for (j=0;j<dimY;j++)
		coordinates[i][j] = m.coordinates[i][j];

	hasDefaultContent = m.hasDefaultContent;
    }

    public boolean isSymmetric(){
	int i,j;
	double m;
	if (dimX != dimY)
	    return false;

	for (i=0;i<dimX-1;i++)
	    for (j=i+1;j<dimY;j++)
		if (Math.abs(coordinates[i][j] - coordinates[j][i]) > Precision_For_Eigensystems){
		    return false;
		}
		else if (coordinates[i][j] != coordinates[j][i]){
		    m = (coordinates[i][j] + coordinates[j][i]) / 2.0;
		    coordinates[i][j] = m;
		    coordinates[j][i] = m;
		}
	return true;
    }

    public void computeSymmetric(){
	int i,j;
	double m;
	boolean found = false;
	if (dimX != dimY)
	    isSymmetric = false;

	for (i=0;i<dimX-1;i++)
	    for (j=i+1;j<dimY;j++)
		if (Math.abs(coordinates[i][j] - coordinates[j][i]) > Precision_For_Eigensystems){
		    isSymmetric = false;
		    found = true;;
		}
		else if (coordinates[i][j] != coordinates[j][i]){
		    m = (coordinates[i][j] + coordinates[j][i]) / 2.0;
		    coordinates[i][j] = m;
		    coordinates[j][i] = m;
		}
	if (found == false)
	    isSymmetric = true;
    }

    public synchronized void toD(Domain d, Matrix W){
	int dim;
	double sum;
	int i,j, nmade = 0, percent ;

	if (dimX != dimY)
	    Matrix.perror("Cannot generate D on a NON SQUARE matrix");

	dim = dimX;

	if ( (dimX != W.dimX) || (dimY != W.dimY) )
	    Matrix.perror("Dimension mismatch between W and the number of selected Genes");

	if (AGCT.Debug) System.out.print("Generating matrix D (" + dim + " x " + dim + ")... ");
	d.myAGCT.myInformationFrame.setTextProgressBar("computing " + name);

	for (i=0;i<dim;i++){
	    sum = 0.0;
	    for (j=0;j<dim;j++){
		percent = (int) (100.0 * ((double) nmade / (double) (dim*dim)));
		d.myAGCT.myInformationFrame.setValueProgressBar(percent);

		sum+=W.coordinates[i][j];
		if (i!=j)
		    coordinates[i][j] = 0.0;
		nmade++;
	    }
	    coordinates[i][i] = sum;
	}

	Matrix.checkSymmetric(this);

	if (AGCT.Debug) System.out.println("ok.");
	d.myAGCT.myInformationFrame.setValueProgressBar(0);
	d.myAGCT.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);
    }
    
    public synchronized void sparsifyWithDistances(AGCT ap){
	//
	if (AGCT.Debug) System.out.print("Filtering W with Distances (" + dimX + " x " + dimY + ")... ");

	String s, currentStep = "", valret;
	Vector ge = null, co;
	int dim = ap.myDomain.numberSelectedGenes;
	int i, j;
	double dist, dum, midu = 1.0, madu = 0.0, fact;
	double nb = 0.0, nbr = 0.0;

	double mindist = 0;

	AGCTCounter cc = new AGCTCounter(ap.myInformationFrame, "Sparsifying with distances " + name + " with Method_D = " + AGCT.Method_D, dim * (dim - 1) / 2);
	for (i=0;i<dim-1;i++){
	    for (j=i+1;j<dim;j++){
		dist = distance( ((Vector) ((Vector) ap.allGenesCoordinates.elementAt(ap.geneToCoord[i])).elementAt(1)), ((Vector) ((Vector) ap.allGenesCoordinates.elementAt(ap.geneToCoord[j])).elementAt(1)));
		nb += 1.0;

		if ( (AGCT.Method_D == 0) && (dist > AGCT.DISTANCE_PARAMETER) ){
		    coordinates[i][j] = 0.0;
		    coordinates[j][i] = 0.0;
		    
		    nbr += 1.0;
		}else if (AGCT.Method_D == 1){
		    dum = coordinates[i][j];
		    fact = Math.exp(-(dist*dist)/AGCT.DISTANCE_PARAMETER);
		    coordinates[i][j] = dum * fact;
		    coordinates[j][i] = dum * fact;
		    if (fact < midu)
			midu = fact;
		    if (fact > madu)
			madu = fact;
		}

		if ( (nbr == 1.0) || (dist < mindist) )
		    mindist = dist;

		cc.increment();
	    }
	}
	cc.end();
	if (AGCT.Method_D == 0)
	    System.out.println("\n -- Zeroed " + nbr + " out of " + nb + " values in W = " + (nbr / nb)*100.0 + " % --- min dist = " + mindist);
	if (AGCT.Method_D == 1)
	    System.out.println("\n -- Dampening factor in interval [" + midu + ", " + madu + "]");

	if (AGCT.Debug) System.out.println("ok.");
    }

    public double distance(Vector ge1, Vector ge2){
	double x1 = ((Double) ge1.elementAt(0)).doubleValue();
	double x2 = ((Double) ge2.elementAt(0)).doubleValue();
	double y1 = ((Double) ge1.elementAt(1)).doubleValue();
	double y2 = ((Double) ge2.elementAt(1)).doubleValue();

	return Math.sqrt( ( (x1 - x2) * (x1 - x2) ) + ( (y1 - y2) * (y1 - y2) ) );
    }

    public synchronized void Natalia_toW(AGCT ag) {// 近さに基づきマトリックスを作る
 
	if (!AGCT.SUBMISSION)
	    System.out.println(" ** Warning : Computation of W Using Natalia's algorithm ** ");

       int dim = dimX, nmade = 0, percent;
        double kappa;
        if (dimX != dimY)
            Matrix.perror("Cannot generate W on a NON SQUARE matrix");


	if (dimX != ag.myDomain.numberSelectedGenes)
	    Matrix.perror("Dimension mismatch between W and the number of selected Genes");

        //if (dimX != ag.data.getMyDomain().numberSelectedGenes)
	//  Matrix.perror("Dimension mismatch between W and the number of selected Genes");

	if (AGCT.Debug) System.out.print("Generating matrix W (" + dim + " x " + dim + ") ... ");
	ag.myInformationFrame.setTextProgressBar("computing " + name);

        Gene[] genes = new Gene[dim];
        for (int i = 0; i < dim; i++) {
            genes[i] = (Gene) ag.myDomain.domainGenes.elementAt(ag.myDomain.selectedGeneNumberToGeneNumber[i]);
	    //gi = (Gene) ag.myDomain.domainGenes.elementAt(ag.myDomain.selectedGeneNumberToGeneNumber[i]);
        }

        if (AGCT.Method_W == 0) {
            for (int i = 0; i < dim; i++) {
                percent = (int) (100.0 * ((double) nmade / (double) (dim * dim)));
 		ag.myInformationFrame.setValueProgressBar(percent);

                for (int j = 0; j < dim; j++) {
                    set(i, j, Gene.getSimilarity(genes[i], genes[j]));
                    nmade++;
                }
            }
        } else if (AGCT.Method_W == 1) {//default
	    AGCTCounter cc = new AGCTCounter(ag.myInformationFrame, "computing " + name, dim);

            for (int i = 0; i < dim; i++) {
                class A implements Comparable<A> {
                    double d;
                    int i;

                    A(double d, int i) {
                        this.d = d;
                        this.i = i;
                    }

                    public int compareTo(A o) {
                        return -(int) Math.signum(d - o.d);
                    }

                    ;

                    @Override
                    public String toString() {
                        return "A{" + "d=" + d + ", i=" + i + '}';
                    }
                }
                A[] as = new A[dim];
                for (int j = 0; j < dim; j++) {
                    as[j] = new A(Gene.getSimilarity(genes[i], genes[j]), j);
                    if (i == j) as[j] = new A(Double.POSITIVE_INFINITY, j);
                }
                Arrays.sort(as);

                for (int j = 0; j < AGCT.Number_Of_Neighbors; j++)
                    if (j + 1 < dim && as[j + 1].i < dimY) {
                        // avoid to Wii > 0
                        set(i, as[j + 1].i, 1.0);
                        set(as[j + 1].i, i, 1.0);
                    } else if (AGCT.Debug)
                        System.out.println("!?");
                cc.increment();
            }

            cc.end();
        }

        if (AGCT.Method_N == 3) {
            ag.myScaling.memorize(this);
            kappa = ag.myScaling.kappa();
            getScaling(kappa);
        }

        hasDefaultContent = false;
        computeSymmetric();
        assert isSymmetric();

        if (AGCT.Debug)
            System.out.println("ok.");
	ag.myInformationFrame.setValueProgressBar(0);
	ag.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);
    }


    public synchronized void toW(AGCT ag){
	Gene gi, gj;
	int dim = dimX, nmade = 0, percent, i, j, k;
	double kappa;
	if (dimX != dimY)
	    Matrix.perror("Cannot generate W on a NON SQUARE matrix");

	if (dimX != ag.myDomain.numberSelectedGenes)
	    Matrix.perror("Dimension mismatch between W and the number of selected Genes");

	if (AGCT.Debug) System.out.print("Generating matrix W (" + dim + " x " + dim + ") ... ");
	ag.myInformationFrame.setTextProgressBar("computing " + name);

	if (AGCT.Method_W == 0){
	    for (i=0;i<dim;i++){
		percent = (int) (100.0 * ((double) nmade / (double) (dim*dim)));
		ag.myInformationFrame.setValueProgressBar(percent);
		gi = (Gene) ag.myDomain.domainGenes.elementAt(ag.myDomain.selectedGeneNumberToGeneNumber[i]);

		for (j=0;j<dim;j++){
		    gj = (Gene) ag.myDomain.domainGenes.elementAt(ag.myDomain.selectedGeneNumberToGeneNumber[j]);
		    coordinates[i][j] = Gene.getSimilarity(gi, gj);

		    if ( (j==i) && (coordinates[i][j] == 0) )
			Matrix.perror("Zero coordinate in diagonal for " + i);

		    nmade++;
		}
	    }
	}else if (AGCT.Method_W == 1){
	    double [] ordsimil = new double[dim];
	    int [] ordindex = new int[dim];
	    int it;
	    double dt; 

	    AGCTCounter cc = new AGCTCounter(ag.myInformationFrame, "computing " + name, dim);


	    for (i=0;i<dim;i++){
		gi = (Gene) ag.myDomain.domainGenes.elementAt(ag.myDomain.selectedGeneNumberToGeneNumber[i]);
		for (j=0;j<dim;j++){
		    gj = (Gene) ag.myDomain.domainGenes.elementAt(ag.myDomain.selectedGeneNumberToGeneNumber[j]);
		    ordsimil[j] =  Gene.getSimilarity(gi, gj);
		    ordindex[j] = j;
		}
		
		QuickSort.quicksort(ordsimil, ordindex);
		//put in decreasing order
		for (j=0;j<dim/2;j++){
		    dt = ordsimil[j];
		    ordsimil[j] = ordsimil[dim-1-j];
		    ordsimil[dim-1-j] = dt;
		    
		    it = ordindex[j];
		    ordindex[j] = ordindex[dim-1-j];
		    ordindex[dim-1-j] = it;
		}

		/*for (j=0;j<dim-1;j++){
		    for (k=j+1;k<dim;k++){
			if (ordsimil[j]<ordsimil[k]){
			    dt = ordsimil[j];
			    ordsimil[j] = ordsimil[k];
			    ordsimil[k] = dt;

			    it = ordindex[j];
			    ordindex[j] = ordindex[k]; 
			    ordindex[k] = it;
			}
		    }
		    //cc.increment();
		    }*/
		
		for (j=0;j<AGCT.Number_Of_Neighbors;j++){
		    coordinates[i][ordindex[j+1]] = 1.0;
		    coordinates[ordindex[j+1]][i] = 1.0;
		}
		cc.increment();
	    }

	    cc.end();
	}

	if (AGCT.Filter_Similarities_With_Distances)
	    sparsifyWithDistances(ag);

	if (AGCT.Method_N == 3){
	    ag.myScaling.memorize(this);
	    kappa = ag.myScaling.kappa();
	    getScaling(kappa);
	}
	//getScaling();
	//last checks

	hasDefaultContent = false;
	computeSymmetric();
	Matrix.checkSymmetric(this);

	if (AGCT.Debug) System.out.println("ok.");
	ag.myInformationFrame.setValueProgressBar(0);
	ag.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);
    }

    public synchronized void toC(Domain d, boolean beforeSelection){
	int dim = dimX, nmade = 0, percent, i, j, k, x, y, z, rx, ry;
	double slx, sly, avex, avey;
	Gene gg;
	if (dimX != dimY)
	    Matrix.perror("Cannot generate C on a NON SQUARE matrix");

	if (dimX != d.myAGCT.dimFeatures)
	    Matrix.perror("Dimension mismatch");

	if (d.domainGenes != null){

	    AGCTCounter cc = new AGCTCounter(d.myAGCT.myInformationFrame, "computing " + name, (2*dimX*d.numberSelectedGenes) + (dimX * (dimX - 1) * d.numberSelectedGenes / 2) );

	    if (AGCT.Debug) System.out.print("Generating matrix C (" + dim + " x " + dim + ") ... ");

	    d.myAGCT.average_features = new double[dimX];
	    d.myAGCT.sigma_features = new double[dimX];

	    for (x=0;x<dimX;x++){
		d.myAGCT.average_features[x]=0.0;
		for (y=0;y<d.numberSelectedGenes;y++){
		    gg = d.getSelectedGene(y); //(Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[y]);
		    d.myAGCT.average_features[x] += gg.getFinalCoordinates(x,beforeSelection);

		    cc.increment();
		}
	    }

	    for (x=0;x<dimX;x++)
		d.myAGCT.average_features[x] /= (double) d.numberSelectedGenes;

	    for (x=0;x<dimX;x++){
		d.myAGCT.sigma_features[x]=0.0;
		for (y=0;y<d.numberSelectedGenes;y++){
		    gg = d.getSelectedGene(y); //(Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[y]);
		    d.myAGCT.sigma_features[x] += ( (gg.getFinalCoordinates(x,beforeSelection) - d.myAGCT.average_features[x]) 
						    * (gg.getFinalCoordinates(x,beforeSelection) - d.myAGCT.average_features[x]) );

		    cc.increment();		}
	    }

	    for (x=0;x<dimX;x++)
		d.myAGCT.sigma_features[x] = Math.sqrt(d.myAGCT.sigma_features[x]);

	    for (x=0;x<dimX;x++){
		for (y=x;y<dimX;y++){
		    coordinates[x][y] = 0.0;
		    if (y!=x)
			coordinates[y][x] = 0.0;
		    for (z=0;z<d.numberSelectedGenes;z++){
			gg = d.getSelectedGene(z); //(Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[z]);
			slx = gg.getFinalCoordinates(x,beforeSelection);
			sly = gg.getFinalCoordinates(y,beforeSelection);
			avex = d.myAGCT.average_features[x];
			avey = d.myAGCT.average_features[y];

			coordinates[x][y] += ( (slx - avex) * (sly - avey ) );

			cc.increment();
		    }

		    coordinates[x][y] /= ( d.myAGCT.sigma_features[x] * d.myAGCT.sigma_features[y] );
		    if ( (d.myAGCT.sigma_features[x] == 0.0) || (d.myAGCT.sigma_features[y] == 0.0) ){
			if (x!=y)
			    coordinates[x][y] = 0.0;
			else
			    coordinates[x][y] = ( (double) d.numberSelectedGenes );
		    }
		    if (y!=x)
			coordinates[y][x] = coordinates[x][y];
		}
	    }

	    Matrix.checkSymmetric(this);
	    cc.end();
	}

	if (AGCT.Debug) System.out.println("ok.");
    }
    

    public synchronized void toE(Domain d, Matrix cop){
	if (dimX != dimY)
	    Matrix.perror("Cannot generate E on a NON SQUARE matrix");

	if (dimX != d.myAGCT.dimFeatures)
	    Matrix.perror("Dimension mismatch");

	d.myAGCT.pca_eigenvalues = new double [dimX];
	double [] e = new double[dimX];

	copyOnThis(cop);

	//System.out.println(this);

	tred2(d,d.myAGCT.pca_eigenvalues,e);
	tqli(d,d.myAGCT.pca_eigenvalues,e);

	eigenSort(d.myAGCT.pca_eigenvalues);
	transposeSquare();

	//checkOrthonormality(d,-1);
	//Matrix.checkEigensystem(d,cop,this,d.myAGCT.pca_eigenvalues,-1);
    }

    public synchronized void toV(Domain d){
	double no, nol, ave, sig, ps, slo, cos, pco;
	int i, x, y;
	Gene gg;

	if (dimX != dimY)
	    Matrix.perror("Cannot generate V on a NON SQUARE matrix");

	if (dimX != d.myAGCT.dimFeatures)
	    Matrix.perror("Dimension mismatch");

	if (d.domainGenes != null){

	    AGCTCounter cc = new AGCTCounter(d.myAGCT.myInformationFrame, "computing " + name, dimX);

	    for (i=0;i<dimX;i++){
		no = 0.0;
		for (x=0;x<d.numberSelectedGenes;x++){
		    gg = (Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[x]);
		    pco = gg.pca_components.coordinates[i];
		    pco *= pco;
		    no += pco;
		}

		no = Math.sqrt(no);

		for (y=0;y<dimX;y++){
		    ave = d.myAGCT.average_features[y];
		    sig = d.myAGCT.sigma_features[y];

		    nol = 0.0;
		    ps = 0.0;

		    for (x=0;x<d.numberSelectedGenes;x++){
			gg = (Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[x]);
			slo = gg.getFinalCoordinates(y,true);

			if (sig > 0.0){
			    nol += ( ( (slo - ave ) / sig ) * ( (slo - ave ) / sig ) );
			    ps += ( (slo - ave ) / sig ) * gg.pca_components.coordinates[i];
			}
		    }

		    if (sig > 0.0){
			nol = Math.sqrt(nol);
			cos = ps/(no * nol);
		    }else{
			cos = 0.0;
		    }

		    coordinates[y][i] = cos;
		}
		cc.increment();
	    }
	    cc.end();
	}
    }

    public synchronized void toVV(Domain d, Clustering cl){
	if (d.domainGenes == null)
	    Matrix.perror("No Genes !");
	
	AGCTCounter cc = new AGCTCounter(d.myAGCT.myInformationFrame, "computing " + name, dimX *  dimY * d.numberSelectedGenes);
	if (AGCT.Debug) System.out.print("Generating matrix VV (" + dimX + " x " + dimY + ") ... ");
	
	int i, j, k;
	Gene gg;
	
	double K = (double) dimX;
	double N = (double) d.numberSelectedGenes;
	Matrix soft = cl.myClusteringAlgorithm.soft_memberships;
	
	for (i=0;i<dimX;i++)
	    for (j=0;j<dimY;j++)
		coordinates[i][j] = 0.0;
	
	for (i=0;i<dimX;i++){
	    for (j=0;j<dimY;j++){
		for (k=0;k<d.numberSelectedGenes;k++){
		    gg = (Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumber[k]);
		    
		    coordinates[i][j] += soft.coordinates[i][k] * gg.pca_components.coordinates[j] * K / N;
		    
		    cc.increment();
		}
	    }
	}
	cc.end();
	if (AGCT.Debug) System.out.println("ok.");
    }


    public synchronized void getScaling(double kappa){
	//Doubly stochastic scaling
	if (AGCT.Debug)
	    System.out.println("Doubly Stochastic Scaling for Matrix " + name + " with value kappa = " + DF.format(kappa) + "...");
	int i, j;
	for (i=0;i<dimX;i++)
	    for (j=0;j<dimY;j++)
		coordinates[i][j] *= kappa;
    }

    public synchronized void toN(AGCT a, Matrix D, Matrix W){
	int dim = dimX, nmade = 0, percent;

	if (dimX != dimY)
	    Matrix.perror("Cannot generate N on a NON SQUARE matrix");

	if (dimX != a.myDomain.numberSelectedGenes)
	    Matrix.perror("Dimension mismatch between N and the number of selected Genes");

	if (dimX != D.dimX)
	    Matrix.perror("Dimension mismatch between N and D");

	if (dimX != W.dimX)
	    Matrix.perror("Dimension mismatch between N and W");

	if (AGCT.Debug) System.out.print("Generating matrix N (" + dim + " x " + dim + ")... ");
	a.myInformationFrame.setTextProgressBar("computing " + name);

	int x,y, max = (dim * dim);
	if (AGCT.Method_N == 1)
	    max *= 2;
	
	if (AGCT.Method_N == 1)
	    a.DWD = new Matrix("DWD", dim, dim);
  
	if ( ( AGCT.Method_N == 0 ) || ( AGCT.Method_N == 1 ) ){
	    for (x=0;x<dim;x++){
		for (y=0;y<dim;y++){
		    if (AGCT.Method_N == 0){
			coordinates[x][y] = W.coordinates[x][y]/(Math.sqrt(D.coordinates[x][x])*Math.sqrt(D.coordinates[y][y]));
		    }else if (AGCT.Method_N == 1){
			a.DWD.coordinates[x][y] = W.coordinates[x][y]/(Math.sqrt(D.coordinates[x][x])*Math.sqrt(D.coordinates[y][y]));
		    }
		    
		    percent = (int) (100.0 * ((double) nmade / (double) (max)));
		    a.myInformationFrame.setValueProgressBar(percent);
		    nmade++;
		}
	    }
    
	    if (AGCT.Method_N == 1){    
		for (x=0;x<dim;x++){
		    if (D.coordinates[x][x] == 0.0)
			Matrix.perror("Empty diagonal in " + x);

		    for (y=0;y<dim;y++){
			coordinates[x][y] = W.coordinates[x][y]/D.coordinates[x][x];
	  
			percent = (int) (100.0 * ((double) nmade / (double) (max)));
			a.myInformationFrame.setValueProgressBar(percent);
			nmade++;
		    }
		}
	    }
	}else if ( (AGCT.Method_N == 2) || (AGCT.Method_N == 3) ){
	    for (x=0;x<dim;x++){
		for (y=0;y<dim;y++){
		    coordinates[x][y] = W.coordinates[x][y];

		    percent = (int) (100.0 * ((double) nmade / (double) (max)));
		    a.myInformationFrame.setValueProgressBar(percent);
		    nmade++;
		}
	    }
	}else if ( (AGCT.Method_N == 4) || (AGCT.Method_N == 5) ){
	    if (AGCT.Method_W != 1)
		Matrix.perror("Matrix.class :: to use Bethe Hessian, please filter W with Symmetric NN");
	    
	    double rc = Math.sqrt(AGCT.Bethe_Hessian_Factor * AGCT.Number_Of_Neighbors);
	    if (AGCT.Method_N == 5)
		rc = -rc;

	    for (x=0;x<dim;x++){
		for (y=0;y<dim;y++){
		    coordinates[x][y] = - (rc * W.coordinates[x][y]);
		    if (x == y){
			coordinates[x][y] += ( ( rc * rc ) - 1.0 ) + D.coordinates[x][x];
		    }

		    //flips the sign so that the search for eigenvalues still picks the largest
		    coordinates[x][y] = -coordinates[x][y];
		}
	    }
	}

	a.myInformationFrame.setValueProgressBar(0);
	a.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);

	if ( (AGCT.Method_N == 2) || (AGCT.Method_N == 3) )
	    fastDoublyStochasticApproximation(a.myDomain);

	if ( (AGCT.Method_N == 0) || (AGCT.Method_N == 4) || (AGCT.Method_N == 5) ){
	    Matrix.checkSymmetric(this);
	}else if (AGCT.Method_N == 1){
	    Matrix.checkRowStochastic(this);
	    Matrix.checkSymmetric(a.DWD);
	}else if ( (AGCT.Method_N == 2) || (AGCT.Method_N == 3) ){
	    Matrix.checkDoublyStochastic(this);
	}

	if (AGCT.Debug) System.out.println("ok.");
    }

    public void sparsifyMatrix(){
	if (AGCT.Debug) System.out.print("Sparsifying matrix (" + (AGCT.Sparsify_Statistic*100) + " % removed)... ");

	int td = dimX * dimY, i, j, ndiff = 0;
	double[] all_values = new double [td];
	int[] cardinals = new int [td];
	double[] values = new double [td];
	double nremove, curc, valr;

	for (i=0;i<td;i++){
	    all_values[i] = 0.0;
	    cardinals[i] = 0;
	    values[i] = 0.0;
	}

	ndiff = 0;
	for (i=0;i<dimX;i++)
	    for (j=0;j<dimY;j++){
		all_values[ndiff] = coordinates[i][j];
		ndiff++;
	    }

	nremove = AGCT.Sparsify_Statistic * (double) ndiff;
	QuickSort.quicksort(all_values);
	ndiff = 0;
	i = 0;
	do{
	    if (i==0){
		values[0] = all_values[0];
		cardinals[0]++;
	    }else{
		if (all_values[i] != values[ndiff]){
		    ndiff++;
		    values[ndiff] = all_values[i];
		}
		cardinals[ndiff]++;
	    }
	    i++;
	}while(i<td);

	curc = 0.0;
	i = 0;
	valr = 0.0;
	do{
	    if ( (i==0) || (curc < nremove) ){
		valr = values[i];
		curc += (double) cardinals[i];

		//System.out.println(" i = " + i + " ; nremove = " + nremove + " ; curc = " + curc + " valr = " + valr);

		i++;
	    }
	}while( curc < nremove );
	valr = values[i];

	for (i=0;i<dimX;i++)
	    for (j=0;j<dimY;j++){
		if (coordinates[i][j] < valr)
		    coordinates[i][j] = 0.0;
	    }

	if (AGCT.Debug) System.out.print("ok.\n");
    }

    public synchronized void fastDoublyStochasticApproximation(Domain d){
	d.myAGCT.myInformationFrame.setTextProgressBar("Fast DSA of " + name );

	if (AGCT.Debug) System.out.print("Computing the fast doubly stochastic approximation for system " + name + " (" + dimX + " x " + dimY + ")... ");
	double lastfrob = 0.0, lastmin = 0.0, fs;
	int niter = 0, niterid = 0, i, j, n = dimX, ratio;

	int lastp = 100 + Inper;

	Matrix.checkSquare(this);

	boolean [] f_Stop = new boolean[1];
	double [] f_Sumlig = new double[n];
	double [] f_Sumcol = new double[n];
	double [] f_Sumgen = new double[1];
	double [] f_Frob = new double[1];
	double [] f_Min = new double[1];
	fs = 0.0;

	f_Stop[0] = false;

	for (i=0;i<n;i++){
	    f_Sumlig[i] = 0.0;
	    f_Sumcol[i] = 0.0;
	}

	for (i=0;i<n;i++){
	    for (j=0;j<n;j++){
		f_Sumlig[i] += coordinates[i][j];
		f_Sumcol[i] += coordinates[j][i];
		fs += coordinates[i][j];
	    }
	}

	f_Sumgen[0] = fs;

	do{
	    fastDsaP1P2(f_Stop, f_Sumlig, f_Sumcol, f_Sumgen, f_Frob, f_Min);


	    //if ( (AGCT.Debug) && (niter % 10 == 0) ) 
	    //System.out.print("\n. Frob = " + f_Frob[0] + " Min = " + f_Min[0] + " (prec. = " + Precision_For_Doubly_Stochastic_Approximation + ", niterid = " + niterid + " ");

	    niter++;

	    //if (niter%100 == 0)
	    //System.out.println("Current niter = " + niter);

	    if ( niter == 1 ){
		lastfrob = f_Frob[0];
		lastmin = f_Min[0];
		niterid = 0;
	    }else{
		if ( (lastfrob > f_Frob[0] - Precision_For_Doubly_Stochastic_Approximation) && 
		     (lastfrob < f_Frob[0] + Precision_For_Doubly_Stochastic_Approximation) && 
		     (lastmin > f_Min[0] - Precision_For_Doubly_Stochastic_Approximation) && 
		     (lastmin < f_Min[0] + Precision_For_Doubly_Stochastic_Approximation) )
		    niterid ++;
		else{
		    if (lastfrob != f_Frob[0]){
			lastfrob = f_Frob[0];
			niterid = 0;
		    }
		    if (lastmin != f_Min[0]){
			lastmin = f_Min[0];
			niterid = 0;
		    }
		}
	    }

	    //System.out.println("Niterid = " + niterid + "; DS_Max = " + Doubly_Stochastic_Iterations_Id_Max + ".");

	    if (niterid == Doubly_Stochastic_Iterations_Id_Max)
		f_Stop[0] = true;

	    d.myAGCT.myInformationFrame.setValueProgressBar( (int) (100.0 * (double) niterid / (double) Doubly_Stochastic_Iterations_Id_Max ) );
	    d.myAGCT.myInformationFrame.setTextProgressString(niterid + " out of " + Doubly_Stochastic_Iterations_Id_Max + " in plateau");

	}while(f_Stop[0] == false);

	d.myAGCT.myInformationFrame.resetTextProgressString();
	d.myAGCT.myInformationFrame.setValueProgressBar(0);
	d.myAGCT.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);
	if (AGCT.Debug) System.out.print(" [#Iterations : " + niter + "] ok.\n");
    }

    public synchronized void fastDsaP1P2(boolean [] f_Stop, double [] f_Sumlig, double [] f_Sumcol, double [] f_Sumgen, double [] f_Frob, double [] f_Min){

	int i,j, n = dimX;
	double nextsumgen = 0.0;
	double v, num, den, dn = (double) n;

	f_Frob[0] = 0.0;
	f_Min[0] = 0.0;

	boolean allpos = true;

	double [] nextsumlig = new double[n];
	double [] nextsumcol = new double[n];

	for (i=0;i<n;i++){
	    nextsumlig[i] = f_Sumlig[i];
	    nextsumcol[i] = f_Sumcol[i];
	}
	nextsumgen = f_Sumgen[0];

	for (i=0;i<n;i++){
	    for (j=0;j<n;j++){
		num = (coordinates[i][j] * dn * dn) + dn + f_Sumgen[0] - ( (f_Sumlig[i] + f_Sumcol[j]) * dn);
		den = dn * dn;

		v = ( num / den );

		f_Frob[0] += ( (v - coordinates[i][j]) * (v - coordinates[i][j]) );

		nextsumgen -= coordinates[i][j];
		nextsumlig[i] -= coordinates[i][j];
		nextsumcol[j] -= coordinates[i][j];

		if (v < 0.0){
		    allpos = false;
		    if (v < f_Min[0])
			f_Min[0] = v;
		    coordinates[i][j] = 0.0;
		}else{
		    coordinates[i][j] = v;

		    nextsumgen += v;
		    nextsumlig[i] += v;
		    nextsumcol[j] += v;
		}
	    }
	}

	if (allpos == true)
	    f_Stop[0] = true;

	for (i=0;i<n;i++){
	    f_Sumlig[i] = nextsumlig[i];
	    f_Sumcol[i] = nextsumcol[i];
	}
	f_Sumgen[0] = nextsumgen;
    }

    public synchronized void tred2(Domain dd, double d[], double e[]){
	AGCTCounter cc = new AGCTCounter(dd.myAGCT.myInformationFrame, "Householder (" + name + ")", 2 * dimX - 2);

	if (AGCT.Debug) System.out.print("Householder tridiagonalization on Matrix " + name + " (" + dimX + " x " + dimY + ") (start)... ");

	Matrix.checkSymmetric(this);

	int l,k,j,i, n = dimX;
	double scale,hh,h,g,f;

	for (i=n;i>=2;i--){
	    cc.increment();

	    l=i-1;
	    h=scale=0.0;
	    if (l > 1){
		for (k=1;k<=l;k++)
		    scale += Math.abs(coordinates[i-1][k-1]);
		if (scale == 0.0)
		    e[i-1] = coordinates[i-1][l-1];
		else {
		    for (k=1;k<=l;k++){
			coordinates[i-1][k-1] /= scale;
			h += coordinates[i-1][k-1]*coordinates[i-1][k-1];
		    }
		    f=coordinates[i-1][l-1];
		    g=(f >= 0.0 ? -Math.sqrt(h) : Math.sqrt(h));
		    e[i-1]=scale*g;
		    h -= f*g;
		    coordinates[i-1][l-1] = f-g;
		    f=0.0;
		    for (j=1;j<=l;j++){
			coordinates[j-1][i-1] = coordinates[i-1][j-1]/h;
			g=0.0;
			for (k=1;k<=j;k++)
			    g += coordinates[j-1][k-1]*coordinates[i-1][k-1];
			for (k=j+1;k<=l;k++)
			    g += coordinates[k-1][j-1]*coordinates[i-1][k-1];
			e[j-1] = g/h;
			f += e[j-1]*coordinates[i-1][j-1];
		    }
		    hh = f/(h+h);
		    for (j=1;j<=l;j++){
			f=coordinates[i-1][j-1];
			e[j-1]=g=e[j-1]-hh*f;
			for (k=1;k<=j;k++)
			    coordinates[j-1][k-1] -= (f*e[k-1] + g*coordinates[i-1][k-1]);
		    }
		}
	    }else
		e[i-1]=coordinates[i-1][l-1];
	    d[i-1]=h;
	}
	
	d[0]=0.0;
	e[0]=0.0;
	for (i=1;i<=n;i++) {
	    cc.increment();
	    
	    l=i-1;
	    if (d[i-1] != 0.0){
		for (j=1;j<=l;j++){
		    g=0.0;
		    for (k=1;k<=l;k++)
			g += coordinates[i-1][k-1]*coordinates[k-1][j-1];
		    for (k=1;k<=l;k++)
			coordinates[k-1][j-1] -= g*coordinates[k-1][i-1];
		}
	    }
	    d[i-1]=coordinates[i-1][i-1];
	    coordinates[i-1][i-1]=1.0;
	    for (j=1;j<=l;j++) coordinates[j-1][i-1]=coordinates[i-1][j-1]=0.0;
	}
	
	if (AGCT.Debug) System.out.print("ok.\n");
	cc.end(); 
    }


    public synchronized void tqli(Domain ddd, double d[], double e[]){
	AGCTCounter cc = new AGCTCounter(ddd.myAGCT.myInformationFrame, "QL-Implicit Shifts (" + name + ")", dimX);

	Matrix.checkSquare(this);

	int m,l,iter,i,k, n = dimX;
	double s,r,p,g,f,dd,c,b,sf;
	for (i=2;i<=n;i++) e[i-2]=e[i-1];
	e[n-1]=0.0;

	for (l=1;l<=n;l++){
	    cc.increment();

	    iter=0;
	    do{
		for (m=l;m<=n-1;m++){
		    dd=Math.abs(d[m-1])+Math.abs(d[m]);
		    sf= (double) (Math.abs(e[m-1])+dd);
		    if ( sf == dd )
			break;
		}
		if (m != l){
		    if (iter++ == Numerical_Iterations_Max) Matrix.perror("Too many iterations in tqli");
		    g=(d[l]-d[l-1])/(2.0*e[l-1]);
		    r=Matrix.pythag(g,1.0);
		    g=d[m-1]-d[l-1]+e[l-1]/(g+Matrix.sign(r,g));
		    s=c=1.0;
		    p=0.0;
		    for (i=m-1;i>=l;i--){
			f=s*e[i-1];
			b=c*e[i-1];
			e[i]=(r=Matrix.pythag(f,g));
			if (r == 0.0) {
			    d[i] -= p;
			    e[m-1]=0.0;
			    break;
			}
			s=f/r;
			c=g/r;
			g=d[i]-p;
			r=(d[i-1]-g)*s+2.0*c*b;
			d[i]=g+(p=s*r);
			g=c*r-b;
			for (k=1;k<=n;k++){
			    f=coordinates[k-1][i];
			    coordinates[k-1][i]=s*coordinates[k-1][i-1]+c*f;
			    coordinates[k-1][i-1]=c*coordinates[k-1][i-1]-s*f;
			}
		    }
		    if (r == 0.0 && i >= l) continue;
		    d[l-1] -= p;
		    e[l-1]=g;
		    e[m-1]=0.0;
		}
	    }while(m != l);
	}
	if (AGCT.Debug) System.out.print("ok.\n");
	cc.end();
    }

    public void eigenSort(double d[]){
	if (AGCT.Debug) System.out.print("Sorting eigenvalues of " + name + " by descending order... ");

	Matrix.checkSquare(this);

	int k,j,i, n = dimX;
	double p;

	for (i=0;i<=n-1;i++){
	    p=d[k=i];
	    for (j=i+1;j<=n-1;j++)
		if (d[j]>=p) p = d[k=j];
	    if (k != i) {
		d[k] = d[i];
		d[i] = p;
		for (j=0;j<=n-1;j++){
		    p = coordinates[j][i];
		    coordinates[j][i] = coordinates[j][k];
		    coordinates[j][k] = p;
		}
	    }
	}

	if (AGCT.Debug) System.out.print("ok.\n");
    }

    public void transpose(Matrix w){
	//puts the transpose in this
	//does not check for all other parameters
	if ( (dimX != w.dimY) || (dimY != w.dimX) )
	    Matrix.perror("Matrix.class :: Dimension mismatch to transpose matrice");

	int i,j;
	for (i=0;i<dimX;i++)
	    for (j=0;j<dimY;j++)
		coordinates[i][j] = w.coordinates[j][i];

	if (dimX == 1)
	    isVector = true;
	else
	    isVector = false;
	hasDefaultContent = w.hasDefaultContent;
	if (dimX == dimY)
	    isSquared = true;
	else
	    isSquared = false;
	isSymmetric = w.isSymmetric;
    }

    public void transposeSquare(){
	//if (AGCT.Debug) System.out.print("Transposing square matrix " + name + "... ");

	Matrix.checkSquare(this);

	int x,y, n = dimX;
	double val;
	for (x=0;x<n;x++){
	    for (y=x+1;y<n;y++){
		val = coordinates[x][y];
		coordinates[x][y] = coordinates[y][x];
		coordinates[y][x] = val;
	    }
	}
	//if (AGCT.Debug) System.out.print("ok.\n");
    }

    public boolean checkOrthonormality(Domain dd, int k){
	dd.myAGCT.myInformationFrame.setTextProgressBar("Checkings (" + name + ")");

	double val;
	int x1,x2,xxmax,i;
	boolean defective = false;
	int n = dimX, vm, nmade = 0, percent;
	
	Matrix.checkSquare(this);
	double [] p1 = new double[n];
	double [] p2 = new double[n];

	if ( (AGCT.Method_N == 0) || (AGCT.Method_N == 2) || (AGCT.Method_N == 3) ){
	    if (AGCT.Debug) System.out.print("Checking Orthonormality for system " + name + " (" + dimX + " x " + dimY + ") ... ");

	    xxmax = (k==-1) ? n : k;
	    vm = xxmax * xxmax;
    
	    for (x1=0;x1<xxmax;x1++){
		for (x2=x1;x2<xxmax;x2++){
		    percent = (int) (100.0 * ((double) nmade / (double) (vm)));
		    dd.myAGCT.myInformationFrame.setValueProgressBar(percent);

		    for (i=0;i<n;i++){
			p1[i] = coordinates[x1][i];
			p2[i] = coordinates[x2][i];
		    }
		    val = dot(p1, p2);

		    if (x1 == x2) val = 1.0 - val;
		    if ( (Math.abs(val)) > Precision_For_Eigensystems ){
			if (AGCT.Debug) System.out.print("\nOrthogonality problem on vectors " + x1 + " and " + x2 + ": (1-)innerprod = " + val + " (limit precision = " + Precision_For_Eigensystems + ")");

			dd.myAGCT.myInformationFrame.setValueProgressBar(0);
			dd.myAGCT.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);

			return false;
		    }  
		    nmade ++;
		}
	    }

	    if (AGCT.Debug) System.out.print("ok.\n");
	}else if (AGCT.Method_N == 1){
	    if (AGCT.Debug) System.out.print("Checking Normality for system " + name + " (" + dimX + " x " + dimY + ")... ");

	    xxmax = (k==-1) ? n : k;
	    vm = xxmax * xxmax;
	    for (x1=0;x1<xxmax;x1++){
		for (x2=x1;x2<xxmax;x2++){
		    percent = (int) (100.0 * ((double) nmade / (double) (vm)));
		    dd.myAGCT.myInformationFrame.setValueProgressBar(percent);

		    for (i=0;i<n;i++){
			p1[i] = coordinates[x1][i];
			p2[i] = coordinates[x2][i];
		    }
		    val = dot(p1, p2);

		    if (x1 == x2){
			val = 1-val;
			if ( (Math.abs(val)) > Precision_For_Eigensystems ){
			    if (AGCT.Debug) System.out.print("\nOrthogonality problem on vectors " + x1 + " and " + x2 + ": (1-)innerprod = " + val + " (limit precision = " + Precision_For_Eigensystems + ")");
			    dd.myAGCT.myInformationFrame.setValueProgressBar(0);
			    dd.myAGCT.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);

			    return false;
			}
		    }
		    else{
			if ( (Math.abs(val)) > Precision_For_Eigensystems ){
			    defective = true;
			}
		    }
		    nmade ++;
		}
	    }
	    if (AGCT.Debug){
		System.out.print("ok");
		if (defective == true)
		    System.out.print(" (but " + name + " is defective)\n");
		else
		    System.out.print(" (and " + name + " is not defective)\n");
	    }
	}

	dd.myAGCT.myInformationFrame.setValueProgressBar(0);
	dd.myAGCT.myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);

	return true;
    }

    public synchronized void ludcmp(int [] indx, double [] d, boolean [] singular){
	if ( (dimX != dimY) || (dimX != indx.length) )
	    Matrix.perror("Dimension mismatch");
	if ( (d.length != 1) || (singular.length != 1) )
	    Matrix.perror("Dimension of d is not 1");
	int i, imax = -1, j, k, n = dimX;
	double big, dum, sum, temp;
	double [] vv = new double [dimX];
	d[0] = 1.0;
	
	singular[0] = false;

	for (i=1;i<=n;i++){
	    big = 0.0;
	    for (j=1;j<=n;j++)
		if ((temp=Math.abs(coordinates[i-1][j-1])) > big)
		    big = temp;
	    if (big == 0.0){
		singular[0] = true;
		return;
	    }
	    vv[i-1] = 1.0/big;
	}
	for (j=1;j<=n;j++){
	    for (i=1;i<j;i++){
		sum = coordinates[i-1][j-1];
		for (k=1;k<i;k++)
		    sum -= coordinates[i-1][k-1] * coordinates[k-1][j-1];
		coordinates[i-1][j-1] = sum;
	    }
	    big = 0.0;
	    for (i=j;i<=n;i++){
		sum = coordinates[i-1][j-1];
		for (k=1;k<j;k++)
		    sum -= coordinates[i-1][k-1] * coordinates[k-1][j-1];
		coordinates[i-1][j-1] = sum;
		if ( (dum=vv[i-1]*Math.abs(sum)) >= big ){
		    big = dum;
		    imax = i;
		}
	    }
	    if (j != imax){
		for (k=1;k<=n;k++){
		    dum = coordinates[imax-1][k-1];
		    coordinates[imax-1][k-1] = coordinates[j-1][k-1];
		    coordinates[j-1][k-1] = dum;
		}
		d[0] = -d[0];
		vv[imax-1] = vv[j-1];
	    }
	    indx[j-1] = imax;
	    if (coordinates[j-1][j-1] == 0.0)
		coordinates[j-1][j-1] = Matrix.Tiny;
	    if (j!= n){
		dum = 1.0 / coordinates[j-1][j-1];
		for (i=j+1;i<=n;i++)
		    coordinates[i-1][j-1] *= dum;
	    }
	}
    }

    public synchronized void determinant(int [] indx, double [] d, boolean [] singular, double [] determine){
	int j;
	double det = 0.0;
	ludcmp(indx, d, singular);
	if (singular[0] == false){
	    det = d[0];
	    for (j=0;j<dimX;j++)
		det *= coordinates[j][j];
	}
	determine[0] = det;
    }


    public void completePositiveFactorization(Domain dd, Matrix DS){
	int i, j, K = dimX, N = dimY, r, s, niter=0, iq;
	double frob = 0.0, lastfrob = 0.0, deltafrob, num, den, sum;
	boolean stop = false;
	  
	Matrix.checkSquare(DS);

	AGCTCounter cc = new AGCTCounter(dd.myAGCT.myInformationFrame, "Pos. fact. of " + name, 100 );

	if (DS.dimX != N)
	    Matrix.perror("Dimension mismatch for CP factorization.");

	if (AGCT.Debug) System.out.print("Computing the complete positive factorization of system " + name + " (" + dimX + " x " + dimY + ") ... ");

	do{
	    for (r=0;r<K;r++){
		for (s=0;s<N;s++){
		    num = 0.0;
		    for (i=0;i<N;i++){
			if (i!=s)
			    num += (coordinates[r][i] * DS.coordinates[s][i]);
		    }
		    num *= coordinates[r][s];

		    den = 0.0;
		    for (j=0;j<K;j++){
			for (i=0;i<N;i++){
			    if (i!=s)
				den += ( coordinates[j][s] * coordinates[j][i] * coordinates[r][i] );
			}
		    }

		    if ( den != 0.0 ){
			coordinates[r][s] = ( num / den ); 
		    }
		    else{
			if (num == 0.0)
			    coordinates[r][s] = 0.0;
			else{
			    Matrix.perror("CP-factorization: num = " + num + ", den = " + den);
			}
		    }
		}
	    }
	    
	    frob = 0.0;
	    for (i=0;i<N;i++){
		for (j=0;j<N;j++){
		    if (j!=i){
			sum = 0.0;
			for (r=0;r<K;r++)
			    sum += (coordinates[r][i] * coordinates[r][j]);
			frob += ( (DS.coordinates[i][j] - sum) * (DS.coordinates[i][j] - sum) );
		    }
		}
	    }

	    if (niter == 0)
		lastfrob = frob;
	    else{
		deltafrob = lastfrob - frob;
		if (deltafrob < 0.0)
		    deltafrob = -deltafrob;
		
		iq = (int) (100.0 * Precision_For_Doubly_Stochastic_Approximation/deltafrob);
		if (iq > 100)
		    iq = 100;
		cc.setCounter(iq);

		if (deltafrob < Precision_For_Doubly_Stochastic_Approximation)
		    stop = true;
		
		lastfrob = frob;
	    }
	    
	    niter++;
	}while (stop == false);
	
	cc.end();
	if (AGCT.Debug) System.out.print("ok.\n");
    }

    public void normalize(){
	int i,j, nlig = dimX, ncol = dimY;
	double sum=0.0;
	
	for (j=0;j<ncol;j++){
	    sum=0.0;
	    for (i=0;i<nlig;i++){
		if (coordinates[i][j] < Precision_For_Eigensystems)
		    coordinates[i][j] = 0.0;
		sum+=coordinates[i][j];
	    }

	    if (sum == 0.0)
		Matrix.perror("The sum, " + sum + ", is not strictly positive");
	    
	    for (i=0;i<nlig;i++)
		coordinates[i][j] /= sum;
	}
    }

    public double logdet(){
	double sum = 0.0;
	if (dimX != dimY)
	    Matrix.perror("Dimension mismatch");
	for (int i=0;i<dimX;i++)
	    sum += Math.log(choleski_decomposition.coordinates[i][i]);
	return 2.0*sum;
    }

    public void elsolve(double [] b, double [] y){
	int i, j, n = dimX;
	double sum;
	if ( (b.length != dimX) || (y.length != dimX) || (dimX != dimY) )
	    Matrix.perror("Dimension mismatch");
	for (i=0;i<n;i++){
	    for (sum=b[i], j=0; j<i ; j++)
		sum -= choleski_decomposition.coordinates[i][j] * y[j];
	    y[i] = sum / choleski_decomposition.coordinates[i][i];
	}
    }

    public void toP(Domain d, Clustering clu, double valP, double valSKK){
	if (d.domainGenes == null)
	    Matrix.perror("No Genes !");
	
	AGCTCounter cc = new AGCTCounter(d.myAGCT.myInformationFrame, "computing " + name, dimX *  dimY);
	if (AGCT.Debug) System.out.print("Generating matrix P (" + dimX + " x " + dimY + ") ... ");
	
	int i, j, k;
	Gene gg;
	
	double K = (double) dimX;
	double N = (double) d.numberSelectedGenes;
	int iN = d.numberSelectedGenes;
	if ( (dimX != iN) || (dimY != iN) )
	    Matrix.perror("Dimension mismatch");

	for (i=0;i<dimX;i++)
	    for (j=0;j<dimY;j++)
		coordinates[i][j] = 0.0;
	
	for (i=0;i<dimX;i++){
	    for (j=0;j<dimY;j++){
		if (i!=j)
		    coordinates[i][j] = - Distortion.Bregman(clu.getRightPoint(i, clu.myClusteringAlgorithm.before),
							     clu.getRightPoint(j, clu.myClusteringAlgorithm.before),
							     4,
							     valP);
		else
		    coordinates[i][j] = valSKK;
		cc.increment();
	    }
	}
	cc.end();
	if (AGCT.Debug) System.out.println("ok.");
    }


    /// INTRODUCTION LANCZOS

    static class Entry implements Comparable<Entry> {// 降順にソート

        int i;
        double d;

        Entry(int i, double d) {
            this.d = d;
            this.i = i;
        }

        @Override
        public int compareTo(Entry o) {
            return (int) Math.signum(o.d - d);
        }
    }



    public static void getEigenSystem(Domain ddd, Matrix mat, double[] eigenValues) {// mat.row(i)が、対応するeigenVectorになる
        if(!mat.isSymmetric())
	    Matrix.perror("Matrix not symmetric");
        int n = mat.dimY;
        double[] beta = new double[n];
        int iter = Math.min(AGCT.MAX_COL_OF_EIGENVECTOR, mat.dimX);	

	CHECK_NaN("MyEigenSystem :: Before Lanczos", mat);

        lanczos(ddd, mat, eigenValues, beta, iter);
        System.gc();

	CHECK_NaN("MyEigenSystem :: Before calcEigenVectorOfTridiagonalizedMatrix", mat);

        calcEigenVectorOfTridiagonalizedMatrix(ddd, mat, eigenValues, beta, iter);
        System.gc();

	CHECK_NaN("MyEigenSystem :: Before lanczosEigenSort", mat);

        lanczosEigenSort(mat, eigenValues, iter);
        System.gc();

	CHECK_NaN("MyEigenSystem :: Before transposeSquare", mat);

	mat.transposeSquare();
        //transpose(mat);
        System.gc();
    }


    public Matrix mul(Matrix mat) {
        assert dimY == mat.dimX;
        int m = dimX, n = dimY, l = mat.dimY;
        Matrix res = new Matrix("Dummy", m, l);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < l; j++) {
                double value = 0;
                for (int k = 0; k < n; k++) {
                    value += get(i, k) * mat.get(k, j);
                }
                res.set(i, j, value);
            }
        }
        return res;
    }

    public MyVector getRow(int x) {
        DenseVector res = new DenseVector(dimY);
        for (int i = 0; i < dimY; i++) {
            res.set(i, get(x, i));
        }
        return res;
    }

    public MyVector mul(MyVector v) {
        MyVector res = new DenseVector(dimX);
        for (int i = 0; i < dimX; i++) {
            res.set(i, v.dot(getRow(i)));
        }
        return res;
    }

    private static MyVector lanczos(Domain ddd, Matrix A, double[] alpha, double[] beta,
                                  int iter) {
        int n = A.dimX;
        Matrix B = new Matrix(A.name,A.dimX,A.dimY);
	B.copyOnThis(A);

	AGCTCounter cc = new AGCTCounter(ddd.myAGCT.myInformationFrame, "Lanczos " + A.name + " with N_Cols = " + AGCT.MAX_COL_OF_EIGENVECTOR, iter);

        double[] random = new double[n];
        Random rand = new Random(14214987);
        double norm = 0;
        for (int i = 0; i < n; i++) {
            random[i] = rand.nextDouble();
            norm += random[i] * random[i];
        }
        norm = Math.sqrt(norm);
        MyVector v0 = new SparseVector(n), v1 = new SparseVector(n);
        MyVector res = new DenseVector(n);
        for (int i = 0; i < n; i++) {
            v1.set(i, random[i] / norm);
            res.set(i, random[i] / norm);
        }
        MyVector w;
        beta[0] = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A.set(i, j, 0);
            }
        }
        for (int j = 0; j < iter; j++) {
            assert eq(v1.norm(), 1);
	    cc.increment();

            for (int i = 0; i < n; i++) {
                A.set(i, j, v1.get(i));
            }
            w = B.mul(v1).sub(v0.mul(beta[j]));
            alpha[j] = w.dot(v1);
            if (j == n - 1)
                break;
            w = w.sub(v1.mul(alpha[j]));
            beta[j + 1] = w.norm();
            if (beta[j + 1] == 0)
                break;
            v0 = v1;
            v1 = w.div(beta[j + 1]);
        }
	cc.end();
        return res;
    }

    private static void calcEigenVectorOfTridiagonalizedMatrix(Domain ddd, Matrix A,
                                                               double[] alpha, double[] beta, int row) {
        int n = A.dimX;
	AGCTCounter cc = new AGCTCounter(ddd.myAGCT.myInformationFrame, "Computing eigenvectors " + A.name, row);

	if (A.dimX != A.dimY)
	    Matrix.perror("Non squared matrix");

        int m, L, iter, i, k;

        double s, r, p, g, f, dd, c, b, sf;
        for (i = 2; i <= n; i++)
            beta[i - 2] = beta[i - 1];
        beta[n - 1] = 0.0;

        for (L = 1; L <= row; L++) {
            cc.increment();

            iter = 0;
            do {
                for (m = L; m <= row - 1 - 1; m++) {
                    dd = Math.abs(alpha[m - 1]) + Math.abs(alpha[m]);
                    sf = (double) (Math.abs(beta[m - 1]) + dd);
                    if (sf == dd)
                        break;
                }
                if (m != L) {
                    if (iter++ == Numerical_Iterations_Max)
                        Matrix.perror("Too many iterations in tqli");
                    g = (alpha[L] - alpha[L - 1]) / (2.0 * beta[L - 1]);
                    r = pythag(g, 1.0);
                    g = alpha[m - 1] - alpha[L - 1] + beta[L - 1] / (g + sign(r, g));
                    s = c = 1.0;
                    p = 0.0;
                    for (i = m - 1; i >= L; i--) {
                        f = s * beta[i - 1];
                        b = c * beta[i - 1];
                        beta[i] = (r = pythag(f, g));
                        if (r == 0.0) {
                            alpha[i] -= p;
                            beta[m - 1] = 0.0;
                            break;
                        }
                        s = f / r;
                        c = g / r;
                        g = alpha[i] - p;
                        r = (alpha[i - 1] - g) * s + 2.0 * c * b;
                        alpha[i] = g + (p = s * r);
                        g = c * r - b;
                        for (k = 1; k <= n; k++) {// dimXまでが必須
                            f = A.get(k - 1, i);
                            A.set(k - 1, i, s * A.get(k - 1, i - 1) + c * f);
                            A.set(k - 1, i - 1, c * A.get(k - 1, i - 1) - s * f);
                        }
                    }
                    if (r == 0.0 && i >= L)
                        continue;
                    alpha[L - 1] -= p;
                    beta[L - 1] = g;
                    beta[m - 1] = 0.0;
                }
            } while (m != L);
        }
        if (AGCT.Debug)
            System.out.print("ok.\n");
        cc.end();
    }


    private static final double EPS = 1e-3;

    private static boolean eq(double a, double b) {
        //		if (Math.abs(a - b) < EPS)
        //			return true;
        if (Math.abs(1 - a / b) < EPS)
            return true;
        if (Math.abs(1 - b / a) < EPS)
            return true;
        return false;
    }

    private static void lanczosEigenSort(Matrix mat, double[] eigen, int iter) {//先頭からiterのみ計算する
 	
	int n = mat.dimX;
	if (mat.dimX != mat.dimY)
	    Matrix.perror("Non squared matrix");


        Entry[] es = new Entry[iter];
        for (int i = 0; i < iter; i++) {
            es[i] = new Entry(i, eigen[i]);
        }
        Arrays.sort(es);
        ArrayList<Integer> same = new ArrayList<Integer>();
        ArrayList<Integer> diff = new ArrayList<Integer>();
        for (int i = 0; i < iter; ) {
            diff.add(es[i].i);
            int j = i + 1;
            while (j < iter && eq(eigen[es[i].i], eigen[es[j].i])) {
                same.add(es[j].i);
                j++;
            }
            i = j;
        }

	Matrix dupM = new Matrix(mat.name, mat.dimX, mat.dimY);
	dupM.copyOnThis(mat);

        double[] dupE = eigen.clone();
        for (int j = 0; j < diff.size(); j++) {
            for (int i = 0; i < n; i++) {
                mat.set(i, j, dupM.get(i, diff.get(j)));
            }
            eigen[j] = dupE[diff.get(j)];
        }
        for (int j = 0; j < same.size(); j++) {
            for (int i = 0; i < n; i++) {
                mat.set(i, j + diff.size(), dupM.get(i, same.get(j)));
            }
            eigen[j + diff.size()] = dupE[same.get(j)];
        }
    }

    // TO REMOVE WHEN TESTS OKAY
    // TEST: eigensystem of Natalia

    public static void CHECK_NaN(String misc, Matrix mat){
	int i, j;
	for (i=0;i<mat.dimX;i++){
	    for (j=0;j<mat.dimY;j++){
		if (Statistics.isNaN(mat.get(i, j)))
		    Matrix.perror("coordinate (" + i + ", " + j + ") in matrix " + mat.name + " is NaN (Step " + misc +")");
	    }
	}
    }

    public static void Natalia_getEigenSystem(Domain ddd, Matrix mat, double[] eigenValues) {// mat.row(i)が、対応するeigenVectorになる

	if (!AGCT.SUBMISSION)
	    System.out.println(" ** Warning : Eigensystem Using Natalia's parameters ** ");
	//System.out.println("Mat 1 = + " + mat.affiche(100));

	if(!mat.isSymmetric())
	    Matrix.perror("Matrix not symmetric");
	
        int n = mat.dimY;
        double[] d = new double[n];

	CHECK_NaN("Before tridiagonalize", mat);

        Natalia_Tridiagonalize(ddd, mat, eigenValues, d);

	CHECK_NaN("Before eigenCalc", mat);

	//System.out.println("Mat 2 = + " + mat.affiche(100));

        Natalia_calcEigenVectorOfTridiagonalizedMatrix(ddd, mat, eigenValues, d);

	//System.out.println("Mat 3 = + " + mat.affiche(100));

	CHECK_NaN("Before eigenSort", mat);

        Natalia_eigenSort(mat, eigenValues);

	//System.out.println("Mat 4 = + " + mat.affiche(100));

	CHECK_NaN("Before transpose", mat);

        Natalia_transpose(mat);

	CHECK_NaN("After transpose", mat);
    }

    static void Natalia_transpose(Matrix mat) {
        assert mat.dimX == mat.dimY;
        int n = mat.dimX;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double d = mat.get(i, j);
                mat.set(i, j, mat.get(j, i));
                mat.set(j, i, d);
            }
        }
    }

    public static void Natalia_Tridiagonalize(Domain ddd, Matrix mat, double d[], double e[]) {
        if (Natalia_useLanczos(d.length)) {
            Natalia_lanczos(ddd, mat, d, e);
            return;
        }

	int cp = 10, ccur = 0, ctot = mat.dimX - 2 + mat.dimX - 1;

        int dimX = mat.dimX;
        int dimY = mat.dimY;

	AGCTCounter cc = new AGCTCounter(ddd.myAGCT.myInformationFrame, "Householder (" + mat.name + ")", ctot);

        if (AGCT.Debug) {
            System.out.print("Householder tridiagonalization on Matrix " + " ("
                    + dimX + " x " + dimY + ") (start)... ");
        }

	if(!mat.isSymmetric())
	    Matrix.perror("Matrix not symmetric");

        int l, k, j, i, n = dimX;
        double scale, hh, h, g, f;

        for (i = n; i >= 2; i--) {
	    cc.increment();

	    /*if (ccur%(ctot/cp) == 0)
		System.out.print( (ccur/(ctot/cp))*cp + " % ");
		ccur++;*/

            l = i - 1;
            h = scale = 0.0;
            if (l > 1) {
                for (k = 1; k <= l; k++) {
                    scale += Math.abs(mat.get(i - 1, k - 1));
                }
                if (scale == 0.0) {
                    e[i - 1] = mat.get(i - 1, l - 1);
                } else {
                    for (k = 1; k <= l; k++) {
                        mat.set(i - 1, k - 1, mat.get(i - 1, k - 1) / scale); // coordinates[i
                        // -
                        // 1][k
                        // - 1]
                        // /=
                        // scale;
                        h += mat.get(i - 1, k - 1) * mat.get(i - 1, k - 1);
                    }
                    f = mat.get(i - 1, l - 1);
                    g = (f >= 0.0 ? -Math.sqrt(h) : Math.sqrt(h));
                    e[i - 1] = scale * g;
                    h -= f * g;
                    mat.set(i - 1, l - 1, f - g);
                    f = 0.0;
                    for (j = 1; j <= l; j++) {
                        mat.set(j - 1, i - 1, mat.get(i - 1, j - 1) / h);
                        g = 0.0;
                        for (k = 1; k <= j; k++) {
                            g += mat.get(j - 1, k - 1) * mat.get(i - 1, k - 1);
                        }
                        for (k = j + 1; k <= l; k++) {
                            g += mat.get(k - 1, j - 1) * mat.get(i - 1, k - 1);
                        }
                        e[j - 1] = g / h;
                        f += e[j - 1] * mat.get(i - 1, j - 1);
                    }
                    hh = f / (h + h);
                    for (j = 1; j <= l; j++) {
                        f = mat.get(i - 1, j - 1);
                        e[j - 1] = g = e[j - 1] - hh * f;
                        for (k = 1; k <= j; k++) {
                            mat.set(j - 1, k - 1, mat.get(j - 1, k - 1)
                                    - (f * e[k - 1] + g * mat.get(i - 1, k - 1)));
                        }
                        // coordinates[j - 1][k - 1] -= (f * e[k - 1] + g *
                        // get(i - 1,k - 1));
                    }
                }
            } else {
                e[i - 1] = mat.get(i - 1, l - 1);
            }
            d[i - 1] = h;
        }

        d[0] = 0.0;
        e[0] = 0.0;
        for (i = 1; i <= n; i++) {
	    cc.increment();

	    /*if (ccur%(ctot/cp) == 0)
		System.out.print( (ccur/(ctot/cp))*cp + " % ");
		ccur++;*/

            l = i - 1;
            if (d[i - 1] != 0.0) {
                for (j = 1; j <= l; j++) {
                    g = 0.0;
                    for (k = 1; k <= l; k++) {
                        g += mat.get(i - 1, k - 1) * mat.get(k - 1, j - 1);
                    }
                    for (k = 1; k <= l; k++) // coordinates[k - 1][j - 1] -= g * get(k - 1,i - 1);
                    {
                        mat.set(k - 1, j - 1, mat.get(k - 1, j - 1) - g
                                * mat.get(k - 1, i - 1));
                    }
                }
            }
            d[i - 1] = mat.get(i - 1, i - 1);
            mat.set(i - 1, i - 1, 1.0);
            for (j = 1; j <= l; j++) {
                mat.set(j - 1, i - 1, 0);
                mat.set(i - 1, j - 1, 0.0);
            }
        }

        if (AGCT.Debug) {
            System.out.print("ok.\n");
        }

	cc.end();
    }


    private static boolean Natalia_useLanczos(int dim) {
        return AGCT.Natalia_lanczos && dim > AGCT.Natalia_MAX_COL_OF_EIGENVECTOR;
    }

    public static void Natalia_lanczos(Domain ddd, Matrix A, double[] alpha, double[] beta) {
        int dimX = A.dimX;
        int dimY = A.dimY;

        assert A.dimX == A.dimY;

	AGCTCounter cc = new AGCTCounter(ddd.myAGCT.myInformationFrame, "Lanczos " + A.name + " with N_Cols = " + AGCT.Natalia_MAX_COL_OF_EIGENVECTOR, A.dimY);

       if(!A.isSymmetric())
	    Matrix.perror("Matrix not symmetric");

       Matrix B = new Matrix(A.name,A.dimX,A.dimY);
       B.copyOnThis(A);

       double[] random = new double[A.dimX];
       Random rand = new Random(0);
       double norm = 0;
       for (int i = 0; i < A.dimX; i++) {
	   random[i] = rand.nextDouble();
	   norm += random[i] * random[i];
       }
       norm = Math.sqrt(norm);
       MyVector v0 = new SparseVector(A.dimX), v1 = new SparseVector(A.dimX);
       for (int i = 0; i < dimX; i++) {
	   v1.set(i, random[i] / norm);
       }
       MyVector w;
       beta[0] = 0;
       for (int i = 0; i < A.dimX; i++) {
	   for (int j = 0; j < A.dimY; j++) {
	       A.set(i, j, 0);
	   }
       }
       for (int j = 0; j < A.dimY; j++) {
	   assert eq(v1.norm(), 1);
	   cc.increment();

	   for (int i = 0; i < A.dimX; i++) {
	       A.set(i, j, v1.get(i));
	   }
	   w = B.mul(v1).sub(v0.mul(beta[j]));
	   alpha[j] = w.dot(v1);
	   if (j == A.dimY - 1) {
	       break;
	   }
	   w = w.sub(v1.mul(alpha[j]));
	   beta[j + 1] = w.norm();
	   if (beta[j + 1] == 0) {
	       break;
	   }
	   v0 = v1;
	   v1 = w.div(beta[j + 1]);
       }
       cc.end();
    }

    public static void Natalia_calcEigenVectorOfTridiagonalizedMatrix(Domain ddd, Matrix mat,
                                                              double[] alpha, double[] beta) {
        int dimX = mat.dimX;

        assert mat.dimX == mat.dimY;

        int m, L, iter, i, k, n = mat.dimX;
        double s, r, p, g, f, dd, c, b, sf;
        for (i = 2; i <= n; i++) {
            beta[i - 2] = beta[i - 1];
        }
        beta[n - 1] = 0.0;

	AGCTCounter cc = new AGCTCounter(ddd.myAGCT.myInformationFrame, "Computing eigenvectors " + mat.name, n);

        for (L = 1; L <= n; L++) {
            cc.increment();
            iter = 0;
            do {
                for (m = L; m <= n - 1; m++) {
                    dd = Math.abs(alpha[m - 1]) + Math.abs(alpha[m]);
                    sf = (double) (Math.abs(beta[m - 1]) + dd);
                    if (sf == dd) {
                        break;
                    }
                }
                if (m != L) {
                    if (iter++ == Natalia_Numerical_Iterations_Max) {
                        Matrix.perror("Too many iterations in tqli");
                    }
                    g = (alpha[L] - alpha[L - 1]) / (2.0 * beta[L - 1]);
                    r = pythag(g, 1.0);
                    g = alpha[m - 1] - alpha[L - 1] + beta[L - 1] / (g + sign(r, g));
                    s = c = 1.0;
                    p = 0.0;
                    for (i = m - 1; i >= L; i--) {
                        f = s * beta[i - 1];
                        b = c * beta[i - 1];
                        beta[i] = (r = pythag(f, g));
                        if (r == 0.0) {
                            alpha[i] -= p;
                            beta[m - 1] = 0.0;
                            break;
                        }
                        s = f / r;
                        c = g / r;
                        g = alpha[i] - p;
                        r = (alpha[i - 1] - g) * s + 2.0 * c * b;
                        alpha[i] = g + (p = s * r);
                        g = c * r - b;
                        for (k = 1; k <= n; k++) {
                            f = mat.get(k - 1, i);
                            mat.set(k - 1, i, s * mat.get(k - 1, i - 1) + c * f);
                            mat.set(k - 1, i - 1, c * mat.get(k - 1, i - 1) - s * f);
                        }
                    }
                    if (r == 0.0 && i >= L) {
                        continue;
                    }
                    alpha[L - 1] -= p;
                    beta[L - 1] = g;
                    beta[m - 1] = 0.0;
                }
            } while (m != L);
        }
        if (AGCT.Debug) {
            System.out.print("ok.\n");
        }
	cc.end();
    }


    public static void Natalia_eigenSort(Matrix mat, double[] eigen) {
        if (Natalia_useLanczos(eigen.length)) {
            Natalia_lnczosEigenSort(mat, eigen);
            return;
        }

        if (AGCT.Debug) {
            System.out.print("Sorting eigenvalues" + " by descending order... ");
        }

        assert mat.dimX == mat.dimY;

        int k, j, i, n = mat.dimX;
        double p;
        for (i = 0; i <= n - 1; i++) {
            p = eigen[k = i];
            for (j = i + 1; j <= n - 1; j++) {
                if (eigen[j] >= p) {
                    p = eigen[k = j];
                }
            }
            if (k != i) {
                eigen[k] = eigen[i];
                eigen[i] = p;
                for (j = 0; j <= n - 1; j++) {
                    p = mat.get(j, i);
                    mat.set(j, i, mat.get(j, k));
                    mat.set(j, k, p);
                }
            }
        }
        if (AGCT.Debug) {
            System.out.print("ok.\n");
        }
    }


    private static void Natalia_lnczosEigenSort(Matrix mat, double[] eigen) {
        int n = mat.dimX;
        assert mat.dimX == mat.dimY;
        Entry[] es = new Entry[n];
        for (int i = 0; i < n; i++) {
            es[i] = new Entry(i, eigen[i]);
        }
        Arrays.sort(es);
        ArrayList<Integer> same = new ArrayList<Integer>();
        ArrayList<Integer> diff = new ArrayList<Integer>();
        for (int i = 0; i < n; ) {
            diff.add(es[i].i);
            int j = i + 1;
            while (j < n && eq(eigen[es[i].i], eigen[es[j].i])) {
                same.add(es[j].i);
                j++;
            }
            i = j;
        }
        Matrix dupM = new Matrix(mat.name,mat.dimX,mat.dimY);
	dupM.copyOnThis(mat);

        double[] dupE = eigen.clone();
        for (int j = 0; j < diff.size(); j++) {
            for (int i = 0; i < n; i++) {
                mat.set(i, j, dupM.get(i, diff.get(j)));
            }
            eigen[j] = dupE[diff.get(j)];
        }
        for (int j = 0; j < same.size(); j++) {
            for (int i = 0; i < n; i++) {
                mat.set(i, j + diff.size(), dupM.get(i, same.get(j)));
            }
            eigen[j + diff.size()] = dupE[same.get(j)];
        }
    }
}

