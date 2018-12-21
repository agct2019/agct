import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import javax.swing.event.*;
import java.util.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;


public class AGCT extends javax.swing.JApplet implements Runnable, Debuggable, ActionListener{

    /***
     * Parameters for submission of code
     ***/

    public static boolean VERBOSE = false;

    public static boolean SUBMISSION = false;//true;

    public static String SUBMISSION_TEXT = "\n//  Supplementary Java code supporting the submission  of the paper:\n// \n"
			   +"//  \"AGCT: A Geometric Clustering Tool for robustly unravelling the inner cluster structures of gene expressions\",\n"
			   +"//  by R. Nock, N. Polouliakh, K. Oka, F. Nielsen and H. Kitano\n// \n"
	+"//  This version of AGCT goes with our Submission\n//  and it may be used with the purpose of the Review process ONLY\n// \n//  *********************\n//  * Do NOT distribute *\n//  *********************\n// \n//  CONSULT THE README FILE\n// \n//  Contact : richard.nock@nicta.com.au or nata@csl.sony.co.jp\n";



    /***
     * Show all options for filtering the similarity matrix with another data involving "distances"
     ***/
    public static boolean SHOW_FILTERING = false;

    /*** 
     * VARIABLE that tells whether the computations are using Natalia's AGCT algorithms, so that we math the output
     * of Natalia's version of AGCT
     * (apparently, the choice of random variables is made differently in our version so the results may differ,
     * furthermore, the parameters in Lanczos are different
     ***/
    public static boolean USE_NATALIA = true;

    /*****
     * Variable to try to gain memory : set to false to prevent bad surprises
     *****/

    public static boolean SAVE_MEMORY = true;
    public static boolean SAVE_MEMORY_NO_ANN = true;

    public static boolean No_GeneSelectionAvailable = true;
    // prevents Gene selection, and avoids a StackOverflow Error on the checkboxes when clicking on a different tab than the Gene selection tab
    
    /**********************************************************************************
     * Static stuff
     *****/

    public static boolean NORMALIZE_UNIT_PCA = true;
    //if true, norms are normalized in PCA so that the largest has unit norm
    //(prevents all-very small or all-very large norms)

    public static int Max_K_Dimension_Estimation = 15;//XXX 50
    public static int Min_K_Dimension_Estimation = 2;
    public static int K_Dimension_Estimation = 10;
    //must be < Maximal_Number_Of_Neighbors_Dimension_Estimation

    public static double MIN_DIMENSION = 0.5;

    public static boolean COMPUTE_DIMENSION = true;
    // allows to compute dimension estimators

    public static boolean LOCAL_DIMENSION_COMPUTED = false;

    public static String Token_Ordered_List_Names_Begin = "@Ordered_List_Names_Begin",
	Token_Ordered_List_Names_End = "@Ordered_List_Names_End",
	Token_Matrix_W_Begin = "@Matrix_W_Begin",
	Token_Matrix_W_End = "@Matrix_W_End",
	Token_Manifold_Eigenvalues_Begin = "@Manifold_Eigenvalues_Begin",
	Token_Manifold_Eigenvalues_End = "@Manifold_Eigenvalues_End",
	Token_Matrix_M_Begin = "@Matrix_M_Begin",
	Token_Matrix_M_End = "@Matrix_M_End",
	Token_Matrix_Just_W = "@Just_W";

    public static String allTokens_Manifold [] = {Token_Ordered_List_Names_Begin,
						  Token_Ordered_List_Names_End,
						  Token_Matrix_W_Begin,
						  Token_Matrix_W_End,
						  Token_Manifold_Eigenvalues_Begin,
						  Token_Manifold_Eigenvalues_End,
						  Token_Matrix_M_Begin,
						  Token_Matrix_M_End,
						  Scenario.allKeywords[35],
						  Token_Matrix_Just_W};

    public static boolean Loading_From_Scenario_Manifold_Begin, Loading_From_Scenario_Manifold_End;
    public static boolean Loading_From_Scenario_Dimension_Begin, Loading_From_Scenario_Dimension_End;

    public static boolean MULTIPLE_ROWS_FOR_A_SAME_GENE_CONSIDERED_DIFFERENT_GENES = false;
    public static boolean MULTIPLE_ROWS_FOR_A_SAME_GENE_DELETED_EXCEPT_FIRST_ONE = true;
    public static boolean CHECK_TRIANGLE_CONSISTENCY = true;

    public static boolean HAAR_WAVELETS_REMOVE_CONSTANT = true;
    //we remove the product wih the constant function

    public static Random RG;
    //global random number generator

    public static double Very_Small = 10E-6;
    
    public static int DEFAULT_Max_Number_Genes = 10;

    public static final String MS = "MS",
	MF = "MF",
	MW = "MW",
	MN = "MN",
	ME = "ME",
	MD = "MD";

    public static String [] Feature_Method = {"Slopes", 
					      "Haar_Wavelets", 
					      "Daubechies_D4_Wavelets_Interval", 
					      "Daubechies_D4_Wavelets", 
					      "Daubechies_D8_Wavelets", 
					      "Daubechies_D10_Wavelets", 
					      "Daubechies_D12_Wavelets", 
					      "Daubechies_D14_Wavelets", 
					      "Daubechies_D16_Wavelets", 
					      "Daubechies_D18_Wavelets",  
					      "Daubechies_D20_Wavelets"};

    public static String [] Feature_Selection_Method = {"None",
							"K_Smallest_Average_Determination_Coefficient"};

    public static String [] Prototype_Selection_Method = {"None",
							  "K_Means_Cluster_Prototype_Aggregation"};

    public static int Number_Of_Clusters, Number_Of_Neighbors = 10, WindowWidth = 800;

    public static double Bethe_Hessian_Factor = 2.0;

    //Lanczos

    public static int MAX_COL_OF_EIGENVECTOR = 50; // WAS 500 BEFORE THE COMPArIsONS WITH NATALIA;
    public static int MAX_MAX_COL_OF_EIGENVECTOR = 500;

    // Natalia's implementation of Lanczos, for comparison.

    public static boolean Natalia_lanczos = false;
    public static int Natalia_MAX_COL_OF_EIGENVECTOR = 50;

    /*********************************
     * Feature computation and selection
     *****/

    public static int MAX_WAVELET_COEFFICIENTS_SAVE = 1024;
    public static boolean COMPLETE_WITH_EMPTY_DATA = false;
    // IF TRUE, when saving wavelet coefficients, ensures that we have MAX_WAVELET_COEFFICIENTS_SAVE columns in saving file; 
    // may be useful for post-processing

    public static int Method_F = 0;
    // computes features
    // 0 : use slopes
    // 1 : use Haar wavelet coefficients
    // 2 : use Daubechies D4 wavelets on the interval (Cohen, Daubechies & Vial, 1993)
    // 3 : use Daubechies D4
    // 4 : use Daubechies D8
    // 5 : use Daubechies D10
    // 6 : use Daubechies D12
    // 7 : use Daubechies D14
    // 8 : use Daubechies D16
    // 9 : use Daubechies D18
    // 10 : use Daubechies D20

    public static int Method_FS = 0;
    // feature selection methods
    // 0 : no selection
    // 1 : uses average determination coefficient (repeatedly removes the variable with the largest value)

    public static int Method_PS = 0;
    // prototype selection methods
    // 0 : no selection
    // 1 : uses k-means (runs k-means, and then keeps the k prototypes that are closest to their cluster center)

    public static int Max_Number_Of_Features = -1;
    // maximum number of features authorized (useful e.g. when using wavelets with large number of coefficients)
    // if -1, uses DEFAULT

    public static int DEFAULT_Max_Number_Of_Features = 50;

    public static int Max_Number_Of_Prototypes = -1;
    // maximum number of prototypes authorized (useful e.g. when extremely large number of prototypes)
    // if -1, uses DEFAULT

    public static int DEFAULT_Max_Number_Of_Prototypes = 1000;

    public static int Number_Of_Wavelet_Stamps = -1;
    // number of time stamps used to compute the wavelet decomposition
    // must be a power of 2, or -1 (the program chooses at best)

    public static int DEFAULT_Number_Of_Wavelet_Stamps = 256;



    // NATURE BIO : 2 -> 1

    public static int Method_S = 1;
    // computes similarity
    // 0 : Heat Kernel
    // 1 : Cosine distance
    // 2 : Absolute Cosine distance (= absolute correlation)

    public static int Method_W = 1;
    // computes similarity matrix W
    // 0 : keep similarity as is
    // 1 : filter = use Symmetric NN (very sparse matrix)
    // 2 : filter = use Natural Neighbors (very sparse matrix)

    public static int Method_N = 1;
    // computes normalized matrix N
    // 0 : N = D^{-1/2} W D^{-1/2}
    // 1 : N = D^{-1} W
    // 2 : N = 2*Stochastic(W)
    // 3 : N = 2*Stochastic(DS Scaling(W))
    // 4 : N = (r^2 - 1) Id - rA + D --- Bethe Hessian (Assortative), Krzakala PNAS, with r = +sqrt(2k) and k = # neighbors
    // 5 : N = (r^2 - 1) Id - rA + D --- Bethe Hessian (Disassortative), Krzakala PNAS, with r = -sqrt(2k) and k = # neighbors
    // for 4 or 5, the similarity matrix W must be filtered via symmetric NN

    public static int Method_E = 1;
    // computes eigensystem
    // 0 : QR
    // 1 : Lanczos, faster but approx

    public static int Method_D = 1;
    // filter with distances
    // 0 : remove all wrt threshold
    // 1 : convolutes similarities with a Gaussian kernel

    //Doubly_Stochastic_Scaling:
    //When true, the matrix is first multiplied by a real T so as to be the closest to
    //a doubly stochastic matrix. This does not change eigenvectors, but greatly reduces
    //the doubly stochastic approximation time
    //When used as similarity matrix, results become more "flattened"

    public static boolean Sort_Depth, No_GeneSelection = No_GeneSelectionAvailable, Use_Shadow, Debug, Warning, Perspective, Only_Referenced, Only_ReferencedEdges, Referenced_Available;
    public static int Max_Reference = -1;

    public static int Show_Correlation = -1;
    //-1 -> no discrimination
    //0 -> show neg
    //1 -> show pos

    public static boolean Show_Only_Genes_With_Visible_Edge = false;
    //0 -> no discrimination
    //1 -> is triangulation computed, shows only genes with at least one visible edge

    public static double T_Heat_Kernel = 2.0;
    //T parameter for the heat kernel;

    public static double T_Heat_Kernel_Similarity_Neighbors = 0.1;
    // T_Heat_Kernel used to compute similarities on the *manifold* (2D)

    public static double Sparsify_Statistic = 0.0;
    // p*100 % of the smallest similarities in W are replaced by 0

    public static boolean Filter_Triangulation_With_Distances = false;
    public static boolean Filter_Similarities_With_Distances = false;
    // if true, then a distance file is used to filter triangulation edges / similarities between farthest entities -> 0
    
    public static boolean Correction_Profiles_By_Average = false;
    // if true, profiles are normalized to averaged growth with distances

    public static int Number_Of_Manifold_Components = 20;
    // number of components on the manifold to keep for each individual

    public static int Number_Of_Triangulation_Dimensions = 3;
    // number of dimensions to consider for triangulation

    public static int indexInTokens_Manifold(String s){
	String ref, t;
	int i;
	for (i=0;i<allTokens_Manifold.length;i++){
	    ref = allTokens_Manifold[i];
	    if (s.length() >= ref.length()){
		t = s.substring(0,ref.length());
		if (t.equals(ref))
		    return i;
	    }
	}
	return -1;
    }

    public static void executeManifold(String ex, AGCT ap){
	StringTokenizer t;
	String s1, s2;
	Vector v;
	int i, nc;
	double dum;
	int dumi;

	if (!ControlProcess.hasTrue("filtersSelected"))
	    Matrix.perror("AGCT.class :: genes not computed");

	//bOrdered_List_Names, bMatrix_W, bManifold_Eigenvalues, bMatrix_M

	//Names must be ordered in the same way as in the scenario; otherwise stops
	if (Scenario.bOrdered_List_Names){
	    t = new StringTokenizer(ex, Scenario.myToken);
	    if (t.countTokens() != ap.myDomain.numberSelectedGenes)
		Matrix.perror("AGCT.class :: bad number of tokens for gene list (" + t.countTokens() + " != " + ap.myDomain.selectedGeneNumberToGeneNumber.length + ")");
	    for (i=0;i<ap.myDomain.numberSelectedGenes;i++){
		s1 = ((Gene) ap.myDomain.domainGenes.elementAt(ap.myDomain.selectedGeneNumberToGeneNumber[i])).name;
		s2 = t.nextToken();

		if (!s1.equals(s2))
		    Matrix.perror("AGCT.class :: bad name correspondence for gene (" + s1 + " != " + s2 + ")");		    
	    }
	}else if (Scenario.bMatrix_W){
	    if (ap.currentIndexW == -1){
		ap.W = new Matrix("W", ap.myDomain.numberSelectedGenes, ap.myDomain.numberSelectedGenes);
		ap.CW = new Matrix("W_Original", ap.myDomain.numberSelectedGenes, ap.myDomain.numberSelectedGenes);
		ap.currentIndexW++;
	    }
	    t = new StringTokenizer(ex, Scenario.myToken);
	    if ( (Scenario.MATRIX_STORING_METHOD == Scenario.METHOD_UAG) && (t.countTokens() != ap.myDomain.numberSelectedGenes) )
		Matrix.perror("AGCT.class :: bad number of tokens for Matrix W (" + t.countTokens() + " != " + ap.myDomain.selectedGeneNumberToGeneNumber.length + ")");
	    if (Scenario.MATRIX_STORING_METHOD == Scenario.METHOD_UAG)
		for (i=0;i<ap.myDomain.numberSelectedGenes;i++){
		    dum = Double.parseDouble(t.nextToken());
		    ap.W.coordinates[ap.currentIndexW][i] = dum;
		    ap.CW.coordinates[ap.currentIndexW][i] = dum;
		}
	    else if (Scenario.MATRIX_STORING_METHOD == Scenario.METHOD_SONY)
		while(t.hasMoreTokens()){
		    dumi = Integer.parseInt(t.nextToken());
		    dum = Double.parseDouble(t.nextToken());
		    ap.W.coordinates[ap.currentIndexW][dumi] = dum;
		    ap.CW.coordinates[ap.currentIndexW][dumi] = dum;
		}
	    else
		Matrix.perror("This storing method, " + Scenario.MATRIX_STORING_METHOD + ", is not recognized in a Scenario");
	    ap.currentIndexW++;
	}else if (Scenario.bManifold_Eigenvalues){
	    t = new StringTokenizer(ex, Scenario.myToken);
	    if (t.countTokens() != ap.myDomain.numberSelectedGenes)
		Matrix.perror("AGCT.class :: bad number of tokens for manifold eigenvalues (" + t.countTokens() + " != " + ap.myDomain.selectedGeneNumberToGeneNumber.length + ")");
	    for (i=0;i<ap.myDomain.numberSelectedGenes;i++){
		if (i==0)
		    ap.manifold_eigenvalues = new double[ap.myDomain.numberSelectedGenes];
		ap.manifold_eigenvalues[i] = Double.parseDouble(t.nextToken());	    
	    }
	}else if (Scenario.bMatrix_M){
	    if (ap.currentIndexM == -1){
		ap.M = new Matrix("Manifold Eigenvectors", ap.myDomain.numberSelectedGenes, ap.myDomain.numberSelectedGenes);
		ap.generate_Neighbors();
		ControlProcess.put("neighborsComputed",true);
		ap.currentIndexM++;
	    }
	    t = new StringTokenizer(ex, Scenario.myToken);
	    if ( (Scenario.MATRIX_STORING_METHOD == Scenario.METHOD_UAG) && (t.countTokens() != ap.myDomain.numberSelectedGenes) )
		Matrix.perror("AGCT.class :: bad number of tokens for Matrix M (" + t.countTokens() + " != " + ap.myDomain.selectedGeneNumberToGeneNumber.length + ")");

	    if (Scenario.MATRIX_STORING_METHOD == Scenario.METHOD_UAG)
		for (i=0;i<ap.myDomain.numberSelectedGenes;i++){
		    ap.M.coordinates[ap.currentIndexM][i] = Double.parseDouble(t.nextToken());
		}
	    else if (Scenario.MATRIX_STORING_METHOD == Scenario.METHOD_SONY)
		while(t.hasMoreTokens()){
		    try{
			dumi = Integer.parseInt(t.nextToken());
			dum = Double.parseDouble(t.nextToken());
			ap.M.coordinates[ap.currentIndexM][dumi] = dum;
		    }catch(java.util.NoSuchElementException e){
			System.out.println("Matrix M ic cut in Scenario at index = " + ap.currentIndexM);
			System.out.println(ex);
			System.exit(0);
		    }
		}
	    else
		Matrix.perror("This storing method, " + Scenario.MATRIX_STORING_METHOD + ", is not recognized in a Scenario");
	    ap.currentIndexM++;
	}else
	    Matrix.perror("AGCT.class :: not in a manifold scenario");
    }

    public static void initDefaultClassVariables(){
	Number_Of_Neighbors = 10;
	Number_Of_Wavelet_Stamps = -1;
	T_Heat_Kernel = 2.0;
	T_Heat_Kernel_Similarity_Neighbors = 0.1;
	Sparsify_Statistic = 0.0;
	Number_Of_Manifold_Components = 20;
	Number_Of_Triangulation_Dimensions = 3;
	Max_Number_Of_Features = -1;
    }

    public static int Max_Number_Of_Genes_Kept = 10;
    // maximum number of genes to prevent storage capacity problems in Strings


    ///////////////////////////////////////////////////////////////////////////////////////////

    static double DISTANCE_PARAMETER = 10.0;
    // in units, distance max between genes or entities when filtered with distance which => similarity = 0
    // or T heat kernel if distances smooth similarities

    static double MAX_DISTANCE_TRIANGULATION = 10.0;
    // Puts edges only when distances < value

    /**************************************************************************************
     * Graphics part
     *****/
    //Buttons and fields that appear on the top buttonbar

    JButton loadButton, saveButton, saveButtonStats, saveButtonDelaunay, saveButtonWeb, binButton, cameraButton, scenarioLoadButton, scenarioSaveButton, scenarioDeleteButton, scenarioRecordButton, scenarioDisplayButton;
    ImageIcon scenarioRecordIdleIcon, scenarioRecordingIcon, scenarioLoadIdleIcon, scenarioRunningIcon;
    JLabel clockLabel;

    //items and alike that appear on the menus

    JCheckBoxMenuItem useDebug, useWarning, sortDepth, useShadow, perspective, onlyReferenced, onlyReferencedEdges, noGeneSelection, saveConstantWavelet, filterSimilaritiesWithDistance, filterTriangulationWithDistance, correctionProfiles;
    JRadioButtonMenuItem 
	methodF_0, methodF_1, methodF_2, methodF_3, methodF_4, methodF_5, methodF_6, methodF_7, methodF_8, methodF_9, methodF_10,
	methodS_0, methodS_1, methodS_2,
	methodW_0, methodW_1, methodW_2,
	methodN_0, methodN_1, methodN_2, methodN_3, methodN_4, methodN_5,
	methodE_0, methodE_1,
	methodD_0, methodD_1;
    JMenuItem numberOfNeighbors, betheHessianFactor, kDimensionEstimation, displayDimension, numberOfManifoldComponents, numberOfTriangulationDimensions, tHeatKernel, tHeatKernelSimilarityNeighbors, sparsifyStatistic, limitp, bregDiv, limitpdelaunay, numberOfProfiles, loadHighlightFile, setWC, filterDistanceFile, chooseDistanceParameterSimilarity, chooseDistanceParameterTriangulation;

    //"pointer" on the master frame
    JFrame myFrame;

    //"pointer" on the JinformationFrame uses to display all informations about processing
    JInformationFrame myInformationFrame;

    //id but for the Annotation Frame
    JAnnotationFrame myAnnotationFrame;

    //id but to display global informations about clusterings
    JClusteringFrame myClusteringFrame;

    //id but to display cluster profiles
    JClusteringProfileFrame myClusteringProfileFrame;

    //tabbedPane
    JMainFrameTabbedPane myTabbedPane;

    //DSScaling
    JDSScaling myScaling;

    boolean timerIsTicking = false, triangulationTimerIsTicking = false;
    AGCTTimer myTickingTimer;
    AGCTTimer myTriangulationTickingTimer;
    AGCTTimer myScenarioTimer;

    /**************************************************************************************
     * filter with distances
     *****/
    File myFilterDistanceFile;

    /**************************************************************************************
     * variables
     *****/

    //variables to get ressource data
    AGCTRessourceFileManager myRFM;

    //int used for Scenario
    int currentIndexW, currentIndexM;

    //FileWriter for convenience
    AGCTFileWriter myAGCTFileWriter;

    //pointer on the Domain of genes
    Domain myDomain;
    File myDomainFile;


    // Used when filtering with distances
    Vector allGenesCoordinates = new Vector();
    int [] geneToCoord;

    int dimElements;
    // dimension of all square matrices below:
    Matrix W, CW, D_from_W, N_from_W, DWD, M;
    // W = Similarity matrix
    // CW = Copy of W for clustering;
    // D_from_W = Diagonal matrix built on W
    // N_from_W = Normalized matrix = f(D, W)
    // DWD = toy normalized matrix
    // M = Manifold eigenvector matrix

    int dimFeatures;
    // dimension of all square matrices below:
    Matrix C, E, V;
    // C = Correlation matrix
    // E = Eigenvectors of C
    // V = Correlation circle variable components (in ROWS)

    Vector feature_pca_Pnt3D;
    Pnt3D min_CorrelationPnt3D;
    Pnt3D max_CorrelationPnt3D;

    Vector feature_names;
    //All names (typically, ligand name + # feature if Haar for example)

    double [] average_features;
    double [] sigma_features;
    // average and sigmas for features

    double [] manifold_eigenvalues;
    double [] pca_eigenvalues;

    int [][] nearest_neighbors;
    // nearest neighbors for each gene

    int [][] edges;
    // edges in the triangulation

    
    /**************************************************************************************
     * clustering
     *****/
    Vector allClusterings;
    int nbClustering;

    // some nice booleans

    boolean rawDomainExists; // yes iff raw data already loaded
    boolean annotationsExist, annotationsExistP, annotationsExistF, annotationsExistC; //yes iff annotations loaded
    boolean isTriangulated; // yes if triangulated
    boolean pca_computed; //yes iff pca computed
    boolean soft_clustering_computed; //yes iff a soft membership is computed
    boolean hard_clustering_computed; //yes iff a hard membership is computed

    AGCT(JBugFreeFrame myF){
	super();
	AGCT.Loading_From_Scenario_Manifold_Begin = AGCT.Loading_From_Scenario_Manifold_End = false;
	AGCT.Loading_From_Scenario_Dimension_Begin = AGCT.Loading_From_Scenario_Dimension_End = false;

	Referenced_Available = false;
	myFrame = myF;
	AGCT.RG = new Random();
	AGCTClustering_CP.init();

	myAGCTFileWriter = new AGCTFileWriter(this);
	myAnnotationFrame = null;
	myClusteringFrame = new JClusteringFrame("AGCT --- Clustering Result Visualization Frame", this);
	myClusteringProfileFrame = new JClusteringProfileFrame("AGCT --- [Clustering|Gene] Profile Visualization Frame", this);
	myScaling = new JDSScaling(this);

	min_CorrelationPnt3D = new Pnt3D();
	max_CorrelationPnt3D = new Pnt3D();

	nbClustering = 0;
	allClusterings = null;
	currentIndexW = currentIndexM  = -1;

	allGenesCoordinates = null;
	geneToCoord = null;

	Statistics.flushAll();
    }

    public Gene getDirectGene(int iGene){
	return (Gene) myDomain.domainGenes.elementAt(iGene);
    }

    public static void setCluster(int v){
	Number_Of_Clusters = v;
    }

    public static void setSort_Depth(boolean b){
	Sort_Depth = b;
    }

    public static void setNo_GeneSelection(boolean b){
	No_GeneSelection = b;
    }

    public static void setUse_Shadow(boolean b){
	Use_Shadow = b;
    }

    public static void setPerspective(boolean b){
	Perspective = b;
    }

    public static void setOnly_Referenced(boolean b){
	Only_Referenced = b;
    }

   public static void setOnly_ReferencedEdges(boolean b){
	Only_ReferencedEdges = b;
    }

    public static void setDebug(boolean b){
	Debug = b;
    }

    public static void setWarning(boolean b){
	Warning = b;
    }

    public static void setSaveConstantWavelet(boolean b){
       
	if ( (HAAR_WAVELETS_REMOVE_CONSTANT) && (b) )
	    MAX_WAVELET_COEFFICIENTS_SAVE++;
	else if ( (!HAAR_WAVELETS_REMOVE_CONSTANT) && (!b) )
	    MAX_WAVELET_COEFFICIENTS_SAVE--;
	
	HAAR_WAVELETS_REMOVE_CONSTANT = !b;
    }

    public static void setFilterSimilaritiesWithDistances(boolean b){
	Filter_Similarities_With_Distances = b;
    }

    public static void setFilterTriangulationWithDistances(boolean b){
	Filter_Triangulation_With_Distances = b;
    }

    public static void setCorrectionProfiles(boolean b){
	Correction_Profiles_By_Average = b;
	if (b)
	    System.out.println("Profile correction on :: WARNING --- Profiles corrected wrt average observed variation");
    }

    public static void setK_Dimension_Estimation(int n){
	K_Dimension_Estimation = n;
    }

    public static void setDisplay_Dimension(int n){
	JAGCTGraphicsPanel.DISPLAY_K_DIMENSION_ESTIMATION = n;
    }

    public static void setNumber_Of_Neighbors(int n){
	Number_Of_Neighbors = n;
    }

    public static void setBethe_Hessian_Factor(double n){
	Bethe_Hessian_Factor = n;
    }

    public static void setNumber_Of_Manifold_Components(int n){
	Number_Of_Manifold_Components = n;
    }

    public static void setNumber_Of_Wavelet_Stamps(int n){
	Number_Of_Wavelet_Stamps = n;
    }

    public static void setWC_Save(int n){
	MAX_WAVELET_COEFFICIENTS_SAVE = n;
    }

    public static void setDistanceParameterSimilarity(double n){
	DISTANCE_PARAMETER = n;
    }

    public static void setDistanceParameterTriangulation(double n){
	MAX_DISTANCE_TRIANGULATION = n;
    }

    public static void setMax_Number_Of_Features(int n){
	Max_Number_Of_Features = n;
    }

    public static void setMax_Number_Of_Prototypes(int n){
	Max_Number_Of_Prototypes = n;
    }

    public static void setNumber_Of_Profiles(int n){
	DEFAULT_Max_Number_Genes = n;
    }

    public static void setNumber_Of_Triangulation_Dimensions(int n){
	Number_Of_Triangulation_Dimensions = n;
    } 

    public void setBregDiv(int n){
	myScaling.indexSelected = n;
    }

    public static void setT_Heat_Kernel(double n){
	T_Heat_Kernel = n;
    }

    public static void setT_Heat_Kernel_Similarity_Neighbors(double n){
	T_Heat_Kernel_Similarity_Neighbors = n;
    }

    public static void setSparsify_Statistic(double n){
	Sparsify_Statistic = n;
    }

    public static void setLimit_P_Chi2(double p){
	Statistics.LIMIT_P_CHI2 = p;
    }

    public static void setLimit_P_Delaunay(double p){
	Statistics.LIMIT_P_DELAUNAY = p;
    }

    public static void setMethod_F(int v){
	Method_F = v;
    }

    public static void setMethod_S(int v){
	Method_S = v;
    }

    public static void setMethod_W(int v){
	Method_W = v;
    }

    public static void setMethod_N(int v){
	Method_N = v;
    }

    public static void setMethod_E(int v){
	Method_E = v;
    }

    public static void setMethod_D(int v){
	Method_D = v;
    }

    public void setMyFrame(JFrame f){
	myFrame = f;
    }


    public void requestModificationMethod_F(int v){
	if (v==0){
	    methodF_0.setSelected(true);
	}else if (v==1){
	    methodF_1.setSelected(true);
	}else if (v==2){
	    methodF_2.setSelected(true);
	}else if (v==3){
	    methodF_3.setSelected(true);
	}else if (v==4){
	    methodF_4.setSelected(true);
	}else if (v==5){
	    methodF_5.setSelected(true);
	}else if (v==6){
	    methodF_6.setSelected(true);
	}else if (v==7){
	    methodF_7.setSelected(true);
	}else if (v==8){
	    methodF_8.setSelected(true);
	}else if (v==9){
	    methodF_9.setSelected(true);
	}else if (v==10){
	    methodF_10.setSelected(true);
	}
	AGCT.setMethod_F(v);
    }

    public void requestMofificationMethod_F(){
	if (! ( (methodF_0.isSelected()) 
		|| (methodF_1.isSelected()) 
		|| (methodF_2.isSelected()) 
		|| (methodF_3.isSelected()) 
		|| (methodF_4.isSelected()) 
		|| (methodF_5.isSelected())
		|| (methodF_6.isSelected()) 
		|| (methodF_7.isSelected()) 
		|| (methodF_8.isSelected()) 
		|| (methodF_9.isSelected()) 
		|| (methodF_10.isSelected()) ) )
	    Matrix.perror("No radio button for Features selected");
	if (methodF_0.isSelected())
	    AGCT.setMethod_F(0);
	if (methodF_1.isSelected())
	    AGCT.setMethod_F(1);
	if (methodF_2.isSelected())
	    AGCT.setMethod_F(2);
	if (methodF_3.isSelected())
	    AGCT.setMethod_F(3);
	if (methodF_4.isSelected())
	    AGCT.setMethod_F(4);
	if (methodF_5.isSelected())
	    AGCT.setMethod_F(5);
 	if (methodF_6.isSelected())
	    AGCT.setMethod_F(6);
	if (methodF_7.isSelected())
	    AGCT.setMethod_F(7);
	if (methodF_8.isSelected())
	    AGCT.setMethod_F(8);
	if (methodF_9.isSelected())
	    AGCT.setMethod_F(9);
	if (methodF_10.isSelected())
	    AGCT.setMethod_F(10);

	Scenario.add("AGCT_Modification_Method_F", AGCT.Method_F);
    }

    public void requestModificationMethod_S(int v){
	if (v==0)
	    methodS_0.setSelected(true);
	else if (v==1)
	    methodS_1.setSelected(true);
	else if (v==2)
	    methodS_2.setSelected(true);
	AGCT.setMethod_S(v);
    }

    public void requestMofificationMethod_S(){
	if (! ( (methodS_0.isSelected()) || (methodS_1.isSelected()) || (methodS_2.isSelected()) ) )
	    Matrix.perror("No radio button for Similarity selected");
	if (methodS_0.isSelected())
	    AGCT.setMethod_S(0);
	if (methodS_1.isSelected())
	    AGCT.setMethod_S(1);
	if (methodS_2.isSelected())
	    AGCT.setMethod_S(2);

	Scenario.add("AGCT_Modification_Method_S", AGCT.Method_S);
    }

    public void requestModificationMethod_W(int v){
	if (v==0)
	    methodW_0.setSelected(true);
	else if (v==1)
	    methodW_1.setSelected(true);
	AGCT.setMethod_W(v);
    }

    public void requestMofificationMethod_W(){
	if (! ( (methodW_0.isSelected()) || (methodW_1.isSelected()) ) )
	    Matrix.perror("No radio button for Matrix W selected");
	if (methodW_0.isSelected())
	    AGCT.setMethod_W(0);
	if (methodW_1.isSelected())
	    AGCT.setMethod_W(1);

	Scenario.add("AGCT_Modification_Method_W", AGCT.Method_W);
    }

    public void requestModificationMethod_N(int v){
	if (v==0)
	    methodN_0.setSelected(true);
	else if (v==1)
	    methodN_1.setSelected(true);
	else if (v==2)
	    methodN_2.setSelected(true);
	else if (v==3)
	    methodN_3.setSelected(true);
	else if (v==4)
	    methodN_4.setSelected(true);
	else if (v==5)
	    methodN_5.setSelected(true);
	AGCT.setMethod_N(v);
    }

   public void requestModificationMethod_E(int v){
	if (v==0)
	    methodE_0.setSelected(true);
	else if (v==1)
	    methodE_1.setSelected(true);
	AGCT.setMethod_E(v);
    }

   public void requestModificationMethod_D(int v){
	if (v==0)
	    methodD_0.setSelected(true);
	else if (v==1)
	    methodD_1.setSelected(true);
	AGCT.setMethod_D(v);
    }

    public void requestMofificationMethod_N(){
	if (! ( (methodN_0.isSelected()) || (methodN_1.isSelected()) || (methodN_2.isSelected()) || (methodN_3.isSelected()) || (methodN_4.isSelected()) || (methodN_5.isSelected()) ) )
	    Matrix.perror("No radio button for Matrix N selected");
	if (methodN_0.isSelected())
	    AGCT.setMethod_N(0);
	if (methodN_1.isSelected())
	    AGCT.setMethod_N(1);
	if (methodN_2.isSelected())
	    AGCT.setMethod_N(2);
	if (methodN_3.isSelected())
	    AGCT.setMethod_N(3);
	if (methodN_4.isSelected())
	    AGCT.setMethod_N(4);
	if (methodN_5.isSelected())
	    AGCT.setMethod_N(5);

	Scenario.add("AGCT_Modification_Method_N", AGCT.Method_N);
    }

    public void requestMofificationMethod_E(){
	if (! ( (methodE_0.isSelected()) || (methodE_1.isSelected()) ) )
	    Matrix.perror("No radio button for Eigenmethod selected");
	if (methodE_0.isSelected())
	    AGCT.setMethod_E(0);
	if (methodE_1.isSelected())
	    AGCT.setMethod_E(1);

	Scenario.add("AGCT_Modification_Method_E", AGCT.Method_E);
    }


    public void requestMofificationMethod_D(){
	if (! ( (methodD_0.isSelected()) || (methodD_1.isSelected()) ) )
	    Matrix.perror("No radio button for Distance filtering selected");
	if (methodD_0.isSelected())
	    AGCT.setMethod_D(0);
	if (methodD_1.isSelected())
	    AGCT.setMethod_D(1);

	Scenario.add("AGCT_Modification_Method_D", AGCT.Method_D);
    }

    public void requestModificationDisplay_Dimension(int v){
	AGCT.setDisplay_Dimension(v);
    }

    public void requestModificationDisplay_Dimension(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Display of local dimension\n0 = no display\n1 = pies\n2 = numbers",(new Integer(JAGCTGraphicsPanel.DISPLAY_K_DIMENSION_ESTIMATION)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = 1;
	    }
	    if ( (i == 0) || (i == 1) || (i == 2) ){
		AGCT.setDisplay_Dimension(i);
		Scenario.add("AGCT_Display_Dimension",i);
		repaint();
	    }
	}
    }


    public void requestModificationK_Dimension_Estimation(int v){
	AGCT.setK_Dimension_Estimation(v);
    }

    public void requestModificationK_Dimension_Estimation(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Number of Neighbors\nused in Dimension estimation,\n in [" + AGCT.Min_K_Dimension_Estimation + ", " + AGCT.Max_K_Dimension_Estimation + "[",(new Integer(AGCT.K_Dimension_Estimation)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = 10;
	    }
	    if ( (i>0) && (i>=AGCT.Min_K_Dimension_Estimation) && (i<AGCT.Max_K_Dimension_Estimation) ){
		AGCT.setK_Dimension_Estimation(i);
		Scenario.add("AGCT_Modification_K_Dimension_Estimation",i);
		repaint();
	    }
	}
    }

    public void requestModificationBethe_Hessian_Factor(double v){
	AGCT.setBethe_Hessian_Factor(v);
    }

    public void requestModificationBethe_Hessian_Factor(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Bethe Hessian factor\n(q in r = +/- sqrt(q k), default is 2.0)",(new Double(AGCT.Bethe_Hessian_Factor)).toString());
	double i;

	if (ret != null){
	    try{
		i = Double.parseDouble(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an double !");
		i = 2.0;
	    }
	    if (i>0){
		AGCT.setBethe_Hessian_Factor(i);
		Scenario.add("AGCT_Modification_Bethe_Hessian_Factor",i);
	    }
	}
    }

    public void requestModificationNumber_Of_Neighbors(int v){
	AGCT.setNumber_Of_Neighbors(v);
    }

    public void requestModificationNumber_Of_Neighbors(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Number of Neighbors\n(Useful when SNN for W chosen)",(new Integer(AGCT.Number_Of_Neighbors)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = 5;
	    }
	    if (i>0){
		AGCT.setNumber_Of_Neighbors(i);
		Scenario.add("AGCT_Modification_Number_Of_Neighbors",i);
	    }
	}
    }

    public void requestModificationNumber_Of_Profiles(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Max # of profiles displayed\n(in clustering profile frame)",(new Integer(AGCT.DEFAULT_Max_Number_Genes)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = 10;
	    }
	    if (i>0){
		AGCT.setNumber_Of_Profiles(i);
	    }
	}
    }

   public void requestModificationMax_Number_Of_Features(int v){
       AGCT.setMax_Number_Of_Features(v);
   }

   public void requestModificationMax_Number_Of_Prototypes(int v){
       AGCT.setMax_Number_Of_Prototypes(v);
   }


    public void requestModificationMax_Number_Of_Features(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Maximum Number of Final Features by Gene:\n* -1 implies no selection, \n* or > 2",(new Integer(AGCT.Max_Number_Of_Features)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = -1;
	    }
	    if ( (i>2) || (i==-1) ){
		AGCT.setMax_Number_Of_Features(i);
		Scenario.add("JAGCTSelectionPane_Max_Number_Of_Features",i);
	    }
	}
    }


    public void requestModificationMax_Number_Of_Prototypes(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Maximum Number of Prototypes:\n* -1 implies no selection, \n* or > 2",(new Integer(AGCT.Max_Number_Of_Prototypes)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = -1;
	    }
	    if ( (i>2) || (i==-1) ){
		AGCT.setMax_Number_Of_Prototypes(i);
		Scenario.add("JAGCTSelectionPane_Max_Number_Of_Prototypes",i);
	    }
	}
    }

    public void defaultNWS(){
	requestModificationNumber_Of_Wavelet_Stamps(AGCT.DEFAULT_Number_Of_Wavelet_Stamps);
    }

   public void requestModificationNumber_Of_Wavelet_Stamps(int v){
       AGCT.setNumber_Of_Wavelet_Stamps(v);
   }

    public void requestModificationNumber_Of_Wavelet_Stamps(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Number of Wavelet Stamps:\n* -1 to let the program choose, \n* or a power of 2 > 1",(new Integer(AGCT.Number_Of_Wavelet_Stamps)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = -1;
	    }
	    if ( ( (i>1) && (Feature.powerOfTwo(i)) ) || ( (i==-1) && (!Feature.checkAllWaveletStampsAutomatic(this, myDomain.domainTimes)) ) ){
		AGCT.setNumber_Of_Wavelet_Stamps(i);
		Scenario.add("AGCT_Modification_Number_Of_Wavelet_Stamps",i);
	    }else if ( (!Feature.powerOfTwo(i)) || (i==-1) ){
		defaultNWS();
		Scenario.add("AGCT_Modification_Number_Of_Wavelet_Stamps", AGCT.Number_Of_Wavelet_Stamps);
	    }
	}
    }

   public void requestModificationNumber_Of_Manifold_Components(int v){
       AGCT.setNumber_Of_Manifold_Components(v);
   }

    public void requestModificationNumber_Of_Manifold_Components(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Number of Manifold Components to keep",(new Integer(AGCT.Number_Of_Manifold_Components)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = 5;
	    }
	    if (i>2){
		AGCT.setNumber_Of_Manifold_Components(i);
		Scenario.add("AGCT_Modification_Number_Of_Manifold_Components",i);
	    }
	}
    }

    public void requestModificationBregDiv(int i){
	myScaling.indexSelected = i;
    }

    public void requestModificationBregDiv(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Bregman divergence for DS Scaling\n" + JDSScaling.getString(),(new Integer(myScaling.indexSelected)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = 0;
	    }
	    if (JDSScaling.isValid(i)){
		myScaling.indexSelected = i;
		Scenario.add("JDSScaling_Modification_Bregman_Divergence",i);
	    }
	}
    }

    public void requestModificationNumber_Of_Triangulation_Dimensions(int v){
	AGCT.setNumber_Of_Triangulation_Dimensions(v);
    }

    public void requestModificationNumber_Of_Triangulation_Dimensions(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Which dimension for the triangulation ?\n(must be <= #Manifold components)",(new Integer(AGCT.Number_Of_Triangulation_Dimensions)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = 3;
	    }
	    if ( (i>2) && (i<=AGCT.Number_Of_Manifold_Components) ){
		AGCT.setNumber_Of_Triangulation_Dimensions(i);
		Scenario.add("AGCT_Modification_Number_Of_Triangulation_Dimensions",i);
	    }
	}
    }


    public void requestChooseDistanceParameterSimilarity(){
	JOptionPane d = new JOptionPane();
	double i;
	String ret = d.showInputDialog(this,"Distance parameter\n(When used with distance filtering, \nput 0 to all similarities for which\n the distance between the entities is >= this value,\n or uses it in Heat kernel to convolute similarities)",(new Double(AGCT.DISTANCE_PARAMETER)).toString());
	
	if (ret != null){
	    try{
		i = Double.parseDouble(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not a double !");
		i = 2.0;
	    }
	    if (i>=0.0){
		AGCT.setDistanceParameterSimilarity(i);
	    }
	}
    }
	
    public void requestChooseDistanceParameterTriangulation(){
	JOptionPane d = new JOptionPane();
	double i;
	String ret = d.showInputDialog(this,"Distance parameter\n(When used in triangulation, \nremoves edges whose object dist > value)",(new Double(AGCT.MAX_DISTANCE_TRIANGULATION)).toString());
	
	if (ret != null){
	    try{
		i = Double.parseDouble(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not a double !");
		i = 2.0;
	    }
	    if (i>=0.0){
		AGCT.setDistanceParameterTriangulation(i);
	    }
	}
    }
	

   public void requestModificationWC_Save(){
	JOptionPane d = new JOptionPane();
	String ret = d.showInputDialog(this,"Max number of Wavelet Coefficients used when Saving ?\n(not too large)",(new Integer(AGCT.MAX_WAVELET_COEFFICIENTS_SAVE)).toString());
	int i;

	if (ret != null){
	    try{
		i = Integer.parseInt(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not an integer !");
		i = 3;
	    }
	    if (i>0){
		AGCT.setWC_Save(i);
	    }
	}
    }

 

    public void requestModificationT_Heat_Kernel(double v){
	AGCT.setT_Heat_Kernel(v);
    }

    public void requestModificationT_Heat_Kernel(){
	JOptionPane d = new JOptionPane();
	double i;
	String ret = d.showInputDialog(this,"T Heat Kernel\n(W)",(new Double(AGCT.T_Heat_Kernel)).toString());
	
	if (ret != null){
	    try{
		i = Double.parseDouble(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not a double !");
		i = 2.0;
	    }
	    if (i>0){
		AGCT.setT_Heat_Kernel(i);
		Scenario.add("AGCT_Modification_T_Heat_Kernel",i);
	    }
	}
    }

    public void requestModificationT_Heat_Kernel_Similarity_Neighbors(double v){
	AGCT.setT_Heat_Kernel_Similarity_Neighbors(v);
    }

    public void requestModificationT_Heat_Kernel_Similarity_Neighbors(){
	JOptionPane d = new JOptionPane();
	double i;
	String ret = d.showInputDialog(this,"T Heat Kernel\n(2D manifold display)",(new Double(AGCT.T_Heat_Kernel_Similarity_Neighbors)).toString());
	
	if (ret != null){
	    try{
		i = Double.parseDouble(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not a double !");
		i = 2.0;
	    }
	    if (i>0){
		AGCT.setT_Heat_Kernel_Similarity_Neighbors(i);
		Scenario.add("AGCT_Modification_T_Heat_Kernel_Similarity_Neighbors",i);
	    }
	}
    }

    public void requestModificationSparsify_Statistic(double v){
	AGCT.setSparsify_Statistic(v);
    }

    public void requestModificationSparsify_Statistic(){
	JOptionPane d = new JOptionPane();
	double i;
	String ret = d.showInputDialog(this,"Sparsify Statistic\n(sparsify 100p% entries of W)",(new Double(AGCT.Sparsify_Statistic)).toString());
	
	if (ret != null){
	    try{
		i = Double.parseDouble(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not a double !");
		i = 2.0;
	    }
	    if ( (i>=0.0) && (i<=1.0) ){
		AGCT.setSparsify_Statistic(i);
		Scenario.add("AGCT_Modification_Sparsify_Statistic",i);
	    }
	}
    }

    public void requestModificationLimit_P_Chi2(double v){
	AGCT.setLimit_P_Chi2(v);
	if (myClusteringFrame != null)
	  myClusteringFrame.graphPanel.update=true;
    }

    public void requestModificationLimit_P_Delaunay(double v){
	AGCT.setLimit_P_Delaunay(v);
	if (myClusteringFrame != null)
	  myClusteringFrame.graphPanel.update=true;
    }

    public void requestModificationLimit_P_Chi2(){
	JOptionPane d = new JOptionPane();
	double i;
	String ret = d.showInputDialog(this,"Limit P-value for Chi2\n(typically < 0.05)",(new Double(Statistics.LIMIT_P_CHI2)).toString());
	
	if (ret != null){
	    try{
		i = Double.parseDouble(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not a double !");
		i = 0.01;
	    }
	    if (i<0.0)
		i = 0.0;
	    if (i>1.0)
		i = 1.0;

	    if (i>0.0){
		AGCT.setLimit_P_Chi2(i);
		Scenario.add("Statistics_Modification_Limit_P",i);
	    }
	    
	    if (myClusteringFrame != null)
		myClusteringFrame.graphPanel.update=true;
	}
    }

    public void requestModificationLimit_P_Delaunay(){
	JOptionPane d = new JOptionPane();
	double i;
	String ret = d.showInputDialog(this,"Limit P-value for plotting Delaunay triangulation\n(not too small)",(new Double(Statistics.LIMIT_P_DELAUNAY)).toString());
	
	if (ret != null){
	    try{
		i = Double.parseDouble(ret);
	    }catch (NumberFormatException n){
		System.out.println("Not a double !");
		i = 0.3;
	    }
	    if (i<0.0)
		i = 0.0;
	    if (i>1.0)
		i = 1.0;

	    if (i>0.0){
		AGCT.setLimit_P_Delaunay(i);
		//Scenario.add("Statistics_Modification_Limit_P_Delaunay",i); TBA
	    }
	    
	    if (myClusteringFrame != null)
		myClusteringFrame.graphPanel.update=true;
	}
    }

    public void requestModificationSort_Depth(boolean v){
	sortDepth.setSelected(v);
	requestModificationSort_Depth();
    }

    public void requestModificationSort_Depth(){
	AGCT.setSort_Depth(sortDepth.isSelected());
	Scenario.add("AGCT_Modification_Sort_Depth",sortDepth.isSelected());
	repaint();
    }

    public void requestModificationNo_GeneSelection(boolean v){
	noGeneSelection.setSelected(v);
	requestModificationNo_GeneSelection();
    }

    public void requestModificationNo_GeneSelection(){
	AGCT.setNo_GeneSelection(noGeneSelection.isSelected());
    }

    public void requestModificationUse_Shadow(boolean v){
	useShadow.setSelected(v);
	requestModificationUse_Shadow();
    }

    public void requestModificationUse_Shadow(){
	AGCT.setUse_Shadow(useShadow.isSelected());
	Scenario.add("AGCT_Modification_Use_Shadow",useShadow.isSelected());
	repaint();
    }

    public void requestModificationPerspective(boolean v){
	perspective.setSelected(v);
	requestModificationPerspective();
    }

    public void requestModificationOnly_Referenced(boolean v){
	onlyReferenced.setSelected(v);
	requestModificationOnly_Referenced();
    }

    public void requestModificationOnly_ReferencedEdges(boolean v){
	onlyReferencedEdges.setSelected(v);
	requestModificationOnly_ReferencedEdges();
    }

    public void requestModificationPerspective(){
	AGCT.setPerspective(perspective.isSelected());	
	Scenario.add("AGCT_Modification_Perspective",perspective.isSelected());
	myTabbedPane.myManifoldPane.visualizationPanel.myView3D.setPerspective(AGCT.Perspective);
	myTabbedPane.myPCAPane.visualizationPanel.myView3D.setPerspective(AGCT.Perspective);
	myTabbedPane.myCorrelationPane.visualizationPanel.myView3D.setPerspective(AGCT.Perspective);
	repaint();
    }

    public void requestModificationOnly_Referenced(){
	AGCT.setOnly_Referenced(onlyReferenced.isSelected());	
	myTabbedPane.myManifoldPane.visualizationPanel.myView3D.setOnlyReferenced(AGCT.Only_Referenced);
	myTabbedPane.myPCAPane.visualizationPanel.myView3D.setOnlyReferenced(AGCT.Only_Referenced);
	myTabbedPane.myCorrelationPane.visualizationPanel.myView3D.setOnlyReferenced(AGCT.Only_Referenced);
	repaint();
    }

    public void requestModificationOnly_ReferencedEdges(){
	AGCT.setOnly_ReferencedEdges(onlyReferencedEdges.isSelected());	
	myTabbedPane.myManifoldPane.visualizationPanel.myView3D.setOnlyReferencedEdges(AGCT.Only_ReferencedEdges);
	myTabbedPane.myPCAPane.visualizationPanel.myView3D.setOnlyReferencedEdges(AGCT.Only_ReferencedEdges);
	myTabbedPane.myCorrelationPane.visualizationPanel.myView3D.setOnlyReferencedEdges(AGCT.Only_ReferencedEdges);
	repaint();
    }

    public void requestModificationUseDebug(boolean v){
	useDebug.setSelected(v);
	requestModificationUseDebug();
    }

    public void requestModificationUseDebug(){
	AGCT.setDebug(useDebug.isSelected());
	Scenario.add("AGCT_Modification_Use_Debug",useDebug.isSelected());
    }

    public void requestModificationUseWarning(boolean v){
	useWarning.setSelected(v);
	requestModificationUseWarning();
    }

    public void requestModificationUseWarning(){
	AGCT.setWarning(useWarning.isSelected());
	Scenario.add("AGCT_Modification_Warning",useWarning.isSelected());
     }

    public void requestModificationSaveConstantWavelet(boolean v){
	saveConstantWavelet.setSelected(v);
	requestModificationSaveConstantWavelet();
    }

    public void requestModificationSaveConstantWavelet(){
	AGCT.setSaveConstantWavelet(saveConstantWavelet.isSelected());
	Scenario.add("AGCT_Save_Constant_Wavelet",saveConstantWavelet.isSelected());
     }

    public void requestModificationFilterSimilaritiesWithDistances(boolean v){
	filterSimilaritiesWithDistance.setSelected(v);
	requestModificationFilterSimilaritiesWithDistances();
    }

    public void requestModificationFilterTriangulationWithDistances(boolean v){
	filterTriangulationWithDistance.setSelected(v);
	requestModificationFilterTriangulationWithDistances();
    }

    public void requestModificationCorrectionProfiles(boolean v){
	correctionProfiles.setSelected(v);
	requestModificationCorrectionProfiles();
    }

    public void requestModificationCorrectionProfiles(){
	AGCT.setCorrectionProfiles(correctionProfiles.isSelected());
	Scenario.add("AGCT_Correction_Profiles",correctionProfiles.isSelected());
     }


    public void requestModificationFilterSimilaritiesWithDistances(){
	AGCT.setFilterSimilaritiesWithDistances(filterSimilaritiesWithDistance.isSelected());
	Scenario.add("AGCT_Filter_Similarities_With_Distance",filterSimilaritiesWithDistance.isSelected());
     }

    public void requestModificationFilterTriangulationWithDistances(){
	AGCT.setFilterTriangulationWithDistances(filterTriangulationWithDistance.isSelected());
	Scenario.add("AGCT_Filter_Triangulation_With_Distance",filterTriangulationWithDistance.isSelected());
     }

    public void requestLoad_Raw_Domain(String fname){
	myInformationFrame.appendText("Opening file: " + fname);
	rawDomainExists = true;
	myDomainFile = new File(fname);

	myDomain = new Domain(myDomainFile, this);
	myTabbedPane.setDomain(myDomain);
	ControlProcess.put("dataLoaded",true);
    }


    public void requestChooseFilterDistanceFile(){
	JFileChooser chooser = new JFileChooser(myRFM.Path_Data_Files);
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("agct");
	filter.setDescription(".agct Files");
	chooser.setFileFilter(filter);
	int returnVal = chooser.showOpenDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    myInformationFrame.appendText("Selected distance file: " +
			       chooser.getSelectedFile().getName());
	    myFilterDistanceFile = chooser.getSelectedFile();
	}
	
    }

    public void requestLoad_Raw_Domain(){
	JFileChooser chooser = new JFileChooser(myRFM.Path_Data_Files);
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("agct");
	filter.setDescription(".agct Files");
	chooser.setFileFilter(filter);
	int returnVal = chooser.showOpenDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    myInformationFrame.appendText("Opening file: " +
			       chooser.getSelectedFile().getName());
	    rawDomainExists = true;
	    myDomainFile = chooser.getSelectedFile();

	    myDomain = new Domain(myDomainFile, this);
	    myTabbedPane.setDomain(myDomain);
	    Scenario.add("Load",chooser.getSelectedFile().getAbsolutePath());
	    ControlProcess.put("dataLoaded",true);
	}
    }

    public void requestLoad_HighlightFile(){
	JFileChooser chooser = new JFileChooser(myRFM.Path_Highlight_Files);
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("txt");
	filter.setDescription(".txt Files");
	chooser.setFileFilter(filter);
	int returnVal = chooser.showOpenDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    myInformationFrame.appendText("Opening file: " +
			       chooser.getSelectedFile().getName());
	    myDomain.loadHighlights(chooser.getSelectedFile());
	}
    }

    public void requestLoadExecute_Scenario(){
	JFileChooser chooser = new JFileChooser(myRFM.Path_Scenario_Files);
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("txt");
	filter.setDescription(".txt Files");
	chooser.setFileFilter(filter);
	chooser.setApproveButtonText("Load Scenario");
	int returnVal = chooser.showOpenDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    myInformationFrame.setText("Opening and executing scenario: " +
				       chooser.getSelectedFile().getName());
	    Scenario.loadExecute(chooser.getSelectedFile(), this);
	    ControlProcess.put("scenarioLoaded",true);
	    myInformationFrame.setText("Opening and executing scenario: " +
				       chooser.getSelectedFile().getName() + " done...");
	}
    }

    public void requestLoad_Clustering(){
	JFileChooser chooser = new JFileChooser(myRFM.Path_Clustering_Files);
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("txt");
	filter.setDescription(".txt Files");
	chooser.setFileFilter(filter);
	chooser.setApproveButtonText("Load/Exec. Clustering");
	int returnVal = chooser.showOpenDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    myInformationFrame.setText("Opening and executing clustering: " +
				       chooser.getSelectedFile().getName());
	    myAGCTFileWriter.loadClustering(chooser.getSelectedFile(), this);
	    myInformationFrame.setText("Opening and executing clustering: " +
				       chooser.getSelectedFile().getName() + " done...");
	}
    }

    public void requestSave_Scenario(){
	if ( (Scenario.allStrings != null) && (Scenario.allStrings.size() > 0) ){
	    JFileChooser chooser = new JFileChooser(myRFM.Path_Scenario_Files);
	    ExampleFileFilter filter = new ExampleFileFilter();
	    filter.addExtension("txt");
	    filter.setDescription(".txt Files");
	    chooser.setFileFilter(filter);
	    chooser.setApproveButtonText("Save Scenario");
	    int returnVal = chooser.showSaveDialog(myFrame);
	    if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
		myInformationFrame.setText("Saving scenario to file: " +
					   chooser.getSelectedFile().getName() + " (please wait)...");
		myAGCTFileWriter.toSavingScenario(chooser.getSelectedFile());
		ControlProcess.put("scenarioSaved",true);
		myInformationFrame.setText("Saving scenario to file: " +
					   chooser.getSelectedFile().getName() + " done...");		
	    }
	}
    }

    public void requestSave_WC(){
	JFileChooser chooser = new JFileChooser(myRFM.Path_Data_Files);
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("csv");
	filter.setDescription(".csv Files");
	chooser.setFileFilter(filter);
	chooser.setApproveButtonText("Save Data's Wavelet Coeff. (Max = " + AGCT.MAX_WAVELET_COEFFICIENTS_SAVE +")");
	int returnVal = chooser.showSaveDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    myInformationFrame.setText("Saving data's wavelet coefficients to file: " +
			       chooser.getSelectedFile().getName());
	    myAGCTFileWriter.toSavingWC(chooser.getSelectedFile());
	    myInformationFrame.setText("Saving data's wavelet coefficients to file: " +
			       chooser.getSelectedFile().getName() + " done...");
	}
    }

    public void requestSave_Data(){
	JFileChooser chooser = new JFileChooser(myRFM.Path_Data_Files);
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("txt");
	filter.setDescription(".txt Files");
	chooser.setFileFilter(filter);
	chooser.setApproveButtonText("Save Data");
	int returnVal = chooser.showSaveDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    myInformationFrame.setText("Saving processed data to file: " +
			       chooser.getSelectedFile().getName());
	    myAGCTFileWriter.toSaving(chooser.getSelectedFile());
	    ControlProcess.put("dataSaved",true);
	    myInformationFrame.setText("Saving processed data to file: " +
			       chooser.getSelectedFile().getName() + " done...");
	}
    }

    public void requestSave_ThisClustering(int nclustering){
	JFileChooser chooser = new JFileChooser(myRFM.Path_Clustering_Files);
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("txt");
	filter.setDescription(".txt Files");
	chooser.setFileFilter(filter);
	chooser.setApproveButtonText("Save The Visible Clustering");
	int returnVal = chooser.showSaveDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    myInformationFrame.setText("Saving clustering to file: " +
			       chooser.getSelectedFile().getName());
	    myAGCTFileWriter.toSavingClustering(chooser.getSelectedFile(), nclustering);
	    ControlProcess.put("dataSaved",true);
	    myInformationFrame.setText("Saving clustering to file: " +
			       chooser.getSelectedFile().getName() + " done...");
	}
    }

    public void requestSave_DataChi2(){
	if (Statistics.allChi2Tests){
	    JFileChooser chooser = new JFileChooser(myRFM.Path_Chi2_Files);
	    ExampleFileFilter filter = new ExampleFileFilter();
	    filter.addExtension("txt");
	    filter.setDescription(".txt Files");
	    chooser.setFileFilter(filter);
	    chooser.setApproveButtonText("Save Data");
	    int returnVal = chooser.showSaveDialog(myFrame);
	    if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
		myInformationFrame.setText("Saving Chi2 data to file: " +
					   chooser.getSelectedFile().getName());
		myAGCTFileWriter.toSavingChi2(chooser.getSelectedFile(), this);
		ControlProcess.put("dataSavedChi2",true);
		myInformationFrame.setText("Saving Chi2 data to file: " +
					   chooser.getSelectedFile().getName() + " done...");
	    }
	}
    }

    public void requestSave_ClusterProfile(){
	JFileChooser chooser = new JFileChooser(myRFM.Path_Clustering_Files);
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("txt");
	filter.setDescription(".txt Files");
	chooser.setFileFilter(filter);
	chooser.setApproveButtonText("Save Data");
	int returnVal = chooser.showSaveDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    myInformationFrame.setText("Saving Cluster profile data to file: " +
				       chooser.getSelectedFile().getName());
	    myAGCTFileWriter.toSavingClusterProfile(chooser.getSelectedFile());
	    ControlProcess.put("dataSavedClusterProfile",true);
	    myInformationFrame.setText("Saving Cluster profile data to file: " +
				       chooser.getSelectedFile().getName() + " done...");
	}
    }


    public void requestSave_DataDelaunay(){
	if (isTriangulated){
	    JFileChooser chooser = new JFileChooser();
	    ExampleFileFilter filter = new ExampleFileFilter();
	    filter.addExtension("txt");
	    filter.setDescription(".txt Files");
	    chooser.setFileFilter(filter);
	    chooser.setApproveButtonText("Save Data");
	    int returnVal = chooser.showSaveDialog(myFrame);
	    if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
		myInformationFrame.setText("Saving Triangulated data to file: " +
					   chooser.getSelectedFile().getName());
		myAGCTFileWriter.toSavingDelaunay(chooser.getSelectedFile());
		ControlProcess.put("dataSavedDelaunay",true);
		myInformationFrame.setText("Saving Triangulated data to file: " +
					   chooser.getSelectedFile().getName() + " done...");
	    }
	}
    }

   public void requestSave_DataWeb(){
       JFileChooser chooser = new JFileChooser();
       ExampleFileFilter filter = new ExampleFileFilter();
       filter.addExtension("txt");
       filter.setDescription(".txt Files");
       chooser.setFileFilter(filter);
       chooser.setApproveButtonText("Save (Web)");
       int returnVal = chooser.showSaveDialog(myFrame);
       if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	   myInformationFrame.setText("Saving manifold to web file: " +
				      chooser.getSelectedFile().getName());
	   myAGCTFileWriter.toSavingWeb(chooser.getSelectedFile());
	   myInformationFrame.setText("Saving manifold to web file: " +
				      chooser.getSelectedFile().getName() + " done...");
       }
   }


    public void captureAndSave(){
	Rectangle rect = myTabbedPane.getCaptureRectangle();
	BufferedImage screencapture = null, framecapture = null;
	try{
	    framecapture = new Robot().createScreenCapture(rect);
	}catch(AWTException e){}
	
	JFileChooser chooser = new JFileChooser();
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("png");
	filter.setDescription(".png Files");
	chooser.setFileFilter(filter);
	chooser.setApproveButtonText("Save");
	int returnVal = chooser.showSaveDialog(myFrame);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    try{
		ImageIO.write(framecapture, "png", chooser.getSelectedFile());
	    }catch(IOException e){}
	    myInformationFrame.setText("Saving processed data to file: " +
				       chooser.getSelectedFile().getName());
	    ControlProcess.put("frameCaptured",true);
	}
    }

    public void setInfoTitle(String s){
	myInformationFrame.setTitle(s);
    }

    public void setInfoText(String s){
	myInformationFrame.setText(s);
    }

    public void appendInfoText(String s){
	myInformationFrame.appendText(s);
    }

    public void displayHTMLGene(Gene gg){
	myInformationFrame.setText(gg.toLightHTMLString());
    }
    
    public void init () {
	dimElements = -1;

        try {SwingUtilities.invokeAndWait(this);}
        catch (Exception e) {System.err.println("Initialization failure");}
    }

    public void initMenuItems(){
	betheHessianFactor.setEnabled(true);

	kDimensionEstimation.setEnabled(false);
	displayDimension.setEnabled(false);

	useDebug.setSelected(true);
	requestModificationUseDebug();

	useWarning.setSelected(true);
	requestModificationUseWarning();

	// NATURE BIO : true -> false
	saveConstantWavelet.setSelected(false);
	requestModificationSaveConstantWavelet();

	filterSimilaritiesWithDistance.setSelected(false);
	requestModificationFilterSimilaritiesWithDistances();

	filterTriangulationWithDistance.setSelected(false);
	requestModificationFilterTriangulationWithDistances();

	correctionProfiles.setSelected(false);
	requestModificationCorrectionProfiles();

	sortDepth.setSelected(false);
	requestModificationSort_Depth();

	useShadow.setSelected(false);
	requestModificationUse_Shadow();

	perspective.setSelected(false);
	requestModificationPerspective();

	onlyReferenced.setEnabled(false);
	onlyReferenced.setSelected(false);
	requestModificationOnly_Referenced();

	onlyReferencedEdges.setEnabled(false);
	onlyReferencedEdges.setSelected(false);
	requestModificationOnly_ReferencedEdges();

	loadHighlightFile.setEnabled(false);

	methodF_0.setSelected(true);
	methodF_1.setSelected(false);
	methodF_2.setSelected(false);
	methodF_3.setSelected(false);
	methodF_4.setSelected(false);
	methodF_5.setSelected(false);
	methodF_6.setSelected(false);
	methodF_7.setSelected(false);
	methodF_8.setSelected(false);
	methodF_9.setSelected(false);
	methodF_10.setSelected(false);


	// NATURE BIO : 2 -> 1
	methodS_0.setSelected(false);
	methodS_1.setSelected(true);
	methodS_2.setSelected(false);

	methodW_0.setSelected(false);
	methodW_1.setSelected(true);
	methodW_2.setSelected(false);
	methodN_0.setSelected(false);
	methodN_1.setSelected(true);
	methodN_2.setSelected(false);
	methodN_3.setSelected(false);
	methodN_4.setSelected(false);
	methodN_5.setSelected(false);

	methodE_0.setSelected(false);
	methodE_1.setSelected(true);

	methodD_0.setSelected(false);
	methodD_1.setSelected(true);

	requestMofificationMethod_F();
	requestMofificationMethod_S();
	requestMofificationMethod_W();
	requestMofificationMethod_N();
	requestMofificationMethod_E();
	requestMofificationMethod_D();
    }

    public String getParameterString(){
	String val;
	
	val = " Feature type: ";
	if (Method_F == 0)
	    val += "slopes\n";
	else if (Method_F == 1)
	    val += "Haar wavelets\n";
	else if (Method_F == 2)
	    val += "Daubechies D4 wavelets on an interval\n";
	else if (Method_F == 3)
	    val += "Daubechies D4 wavelets\n";
	else if (Method_F == 4)
	    val += "Daubechies D8 wavelets\n";
	else if (Method_F == 5)
	    val += "Daubechies D10 wavelets\n";
	else if (Method_F == 6)
	    val += "Daubechies D12 wavelets\n";
	else if (Method_F == 7)
	    val += "Daubechies D14 wavelets\n";
	else if (Method_F == 8)
	    val += "Daubechies D16 wavelets\n";
	else if (Method_F == 9)
	    val += "Daubechies D18 wavelets\n";
	else if (Method_F == 10)
	    val += "Daubechies D20 wavelets\n";

	val += " Similarity measure: ";
	if (Method_S == 0)
	    val += "Heat kernel with K = " + AGCT.T_Heat_Kernel + "\n";
	else if (Method_S == 1)
	    val += "Cosine similarity\n";

	val += " Similarity matrix W: ";
	if (Method_W == 0)
	    val += "used similarity as is\n";
	else if (Method_W == 1)
	    val += "filtered with K=" + AGCT.Number_Of_Neighbors + " symmetric neighbors\n";

	val += " Sparsification of W: " + (100*AGCT.Sparsify_Statistic) + "% of the smallest values become 0.0\n";

	val += " Normalization of W: ";
	if (Method_N == 0)
	    val += "Normalized Graph Laplacian, N = D^{-1/2} W D^{-1/2}\n";
	else if (Method_N == 1)
	    val += "Row-Stochastic Normalization, N = D^{-1} W\n";
	else if (Method_N == 2)
	    val += "Doubly-Stochastic approximation of W\n";
	else if (Method_N == 3)
	    val += "Doubly-Stochastic approximation of (Doubly-Stochastic Scaling of W)\n";
	else if (Method_N == 4)
	    val += "Bethe Hessian (Assortative) with r = sqrt(" + (2.0 * AGCT.Number_Of_Neighbors) + ") \n";
	else if (Method_N == 5)
	    val += "Bethe Hessian (Disassortative) with r = -sqrt(" + (2.0 * AGCT.Number_Of_Neighbors) + ") \n";
	
	if (Method_N == 3)
	    val += "Divergence chosen : " + JDSScaling.divergences[myScaling.indexSelected] + "\n";

	val += " Number of Manifold components kept for each gene (saves space): " + AGCT.Number_Of_Manifold_Components + "\n";
	val += " Number of dimensions used for triangulation: " + AGCT.Number_Of_Triangulation_Dimensions + "\n";

	val += "\n Plotting options at saving time: ";

	if (AGCT.Sort_Depth)
	    val += " Genes sorted wrt depth, ";
	else
	    val += " Depth not used, ";

	if (AGCT.Use_Shadow)
	    val += " Shadow plots, ";
	else
	    val += " Shadow not used, ";

	if (AGCT.Perspective)
	    val += " Perspective used.";
	else
	    val += " Perspective not used.";

	if (AGCT.Only_Referenced)
	    val += " Only referenced genes used.";
	else
	    val += " All genes used.";

	val += "\n";
	
	return val;
    }

    public JMenuBar menuAGCT(){
	useDebug = new JCheckBoxMenuItem("Debug");
	useDebug.setSelected(true);
	requestModificationUseDebug();
	useDebug.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationUseDebug();
		}
	    } );

	saveConstantWavelet = new JCheckBoxMenuItem("Keep Scaling Wavelet's Coeff.");
	saveConstantWavelet.setSelected(true);
	requestModificationSaveConstantWavelet();
	saveConstantWavelet.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationSaveConstantWavelet();
		}
	    } );

	filterSimilaritiesWithDistance = new JCheckBoxMenuItem("Filter Similarities with Distances");
	filterSimilaritiesWithDistance.setSelected(true);
	requestModificationFilterSimilaritiesWithDistances();
	filterSimilaritiesWithDistance.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationFilterSimilaritiesWithDistances();
		}
	    } );

	filterTriangulationWithDistance = new JCheckBoxMenuItem("Filter Triangulation with Distances");
	filterTriangulationWithDistance.setSelected(true);
	requestModificationFilterTriangulationWithDistances();
	filterTriangulationWithDistance.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationFilterTriangulationWithDistances();
		}
	    } );

	correctionProfiles = new JCheckBoxMenuItem("Correct Profiles with Average Vars");
	correctionProfiles.setSelected(false);
	requestModificationCorrectionProfiles();
	correctionProfiles.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationCorrectionProfiles();
		}
	    } );

	useWarning = new JCheckBoxMenuItem("Warning");
	useWarning.setSelected(true);
	requestModificationUseWarning();
	useWarning.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationUseWarning();
		}
	    } );

	sortDepth = new JCheckBoxMenuItem("Use depth (3D)");
	sortDepth.setSelected(false);
	requestModificationSort_Depth();
	sortDepth.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationSort_Depth();
		}
	    } );

	noGeneSelection = new JCheckBoxMenuItem("No gene selection");
	noGeneSelection.setSelected(AGCT.No_GeneSelection);
	requestModificationNo_GeneSelection();
	noGeneSelection.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationNo_GeneSelection();
		}
	    } );

	useShadow = new JCheckBoxMenuItem("Use shadow (3D)");
	useShadow.setSelected(false);
	requestModificationSort_Depth();
	useShadow.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationUse_Shadow();
		}
	    } );

	numberOfProfiles = new JMenuItem("Max #Profiles Displayed");
	numberOfProfiles.setMnemonic(KeyEvent.VK_D);
	numberOfProfiles.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationNumber_Of_Profiles();
		}
	    } );


	perspective = new JCheckBoxMenuItem("Perspective");
	perspective.setSelected(false);
	requestModificationPerspective();
	perspective.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationPerspective();
		}
	    } );


	displayDimension = new JMenuItem("Dimension display");
	displayDimension.setMnemonic(KeyEvent.VK_Y);
	displayDimension.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationDisplay_Dimension();
		}
	    } );


	onlyReferenced = new JCheckBoxMenuItem("Highlight: referenced genes");
	onlyReferenced.setSelected(false);
	requestModificationOnly_Referenced();
	onlyReferenced.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationOnly_Referenced();
		}
	    } );

	onlyReferencedEdges = new JCheckBoxMenuItem("Highlight: referenced edges");
	onlyReferencedEdges.setSelected(false);
	requestModificationOnly_ReferencedEdges();
	onlyReferencedEdges.addItemListener( new ItemListener() {
		public void itemStateChanged( ItemEvent e )
		{
		    requestModificationOnly_ReferencedEdges();
		}
	    } );
	
	loadHighlightFile = new JMenuItem("Load: highlight file");
	loadHighlightFile.setMnemonic(KeyEvent.VK_L);
	loadHighlightFile.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestLoad_HighlightFile();
		}
	    } );

	methodF_0 = new JRadioButtonMenuItem("Slopes");
	methodF_0.setActionCommand(AGCT.MF);

	methodF_1 = new JRadioButtonMenuItem("Haar wavelet coefficients");
	methodF_1.setActionCommand(AGCT.MF);

	methodF_2 = new JRadioButtonMenuItem("D4i wavelet coefficients");
	methodF_2.setActionCommand(AGCT.MF);

	methodF_3 = new JRadioButtonMenuItem("D4 wavelet coefficients");
	methodF_3.setActionCommand(AGCT.MF);

	methodF_4 = new JRadioButtonMenuItem("D8 wavelet coefficients");
	methodF_4.setActionCommand(AGCT.MF);

	methodF_5 = new JRadioButtonMenuItem("D10 wavelet coefficients");
	methodF_5.setActionCommand(AGCT.MF);

	methodF_6 = new JRadioButtonMenuItem("D12 wavelet coefficients");
	methodF_6.setActionCommand(AGCT.MF);

	methodF_7 = new JRadioButtonMenuItem("D14 wavelet coefficients");
	methodF_7.setActionCommand(AGCT.MF);

	methodF_8 = new JRadioButtonMenuItem("D16 wavelet coefficients");
	methodF_8.setActionCommand(AGCT.MF);

	methodF_9 = new JRadioButtonMenuItem("D18 wavelet coefficients");
	methodF_9.setActionCommand(AGCT.MF);

	methodF_10 = new JRadioButtonMenuItem("D20 wavelet coefficients");
	methodF_10.setActionCommand(AGCT.MF);

	ButtonGroup groupF = new ButtonGroup();
	groupF.add(methodF_0);
	groupF.add(methodF_1);
	groupF.add(methodF_2);
	groupF.add(methodF_3);
	groupF.add(methodF_4);
	groupF.add(methodF_5);
	groupF.add(methodF_6);
	groupF.add(methodF_7);
	groupF.add(methodF_8);
	groupF.add(methodF_9);
	groupF.add(methodF_10);

	methodS_0 = new JRadioButtonMenuItem("Heat kernel");
	methodS_0.setActionCommand(AGCT.MS);

	methodS_1 = new JRadioButtonMenuItem("Cosine similarity");
	methodS_1.setActionCommand(AGCT.MS);

	methodS_2 = new JRadioButtonMenuItem("Absolute Cosine similarity");
	methodS_2.setActionCommand(AGCT.MS);

	ButtonGroup groupS = new ButtonGroup();
	groupS.add(methodS_0);
	groupS.add(methodS_1);
	groupS.add(methodS_2);

	methodW_0 = new JRadioButtonMenuItem("Keep similarity as is");
	methodW_0.setActionCommand(AGCT.MW);

	methodW_1 = new JRadioButtonMenuItem("Filter via Symmetric NN");
	methodW_1.setActionCommand(AGCT.MW);

	methodW_2 = new JRadioButtonMenuItem("Filter via Natural Neighbors");
	methodW_2.setActionCommand(AGCT.MW);

	ButtonGroup groupW = new ButtonGroup();
	groupW.add(methodW_0);
	groupW.add(methodW_1);
	groupW.add(methodW_2);

	methodN_0 = new JRadioButtonMenuItem("Normalized Graph Laplacian(W)");
	methodN_0.setActionCommand(AGCT.MN);

	methodN_1 = new JRadioButtonMenuItem("Row Stochastic Normalization(W)");
	methodN_1.setActionCommand(AGCT.MN);

	methodN_2 = new JRadioButtonMenuItem("Doubly Stochastic Approximation(W)");
	methodN_2.setActionCommand(AGCT.MN);

	methodN_3 = new JRadioButtonMenuItem("DS Approximation(DS Scaling(W))");
	methodN_3.setActionCommand(AGCT.MN);

	methodN_4 = new JRadioButtonMenuItem("Bethe Hessian (Assortative) (W)");
	methodN_4.setActionCommand(AGCT.MN);

	methodN_5 = new JRadioButtonMenuItem("Bethe Hessian (Disassortative) (W)");
	methodN_5.setActionCommand(AGCT.MN);

	ButtonGroup groupN = new ButtonGroup();
	groupN.add(methodN_0);
	groupN.add(methodN_1);
	groupN.add(methodN_2);
	groupN.add(methodN_3);
	groupN.add(methodN_4);
	groupN.add(methodN_5);

	methodE_0 = new JRadioButtonMenuItem("QR");
	methodE_0.setActionCommand(AGCT.ME);

	methodE_1 = new JRadioButtonMenuItem("Lanczos");
	methodE_1.setActionCommand(AGCT.ME);

	ButtonGroup groupE = new ButtonGroup();
	groupE.add(methodE_0);
	groupE.add(methodE_1);


	methodD_0 = new JRadioButtonMenuItem("Hard threshold / distances");
	methodD_0.setActionCommand(AGCT.MD);

	methodD_1 = new JRadioButtonMenuItem("Heat Kernel convolution / distances");
	methodD_1.setActionCommand(AGCT.MD);

	ButtonGroup groupD = new ButtonGroup();
	groupD.add(methodD_0);
	groupD.add(methodD_1);


	JMenu tete = new JMenu("Comput. & Display Parameters");
	tete.setMnemonic(KeyEvent.VK_T);

	JMenu proc = new JMenu("Processing Parameters");
	proc.setMnemonic(KeyEvent.VK_P);

	JMenu about = new JMenu("About");
	about.setMnemonic(KeyEvent.VK_A);

	JMenu subF = new JMenu("Feature Types");
	subF.setMnemonic(KeyEvent.VK_F);
	subF.add(methodF_0);
	subF.add(methodF_1);
	subF.add(methodF_2);
	subF.add(methodF_3);
	subF.add(methodF_4);
	subF.add(methodF_5);
	subF.add(methodF_6);
	subF.add(methodF_7);
	subF.add(methodF_8);
	subF.add(methodF_9);
	subF.add(methodF_10);
	
	JMenu subS = new JMenu("Similarity Measure");
	subS.setMnemonic(KeyEvent.VK_S);
	subS.add(methodS_0);
	subS.add(methodS_1);
	subS.add(methodS_2);

	JMenu subW = new JMenu("Similarity Matrix W");
	subW.setMnemonic(KeyEvent.VK_S);
	subW.add(methodW_0);
	subW.add(methodW_1);
	
	JMenu subE = new JMenu("Computation of Eigensystem");
	subE.setMnemonic(KeyEvent.VK_C);
	subE.add(methodE_0);
	subE.add(methodE_1);
	
	JMenu subD = new JMenu("Distance Filtering method");
	subD.setMnemonic(KeyEvent.VK_F);
	subD.add(methodD_0);
	subD.add(methodD_1);
	
	JMenu subN = new JMenu("Normalized Matrix N");
	subN.setMnemonic(KeyEvent.VK_N);
	subN.add(methodN_0);
	subN.add(methodN_1);
	subN.add(methodN_2);
	subN.add(methodN_3);
	subN.add(methodN_4);
	subN.add(methodN_5);

	methodF_0.addActionListener(this);
	methodF_1.addActionListener(this);
	methodF_2.addActionListener(this);
	methodF_3.addActionListener(this);
	methodF_4.addActionListener(this);
	methodF_5.addActionListener(this);
	methodF_6.addActionListener(this);
	methodF_7.addActionListener(this);
	methodF_8.addActionListener(this);
	methodF_9.addActionListener(this);
	methodF_10.addActionListener(this);
	methodS_0.addActionListener(this);
	methodS_1.addActionListener(this);
	methodS_2.addActionListener(this);
	methodW_0.addActionListener(this);
	methodW_1.addActionListener(this);
	methodN_0.addActionListener(this);
	methodN_1.addActionListener(this);
	methodN_2.addActionListener(this);
	methodN_3.addActionListener(this);
	methodN_4.addActionListener(this);
	methodN_5.addActionListener(this);
	methodE_0.addActionListener(this);
	methodE_1.addActionListener(this);
	methodD_0.addActionListener(this);
	methodD_1.addActionListener(this);

	numberOfNeighbors = new JMenuItem("# Neighbors (SNN)");
	numberOfNeighbors.setMnemonic(KeyEvent.VK_O);
	numberOfNeighbors.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationNumber_Of_Neighbors();
		}
	    } );

	betheHessianFactor = new JMenuItem("Bethe Hessian factor (r)");
	betheHessianFactor.setMnemonic(KeyEvent.VK_B);
	betheHessianFactor.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationBethe_Hessian_Factor();
		}
	    } );

	kDimensionEstimation = new JMenuItem("# Neighbors (Dimension)");
	kDimensionEstimation.setMnemonic(KeyEvent.VK_D);
	kDimensionEstimation.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationK_Dimension_Estimation();
		}
	    } );

	numberOfManifoldComponents = new JMenuItem("Number Of Manifold Components");
	numberOfManifoldComponents.setMnemonic(KeyEvent.VK_M);
	numberOfManifoldComponents.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationNumber_Of_Manifold_Components();
		}
	    } );

	numberOfTriangulationDimensions = new JMenuItem("Dimension of the Triangulation");
	numberOfTriangulationDimensions.setMnemonic(KeyEvent.VK_T);
	numberOfTriangulationDimensions.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationNumber_Of_Triangulation_Dimensions();
		}
	    } );

	tHeatKernel = new JMenuItem("T Heat Kernel (W)");
	tHeatKernel.setMnemonic(KeyEvent.VK_T);
	tHeatKernel.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationT_Heat_Kernel();
		}
	    } );

	tHeatKernelSimilarityNeighbors = new JMenuItem("T Heat Kernel (Manifold)");
	tHeatKernelSimilarityNeighbors.setMnemonic(KeyEvent.VK_H);
	tHeatKernelSimilarityNeighbors.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationT_Heat_Kernel_Similarity_Neighbors();
		}
	    } );

	sparsifyStatistic = new JMenuItem("Sparsify value");
	sparsifyStatistic.setMnemonic(KeyEvent.VK_S);
	sparsifyStatistic.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationSparsify_Statistic();
		}
	    } );

	limitp = new JMenuItem("Limit P (Chi2)");
	limitp.setMnemonic(KeyEvent.VK_P);
	limitp.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationLimit_P_Chi2();
		}
	    } );

	limitpdelaunay = new JMenuItem("Limit P (Delaunay)");
	limitpdelaunay.setMnemonic(KeyEvent.VK_D);
	limitpdelaunay.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationLimit_P_Delaunay();
		}
	    } );

	setWC = new JMenuItem("Set max. #Wavelets_Coeff. for saving");
	setWC.setMnemonic(KeyEvent.VK_W);
	setWC.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationWC_Save();
		}
	    } );

	filterDistanceFile = new JMenuItem("Distance File");
	filterDistanceFile.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestChooseFilterDistanceFile();
		}
	    } );

	chooseDistanceParameterSimilarity = new JMenuItem("Distance Parameter (Similarities)");
	chooseDistanceParameterSimilarity.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestChooseDistanceParameterSimilarity();
		}
	    } );

	chooseDistanceParameterTriangulation = new JMenuItem("Distance Parameter (Triangulation)");
	chooseDistanceParameterTriangulation.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestChooseDistanceParameterTriangulation();
		}
	    } );


	bregDiv = new JMenuItem("Bregman Divergence (DS*)");
	bregDiv.setMnemonic(KeyEvent.VK_B);
	bregDiv.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    requestModificationBregDiv();
		}
	    } );

	JMenuItem aboutItem = new JMenuItem("About AGCT");
	aboutItem.setMnemonic(KeyEvent.VK_A);
	aboutItem.addActionListener( new ActionListener() {
		public void actionPerformed( ActionEvent e )
		{
		    showAbout();
		}
	    } );
	about.add(aboutItem);

	JMenuBar menuBar = new JMenuBar();
	menuBar.add(proc);
	menuBar.add(tete);
	menuBar.add(about);


	if (!AGCT.SUBMISSION){
	    tete.add(useDebug);
	    tete.add(useWarning);
	    tete.addSeparator();
	}

	if (!AGCT.SUBMISSION){
	    tete.add(noGeneSelection);	
	    tete.addSeparator();
	}

	tete.add(displayDimension);		
	tete.addSeparator();
	if (!AGCT.SUBMISSION)
	    tete.add(perspective);
	
	tete.add(sortDepth);
	tete.add(useShadow);
	tete.addSeparator();

	if (!AGCT.SUBMISSION){
	    tete.add(numberOfProfiles);
	    tete.addSeparator();
	}

	tete.add(loadHighlightFile);	
	tete.add(onlyReferenced);	
	tete.add(onlyReferencedEdges);	

	//proc.add(subF); Only accessible via the PANEL
	proc.add(subS);
	proc.add(subW);
	proc.add(subN);

	if (!AGCT.SUBMISSION)
	    proc.add(subE);
	proc.addSeparator();
	proc.add(betheHessianFactor);
	proc.addSeparator();
	proc.add(numberOfNeighbors);
	proc.add(kDimensionEstimation);
	proc.add(numberOfManifoldComponents);
	proc.add(numberOfTriangulationDimensions);
	proc.add(tHeatKernel);
	proc.add(tHeatKernelSimilarityNeighbors);
	proc.add(sparsifyStatistic);
	proc.add(bregDiv);
	proc.addSeparator();
	proc.add(limitp);
	proc.add(limitpdelaunay);
	proc.addSeparator();
	if (!AGCT.SUBMISSION)
	    proc.add(saveConstantWavelet);
	proc.add(setWC);

	if (AGCT.SHOW_FILTERING){
	    proc.addSeparator();
	    proc.add(correctionProfiles);
	    proc.addSeparator();
	    proc.add(filterSimilaritiesWithDistance);
	    proc.add(filterTriangulationWithDistance);
	    proc.add(filterDistanceFile);
	    proc.add(chooseDistanceParameterSimilarity);
	    proc.add(chooseDistanceParameterTriangulation);
	    proc.add(subD);
	}

	initMenuItems();

	return menuBar;
    }

    public void showAbout(){
	JOptionPane.showMessageDialog(this, History.getLatestVersion() + ", " + History.getLatestDate() + "\n Programming: R. Nock" + "\n Specific Routines: R. Nock, F. Nielsen, K. Oka\n Testing: N. Polouliakh", "About AGCT", JOptionPane.INFORMATION_MESSAGE, new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/information.png"))));
    }

    public void timerTickOn(){
	timerIsTicking = true;
    }

    public void timerTickOff(){
	timerIsTicking = false;
	clockLabel.setEnabled(false);
    }

    public void triangulationTimerTickOn(){
	triangulationTimerIsTicking = true;
    }

    public void triangulationTimerTickOff(){
	triangulationTimerIsTicking = false;
	myTabbedPane.myManifoldPane.delauButton.setIcon(myTabbedPane.myManifoldPane.triangulationIdleIcon);
    }

    public void setDomain(String v){

    }

    public void run () {
	myRFM = new AGCTRessourceFileManager();
	myRFM.loadRessource();

	Number_Of_Clusters = 3;
	isTriangulated = false;
	rawDomainExists = false;
	Scenario.tickOff();

        Container pane = getContentPane();
        pane.setLayout(new BorderLayout());

	loadButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/folder_go.png")))); 
	loadButton.setToolTipText("load domain (raw data)");
        loadButton.setActionCommand("load_raw_domain");
	loadButton.setEnabled(true);

	saveButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/disk.png")))); 
	saveButton.setToolTipText("save current results");
        saveButton.setActionCommand("save_data");
	saveButton.setEnabled(true);

	saveButtonStats = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/diskChi2.png")))); 
	saveButtonStats.setToolTipText("save Chi2 results");
        saveButtonStats.setActionCommand("save_data_chi2");
	saveButtonStats.setEnabled(false);

	saveButtonDelaunay = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/diskDel.png")))); 
	saveButtonDelaunay.setToolTipText("save the *filtered* Delaunay triangulation");
        saveButtonDelaunay.setActionCommand("save_data_delaunay");
	saveButtonDelaunay.setEnabled(false);

	saveButtonWeb = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/diskWeb.png")))); 
	saveButtonWeb.setToolTipText("save the *visible* pane coordinates");
        saveButtonWeb.setActionCommand("save_data_web");
	saveButtonWeb.setEnabled(false);

	binButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/bin_closed.png"))));
	binButton.setToolTipText("flush results");
        binButton.setActionCommand("flush");
	//loadButton.setEnabled(false);

	cameraButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/camera.png"))));
	cameraButton.setToolTipText("capture the visible graphics panel");
        cameraButton.setActionCommand("capture");
	//loadButton.setEnabled(false);

	clockLabel = new JLabel(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/time_add.png"))));
	clockLabel.setEnabled(false);

	myTickingTimer = new AGCTTimer(1000, this, AGCTTimer.GeneralTimerString);
	myTickingTimer.setInitialDelay(0);
	myTickingTimer.setCoalesce(true);
	myTickingTimer.start();

	myScenarioTimer = new AGCTTimer(1000, this, AGCTTimer.ScenarioTimerString);
	myScenarioTimer.setInitialDelay(0);
	myScenarioTimer.setCoalesce(true);
	myScenarioTimer.start();

	scenarioRecordIdleIcon = new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/script_add.png")));
	scenarioRecordIdleIcon.setDescription("Idle");
	scenarioRecordingIcon = new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/script_go.png")));
	scenarioRecordingIcon.setDescription("Recording");

	scenarioLoadIdleIcon = new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/folder_go.png")));
	scenarioLoadIdleIcon.setDescription("Idle");
	scenarioRunningIcon = new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/cog.png")));
	scenarioRunningIcon.setDescription("Running");

	scenarioLoadButton = new JButton(scenarioLoadIdleIcon); 
	scenarioLoadButton.setToolTipText("load and execute scenario");
        scenarioLoadButton.setActionCommand("load_scenario");

	scenarioSaveButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/disk.png")))); 
	scenarioSaveButton.setToolTipText("save currently recorded scenario");
        scenarioSaveButton.setActionCommand("save_scenario");

	scenarioDeleteButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/script_delete.png")))); 
	scenarioDeleteButton.setToolTipText("erase currently recorded scenario");
        scenarioDeleteButton.setActionCommand("delete_scenario");

	scenarioRecordButton = new JButton(scenarioRecordIdleIcon); 
	scenarioRecordButton.setToolTipText("record/stop recording scenario");
        scenarioRecordButton.setActionCommand("record_scenario");

	scenarioDisplayButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/magnifier.png")))); 
	scenarioDisplayButton.setToolTipText("display current scenario in information frame");
        scenarioDisplayButton.setActionCommand("display_scenario");

	myTriangulationTickingTimer = new AGCTTimer(1000, this, AGCTTimer.TriangulationTimerString);
	myTriangulationTickingTimer.setInitialDelay(0);
	myTriangulationTickingTimer.setCoalesce(true);
	myTriangulationTickingTimer.start();

	setDomain("");

	JToolBar mainBoxButtons = new JToolBar();
	mainBoxButtons.setFloatable(false);
	mainBoxButtons.setBorder(BorderFactory.createTitledBorder(""));

	JToolBar scenarioButtons = new JToolBar();
	scenarioButtons.setFloatable(false);
	scenarioButtons.setBorder(BorderFactory.createTitledBorder(""));

	JToolBar miscButtons = new JToolBar();
	miscButtons.setFloatable(false);
	miscButtons.setBorder(BorderFactory.createTitledBorder(""));

	JToolBar boxButtons = new JToolBar();

	mainBoxButtons.add(new JLabel("Data:"));
	mainBoxButtons.add(loadButton);
	mainBoxButtons.add(saveButton);
	mainBoxButtons.add(saveButtonStats);
	mainBoxButtons.add(saveButtonDelaunay);
	mainBoxButtons.addSeparator();
	mainBoxButtons.add(saveButtonWeb);

	scenarioButtons.add(new JLabel("Scenario:"));
	scenarioButtons.add(scenarioLoadButton);
	scenarioButtons.add(scenarioSaveButton);
	scenarioButtons.addSeparator();
	scenarioButtons.add(scenarioRecordButton);
	scenarioButtons.add(scenarioDisplayButton);
	scenarioButtons.addSeparator();
	scenarioButtons.add(scenarioDeleteButton);

	miscButtons.add(new JLabel("Misc:"));
	//miscButtons.add(cameraButton);
	//miscButtons.addSeparator();
	miscButtons.add(binButton);
	miscButtons.add(clockLabel);

	boxButtons.add(mainBoxButtons);
	boxButtons.add(Box.createHorizontalGlue());
	boxButtons.add(scenarioButtons);
	boxButtons.add(miscButtons);

	pane.add(boxButtons, BorderLayout.NORTH);

	myTabbedPane = new JMainFrameTabbedPane(this, myDomain);
	pane.add(myTabbedPane, BorderLayout.CENTER);

	//pane.add(lowPane, BorderLayout.SOUTH);

	setJMenuBar(menuAGCT());

	Image img = Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/chart_line.png"));
	myFrame.setIconImage( img );

        loadButton.addActionListener(this);
        saveButton.addActionListener(this);
        saveButtonStats.addActionListener(this);
        saveButtonDelaunay.addActionListener(this);
        saveButtonWeb.addActionListener(this);
        binButton.addActionListener(this);
        cameraButton.addActionListener(this);

        scenarioLoadButton.addActionListener(this);
        scenarioSaveButton.addActionListener(this);
        scenarioDeleteButton.addActionListener(this);
        scenarioRecordButton.addActionListener(this);
	scenarioDisplayButton.addActionListener(this);

	myInformationFrame = new JInformationFrame("AGCT --- Information Frame", this);
	myInformationFrame.setSize(AGCT.WindowWidth,400);

	if (AGCT.SUBMISSION)
	    myInformationFrame.setText(AGCT.SUBMISSION_TEXT);
    }

    

    public static void go(){
	Prototype.flushAll();
	JBugFreeFrame f = null;
	AGCT applet = null;
	
	f=new JBugFreeFrame("AGCT --- A Geometric Clustering Tool");
	applet = new AGCT(f);
	applet.init();
	
	
	f.setSize(AGCT.WindowWidth,700);
	f.getContentPane().setLayout(new BorderLayout()); 
	f.getContentPane().add(applet, "Center");
	
	f.addWindowListener(new FermetureListener("Closing AGCT\n"));
	
	f.setVisible(true);
	    
	//f.setAlwaysOnTop(true);
	//f.pack();
	//f.show();
    }


    public static void affiche_warning(){
	System.out.println(SUBMISSION_TEXT);
    }

    public static void main(String [] argv){
	if (SUBMISSION)
	    affiche_warning();
	go();
    }


    public void actionPerformed (ActionEvent e) {
	if (e == null)
	    return;

	String command = e.getActionCommand();
	String eClass = e.getSource().getClass().getName();
	String idleString;

	if (eClass == "AGCTTimer"){
	    if ( ((AGCTTimer) e.getSource()).name.equals(AGCTTimer.GeneralTimerString)){
		if (timerIsTicking == true){
		    if (clockLabel.isEnabled())
			clockLabel.setEnabled(false);
		    else
			clockLabel.setEnabled(true);
		}
	    }
	    else if ( ((AGCTTimer) e.getSource()).name.equals(AGCTTimer.ScenarioTimerString)){
		if (Scenario.isTicking()){
		    idleString = ((ImageIcon) scenarioRecordButton.getIcon()).getDescription();
		    if (idleString == "Idle")
			scenarioRecordButton.setIcon(scenarioRecordingIcon);
		    else
			scenarioRecordButton.setIcon(scenarioRecordIdleIcon);
		}
		if (Scenario.isRunning()){
		    idleString = ((ImageIcon) scenarioLoadButton.getIcon()).getDescription();
		    if (idleString == "Idle")
			scenarioLoadButton.setIcon(scenarioRunningIcon);
		    else
			scenarioLoadButton.setIcon(scenarioLoadIdleIcon);
		}
	    }
	    else if ( ( (AGCTTimer) e.getSource() ).name.equals(AGCTTimer.TriangulationTimerString)){
		if (triangulationTimerIsTicking == true){
		    idleString = ((ImageIcon) myTabbedPane.myManifoldPane.delauButton.getIcon()).getDescription();
		    if (idleString == "Idle")
			myTabbedPane.myManifoldPane.delauButton.setIcon(myTabbedPane.myManifoldPane.triangulationRunningIcon);
		    else
			myTabbedPane.myManifoldPane.delauButton.setIcon(myTabbedPane.myManifoldPane.triangulationIdleIcon);
		}
	    }
	}
	else if (command != null){
	    if (command.equals("load_raw_domain")){
		goLoadDomain();
	    }else if (command.equals("save_data")){
		requestSave_Data();
	    }else if (command.equals("save_data_chi2")){
		requestSave_DataChi2();
	    }else if (command.equals("save_data_delaunay")){
		requestSave_DataDelaunay();
	    }else if (command.equals("save_data_web")){
		requestSave_DataWeb();
	    }else if (command.equals("flush")){
		flushAll();
	    }else if (command.equals("capture")){
		captureAndSave();
	    }else if (command.equals("record_scenario")){
		Scenario.tickSwitch();
		if (!Scenario.isTicking())
		    scenarioRecordButton.setIcon(scenarioRecordIdleIcon);
	    }else if (command.equals("save_scenario")){
		requestSave_Scenario();
	    }else if (command.equals("load_scenario")){
		requestLoadExecute_Scenario();
	    }else if (command.equals("display_scenario")){
		Scenario.displayInformationFrame(this);
	    }else if (command.equals("delete_scenario")){
		if (confirmEraseScenario()){
		    Scenario.flushAll();
		    scenarioRecordButton.setIcon(scenarioRecordIdleIcon);
		}
	    }else if (command.equals(AGCT.MF)){
		requestMofificationMethod_F();
	    }else if (command.equals(AGCT.MW)){
		requestMofificationMethod_W();
	    }else if (command.equals(AGCT.MS)){
		requestMofificationMethod_S();
	    }else if (command.equals(AGCT.MN)){
		requestMofificationMethod_N();
	    }else if (command.equals(AGCT.ME)){
		requestMofificationMethod_E();
	    }else if (command.equals(AGCT.MD)){
		requestMofificationMethod_D();
	    }
	}

	//requestFocus();
	//(myTabbedPane.getSelectedComponent()).requestFocus();
    }
    
    public boolean confirmEraseScenario(){
	int i = JOptionPane.showConfirmDialog(this, "Confirm the erasing of the scenario");

	if (i==JOptionPane.YES_OPTION)
	    return true;
	return false;
    }


    public void goLoadDomain(){
	Thread t = new Thread(){
		public void run(){
		    requestLoad_Raw_Domain();
		    if (rawDomainExists == true){
			myTabbedPane.filterData();
			dimElements = myDomain.numberSelectedGenes;
		    }
		}
	    };
	t.start();
    }

    public void goLoadDomain(String s){
	final String s2 = s;
	Thread t = new Thread(){
		public void run(){
		    requestLoad_Raw_Domain(s2);
		    if (rawDomainExists == true){
			myTabbedPane.filterData();
			dimElements = myDomain.numberSelectedGenes;
		    }
		}
	    };
	t.start();
    }

    public void flushAll(){
	AGCT.initDefaultClassVariables();

	ControlProcess.removeAll();
	if (myAnnotationFrame != null)
	    myAnnotationFrame.dispose();
	myAnnotationFrame = null;
	annotationsExist = annotationsExistF = annotationsExistP = annotationsExistC = false;

	myClusteringFrame.setVisible(false);
	myClusteringFrame.flushAll();
	myClusteringProfileFrame.setVisible(false);
	myClusteringProfileFrame.flushAll();
	myScaling.flushAll();

	setDomain("");
	myDomain = null;
	myDomainFile = null;
	myTabbedPane.setDomain(myDomain);
	Gene.flushAll();
	Feature.flushAll();
	Prototype.flushAll();

	myTabbedPane.mySelectionPane.toTrash();
	myTabbedPane.myClusteringPane.toTrash();
	myTabbedPane.myPCAPane.toTrash();
	myTabbedPane.myManifoldPane.toTrash();
	myTabbedPane.myCorrelationPane.toTrash();

	myTabbedPane.setSelectedIndex(0);

	W = CW = D_from_W = N_from_W = DWD = M = C = E = V = null;
	AGCTClustering_CP.init();
	AGCT.Loading_From_Scenario_Manifold_Begin = AGCT.Loading_From_Scenario_Manifold_End = false;
	AGCT.Show_Correlation = -1;
	AGCT.Show_Only_Genes_With_Visible_Edge = false;
	AGCT.Referenced_Available = false;
	AGCT.Max_Reference = -1;
	isTriangulated = false;
	saveButtonDelaunay.setEnabled(false);
	saveButtonWeb.setEnabled(false);
	if (myTabbedPane.myManifoldPane.delauConsistencyButton != null)
	    myTabbedPane.myManifoldPane.delauConsistencyButton.setEnabled(false);
	if (myTabbedPane.myManifoldPane.saveClusteringButton != null)
	    myTabbedPane.myManifoldPane.saveClusteringButton.setEnabled(false);
	if (myTabbedPane.myManifoldPane.loadClusteringButton != null)
	    myTabbedPane.myManifoldPane.loadClusteringButton.setEnabled(false);
	if (myTabbedPane.myPCAPane.saveClusteringButton != null)
	    myTabbedPane.myPCAPane.saveClusteringButton.setEnabled(false);
	if (myTabbedPane.myPCAPane.loadClusteringButton != null)
	    myTabbedPane.myPCAPane.loadClusteringButton.setEnabled(false);

	initMenuItems();
	allClusterings = null;
	nbClustering = 0;
	currentIndexW = currentIndexM = -1;

	myInformationFrame.setText("Flushed all data");
	myInformationFrame.setTitle(JInformationFrame.defaultProgressString);

	rawDomainExists = false;
	Statistics.flushAll(this);
	Scenario.flushAll();

	repaint();
    }


    ////////////////////////////////////////////////////////////////////////////////////

    public void computeFeatureNames(){
	feature_names = new Vector();
	int i, j, rx;
	String ln;

	rx = 0;
	Gene gg = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumberBefore[0]);

	for (i=0;i<myDomain.numberLigands;i++){
	    ln = (String) ((Vector) myDomain.domainLigands.elementAt(i)).elementAt(0);
	    if (myDomain.selectedLigand[i]){
		if (gg.finalCoordinates[rx].length == 1)
		    feature_names.addElement(ln);
		else{
		    for (j=0;j<gg.finalCoordinates[rx].length;j++)
			feature_names.addElement(ln + (j+1) );
		}
		rx++;
	    }
	}

	//printAllFeatures();
    }

    public void printAllFeatures(){
	int i, j, rx = 0, ry, index = 0;
	Gene gg = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[0]);
	for (i=0;i<myDomain.numberLigands;i++){
	    if (myDomain.selectedLigand[i]){
		for (j=0;j<gg.finalCoordinates[rx].length;j++){
		    System.out.print("[" + i + ", " + j + "] Feature " + feature_names.elementAt(index) + " : " + Feature.trashed[index] + " -- ");
		    index++;
		}
		System.out.println();
		rx++;
	    }
	}
    }

    public String getFeatureName(int i){
	return (String) feature_names.elementAt( ( (Integer) Feature.featureIndexToNonTrashedIndex.get(new Integer(i)) ).intValue());
    }

    public void justBeforeManifold(){
	myTabbedPane.myManifoldPane.initDepthOrderSelectedGenes();
	myInformationFrame.setText("...");
    }

    public void dimension(){
	ControlProcess.put("justBeforeDimension", true);

	Thread t = new Thread(){
		public void run(){
		    timerTickOn();
		    myDomain.toNeighborsDimension();
		    ControlProcess.put("dimensionProcessed",true);
		    timerTickOff();
		}
	    };
	t.start();
    }

    public void manifoldEigensystem_W(){
	generate_W();

	// BEWARE : if W was prefiltered using distances, it may be already very very sparse
	if (AGCT.Sparsify_Statistic > 0.0)
	    W.sparsifyMatrix();

	myInformationFrame.setText("Matrix W (10 x 10 upperleft block):\n" + W.affiche(10));
    }

    public void manifoldEigensystem(boolean includeW){
	if ( (AGCT.SAVE_MEMORY) && (!AGCT.SUBMISSION) )
	    System.out.println("\n*************************************************************************************\nWarning : trying to save memory --- This may prevent clustering or several algorithms\n*************************************************************************************\n");

	if (includeW){
	    justBeforeManifold();
	}
	manifold_eigenvalues = new double[dimElements];
	ControlProcess.put("justBeforeManifold", true);
	M = new Matrix("Manifold Eigenvectors", dimElements, dimElements);

	final boolean includeWf = includeW;

	Thread t = new Thread(){
		public void run(){
		    
		    if ( (AGCT.Filter_Similarities_With_Distances) || (AGCT.Filter_Triangulation_With_Distances) )
			loadDistances();

		    int x, y;
		    double [] e = new double[dimElements];
		    double norm;

		    timerTickOn();

		    if (includeWf)
			manifoldEigensystem_W();

		    generate_D();

		    myInformationFrame.setText("Matrix D (10 x 10 upperleft block):\n" + D_from_W.affiche(10));

		    generate_N(); 
		    
		    myInformationFrame.setText("Matrix N (10 x 10 upperleft block):\n" + N_from_W.affiche(10));

		    if ( (Method_N == 0) || (Method_N == 2) || (Method_N == 3) || (Method_N == 4) || (Method_N == 5) )
			M.copyOnThis(N_from_W);
		    else if (Method_N == 1)
			M.copyOnThis(DWD);
		    

		    if (AGCT.SAVE_MEMORY){
			if (!AGCT.SUBMISSION)
			    System.out.println(" * Saving memory --- erasing DWD & N_from_W");
			DWD = null;
			N_from_W = null;
		    }

		    if (AGCT.Method_E == 0){
			M.tred2(myDomain, manifold_eigenvalues,e);
			M.tqli(myDomain, manifold_eigenvalues,e);

			M.eigenSort(manifold_eigenvalues);
			M.transposeSquare();
		    }else if (AGCT.Method_E == 1){
			AGCT.MAX_COL_OF_EIGENVECTOR = M.dimX / 5;
			if (AGCT.MAX_COL_OF_EIGENVECTOR > AGCT.MAX_MAX_COL_OF_EIGENVECTOR)
			    AGCT.MAX_COL_OF_EIGENVECTOR = AGCT.MAX_MAX_COL_OF_EIGENVECTOR;

			//System.out.println("Dim M = " + M.dimX + ", " + M.dimY); 
			//System.out.println("Mat 0 = + " + M.affiche(100));

			if (AGCT.USE_NATALIA)
			    Matrix.Natalia_getEigenSystem(myDomain, M, manifold_eigenvalues);
			else
			    Matrix.getEigenSystem(myDomain, M, manifold_eigenvalues);
		    }else
			Matrix.perror("Method_E : invalid parameter " + AGCT.Method_E);
		    
		    if ( (Method_N == 4) || (Method_N == 5) ){
			int i, sum = 0;
			for (i=0;i<manifold_eigenvalues.length;i++)
			    if (manifold_eigenvalues[i] > 0.0)
				sum++;
			System.out.println(sum + " < 0 eigenvalues ");
			if (sum < 3){
			    System.out.println(" ** Warning : too few negative eigenvalues for Bethe Hessian with AGCT.Method_N = " + AGCT.Method_N);
			    System.out.println("List of Eigenvalues");
			    for (i=0;i<manifold_eigenvalues.length;i++)
				System.out.println(i + " --- " + (-manifold_eigenvalues[i]));
			}
		    }
		    

		    /*
		    //Test Lanczos
		    Matrix M2 = new Matrix("Manifold Eigenvectors", M.dimX, M.dimY);
		    M2.copyOnThis(M);

		    double [] manifold_eigenvalues_2 = new double[M.dimX];
		    Matrix.getEigenSystem(myDomain, M2, manifold_eigenvalues_2);

		    //M2.times(-1.0);
		    // FIN test

		    M.tred2(myDomain, manifold_eigenvalues,e);
		    M.tqli(myDomain, manifold_eigenvalues,e);

		    myInformationFrame.setText("Matrix M (10 x 10 upperleft block):\n" + M.affiche(10));

		    M.eigenSort(manifold_eigenvalues);
		    M.transposeSquare();
		    */


		    /*M2.times(-1,0);
		    for (x=0;x<manifold_eigenvalues.length;x++)
			System.out.println(manifold_eigenvalues[x] + " " + manifold_eigenvalues_2[x]);
		    
		    System.out.println(M2.afficheComplet(10));
		    System.out.println("\n\n" + M.afficheComplet(10));

		    System.out.println("\n\n Frob = " + Matrix.FrobAbs(M, M2, 4) + "; Linfty = " + Matrix.LinftyAbs(M, M2, 4) );*/

		    if (Method_N == 1){
			for (x=0;x<dimElements;x++){
			    norm=0.0;
			    for (y=0;y<dimElements;y++){
				M.coordinates[x][y] = M.coordinates[x][y] / Math.sqrt(D_from_W.coordinates[y][y]);
				norm += (M.coordinates[x][y]*M.coordinates[x][y]);
			    }
			    for (y=0;y<dimElements;y++){
				M.coordinates[x][y] = M.coordinates[x][y] / Math.sqrt(norm);
			    }
			}
		    }
		    


		    if (AGCT.SAVE_MEMORY){
			if (!AGCT.SUBMISSION)
			    System.out.println(" * Saving memory --- erasing D_from_W");
			D_from_W = null;
		    }

		    //M.checkOrthonormality(myDomain, -1);
		    //Matrix.checkEigensystem(myDomain, N,M,manifold_eigenvalues,-1);
		    
		    //Scenario : save M from here.
		    justAfterManifold();

		    if (myTabbedPane.myManifoldPane.keepManifold.isSelected())
			Scenario.add("JAGCTManifoldPane_Keep_Manifold");

		    timerTickOff();

		}
	    };
	t.start();
    }

    public void justAfterManifold(){

	myDomain.fillManifoldComponents(M);
	
	myTabbedPane.myManifoldPane.setGenesButtons(true);
	myTabbedPane.myManifoldPane.setHiButtons(true);
	
	//propagating events to Swing units
	
	myTabbedPane.myManifoldPane.visualizationPanel.plotAvailable = true;
	
	myTabbedPane.myManifoldPane.initManifold3D();
	myTabbedPane.myManifoldPane.visualizationPanel.initAll(true, true, true);
	myTabbedPane.repaint();
	myTabbedPane.myManifoldPane.xText.setEnabled(true);
	myTabbedPane.myManifoldPane.yText.setEnabled(true);
	myTabbedPane.myManifoldPane.zText.setEnabled(true);
	myTabbedPane.myManifoldPane.delauButton.setEnabled(true);
	myTabbedPane.myManifoldPane.searchText.setEnabled(true);	
	myTabbedPane.myManifoldPane.viewSelect.setEnabled(true);
	
	myTabbedPane.myManifoldPane.keepManifold.setEnabled(false);
	

	ControlProcess.put("manifoldProcessed",true);
	checkClusteringAvailable();
	if (myTabbedPane.getSelectedIndex() == 1)
	    saveButtonWeb.setEnabled(true);
   }

    public void checkClusteringAvailable(){
	if ( (ControlProcess.hasTrue("manifoldProcessed")) && (ControlProcess.hasTrue("pcaProcessed")) ){
	    //Clustering becomes available !
	    myTabbedPane.toClustering();
	    myTabbedPane.myClusteringPane.activate();
	    ControlProcess.put("clusteringAvailable",true);
	}
    }

    public void pcaEigensystem(){
	myInformationFrame.setText("...");
	Thread t = new Thread(){
		public void run(){
		    timerTickOn();

		    generate_C(true);

		    myInformationFrame.setText("Matrix C (10 x 10 upperleft block):\n" + C.affiche(10));

		    generate_E();

		    myInformationFrame.setText("Matrix E (10 x 10 upperleft block):\n" + E.affiche(10));

		    myDomain.fillPCAComponents(E,average_features,sigma_features);

		    generate_V();

		    myInformationFrame.setText("Matrix V (10 x 10 upperleft block):\n" + V.affiche(10));

		    myTabbedPane.myPCAPane.setGenesButtons(true);
		    myTabbedPane.myPCAPane.setHiButtons(true);

		    myTabbedPane.myPCAPane.visualizationPanel.plotAvailable = true;
		    myTabbedPane.myPCAPane.initPCA3D();
		    myTabbedPane.myPCAPane.visualizationPanel.initAll(true, true, true);
		    myTabbedPane.repaint();
		    myTabbedPane.myPCAPane.xText.setEnabled(true);
		    myTabbedPane.myPCAPane.yText.setEnabled(true);
		    myTabbedPane.myPCAPane.zText.setEnabled(true);
		    myTabbedPane.myPCAPane.searchText.setEnabled(true);	
		    myTabbedPane.myPCAPane.viewSelect.setEnabled(true);

		    fillCorrelationComponents();
		    myTabbedPane.myCorrelationPane.visualizationPanel.plotAvailable = true;
		    myTabbedPane.myCorrelationPane.initCorrelation3D();
		    myTabbedPane.myCorrelationPane.visualizationPanel.initAll(true, true, true);
		    myTabbedPane.repaint();
		    myTabbedPane.myCorrelationPane.xText.setEnabled(true);
		    myTabbedPane.myCorrelationPane.yText.setEnabled(true);
		    myTabbedPane.myCorrelationPane.zText.setEnabled(true);
		    myTabbedPane.myCorrelationPane.searchText.setEnabled(true);	
		    myTabbedPane.myCorrelationPane.viewSelect.setEnabled(true);

		    ControlProcess.put("pcaProcessed",true);
		    checkClusteringAvailable();
		    myTabbedPane.myManifoldPane.loadClusteringButton.setEnabled(true);
		    myTabbedPane.myPCAPane.loadClusteringButton.setEnabled(true);
		    if (myTabbedPane.getSelectedIndex() == 2)
			saveButtonWeb.setEnabled(true);

		    timerTickOff();

		}
	    };
	t.start();
    }

    public void generate_E(){
	E = new Matrix("E", dimFeatures, dimFeatures);
	E.toE(myDomain, C);
    }
    
    public void generate_C(boolean after){
	C = new Matrix("C", dimFeatures, dimFeatures);
	C.toC(myDomain,after);
    }

    public void generate_V(){
	V = new Matrix("V", dimFeatures, dimFeatures);
	V.toV(myDomain);
    }
    

    public void computeMinMax_CorrelationPnt3D(){
	//Makes the assumption that V is computed

	int i, j;
	Pnt3D pt;
	for (i=0;i<dimFeatures;i++){
	    pt = (Pnt3D) feature_pca_Pnt3D.elementAt(i);
	    for (j=0;j<3;j++){
		if ( (i==0) || (pt.coordinates[j] < min_CorrelationPnt3D.coordinates[j]) )
		    min_CorrelationPnt3D.coordinates[j] = pt.coordinates[j];
		if ( (i==0) || (pt.coordinates[j] > max_CorrelationPnt3D.coordinates[j]) )
		    max_CorrelationPnt3D.coordinates[j] = pt.coordinates[j];
	    }
	}
    }

    public void fillCorrelationComponents(){
	int i;
	feature_pca_Pnt3D = new Vector();
	for (i=0;i<dimFeatures;i++){
	    feature_pca_Pnt3D.addElement(new Pnt3D(V.coordinates[i][0],
						   V.coordinates[i][1],
						   V.coordinates[i][2]));
	}
    }

    public void updateCorrelationComponents(int xAxis, int yAxis, int zAxis){
	int i;
	Pnt3D p;
	if (feature_pca_Pnt3D == null)
	    Matrix.perror("No features Pnt 3D !");
	for (i=0;i<dimFeatures;i++){
	    p = (Pnt3D) feature_pca_Pnt3D.elementAt(i);
	    p.coordinates[0] = V.coordinates[i][xAxis];
	    p.coordinates[1] = V.coordinates[i][yAxis];
	    p.coordinates[2] = V.coordinates[i][zAxis];
	}
    }

    public void initCorrelation3D(int xAxis, int yAxis, int zAxis){
	int i;
	Pnt3D pt;
	for (i=0;i<dimFeatures;i++){
	    pt = (Pnt3D) feature_pca_Pnt3D.elementAt(i);
	    pt.coordinates[0] = V.coordinates[i][xAxis];
	    pt.coordinates[1] = V.coordinates[i][yAxis];
	    pt.coordinates[2] = V.coordinates[i][zAxis];
	}
    }

    public synchronized void loadDistances(){
	//
	
	allGenesCoordinates = new Vector();

	System.out.println("Opening " + myFilterDistanceFile);
	
	File rf = myFilterDistanceFile;
	FileReader e;
	BufferedReader br;
	StringTokenizer t;
	String s, currentStep = "", valret;
	Vector v, v2, vt;
	boolean collectData = false;
	Gene g;
	Vector ge = null, co;
	int dim = myDomain.numberSelectedGenes;
	int i, j;
	String gin, gjn;
	boolean trouve;
	double dist;
	Gene gi, gj;
	double nb = 0.0, nbr = 0.0;

	try{
	    e = new FileReader(rf);
	}catch(FileNotFoundException ex){ 
	    myInformationFrame.setText("The file you try to load does not exist");
	    rawDomainExists = false;
	    return;
	}
	br = new BufferedReader(e);
	try{
	    while ( (s=br.readLine()) != null){
		if ( (s.length()>=AGCTFileWriter.DATA_Data.length()) && (s.equals(AGCTFileWriter.DATA_Data)) )
		    collectData = true;
		else if ( (collectData) && (s.length()>1) && (!s.substring(0,AGCTFileWriter.DATA_Comments.length()).equals(AGCTFileWriter.DATA_Comments))){
		    t = new StringTokenizer(s.replace('\t',' '), " ");

		    if (t.countTokens()>0) {
			s = t.nextToken();
			ge = new Vector();
			co = new Vector();
			ge.addElement(s);
			s = t.nextToken();
			co.addElement(new Double(Double.parseDouble(s)));
			s = t.nextToken();
			co.addElement(new Double(Double.parseDouble(s)));
			ge.addElement(co);
			allGenesCoordinates.addElement(ge);
		    }
		}
	    }
	    e.close();
	}catch(IOException ex){
	    myInformationFrame.setText("IOError when loading domain");
	    rawDomainExists = false;
	    return;
	}
	e = null;

	AGCTCounter cc = new AGCTCounter(myInformationFrame, "Searching genes in File " + myFilterDistanceFile.getName(), dim);
	geneToCoord = new int [dim];
	for (i=0;i<dim;i++){
	    j=0;
	    trouve = false;
	    gi = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[i]);
	    gin = gi.name;
	    do{
		gjn = (String) ((Vector) allGenesCoordinates.elementAt(j)).elementAt(0);
		if (gin.equals(gjn))
		    trouve = true;
		else
		    j++;
	    }while ( (!trouve) && (!gin.equals(gjn)) );
	    if (!trouve)
		Matrix.perror("Gene " + gin + " not found in distances");
	    else
		geneToCoord[i] = j;
	    cc.increment();
	}
	cc.end();
    }

    public void generate_W(){
	W = new Matrix("W", dimElements, dimElements);


	if (AGCT.USE_NATALIA)
	    W.Natalia_toW(this);
	else
	    W.toW(this);

	//System.out.println("W = " + W.affiche(10));

	if (AGCT.SAVE_MEMORY_NO_ANN){
	    CW = new Matrix("W_Original", dimElements, dimElements);
	    CW.copyOnThis(W);
	    generate_Neighbors();

	    if (AGCT.SAVE_MEMORY){
		if (!AGCT.SUBMISSION)
		    System.out.println(" * Saving memory --- erasing CW");
		CW = null;
	    }
	}

	ControlProcess.put("neighborsComputed",true);
    }
    
    public void generate_D(){
	D_from_W = new Matrix("D", dimElements, dimElements);
	D_from_W.toD(myDomain, W);
    }
    
    public void generate_N(){
	N_from_W = new Matrix("N", dimElements, dimElements);
	N_from_W.toN(this, D_from_W, W);
    }

    public void generate_Neighbors(){
	if (AGCT.Debug) System.out.print("Computing All Nearest Neighbors... ");
	myInformationFrame.setTextProgressBar("ANN");
	int i, j, k, index, ij, ik, dum, percent, max, nmade;
	double [] sim = new double [dimElements - 1];
	double ddum;

	max = dimElements ;
	nmade = 0;

	nearest_neighbors = new int [dimElements][];
	for (i=0;i<dimElements;i++){
	    nearest_neighbors[i] = new int [dimElements-1];
	    index = 0;
	    for (j=0;j<dimElements;j++)
		if (j!=i){
		    nearest_neighbors[i][index] = j;
		    sim[index] = CW.coordinates[i][j];
		    index++;
		}
	    QuickSort.quicksort(sim,nearest_neighbors[i]);
	    for (j=0;j<(dimElements-2)/2;j++){
		dum = nearest_neighbors[i][j];
		nearest_neighbors[i][j] = nearest_neighbors[i][dimElements - 2 - j];
		nearest_neighbors[i][dimElements - 2 - j] = dum;
		
		ddum = sim[j];
		sim[j] = sim[dimElements - 2 - j];
		sim[dimElements - 2 - j] = ddum;
	    }

	    percent = (int) (100.0 * ((double) nmade / (double) (max)));
	    myInformationFrame.setValueProgressBar(percent);
	    nmade++;
	}
	if (AGCT.Debug) System.out.println("ok.");
	myInformationFrame.setValueProgressBar(0);
	myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);
    }

    public void triangulation(){

	Thread t = new Thread(){
		public void run(){
		    myInformationFrame.setTextProgressBar("Delaunay Triang.");
		    myInformationFrame.setText("Batch processing of Delaunay Triangulation... ");
		    int i, j, d = AGCT.Number_Of_Triangulation_Dimensions, nmade = 0, percent;
		    Gene g;

		    timerTickOn();
		    triangulationTimerTickOn();

		    double [] min = new double [d], max = new double [d], delta = new double [d], extreme = new double [d], comp;
		    Pnt [] big;
		    double dum;
		    Simplex ori;
		    DelaunayTriangulation dt;

		    for (i=0;i<myDomain.numberSelectedGenes;i++){
			g = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[i]);
			comp = g.manifold_components_triangulation.coordinates;

			for (j=0;j<d;j++){
			    if ( (i==0) || (comp[j] < min[j]) )
				min[j] = comp[j];

			    if ( (i==0) || (comp[j] > max[j]) )
				max[j] = comp[j];
			}
		    }

		    for (i=0;i<d;i++)
			delta[i] = max[i] - min[i];
		    
		    big = new Pnt[d+1];
		    comp = new double[d];
		    for (i=0;i<d;i++)
			comp[i] = min[i] - delta[i];
		    big[0] = new Pnt(comp, -1);

		    for (i=0;i<d;i++){
			extreme[i] = max[i] + ( (double) d * 2)*delta[i];

			comp = new double[d];
			for (j=0;j<d;j++){
			    if (j==i)
				comp[j] = max[i] + ( (double) d * 2)*delta[i];
			    else
				comp[j] = 0.0;
			}
			big[i+1] = new Pnt(comp, -(i+1));
		    }

		    for (i=0;i<myDomain.numberSelectedGenes;i++){
			g = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[i]);
			comp = g.manifold_components_triangulation.coordinates;
			dum = 0.0;
			for (j=0;j<d;j++)
			    dum += ( comp[j] / extreme[j] ) ;
			if (dum >= 1.0)
			    Matrix.perror("Gene " + i + " is outside big simplex");
		    }

		    ori = new Simplex(big);
		    dt = new DelaunayTriangulation(ori);

		    for (i=0;i<myDomain.numberSelectedGenes;i++){
			percent = (int) (100.0 * ((double) nmade / (double) (myDomain.numberSelectedGenes)));
			myInformationFrame.setValueProgressBar(percent);
			nmade++;

			g = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[i]);
			comp = g.manifold_components_triangulation.coordinates;
			dt.delaunayPlace(new Pnt(comp, i));
		    }
		    myDomain.myDelaunayTriangulation = dt;

		    myInformationFrame.setText("Batch processing of Delaunay Triangulation... ok.");
		    myInformationFrame.setValueProgressBar(0);
		    myInformationFrame.setTextProgressBar(JInformationFrame.defaultProgressString);

		    myTabbedPane.myManifoldPane.setManifoldButtons(true);
		    myTabbedPane.myPCAPane.setManifoldButtons(true);

		    myDomain.toEdges();

		    triangulationTimerTickOff();
		    timerTickOff();

		    isTriangulated = true;
		    saveButtonDelaunay.setEnabled(true);
		    myTabbedPane.myManifoldPane.delauConsistencyButton.setEnabled(true);

		    ControlProcess.put("manifoldTriangulated",true);
		}
	    };
	t.start();
    }

    public void triangleConsistency(){
	Thread t = new Thread(){
		public void run(){
		    int allP = 0, oneN = 0, twoN = 0, allN = 0, tot = 0, nid, corel;
		    int i, j, k, xi1, xi2, xj1, xj2, xk1, xk2, ci, cj, ck, cij, cia, cja, cs, mi, ma, mi2, ma2, tmp, nPos = 0, nNeg = 0;
		    Vector eli, elj, elk;
		    Gene gg, ggn;
		    double nangle, maxl = -1.0, curl, maxltot = -1.0, curltot, sizetot = 0.0, avg = 0.0, avgtot = 0.0, sizequant = 0.0, sloc = 0.0, totloc = 0.0, stot = 0.0, tottot = 0.0;
		    boolean okPlot;
		    AGCTCounter cc;
		    
		    Vector listOfEdges = new Vector();
		    Vector element, element2;
		    String saff = "";
		    
		    cc = new AGCTCounter(myInformationFrame, "[Lengths] Step 1 (largest visible)", myDomain.numberSelectedGenes);
		    for (i=0;i<myDomain.numberSelectedGenes;i++){
			gg = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[i]);
			if ( (gg.neighbors != null) && (gg.neighbors.size()>0) ){
			    for (j=0;j<gg.neighbors.size();j++){
				nid = ((Integer) gg.neighbors.elementAt(j)).intValue();
				nangle = ((Double) gg.neighbor_angles.elementAt(j)).doubleValue();
				corel = 0;
				okPlot = false;
				if (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0){
				    okPlot = true;
				    corel = 1;
				    nPos++;
				}else if (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0))){
				    okPlot = true;
				    corel = -1;
				    nNeg++;
				}
				if (okPlot){
				    ggn = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[((Integer) gg.neighbors.elementAt(j)).intValue()]);
				    curl = Distortion.distortion_l22(gg.manifold_components_triangulation, ggn.manifold_components_triangulation);
				    sloc += 1.0;
				    totloc += curl;

				    if (curl>maxl)
					maxl = curl;
				}
			    }
			}
			cc.increment();
		    }
		    cc.end();

		    if (sloc > 0.0){
			avg = totloc / sloc;

			cc = new AGCTCounter(myInformationFrame, "[Lengths] Step 2 (quantile)", myDomain.numberSelectedGenes);
			for (i=0;i<myDomain.numberSelectedGenes;i++){
			    gg = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[i]);
			    if ( (gg.neighbors != null) && (gg.neighbors.size()>0) ){
				for (j=0;j<gg.neighbors.size();j++){
				    sizetot+=1.0;
				    ggn = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[((Integer) gg.neighbors.elementAt(j)).intValue()]);
				    curltot = Distortion.distortion_l22(gg.manifold_components_triangulation, ggn.manifold_components_triangulation);
				    tottot += curltot;
				    if (curltot > maxltot)
					maxltot = curltot;
				    if (curltot < maxl)
					sizequant += 1.0;
				}
			    }
			    cc.increment();
			}
			cc.end();
			avgtot = tottot / sizetot;

			saff += "Total observable edges / # >0 Corr. / # <0 Corr. = " + (nPos + nNeg) + " / " + nPos + " ( = " + ( 100.0 * (double) nPos / (nPos + nNeg) ) + "%) / " + nNeg + " ( = " + ( 100.0 * (double) nNeg / (nPos + nNeg) ) + "%)\n";
			saff += "Average observable length / Average length = " + avg + " / " + avgtot + " = " + DF.format((avg / avgtot)) + "\n";
			saff += "Max obervable length / Max length = " + maxl + " / " + maxltot + " = " + DF.format((maxl / maxltot)) + "\n";
			
			cc = new AGCTCounter(myInformationFrame, "[Triangles] Step 1 (edges)", myDomain.numberSelectedGenes);
			for (i=0;i<myDomain.numberSelectedGenes;i++){
			    gg = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[i]);
			    if ( (gg.neighbors != null) && (gg.neighbors.size()>0) ){
				for (j=0;j<gg.neighbors.size();j++){
				    nid = ((Integer) gg.neighbors.elementAt(j)).intValue();
				    nangle = ((Double) gg.neighbor_angles.elementAt(j)).doubleValue();
				    corel = 0;
				    okPlot = false;
				    if (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0){
					okPlot = true;
					corel = 1;
				    }else if (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0))){
					okPlot = true;
					corel = -1;
				    }
				    
				    if (okPlot){
					mi = myDomain.selectedGeneNumberToGeneNumber[i];
					ma = myDomain.selectedGeneNumberToGeneNumber[((Integer) gg.neighbors.elementAt(j)).intValue()];
					if (mi > ma){
					    tmp = mi;
					    mi = ma;
					    ma = tmp;
					}
					
					element = new Vector();
					element.addElement(new Integer(mi));
					element.addElement(new Integer(ma));
					element.addElement(new Integer(corel));
					listOfEdges.addElement(element);
				    }
				}
			    }
			    cc.increment();
			}
			cc.end();
			
			QuickSort.quicksortVectorInteger(listOfEdges,0);
			i=j=0;
			do{
			    element = (Vector) listOfEdges.elementAt(i);
			    mi = ((Integer) element.elementAt(0)).intValue();
			    ma = ((Integer) element.elementAt(1)).intValue();
			    j=i+1;
			    do{
				element2 = (Vector) listOfEdges.elementAt(j);
				mi2 = ((Integer) element2.elementAt(0)).intValue();
				ma2 = ((Integer) element2.elementAt(1)).intValue();
				if ( (mi == mi2) && (ma == ma2) )
				    listOfEdges.removeElementAt(j);
				else
				    j++;
			    }while(j<listOfEdges.size());
			    i++;
			}while(i<listOfEdges.size()-1);
			
			cc = new AGCTCounter(myInformationFrame, "[Triangles] Step 2 (checking #" + listOfEdges.size() + " edges)", listOfEdges.size()-2);
			for (i=0;i<listOfEdges.size()-2;i++){
			    eli = (Vector) listOfEdges.elementAt(i);
			    xi1 = ((Integer) eli.elementAt(0)).intValue();
			    xi2 = ((Integer) eli.elementAt(1)).intValue();
			    ci = ((Integer) eli.elementAt(2)).intValue();
			    for (j=i+1;j<listOfEdges.size()-1;j++){
				elj = (Vector) listOfEdges.elementAt(j);
				xj1 = ((Integer) elj.elementAt(0)).intValue();
				xj2 = ((Integer) elj.elementAt(1)).intValue();
				cj = ((Integer) elj.elementAt(2)).intValue();
				if ( (xi1 == xj1) || (xi1 == xj2) || (xi2 == xj1) || (xi2 == xj2) ){
				    if ( (xi1 == xj1) || (xi1 == xj2) ){
					cij = xi1;
					cia = xi2;
				    }else{
					cij = xi2;
					cia = xi1;
				    }
				    if (cij == xj1)
					cja = xj2;
				    else
					cja = xj1;
				    
				    for (k=j+1;k<listOfEdges.size();k++){
					elk = (Vector) listOfEdges.elementAt(k);
					xk1 = ((Integer) elk.elementAt(0)).intValue();
					xk2 = ((Integer) elk.elementAt(1)).intValue();		    
					ck = ((Integer) elk.elementAt(2)).intValue();
					
					cs = ci+cj+ck;
					if ( ( (xk1 == cia) && (xk2 == cja) )
					     || ( (xk2 == cia) && (xk1 == cja) ) ){
					    //System.out.println("Found triangle between cij = " + cij +", xk1 = " + xk1 + ", xk2 = " + xk2 + " (Sigma = " + cs + ")");
					    
					    if (cs == 3)
						allP++;
					    else if (cs == 1)
						oneN++;
					    else if (cs == -1)
						twoN++;
					    else if (cs == -3)
						allN++;
					}
				    }
				}
			    }
			    cc.increment();
			}
			cc.end();
			tot = allP + oneN + twoN + allN;
			saff += "Triangle consistency: found " + tot + " triangles, (all+) pattern = " + allP + ", (two-) pattern = " + twoN + " (one- pattern = " + oneN + ", all- pattern = " + allN + ")\n";

			myInformationFrame.setText("Checking Delaunay properties for P-Delaunay = " + Statistics.LIMIT_P_DELAUNAY + ":\n" + saff);
		    }else
			myInformationFrame.setText("No visible edge");
		}
	    };
	t.start();
    }

    public void nonThreadClustering(final String lineu){
	String line = lineu;
	int ii, rep = 1;
	if (Clustering.containsRepeat(line)){
	    rep = Clustering.getRepeatTimes(line);
	    line = Clustering.flushRepeatString(line);
	}
	
	for (ii=0;ii<rep;ii++)
	    if (Clustering.containsLoop(line)){
		int [] bds = Clustering.getBoundsNclusters(line);
		String scur;
		int i;
		for (i=bds[0];i<=bds[1];i++){
		    scur = Clustering.getReplacementString(i, line);
		    synchronized (this){
			performClustering(scur);
		    }
		}
	    }else
		performClustering(line);
    }
    
    public void clustering(final String lineu){

	Thread t = new Thread(){
		public void run(){
		    nonThreadClustering(lineu);
		}
	    };
	t.start();
    }

    public void clustersToDimension(){
	int i;
	if (allClusterings!=null){
	    for (i=0;i<allClusterings.size();i++)
		((Clustering) allClusterings.elementAt(i)).myClusteringAlgorithm.compute_dimension();
	}
    }

    public void performClustering(String line){
	nbClustering++;
	final Clustering newClustering;
	newClustering = new Clustering(this,nbClustering-1);
	newClustering.getOptions(line);
	myInformationFrame.setTextProgressBar("Clustering");
	myInformationFrame.setText("Processing clustering : " + newClustering);
	timerTickOn();
	
	newClustering.toClustering();
	newClustering.saveSpace();
	
	if (nbClustering == 1){
	    allClusterings = new Vector();
	    myTabbedPane.myPCAPane.currentClustering = 0;
	    myTabbedPane.myManifoldPane.currentClustering = 0;
	    myTabbedPane.myClusteringPane.activateSearch();
	}
	allClusterings.add(newClustering);
	myClusteringProfileFrame.addClusteringLF(newClustering, nbClustering-1);
	
	updateClusteringLF(JAGCTVisualizationPane.stringClustering(nbClustering-1));//(newClustering.getName() + "(" + newClustering.nclusters + ")");
	
	myInformationFrame.setText("Processing clustering : " + newClustering + "... done");
	myTabbedPane.myClusteringPane.listClusteringPane.setText(myTabbedPane.myClusteringPane.clusteringList());
	timerTickOff();    
    }

    public void updateClusteringLF(String s){
	if (nbClustering == 1){
	    myTabbedPane.myPCAPane.clusteringSelect.removeAllItems();
	    myTabbedPane.myManifoldPane.clusteringSelect.removeAllItems();
	    myTabbedPane.myCorrelationPane.clusteringSelect.removeAllItems();

	    myTabbedPane.myPCAPane.clusteringSelect.addItem(s);
	    myTabbedPane.myManifoldPane.clusteringSelect.addItem(s);
	    myTabbedPane.myCorrelationPane.clusteringSelect.addItem(s);

	    myTabbedPane.myPCAPane.clusteringSelect.setEnabled(true);
	    myTabbedPane.myManifoldPane.clusteringSelect.setEnabled(true);
	    myTabbedPane.myCorrelationPane.clusteringSelect.setEnabled(true);


	    myTabbedPane.myPCAPane.clusterChoice.removeAllItems();
	    myTabbedPane.myManifoldPane.clusterChoice.removeAllItems();

	    int i;
	    Clustering cc = (Clustering) allClusterings.elementAt(0);
	    String [] col = cc.myClusteringAlgorithm.clusterChoices;
	    for(i=0;i<col.length;i++){
		myTabbedPane.myPCAPane.clusterChoice.addItem(col[i]);
		myTabbedPane.myManifoldPane.clusterChoice.addItem(col[i]);
	    }

	    myTabbedPane.myPCAPane.clusterChoice.setEnabled(true);
	    myTabbedPane.myManifoldPane.clusterChoice.setEnabled(true);


	}else{
	    if (nbClustering < 1)
		Matrix.perror("0 clusterings to display !");
	    myTabbedPane.myPCAPane.clusteringSelect.addItem(s);
	    myTabbedPane.myManifoldPane.clusteringSelect.addItem(s);
	    myTabbedPane.myCorrelationPane.clusteringSelect.addItem(s);
	}

	myTabbedPane.myPCAPane.updateStatButtons();
	myTabbedPane.myManifoldPane.updateStatButtons();

	myClusteringFrame.updateClusteringLF();
    }

    public Clustering getClustering(int num){

	if (num < 0)
	    Matrix.perror("Negative value for clustering number !");
	if (num > allClusterings.size() - 1)
	    Matrix.perror("Clustering number out of range");

	return (Clustering) allClusterings.elementAt(num);
    }

}
