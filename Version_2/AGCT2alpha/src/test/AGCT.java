package test;

import cellDesignerPlugin.CellDesignerInput;
import clusteringProfile.ClusteringProfileFile;
import clusteringProfile.JClusteringProfileFrame;
import combinationRun.IterativeRun;
import combinationRun.IterativeRun.SortMode;
import filtering.JLeftIndicator;
import forDebug.Debug;
import garudaBinding.AGCTPropertyChangeListener;
import garudaBinding.AGCTgaruda;
import gene.Gene;
import gene.GeneList;
import gene.GeneListSaver;
import gui.ClosingDialog;
import gui.JNormalizationMenu;
import ligand.Ligand;
import matrix.MatrixConfig;
import matrix.MatrixMath;
import normalizationData.NormalizationData;
import timer.Timer;
import timer.Timer.TimerKey;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;

public class AGCT extends JApplet implements Runnable, Debuggable, ActionListener {

    public static String SUFFIX = System.currentTimeMillis() + "";
    private Timer timer = new Timer();

    public Timer getTimer() {
        return timer;
    }

    /**
     * ファイルをロードしたら、デフォルトのセッティングで、自動的に最後まで解析を行う。終了したら、時間を表示する
     */
    public static boolean AUTO = true;
    private static AGCT singleton;

    public static AGCT getInstance() {
        return singleton;
    }

    public static boolean highlightOther = true;
    private final NormalizationData normalizationData = new NormalizationData();// 2回目の呼び出し？
    public static int HighlightMask = 0;

    public NormalizationData getNormalizationData() {
        return normalizationData;
    }

    public static final boolean MYDEBUG = true;
    public static final boolean LOAD_NEW_ANNOTATION_FILE = true;

    public static void debug(Object... os) {
        System.err.println(Arrays.deepToString(os));
    }

    /**
     * *******************************************************************************
     * Static stuff
     * ***
     */
    public static final String Token_Ordered_List_Names_Begin = "@Ordered_List_Names_Begin", Token_Ordered_List_Names_End = "@Ordered_List_Names_End",
            Token_Matrix_W_Begin = "@Matrix_W_Begin", Token_Matrix_W_End = "@Matrix_W_End", Token_Manifold_Eigenvalues_Begin = "@Manifold_Eigenvalues_Begin",
            Token_Manifold_Eigenvalues_End = "@Manifold_Eigenvalues_End", Token_Matrix_M_Begin = "@Matrix_M_Begin", Token_Matrix_M_End = "@Matrix_M_End";
    /**
     * *******************************************************************************
     * Static stuff
     * ***
     */
    public static boolean Loading_From_Scenario_Manifold_Begin, Loading_From_Scenario_Manifold_End;
    public static boolean MULTIPLE_ROWS_FOR_A_SAME_GENE_CONSIDERED_DIFFERENT_GENES = false;
    public static boolean MULTIPLE_ROWS_FOR_A_SAME_GENE_DELETED_EXCEPT_FIRST_ONE = false;
    public static boolean CHECK_TRIANGLE_CONSISTENCY = true;
    public static final boolean HAAR_WAVELETS_REMOVE_CONSTANT = true;
    // we remove the product wih the constant function
    public static Random RandomGenerator;
    // global random number generator
    public static final double Very_Small = 10E-6;
    public static int DEFAULT_Max_Number_Genes = 10;
    public static long DEFAULT_RANDOMSEED = new Random(124124).nextLong();
    public static final String MS = "MS", MF = "MF", MW = "MW", MN = "MN";
    public static final String[] Feature_Method = {"Slopes", "Haar_Wavelets", "Daubechies_D4_Wavelets_Interval", "Daubechies_D4_Wavelets", "Daubechies_D8_Wavelets",
            "Daubechies_D10_Wavelets", "Daubechies_D12_Wavelets", "Daubechies_D14_Wavelets", "Daubechies_D16_Wavelets", "Daubechies_D18_Wavelets",
            "Daubechies_D20_Wavelets"};
    public static final String[] Feature_Selection_Method = {"None", "K_Smallest_Average_Determination_Coefficient"};
    public static int Number_Of_Clusters, Number_Of_Neighbors = 10, WindowWidth = 800;
    public static int Method_F = 0;
    /**
     * // computes features // 0 : use slopes // 1 : use Haar wavelet
     * coefficients // 2 : use Daubechies D4 wavelets on the interval (Cohen,
     * Daubechies & Vial, // 1993) // 3 : use Daubechies D4 // 4 : use
     * Daubechies D8 // 5 : use Daubechies D10 // 6 : use Daubechies D12 // 7 :
     * use Daubechies D14 // 8 : use Daubechies D16 // 9 : use Daubechies D18 //
     * 10 : use Daubechies D20
     */
    public static int Method_FS = 0;
    /**
     * // feature selection methods // 0 : no selection // 1 : uses average
     * determination coefficient (repeatedly removes the // variable with the
     * largest value)
     */
    public static int Method_PS = 0;
    /**
     * // prototype selection methods // 0 : no selection // 1 : uses k-means
     * (runs k-means, and then keeps the k prototypes that are // closest to
     * their cluster center)
     */
    public static int Max_Number_Of_Features = -1;
    /**
     * // maximum number of features authorized (useful e.g. when using wavelets
     * // with large number of coefficients) // if -1, uses DEFAULT
     */
    public static final int DEFAULT_Max_Number_Of_Features = 50;
    public static int Max_Number_Of_Prototypes = -1;
    /**
     * // maximum number of prototypes authorized (useful e.g. when extremely
     * large // number of prototypes) // if -1, uses DEFAULT
     */
    public static final int DEFAULT_Max_Number_Of_Prototypes = 1000;
    public static int Number_Of_Wavelet_Stamps = -1;
    /**
     * // number of time stamps used to compute the wavelet decomposition //
     * must be a power of 2, or -1 (the program chooses at best)
     */
    public static final int DEFAULT_Number_Of_Wavelet_Stamps = 256;
    /**
     * // computes similarity // 0 : Heat Kernel // 1 : Cosine distance // 2 :
     * Absolute Cosine distance (= absolute correlation)
     */
    public static int Method_S = 1;
    /**
     * // computes similarity matrix W // 0 : keep similarity as is // 1 :
     * filter = use Symmetric NN (very sparse matrix) // 2 : filter = use
     * Natural Neighbors (very sparse matrix)
     */
    public static int Method_W = 1;
    /**
     * // computes normalized matrix N // 0 : N = D^{-1/2} W D^{-1/2} // 1 : N =
     * D^{-1} W // 2 : N = 2*Stochastic(W)
     * 3 : N = 2*Stochastic(DS Scaling(W)) // Doubly_Stochastic_Scaling: // When true, the matrix is
     * first multiplied by a real T so as to be the // closest to // a doubly
     * stochastic matrix. This does not change eigenvectors, but // greatly
     * reduces // the doubly stochastic approximation time // When used as
     * similarity matrix, results become more "flattened"
     */
    public static int Method_N = 1;
    public static boolean Sort_Depth, No_GeneSelection, Use_Shadow, DoDebug,
            Warning, Perspective, Only_Referenced, Only_ReferencedEdges, Referenced_Available;
    public static int Max_Reference = -1;
    /**
     * // -1 -> no discrimination // 0 -> show neg // 1 -> show pos
     */
    public static int Show_Correlation = -1;
    /**
     * // 0 -> no discrimination // 1 -> is triangulation computed, shows only
     * genes with at least one // visible edge
     */
    public static boolean Show_Only_Genes_With_Visible_Edge = false;
    /**
     * // T parameter for the heat kernel;
     */
    public static double T_Heat_Kernel = 2.0;
    /**
     * // T_Heat_Kernel used to compute similarities on the *manifold* (2D)
     */
    public static double T_Heat_Kernel_Similarity_Neighbors = 0.1;
    /**
     * // p*100 % of the smallest similarities in W are replaced by 0
     */
    public static double Sparsify_Statistic = 0.01;
    /**
     * // number of components on the manifold to keep for each individual
     */
    public static int Number_Of_Manifold_Components = 20;
    /**
     * // number of dimensions to consider for triangulation
     */
    public static int Number_Of_Triangulation_Dimensions = 3;
    JNormalizationMenu jNormalizationMenu;
    /**
     * GarudaOptionMenu (add 20130122)
     */
    /**
     * **************************************
     * garuda section
     * (JMenu is implemented together with the others)
     * 2013.2.20 garuda9
     * coded by TM
     * <p/>
     * <p/>
     * ***************************************
     */
    private AGCTPropertyChangeListener apcl;

    /**
     *
     */
    private void startGaruda(AGCT singleton) {
        apcl = new AGCTPropertyChangeListener(singleton);
        System.out.println("connected to garuda");
    }

    ////////////////// END OF GARUDA //////////////////////////////////


    /**
     * @param ex
     * @param ap
     */
    public static void executeManifold(String ex, AGCT ap) {
        System.out.println("AGCT.executeManifold()");
        StringTokenizer t;
        String s1, s2;
        int i;

        if (!ControlProcess.hasTrue("filtersSelected")) {
            Matrix.perror("AGCT.class :: genes not computed");
        }
        if (Scenario.bOrdered_List_Names) {
            //checkしているだけ
            t = new StringTokenizer(ex, Scenario.myToken);
            if (t.countTokens() != ap.data.getMyDomain().numberSelectedGenes) {
                Matrix.perror("AGCT.class :: bad number of tokens for gene list (" + t.countTokens() + " != " + Arrays.toString(ap.data.getMyDomain().selectedGeneNumberToGeneNumber) + ")");
            }
            for (i = 0; i < ap.data.getMyDomain().numberSelectedGenes; i++) {
                s1 = ap.data.getMyDomain().getGenes().get(ap.data.getMyDomain().selectedGeneNumberToGeneNumber[i]).getName();
                s2 = t.nextToken();

                // System.out.println("Bluk names : gene = " + s1 +
                // ", fileName = " + s2);

                if (!s1.equals(s2)) {
                    Matrix.perror("AGCT.class :: bad name correspondence for gene (" + s1 + " != " + s2 + ")");
                }
            }
        } else if (Scenario.bMatrix_W) {
            // 0でない部分に関して、id|value という形式で書かれている
            if (ap.data.getCurrentIndexW() == -1) {
                ap.data.setW(new Matrix("W", ap.data.getMyDomain().numberSelectedGenes, ap.data.getMyDomain().numberSelectedGenes));
                ap.data.setCW(new Matrix("W_Original", ap.data.getMyDomain().numberSelectedGenes, ap.data.getMyDomain().numberSelectedGenes));
                ap.data.setCurrentIndexW(ap.data.getCurrentIndexW() + 1);
            }
            t = new StringTokenizer(ex, Scenario.myToken);
            while (t.hasMoreTokens()) {
                int j = Integer.valueOf(t.nextToken());
                double val = Double.valueOf(t.nextToken());
                ap.data.getW().set(ap.data.getCurrentIndexW(), j, val);
                ap.data.getCopyOfSimilarityMatrix().set(ap.data.getCurrentIndexW(), j, val);
            }
            ap.data.setCurrentIndexW(ap.data.getCurrentIndexW() + 1);
        } else if (Scenario.bManifold_Eigenvalues) {
            t = new StringTokenizer(ex, Scenario.myToken);
            if (t.countTokens() != ap.data.getMyDomain().numberSelectedGenes) {
                Matrix.perror("AGCT.class :: bad number of tokens for manifold eigenvalues (" + t.countTokens() + " != "
                        + ap.data.getMyDomain().selectedGeneNumberToGeneNumber + ")");
            }
            for (i = 0; i < ap.data.getMyDomain().numberSelectedGenes; i++) {
                if (i == 0) {
                    ap.data.setManifold_eigenvalues(new double[ap.data.getMyDomain().numberSelectedGenes]);
                }
                ap.data.getManifold_eigenvalues()[i] = Double.parseDouble(t.nextToken());
            }
        } else if (Scenario.bMatrix_M) {
            if (ap.data.getCurrentIndexM() == -1) {
                ap.data.setM(new Matrix("Manifold Eigenvectors", ap.data.getMyDomain().numberSelectedGenes, ap.data.getMyDomain().numberSelectedGenes));
                ap.generate_Neighbors();
                ControlProcess.put("neighborsComputed", true);
                ap.data.setCurrentIndexM(ap.data.getCurrentIndexM() + 1);
            }
            t = new StringTokenizer(ex, Scenario.myToken);
            while (t.hasMoreTokens()) {
                int j = Integer.valueOf(t.nextToken());
                double val = Double.valueOf(t.nextToken());
                ap.data.getM().set(ap.data.getCurrentIndexM(), j, val);
            }
            ap.data.setCurrentIndexM(ap.data.getCurrentIndexM() + 1);
        } else {
            Matrix.perror("AGCT.class :: not in a manifold scenario");
        }
    }

    public static void initDefaultClassVariables() {
        Number_Of_Neighbors = 10;
        Number_Of_Wavelet_Stamps = -1;
        T_Heat_Kernel = 2.0;
        T_Heat_Kernel_Similarity_Neighbors = 0.1;
        Sparsify_Statistic = 0.01;
        Number_Of_Manifold_Components = 20;
        Number_Of_Triangulation_Dimensions = 3;
        Max_Number_Of_Features = -1;
    }

    public static int Max_Number_Of_Genes_Kept = 10;
    /**
     * // maximum number of genes to prevent storage capacity problems in
     * Strings
     */
    /**************************************************************************************
     * Graphics part
     *****/
    /**
     * // Buttons and fields that appear on the top buttonbar
     */
    JButton loadButton, saveButton, saveManifold, showAxis, performReflection, saveButtonStats, saveButtonDelaunay, saveButtonWeb, binButton, saveClusteringButton,
            cameraButton,
            scenarioLoadButton, scenarioSaveButton, scenarioDeleteButton, scenarioRecordButton, scenarioDisplayButton;
    ImageIcon scenarioRecordIdleIcon, scenarioRecordingIcon, scenarioLoadIdleIcon, scenarioRunningIcon;
    JLabel clockLabel;
    /**
     * items and alike that appear on the menus
     */
    JCheckBoxMenuItem useDebug, useWarning, sortDepth, useShadow, perspective, highLightGene, onlyReferencedEdges, noGeneSelection;
    JRadioButtonMenuItem methodF_0, methodF_1, methodF_2, methodF_3, methodF_4, methodF_5, methodF_6, methodF_7, methodF_8, methodF_9, methodF_10, methodS_0, methodS_1,
            methodS_2, methodW_0, methodW_1, methodW_2, methodN_0, methodN_1, methodN_2, methodN_3;
    JMenuItem numberOfNeighbors, numberOfManifoldComponents, numberOfTriangulationDimensions, tHeatKernel, tHeatKernelSimilarityNeighbors, sparsifyStatistic, limitp,
            bregDiv, limitpdelaunay, limitpmagnitude, numberOfProfiles, loadHighlightFile, randomSeed;
    JMenuItem goFiltering;
    public JMenuItem loadResponsiveFile;
    /**
     * "pointer" on the master frame
     */
    public JFrame myFrame;
    /**
     * // "pointer" on the JinformationFrame uses to display all informations
     * about // processing
     */
    // public JInformationFrame JInformationFrame.getInstance();
    /**
     * id but for the Annotation Frame
     */
    JAnnotationFrame myAnnotationFrame;
    /**
     * id but to display global informations about clusterings
     */
    JClusteringFrame myClusteringFrame;
    /**
     * id but to display cluster profiles
     */

    public JClusteringProfileFrame myClusteringProfileFrame;
    /**
     * tabbedPane
     */
    public JMainFrameTabbedPane myTabbedPane;
    /**
     * DSScaling
     */
    JDSScaling myScaling;
    public AGCTData data = new AGCTData(false, false, false, false, false, false);

    /**
     * yes iff a hard membership is computed
     */
    private AGCT(JBugFreeFrame myF) {
        super();
        if (AGCT.MYDEBUG) {
            AGCT.debug("AGCT.AGCT(jBugFreeFrame)");
        }
        AGCT.Loading_From_Scenario_Manifold_Begin = AGCT.Loading_From_Scenario_Manifold_End = false;

        Referenced_Available = false;
        myFrame = myF;
        myFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        myFrame.addWindowListener(new WindowAdapter() {

            public void windowClosing(WindowEvent ev) {
                final Frame frame = (Frame) ev.getSource();
                new ClosingDialog(frame);
            }
        });
        AGCT.RandomGenerator = new Random(AGCT.DEFAULT_RANDOMSEED);
        AGCTClustering_CP.init();

        data.setMyAGCTFileWriter(new AGCTFileWriter(this));
        myAnnotationFrame = null;
        // Here Clustering Reslut Vis... starting.
        myClusteringFrame = new JClusteringFrame("AGCT --- Clustering Result Visualization Frame", this);
        myClusteringProfileFrame = (new JClusteringProfileFrame("AGCT --- [Clustering|Gene] Profile Visualization Frame", this));// scroll
        myScaling = new JDSScaling(this);

        data.setMin_CorrelationPnt3D(new Pnt3D());
        data.setMax_CorrelationPnt3D(new Pnt3D());

        data.setNbClustering(0);
        data.setAllClusterings(null);
        data.setCurrentIndexW(-1);
        data.setCurrentIndexM(-1);
        Statistics.flushAll();
    }

    public Gene getDirectGene(int iGene) {
        return data.getMyDomain().getGenes().get(iGene);
    }

    public static void setCluster(int v) {
        Number_Of_Clusters = v;
    }

    public static void setSort_Depth(boolean b) {
        Sort_Depth = b;
    }

    public static void setNo_GeneSelection(boolean b) {
        No_GeneSelection = b;
    }

    public static void setUse_Shadow(boolean b) {
        Use_Shadow = b;
    }

    public static void setPerspective(boolean b) {
        Perspective = b;
    }

    public static void setOnly_Referenced(boolean b) {
        Only_Referenced = b;
    }

    public static void setOnly_ReferencedEdges(boolean b) {
        Only_ReferencedEdges = b;
    }

    public static void setDebug(boolean b) {
        DoDebug = b;
    }

    public static void setWarning(boolean b) {
        Warning = b;
    }

    public static void setNumber_Of_Neighbors(int n) {
        Number_Of_Neighbors = n;
    }

    public static void setNumber_Of_Manifold_Components(int n) {
        Number_Of_Manifold_Components = n;
    }

    public static void setNumber_Of_Wavelet_Stamps(int n) {
        System.out.println("AGCT.setNumber_Of_Wavelet_Stamps(" + n + ")");
        Number_Of_Wavelet_Stamps = n;
    }

    public static void setMax_Number_Of_Features(int n) {
        System.out.println("AGCT.setMax_Number_Of_Features(" + n + ")");
        Max_Number_Of_Features = n;
    }

    public static void setMax_Number_Of_Prototypes(int n) {
        System.out.println("AGCT.setMax_Number_Of_Prototypes(" + n + ")");
        Max_Number_Of_Prototypes = n;
    }

    public static void setNumber_Of_Profiles(int n) {
        DEFAULT_Max_Number_Genes = n;
    }

    public static void setRandomSeed(int n) {
        DEFAULT_RANDOMSEED = n;
    }

    public static void setNumber_Of_Triangulation_Dimensions(int n) {
        Number_Of_Triangulation_Dimensions = n;
    }

    public void setBregDiv(int n) {
        myScaling.indexSelected = n;
    }

    public static void setT_Heat_Kernel(double n) {
        T_Heat_Kernel = n;
    }

    public static void setT_Heat_Kernel_Similarity_Neighbors(double n) {
        T_Heat_Kernel_Similarity_Neighbors = n;
    }

    public static void setSparsify_Statistic(double n) {
        Sparsify_Statistic = n;
    }

    public static void setLimit_P_Chi2(double p) {
        Statistics.LIMIT_P_CHI2 = p;
    }

    public static void setLimit_P_Magnitude(double p) {
        Statistics.LIMIT_P_MAGNITUDE = p;
    }

    public static void setMethod_F(int v) {
        Method_F = v;
    }

    public static void setMethod_S(int v) {
        Method_S = v;
    }

    public static void setMethod_W(int v) {
        Method_W = v;
    }

    public static void setMethod_N(int v) {
        Method_N = v;
    }

    public void setMyFrame(JFrame f) {
        myFrame = f;
    }


    public void requestModificationNumber_Of_Matrx_Max_Iteration() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Max # of iterations on eigenvector calculation", Integer.toString(MatrixMath.getNumericalIterationsMax()));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = MatrixMath.getNumericalIterationsMax();
            }
            if (i > 0) {
                MatrixMath.setNumerical_Iterations_Max(i);
            }
        }
    }

    public void requestModificationMethod_F(int v) {
        if (v == 0) {
            methodF_0.setSelected(true);
        } else if (v == 1) {
            methodF_1.setSelected(true);
        } else if (v == 2) {
            methodF_2.setSelected(true);
        } else if (v == 3) {
            methodF_3.setSelected(true);
        } else if (v == 4) {
            methodF_4.setSelected(true);
        } else if (v == 5) {
            methodF_5.setSelected(true);
        } else if (v == 6) {
            methodF_6.setSelected(true);
        } else if (v == 7) {
            methodF_7.setSelected(true);
        } else if (v == 8) {
            methodF_8.setSelected(true);
        } else if (v == 9) {
            methodF_9.setSelected(true);
        } else if (v == 10) {
            methodF_10.setSelected(true);
        }
        AGCT.setMethod_F(v);
    }

    public void requestMofificationMethod_F() {
        if (!(methodF_0.isSelected() || methodF_1.isSelected() || methodF_2.isSelected() || methodF_3.isSelected() || methodF_4.isSelected() || methodF_5.isSelected()
                || methodF_6.isSelected() || methodF_7.isSelected() || methodF_8.isSelected() || methodF_9.isSelected() || methodF_10.isSelected())) {
            Matrix.perror("No radio button for Features selected");
        }
        if (methodF_0.isSelected()) {
            AGCT.setMethod_F(0);
        }
        if (methodF_1.isSelected()) {
            AGCT.setMethod_F(1);
        }
        if (methodF_2.isSelected()) {
            AGCT.setMethod_F(2);
        }
        if (methodF_3.isSelected()) {
            AGCT.setMethod_F(3);
        }
        if (methodF_4.isSelected()) {
            AGCT.setMethod_F(4);
        }
        if (methodF_5.isSelected()) {
            AGCT.setMethod_F(5);
        }
        if (methodF_6.isSelected()) {
            AGCT.setMethod_F(6);
        }
        if (methodF_7.isSelected()) {
            AGCT.setMethod_F(7);
        }
        if (methodF_8.isSelected()) {
            AGCT.setMethod_F(8);
        }
        if (methodF_9.isSelected()) {
            AGCT.setMethod_F(9);
        }
        if (methodF_10.isSelected()) {
            AGCT.setMethod_F(10);
        }

        Scenario.add("AGCT_Modification_Method_F", AGCT.Method_F);
    }

    public void requestModificationMethod_S(int v) {
        if (v == 0) {
            methodS_0.setSelected(true);
        } else if (v == 1) {
            methodS_1.setSelected(true);
        } else if (v == 2) {
            methodS_2.setSelected(true);
        }
        AGCT.setMethod_S(v);
    }

    public void requestModificationMethod_S() {
        if (!(methodS_0.isSelected() || methodS_1.isSelected() || methodS_2.isSelected())) {
            Matrix.perror("No radio button for Similarity selected");
        }
        if (methodS_0.isSelected()) {
            AGCT.setMethod_S(0);
        }
        if (methodS_1.isSelected()) {
            AGCT.setMethod_S(1);
        }
        if (methodS_2.isSelected()) {
            AGCT.setMethod_S(2);
        }

        Scenario.add("AGCT_Modification_Method_S", AGCT.Method_S);
    }

    public void requestModificationMethod_W(int v) {
        if (v == 0) {
            methodW_0.setSelected(true);
        } else if (v == 1) {
            methodW_1.setSelected(true);
        }
        AGCT.setMethod_W(v);
    }

    public void requestMofificationMethod_W() {
        if (!(methodW_0.isSelected() || methodW_1.isSelected())) {
            Matrix.perror("No radio button for Matrix W selected");
        }
        if (methodW_0.isSelected()) {
            AGCT.setMethod_W(0);
        }
        if (methodW_1.isSelected()) {
            AGCT.setMethod_W(1);
        }

        Scenario.add("AGCT_Modification_Method_W", AGCT.Method_W);
    }

    public void requestModificationMethod_N(int v) {
        if (v == 0) {
            methodN_0.setSelected(true);
        } else if (v == 1) {
            methodN_1.setSelected(true);
        } else if (v == 2) {
            methodN_2.setSelected(true);
        } else if (v == 3) {
            methodN_3.setSelected(true);
        }
        AGCT.setMethod_N(v);
    }

    public void requestMofificationMethod_N() {
        if (!(methodN_0.isSelected() || methodN_1.isSelected() || methodN_2.isSelected() || methodN_3.isSelected())) {
            Matrix.perror("No radio button for Matrix N selected");
        }
        if (methodN_0.isSelected()) {
            AGCT.setMethod_N(0);
        }
        if (methodN_1.isSelected()) {
            AGCT.setMethod_N(1);
        }
        if (methodN_2.isSelected()) {
            AGCT.setMethod_N(2);
        }
        if (methodN_3.isSelected()) {
            AGCT.setMethod_N(3);
        }

        Scenario.add("AGCT_Modification_Method_N", AGCT.Method_N);
    }

    public void requestModificationNumber_Of_Neighbors(int v) {
        AGCT.setNumber_Of_Neighbors(v);
    }

    private void requestModificationNumber_Of_Neighbors() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Number of Neighbors\n(Useful when SNN for W chosen)", Integer.toString(AGCT.Number_Of_Neighbors));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = 5;
            }
            if (i > 0) {
                AGCT.setNumber_Of_Neighbors(i);
                Scenario.add("AGCT_Modification_Number_Of_Neighbors", i);
            }
        }
    }

    public void requestModificationNumber_Of_Profiles() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Max # of profiles displayed\n(in clustering profile frame)", Integer.toString(AGCT.DEFAULT_Max_Number_Genes));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = 10;
            }
            if (i > 0) {
                AGCT.setNumber_Of_Profiles(i);
            }
        }
    }

    public void requestModificationRandomSeed() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Random seed:", Long.toString(AGCT.DEFAULT_RANDOMSEED));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = 10;
            }
            if (i > 0) {
                AGCT.setRandomSeed(i);
            }
        }
    }

    public void requestModificationMax_Number_Of_Features(int v) {
        AGCT.setMax_Number_Of_Features(v);
    }

    public void requestModificationMax_Number_Of_Prototypes(int v) {
        AGCT.setMax_Number_Of_Prototypes(v);
    }

    public void requestModificationMax_Number_Of_Features() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Maximum Number of Final Features by Gene:\n* -1 implies no selection, \n* or > 2", Integer.toString(AGCT.Max_Number_Of_Features));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = -1;
            }
            if (i > 2 || i == -1) {
                AGCT.setMax_Number_Of_Features(i);
                Scenario.add("JAGCTSelectionPane_Max_Number_Of_Features", i);
            }
        }
    }

    public void requestModificationMax_Number_Of_Prototypes() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Maximum Number of Prototypes:\n* -1 implies no selection, \n* or > 2", Integer.toString(AGCT.Max_Number_Of_Prototypes));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = -1;
            }
            if (i > 2 || i == -1) {
                AGCT.setMax_Number_Of_Prototypes(i);
                Scenario.add("JAGCTSelectionPane_Max_Number_Of_Prototypes", i);
            }
        }
    }

    public void defaultNWS() {
        requestModificationNumber_Of_Wavelet_Stamps(AGCT.DEFAULT_Number_Of_Wavelet_Stamps);
    }

    public void requestModificationNumber_Of_Wavelet_Stamps(int v) {
        AGCT.setNumber_Of_Wavelet_Stamps(v);
    }

    public void requestModificationNumber_Of_Wavelet_Stamps() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Number of Wavelet Stamps:\n* -1 to let the program choose, \n* or a power of 2 > 1", Integer.toString(AGCT.Number_Of_Wavelet_Stamps));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = -1;
            }
            if (i > 1 && Util_Feature.powerOfTwo(i) || i == -1 && !Util_Feature.checkAllWaveletStampsAutomatic(this)) {
                AGCT.setNumber_Of_Wavelet_Stamps(i);
                Scenario.add("AGCT_Modification_Number_Of_Wavelet_Stamps", i);
            } else if (!Util_Feature.powerOfTwo(i) || i == -1) {
                defaultNWS();
                Scenario.add("AGCT_Modification_Number_Of_Wavelet_Stamps", AGCT.Number_Of_Wavelet_Stamps);
            }
        }
    }

    public void requestModificationNumber_Of_Manifold_Components(int v) {
        AGCT.setNumber_Of_Manifold_Components(v);
    }

    public void requestModificationNumber_Of_Manifold_Components() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Number of Manifold Components to keep", Integer.toString(AGCT.Number_Of_Manifold_Components));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = 5;
            }
            if (i > 2) {
                AGCT.setNumber_Of_Manifold_Components(i);
                Scenario.add("AGCT_Modification_Number_Of_Manifold_Components", i);
            }
        }
    }

    public void requestModificationBregDiv(int i) {
        myScaling.indexSelected = i;
    }

    public void requestModificationBregDiv() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Bregman divergence for DS Scaling\n" + JDSScaling.getString(), Integer.toString(myScaling.indexSelected));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = 0;
            }
            if (JDSScaling.isValid(i)) {
                myScaling.indexSelected = i;
                Scenario.add("JDSScaling_Modification_Bregman_Divergence", i);
            }
        }
    }

    public void requestModificationNumber_Of_Triangulation_Dimensions(int v) {
        AGCT.setNumber_Of_Triangulation_Dimensions(v);
    }

    public void requestModificationNumber_Of_Triangulation_Dimensions() {
        final JOptionPane d = new JOptionPane();
        final String ret = d.showInputDialog(this, "Which dimension for the triangulation ?\n(must be <= #Manifold components)", Integer.toString(AGCT.Number_Of_Triangulation_Dimensions));
        int i;

        if (ret != null) {
            try {
                i = Integer.parseInt(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not an integer !");
                i = 3;
            }
            if (i > 2 && i <= AGCT.Number_Of_Manifold_Components) {
                AGCT.setNumber_Of_Triangulation_Dimensions(i);
                Scenario.add("AGCT_Modification_Number_Of_Triangulation_Dimensions", i);
            }
        }
    }

    public void requestModificationT_Heat_Kernel(double v) {
        AGCT.setT_Heat_Kernel(v);
    }

    public void requestModificationT_Heat_Kernel() {
        final JOptionPane d = new JOptionPane();
        double i;
        final String ret = d.showInputDialog(this, "T Heat Kernel\n(W)", Double.toString(AGCT.T_Heat_Kernel));

        if (ret != null) {
            try {
                i = Double.parseDouble(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not a double !");
                i = 2.0;
            }
            if (i > 0) {
                AGCT.setT_Heat_Kernel(i);
                Scenario.add("AGCT_Modification_T_Heat_Kernel", i);
            }
        }
    }

    public void requestModificationT_Heat_Kernel_Similarity_Neighbors(double v) {
        AGCT.setT_Heat_Kernel_Similarity_Neighbors(v);
    }

    public void requestModificationT_Heat_Kernel_Similarity_Neighbors() {
        final JOptionPane d = new JOptionPane();
        double i;
        final String ret = d.showInputDialog(this, "T Heat Kernel\n(2D manifold display)", Double.toString(AGCT.T_Heat_Kernel_Similarity_Neighbors));

        if (ret != null) {
            try {
                i = Double.parseDouble(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not a double !");
                i = 2.0;
            }
            if (i > 0) {
                AGCT.setT_Heat_Kernel_Similarity_Neighbors(i);
                Scenario.add("AGCT_Modification_T_Heat_Kernel_Similarity_Neighbors", i);
            }
        }
    }

    public void requestModificationSparsify_Statistic(double v) {
        AGCT.setSparsify_Statistic(v);
    }

    public void requestModificationSparsify_Statistic() {
        final JOptionPane d = new JOptionPane();
        double i;
        final String ret = d.showInputDialog(this, "Sparsify Statistic\n(sparsify 100p% entries of W)", Double.toString(AGCT.Sparsify_Statistic));

        if (ret != null) {
            try {
                i = Double.parseDouble(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not a double !");
                i = 2.0;
            }
            if (i >= 0.0 && i <= 1.0) {
                AGCT.setSparsify_Statistic(i);
                Scenario.add("AGCT_Modification_Sparsify_Statistic", i);
            }
        }
    }

    public void requestModificationLimit_P_Chi2(double v) {
        AGCT.setLimit_P_Chi2(v);
        if (myClusteringFrame != null) {
            myClusteringFrame.graphPanel.update = true;
        }
    }

    public void requestModificationLimit_P_Delaunay(double v) {
        Statistics.setLIMIT_P_DELAUNAY(v);
        if (myClusteringFrame != null) {
            myClusteringFrame.graphPanel.update = true;
        }
    }

    public void requestModificationLimit_P_Chi2() {
        final JOptionPane d = new JOptionPane();
        double i;
        final String ret = d.showInputDialog(this, "Limit P-value for Chi2\n(typically < 0.05)", Double.toString(Statistics.LIMIT_P_CHI2));

        if (ret != null) {
            try {
                i = Double.parseDouble(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not a double !");
                i = 0.01;
            }
            if (i < 0.0) {
                i = 0.0;
            }
            if (i > 1.0) {
                i = 1.0;
            }

            if (i > 0.0) {
                AGCT.setLimit_P_Chi2(i);
                Scenario.add("Statistics_Modification_Limit_P", i);
            }

            if (myClusteringFrame != null) {
                myClusteringFrame.graphPanel.update = true;
            }
        }
    }

    public void requestModificationLimit_P_Delaunay() {
        final JOptionPane d = new JOptionPane();
        double i;
        final String ret = d.showInputDialog(this, "Limit P-value for plotting Delaunay triangulation\n(not too small)", Double.toString(Statistics.getLIMIT_P_DELAUNAY()));

        if (ret != null) {
            try {
                i = Double.parseDouble(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not a double !");
                i = 0.2;
            }
            if (i < 0.0) {
                i = 0.0;
            }
            if (i > 1.0) {
                i = 1.0;
            }

            if (i > 0.0) {
                Statistics.setLIMIT_P_DELAUNAY(i);
            }

            if (myClusteringFrame != null) {
                myClusteringFrame.graphPanel.update = true;
            }
        }
    }

    public void requestModificationLimit_P_Magnitude() {
        final JOptionPane d = new JOptionPane();
        double i;
        final String ret = d.showInputDialog(this, "Limit P-value for TotalDistance\n", Double.toString(Statistics.LIMIT_P_MAGNITUDE));

        if (ret != null) {
            try {
                i = Double.parseDouble(ret);
            } catch (final NumberFormatException n) {
                System.out.println("Not a double !");
                i = 0.2;
            }
            if (i < 0.0) {
                i = 0.0;
            }
            if (i > 1.0) {
                i = 1.0;
            }

            if (i > 0.0) {
                AGCT.setLimit_P_Magnitude(i);
            }

            if (myClusteringFrame != null) {
                myClusteringFrame.graphPanel.update = true;
            }
        }
    }

    public void requestModificationSort_Depth(boolean v) {
        sortDepth.setSelected(v);
        requestModificationSort_Depth();
    }

    public void requestModificationSort_Depth() {
        AGCT.setSort_Depth(sortDepth.isSelected());
        Scenario.add("AGCT_Modification_Sort_Depth", sortDepth.isSelected());
        repaint();
    }

    public void requestModificationNo_GeneSelection(boolean v) {
        noGeneSelection.setSelected(v);
        requestModificationNo_GeneSelection();
    }

    public void requestModificationNo_GeneSelection() {
        AGCT.setNo_GeneSelection(noGeneSelection.isSelected());
    }

    public void requestModificationUse_Shadow(boolean v) {
        useShadow.setSelected(v);
        requestModificationUse_Shadow();
    }

    public void requestModificationUse_Shadow() {
        AGCT.setUse_Shadow(useShadow.isSelected());
        Scenario.add("AGCT_Modification_Use_Shadow", useShadow.isSelected());
        repaint();
    }

    public void requestModificationPerspective(boolean v) {
        perspective.setSelected(v);
        requestModificationPerspective();
    }

    public void requestModificationOnly_Referenced(boolean v) {
        highLightGene.setSelected(v);
        requestModificationOnly_Referenced();
    }

    public void requestModificationOnly_ReferencedEdges(boolean v) {
        onlyReferencedEdges.setSelected(v);
        requestModificationOnly_ReferencedEdges();
    }

    public void requestModificationPerspective() {
        AGCT.setPerspective(perspective.isSelected());
        Scenario.add("AGCT_Modification_Perspective", perspective.isSelected());
        myTabbedPane.myManifoldPane.visualizationPanel.myView3D.setPerspective(AGCT.Perspective);
        myTabbedPane.myPCAPane.visualizationPanel.myView3D.setPerspective(AGCT.Perspective);
        myTabbedPane.myCorrelationPane.visualizationPanel.myView3D.setPerspective(AGCT.Perspective);
        repaint();
    }

    public void requestModificationOnly_Referenced() {
        AGCT.setOnly_Referenced(highLightGene.isSelected());
        myTabbedPane.myManifoldPane.visualizationPanel.myView3D.setHighlight(AGCT.Only_Referenced);
        myTabbedPane.myPCAPane.visualizationPanel.myView3D.setHighlight(AGCT.Only_Referenced);
        myTabbedPane.myCorrelationPane.visualizationPanel.myView3D.setHighlight(AGCT.Only_Referenced);
        repaint();
    }

    public void requestModificationOnly_ReferencedEdges() {
        AGCT.setOnly_ReferencedEdges(onlyReferencedEdges.isSelected());
        myTabbedPane.myManifoldPane.visualizationPanel.myView3D.setOnlyReferencedEdges(AGCT.Only_ReferencedEdges);
        myTabbedPane.myPCAPane.visualizationPanel.myView3D.setOnlyReferencedEdges(AGCT.Only_ReferencedEdges);
        myTabbedPane.myCorrelationPane.visualizationPanel.myView3D.setOnlyReferencedEdges(AGCT.Only_ReferencedEdges);
        repaint();
    }

    public void requestModificationUseDebug(boolean v) {
        useDebug.setSelected(v);
        requestModificationUseDebug();
    }

    public void requestModificationUseDebug() {
        AGCT.setDebug(useDebug.isSelected());
        Scenario.add("AGCT_Modification_Use_Debug", useDebug.isSelected());
    }

    public void requestModificationUseWarning(boolean v) {
        useWarning.setSelected(v);
        requestModificationUseWarning();
    }

    public void requestModificationUseWarning() {
        AGCT.setWarning(useWarning.isSelected());
        Scenario.add("AGCT_Modification_Warning", useWarning.isSelected());
    }

    /**
     * Scenario以外の場合のデータロードトリガです
     *
     * @param fname
     */
    public void requestLoad_Raw_Domain(String fname) {
        Debug.debug("AGCT # requestLoad_Raw_Domain()", fname);


        JInformationFrame.getInstance().appendText("Opening file: " + fname);
        data.setRawDomainExists(true);
        data.setMyDomainFile(new File(fname));


        Domain.getInstance().init(data.getMyDomainFile());
        data.setMyDomain(Domain.getInstance());
        myTabbedPane.setDomain(data.getMyDomain());
        ControlProcess.put("dataLoaded", true);
    }

    /**
     * Scenarioの方
     */
    private void requestLoad_Raw_Domain() {
        System.out.println("AGCT.requestLoad_Raw_Domain()<Scenario>");
        if (AGCT.MYDEBUG) {
            AGCT.debug("AGCT.requestLoad_Raw_Domain()");
        }
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        final int returnVal = chooser.showOpenDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().appendText("Opening file: " + chooser.getSelectedFile().getName());
            data.setRawDomainExists(true);
            data.setMyDomainFile(chooser.getSelectedFile());

            Domain.getInstance().init(data.getMyDomainFile());
            data.setMyDomain(Domain.getInstance());// ここに注目、Domain作成開始点
            // domain.setinformationFrame(JInformationFrame.getInstance());
            myTabbedPane.setDomain(data.getMyDomain());
            Scenario.add("Load", chooser.getSelectedFile().getAbsolutePath());
            ControlProcess.put("dataLoaded", true);
        }
    }

    public void requestLoadResponsiveFile() {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        chooser.setMultiSelectionEnabled(true);
        final int ret = chooser.showOpenDialog(myFrame);
        if (chooser.getSelectedFile() != null && ret == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().appendText("Opening file: " + chooser.getSelectedFile().getName());
            // data.getMyDomain().getGenes().loadResponsiveGene(chooser.getSelectedFile());
            JLeftIndicator.getInstance().readResponsiveGeneFile(chooser.getSelectedFiles());
        }
    }

    public void requestSaveHighlightFile() {
        final char[] mask = JOptionPane.showInputDialog(this, "Bitmask of ligands", "").toCharArray();
        final BitSet bitmask = new BitSet();
        for (int i = 0; i < mask.length; i++) {
            if (mask[i] == '1') {
                bitmask.set(i);
            }
        }
        final double percentile = Double.parseDouble(JOptionPane.showInputDialog(this, "percentile", "")) / 100.;

        final JFileChooser chooser = new JFileChooser();
        chooser.setApproveButtonText("Save as:");
        final int ret = chooser.showSaveDialog(myFrame);
        if (chooser.getSelectedFile() != null && ret == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().appendText("Opening file: " + chooser.getSelectedFile().getName());
            data.getMyDomain().getGenes().saveAnnotationFile(chooser.getSelectedFile(), bitmask, percentile);
        }
    }

    public void requestLoad_HighlightFile() {
        final JFileChooser chooser = new JFileChooser();
        final int returnVal = chooser.showOpenDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().appendText("Opening file: " + chooser.getSelectedFile().getName());
            data.getMyDomain().loadHighlights(chooser.getSelectedFile());
            if (!highLightGene.isSelected()) {
                highLightGene.doClick();
            }
        }
    }
    /* scenario実行*/

    public void requestLoadExecute_Scenario() {
        final JFileChooser chooser = new JFileChooser();
        final ExampleFileFilter filter = new ExampleFileFilter();
        filter.addExtension("txt");
        filter.setDescription(".txt Files");
        chooser.setFileFilter(filter);
        chooser.setApproveButtonText("Load Scenario");
        final int returnVal = chooser.showOpenDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().setText("Opening and executing scenario: " + chooser.getSelectedFile().getName());
            Scenario.loadExecute(chooser.getSelectedFile(), this);
            ControlProcess.put("scenarioLoaded", true);
            JInformationFrame.getInstance().setText("Opening and executing scenario: " + chooser.getSelectedFile().getName() + " done...");
        }
    }

    public void requestLoad_Clustering() {
        final JFileChooser chooser = new JFileChooser();
        final ExampleFileFilter filter = new ExampleFileFilter();
        filter.addExtension("txt");
        filter.setDescription(".txt Files");
        chooser.setFileFilter(filter);
        chooser.setApproveButtonText("Load/Exec. Clustering");
        final int returnVal = chooser.showOpenDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().setText("Opening and executing clustering: " + chooser.getSelectedFile().getName());
            data.getMyAGCTFileWriter().loadClustering(chooser.getSelectedFile(), this);
            JInformationFrame.getInstance().setText("Opening and executing clustering: " + chooser.getSelectedFile().getName() + " done...");
        }
    }

    public void requestSave_Scenario() {
        if (Scenario.allStrings != null && Scenario.allStrings.size() > 0) {
            final JFileChooser chooser = new JFileChooser();
            final ExampleFileFilter filter = new ExampleFileFilter();
            filter.addExtension("txt");
            filter.setDescription(".txt Files");
            chooser.setFileFilter(filter);
            chooser.setApproveButtonText("Save Scenario");
            final int returnVal = chooser.showSaveDialog(myFrame);
            if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
                JInformationFrame.getInstance().setText("Saving scenario to file: " + chooser.getSelectedFile().getName() + " (please wait)...");
                data.getMyAGCTFileWriter().toSavingScenario(chooser.getSelectedFile());
                ControlProcess.put("scenarioSaved", true);
                JInformationFrame.getInstance().setText("Saving scenario to file: " + chooser.getSelectedFile().getName() + " done...");
            }
        }
    }

    public void requestSave_Data() {
        if (AGCT.MYDEBUG) {
            AGCT.debug("AGCT.requestSave_Data()");
        }

        final JFileChooser chooser = new JFileChooser();
        final ExampleFileFilter filter = new ExampleFileFilter();
        filter.addExtension("txt");
        filter.setDescription(".txt Files");
        chooser.setFileFilter(filter);
        chooser.setApproveButtonText("Save Data");
        final int returnVal = chooser.showSaveDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().setText("Saving processed data to file: " + chooser.getSelectedFile().getName());
            //						data.getMyAGCTFileWriter().toSaving(chooser.getSelectedFile());
            CellDesignerInput.save(chooser.getSelectedFile());
            ControlProcess.put("dataSaved", true);
            JInformationFrame.getInstance().setText("Saving processed data to file: " + chooser.getSelectedFile().getName() + " done...");
        }
    }

    // frank
    public void requestReflection() {

        System.out.println("Performing point set reflexions, axis Z");
        int i;
        Gene g;

        for (i = 0; i < data.getMyDomain().numberSelectedGenes; i++) {
            g = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[i]);

            final double z = g.manifold_point3D.coordinates[2];
            g.manifold_point3D.coordinates[2] = -z;
            System.out.println(z + "-->" + -z);

            data.getMyDomain().getGenes().set(data.getMyDomain().selectedGeneNumberToGeneNumber[i], g);

        }

    }

    // frank
    public void requestSave_Manifold() {
        final JFileChooser chooser = new JFileChooser();
        final ExampleFileFilter filter = new ExampleFileFilter();
        filter.addExtension("txt");
        filter.setDescription(".txt Files");
        chooser.setFileFilter(filter);
        chooser.setApproveButtonText("Save 3D manifold data");
        final int returnVal = chooser.showSaveDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().setText("Saving manifold data to file: " + chooser.getSelectedFile().getName());

            data.getMyAGCTFileWriter().toSavingManifold(chooser.getSelectedFile());

            ControlProcess.put("dataSaved", true);
            JInformationFrame.getInstance().setText("Saving manifold data to file: " + chooser.getSelectedFile().getName() + " done...");
        }
    }

    public void requestSave_ThisClustering(int nclustering) {
        final JFileChooser chooser = new JFileChooser();
        final ExampleFileFilter filter = new ExampleFileFilter();
        filter.addExtension("txt");
        filter.setDescription(".txt Files");
        chooser.setFileFilter(filter);
        chooser.setApproveButtonText("Save The Visible Clustering");
        final int returnVal = chooser.showSaveDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().setText("Saving clustering to file: " + chooser.getSelectedFile().getName());
            data.getMyAGCTFileWriter().toSavingClustering(chooser.getSelectedFile(), nclustering);
            ControlProcess.put("dataSaved", true);
            JInformationFrame.getInstance().setText("Saving clustering to file: " + chooser.getSelectedFile().getName() + " done...");
        }
    }

    public void requestSave_DataChi2() {
        if (Statistics.allChi2Tests) {
            final JFileChooser chooser = new JFileChooser();
            final ExampleFileFilter filter = new ExampleFileFilter();
            filter.addExtension("txt");
            filter.setDescription(".txt Files");
            chooser.setFileFilter(filter);
            chooser.setApproveButtonText("Save Data");
            final int returnVal = chooser.showSaveDialog(myFrame);
            if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
                JInformationFrame.getInstance().setText("Saving Chi2 data to file: " + chooser.getSelectedFile().getName());
                data.getMyAGCTFileWriter().toSavingChi2(chooser.getSelectedFile());
                ControlProcess.put("dataSavedChi2", true);
                JInformationFrame.getInstance().setText("Saving Chi2 data to file: " + chooser.getSelectedFile().getName() + " done...");
            }
        }
    }

    public void requestSave_ClusterProfile() {
        Debug.debug("AGCT # requestSave_ClusterProfile");
        final JFileChooser chooser = new JFileChooser();
        final ExampleFileFilter filter = new ExampleFileFilter();
        filter.addExtension("txt");
        filter.setDescription(".txt Files");
        chooser.setFileFilter(filter);
        chooser.setApproveButtonText("Save Data");
        final int returnVal = chooser.showSaveDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().setText("Saving Cluster profile data to file: " + chooser.getSelectedFile().getName());
            ClusteringProfileFile.toSavingClusterProfile(chooser.getSelectedFile());
            ControlProcess.put("dataSavedClusterProfile", true);
            JInformationFrame.getInstance().setText("Saving Cluster profile data to file: " + chooser.getSelectedFile().getName() + " done...");
        }
    }

    public void requestSave_DataDelaunay() {
        if (data.isTriangulated()) {
            final JFileChooser chooser = new JFileChooser();
            final ExampleFileFilter filter = new ExampleFileFilter();
            filter.addExtension("txt");
            filter.setDescription(".txt Files");
            chooser.setFileFilter(filter);
            chooser.setApproveButtonText("Save Data");
            final int returnVal = chooser.showSaveDialog(myFrame);
            if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
                JInformationFrame.getInstance().setText("Saving Triangulated data to file: " + chooser.getSelectedFile().getName());
                data.getMyAGCTFileWriter().toSavingDelaunay(chooser.getSelectedFile());
                ControlProcess.put("dataSavedChi2", true);
                JInformationFrame.getInstance().setText("Saving Triangulated data to file: " + chooser.getSelectedFile().getName() + " done...");
            }
        }
    }

    public void requestSave_DataWeb() {
        final JFileChooser chooser = new JFileChooser();
        final ExampleFileFilter filter = new ExampleFileFilter();
        filter.addExtension("txt");
        filter.setDescription(".txt Files");
        chooser.setFileFilter(filter);
        chooser.setApproveButtonText("Save (Web)");
        final int returnVal = chooser.showSaveDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            JInformationFrame.getInstance().setText("Saving manifold to web file: " + chooser.getSelectedFile().getName());
            data.getMyAGCTFileWriter().toSavingWeb(chooser.getSelectedFile());
            JInformationFrame.getInstance().setText("Saving manifold to web file: " + chooser.getSelectedFile().getName() + " done...");
        }
    }

    public void captureAndSave() {
        final Rectangle rect = myTabbedPane.getCaptureRectangle();
        BufferedImage framecapture = null;
        try {
            framecapture = new Robot().createScreenCapture(rect);
        } catch (final AWTException e) {
            Debug.debug("Error on AGCT # captureAndSave()");
        }

        final JFileChooser chooser = new JFileChooser();
        final ExampleFileFilter filter = new ExampleFileFilter();
        filter.addExtension("png");
        filter.setDescription(".png Files");
        chooser.setFileFilter(filter);
        chooser.setApproveButtonText("Save");
        final int returnVal = chooser.showSaveDialog(myFrame);
        if (chooser.getSelectedFile() != null && returnVal == JFileChooser.APPROVE_OPTION) {
            try {
                ImageIO.write(framecapture, "png", chooser.getSelectedFile());
            } catch (final IOException e) {
                Debug.debug("error on AGCT # capcureAndSave");
            }
            JInformationFrame.getInstance().setText("Saving processed data to file: " + chooser.getSelectedFile().getName());
            ControlProcess.put("frameCaptured", true);
        }
    }

    public void setInfoTitle(String s) {
        JInformationFrame.getInstance().setTitle(s);
    }

    public void setInfoText(String s) {
        JInformationFrame.getInstance().setText(s);
    }

    public void appendInfoText(String s) {
        JInformationFrame.getInstance().appendText(s);
    }

    void displayHTMLGene(Gene gg) {
        JInformationFrame.getInstance().setText(gg.toLightHTMLString());
    }

    public void init() {
        if (AGCT.MYDEBUG) {
            AGCT.debug("AGCT.init()");
        }
        data.setDimElements(-1);

        try {
            SwingUtilities.invokeAndWait(this);
        } catch (final Exception e) {
            System.err.println("Initialization failure");
        }
    }

    public void initMenuItems() {
        useDebug.setSelected(true);
        requestModificationUseDebug();

        useWarning.setSelected(true);
        requestModificationUseWarning();

        sortDepth.setSelected(false);
        requestModificationSort_Depth();

        useShadow.setSelected(false);
        requestModificationUse_Shadow();

        perspective.setSelected(false);
        requestModificationPerspective();

        highLightGene.setEnabled(false);
        highLightGene.setSelected(false);
        requestModificationOnly_Referenced();

        onlyReferencedEdges.setEnabled(false);
        onlyReferencedEdges.setSelected(false);
        requestModificationOnly_ReferencedEdges();

        // save20100216.setEnabled(false);

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

        requestMofificationMethod_F();
        requestModificationMethod_S();
        requestMofificationMethod_W();
        requestMofificationMethod_N();
    }

    public String getParameterString() {
        String val;

        val = " Feature type: ";
        if (Method_F == 0) {
            val += "slopes\n";
        } else if (Method_F == 1) {
            val += "Haar wavelets\n";
        } else if (Method_F == 2) {
            val += "Daubechies D4 wavelets on an interval\n";
        } else if (Method_F == 3) {
            val += "Daubechies D4 wavelets\n";
        } else if (Method_F == 4) {
            val += "Daubechies D8 wavelets\n";
        } else if (Method_F == 5) {
            val += "Daubechies D10 wavelets\n";
        } else if (Method_F == 6) {
            val += "Daubechies D12 wavelets\n";
        } else if (Method_F == 7) {
            val += "Daubechies D14 wavelets\n";
        } else if (Method_F == 8) {
            val += "Daubechies D16 wavelets\n";
        } else if (Method_F == 9) {
            val += "Daubechies D18 wavelets\n";
        } else if (Method_F == 10) {
            val += "Daubechies D20 wavelets\n";
        }

        val += " Similarity measure: ";
        if (Method_S == 0) {
            val += "Heat kernel with K = " + AGCT.T_Heat_Kernel + "\n";
        } else if (Method_S == 1) {
            val += "Cosine similarity\n";
        }

        val += " Similarity matrix W: ";
        if (Method_W == 0) {
            val += "used similarity as is\n";
        } else if (Method_W == 1) {
            val += "filtered with K=" + AGCT.Number_Of_Neighbors + " symmetric neighbors\n";
        }

        val += " Sparsification of W: " + 100 * AGCT.Sparsify_Statistic + "% of the smallest values become 0.0\n";

        val += " Normalization of W: ";
        if (Method_N == 0) {
            val += "Normalized Graph Laplacian, N = D^{-1/2} W D^{-1/2}\n";
        } else if (Method_N == 1) {
            val += "Row-Stochastic Normalization, N = D^{-1} W\n";
        } else if (Method_N == 2) {
            val += "Doubly-Stochastic approximation of W\n";
        } else if (Method_N == 3) {
            val += "Doubly-Stochastic approximation of (Doubly-Stochastic Scaling of W)\n";
        }

        if (Method_N == 3) {
            val += "Divergence chosen : " + JDSScaling.divergences[myScaling.indexSelected] + "\n";
        }

        val += " Number of Manifold components kept for each gene (saves space): " + AGCT.Number_Of_Manifold_Components + "\n";
        val += " Number of dimensions used for triangulation: " + AGCT.Number_Of_Triangulation_Dimensions + "\n";

        val += "\n Plotting options at saving time: ";

        if (AGCT.Sort_Depth) {
            val += " Genes sorted wrt depth, ";
        } else {
            val += " Depth not used, ";
        }

        if (AGCT.Use_Shadow) {
            val += " Shadow plots, ";
        } else {
            val += " Shadow not used, ";
        }

        if (AGCT.Perspective) {
            val += " Perspective used.";
        } else {
            val += " Perspective not used.";
        }

        if (AGCT.Only_Referenced) {
            val += " Only referenced genes used.";
        } else {
            val += " All genes used.";
        }

        val += "\n";

        return val;
    }

    public JMenuBar menuAGCT() {
        useDebug = new JCheckBoxMenuItem("Debug");
        useDebug.setSelected(true);
        requestModificationUseDebug();
        useDebug.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent e) {
                requestModificationUseDebug();
            }
        });

        useWarning = new JCheckBoxMenuItem("Warning");
        useWarning.setSelected(true);
        requestModificationUseWarning();
        useWarning.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent e) {
                requestModificationUseWarning();
            }
        });

        sortDepth = new JCheckBoxMenuItem("Use depth");
        sortDepth.setSelected(false);
        requestModificationSort_Depth();
        sortDepth.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent e) {
                requestModificationSort_Depth();
            }
        });

        noGeneSelection = new JCheckBoxMenuItem("No gene selection");
        noGeneSelection.setSelected(false);
        requestModificationNo_GeneSelection();
        noGeneSelection.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent e) {
                requestModificationNo_GeneSelection();
            }
        });

        useShadow = new JCheckBoxMenuItem("Use shadow");
        useShadow.setSelected(false);
        requestModificationSort_Depth();
        useShadow.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent e) {
                requestModificationUse_Shadow();
            }
        });

        numberOfProfiles = new JMenuItem("Max #Profiles Displayed");
        numberOfProfiles.setMnemonic(KeyEvent.VK_D);
        numberOfProfiles.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationNumber_Of_Profiles();
            }
        });

        // Frank
        randomSeed = new JMenuItem("Random seed");
        randomSeed.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationRandomSeed();
            }
        });

        //

        perspective = new JCheckBoxMenuItem("Perspective");
        perspective.setSelected(false);
        requestModificationPerspective();
        perspective.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent e) {
                requestModificationPerspective();
            }
        });

        highLightGene = new JCheckBoxMenuItem("Highlight: referenced genes");
        highLightGene.setSelected(false);
        requestModificationOnly_Referenced();
        highLightGene.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent e) {
                requestModificationOnly_Referenced();
            }
        });

        onlyReferencedEdges = new JCheckBoxMenuItem("Filter: referenced edges");
        onlyReferencedEdges.setSelected(false);
        requestModificationOnly_ReferencedEdges();
        onlyReferencedEdges.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent e) {
                requestModificationOnly_ReferencedEdges();
            }
        });

        JMenuItem numericalIterationSettingMenuItem = new JMenuItem("Number of Numerical Iteration (EigenVector)");
        numericalIterationSettingMenuItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                requestModificationNumber_Of_Matrx_Max_Iteration();
            }
        });


        loadResponsiveFile = new JMenuItem("Load: responsive file");
        loadResponsiveFile.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent arg0) {
                requestLoadResponsiveFile();
            }
        });

        JMenuItem saveHighlightFile = new JMenuItem("Save: highlight file(same peak. choose usee ligands & percentile)");
        saveHighlightFile.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent arg0) {
                requestSaveHighlightFile();

                // new Thread(){

                // };
                // }.start();

            }
        });

        loadHighlightFile = new JMenuItem("Load: highlight file");
        loadHighlightFile.setMnemonic(KeyEvent.VK_L);
        loadHighlightFile.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestLoad_HighlightFile();
            }
        });

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

        final ButtonGroup groupF = new ButtonGroup();
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

        // meanMode_ArithmeticMean = new
        // JRadioButtonMenuItem("Arithmetic Mean");
        // meanMode_ArithmeticMean.setActionCommand(AGCT.MS);
        // meanMode_GeometricMean = new JRadioButtonMenuItem("Geometric Mean");
        // meanMode_GeometricMean.setActionCommand(AGCT.MS);
        // groupMeanMode.add(meanMode_ArithmeticMean);
        // groupMeanMode.add(meanMode_GeometricMean);

        methodS_0 = new JRadioButtonMenuItem("Heat kernel");
        methodS_0.setActionCommand(AGCT.MS);

        methodS_1 = new JRadioButtonMenuItem("Cosine similarity");
        methodS_1.setActionCommand(AGCT.MS);

        methodS_2 = new JRadioButtonMenuItem("Absolute Cosine similarity");
        methodS_2.setActionCommand(AGCT.MS);

        final ButtonGroup groupS = new ButtonGroup();
        groupS.add(methodS_0);
        groupS.add(methodS_1);
        groupS.add(methodS_2);

        methodW_0 = new JRadioButtonMenuItem("Keep similarity as is");
        methodW_0.setActionCommand(AGCT.MW);

        methodW_1 = new JRadioButtonMenuItem("Filter via Symmetric NN");
        methodW_1.setActionCommand(AGCT.MW);

        methodW_2 = new JRadioButtonMenuItem("Filter via Natural Neighbors");
        methodW_2.setActionCommand(AGCT.MW);

        final ButtonGroup groupW = new ButtonGroup();
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

        final ButtonGroup groupN = new ButtonGroup();
        groupN.add(methodN_0);
        groupN.add(methodN_1);
        groupN.add(methodN_2);
        groupN.add(methodN_3);
        final JMenu jFiltering = new JMenu("Filtering");
//        final JMenu jIterativeFiltering = new JMenu("Iterative Filtering");
//        final JMenuItem jFilterByClusterIds = new JMenuItem("Filter by clusterId ");
        final JMenuItem jFilterByRangeOfValue = new JMenuItem("Filter by Range of Value");
        final JMenuItem jFilterByPeakTime = new JMenuItem("Filter by Peak Time");
        final JMenuItem jFilterByMagnitude = new JMenuItem("Filter by Magnitude");
        final JMenuItem jFilterByDelaunay = new JMenuItem("Filter by Delaunay");

        final JMenu jAddFilterUsingThisView = new JMenu("add Filter using this View");

        final JMenuItem jAddFilterRangeOfMagnitude = new JMenuItem("range of Magnitude");
        final JMenuItem jAddFilterRangeOfDelaunay = new JMenuItem("range of delaunay");
        final JMenuItem jAddFilterRangeOfPeakValue = new JMenuItem("range of peak value for all doses");
        jAddFilterUsingThisView.add(jAddFilterRangeOfMagnitude);
        jAddFilterUsingThisView.add(jAddFilterRangeOfDelaunay);
        jAddFilterUsingThisView.add(jAddFilterRangeOfPeakValue);

//        jFiltering.add(jIterativeFiltering);
//        jFiltering.add(jFilterByClusterIds);
        jFiltering.add(jFilterByRangeOfValue);
        jFiltering.add(jFilterByPeakTime);
        jFiltering.add(jFilterByMagnitude);
        jFiltering.add(jFilterByDelaunay);
        jFiltering.add(jAddFilterUsingThisView);

        final JMenu jComputeAndDisplayParameters = new JMenu("Comput. & Display Parameters");
        jComputeAndDisplayParameters.setMnemonic(KeyEvent.VK_T);

        jNormalizationMenu = new JNormalizationMenu();
        jNormalizationMenu.jReplicateNormalizationMenu.setNormalizationData(normalizationData);
        jNormalizationMenu.addSDThresholdChangeListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                AGCT.getInstance().myTabbedPane.mySelectionPane.computeLabel();
            }
        });

        final JMenu proc = new JMenu("Processing Parameters");
        proc.setMnemonic(KeyEvent.VK_P);

        final JMenu about = new JMenu("About");
        about.setMnemonic(KeyEvent.VK_A);

        final JMenu subF = new JMenu("Feature Types");
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

        // JMenu meanMode = new JMenu("Mean Mode(set before loading file.)");
        // meanMode.add(meanMode_ArithmeticMean);
        // meanMode.add(meanMode_GeometricMean);

        final JMenu similarityMeasureMenu = new JMenu("Similarity Measure");
        similarityMeasureMenu.setMnemonic(KeyEvent.VK_S);
        similarityMeasureMenu.add(methodS_0);
        similarityMeasureMenu.add(methodS_1);
        similarityMeasureMenu.add(methodS_2);

        final JMenu subW = new JMenu("Similarity Matrix W");
        subW.setMnemonic(KeyEvent.VK_S);
        subW.add(methodW_0);
        subW.add(methodW_1);

        final JMenu subN = new JMenu("Normalized Matrix N");
        subN.setMnemonic(KeyEvent.VK_N);
        subN.add(methodN_0);
        subN.add(methodN_1);
        subN.add(methodN_2);
        subN.add(methodN_3);

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
        // meanMode_ArithmeticMean.addActionListener(this);
        // meanMode_GeometricMean.addActionListener(this);
        methodS_0.addActionListener(this);
        methodS_1.addActionListener(this);
        methodS_2.addActionListener(this);
        methodW_0.addActionListener(this);
        methodW_1.addActionListener(this);
        methodN_0.addActionListener(this);
        methodN_1.addActionListener(this);
        methodN_2.addActionListener(this);
        methodN_3.addActionListener(this);

        numberOfNeighbors = new JMenuItem("Number Of Neighbors");
        numberOfNeighbors.setMnemonic(KeyEvent.VK_O);
        numberOfNeighbors.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationNumber_Of_Neighbors();
            }
        });

        numberOfManifoldComponents = new JMenuItem("Number Of Manifold Components");
        numberOfManifoldComponents.setMnemonic(KeyEvent.VK_M);
        numberOfManifoldComponents.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationNumber_Of_Manifold_Components();
            }
        });

        numberOfTriangulationDimensions = new JMenuItem("Dimension of the Triangulation");
        numberOfTriangulationDimensions.setMnemonic(KeyEvent.VK_T);
        numberOfTriangulationDimensions.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationNumber_Of_Triangulation_Dimensions();
            }
        });

        tHeatKernel = new JMenuItem("T Heat Kernel (W)");
        tHeatKernel.setMnemonic(KeyEvent.VK_T);
        tHeatKernel.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationT_Heat_Kernel();
            }
        });

        tHeatKernelSimilarityNeighbors = new JMenuItem("T Heat Kernel (Manifold)");
        tHeatKernelSimilarityNeighbors.setMnemonic(KeyEvent.VK_H);
        tHeatKernelSimilarityNeighbors.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationT_Heat_Kernel_Similarity_Neighbors();
            }
        });

        sparsifyStatistic = new JMenuItem("Sparsify value");
        sparsifyStatistic.setMnemonic(KeyEvent.VK_S);
        sparsifyStatistic.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationSparsify_Statistic();
            }
        });

        limitp = new JMenuItem("Limit P (Chi2)");
        limitp.setMnemonic(KeyEvent.VK_P);
        limitp.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationLimit_P_Chi2();
            }
        });

        limitpdelaunay = new JMenuItem("Limit P (Delaunay)");
        limitpdelaunay.setMnemonic(KeyEvent.VK_D);
        limitpdelaunay.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationLimit_P_Delaunay();
            }
        });

        limitpmagnitude = new JMenuItem("Limit P (Magnitiude)");
        limitpmagnitude.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationLimit_P_Magnitude();
            }
        });

        bregDiv = new JMenuItem("Bregman Divergence (DS*)");
        bregDiv.setMnemonic(KeyEvent.VK_B);
        bregDiv.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                requestModificationBregDiv();
            }
        });

        final JMenuItem aboutItem = new JMenuItem("About AGCT");
        aboutItem.setMnemonic(KeyEvent.VK_A);
        aboutItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                showAbout();
            }
        });
        about.add(aboutItem);

        final JMenuBar menuBar = new JMenuBar();
        jNormalizationMenu.addSeparator();
        menuBar.add(jNormalizationMenu);
        // menuBar.add(jNormalization);
        menuBar.add(proc);
        menuBar.add(jFiltering);
        menuBar.add(jComputeAndDisplayParameters);
//        menuBar.add(jGarudaOptionMenu);
        menuBar.add(about);

        final JMenu jSortMode = new JMenu("sort by");
        final JRadioButtonMenuItem commonGene = new JRadioButtonMenuItem("common genes");
        final JRadioButtonMenuItem KBCoverage = new JRadioButtonMenuItem("KB coverage");
        final ButtonGroup group = new ButtonGroup();
        group.add(commonGene);
        group.add(KBCoverage);
        jSortMode.add(commonGene);
        jSortMode.add(KBCoverage);
        KBCoverage.setSelected(true);
//        jIterativeFiltering.add(jSortMode);
        goFiltering = new JMenuItem("Go");
        goFiltering.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent arg0) {
                new Thread() {

                    public void run() {
                        final IterativeRun it = new IterativeRun();
                        it.sortMode = commonGene.isSelected() ? SortMode.CommonGenes : SortMode.KBCoverage;
                        it.setDomain(data.getMyDomain());// TODO FIX!
                        if (data.getMyDomain().getSelectedGenes() == null) {
                            it.initialGenes = data.getMyDomain().getGenes();
                        } else {
                            it.initialGenes = data.getMyDomain().getSelectedGenes();
                        }
                        it.iterativeRun();
                    }

                }

                        .start();
            }
        });
//        jIterativeFiltering.add(goFiltering);

//        jFilterByClusterIds.addActionListener(new ActionListener() {
//
//            @Override
//            public void actionPerformed(ActionEvent arg0) {
//                final String _clusterIds = JOptionPane.showInputDialog(null, "clusterId (punctuate with ,)");
//                ArrayList<Integer> list = new ArrayList<Integer>();
//                try {
//                    for (String s : _clusterIds.split(",")) {
//                        list.add(Integer.valueOf(s));
//                    }
//                } catch (Exception e) {
//                    JOptionPane.showMessageDialog(null, "Invalid input.");
//                    return;
//                }
//                int[] is = new int[list.size()];
//                for (int i = 0; i < is.length; i++) {
//                    is[i] = list.get(i);
//                }
//                JLeftIndicator.getInstance().setCluster(is);
//            }
//        });

        jFilterByRangeOfValue.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                final String _doseId = JOptionPane.showInputDialog(null, "doseId(0-origin)", 0 + "");
                int doseId;
                try {
                    doseId = Integer.valueOf(_doseId);
                } catch (final NumberFormatException nfe) {
                    JOptionPane.showMessageDialog(null, "Invalid input.");
                    return;
                }
                double from = 0, to = Integer.MAX_VALUE;
                final String fromto = JOptionPane.showInputDialog(null, "from-to", from + "-" + to);
                try {
                    final String[] ss = fromto.split("-");
                    from = Double.valueOf(ss[0]);
                    to = Double.valueOf(ss[1]);
                } catch (final NumberFormatException nfe) {
                    JOptionPane.showMessageDialog(null, "Invalid input.");
                    return;
                }

                JLeftIndicator.getInstance().setRange(doseId, from, to);

            }
        });

        jFilterByPeakTime.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                final String _doseId = JOptionPane.showInputDialog(null, "doseId(0-origin)", 0 + "");
                int doseId;
                try {
                    doseId = Integer.valueOf(_doseId);
                } catch (final NumberFormatException nfe) {
                    JOptionPane.showMessageDialog(null, "Invalid input.");
                    return;
                }
                final String _peakTimeId = JOptionPane.showInputDialog(null, "timeId(0-origin)", 0 + "");
                int peakTimeId = -1;
                try {
                    peakTimeId = Integer.valueOf(_peakTimeId);
                } catch (final NumberFormatException nfe) {
                    JOptionPane.showMessageDialog(null, "Invalid input.");// spell を直しただけ(tk)
                }
                JLeftIndicator.getInstance().setPeakTime(doseId, peakTimeId);
            }
        });

        jFilterByMagnitude.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                double from = 0, to = Integer.MAX_VALUE;
                final String fromto = JOptionPane.showInputDialog(null, "min-max", from + "-" + to);
                try {
                    final String[] ss = fromto.split("-");
                    from = Double.valueOf(ss[0]);
                    to = Double.valueOf(ss[1]);
                } catch (final NumberFormatException nfe) {
                    JOptionPane.showMessageDialog(null, "Invalid input.");
                    return;
                }
                JLeftIndicator.getInstance().setRangeOfMagnitude(from - 1e-3, to + 1e-3);
            }
        });

        jFilterByDelaunay.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                double from = 0, to = Integer.MAX_VALUE;
                final String fromto = JOptionPane.showInputDialog(null, "min-max", from + "-" + to);
                if (fromto == null) {
                    return;
                }
                try {
                    final String[] ss = fromto.split("-");
                    from = Double.valueOf(ss[0]);
                    to = Double.valueOf(ss[1]);
                } catch (final NumberFormatException nfe) {
                    JOptionPane.showMessageDialog(null, "Invalid input.");
                    return;
                }
                JLeftIndicator.getInstance().setRangeOfDelaunay(from - 1e-3, to + 1e-3);
            }
        });

        jAddFilterRangeOfMagnitude.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                boolean ok = false;
                double min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;
                for (Gene gene : Gene.getGeneList()) {
                    if (gene.isVisible()) {
                        ok = true;
                        double val = gene.calcTotalDistance();
                        min = Math.min(min, val);
                        max = Math.max(max, val);
                    }
                }
                if (!ok) {
                    JOptionPane.showMessageDialog(null, "No gene selected.");
                    return;
                }
                String[] msg = {"range of magnitude is", String.format("%.2f-%.2f", min, max), "add filter (from-to)"};
                String fromto = JOptionPane.showInputDialog(null, msg, String.format("%f-%f", min - 1e-3, max + 1e-3));
                if (fromto != null) {
                    double from, to;
                    try {
                        final String[] ss = fromto.split("-");
                        from = Double.valueOf(ss[0]);
                        to = Double.valueOf(ss[1]);
                    } catch (final NumberFormatException nfe) {
                        JOptionPane.showMessageDialog(null, "Invalid input.");
                        return;
                    }
                    JLeftIndicator.getInstance().setRangeOfMagnitude(from, to);
                }
            }
        });

        jAddFilterRangeOfDelaunay.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                boolean ok = false;
                double min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;
                for (Gene gene : Gene.getGeneList()) {
                    if (gene.isVisible()) {
                        ok = true;
                        double val = gene.getMaxDelaunay();
                        if (!Double.isNaN(val)) {
                            min = Math.min(min, val);
                            max = Math.max(max, val);
                        }
                    }
                }
                if (!ok) {
                    JOptionPane.showMessageDialog(null, "No gene selected.");
                    return;
                }
                if (min == Integer.MAX_VALUE) {
                    JOptionPane.showMessageDialog(null, "No valid value.");
                    return;
                }
                String[] msg = {"range of delaunay is", String.format("%.2f-%.2f", min, max), "add filter (from-to)"};
                String fromto = JOptionPane.showInputDialog(null, msg, String.format("%f-%f", min, max));
                if (fromto != null) {
                    double from, to;
                    try {
                        final String[] ss = fromto.split("-");
                        from = Double.valueOf(ss[0]);
                        to = Double.valueOf(ss[1]);
                    } catch (final NumberFormatException nfe) {
                        JOptionPane.showMessageDialog(null, "Invalid input.");
                        return;
                    }
                    JLeftIndicator.getInstance().setRangeOfDelaunay(from - 1E-3, to + 1E-3);
                }
            }
        });

        jAddFilterRangeOfPeakValue.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                boolean ok = false;
                for (Gene gene : Gene.getGeneList()) {
                    if (gene.isVisible()) {
                        ok = true;
                        break;
                    }
                }
                if (!ok) {
                    JOptionPane.showMessageDialog(null, "No gene selected.");
                    return;
                }
                int n = Gene.getGeneList().getNumberOfDoses();
                boolean[] bs = new boolean[n];
                double[] froms = new double[n];
                double[] tos = new double[n];
                for (int doseId = 0; doseId < n; doseId++) {
                    double min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;
                    for (Gene gene : Gene.getGeneList()) {
                        if (gene.isVisible()) {
                            double val = gene.peakInitialValue(doseId);
                            min = Math.min(min, val);
                            max = Math.max(max, val);
                        }
                    }
                    String[] msg = {String.format("dose id: %d", doseId), "range of peak value is", String.format("%.2f-%.2f", min, max), "add filter (from-to)"};
                    String fromto = JOptionPane.showInputDialog(null, msg, String.format("%f-%f", min, max));
                    if (fromto != null) {
                        double from, to;
                        try {
                            final String[] ss = fromto.split("-");
                            from = Double.valueOf(ss[0]);
                            to = Double.valueOf(ss[1]);
                        } catch (final NumberFormatException nfe) {
                            JOptionPane.showMessageDialog(null, "Invalid input.");
                            return;
                        }
                        froms[doseId] = from - 1e-3;
                        tos[doseId] = to + 1e-3;
                        bs[doseId] = true;
                    }
                }
                for (int i = 0; i < n; i++) {
                    if (bs[i]) {
                        JLeftIndicator.getInstance().setRange(i, froms[i], tos[i]);
                    }
                }
            }
        });

        jComputeAndDisplayParameters.add(useDebug);
        jComputeAndDisplayParameters.add(useWarning);
        jComputeAndDisplayParameters.addSeparator();
        jComputeAndDisplayParameters.add(noGeneSelection);
        jComputeAndDisplayParameters.addSeparator();
        jComputeAndDisplayParameters.add(perspective);
        jComputeAndDisplayParameters.add(sortDepth);
        jComputeAndDisplayParameters.add(useShadow);
        jComputeAndDisplayParameters.addSeparator();
        jComputeAndDisplayParameters.add(numberOfProfiles);
        jComputeAndDisplayParameters.add(randomSeed);
        jComputeAndDisplayParameters.addSeparator();

        jComputeAndDisplayParameters.add(saveHighlightFile);
        jComputeAndDisplayParameters.add(loadResponsiveFile);
        jComputeAndDisplayParameters.addSeparator();

        jComputeAndDisplayParameters.add(loadHighlightFile);
        jComputeAndDisplayParameters.add(highLightGene);
        jComputeAndDisplayParameters.add(onlyReferencedEdges);
        jComputeAndDisplayParameters.addSeparator();
        jComputeAndDisplayParameters.add(numericalIterationSettingMenuItem);

        proc.add(similarityMeasureMenu);
        proc.add(subW);
        proc.add(subN);
        proc.addSeparator();
        proc.add(numberOfNeighbors);
        proc.add(numberOfManifoldComponents);
        proc.add(numberOfTriangulationDimensions);
        proc.add(tHeatKernel);
        proc.add(tHeatKernelSimilarityNeighbors);
        proc.add(sparsifyStatistic);
        proc.add(bregDiv);
        proc.addSeparator();
        proc.add(limitp);
        proc.add(limitpdelaunay);
        proc.add(limitpmagnitude);
        initMenuItems();

        return menuBar;
    }

    public void showAbout() {
        JOptionPane.showMessageDialog(this, History.getLatestVersion() + ", " + History.getLatestDate() + "  was developed in collaboration of \n" +
                " Natalia Polouliakh, Frank Nielsen, Hiroaki Kitano (Sony CSL)\n" +
                " Richard Nock (University of Martinique)\n" +
                " and Keigo Oka (University of Tokyo).", "About AGCT", JOptionPane.INFORMATION_MESSAGE, new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/information.png"))));
    }

    public void timerTickOn() {
        data.setTimerIsTicking(true);
    }

    public void timerTickOff() {
        data.setTimerIsTicking(false);
        clockLabel.setEnabled(false);
    }

    public void triangulationTimerTickOn() {
        data.setTriangulationTimerIsTicking(true);
    }

    public void triangulationTimerTickOff() {
        data.setTriangulationTimerIsTicking(false);
        myTabbedPane.myManifoldPane.delauButton.setIcon(myTabbedPane.myManifoldPane.triangulationIdleIcon);
    }

    static ImageIcon myImageIcon(String pathname) {
        System.out.print("My image icon");
        java.net.URL imageURL;
        imageURL = AGCT.class.getResource(pathname);
        if (imageURL != null) {
            System.out.println(" done!");
            return new ImageIcon(imageURL);
        } else {
            return null;
        }
    }

    static JButton ButtonIcon(String pathname) {
        System.out.print("Button icon");
        final java.net.URL imageURL = AGCT.class.getResource(pathname);
        if (imageURL != null) {
            System.out.println(" done!");
            return new JButton(new ImageIcon(imageURL));
        } else {
            return null;
        }
    }

    public void run() {
        if (AGCT.MYDEBUG) {
            AGCT.debug("AGCT.run()");
        }
        Number_Of_Clusters = 3;
        data.setTriangulated(false);
        data.setRawDomainExists(false);
        Scenario.tickOff();

        final Container pane = getContentPane();
        pane.setLayout(new BorderLayout());

        loadButton = ButtonIcon("/Images/folder_go.png");
        loadButton.setToolTipText("load domain (raw data)");
        loadButton.setActionCommand("load_raw_domain");
        loadButton.setEnabled(true);

        saveButton = ButtonIcon("/Images/disk.png");
        // new JButton(new
        // ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource())));
        saveButton.setToolTipText("save current results");
        saveButton.setActionCommand("save_data");
        saveButton.setEnabled(true);

        // Frank

        saveManifold = ButtonIcon("/Images/manifold.png");
        // new JButton(new
        // ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("/Images/manifold.png"))));
        saveManifold.setToolTipText("save manifold");
        saveManifold.setActionCommand("save_manifold");
        saveManifold.setEnabled(true);
        //

        // Frank
        showAxis = ButtonIcon("/Images/axis.png");
        // new JButton(new
        // ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("/Images/axis.png"))));
        showAxis.setToolTipText("show reference");
        showAxis.setActionCommand("show_axis");
        showAxis.setEnabled(true);
        //

        performReflection = ButtonIcon("/Images/reflection.png");
        performReflection.setToolTipText("perform reflection");
        performReflection.setActionCommand("perform_reflection");
        performReflection.setEnabled(true);

        saveButtonStats = ButtonIcon("/Images/diskChi2.png");
        saveButtonStats.setToolTipText("save Chi2 results");
        saveButtonStats.setActionCommand("save_data_chi2");
        saveButtonStats.setEnabled(false);

        saveButtonDelaunay = ButtonIcon("/Images/diskDel.png");
        saveButtonDelaunay.setToolTipText("save the *filtered* Delaunay triangulation");
        saveButtonDelaunay.setActionCommand("save_data_delaunay");
        saveButtonDelaunay.setEnabled(false);

        saveButtonWeb = ButtonIcon("/Images/diskWeb.png");
        saveButtonWeb.setToolTipText("save the *visible* pane coordinates");
        saveButtonWeb.setActionCommand("save_data_web");
        saveButtonWeb.setEnabled(false);

        binButton = ButtonIcon("/Images/bin_closed.png");
        binButton.setToolTipText("flush results");
        binButton.setActionCommand("flush");

        saveClusteringButton = ButtonIcon("/Images/diskClustering.png");
        saveClusteringButton.setToolTipText("save current Clusters");
        saveClusteringButton.setActionCommand("save_clustering");

        cameraButton = ButtonIcon("/Images/camera.png");

        cameraButton.setToolTipText("capture the visible graphics panel");
        cameraButton.setActionCommand("capture");

        clockLabel = new JLabel(myImageIcon("/Images/time_add.png"));
        clockLabel.setEnabled(false);

        data.setMyTickingTimer(new AGCTTimer(1000, this, AGCTTimer.GeneralTimerString));
        data.getMyTickingTimer().setInitialDelay(0);
        data.getMyTickingTimer().setCoalesce(true);
        data.getMyTickingTimer().start();

        data.setMyScenarioTimer(new AGCTTimer(1000, this, AGCTTimer.ScenarioTimerString));
        data.getMyScenarioTimer().setInitialDelay(0);
        data.getMyScenarioTimer().setCoalesce(true);
        data.getMyScenarioTimer().start();

        scenarioRecordIdleIcon = myImageIcon("/Images/script_add.png");
        scenarioRecordIdleIcon.setDescription("Idle");

        scenarioRecordingIcon = myImageIcon("/Images/script_go.png");
        scenarioRecordingIcon.setDescription("Recording");

        scenarioLoadIdleIcon = myImageIcon("/Images/folder_go.png");
        scenarioLoadIdleIcon.setDescription("Idle");

        scenarioRunningIcon = myImageIcon("/Images/cog.png");
        scenarioRunningIcon.setDescription("Running");

        scenarioLoadButton = new JButton(scenarioLoadIdleIcon);
        scenarioLoadButton.setToolTipText("load and execute scenario");
        scenarioLoadButton.setActionCommand("load_scenario");

        scenarioSaveButton = ButtonIcon("/Images/disk.png");
        scenarioSaveButton.setToolTipText("save currently recorded scenario");
        scenarioSaveButton.setActionCommand("save_scenario");

        scenarioDeleteButton = ButtonIcon("/Images/script_delete.png");
        scenarioDeleteButton.setToolTipText("erase currently recorded scenario");
        scenarioDeleteButton.setActionCommand("delete_scenario");

        scenarioRecordButton = new JButton(scenarioRecordIdleIcon);
        scenarioRecordButton.setToolTipText("record/stop recording scenario");
        scenarioRecordButton.setActionCommand("record_scenario");
        scenarioDisplayButton = ButtonIcon("/Images/magnifier.png");
        scenarioDisplayButton.setToolTipText("display current scenario in information frame");
        scenarioDisplayButton.setActionCommand("display_scenario");

        data.setMyTriangulationTickingTimer(new AGCTTimer(1000, this, AGCTTimer.TriangulationTimerString));
        data.getMyTriangulationTickingTimer().setInitialDelay(0);
        data.getMyTriangulationTickingTimer().setCoalesce(true);
        data.getMyTriangulationTickingTimer().start();

        final JToolBar mainBoxButtons = new JToolBar();
        mainBoxButtons.setFloatable(false);
        mainBoxButtons.setBorder(BorderFactory.createTitledBorder(""));

        final JToolBar scenarioButtons = new JToolBar();
        scenarioButtons.setFloatable(false);
        scenarioButtons.setBorder(BorderFactory.createTitledBorder(""));

        final JToolBar miscButtons = new JToolBar();
        miscButtons.setFloatable(false);
        miscButtons.setBorder(BorderFactory.createTitledBorder(""));

        final JToolBar boxButtons = new JToolBar();

        mainBoxButtons.add(new JLabel("Data:"));
        mainBoxButtons.add(loadButton);
        mainBoxButtons.add(saveButton);
        // Frank
        mainBoxButtons.add(saveManifold);
        mainBoxButtons.add(performReflection);
        mainBoxButtons.add(showAxis);

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
        // miscButtons.add(cameraButton);
        // miscButtons.addSeparator();
        miscButtons.add(binButton);
        miscButtons.add(saveClusteringButton);
        miscButtons.add(clockLabel);

        boxButtons.add(mainBoxButtons);
        boxButtons.add(Box.createHorizontalGlue());
        boxButtons.add(scenarioButtons);
        boxButtons.add(miscButtons);

        pane.add(boxButtons, BorderLayout.NORTH);

        myTabbedPane = new JMainFrameTabbedPane(this, data.getMyDomain());
        pane.add(myTabbedPane, BorderLayout.CENTER);

        // pane.add(lowPane, BorderLayout.SOUTH);

        setJMenuBar(menuAGCT());

        final Image img = Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/chart_line.png"));
        myFrame.setIconImage(img);

        loadButton.addActionListener(this);
        saveButton.addActionListener(this);

        // Frank
        saveManifold.addActionListener(this);
        showAxis.addActionListener(this);
        performReflection.addActionListener(this);

        saveButtonStats.addActionListener(this);
        saveButtonDelaunay.addActionListener(this);
        saveButtonWeb.addActionListener(this);
        binButton.addActionListener(this);
        saveClusteringButton.addActionListener(this);
        cameraButton.addActionListener(this);

        scenarioLoadButton.addActionListener(this);
        scenarioSaveButton.addActionListener(this);
        scenarioDeleteButton.addActionListener(this);
        scenarioRecordButton.addActionListener(this);
        scenarioDisplayButton.addActionListener(this);

        // myInformationFrame = new
        // JInformationFrame("AGCT --- Information Frame", this);
        JInformationFrame.getInstance().setSize(AGCT.WindowWidth, 400);
    }


    public static void main(String[] argv) {
        if (AGCT.MYDEBUG) {
            AGCT.debug("AGCT.main(argv)");
        }
        Prototype.flushAll();

        final JBugFreeFrame f = initAGCT();

        f.setSize(AGCT.WindowWidth, 700);
        f.getContentPane().setLayout(new BorderLayout());
        f.getContentPane().add(singleton, "Center");

        f.setVisible(true);
        singleton.startGaruda(singleton);
    }

    public static JBugFreeFrame initAGCT() {
        final JBugFreeFrame f = new JBugFreeFrame("AGCT --- A Geometric Clustering Tool");
        singleton = new AGCT(f);
        singleton.init();
        return f;
    }

    public void actionPerformed(ActionEvent e) {
        final String command = e.getActionCommand();
        if (command != null) {
            Debug.debug("AGCT # actionPerformed ", command);
        }
        //		 if (AGCT.MYDEBUG)
        //		 AGCT.debug("ActionPerformed",command);
        final String eClass = e.getSource().getClass().getName();
        //		debug("eClass",eClass);
        String idleString;

        if (eClass.equals("test.AGCTTimer")) {
            //			debug("name : ",(((AGCTTimer)e.getSource()).name));
            if (((AGCTTimer) e.getSource()).name.equals(AGCTTimer.GeneralTimerString)) {
                if (data.isTimerIsTicking()) {
                    if (clockLabel.isEnabled()) {
                        clockLabel.setEnabled(false);
                    } else {
                        clockLabel.setEnabled(true);
                    }
                }
            } else if (((AGCTTimer) e.getSource()).name.equals(AGCTTimer.ScenarioTimerString)) {
                if (Scenario.isTicking()) {
                    idleString = ((ImageIcon) scenarioRecordButton.getIcon()).getDescription();
                    if (idleString.equals("Idle")) {
                        scenarioRecordButton.setIcon(scenarioRecordingIcon);
                    } else {
                        scenarioRecordButton.setIcon(scenarioRecordIdleIcon);
                    }
                }
                if (Scenario.isRunning()) {
                    idleString = ((ImageIcon) scenarioLoadButton.getIcon()).getDescription();
                    if (idleString.equals("Idle")) {
                        scenarioLoadButton.setIcon(scenarioRunningIcon);
                    } else {
                        scenarioLoadButton.setIcon(scenarioLoadIdleIcon);
                    }
                }
            } else if (((AGCTTimer) e.getSource()).name.equals(AGCTTimer.TriangulationTimerString)) {
                if (data.isTriangulationTimerIsTicking()) {
                    idleString = ((ImageIcon) myTabbedPane.myManifoldPane.delauButton.getIcon()).getDescription();
                    if (idleString.equals("Idle")) {
                        myTabbedPane.myManifoldPane.delauButton.setIcon(myTabbedPane.myManifoldPane.triangulationRunningIcon);
                    } else {
                        myTabbedPane.myManifoldPane.delauButton.setIcon(myTabbedPane.myManifoldPane.triangulationIdleIcon);
                    }
                }
            }
        } else if (command == null) {
            // TODO handle The exception.
            forDebug.Debug.debug("AGCT : command == null");
            System.exit(1);
        } else if (command.equals("load_raw_domain")) {
            goLoadDomain();
        } else if (command.equals("save_data")) {
            requestSave_Data();
        } else if (command.equals("perform_reflection")) {
            System.out.print("Frank: Perform reflection...");
            // requestReflection();
            data.setReflectionZ(!data.isReflectionZ());
            myTabbedPane.repaint();
        } else if (command.equals("save_manifold")) {
            System.out.print("Frank: Save Manifold...");
            requestSave_Manifold();
            System.out.println("done");
        } else if (command.equals("show_axis")) {
            System.out.print("Frank: Toggle show axis...");
            data.setDrawFrame(!data.isDrawFrame());
            myTabbedPane.repaint();
        } else if (command.equals("save_data_chi2")) {
            requestSave_DataChi2();
        } else if (command.equals("save_data_delaunay")) {
            requestSave_DataDelaunay();
        } else if (command.equals("save_data_web")) {
            requestSave_DataWeb();
        } else if (command.equals("flush")) {
            flushAll();
        } else if (command.equals("save_clustering")) {
            saveClusteringProfile();
        } else if (command.equals("capture")) {
            captureAndSave();
        } else if (command.equals("record_scenario")) {
            Scenario.tickSwitch();
            if (!Scenario.isTicking()) {
                scenarioRecordButton.setIcon(scenarioRecordIdleIcon);
            } else {
                scenarioRecordButton.setIcon(scenarioRecordingIcon);
            }
        } else if (command.equals("save_scenario")) {
            requestSave_Scenario();
        } else if (command.equals("load_scenario")) {
            requestLoadExecute_Scenario();
        } else if (command.equals("display_scenario")) {
            Scenario.displayInformationFrame(this);
        } else if (command.equals("delete_scenario")) {
            if (confirmEraseScenario()) {
                Scenario.flushAll();
                scenarioRecordButton.setIcon(scenarioRecordIdleIcon);
            }
        } else if (command.equals(AGCT.MF)) {
            requestMofificationMethod_F();
        } else if (command.equals(AGCT.MW)) {
            requestMofificationMethod_W();
        } else if (command.equals(AGCT.MS)) {
            requestModificationMethod_S();
        } else if (command.equals(AGCT.MN)) {
            requestMofificationMethod_N();
        }
    }

    public boolean confirmEraseScenario() {
        final int i = JOptionPane.showConfirmDialog(this, "Confirm the erasing of the scenario");

        return i == JOptionPane.YES_OPTION;
    }

    /**
     * fileをロードする。
     */
    private void goLoadDomain() {
        System.out.println("AGCT.goLoadDomain()");

        if (AGCT.MYDEBUG) {
            Debug.debug("AGCT.goLoadDomain");
            Debug.debug("stacktrace: ", Thread.currentThread().getStackTrace());
        }

        final Thread t = new Thread() {

            public void run() {

                AGCT.getInstance().timer.start(Timer.TimerKey.LoadInitialFile);

                requestLoad_Raw_Domain();

                AGCT.getInstance().timer.finish(Timer.TimerKey.LoadInitialFile);

                if (data.isRawDomainExists()) {
                    myTabbedPane.filterData();
                    data.setDimElements(data.getMyDomain().numberSelectedGenes);
                }
            }
        };
        t.start();
    }

    public void goLoadDomain(String s) {
        Debug.debug("AGCT # goLoadDomain()", s);
        final String s2 = s;
        final Thread t = new Thread() {

            public void run() {

                AGCT.getInstance().timer.start(Timer.TimerKey.LoadInitialFile);

                requestLoad_Raw_Domain(s2);//

                AGCT.getInstance().timer.start(Timer.TimerKey.LoadInitialFile);

                if (data.isRawDomainExists()) {
                    myTabbedPane.filterData();
                    data.setDimElements(data.getMyDomain().numberSelectedGenes);
                }
            }
        };
        t.start();
    }

    public void flushAll() {
        AGCT.initDefaultClassVariables();

        ControlProcess.removeAll();
        if (myAnnotationFrame != null) {
            myAnnotationFrame.dispose();
        }
        myAnnotationFrame = null;
        data.setAnnotationsExist(false);
        data.setAnnotationsExistF(false);
        data.setAnnotationsExistP(false);
        data.setAnnotationsExistC(false);

        myClusteringFrame.setVisible(false);
        myClusteringFrame.flushAll();
        myClusteringProfileFrame.setVisible(false);
        myClusteringProfileFrame.flushAll();
        myScaling.flushAll();

        // setDomain("");
        data.setMyDomain(null);
        data.setMyDomainFile(null);
        myTabbedPane.setDomain(data.getMyDomain());
        Gene.flushAll();
        Util_Feature.flushAll();
        Prototype.flushAll();

        myTabbedPane.mySelectionPane.toTrash();
        myTabbedPane.myClusteringPane.toTrash();
        myTabbedPane.myPCAPane.toTrash();
        myTabbedPane.myManifoldPane.toTrash();
        myTabbedPane.myCorrelationPane.toTrash();

        myTabbedPane.setSelectedIndex(0);

        data.setW(null);
        data.setCW(null);
        data.setD_from_W(null);
        data.setN_from_W(null);
        data.setDWD(null);
        data.setM(null);
        data.setC(null);
        data.setE(null);
        data.setV(null);
        AGCTClustering_CP.init();
        AGCT.Loading_From_Scenario_Manifold_Begin = AGCT.Loading_From_Scenario_Manifold_End = false;
        AGCT.Show_Correlation = -1;
        AGCT.Show_Only_Genes_With_Visible_Edge = false;
        AGCT.Referenced_Available = false;
        AGCT.Max_Reference = -1;
        data.setTriangulated(false);
        saveButtonDelaunay.setEnabled(false);
        saveButtonWeb.setEnabled(false);
        myTabbedPane.myManifoldPane.delauConsistencyButton.setEnabled(false);
        myTabbedPane.myManifoldPane.saveClusteringButton.setEnabled(true);
        myTabbedPane.myManifoldPane.loadClusteringButton.setEnabled(false);
        myTabbedPane.myPCAPane.saveClusteringButton.setEnabled(true);
        myTabbedPane.myPCAPane.loadClusteringButton.setEnabled(false);

        initMenuItems();
        data.setAllClusterings(null);
        data.setNbClustering(0);
        data.setCurrentIndexW(-1);
        data.setCurrentIndexM(-1);

        JInformationFrame.getInstance().setText("Flushed all data");
        JInformationFrame.getInstance().setTitle(JInformationFrame.defaultProgressString);

        data.setRawDomainExists(false);
        Statistics.flushAll(this);
        Scenario.flushAll();

        repaint();
    }

    // //////////////////////////////////////////////////////////////////////////////////
    public void computeFeatureNames() {
        data.setFeature_names(new Vector());
        int rx = 0;
        // Gene gg = (Gene)
        // data.getMyDomain().genes.get(data.getMyDomain().selectedGeneNumberToGeneNumberBefore[0]);
        final Gene gg = data.getMyDomain().getSelectedGenes().get(0);
        for (final Ligand ligand : data.getMyDomain().getLigands()) {
            if (ligand.isChecked()) {
                final String ln = ligand.getName();

                if (gg.finalCoordinates[rx].length == 1) {
                    data.getFeature_names().addElement(ln);
                } else {
                    for (int j = 0; j < gg.finalCoordinates[rx].length; j++) {
                        data.getFeature_names().addElement(ln + (j + 1));
                    }
                }
                rx++;
            }
        }
//	printAllFeatures();
    }

    public void printAllFeatures() {
        int i, j, rx = 0;
        int index = 0;
        final Gene gg = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[0]);
        for (i = 0; i < data.getMyDomain().getLigands().size(); i++) {
            if (data.getMyDomain().getLigands().get(i).isChecked()) {
                for (j = 0; j < gg.finalCoordinates[rx].length; j++) {
                    System.out.print("[" + i + ", " + j + "] Feature " + data.getFeature_names().elementAt(index) + " : " + Util_Feature.trashed[index] + " -- ");
                    index++;
                }
                System.out.println();
                rx++;
            }
        }
    }

    public String getFeatureName(int i) {
        return (String) data.getFeature_names().elementAt((Integer) Util_Feature.featureIndexToNonTrashedIndex.get(new Integer(i)));
    }

    /**
     * たいしたことやってない
     */
    public void justBeforeManifold() {
        myTabbedPane.myManifoldPane.initDepthOrderSelectedGenes();
        JInformationFrame.getInstance().setText("...");
        data.setDimElements(data.getMyDomain().numberSelectedGenes);
    }

    public void manifoldEigensystem() {
        {
            int num = data.getMyDomain().selectedGeneNumberToGeneNumber.length;
            Gene[] genes = new Gene[num];
            for (int i = 0; i < num; i++) {
                genes[i] = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
            }
            int dim = Util_Feature.getNumberOfFeaturesAfterSelection();
            Debug.debug("num of genes and features", num, dim);
            double[][] features = new double[num][dim];
            for (int i = 0; i < num; i++) {
                for (int j = 0; j < dim; j++) {
                    features[i][j] = genes[i].getFinalCoordinates(j, true);
                }
                String name = genes[i].getAnnotation(Gene.ID);
                Debug.fdebug("tmp/genePoints.txt", name + "\t" + Arrays.toString(features[i]).replaceAll(",", "").replaceAll(" ", "\t"));
            }
        }

        justBeforeManifold();
        data.setManifold_eigenvalues(new double[data.getDimElements()]);
        ControlProcess.put("justBeforeManifold", true);
        data.setM(new Matrix("Manifold Eigenvectors", data.getDimElements(), data.getDimElements()));

        int x, y;
        double norm;

        timerTickOn();
        generate_W();

        Debug.fdebug("tmp/Manifold/W_init.txt", data.getW());

        if (AGCT.Sparsify_Statistic > 0.0) {
            data.getW().sparsifyMatrix();
        }
        Debug.fdebug("tmp/Manifold/W.txt", data.getW());

        JInformationFrame.getInstance().setText("Matrix W (10 x 10 upperleft block):\n" + data.getW().affiche(10));

        generate_D();
        debug("density of D : " + data.getD_from_W().density());

        JInformationFrame.getInstance().setText("Matrix D (10 x 10 upperleft block):\n" + data.getD_from_W().affiche(10));

        generate_N();
        debug("density of N : " + data.getN_from_W().density());

        JInformationFrame.getInstance().setText("Matrix N (10 x 10 upperleft block):\n" + data.getN_from_W().affiche(10));

        if (Method_N == 0 || Method_N == 2 || Method_N == 3) {
            data.getM().copyOnThis(data.getN_from_W());
        } else if (Method_N == 1) {
            data.getM().copyOnThis(data.getDWD());// Method_N=0で作ったのと同じ行列
        }

        Debug.fdebug("tmp/Manifold/N.txt", data.getN_from_W());
        Debug.fdebug("tmp/Manifold/L.txt", data.getM());

        MatrixMath.getEigenSystem(data.getM(), data.getManifold_eigenvalues());

        Debug.fdebug("tmp/Manifold/M_init.txt", data.getM());

        forDebug.Debug.debug("M.density()", data.getM().density());
        //		data.getM().writeTo("M.txt");

        if (Method_N == 1) {
            final int n = data.getDimElements();
            final int to = Math.min(n, MatrixConfig.MAX_COL_OF_EIGENVECTOR);
            for (x = 0; x < to; x++) {
                norm = 0.0;
                for (y = 0; y < n; y++) {
                    data.getM().set(x, y, data.getM().get(x, y) / Math.sqrt(data.getD_from_W().get(y, y)));
                    norm += data.getM().get(x, y) * data.getM().get(x, y);
                }
                norm = Math.sqrt(norm);
                for (y = 0; y < n; y++) {
                    data.getM().set(x, y, data.getM().get(x, y) / norm);
                }
            }
        }

        // Scenario : save M from here.
        Debug.fdebug("tmp/Manifold/M.txt", data.getM());
        Debug.fdebug("tmp/Manifold/eigenValues.txt", data.getManifold_eigenvalues());

        justAfterManifold();

        if (myTabbedPane.myManifoldPane.keepManifold.isSelected()) {
            Scenario.add("JAGCTManifoldPane_Keep_Manifold");
        }
        timerTickOff();
    }

    public void justAfterManifold() {
        System.out.println("AGCT.justAfterManifold()");

        data.getMyDomain().fillManifoldComponents(data.getM());

        myTabbedPane.myManifoldPane.setGenesButtons(true);
        myTabbedPane.myManifoldPane.setHiButtons(true);

        // propagating events to Swing units

        myTabbedPane.myManifoldPane.visualizationPanel.plotAvailable = true;
        myTabbedPane.myManifoldPane.initManifold3D();
        myTabbedPane.myManifoldPane.visualizationPanel.initAll(true);
        myTabbedPane.repaint();
        myTabbedPane.myManifoldPane.xText.setEnabled(true);
        myTabbedPane.myManifoldPane.yText.setEnabled(true);
        myTabbedPane.myManifoldPane.zText.setEnabled(true);
        myTabbedPane.myManifoldPane.delauButton.setEnabled(true);
        myTabbedPane.myManifoldPane.searchText.setEnabled(true);
        myTabbedPane.myManifoldPane.viewSelect.setEnabled(true);

        myTabbedPane.myManifoldPane.keepManifold.setEnabled(false);

        ControlProcess.put("manifoldProcessed", true);
        checkClusteringAvailable();
        if (myTabbedPane.getSelectedIndex() == 1) {
            saveButtonWeb.setEnabled(true);
        }
        System.out.println("END_AGCT.justAfterManifold()");
    }

    synchronized public void checkClusteringAvailable() {
        System.out.println("AGCT.checkClusteringAvailable()");
        forDebug.Debug.debug("manifoldProcessed", ControlProcess.hasTrue("manifoldProcessed"), "pcaProcessed", ControlProcess.hasTrue("pcaProcessed"));
        if (ControlProcess.hasTrue("manifoldProcessed") && ControlProcess.hasTrue("pcaProcessed")) {
            // Clustering becomes available !
            myTabbedPane.toClustering();
            myTabbedPane.myClusteringPane.activate();
            ControlProcess.put("clusteringAvailable", true);
        }
    }

    public void pcaEigensystem() {
        JInformationFrame.getInstance().setText("...");
        //		final Thread t=new Thread(){
        //			public void run(){
        timerTickOn();

        generate_C(true);

        JInformationFrame.getInstance().setText("Matrix C (10 x 10 upperleft block):\n" + data.getC().affiche(10));

        generate_E();

        JInformationFrame.getInstance().setText("Matrix E (10 x 10 upperleft block):\n" + data.getE().affiche(10));

        data.getMyDomain().fillPCAComponents(data.getE(), data.getAverage_features(), data.getSigma_features());

        generate_V();

        JInformationFrame.getInstance().setText("Matrix V (10 x 10 upperleft block):\n" + data.getV().affiche(10));

        myTabbedPane.myPCAPane.setGenesButtons(true);
        myTabbedPane.myPCAPane.setHiButtons(true);

        myTabbedPane.myPCAPane.visualizationPanel.plotAvailable = true;
        myTabbedPane.myPCAPane.initPCA3D();
        myTabbedPane.myPCAPane.visualizationPanel.initAll(true);
        myTabbedPane.repaint();
        myTabbedPane.myPCAPane.xText.setEnabled(true);
        myTabbedPane.myPCAPane.yText.setEnabled(true);
        myTabbedPane.myPCAPane.zText.setEnabled(true);
        myTabbedPane.myPCAPane.searchText.setEnabled(true);
        myTabbedPane.myPCAPane.viewSelect.setEnabled(true);

        fillCorrelationComponents();
        myTabbedPane.myCorrelationPane.visualizationPanel.plotAvailable = true;
        myTabbedPane.myCorrelationPane.initCorrelation3D();
        myTabbedPane.myCorrelationPane.visualizationPanel.initAll(true);
        myTabbedPane.repaint();
        while (myTabbedPane.myCorrelationPane == null || myTabbedPane.myCorrelationPane.xText == null) {
        }
        myTabbedPane.myCorrelationPane.xText.setEnabled(true);
        myTabbedPane.myCorrelationPane.yText.setEnabled(true);
        myTabbedPane.myCorrelationPane.zText.setEnabled(true);
        myTabbedPane.myCorrelationPane.searchText.setEnabled(true);
        myTabbedPane.myCorrelationPane.viewSelect.setEnabled(true);

        ControlProcess.put("pcaProcessed", true);
        checkClusteringAvailable();
        myTabbedPane.myManifoldPane.loadClusteringButton.setEnabled(true);
        myTabbedPane.myPCAPane.loadClusteringButton.setEnabled(true);
        if (myTabbedPane.getSelectedIndex() == 2) {
            saveButtonWeb.setEnabled(true);
        }

        timerTickOff();

        //			}
        //		};
        //		t.start();
    }

    public void generate_E() {
        data.setE(new Matrix("E", data.getDimFeatures(), data.getDimFeatures()));
        data.getE().toE(data.getMyDomain(), data.getC());
    }

    /**
     * @param after
     */
    public void generate_C(boolean after) {
        data.setC(new Matrix("C", data.getDimFeatures(), data.getDimFeatures()));
        data.getC().toC(data.getMyDomain(), after);
    }

    public void generate_V() {
        data.setV(new Matrix("V", data.getDimFeatures(), data.getDimFeatures()));
        data.getV().toV(data.getMyDomain());
    }

    public void computeMinMax_CorrelationPnt3D() {
        // Makes the assumption that V is computed

        int i, j;
        Pnt3D pt;
        for (i = 0; i < data.getDimFeatures(); i++) {
            pt = (Pnt3D) data.getFeature_pca_Pnt3D().elementAt(i);
            for (j = 0; j < 3; j++) {
                if (i == 0 || pt.coordinates[j] < data.getMin_CorrelationPnt3D().coordinates[j]) {
                    data.getMin_CorrelationPnt3D().coordinates[j] = pt.coordinates[j];
                }
                if (i == 0 || pt.coordinates[j] > data.getMax_CorrelationPnt3D().coordinates[j]) {
                    data.getMax_CorrelationPnt3D().coordinates[j] = pt.coordinates[j];
                }
            }
        }
    }

    public void fillCorrelationComponents() {
        int i;
        data.setFeature_pca_Pnt3D(new Vector());
        for (i = 0; i < data.getDimFeatures(); i++) {
            data.getFeature_pca_Pnt3D().addElement(new Pnt3D(data.getV().get(i, 0), data.getV().get(i, 1), data.getV().get(i, 2)));
        }
    }

    public void updateCorrelationComponents(int xAxis, int yAxis, int zAxis) {
        int i;
        Pnt3D p;
        if (data.getFeature_pca_Pnt3D() == null) {
            Matrix.perror("No features Pnt 3D !");
        }
        for (i = 0; i < data.getDimFeatures(); i++) {
            p = (Pnt3D) data.getFeature_pca_Pnt3D().elementAt(i);
            p.coordinates[0] = data.getV().get(i, xAxis);
            p.coordinates[1] = data.getV().get(i, yAxis);
            p.coordinates[2] = data.getV().get(i, zAxis);
        }
    }

    public void initCorrelation3D(int xAxis, int yAxis, int zAxis) {
        // int i;
        Pnt3D pt;
        for (int i = 0; i < data.getDimFeatures(); i++) {
            pt = (Pnt3D) data.getFeature_pca_Pnt3D().elementAt(i);
            pt.coordinates[0] = data.getV().get(i, xAxis);
            pt.coordinates[1] = data.getV().get(i, yAxis);
            pt.coordinates[2] = data.getV().get(i, zAxis);
        }
    }

    public void generate_W() {
        data.setW(new Matrix("W", data.getDimElements(), data.getDimElements()));
        data.getW().toW(this);
        data.setCW(new Matrix("W_Original", data.getDimElements(), data.getDimElements()));
        data.getCopyOfSimilarityMatrix().copyOnThis(data.getW());
        generate_Neighbors();// ?
        ControlProcess.put("neighborsComputed", true);
    }

    public void generate_D() {
        data.setD_from_W(new Matrix("D", data.getDimElements(), data.getDimElements()));
        data.getD_from_W().toD(data.getMyDomain(), data.getW());
    }

    public void generate_N() {
        data.setN_from_W(new Matrix("N", data.getDimElements(), data.getDimElements()));
        data.getN_from_W().toN(this, data.getD_from_W(), data.getW());
    }

    public void generate_Neighbors() {
        if (AGCT.DoDebug) {
            System.out.print("Computing All Nearest Neighbors... ");
        }
        JInformationFrame.getInstance().setTextProgressBar("ANN");
        int i, j, index, dum, percent, max, nmade;
        final double[] sim = new double[data.getDimElements() - 1];
        double ddum;

        max = data.getDimElements();
        nmade = 0;

        //		darta.setNearest_neighbors(new int[data.getDimElements()][]);
        for (i = 0; i < data.getDimElements(); i++) {
            final int[] nnei = new int[data.getDimElements() - 1];
            //						data.getNearest_neighbors()[i] = new int[data.getDimElements() - 1];
            index = 0;
            for (j = 0; j < data.getDimElements(); j++) {
                if (j != i) {
                    nnei[index] = j;
                    //					data.getNearest_neighbors()[i][index] = j;
                    sim[index] = data.getCopyOfSimilarityMatrix().get(i, j);
                    index++;
                }
            }
            //			QuickSort.quicksort(sim, data.getNearest_neighbors()[i]);
            QuickSort.quicksort(sim, nnei);
            for (j = 0; j < (data.getDimElements() - 2) / 2; j++) {
                //				dum = data.getNearest_neighbors()[i][j];
                dum = nnei[j];
                //				data.getNearest_neighbors()[i][j] = data.getNearest_neighbors()[i][data
                //						.getDimElements()
                //						- 2 - j];
                nnei[j] = nnei[data.getDimElements() - 2 - j];
                //				data.getNearest_neighbors()[i][data.getDimElements() - 2 - j] = dum;
                nnei[data.getDimElements() - 2 - j] = dum;

                ddum = sim[j];
                sim[j] = sim[data.getDimElements() - 2 - j];
                sim[data.getDimElements() - 2 - j] = ddum;
            }

            percent = (int) (100.0 * (double) nmade / (double) max);
            JInformationFrame.getInstance().setValueProgressBar(percent);
            nmade++;
        }
        //		for(int[] is:data.getNearest_neighbors()){
        //			debug(is);
        //		}

        if (AGCT.DoDebug) {
            System.out.println("ok.");
        }
        JInformationFrame.getInstance().setValueProgressBar(0);
        JInformationFrame.getInstance().setTextProgressBar(JInformationFrame.defaultProgressString);
    }

    public void triangulation() {
        AGCT.getInstance().getTimer().start(TimerKey.Triangulation);

        JInformationFrame.getInstance().setTextProgressBar("Delaunay Triang.");
        JInformationFrame.getInstance().setText("Batch processing of Delaunay Triangulation... ");
        final int d = AGCT.Number_Of_Triangulation_Dimensions;
        int nmade = 0, percent;

        timerTickOn();
        triangulationTimerTickOn();

        final double[] min = new double[d], max = new double[d];//その軸に置ける最小値、最大値
        final double[] delta = new double[d];//max-min
        final double[] extreme = new double[d];
        double[] comp;
        double dum;

        for (int i = 0; i < data.getMyDomain().numberSelectedGenes; i++) {
            final Gene g = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[i]);//選ばれたgeneを列挙
            comp = g.manifold_components_triangulation.coordinates;//fillManifoldComponentsで埋められる。manifoldの座標の先頭
            // triangulationに使われる次元(d) をとってきたもの。

            for (int j = 0; j < d; j++) {
                if (i == 0 || comp[j] < min[j]) {
                    min[j] = comp[j];
                }

                if (i == 0 || comp[j] > max[j]) {
                    max[j] = comp[j];
                }
            }
        }

        for (int i = 0; i < d; i++) {
            delta[i] = max[i] - min[i];
        }

        Pnt[] big = new Pnt[d + 1];
        comp = new double[d];
        for (int i = 0; i < d; i++) {
            comp[i] = min[i] - delta[i];
        }
        big[0] = new Pnt(comp, -1);

        for (int i = 0; i < d; i++) {
            extreme[i] = max[i] + (double) d * 2 * delta[i];

            comp = new double[d];
            for (int j = 0; j < d; j++) {
                if (j == i) {
                    comp[j] = max[i] + (double) d * 2 * delta[i];
                } else {
                    comp[j] = 0.0;
                }
            }
            big[i + 1] = new Pnt(comp, -(i + 1));
        }
        /*
         * delta = max-min
         * big[0] = min-delta
         * big[1+i] = i番目が、max[i] + d*2+dleta[i]
         *
         * extreme = max + 2*d * delta
         */

        for (int i = 0; i < data.getMyDomain().numberSelectedGenes; i++) {
            final Gene g = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
            comp = g.manifold_components_triangulation.coordinates;
            dum = 0.0;
            for (int j = 0; j < d; j++) {
                dum += comp[j] / extreme[j];
            }
            if (dum >= 1.0) {
                Matrix.perror("Gene " + i + " is outside big simplex");
            }
        }

        Simplex ori = new Simplex(big);/*すべての点を含むような大きいsimplex*/
        DelaunayTriangulation dt = new DelaunayTriangulation(ori);

        for (int i = 0; i < data.getMyDomain().numberSelectedGenes; i++) {
            percent = (int) (100.0 * (double) nmade / (double) data.getMyDomain().numberSelectedGenes);
            JInformationFrame.getInstance().setValueProgressBar(percent);
            nmade++;

            final Gene g = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
            comp = g.manifold_components_triangulation.coordinates;
            dt.delaunayPlace(new Pnt(comp, i));
        }
        data.getMyDomain().myDelaunayTriangulation = dt;

        JInformationFrame.getInstance().setText("Batch processing of Delaunay Triangulation... ok.");
        JInformationFrame.getInstance().setValueProgressBar(0);
        JInformationFrame.getInstance().setTextProgressBar(JInformationFrame.defaultProgressString);

        myTabbedPane.myManifoldPane.setManifoldButtons(true);
        myTabbedPane.myPCAPane.setManifoldButtons(true);

        data.getMyDomain().toEdges();

        triangulationTimerTickOff();
        timerTickOff();

        data.setTriangulated(true);
        saveButtonDelaunay.setEnabled(true);
        myTabbedPane.myManifoldPane.delauConsistencyButton.setEnabled(true);

        ControlProcess.put("manifoldTriangulated", true);

        AGCT.getInstance().getTimer().finish(TimerKey.Triangulation);

        if (AGCT.AUTO) {
            AGCT.getInstance().myTabbedPane.myPCAPane.pcaEigenButton.doClick();
        }
    }

    public void triangleConsistency() {
        final Thread t = new Thread() {

            public void run() {
                int allP = 0, oneN = 0, twoN = 0, allN = 0, tot, corel;
                int i, j, k, xi1, xi2, xj1, xj2, xk1, xk2, ci, cj, ck, cij, cia, cja, cs, mi, ma, mi2, ma2, tmp, nPos = 0, nNeg = 0;
                Vector<Integer> eli;
                Vector<Integer> elj;
                Vector<Integer> elk;
                Gene gg, ggn;
                double maxl = -1.0, curl, maxltot = -1.0, curltot, sizetot = 0.0, avg, avgtot, sloc = 0.0, totloc = 0.0;
                double tottot = 0.0;
                boolean okPlot;
                AGCTCounter cc;

                final Vector<Vector<Integer>> listOfEdges = new Vector<Vector<Integer>>();
                Vector<Integer> element;
                Vector<Integer> element2;
                String saff = "";

                cc = new AGCTCounter("[Lengths] Step 1 (largest visible)", data.getMyDomain().numberSelectedGenes);
                for (i = 0; i < data.getMyDomain().numberSelectedGenes; i++) {
                    gg = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
                    if (gg.neighbors != null && gg.neighbors.size() > 0) {
                        for (j = 0; j < gg.neighbors.size(); j++) {
                            okPlot = false;
                            if (Statistics.isPositiveDelauney(gg, j)) {
                                okPlot = true;
                                nPos++;
                            } else if (Statistics.isNegativeDelauney(gg, j)) {
                                okPlot = true;
                                nNeg++;
                            }
                            if (okPlot) {
                                ggn = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[((Integer) gg.neighbors.elementAt(j))]);
                                curl = Distortion.distortion_l22(gg.manifold_components_triangulation, ggn.manifold_components_triangulation);
                                sloc += 1.0;
                                totloc += curl;

                                if (curl > maxl) {
                                    maxl = curl;
                                }
                            }
                        }
                    }
                    cc.increment();
                }
                cc.end();

                if (sloc > 0.0) {
                    avg = totloc / sloc;

                    cc = new AGCTCounter("[Lengths] Step 2 (quantile)", data.getMyDomain().numberSelectedGenes);
                    for (i = 0; i < data.getMyDomain().numberSelectedGenes; i++) {
                        gg = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
                        if (gg.neighbors != null && gg.neighbors.size() > 0) {
                            for (j = 0; j < gg.neighbors.size(); j++) {
                                sizetot += 1.0;
                                ggn = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[((Integer) gg.neighbors.elementAt(j))]);
                                curltot = Distortion.distortion_l22(gg.manifold_components_triangulation, ggn.manifold_components_triangulation);
                                tottot += curltot;
                                if (curltot > maxltot) {
                                    maxltot = curltot;
                                }
                            }
                        }
                        cc.increment();
                    }
                    cc.end();
                    avgtot = tottot / sizetot;

                    saff += "Total observable edges / # >0 Corr. / # <0 Corr. = " + (nPos + nNeg) + " / " + nPos + " ( = " + 100.0 * (double) nPos / (nPos + nNeg) + "%) / " + nNeg + " ( = " + 100.0
                            * (double) nNeg / (nPos + nNeg) + "%)\n";
                    saff += "Average observable length / Average length = " + avg + " / " + avgtot + " = " + DF.format((avg / avgtot)) + "\n";
                    saff += "Max obervable length / Max length = " + maxl + " / " + maxltot + " = " + DF.format((maxl / maxltot)) + "\n";

                    cc = new AGCTCounter("[Triangles] Step 1 (edges)", data.getMyDomain().numberSelectedGenes);
                    for (i = 0; i < data.getMyDomain().numberSelectedGenes; i++) {
                        gg = data.getMyDomain().getGenes().get(data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
                        if (gg.neighbors != null && gg.neighbors.size() > 0) {
                            for (j = 0; j < gg.neighbors.size(); j++) {
                                corel = 0;
                                okPlot = false;
                                if (Statistics.isPositiveDelauney(gg, j)) {
                                    okPlot = true;
                                    corel = 1;
                                } else if (Statistics.isNegativeDelauney(gg, j)) {
                                    okPlot = true;
                                    corel = -1;
                                }

                                if (okPlot) {
                                    mi = data.getMyDomain().selectedGeneNumberToGeneNumber[i];
                                    ma = data.getMyDomain().selectedGeneNumberToGeneNumber[((Integer) gg.neighbors.elementAt(j))];
                                    if (mi > ma) {
                                        tmp = mi;
                                        mi = ma;
                                        ma = tmp;
                                    }

                                    element = new Vector<Integer>();
                                    element.addElement(mi);
                                    element.addElement(ma);
                                    element.addElement(corel);
                                    listOfEdges.addElement(element);
                                }
                            }
                        }
                        cc.increment();
                    }
                    cc.end();

                    QuickSort.quicksortVectorInteger(listOfEdges, 0);
                    i = 0;
                    do {
                        element = listOfEdges.elementAt(i);
                        mi = element.elementAt(0);
                        ma = element.elementAt(1);
                        j = i + 1;
                        do {
                            element2 = listOfEdges.elementAt(j);
                            mi2 = element2.elementAt(0);
                            ma2 = element2.elementAt(1);
                            if (mi == mi2 && ma == ma2) {
                                listOfEdges.removeElementAt(j);
                            } else {
                                j++;
                            }
                        } while (j < listOfEdges.size());
                        i++;
                    } while (i < listOfEdges.size() - 1);

                    cc = new AGCTCounter("[Triangles] Step 2 (checking #" + listOfEdges.size() + " edges)", listOfEdges.size() - 2);
                    for (i = 0; i < listOfEdges.size() - 2; i++) {
                        eli = listOfEdges.elementAt(i);
                        xi1 = eli.elementAt(0);
                        xi2 = eli.elementAt(1);
                        ci = eli.elementAt(2);
                        for (j = i + 1; j < listOfEdges.size() - 1; j++) {
                            elj = listOfEdges.elementAt(j);
                            xj1 = elj.elementAt(0);
                            xj2 = elj.elementAt(1);
                            cj = elj.elementAt(2);
                            if (xi1 == xj1 || xi1 == xj2 || xi2 == xj1 || xi2 == xj2) {
                                if (xi1 == xj1 || xi1 == xj2) {
                                    cij = xi1;
                                    cia = xi2;
                                } else {
                                    cij = xi2;
                                    cia = xi1;
                                }
                                if (cij == xj1) {
                                    cja = xj2;
                                } else {
                                    cja = xj1;
                                }

                                for (k = j + 1; k < listOfEdges.size(); k++) {
                                    elk = listOfEdges.elementAt(k);
                                    xk1 = elk.elementAt(0);
                                    xk2 = elk.elementAt(1);
                                    ck = elk.elementAt(2);

                                    cs = ci + cj + ck;
                                    if (xk1 == cia && xk2 == cja || xk2 == cia && xk1 == cja) {
                                        // System.out.println("Found triangle between cij = "
                                        // + cij +", xk1 = " + xk1 + ", xk2 = "
                                        // + xk2 + " (Sigma = " + cs + ")");

                                        if (cs == 3) {
                                            allP++;
                                        } else if (cs == 1) {
                                            oneN++;
                                        } else if (cs == -1) {
                                            twoN++;
                                        } else if (cs == -3) {
                                            allN++;
                                        }
                                    }
                                }
                            }
                        }
                        cc.increment();
                    }
                    cc.end();
                    tot = allP + oneN + twoN + allN;
                    saff += "Triangle consistency: found " + tot + " triangles, (all+) pattern = " + allP + ", (two-) pattern = " + twoN + " (one- pattern = " + oneN + ", all- pattern = "
                            + allN + ")\n";

                    JInformationFrame.getInstance().setText("Checking Delaunay properties for P-Delaunay = " + Statistics.getLIMIT_P_DELAUNAY() + ":\n" + saff);
                } else {
                    JInformationFrame.getInstance().setText("No visible edge");
                }
            }
        };
        t.start();
    }

    public void nonThreadClustering(final String commandLine) {

        AGCT.getInstance().getTimer().start(TimerKey.Clustering);

        String line = commandLine;
        int rep = 1;
        if (Clustering.containsRepeat(line)) {
            rep = Clustering.getRepeatTimes(line);
            line = Clustering.flushRepeatString(line);
        }

        for (int ii = 0; ii < rep; ii++) {
            if (Clustering.containsLoop(line)) {
                final int[] bds = Clustering.getBoundsNclusters(line);
                String scur;
                int i;
                for (i = bds[0]; i <= bds[1]; i++) {
                    scur = Clustering.getReplacementString(i, line);
                    synchronized (this) {
                        performClustering(scur);
                    }
                }
            } else {
                performClustering(line);
            }
        }

        AGCT.getInstance().getTimer().finish(TimerKey.Clustering);
        if (AGCT.AUTO) {
            AGCT.getInstance().getTimer().start(TimerKey.ALL_CHI_2);
            AGCT.getInstance().myTabbedPane.myManifoldPane.doit();
            AGCT.getInstance().getTimer().finish(TimerKey.ALL_CHI_2);

            GeneList genes = Domain.getInstance().getGenes();
            Gene.showThis(genes.toArray(new Gene[genes.size()]));
        }
        //				AGCT.getInstance().myTabbedPane.setSelectedIndex(1);
        if (Scenario.isRunning) {
            data.getMyAGCTFileWriter().toSavingClustering(new File(".clustering_data" + SUFFIX), 0);
        }

        //				forDebug.DoDebug.debug(JLeftIndicator.getInstance().getCheckBoxes());
        //				ArrayList<JCheckBox> boxes = JLeftIndicator.getInstance().getCheckBoxes();

        if (JLeftIndicator.getInstance().getCheckBoxes().size() > 0) {
            saveClusteringProfile();
            saveButton.doClick();
        }

        //		}

        AGCT.getInstance().getTimer().showResult();
    }


    public void saveClusteringProfile() {
        int experimentId = getCurrentExperimentId();
        forDebug.Debug.debug("AGCT : saveClusteringProfile()");
        int nc = getNumberOfClusters();
        forDebug.Debug.debug("number of clusters = " + nc);
        GeneList[] lists = new GeneList[nc];
        for (int i = 0; i < lists.length; i++) {
            lists[i] = new GeneList();
        }
        for (Gene gene : Domain.getInstance().getSelectedGenes()) {
            if (gene.isVisible()) {
                lists[gene.getClusterNumber(experimentId)].add(gene);
            }
        }

        PrintWriter pw;
        try {
            pw = new PrintWriter(new File(Domain.getInstance().getFileName() + "_clusterProfile.txt"));
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        int id = 0;
        List<File> fileList = new ArrayList<File>();
        for (GeneList list : lists) {
            if (list.size() == 0) {
                continue;
            }
            id++;
            final Gene grav = Gene.center(list);
            Collections.sort(list, new Comparator<Gene>() {

                @Override
                public int compare(Gene o1, Gene o2) {
                    return (int) Math.signum(Gene.getSimilarity(o1, grav) - Gene.getSimilarity(o2, grav));
                }
            });
            Gene centroid = list.get(0);
            pw.print("Cluster" + id + " : " + list.size() + " genes");
            for (int j = 0; j < Domain.getInstance().getTimes().length; j++) {
                pw.print("\ttime" + (Domain.getInstance().getTimes())[j]);
            }
            pw.println();

            for (int j = 0; j < centroid.getNumberOfLigand(); j++) {
                pw.print(getFeatureName(j));
                int fid = (Integer) Util_Feature.featureIndexToNonTrashedIndex.get(new Integer(j));
                for (int k = 0; k < centroid.getInitialCoordinate(fid).length; k++) {
                    double val = centroid.getInitialCoordinate(fid)[k];
                    pw.print("\t" + (Double.isNaN(val) ? "Not_defined" : ("" + val)));
                }
                pw.println();
            }
            pw.println();
            GeneListSaver saver = new GeneListSaver(list);
            File file = new File(Domain.getInstance().getFileName() + "_cluster_" + id + ".txt");
            fileList.add(file);
            saver.setFile(file);
            saver.execute();
        }
        pw.flush();
        pw.close();

        AGCTgaruda garudaBackend = AGCT.getInstance().getGarudaBackend();
        garudaBackend.setCandidateSendingFiles(fileList.toArray(new File[fileList.size()]));
        File file = apcl.agctbackend.getFileSentToGaruda(apcl.agctbackend.getCandidateSendingFiles(), this);
        if (file != null) {
            AGCT.getInstance().sendGaruda(file);
        }
    }

    private AGCTgaruda getGarudaBackend() {
        return apcl.agctbackend;
    }


    public int getCurrentExperimentId() {
        if (JAGCTVisualizationPane.instances.isEmpty()) return -1;
        return JAGCTVisualizationPane.instances.get(0).currentClustering;
    }

    private Clustering getCurrentShownClustering() {
        int experimentId = getCurrentExperimentId();
        if (experimentId == -1) return null;
        return (Clustering) data.getAllClusterings().get(experimentId);
    }

    private int getNumberOfClusters() {
        Clustering clustering = getCurrentShownClustering();
        if (clustering == null) return 0;
        return clustering.getNumberOfClusters();
    }

    public void clustering(final String lineu) {

        final Thread t = new Thread() {

            public void run() {
                nonThreadClustering(lineu);
            }
        };
        t.start();
    }

    public void performClustering(String line) {
        forDebug.Debug.debug("performClustering", line);
        data.setNbClustering(data.getNbClustering() + 1);
        final Clustering newClustering;
        newClustering = new Clustering(this, data.getNbClustering() - 1);
        newClustering.getOptions(line, false);
        if (newClustering.myClusteringAlgorithm == null) {
            return;
        }
        JInformationFrame.getInstance().setTextProgressBar("Clustering");
        JInformationFrame.getInstance().setText("Processing clustering : " + newClustering);
        timerTickOn();

        newClustering.toClustering();
        newClustering.saveSpace();

        if (data.getNbClustering() == 1) {
            data.setAllClusterings(new Vector());
            myTabbedPane.myPCAPane.currentClustering = 0;
            myTabbedPane.myManifoldPane.currentClustering = 0;
            myTabbedPane.myClusteringPane.activateSearch();
        }

        data.getAllClusterings().add(newClustering);// nullPointerException
        myClusteringProfileFrame.addClusteringLF(newClustering, data.getNbClustering() - 1);

        updateClusteringLF(JAGCTVisualizationPane.stringClustering(data.getNbClustering() - 1));

        JInformationFrame.getInstance().setText("Processing clustering : " + newClustering + "... done");
        myTabbedPane.myClusteringPane.listClusteringPane.setText(myTabbedPane.myClusteringPane.clusteringList());
        timerTickOff();
    }

    public void updateClusteringLF(String s) {
        if (data.getNbClustering() == 1) {
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
            final Clustering cc = (Clustering) data.getAllClusterings().elementAt(0);
            final String[] col = cc.myClusteringAlgorithm.getClusterChoices();
            for (i = 0; i < col.length; i++) {
                myTabbedPane.myPCAPane.clusterChoice.addItem(col[i]);
                myTabbedPane.myManifoldPane.clusterChoice.addItem(col[i]);
            }

            myTabbedPane.myPCAPane.clusterChoice.setEnabled(true);
            myTabbedPane.myManifoldPane.clusterChoice.setEnabled(true);

        } else {
            if (data.getNbClustering() < 1) {
                Matrix.perror("0 clusterings to display !");
            }
            myTabbedPane.myPCAPane.clusteringSelect.addItem(s);
            myTabbedPane.myManifoldPane.clusteringSelect.addItem(s);
            myTabbedPane.myCorrelationPane.clusteringSelect.addItem(s);
        }

        myTabbedPane.myPCAPane.updateStatButtons();
        myTabbedPane.myManifoldPane.updateStatButtons();

        myClusteringFrame.updateClusteringLF();
    }

    public Clustering getClustering(int num) {
        if (data.getAllClusterings() == null) {
            return null;
        }
        if (num < 0) {
            Matrix.perror("Negative value for clustering number !");
        }
        if (num > data.getAllClusterings().size() - 1) {
            Matrix.perror("Clustering number out of range");
        }

        return (Clustering) data.getAllClusterings().elementAt(num);
    }

    public void sendGaruda(File file) {
        apcl.agctbackend.sendFileToGaruda(file);
    }

    public void setSDThreshold(double sd){
        normalizationData.setSDThreshold(sd);
    }
}
