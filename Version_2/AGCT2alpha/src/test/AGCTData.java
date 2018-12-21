package test;

import forDebug.Debug;

import java.io.File;
import java.util.Vector;

public class AGCTData {

    /**
     * ******************************
     * Feature cmputation and selection
     * ***
     */
    private boolean drawFrame;
    private boolean reflectionZ;
    private boolean reflectionX;
    private boolean reflectionY;
    private boolean timerIsTicking;
    private boolean triangulationTimerIsTicking;
    private AGCTTimer myTickingTimer;
    private AGCTTimer myTriangulationTickingTimer;
    private AGCTTimer myScenarioTimer;
    /**
     * int used for Scenario
     */
    private int currentIndexW;
    /**
     * int used for Scenario
     */
    private int currentIndexM;
    /**
     * FileWriter for convenience
     */
    private AGCTFileWriter myAGCTFileWriter;
    /**
     * pointer on the Domain of genes
     */
    private Domain myDomain;// 
    private File myDomainFile;
    private int dimElements;
    /**
     * dimension of all square matrices below:
     */
    /* similarity matrix */
    private Matrix W;
    /**
     * dimension of all square matrices below:
     */
    /* copy of W for clustering*/
    private Matrix CW;
    /**
     * dimension of all square matrices below:
     */
    /* Diagonal matrix build on W*/
    private Matrix D_from_W;
    /**
     * dimension of all square matrices below:
     */
    /* Normalized matrix = f(D, W) */
    private Matrix N_from_W;
    /**
     * dimension of all square matrices below:
     */
    /* toy normalized matrix */
    private Matrix DWD;
    /**
     * dimension of all square matrices below:
     */
    /* M = Manifold eigenvector matrix */
    private Matrix M;
    /**
     * // W = Similarity matrix // CW = Copy of W for clustering; // D_from_W =
     * Diagonal matrix built on W // N_from_W = Normalized matrix = f(D, W) //
     * DWD = toy normalized matrix // M = Manifold eigenvector matrix
     */
    private int dimFeatures;
    /**
     * dimension of all square matrices below:
     */
    private Matrix C;
    /**
     * dimension of all square matrices below:
     */
    private Matrix E;
    /**
     * dimension of all square matrices below:
     */
    private Matrix V;
    /**
     * // C = Correlation matrix // E = Eigenvectors of C // V = Correlation
     * circle variable components (in ROWS)
     */
    private Vector feature_pca_Pnt3D;
    private Pnt3D min_CorrelationPnt3D;
    private Pnt3D max_CorrelationPnt3D;
    private Vector feature_names;
    /**
     * All names (typically, ligand name + # feature if Haar for example)
     */
    private double[] average_features;
    private double[] sigma_features;
    /**
     * average and sigmas for features
     */
    private double[] manifold_eigenvalues;
    private double[] pca_eigenvalues;
//	private int[][] nearest_neighbors;
    /**
     * nearest neighbors for each gene
     */
    private int[][] edges;
    /**
     * ***********************************************************************************
     * clustering
     * ***
     */
    private Vector allClusterings;
    private int nbClustering;
    /**
     * some nice booleans
     */
    private boolean rawDomainExists;
    /**
     * yes iff raw data already loaded
     */
    private boolean annotationsExist;
    /**
     * yes iff raw data already loaded
     */
    private boolean annotationsExistP;
    /**
     * yes iff raw data already loaded
     */
    private boolean annotationsExistF;
    /**
     * yes iff raw data already loaded
     */
    private boolean annotationsExistC;
    /**
     * yes // iff // annotations // loaded
     */
    private boolean isTriangulated;
    /**
     * yes if triangulated
     */
    private boolean pca_computed;
    /**
     * yes iff pca computed
     */
    private boolean soft_clustering_computed;
    /**
     * yes iff a soft membership is computed
     */
    private boolean hard_clustering_computed;

    public AGCTData(boolean drawFrame, boolean reflectionZ,
                    boolean reflectionX, boolean reflectionY, boolean timerIsTicking,
                    boolean triangulationTimerIsTicking) {
        this.drawFrame = drawFrame;
        this.reflectionZ = reflectionZ;
        this.reflectionX = reflectionX;
        this.reflectionY = reflectionY;
        this.timerIsTicking = timerIsTicking;
        this.triangulationTimerIsTicking = triangulationTimerIsTicking;
    }

    public boolean isDrawFrame() {
        return drawFrame;
    }

    public void setDrawFrame(boolean drawFrame) {
        this.drawFrame = drawFrame;
    }

    public boolean isReflectionZ() {
        return reflectionZ;
    }

    public void setReflectionZ(boolean reflectionZ) {
        this.reflectionZ = reflectionZ;
    }

    public boolean isReflectionX() {
        return reflectionX;
    }

    public void setReflectionX(boolean reflectionX) {
        this.reflectionX = reflectionX;
    }

    public boolean isReflectionY() {
        return reflectionY;
    }

    public void setReflectionY(boolean reflectionY) {
        this.reflectionY = reflectionY;
    }

    public boolean isTimerIsTicking() {
        return timerIsTicking;
    }

    public void setTimerIsTicking(boolean timerIsTicking) {
        this.timerIsTicking = timerIsTicking;
    }

    public boolean isTriangulationTimerIsTicking() {
        return triangulationTimerIsTicking;
    }

    public void setTriangulationTimerIsTicking(
            boolean triangulationTimerIsTicking) {
        this.triangulationTimerIsTicking = triangulationTimerIsTicking;
    }

    public AGCTTimer getMyTickingTimer() {
        return myTickingTimer;
    }

    public void setMyTickingTimer(AGCTTimer myTickingTimer) {
        this.myTickingTimer = myTickingTimer;
    }

    public AGCTTimer getMyTriangulationTickingTimer() {
        return myTriangulationTickingTimer;
    }

    public void setMyTriangulationTickingTimer(
            AGCTTimer myTriangulationTickingTimer) {
        this.myTriangulationTickingTimer = myTriangulationTickingTimer;
    }

    public AGCTTimer getMyScenarioTimer() {
        return myScenarioTimer;
    }

    public void setMyScenarioTimer(AGCTTimer myScenarioTimer) {
        this.myScenarioTimer = myScenarioTimer;
    }

    public int getCurrentIndexW() {
        return currentIndexW;
    }

    public void setCurrentIndexW(int currentIndexW) {
        this.currentIndexW = currentIndexW;
    }

    public int getCurrentIndexM() {
        return currentIndexM;
    }

    public void setCurrentIndexM(int currentIndexM) {
        this.currentIndexM = currentIndexM;
    }

    public AGCTFileWriter getMyAGCTFileWriter() {
        return myAGCTFileWriter;
    }

    public void setMyAGCTFileWriter(AGCTFileWriter myAGCTFileWriter) {
        this.myAGCTFileWriter = myAGCTFileWriter;
    }

    public Domain getMyDomain() {
        return myDomain;
    }

    public void setMyDomain(Domain myDomain) {
        this.myDomain = myDomain;
    }

    public File getMyDomainFile() {
        return myDomainFile;
    }

    public void setMyDomainFile(File myDomainFile) {
        this.myDomainFile = myDomainFile;
    }

    /**
     * 選ばれたgene の数と等しい
     *
     * @return
     */
    public int getDimElements() {
        return dimElements;
    }

    public void setDimElements(int dimElements) {
        this.dimElements = dimElements;
    }

    public Matrix getW() {
        return W;
    }

    public void setW(Matrix w) {
        W = w;
    }

    public Matrix getCopyOfSimilarityMatrix() {
        return CW;
    }

    public void setCW(Matrix cW) {
        CW = cW;
    }

    public Matrix getD_from_W() {
        return D_from_W;
    }

    public void setD_from_W(Matrix dFromW) {
        D_from_W = dFromW;
    }

    public Matrix getN_from_W() {
        return N_from_W;
    }

    public void setN_from_W(Matrix nFromW) {
        N_from_W = nFromW;
    }

    public Matrix getDWD() {
        return DWD;
    }

    public void setDWD(Matrix dWD) {
        DWD = dWD;
    }

    public Matrix getM() {
        return M;
    }

    public void setM(Matrix m) {
        M = m;
    }

    public int getDimFeatures() {
        return dimFeatures;
    }

    public void setDimFeatures(int dimFeatures) {
        this.dimFeatures = dimFeatures;
    }

    public Matrix getC() {
        return C;
    }

    public void setC(Matrix c) {
        C = c;
    }

    public Matrix getE() {
        return E;
    }

    public void setE(Matrix e) {
        E = e;
    }

    public Matrix getV() {
        return V;
    }

    public void setV(Matrix v) {
        V = v;
    }

    public Vector getFeature_pca_Pnt3D() {
        return feature_pca_Pnt3D;
    }

    public void setFeature_pca_Pnt3D(Vector featurePcaPnt3D) {
        feature_pca_Pnt3D = featurePcaPnt3D;
    }

    public Pnt3D getMin_CorrelationPnt3D() {
        return min_CorrelationPnt3D;
    }

    public void setMin_CorrelationPnt3D(Pnt3D minCorrelationPnt3D) {
        min_CorrelationPnt3D = minCorrelationPnt3D;
    }

    public Pnt3D getMax_CorrelationPnt3D() {
        return max_CorrelationPnt3D;
    }

    public void setMax_CorrelationPnt3D(Pnt3D maxCorrelationPnt3D) {
        max_CorrelationPnt3D = maxCorrelationPnt3D;
    }

    public Vector getFeature_names() {
        return feature_names;
    }

    public void setFeature_names(Vector featureNames) {
        feature_names = featureNames;
    }

    public double[] getAverage_features() {
        return average_features;
    }

    public void setAverage_features(double[] averageFeatures) {
        average_features = averageFeatures;
    }

    public double[] getSigma_features() {
        return sigma_features;
    }

    public void setSigma_features(double[] sigmaFeatures) {
        sigma_features = sigmaFeatures;
    }

    public double[] getManifold_eigenvalues() {
        return manifold_eigenvalues;
    }

    public void setManifold_eigenvalues(double[] manifoldEigenvalues) {
        manifold_eigenvalues = manifoldEigenvalues;
    }

    public double[] getPca_eigenvalues() {
        return pca_eigenvalues;
    }

    public void setPca_eigenvalues(double[] pcaEigenvalues) {
        pca_eigenvalues = pcaEigenvalues;
    }

    /**
     * nearest_neighbors for each gene.
     *
     * @return
     */
//	public int[][] getNearest_neighbors() {
//		return nearest_neighbors;
//	}
//	public void setNearest_neighbors(int[][] nearestNeighbors) {
//		nearest_neighbors = nearestNeighbors;
//	}
    public int[][] getEdges() {
        return edges;
    }

    public void setEdges(int[][] edges) {
        this.edges = edges;
    }

    public Vector getAllClusterings() {
        return allClusterings;
    }

    public void setAllClusterings(Vector allClusterings) {
        this.allClusterings = allClusterings;
    }

    public int getNbClustering() {
        return nbClustering;
    }

    public void setNbClustering(int nbClustering) {
        this.nbClustering = nbClustering;
    }

    public boolean isRawDomainExists() {
        return rawDomainExists;
    }

    public void setRawDomainExists(boolean rawDomainExists) {
        this.rawDomainExists = rawDomainExists;
    }

    public boolean isAnnotationsExist() {
        return annotationsExist;
    }

    public void setAnnotationsExist(boolean annotationsExist) {
        Debug.debug("AGCTData # setAnnotationsExist", annotationsExist);
        this.annotationsExist = annotationsExist;
    }

    public boolean isAnnotationsExistP() {
        return annotationsExistP;
    }

    public void setAnnotationsExistP(boolean annotationsExistP) {
        this.annotationsExistP = annotationsExistP;
    }

    public boolean isAnnotationsExistF() {
        return annotationsExistF;
    }

    public void setAnnotationsExistF(boolean annotationsExistF) {
        this.annotationsExistF = annotationsExistF;
    }

    public boolean isAnnotationsExistC() {
        return annotationsExistC;
    }

    public void setAnnotationsExistC(boolean annotationsExistC) {
        this.annotationsExistC = annotationsExistC;
    }

    public boolean isTriangulated() {
        return isTriangulated;
    }

    public void setTriangulated(boolean isTriangulated) {
        this.isTriangulated = isTriangulated;
    }

    public boolean isPca_computed() {
        return pca_computed;
    }

    public void setPca_computed(boolean pcaComputed) {
        pca_computed = pcaComputed;
    }

    public boolean isSoft_clustering_computed() {
        return soft_clustering_computed;
    }

    public void setSoft_clustering_computed(boolean softClusteringComputed) {
        soft_clustering_computed = softClusteringComputed;
    }

    public boolean isHard_clustering_computed() {
        return hard_clustering_computed;
    }

    public void setHard_clustering_computed(boolean hardClusteringComputed) {
        hard_clustering_computed = hardClusteringComputed;
    }


    public boolean getAnnotationsExist() {
        return annotationsExist;
    }
}