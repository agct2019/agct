package gene;

import forDebug.Debug;
import ligand.LigandList;
import matrix.MatrixConfig;
import test.*;

import java.util.*;

public class Gene implements Debuggable, Comparable<Gene> {

    public static final String Not_defined = "Not_defined";
    public static final String Gene$Ontology$Molecular$Function = "Gene Ontology Molecular Function";
    public static final String Gene$Ontology$Cellular$Component = "Gene Ontology Cellular Component";
    public static final String Gene$Ontology$Biological$Process = "Gene Ontology Biological Process";
    public static final String RefSeq$Transcript$ID = "RefSeq Transcript ID";
    public static final String ENTREZ_GENE_ID = "ENTREZ_GENE_ID";
    public static final String Gene$Symbol = "Gene Symbol";
    public static final String Gene$Title = "Gene Title";
    public static final String Representative$Public$ID = "Representative Public ID";
    public static final String Target$Description = "Target Description";
    public static final String Sequence$Source = "Sequence Source";
    public static final String Sequence$Type = "Sequence Type";
    public static final String Annotation$Date = "Annotation Date";
    public static final String Species$Scientific$Name = "Species Scientific Name";
    public static final String SPOT_ID = "SPOT_ID";
    public static final String GB_ACC = "GB_ACC";
    public static final String ID = "ID";
    public static final List<String> ANNO_LIST = Collections.unmodifiableList(Arrays.asList(new String[]{
            Gene.ID,
            Gene.GB_ACC,
            Gene.SPOT_ID,
            Gene.Species$Scientific$Name,
            Gene.Annotation$Date,
            Gene.Sequence$Type,
            Gene.Sequence$Source,
            Gene.Target$Description,
            Gene.Representative$Public$ID,
            Gene.Gene$Title,
            Gene.Gene$Symbol,
            Gene.ENTREZ_GENE_ID,
            Gene.RefSeq$Transcript$ID,
            Gene.Gene$Ontology$Biological$Process,
            Gene.Gene$Ontology$Cellular$Component,
            Gene.Gene$Ontology$Molecular$Function
    }));
    private static GeneList geneList;

    public static GeneList getGeneList() {
        return geneList;
    }

    public static void setGeneList(GeneList geneList) {
        Gene.geneList = geneList;
    }

    private int idNumber;
    //	private int clusterNumber;
    private int nonVisCount = 0;// 0 iff visible

    public void setVisible(boolean flag) {
        if (flag) {
            nonVisCount++;
        } else {
            nonVisCount--;
        }
//        if(nonVisCount < 0)throw new AssertionError();
    }

    private int visCount = 0;// visible iff >0.

    public void setVisibleForOr(boolean b) {
        if (b) {
            visCount++;
        } else {
            visCount--;
        }
//        if(visCount < 0)throw new AssertionError();
    }

    private static boolean andMode = true;

    public static void setAndMode(boolean andMode) {
        Gene.andMode = andMode;
    }

    static int callIsVisible = 0;
    public boolean isVisible() {
        callIsVisible++;
        if(callIsVisible >= 1e7){
            Debug.debug("Too many call of \"Gene#isVisible\". BEGIN StackTrace:");
            Debug.debug(Thread.currentThread().getStackTrace());
            Debug.debug("END StackTrace");
            callIsVisible = 0;
        }
//        Debug.debug("Gene#isVisible", andMode, nonVisCount, visCount);
        if (andMode) {
            return nonVisCount == 0;
        } else {
            return visCount > 0;
        }
    }

    public int peakTimeId(int doseId) {
        return features.get(doseId).calcPeakTimeId();
    }

    public double peakInitialValue(int doseId) {
        int timeId = features.get(doseId).calcPeakTimeId();
        return features.get(doseId).getInitialActivation(timeId);
    }

    public int getClusterNumber(int experimentId) {
        int res = -1;
        double max = -1;
        try {
            Double[] Ds = getClusterMemberships(experimentId);
            for (int i = 0; i < Ds.length; i++) {
                if (Ds[i] > max) {
                    max = Ds[i];
                    res = i;
                }
            }
        } catch (NullPointerException e) {
            res = -1;
        }
//		DoDebug.debug("Gene :: getClusterNumber() :: end",res);
        return res;
    }

//    public int getClusterNumber() {
////		DoDebug.debug("Gene :: getClusterNumber() :: start");
//        int res = -1;
//        double max = -1;
//        try {
//            Double[] Ds = getClusterMemberships(0);
//            for (int i = 0; i < Ds.length; i++) {
//                if (Ds[i] > max) {
//                    max = Ds[i];
//                    res = i;
//                }
//            }
//        } catch (NullPointerException e) {
//            res = -1;
//        }
////		DoDebug.debug("Gene :: getClusterNumber() :: end",res);
//        return res;
//    }

    //	public void setClusterNumber(int _clusterNumber) {
    //		this.clusterNumber = _clusterNumber;
    //	}
    public void setMyAGCT(AGCT agct) {
        myAGCT = agct;
    }

    @Override
    public int hashCode() {
        // TODO Auto-generated method stub
        return myID.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        // TODO Auto-generated method stub
        return myID.equals(((Gene) obj).myID);
    }

    public int getIdNumber() {
        return idNumber;
    }

    public void setIdNumber(int idNumber) {
        this.idNumber = idNumber;
    }

    private LigandList ligands;

    public void setLigands(LigandList ligands) {
        this.ligands = ligands;
    }

    // private int responsiveCount = 0;
    // private boolean isResponsible = false;
    private boolean isSelected = true;
    private boolean isEnabled = true;
    public final String myID;

    public boolean isEnabled() {
        return isEnabled;
    }

    public void setEnabled(boolean isEnabled) {
        this.isEnabled = isEnabled;
    }

    public boolean isSelected() {
        return isSelected;
    }

    public void setSelected(boolean isSelected) {
        this.isSelected = isSelected;
    }

    public boolean isResponsive() {
        return responsiveMask > 0;
    }

    public int getResponsiveFileId() {
        for (int i = 0; i < 32; i++) {
            if (((responsiveMask >> i) & 1) == 1) {
                return i;
            }
        }
        if (AGCT.MYDEBUG) {
            AGCT.debug("ERROR at Gene.getResponsiveFileId");
        }
        throw new RuntimeException();
    }

    private int responsiveMask = 0;

    public void setResponsive(int id, boolean set) {
        assert 0 <= id && id < 32;
        if (set) {
            responsiveMask |= 1 << id;
        } else {
            responsiveMask &= ~(1 << id);
        }
    }

    public int getNumberOfDoses() {
        return features.size();
    }

    public boolean hasSamePeak(BitSet bitmask) {
        HashSet<Integer> set = new HashSet<Integer>();
        for (int i = 0; i < features.size(); i++) {
            if (bitmask.get(i)) {
                set.add(features.get(i).calcPeakTimeId());
            }
        }
        return set.size() <= 1;
    }

    public int getNumberOfLigand() {
        return rawCoordinates.size();
//		return initialCoordinates.size();
    }

    public static Hashtable FinalIndex;
    // returns an Integer [][] with Integers X et Y
    // Keys are Integers = indexes
    private static boolean FinalIndex_computed = false;
    public AGCT myAGCT;
    static Domain myDomain;

    public static void setMyDomain(Domain d) {
        myDomain = d;
    }

    // public void setName(String name){
    // this.name = name;
    // }
    public String getName() {
        return myID;
    }

    String lightAnnotation;
    public String asciiName;
    public Vector annotationsF;
    public Vector annotationsP;
    public Vector annotationsC;
    public boolean annotationsFilled;
    private boolean referenced; // highlighted.

    public boolean isReferenced() {
        return referenced && typeReferenced >= 0 && ((AGCT.HighlightMask >> typeReferenced) & 1) == 1;
    }

    public void setReferenced(boolean referenced) {
        System.out.println("Gene.setReferenced() " + referenced);
        this.referenced = referenced;
    }

    private int typeReferenced; // group.
    public ArrayList<Feature> features;

    public void makeFeatures() {
        features = new ArrayList<Feature>();
        for (int i = 0; i < rawCoordinates.size(); i++) {
            double[][] ds = input.get(i);
            double[] initialActivation = new double[ds.length];
            for (int j = 0; j < ds.length; j++) {
                initialActivation[j] = Util_Math.calcMean(ds[j], AGCT.getInstance().getNormalizationData().getNormalizeMode());
            }
            Feature feature = new Feature(rawCoordinates.get(i), initialActivation);
            features.add(feature);
        }
    }

    public int getTypeReferenced() {
        return typeReferenced;
    }

    public void setTypeReferenced(int typeReferenced) {
        this.typeReferenced = typeReferenced;
    }

    /**
     * ligandId,timeId,replicateId -> activation.
     */
    private ArrayList<double[][]> input = new ArrayList<double[][]>();
    // ligandId,timeId,replicateId -> activation.
    /**
     * normalize後のcordinates
     */
    private ArrayList<double[]> rawCoordinates = new ArrayList<double[]>();
    /**
     * 最も始めのcondition
     */
    private ArrayList<double[]> initialCoordinates = new ArrayList<double[]>();

    // rawCoordinates.get(ligandId)[timeId] -> activation.
    public void addInput(double[][] vals) {
        input.add(vals);
    }

    /**
     * double[][] : time,rep -> activation.
     *
     * @param ligandId
     * @return time, rep -> activation
     */
    public double[][] getInput(int ligandId) {
        return input.get(ligandId);
    }

    public void addRawCoordinate(double[] vals) {
        rawCoordinates.add(vals);
    }

    public void addInitialCoordinate(double[] vals) {
        initialCoordinates.add(vals);
    }

    public double[] getRawCoordinate(int ligandId) {
        return rawCoordinates.get(ligandId);
    }

    public double[] getInitialCoordinate(int ligandId) {
        return initialCoordinates.get(ligandId);
    }

    // private boolean[][] undeterminedCoordinates;
    private ArrayList<boolean[]> undeterminedCoordinates = new ArrayList<boolean[]>();

    // yes if raw coordinate is "?" for (X, Y)
    public void addUndeterminedCoordinate(boolean[] vals) {
        undeterminedCoordinates.add(vals);
    }
    public boolean[] getUndeterminedCoordinate(int id) {
        /* This function always returns an array of false. Try this:
    	for (boolean b : undeterminedCoordinates.get(id)) {
    		if (b) {
            	System.err.println("UndeterminedCoordinate wasn't false.");
            	System.exit(1);
    		}
    	}
    	*/
        return undeterminedCoordinates.get(id);
    }

    // coordinates used after processing with AGCT.Method_F
    public double[][] finalCoordinates;
    public Pnt finalCoordinates_1D;
    Vector cluster_memberships;
    // each element = Double[] with the cluster memberships

    /*
     * sum of absolute difference
     */
    public double totalDistance;

    public double calcTotalDistance() {
        totalDistance = 0;
        for (double[] ds : initialCoordinates) {
            for (int i = 1; i < ds.length; i++) {
                if (!Double.isNaN(ds[i - 1]) && !Double.isNaN(ds[i])) {
                    totalDistance += Math.abs(ds[i] - ds[i - 1]);
                }
            }
        }
        return totalDistance;
    }

    /**
     * 　copy per cellの　geometric mean を求める
     *
     * @return　geometric mean
     */
    public double calcAverageCopyPerCell() {
        double sum = 1;
        int cnt = 0;
        ArrayList<Double> list = new ArrayList<Double>();
        for (double[] ds : initialCoordinates) {
            for (int t = 0; t < ds.length; t++) {
                if (!Double.isNaN(ds[t])) {
                    sum *= ds[t];
                    list.add(ds[t]);
                    cnt++;
                }
            }
        }
        //		DoDebug.debug("Gene : copy per cell sum",list,sum,cnt,Math.pow(sum,1./cnt));
        return Math.pow(sum, 1. / cnt);
    }

    /**
     * doseId 3,4 (any time) の中で最小のcopy per cellを求める
     *
     * @return
     */
    public double min34CopyPerCell() {
        double res = 0;
        for (int i = 0; i < initialCoordinates.size(); i++) {//
            if (i == 3 || i == 4) {
                double[] ds = initialCoordinates.get(i);
                for (int t = 0; t < ds.length; i++) {
                    if (!Double.isNaN(ds[t])) {
                        res = Math.max(res, ds[t]);
                    }
                }
            }
        }
        return res;
    }

    public boolean doContainsSeriesOfAbove1() {
        for (double[] ds : rawCoordinates) {
            boolean befAbove1 = false;
            for (double d : ds) {
                if (d >= 1) {
                    if (befAbove1) {
                        return true;
                    }
                    befAbove1 = true;
                } else {
                    befAbove1 = false;
                }
            }
        }
        return false;
    }

    /**
     * *******************************************************************************
     * Components after processing
     * ***
     */
    /*
     *  Contains Integers = index among selected Genes of neighbors in the
     *  triangulation
     */
    public Vector neighbors = new Vector();
    /*Contains the angle wrt neighbors*/
    public Vector neighbor_angles;
    // components on the manifold for the Gene
    public Pnt manifold_components_total;
    /**
     * fillManifoldComponentsで埋められる。manifoldの座標の先頭
     * triangulationに使われる次元(d) をとってきたもの。
     */
    public Pnt manifold_components_triangulation;
    // components used for the triangulation
    public Pnt manifold_point3D;// ４次元のPoint
    // idem (in 4D), but with a Pnt3D, should be deprecated (no more scaling)
    public Pnt3D manifold_Pnt3D;// Pnt3Dのインスタンス
    public Pnt pca_components;
    // principal components for the Gene
    public Pnt pca_point3D;
    // idem, but with a Pnt3D
    public Pnt3D pca_Pnt3D;
    boolean selected, manifold_components_computed, pca_components_computed;
    public String geneTitle;
    public String refSeq;

    //	public String Symbol;
    public static void test_computeFinalIndex(double[][] finalC) {
        if (!FinalIndex_computed) {
            computeFinalIndex(finalC);
        }
    }

    public static void flushAll() {
        FinalIndex = null;
        FinalIndex_computed = false;
    }

    public static void computeFinalIndex(double[][] finalC) {
        // on the basis of finalC, computes the Hashtables

        Integer[] I;
        int rx, ry, i = 0;
        boolean arret = false;

        FinalIndex_computed = false;
        FinalIndex = null;
        FinalIndex = new Hashtable();

        rx = ry = 0;
        do {
            I = new Integer[2];
            I[0] = new Integer(rx);
            I[1] = new Integer(ry);

            FinalIndex.put(new Integer(i), I);

            ry++;
            if (ry == finalC[rx].length) {
                rx++;
                ry = 0;

                if (rx == finalC.length) {
                    arret = true;
                }
            }
            i++;

        } while (!arret);

        FinalIndex_computed = true;
    }

    public static void printFinalIndex() {
        Enumeration extensions = FinalIndex.keys();
        Integer I;
        if (extensions != null) {
            while (extensions.hasMoreElements()) {
                I = (Integer) extensions.nextElement();
                System.out.print(I + " -> " + ((Integer[]) FinalIndex.get(I))[0] + ", " + ((Integer[]) FinalIndex.get(I))[1] + "; ");
            }
        }
    }

    public String toLightHTMLString() {
        String val = "<HTML>", sc, sn;

        if (myID != null) {
            val += myID;
            if ((typeReferenced == -1) || (myAGCT.data.getMyDomain().orderedGroups == null)) {
                if (typeReferenced != -1) {
                    val += "(referenced as " + typeReferenced + ")";
                }
            } else {
                sc = Integer.toHexString(JAGCTGraphicsPanel.referencedColors[typeReferenced].getRGB());
                sn = sc.substring(2, sc.length());
                val += "(referenced as <FONT COLOR=#" + sn + ">";
                // val +=
                // myAGCT.data.getMyDomain().orderedGroups.get(typeReferenced);
                val += "</FONT>)";
            }
        }
        if ((asciiName != null) && (asciiName != "")) {
            val += "_" + asciiName + "";
        }
        val += "</HTML>";
        return val;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(myID);
        sb.append("\t");
        if (annotationsF != null) {
//            sb.append("F:");
            for (int i = 0; i < annotationsF.size(); i++) {
                sb.append(" " + ((String[]) annotationsF.get(i))[0] + "  " + ((String[]) annotationsF.get(i))[1]);
            }
        }
        sb.append("\t");
        if (annotationsP != null) {
//            sb.append("P:");
            for (int i = 0; i < annotationsP.size(); i++) {
                sb.append(" " + ((String[]) annotationsP.get(i))[0] + "  " + ((String[]) annotationsP.get(i))[1]);
            }
        }
        sb.append("\t");
        if (annotationsC != null) {
//            sb.append("C:");
            for (int i = 0; i < annotationsC.size(); i++) {
                sb.append(" " + ((String[]) annotationsC.get(i))[0] + "  " + ((String[]) annotationsC.get(i))[1]);
            }
        }
        sb.append("\t");
        if ((lightAnnotation != null) && (lightAnnotation != "")) {
            sb.append("light Annotation (in gene file): " + lightAnnotation);
        }
        return sb.toString();
    }

    public String toString2() {
        String val = "";
        int i, j;
        Gene nn;
        Double[] cc;

        if (myID != null) {
            val += myID + "(" + referenced + "|" + typeReferenced + ")";
        }
        if ((asciiName != null) && (asciiName != "")) {
            val += "_" + asciiName + "";
        }
        val += "\n";

        if (annotationsF != null) {
            val += " -- Annotations (F):\n";
            for (i = 0; i < annotationsF.size(); i++) {
                val += " " + ((String[]) annotationsF.get(i))[0] + "  " + ((String[]) annotationsF.get(i))[1] + "\n";
            }
            val += "\n";
        }
        if (annotationsP != null) {
            val += " -- Annotations (P):\n";
            for (i = 0; i < annotationsP.size(); i++) {
                val += " " + ((String[]) annotationsP.get(i))[0] + "  " + ((String[]) annotationsP.get(i))[1] + "\n";
            }
            val += "\n";
        }
        if (annotationsC != null) {
            val += " -- Annotations (C):\n";
            for (i = 0; i < annotationsC.size(); i++) {
                val += " " + ((String[]) annotationsC.get(i))[0] + "  " + ((String[]) annotationsC.get(i))[1] + "\n";
            }
            val += "\n";
        }

        if (neighbors != null) {
            val += " -- Natural Neighbors:\n";
            for (i = 0; i < neighbors.size(); i++) {
                nn = (Gene) myDomain.getGenes().get(myDomain.selectedGeneNumberToGeneNumber[((Integer) neighbors.get(i)).intValue()]);
                val += " " + nn.myID;
                if ((nn.asciiName != null) && (nn.asciiName != "")) {
                    val += "_" + nn.asciiName + "";
                }
                val += "\n";
            }
        }

        if (cluster_memberships != null) {
            val += " -- Cluster memberships:\n";
            for (i = 0; i < cluster_memberships.size(); i++) {
                cc = getClusterMemberships(i);
                val += " - Memberships [Clustering #" + i + "]: ";
                for (j = 0; j < cc.length; j++) {
                    val += " Cluster_" + j + " : " + DF.format(getClusterMemberships(i, j)) + " --";
                }
                val += "\n";
            }
        }

        if ((lightAnnotation != null) && (lightAnnotation != "")) {
            val += "\n -- light Annotation (in gene file): " + lightAnnotation + "\n";
        }

        val += "\n\n";

        return val;
    }

    public Double[] getClusterMemberships(int i) {
        return (Double[]) cluster_memberships.get(i);
    }

    public double getClusterMemberships(int i, int j) {
        return getClusterMemberships(i)[j].doubleValue();
    }

    public void addClusterMemberships(Double[] cm) {
        if (cluster_memberships == null) {
            cluster_memberships = new Vector();
        }
        cluster_memberships.add(cm);
    }

    public void putReferenced() {
        referenced = true;
    }

    public boolean geneOrNeighborsIsReferencedDelaunay() {
        int j;
        Gene nn;
        if (neighbors != null) {
            for (j = 0; j < neighbors.size(); j++) {
                nn = (Gene) myDomain.getGenes().get(myDomain.selectedGeneNumberToGeneNumber[((Integer) neighbors.get(j)).intValue()]);
                if (nn.referenced) {
                    return true;
                }
            }
        }
        return referenced;
    }

    public double getMaxDelaunay() {
        if (neighbor_angles == null) {
            return Double.NaN;
        }
        double res = -100;
        for (int i = 0; i < neighbor_angles.size(); i++) {
            double nangle = ((Double) neighbor_angles.get(i)).doubleValue();
            res = Math.max(res, nangle);
        }
        return res;
    }

    Gene(String ID) {
        this.myID = ID;
        setAnnotation(Gene.ID, this.myID);
        //		this.Symbol=this.myID;
        setAnnotation(Gene.Gene$Symbol, this.myID);
        typeReferenced = -1;
        lightAnnotation = "";
        manifold_components_triangulation = new Pnt(AGCT.Number_Of_Triangulation_Dimensions);
        manifold_components_total = new Pnt(AGCT.Number_Of_Manifold_Components);
        manifold_point3D = new Pnt(4);
        pca_point3D = new Pnt(4);
        manifold_Pnt3D = new Pnt3D();
        pca_Pnt3D = new Pnt3D();
    }

    public static boolean dimCoincide(Gene a, Gene b) {
        return (a.finalCoordinates.length == b.finalCoordinates.length);
    }

    /**
     * || a - b ||^2
     * ユークリッド距離の二乗
     *
     * @param a
     * @param b
     * @return
     */
    public static double squaredMagnitudeDifference(Gene a, Gene b) {
        double v = 0.0, fca, fcb;
        int i, dim = Util_Feature.getNumberOfFeaturesAfterSelection();
        if (!dimCoincide(a, b)) {
            Matrix.perror("Dimension mismatch between genes " + a.myID + " and " + b.myID);
        }

        for (i = 0; i < dim; i++) {
            fca = a.getFinalCoordinates(i, true);
            fcb = b.getFinalCoordinates(i, true);
            v += ((fca - fcb) * (fca - fcb));
        }
        return v;
    }

    private double norm2() {
        // returns || a ||^2
        double v = 0.0;
        int i, dim = Util_Feature.getNumberOfFeaturesAfterSelection();
        for (i = 0; i < dim; i++) {
            v += (getFinalCoordinates(i, true) * getFinalCoordinates(i, true));
        }
        return v;
    }


    /**
     * return s <a, b>
     *
     * @param b
     * @return
     */
    public double dot(Gene b) {
        double v = 0.0;
        int i, dim = Util_Feature.getNumberOfFeaturesAfterSelection();
        if (!dimCoincide(this, b)) {
            Matrix.perror("Dimension mismatch between genes " + myID + " and " + b.myID);
        }

        // calc step
        for (i = 0; i < dim; i++) {
            v += (getFinalCoordinates(i, true) * b.getFinalCoordinates(i, true));
        }
        return v;
    }

    public static Gene center(GeneList genes) {
        int dim = Util_Feature.getNumberOfFeaturesAfterSelection();
        Gene res = new Gene("");
        res.finalCoordinates = new double[genes.get(0).finalCoordinates.length][genes.get(0).finalCoordinates[0].length];
        for (int i = 0; i < dim; i++) {
            double sum = 0;
            for (Gene gene : genes) {
                sum += gene.getFinalCoordinatesAfterSelection(i);
            }
            res.setFinalCoordinatesAfterSelection(i, sum / geneList.size());
        }
        return res;
    }

    /**
     * | (a,b)/(||a||*||b||)) |
     */
    public static double getAbsoluteCosineSimilarity(Gene a, Gene b) {
        if (!dimCoincide(a, b)) {
            Matrix.perror("Dimension mismatch between genes " + a.myID + " and " + b.myID);
        }

        double ps = a.dot(b);
        double dt1 = Math.sqrt(a.norm2());
        double dt2 = Math.sqrt(b.norm2());
        double val = 0.0;
        double linf = 0.0;
        double lsup = 1.0;

        if ((dt1 == 0.0) || (dt2 == 0.0)) {
            val = 0.0;
        } else {
            val = Math.abs(((ps) / (dt1 * dt2)));

            if ((dt1 < 0.0) || (dt2 < 0.0)) {
                Matrix.perror("Negative norms !");
            }
        }

        if ((val < linf) && (val > -MatrixConfig.Precision_For_Eigensystems)) {
            val = linf;
        }
        if ((val > lsup) && (val < lsup + MatrixConfig.Precision_For_Eigensystems)) {
            val = lsup;
        }
        if ((val < linf) || (val > lsup)) {
            Matrix.perror("Absolute cosine similarity outside bounds");
        }
        return val;
    }

    public static double cosine(Gene a, Gene b) {
        if (!dimCoincide(a, b)) {
            Matrix.perror("Dimension mismatch between genes " + a.myID + " and " + b.myID);
        }

        double ps = a.dot(b);
        double dt1 = Math.sqrt(a.norm2());
        double dt2 = Math.sqrt(b.norm2());
        double val = 0.0;
        double linf = -1.0;
        double lsup = 1.0;

        if ((dt1 == 0.0) || (dt2 == 0.0)) {
            val = 0.0;
        } else {
            val = ((ps) / (dt1 * dt2));

            if ((dt1 < 0.0) || (dt2 < 0.0)) {
                Matrix.perror("Negative norms !");
            }
        }

        if ((val < linf) && (val > linf - MatrixConfig.Precision_For_Cosine)) {
            val = linf;
        }
        if ((val > lsup) && (val < lsup + MatrixConfig.Precision_For_Cosine)) {
            val = lsup;
        }
        if ((val < linf) || (val > lsup)) {
            Matrix.perror("cosine " + val + " outside bounds");
        }
        return val;
    }

    /**
     * (1+(a,b)/(||a||*||b||))/2
     */
    public static double getCosineSimilarity(Gene a, Gene b) {
        if (!dimCoincide(a, b)) {
            Matrix.perror("Dimension mismatch between genes " + a.myID + " and " + b.myID);
        }

        double ps = a.dot(b);
        double dt1 = Math.sqrt(a.norm2());
        double dt2 = Math.sqrt(b.norm2());
        double val = 0.0;
        double linf = 0.0;
        double lsup = 1.0;


        if ((dt1 == 0.0) || (dt2 == 0.0)) {
            val = 0.0;
        } else {
            val = (1.0 + ((ps) / (dt1 * dt2))) / 2.0;
            assert dt1 >= 0 && dt2 >= 0 : "Negative norms !";
        }

        if ((val < linf) && (val > -MatrixConfig.Precision_For_Eigensystems)) {
            val = 0.0;
        }
        if ((val > lsup) && (val < lsup + MatrixConfig.Precision_For_Eigensystems)) {
            val = lsup;
        }
        assert 0 <= val && val <= 1 : "outside bounds.";
        return val;
    }

    public static double getSimilarity(Gene a, Gene b) {
        // double v = 0.0;
        if (AGCT.Method_S == 0) // e^( - ||(a-b)||^2/T )
        {
            return Math.exp(-squaredMagnitudeDifference(a, b) / AGCT.T_Heat_Kernel);
        } else if (AGCT.Method_S == 1) {
            return getCosineSimilarity(a, b);
        } else if (AGCT.Method_S == 2) {// default
            return getAbsoluteCosineSimilarity(a, b);
        } else {
            throw new IllegalArgumentException();
        }
    }

    public String coordinatesString() {
        // we assume rawCoordinates, undeterminedCoordinates AND
        // finalCoordinates computed
        // we also assume ligands selected ARE computed
        String val = " Gene " + myID + " ---\n";
        int i, j;
        val += "Raw Coordinates by ligand : \n";
        for (i = 0; i < rawCoordinates.size(); i++) {
            val += " [ " + i + " : ";
            for (j = 0; j < rawCoordinates.get(i).length; j++) {
                if (getUndeterminedCoordinate(i)[j] == false) {
                    val += Gene.DF.format(getRawCoordinate(i)[j]);
                } else {
                    val += "(?)";
                }
                if (j < rawCoordinates.get(i).length - 1) {
                    val += ", ";
                }
            }
            val += "]\n";
        }
        val += "\n\nFinal Coordinates by ligand : \n";
        for (i = 0; i < finalCoordinates.length; i++) {
            val += " [ " + i + " : ";
            for (j = 0; j < finalCoordinates[i].length; j++) {
                val += Gene.DF.format(finalCoordinates[i][j]);
                if (j < finalCoordinates[i].length - 1) {
                    val += ", ";
                }
            }
            val += "]\n";
        }

        return val;
    }

    double[] finalCoordinatesAfterSelectionMemo;
    double[] finalCoordinatesBeforeSelectionMemo;


    /**
     * @param x
     * @param after
     * @return
     */
    public double getFinalCoordinates(int x, boolean after) {
        // TODO afterが何か調べる。

        if (after) {
            if (finalCoordinatesAfterSelectionMemo == null) {
                int n = Util_Feature.featureIndexToNonTrashedIndex.size() + 1;
                finalCoordinatesAfterSelectionMemo = new double[n];
                Arrays.fill(finalCoordinatesAfterSelectionMemo, 0xDEADBEEF);
            }
            if (finalCoordinatesAfterSelectionMemo[x] != 0xDEADBEEF) {
                return finalCoordinatesAfterSelectionMemo[x];
            }
            return finalCoordinatesAfterSelectionMemo[x] = getFinalCoordinatesAfterSelection(x);
        } else {
            if (finalCoordinatesBeforeSelectionMemo == null) {
                int n = Util_Feature.featureIndexToNonTrashedIndex.size() + 1;
                finalCoordinatesBeforeSelectionMemo = new double[n];
                Arrays.fill(finalCoordinatesBeforeSelectionMemo, 0xDEADBEEF);
            }
            if (finalCoordinatesBeforeSelectionMemo[x] != 0xDEADBEEF) {
                return finalCoordinatesBeforeSelectionMemo[x];
            }
            return finalCoordinatesBeforeSelectionMemo[x] = getFinalCoordinatesBeforeSelection(x);// この構文初めて見た
        }
    }

    /**
     * Pnt finalCoordinates_1Dを計算している。
     * ligandの数の次元で、各値は、featureのデータを表現している。
     */
    public void computeFinalCoordinates_1D() {
        //		System.out.println("Gene.computeFinalCoordinates_1D()");
        // if (AGCT.MYDEBUG)
        // AGCT.debug("Gene.computeFinalCoordinates_1D()");
        int i, dim = Util_Feature.getNumberOfFeaturesAfterSelection();
        finalCoordinates_1D = new Pnt(dim);
        for (i = 0; i < dim; i++) {
            finalCoordinates_1D.coordinates[i] = getFinalCoordinates(i, true);
        }
    }

    /**
     * @param x
     * @return
     */
    public double getFinalCoordinatesBeforeSelection(int x) {
        // new, uses Hashtable FinalIndex
        Integer i = new Integer(x);
        if (!Gene.FinalIndex_computed) {
            Matrix.perror("FinalIndex not computed in Gene (perhaps a Thread problem ?)");
        }

        return finalCoordinates[((Integer[]) Gene.FinalIndex.get(i))[0]][((Integer[]) Gene.FinalIndex.get(i))[1]];
    }

    public double getFinalCoordinatesAfterSelection(int x) {
        // new, uses Hashtable FinalIndex
        Integer i = (Integer) Util_Feature.featureIndexToNonTrashedIndex.get(new Integer(x));
        if (!Gene.FinalIndex_computed) {
            Matrix.perror("FinalIndex not computed in Gene (perhaps a Thread problem ?)");
        }

        if (x > Util_Feature.featureIndexToNonTrashedIndex.size()) {
            Matrix.perror("feature coordinate outside bounds");
        }

        return finalCoordinates[((Integer[]) Gene.FinalIndex.get(i))[0]][((Integer[]) Gene.FinalIndex.get(i))[1]];
    }

    /**
     * @param j
     * @param val
     */
    private void setFinalCoordinatesAfterSelection(int j, double val) {
        Integer i = (Integer) Util_Feature.featureIndexToNonTrashedIndex.get(new Integer(j));
        finalCoordinates[((Integer[]) Gene.FinalIndex.get(i))[0]][((Integer[]) Gene.FinalIndex.get(i))[1]] = val;
    }

    /**
     * myAGCT = a;
     * myDomain = d; 	とするだけ。
     */
    public void finishVariables(AGCT a, Domain d) {
        myAGCT = a;
        myDomain = d;
    }

    public void putLightAnnotation(String s) {
        lightAnnotation = s;
    }

    public void addLightAnnotation(String s) {
        lightAnnotation += s;
    }

    public void addAnnotation(boolean f, boolean p, boolean c, String s) {
        String[] an;

        if ((f && p) || (f && c) || (p && c)) {
            Matrix.perror("Gene.class :: at least two of F, P and C annotations are true for " + s);
        }
        if ((!f) && (!p) && (!c)) {
            Matrix.perror("Gene.class :: F, P and C annotations are false for " + s);
        } else {
            if ((f) && (annotationsF == null)) {
                annotationsF = new Vector();
            } else if ((p) && (annotationsP == null)) {
                annotationsP = new Vector();
            } else if ((c) && (annotationsC == null)) {
                annotationsC = new Vector();
            }

            an = new String[2];
            if (f) {
                an[0] = Domain.ANNOTATION_F;
            } else if (p) {
                an[0] = Domain.ANNOTATION_P;
            } else {
                an[0] = Domain.ANNOTATION_C;
            }

            an[1] = s;
            if (f) {
                annotationsF.add(an);
            } else if (p) {
                annotationsP.add(an);
            }
            if (c) {
                annotationsC.add(an);
            }

            annotationsFilled = true;
        }
    }

    /**
     * 単になくなったリガンドの部分をつめているだけ。
     */
    public void recomputeSelectedFeatures() {
        selected = true;

        ArrayList<double[]> newRawCoordinates = new ArrayList<double[]>();
        // double[][] rc2 = new double[ligands.numberOfSelectedLigands()][];
        ArrayList<boolean[]> newUndeterminedCoordinates = new ArrayList<boolean[]>();
        // boolean[][] uc2 = new boolean[ligands.numberOfSelectedLigands()][];

        // int index = 0;
        for (int i = 0; i < myDomain.getLigands().size(); i++) {
            if (myDomain.getLigands().get(i).isChecked()) {
                // rc2[index] = new double[myDomain.times.length];
                // uc2[index] = new boolean[myDomain.times.length];
                newRawCoordinates.add(getRawCoordinate(i));
                newUndeterminedCoordinates.add(getUndeterminedCoordinate(i));
            }
        }

        rawCoordinates = newRawCoordinates;
        undeterminedCoordinates = newUndeterminedCoordinates;
        finalCoordinates = new double[ligands.numberOfSelectedLigands()][];
    }

    public void fillManifoldComponents(Matrix m, int col) {
        int i;
        //		DoDebug.debug("Gene.fillManifoldComponents()");
        //		DoDebug.debug("manifold_components_total.coordinates.length",
        //				manifold_components_total.coordinates.length);
        for (i = 0; i < manifold_components_total.coordinates.length; i++) {
            if (i + 1 < m.rowCount()) {
                manifold_components_total.coordinates[i] = m.get(i + 1, col);
            }
            if (i < 3) {
                manifold_point3D.coordinates[i] = manifold_components_total.coordinates[i];
                manifold_Pnt3D.coordinates[i] = manifold_components_total.coordinates[i];
            }
            if (i < AGCT.Number_Of_Triangulation_Dimensions) {
                manifold_components_triangulation.coordinates[i] = manifold_components_total.coordinates[i];
            }
        }
        manifold_point3D.coordinates[3] = 1.0;
        manifold_components_computed = true;
        Debug.fdebug("tmp/Gene/maniPoints", manifold_components_total.coordinates);

    }

    public void fillPCAComponents(Matrix E, double[] average, double[] sigma) {
        double acp, slo, ave, sig, curc;
        int y, z;

        pca_components = new Pnt(myAGCT.data.getDimFeatures());

        for (y = 0; y < myAGCT.data.getDimFeatures(); y++) {
            pca_components.coordinates[y] = 0.0;
        }

        for (y = 0; y < myAGCT.data.getDimFeatures(); y++) {
            for (z = 0; z < myAGCT.data.getDimFeatures(); z++) {
                acp = E.get(y, z);
                slo = getFinalCoordinates(z, true);
                ave = average[z];
                sig = sigma[z];
                curc = pca_components.coordinates[y];
                if (sig > 0.0) {
                    curc += (acp * ((slo - ave) / sig));
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

    public void updateManifold_Pnt3D(int xAxis, int yAxis, int zAxis) {
        manifold_Pnt3D.coordinates[0] = manifold_components_total.coordinates[xAxis];
        manifold_Pnt3D.coordinates[1] = manifold_components_total.coordinates[yAxis];
        manifold_Pnt3D.coordinates[2] = manifold_components_total.coordinates[zAxis];
    }

    public void updatePCA_Pnt3D(int xAxis, int yAxis, int zAxis) {
        pca_Pnt3D.coordinates[0] = pca_components.coordinates[xAxis];
        pca_Pnt3D.coordinates[1] = pca_components.coordinates[yAxis];
        pca_Pnt3D.coordinates[2] = pca_components.coordinates[zAxis];
    }

    public void initManifold3D(int xAxis, int yAxis, int zAxis) {
        manifold_point3D.coordinates[0] = manifold_components_total.coordinates[xAxis];
        manifold_point3D.coordinates[1] = manifold_components_total.coordinates[yAxis];
        manifold_point3D.coordinates[2] = manifold_components_total.coordinates[zAxis];
        manifold_point3D.coordinates[3] = 1.0;
    }

    public void initPCA3D(int xAxis, int yAxis, int zAxis) {
        pca_point3D.coordinates[0] = pca_components.coordinates[xAxis];
        pca_point3D.coordinates[1] = pca_components.coordinates[yAxis];
        pca_point3D.coordinates[2] = pca_components.coordinates[zAxis];
        pca_point3D.coordinates[3] = 1.0;
    }

    public void initManifold_Pnt3D() {
        // Makes the assumption that manifold_point3D has been updated !

        int j;
        double dum;

        for (j = 0; j < 3; j++) {
            dum = manifold_point3D.coordinates[j];
            manifold_Pnt3D.coordinates[j] = (2.0 * (dum - myDomain.min_ManifoldPnt3D.coordinates[j]) / (myDomain.max_ManifoldPnt3D.coordinates[j] - myDomain.min_ManifoldPnt3D.coordinates[j])) - 1.0;
        }
    }

    public void putNeighbor(Integer u) {
        // if (name.equals(new String("658")))
        // System.out.println("Neighbor : " + ((Gene)
        // myDomain.domainGenes.get(myDomain.selectedGeneNumberToGeneNumber[u.intValue()])).name);
        double ac;

        if (neighbors == null) {
            neighbors = new Vector();
        }
        if (neighbor_angles == null) {
            neighbor_angles = new Vector();
        }
        if (!neighbors.contains(u)) {
            neighbors.add(u);
            ac = Math.acos(Gene.cosine(this, ((Gene) myDomain.getGenes().get(myDomain.selectedGeneNumberToGeneNumber[u.intValue()]))));
            neighbor_angles.add(new Double(ac));
            // System.out.println(ac + " --- " + Math.PI *
            // Statistics.LIMIT_P/2.0 + ", " + Math.PI * (1.0 -
            // (Statistics.LIMIT_P/2.0)));
        }

        if (neighbors.size() != neighbor_angles.size()) {
            Matrix.perror("Size mismatch");
        }
    }

    public boolean isEligible() {
    	if (AGCT.Method_F < 0 || AGCT.Method_F > 11) {
    		Matrix.perror("Feature method invalid");
    	}
        return geneExpressionContainsEmptyData(2);
    }

    public boolean geneExpressionContainsEmptyData(int n) {
        // returns true iff at least b values NON (?) for EACH selected LIGAND
        // if (ControlProcess.assertTrue("selectedLigands"))
        // Matrix.perror
        // ("You need to validate ligands before processing genes");
        //Debug.debug("Gene :: geneExpressionContainsEmptyData",ligands);
        for (int i = 0; i < ligands.size(); i++) {
            if (ligands.get(i).isChecked() && ligands.get(i).isEnabled()) {
                int nb = 0;
                for (int j = 0; j < getUndeterminedCoordinate(0).length; j++) {
                    if (getUndeterminedCoordinate(i)[j] == false) {
                        nb++;
                    }
                }
                if (nb < n) {
                    return false;
                }
            }
        }
        return true;
    }

    private HashMap<String, String> annotation = new HashMap<String, String>();

    public void setAnnotation(String prop, String value) {
        if (!ANNO_LIST.contains(prop)) {
            throw new IllegalArgumentException();
        }
        annotation.put(prop, value);
    }

    public String getAnnotation(String prop) {
        assert ANNO_LIST.contains(prop);
        //		assert annotation.containsKey(prop);
        return annotation.get(prop);
    }

    /**
     * 遺伝子の比較用の値.
     * save 時に、この値の照準で出力される
     */
    private double comparingValue;

    public void setComparingValue(double comparingValue) {
        this.comparingValue = comparingValue;
    }

    public double getComparingValue() {
        return comparingValue;
    }

    @Override
    public int compareTo(Gene o) {
        return (int) Math.signum(getComparingValue() - o.getComparingValue());
    }

    public boolean anyResponsive = false;

    public void setAnyResponsive(boolean b) {
        anyResponsive = b;
    }


    /**
     * なぜ消去されている？tk
     *
     * @param _genes
     */
    public static void showThis(Gene[] _genes) {
//        ArrayList<Gene> geneList = new ArrayList<Gene>();
//        for (Gene gene : _genes) {
//            if (gene.isEnabled) {
//                geneList.add(gene);
//            }
//        }
//        DoDebug.debug("showTime#geneList.size()", geneList.size());
//        Gene[] genes = geneList.toArray(new Gene[0]);
//
//        NormalizedGeneProfile normalizedGeneProfile = makeNormalizedGeneProfile(genes);
//        GeneClustering geneClustering = makeGeneClustering(genes);
//        try {
//            NormalizedGeneProfile.write(normalizedGeneProfile, "normalizedGeneProfile.txt");
//            GeneClustering.write(geneClustering, "geneClustering.txt");
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(0);
//        }
//        AllClustersGraphVisualizer.visualize(normalizedGeneProfile, geneClustering);
    }

    public int getClusterNumber() {
        int experimentId = AGCT.getInstance().getCurrentExperimentId();
        return getClusterNumber(experimentId);
    }
}
