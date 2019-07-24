package test;

import filtering.JLeftIndicator;
import forDebug.Debug;
import gene.Gene;
import gene.GeneList;
import gene.GeneListSaver;
import ligand.Ligand;
import ligand.LigandGroup;
import ligand.LigandList;
import parser.FileParseException;

import javax.swing.*;
import java.io.*;
import java.util.*;

public class Domain implements Debuggable {

    private static Domain singleton = new Domain();

    private Domain() {
    }

    public static Domain getInstance() {
        return singleton;
    }

    public final String domainName = "Domain";
    private final ArrayList<String> fileNames = new ArrayList<String>();

    public int getFileCount() {
        return fileNames.size();
    }

    public String getFileName() {
        String res = "";
        for (int i = 0; i < fileNames.size(); i++) {
            if (i > 0) {
                res += "_";
            }
            res += fileNames.get(i);
        }
        return res;
    }

    /**
     * ** LIGANDS ****
     */
    private LigandList _ligands;

    public LigandList getLigands() {
        return _ligands;
    }

    /**
     * ** TIMES ****
     */
    private double[] times;// times.lengthとは何を意味するのか？

    public double[] getTimes() {
        return times;
    }

    private int Replication;
    public static final String ANNOTATION_F = "F:";
    public static final String ANNOTATION_P = "P:";
    public static final String ANNOTATION_C = "C:";
    private String geneAnnotationsName = "";
    int maxHierarchy;
    /***** GENES *****/
    /**
     * ここで読み込まれたgeneのリスト
     */
    private GeneList genes;

    public GeneList getGenes() {
        return genes;
    }

    /**
     * ここで選ばれたgeneのリスト。 つまり、SDThresholdではじかれなかったもの
     */
    private GeneList selectedGenes;

    public GeneList getSelectedGenes() {
        return selectedGenes;
    }

    /***** FILTERING VARIABLES *****/
    // Following : Prototype defined variables (Gene name kept for simplicity)
    /**
     * 解析に使われるgeneの数
     */
    public int numberSelectedGenes;
    public int[] selectedGeneNumberToGeneNumber;
    // maps the index of the selected to the index of the corresponding gene in
    // domainGenes
    // private Hashtable accessToSelectedGeneNumberToGeneNumber;
    // does the same as selectedGeneNumberToGeneNumber, but with a Hashtable
    // A CHANGER : : avec Prototypes
    Vector<String> selectedGenesAllAnnotationsF;
    Vector<String> selectedGenesAllAnnotationsP;
    Vector<String> selectedGenesAllAnnotationsC;
    Hashtable<String, Integer> accessToSelectedGenesAllAnnotationsF;
    Hashtable<String, Integer> accessToSelectedGenesAllAnnotationsP;
    Hashtable<String, Integer> accessToSelectedGenesAllAnnotationsC;
    // Keys are annotations strings, returns the index as Integer for the
    // feature in the list of checkboxes of the JAnnotationFrame
    Hashtable<String, Vector<Integer>> annotationFToGenes;
    Hashtable<String, Vector<Integer>> annotationPToGenes;
    Hashtable<String, Vector<Integer>> annotationCToGenes;
    // Keys are annotations, returns an unordered list (Vector) of gene indexes
    // having the annotation
    boolean someAnnotationsF, someAnnotationsP, someAnnotationsC;
    public Pnt3D min_ManifoldPnt3D;
    Pnt3D min_PCAPnt3D;
    public Pnt3D max_ManifoldPnt3D;
    Pnt3D max_PCAPnt3D;
    /**
     * ** TRIANGULATION EDGES ****
     */
    // List of couples of genes
    DelaunayTriangulation myDelaunayTriangulation;
    Vector<Integer[]> myEdges;

    void readFile(File file) throws FileParseException, FileNotFoundException {
        new RawFileReader(file).readFileWithException();
    }

    static boolean READ = false;

    /**
     * 初期化。
     *
     * @param rf - raw data file
     */
    public void init(File rf) {
        READ = true;
        if (AGCT.MYDEBUG) {
            AGCT.debug("Domain # init", rf);
            AGCT.debug("stackTrace: ", Thread.currentThread().getStackTrace());
        }
        initMyValiables();
        try {
            readFile(rf);
        } catch (FileNotFoundException e) {
            JOptionPane.showMessageDialog(null, e);
            return;
        } catch (FileParseException ex) {
            StringBuilder sb = new StringBuilder();
            File file = ex.getFile();
            if (file == null) throw new AssertionError();
            String fileName = file.getName();
            int lineIndex = ex.getLineIndex();
            String line = ex.getLine();
            String message = ex.getMessage();
            sb.append("Parse error on line " + lineIndex + " of the file " + fileName + ".\n");
            sb.append("line:\n");
            sb.append(line + "\n");
            sb.append("reason:\n");
            sb.append(message);
            JOptionPane.showMessageDialog(null, sb.toString().replace('\t', ' '));
            return;
        }
        if (AGCT.getInstance().jNormalizationMenu.jProfileNormalizationMenu.kanoWay
                .isSelected()) {// かぶっている。動作しない設定なのか？？？
            // 条件分岐の草稿
            // if (judgeKanoNormalizable()){
            // FIXME デバッグ中, とりあえず元にもどしてみた。
            genes.kanoWay();// TODO ここ。time 0があるかどうかの条件分岐が行えていない。
        }
        String valret = "File opened: " + rf.getName() + " --- statistics:\n";
        Gene g;
        valret += "\nDomain " + domainName + ": \n";
        valret += "#Groups: ";
        valret += "\n";
        valret += "Groups composition:";
        valret += "#Ligands: " + _ligands.size()
                + " --- Ligand names [Groups: Times]: ";
        for (int i = 0; i < _ligands.size(); i++) {
            valret += _ligands.get(i).getName() + " ";
            valret += "[";
            try {
                valret += _ligands.get(i).getGroup();
            } catch (Exception e) {
                e.printStackTrace();
            }
            valret += ": ";
            for (int j = 0; j < times.length; j++) {
                valret += times[j];
                if (j < times.length - 1) {
                    valret += ", ";
                }
            }
            valret += "] ";
            if (i < _ligands.size() - 1) {
                valret += ", ";
            } else {
                valret += "\n";
            }
        }
        valret += "#Genes: " + getGenes().size() + " --- First genes: ";
        for (int i = 0; i < Math.min(getGenes().size(), 3); i++) {
            g = getGenes().get(i);
            valret += g.toString();
        }
        AGCT.getInstance().setInfoText(valret);

        if (!geneAnnotationsName.equals("")) {
            try {
                if (AGCT.LOAD_NEW_ANNOTATION_FILE) {
                    loadAnnotations_new(rf.getCanonicalPath().replaceAll(
                            rf.getName(), ""));
                } else {
                    loadAnnotations(rf.getCanonicalPath().replaceAll(
                            rf.getName(), ""));
                }
            } catch (IOException ee) {
                System.out.println("No canonical path for the file !");
            }
        } else {
            AGCT.getInstance().data.setAnnotationsExist(false);
            AGCT.getInstance().data.setAnnotationsExistP(false);
            AGCT.getInstance().data.setAnnotationsExistF(false);
            AGCT.getInstance().data.setAnnotationsExistC(false);
        }

        ControlProcess.put("annotationsLoaded", true);

        if (!AGCT.Referenced_Available) {
            AGCT.getInstance().highLightGene.setEnabled(false);
            AGCT.getInstance().onlyReferencedEdges.setEnabled(false);
        } else {
            AGCT.getInstance().highLightGene.setEnabled(true);
            AGCT.getInstance().onlyReferencedEdges.setEnabled(true);
        }
        AGCT.getInstance().loadHighlightFile.setEnabled(true);
        Statistics.flushAll();
    }

    private void initMyValiables() {
        maxHierarchy = -1;

        min_ManifoldPnt3D = new Pnt3D();
        max_ManifoldPnt3D = new Pnt3D();

        min_PCAPnt3D = new Pnt3D();
        max_PCAPnt3D = new Pnt3D();
    }

    public void removeDuplicateStringTables(Vector v) {
        int i = 0, j;
        String ti, tj;
        while (i < v.size() - 1) {
            j = i + 1;
            ti = ((String[]) v.elementAt(i))[1];
            do {
                tj = ((String[]) v.elementAt(j))[1];
                if (ti.equals(tj)) {
                    v.removeElementAt(j);
                } else {
                    j++;
                }
            } while (j < v.size());
            i++;
        }
    }

    public void loadAnnotations(String path) {
        if (AGCT.MYDEBUG) {
            AGCT.debug("Domain.loadAnnotations(" + path + ")");
        }
        FileReader e;
        StringTokenizer t;
        String s;
        BufferedReader br;
        String emptyString = "";
        String fullName = path + geneAnnotationsName;
        boolean collectData = false, foundFormat = false;
        Vector<String> variableTypes = null;
        int typeFormat = -1;

        AGCT.getInstance().setInfoText(
                "found annotations file : " + geneAnnotationsName
                        + " ... loading " + fullName);

        if (AGCT.DoDebug) {
            System.out.println("found annotations file : "
                    + geneAnnotationsName + " ... loading " + fullName);
        }

        try {
            e = new FileReader(new File(fullName));
        } catch (FileNotFoundException ex) {
            JInformationFrame.getInstance().setText(
                    "The annotation file you try to load does not exist");
            AGCT.getInstance().data.setAnnotationsExist(false);
            AGCT.getInstance().data.setAnnotationsExistP(false);
            AGCT.getInstance().data.setAnnotationsExistF(false);
            AGCT.getInstance().data.setAnnotationsExistC(false);
            return;
        }
        br = new BufferedReader(e);
        try {
            while ((s = br.readLine()) != null) {
                if (!collectData) {
                    t = new StringTokenizer(s, "\t");
                    if (t.countTokens() > 0) {
                        s = t.nextToken();

                        if (s.equals(AGCTFileWriter.ANNO_Data)) {
                            collectData = true;
                        } else {
                            if (s.equals(AGCTFileWriter.ANNO_Format)) {
                                s = t.nextToken();
                                if (s.equals(AGCTFileWriter.ANNO_Format_AGCT)) {
                                    typeFormat = 0;
                                    foundFormat = true;
                                } else if (s
                                        .equals(AGCTFileWriter.ANNO_Format_GEO)) {
                                    typeFormat = 1;
                                    foundFormat = true;
                                }
                            } else if (s.equals(AGCTFileWriter.ANNO_Empty)) {
                                emptyString = t.nextToken();
                            } else if (s
                                    .equals(AGCTFileWriter.ANNO_Variable_Definition)) {
                                variableTypes = new Vector<String>();
                                while (t.hasMoreTokens()) {
                                    s = t.nextToken();
                                    if ((!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Dum))
                                            && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_ID))
                                            && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Ascii_Names))
                                            && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_F))
                                            && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_P))
                                            && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_C))
                                            && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Referenced))) {
                                        Matrix.perror("Domain.class :: unknown variable tag while scanning file "
                                                + path);
                                    }
                                    variableTypes.addElement(s);
                                }
                            }
                        }
                    }
                } else {
                    if (!foundFormat) {
                        Matrix.perror("Domain.class :: no file format found inside Annotation file "
                                + path);
                    }

                    if (typeFormat == 0) {
                        loadGeneAnnotationAGCT_Format(s);
                    } else if (typeFormat == 1) {
                        if (variableTypes == null) {
                            Matrix.perror("Domain.class :: no variable format found inside Annotation file "
                                    + path);
                        }
                        if (emptyString.equals("")) {
                            Matrix.perror("Domain.class :: no 'Not_Defined' String defined inside Annotation file "
                                    + path);
                        }
                        loadGeneAnnotationGEO_Format(s, variableTypes,
                                emptyString);
                    } else {
                        Matrix.perror("Domain.class :: bad format value inside Annotation file "
                                + path);
                    }
                }
            }
            e.close();
        } catch (IOException ex) {
            JInformationFrame.getInstance().setText(
                    "IOError when loading annotation file");
            AGCT.getInstance().data.setAnnotationsExist(false);
            AGCT.getInstance().data.setAnnotationsExistP(false);
            AGCT.getInstance().data.setAnnotationsExistF(false);
            AGCT.getInstance().data.setAnnotationsExistC(false);
            return;
        }

        AGCT.getInstance().data.setAnnotationsExist(true);
        AGCT.getInstance().setInfoText(
                "found annotations file : " + geneAnnotationsName
                        + " ... loading... ok.");
    }

    private void loadGeneAnnotationAGCT_Format(String sRef) {
        StringTokenizer t = new StringTokenizer(sRef.replace('\t', ' '), " ");
        boolean readingF, readingP, readingC;

        if (t.countTokens() > 0) {
            String geneName = t.nextToken();
            if (!geneName.substring(0, AGCTFileWriter.DATA_Comments.length())
                    .equals(AGCTFileWriter.DATA_Comments)) {
                if (getGenes().contains(geneName)) {
                    readingF = readingP = readingC = false;
                    Gene gene = getGenes().get(geneName);
                    while (t.hasMoreTokens()) {
                        geneName = t.nextToken();
                        if (geneName.equals(Domain.ANNOTATION_F)) {
                            readingF = true;
                        } else if (geneName.equals(Domain.ANNOTATION_P)) {
                            readingP = true;
                        } else if (geneName.equals(Domain.ANNOTATION_C)) {
                            readingC = true;
                        } else {
                            if (!geneName.equals("")) {
                                if (readingF || readingP || readingC) {
                                    if ((readingF && readingP)
                                            || (readingF && readingC)
                                            || (readingP && readingC)) {
                                        System.out.print(" Gene "
                                                + gene.getName());
                                    }

                                    gene.addAnnotation(readingF, readingP,
                                            readingC, geneName);
                                    if (readingF) {
                                        AGCT.getInstance().data
                                                .setAnnotationsExistF(true);
                                    } else if (readingP) {
                                        AGCT.getInstance().data
                                                .setAnnotationsExistP(true);
                                    } else if (readingC) {
                                        AGCT.getInstance().data
                                                .setAnnotationsExistC(true);
                                    }

                                    readingF = readingP = readingC = false;
                                    gene.annotationsFilled = true;
                                } else {
                                    if (gene.asciiName == null) {
                                        gene.asciiName = geneName;
                                    }
                                }
                            }
                        }

                        if ((gene.annotationsF != null)
                                && (gene.annotationsF.size() > maxHierarchy)) {
                            maxHierarchy = gene.annotationsF.size();
                        }
                        if ((gene.annotationsP != null)
                                && (gene.annotationsP.size() > maxHierarchy)) {
                            maxHierarchy = gene.annotationsP.size();
                        }
                        if ((gene.annotationsC != null)
                                && (gene.annotationsC.size() > maxHierarchy)) {
                            maxHierarchy = gene.annotationsC.size();
                        }

                    }
                }
            }
        }
    }

    private void loadGeneAnnotationGEO_Format(String sRef, Vector<String> varTypes,
                                              String eString) {
        StringTokenizer t = new StringTokenizer(sRef, "\t");
        String geneID = "", gName = "", gF = "", gP = "", gC = "", gRef = "", dumS, refS;
        boolean readingF, readingP, readingC;
        Gene gg;
        int i, refe;
        String[] totalS1;
        String[] totalS2;

        if (t.countTokens() != varTypes.size()) {
            Matrix.perror("1Domain.class :: " + t.countTokens()
                    + " tokens, but " + varTypes.size()
                    + " variables on line \n" + sRef
                    + "\n(perhaps empty lines in .agct file ?)");
        }
        i = 0;
        while (t.hasMoreTokens()) {
            dumS = t.nextToken();
            refS = varTypes.elementAt(i);
            if (refS.equals(AGCTFileWriter.ANNO_Variable_Definition_ID)) {
                geneID = dumS;
            } else if (refS
                    .equals(AGCTFileWriter.ANNO_Variable_Definition_Ascii_Names)) {
                gName = dumS;
            } else if (refS
                    .equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_F)) {
                gF = dumS;
            } else if (refS
                    .equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_P)) {
                gP = dumS;
            } else if (refS
                    .equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_C)) {
                gC = dumS;
            } else if (refS
                    .equals(AGCTFileWriter.ANNO_Variable_Definition_Referenced)) {
                gRef = dumS;
            }
            i++;
        }

        if ((geneID.equals("")) || (geneID.equals(eString))) {
            Matrix.perror("Domain.class :: gene name is empty or dummy");
        }
        if (getGenes().contains(geneID)) {
            gg = getGenes().get(geneID);
            if (!gRef.equals("")) {
                refe = Integer.parseInt(gRef);
                if (refe > -1) {
                    gg.putReferenced();
                    gg.setTypeReferenced(refe);
                    if (refe > JAGCTGraphicsPanel.referencedColors.length) {
                        Matrix.perror("Domain.class :: referenced number "
                                + refe + " > "
                                + JAGCTGraphicsPanel.referencedColors.length);
                    }
                    AGCT.Referenced_Available = true;
                    if (refe > AGCT.Max_Reference) {
                        AGCT.Max_Reference = refe;
                    }
                }
            }

            if (!gg.annotationsFilled) {
                // F
                if (!gF.equals(eString)) {
                    totalS1 = gF.split(AGCTFileWriter.REPLACE_Old_New[0][0]);
                    for (i = 0; i < totalS1.length; i++) {
                        totalS2 = totalS1[i]
                                .split(AGCTFileWriter.REPLACE_Old_New[1][0]);
                        if (totalS2.length != 3) {
                            Matrix.perror("1Domain.class :: bad number of tokens for F in line\n"
                                    + sRef);
                        }
                        dumS = totalS2[1]; // take the central token
                        readingF = true;
                        readingP = readingC = false;
                        gg.addAnnotation(readingF, readingP, readingC,
                                dumS.replace(' ', '_'));

                        AGCT.getInstance().data.setAnnotationsExistF(true);
                    }
                }

                // P
                if (!gP.equals(eString)) {
                    totalS1 = gP.split(AGCTFileWriter.REPLACE_Old_New[0][0]);
                    for (i = 0; i < totalS1.length; i++) {
                        totalS2 = totalS1[i]
                                .split(AGCTFileWriter.REPLACE_Old_New[1][0]);
                        if (totalS2.length != 3) {
                            Matrix.perror("Domain.class :: bad number of tokens for P in line\n"
                                    + sRef);
                        }
                        dumS = totalS2[1]; // take the central token
                        readingP = true;
                        readingF = readingC = false;
                        gg.addAnnotation(readingF, readingP, readingC,
                                dumS.replace(' ', '_'));
                        AGCT.getInstance().data.setAnnotationsExistP(true);
                    }
                }

                // C
                if (!gC.equals(eString)) {
                    totalS1 = gC.split(AGCTFileWriter.REPLACE_Old_New[0][0]);
                    for (i = 0; i < totalS1.length; i++) {
                        totalS2 = totalS1[i]
                                .split(AGCTFileWriter.REPLACE_Old_New[1][0]);
                        if (totalS2.length != 3) {
                            Matrix.perror("Domain.class :: bad number of tokens for C in line\n"
                                    + sRef);
                        }
                        dumS = totalS2[1]; // take the central token
                        readingC = true;
                        readingF = readingP = false;
                        gg.addAnnotation(readingF, readingP, readingC,
                                dumS.replace(' ', '_'));
                        AGCT.getInstance().data.setAnnotationsExistC(true);
                    }
                }

                // Name
                if (!gName.equals(eString)) {
                    gg.asciiName = gName;
                }

                if (gg.annotationsF != null) {
                    removeDuplicateStringTables(gg.annotationsF);
                }
                if (gg.annotationsP != null) {
                    removeDuplicateStringTables(gg.annotationsP);
                }
                if (gg.annotationsC != null) {
                    removeDuplicateStringTables(gg.annotationsC);
                }

                if ((gg.annotationsF != null)
                        && (gg.annotationsF.size() > maxHierarchy)) {
                    maxHierarchy = gg.annotationsF.size();
                }
                if ((gg.annotationsP != null)
                        && (gg.annotationsP.size() > maxHierarchy)) {
                    maxHierarchy = gg.annotationsP.size();
                }
                if ((gg.annotationsC != null)
                        && (gg.annotationsC.size() > maxHierarchy)) {
                    maxHierarchy = gg.annotationsC.size();
                }

                gg.annotationsFilled = true;
            }
        }
    }

    private void loadAnnotations_new(String path) {
        if (AGCT.MYDEBUG) {
            AGCT.debug("Domain.loadAnnotations_new(" + path
                    + geneAnnotationsName + ")");
        }
        FileReader e;
        BufferedReader br;
        String fullName = path + geneAnnotationsName;
        AGCT.getInstance().setInfoText(
                "found annotations file : " + geneAnnotationsName
                        + " ... loading " + fullName);

        if (AGCT.MYDEBUG) {
            AGCT.debug("found annotations file : " + geneAnnotationsName
                    + " ... loading " + fullName);
        }

        try {
            e = new FileReader(new File(fullName));
        } catch (FileNotFoundException ex) {
            if (AGCT.MYDEBUG) {
                Debug.debug("The annotation file you try to load does not exist");
            }
            JInformationFrame.getInstance().setText(
                    "The annotation file you try to load does not exist");
            AGCT.getInstance().data.setAnnotationsExist(false);
            AGCT.getInstance().data.setAnnotationsExistP(false);
            AGCT.getInstance().data.setAnnotationsExistF(false);
            AGCT.getInstance().data.setAnnotationsExistC(false);
            return;
        }
        br = new BufferedReader(e);

        try {
            String s;
            while (br.readLine() != null) {
                while ((s = br.readLine()) != null) {
                    loadGeneAnnotationGPL_Format(s);
                }
            }
            e.close();
        } catch (IOException ex) {
            JInformationFrame.getInstance().setText(
                    "IOError when loading annotation file");
            AGCT.getInstance().data.setAnnotationsExist(false);
            AGCT.getInstance().data.setAnnotationsExistP(false);
            AGCT.getInstance().data.setAnnotationsExistF(false);
            AGCT.getInstance().data.setAnnotationsExistC(false);
            return;
        }
        AGCT.getInstance().data.setAnnotationsExist(true);
        AGCT.getInstance().setInfoText(
                "found annotations file : " + geneAnnotationsName
                        + " ... loading... ok.");
    }

    private String makeName(String name) {
        String[] ss = name.split("///");
        String res = null;
        int best = 1 << 28;
        for (String s : ss) {
            char[] cs = s.toCharArray();
            int t = 0;
            for (char c : cs) {
                if (Character.isDigit(c)) {
                    t++;
                }
            }
            if (res == null || t < best) {
                best = t;
                res = s;
            }
        }
        return res;
    }

    private void loadGeneAnnotationGPL_Format(String sRef) {
        StringTokenizer t = new StringTokenizer(sRef, "\t");
        String geneSymbol = null, gF = null, gP = null, gC = null, gRef = "", dumS, refS, geneTitle = null, refSeq = null;
        boolean readingF, readingP, readingC;
        int i, refE;
        String[] totalS1;
        String[] totalS2;

        if (t.countTokens() != Gene.ANNO_LIST.size()) {
            return;
        }
        i = 0;
        HashSet<String> gIds = new HashSet<String>();
        HashMap<String, String> map = new HashMap<String, String>();
        while (t.hasMoreTokens()) {
            dumS = t.nextToken();
            refS = Gene.ANNO_LIST.get(i);
            map.put(refS, dumS);
            if (refS.equals("ID") || refS.equals("GB_ACC")
                    || refS.equals("Representative Public ID")) {
                gIds.add(dumS);
            } else if (refS.equals("Target Description")) {// 20100824追加
                String[] ss = dumS.split("gb:");
                for (String s : ss) {
                    if (s.contains(" ")) {
                        String u = s.split(" ")[0];
                        if (u.contains(".")) {
                            String v = u.split("\\.")[0];
                            gIds.add(v);
                        } else {
                            gIds.add(u);
                        }
                    }
                }
            } else if (refS.equals("Gene Symbol")) {
                geneSymbol = makeName(dumS);
            } else if (refS.equals("Gene Ontology Molecular Function")) {
                gF = dumS;
            } else if (refS.equals("Gene Ontology Biological Process")) {
                gP = dumS;
            } else if (refS.equals("Gene Ontology Cellular Component")) {
                gC = dumS;
            } else if (refS.equals("Gene Title")) {
                geneTitle = dumS;
            } else if (refS.equals("RefSeq Transcript ID")) {
                refSeq = dumS;
            }
            i++;
        }

        if ((gIds.size() == 0)) {
            return;
        }
        for (String id : gIds) {
            if (getGenes().contains(id)) {
                Gene gene = getGenes().get(id);

                for (Map.Entry<String, String> entry : map.entrySet()) {
                    gene.setAnnotation(entry.getKey(), entry.getValue());
                }

                gene.geneTitle = geneTitle;
                gene.refSeq = refSeq;
                gene.setAnnotation(Gene.Gene$Symbol, geneSymbol);
                if (!gRef.equals("")) {
                    refE = Integer.parseInt(gRef);
                    if (refE > -1) {
                        gene.putReferenced();
                        gene.setTypeReferenced(refE);
                        if (refE > JAGCTGraphicsPanel.referencedColors.length) {
                            Matrix.perror("Domain.class :: referenced number "
                                    + refE
                                    + " > "
                                    + JAGCTGraphicsPanel.referencedColors.length);
                        }
                        AGCT.Referenced_Available = true;
                        if (refE > AGCT.Max_Reference) {
                            AGCT.Max_Reference = refE;
                        }
                    }
                }

                if (!gene.annotationsFilled) {
                    // F
                    if (!gF.equals(Not_defined)) {
                        totalS1 = gF
                                .split(AGCTFileWriter.REPLACE_Old_New[0][0]);
                        for (i = 0; i < totalS1.length; i++) {
                            totalS2 = totalS1[i]
                                    .split(AGCTFileWriter.REPLACE_Old_New[1][0]);
                            if (totalS2.length != 3) {
                                Matrix.perror("Domain.class :: bad number of tokens for F in line\n"
                                        + sRef);
                            }
                            dumS = totalS2[1]; // take the central token
                            readingF = true;
                            readingP = readingC = false;
                            gene.addAnnotation(readingF, readingP, readingC,
                                    dumS.replace(' ', '_'));
                            AGCT.getInstance().data.setAnnotationsExistF(true);
                        }
                    }

                    // P
                    if (!gP.equals(Not_defined)) {
                        totalS1 = gP
                                .split(AGCTFileWriter.REPLACE_Old_New[0][0]);
                        for (i = 0; i < totalS1.length; i++) {
                            totalS2 = totalS1[i]
                                    .split(AGCTFileWriter.REPLACE_Old_New[1][0]);
                            if (totalS2.length != 3) {
                                Matrix.perror("Domain.class :: bad number of tokens for P in line\n"
                                        + sRef);
                            }
                            dumS = totalS2[1]; // take the central token
                            readingP = true;
                            readingF = readingC = false;
                            gene.addAnnotation(readingF, readingP, readingC,
                                    dumS.replace(' ', '_'));
                            AGCT.getInstance().data.setAnnotationsExistP(true);
                        }
                    }

                    // C
                    if (!gC.equals(Not_defined)) {
                        totalS1 = gC
                                .split(AGCTFileWriter.REPLACE_Old_New[0][0]);
                        for (i = 0; i < totalS1.length; i++) {
                            totalS2 = totalS1[i]
                                    .split(AGCTFileWriter.REPLACE_Old_New[1][0]);
                            if (totalS2.length != 3) {
                                Matrix.perror("Domain.class :: bad number of tokens for C in line\n"
                                        + sRef);
                            }
                            dumS = totalS2[1]; // take the central token
                            readingC = true;
                            readingF = readingP = false;
                            gene.addAnnotation(readingF, readingP, readingC,
                                    dumS.replace(' ', '_'));
                            AGCT.getInstance().data.setAnnotationsExistC(true);
                        }
                    }

                    // Name
                    if (!geneSymbol.equals(Not_defined)) {
                        gene.asciiName = geneSymbol;
                    }

                    if (gene.annotationsF != null) {
                        removeDuplicateStringTables(gene.annotationsF);
                    }
                    if (gene.annotationsP != null) {
                        removeDuplicateStringTables(gene.annotationsP);
                    }
                    if (gene.annotationsC != null) {
                        removeDuplicateStringTables(gene.annotationsC);
                    }

                    if ((gene.annotationsF != null)
                            && (gene.annotationsF.size() > maxHierarchy)) {
                        maxHierarchy = gene.annotationsF.size();
                    }
                    if ((gene.annotationsP != null)
                            && (gene.annotationsP.size() > maxHierarchy)) {
                        maxHierarchy = gene.annotationsP.size();
                    }
                    if ((gene.annotationsC != null)
                            && (gene.annotationsC.size() > maxHierarchy)) {
                        maxHierarchy = gene.annotationsC.size();
                    }

                    gene.annotationsFilled = true;
                }
            }
        }
    }

    String Not_defined = "Not_defined";
    public ArrayList<String> orderedGroups;

    // GeneList highlights;
    public void loadHighlights(File rf) {
        if (AGCT.MYDEBUG) {
            AGCT.debug("Domain.loadHighlights");
        }
        boolean someHighlights = false;
        FileReader e;
        Gene gg;
        BufferedReader br;
        String s;
        StringTokenizer t;
        Integer il, ir, is;
        orderedGroups = null;

        // erase
        AGCT.Referenced_Available = false;
        AGCT.Max_Reference = -1;
        AGCT.getInstance().highLightGene.setEnabled(false);
        AGCT.getInstance().onlyReferencedEdges.setEnabled(false);
        for (int i = 0; i < getGenes().size(); i++) {
            if (getGenes().get(i).isSelected()) {
                gg = getGenes().get(i);
                gg.setReferenced(false);
                gg.setTypeReferenced(-1);
            }
        }
        try {
            e = new FileReader(rf);
        } catch (FileNotFoundException ex) {
            JInformationFrame.getInstance().setText(
                    "The file you try to load does not exist");
            return;
        }
        br = new BufferedReader(e);
        /*
      ** Highlighted genes ****
     */
        HashMap<String, Integer> nameToId;
        try {
            orderedGroups = new ArrayList<String>();
            nameToId = new HashMap<String, Integer>();
            while ((s = br.readLine()) != null) {
                if ((s.length() > 1)
                        && (!s.substring(0,
                        AGCTFileWriter.DATA_Comments.length()).equals(
                        AGCTFileWriter.DATA_Comments))) {
                    t = new StringTokenizer(s, ",");
                    if (t.countTokens() > 0) {
                        String name = t.nextToken();
                        String group = t.nextToken();
                        if (!orderedGroups.contains(group)) {
                            orderedGroups.add(group);
                        }
                        // TODO 文字列の場合に対応
                        // nameToId.put(name, orderedGroups.indexOf(group));
                        nameToId.put(name, Integer.valueOf(group));
                    }
                }
            }
            e.close();
        } catch (IOException ex) {
            ex.printStackTrace();
            JInformationFrame.getInstance().setText(
                    "IOError when loading domain(3)");
            return;
        }

        for (int i = 0; i < getGenes().size(); i++) {
            if (getGenes().get(i).isSelected()) {
                gg = getGenes().get(i);
                String name = gg.getName();
                String reference = gg.asciiName;

                il = nameToId.get(name);
                ir = null;
                if (reference != null) {
                    ir = nameToId.get(reference);
                }

                if ((il != null) && (ir != null)) {
                    Matrix.perror("Domain.class :: both il and ir are null");
                }
                if ((il != null) || (ir != null)) {
                    if (ir != null) {
                        is = ir;
                    } else {
                        is = il;
                    }
                    gg.setReferenced(true);
                    gg.setTypeReferenced(is);
                    if (AGCT.MYDEBUG) {
                        AGCT.debug(is);
                    }
                    someHighlights = true;
                }
            }
        }
        JLeftIndicator.getInstance().setHighlightIndicator(times);
        if (!AGCT.getInstance().highLightGene.isSelected()
                || AGCT.getInstance().highLightGene.isEnabled()) {
            AGCT.getInstance().highLightGene.setEnabled(true);
            AGCT.getInstance().highLightGene.doClick();
        }

        if (someHighlights) {
            AGCT.Referenced_Available = true;
            AGCT.Max_Reference = orderedGroups.size() - 1;
            AGCT.getInstance().highLightGene.setEnabled(true);
            AGCT.getInstance().onlyReferencedEdges.setEnabled(true);
        }
    }

    /**
     * 単に選ばれたgeneにAnnotationをつけているだけ。
     */
    public void fillSelectedGenesAllAnnotations() {
        Thread t = new Thread() {

            public void run() {

                selectedGenesAllAnnotationsF = new Vector<String>();
                selectedGenesAllAnnotationsP = new Vector<String>();
                selectedGenesAllAnnotationsC = new Vector<String>();

                annotationFToGenes = annotationPToGenes = annotationCToGenes = null;

                String s, annotation;
                AGCTCounter cc = new AGCTCounter("Filling Annotations [1/3]",
                        selectedGenes.size());
                for (Gene selectedGene : selectedGenes) {
                    if (selectedGene.annotationsF != null) {
                        if (annotationFToGenes == null) {
                            annotationFToGenes = new Hashtable<String, Vector<Integer>>();
                        }
                        for (int j = 0; j < selectedGene.annotationsF.size(); j++) {
                            s = ((String[]) selectedGene.annotationsF.elementAt(j))[0];
                            annotation = ((String[]) selectedGene.annotationsF
                                    .elementAt(j))[1];
                            if (!s.equals(Domain.ANNOTATION_F)) {
                                Matrix.perror("Domain.class :: annotation != F in F annotations");
                            }
                            if (!annotationFToGenes.containsKey(annotation)) {
                                selectedGenesAllAnnotationsF.add(annotation);
                                annotationFToGenes.put(annotation,
                                        new Vector<Integer>());
                            }
                            annotationFToGenes.get(annotation)
                                    .add(selectedGene
                                            .getIdNumber());
                        }
                    }
                    cc.increment();
                }
                cc.end();

                cc = new AGCTCounter("Filling Annotations [2/3]",
                        selectedGenes.size());
                for (Gene selectedGene : selectedGenes) {
                    if (selectedGene.annotationsP != null) {
                        if (annotationPToGenes == null) {
                            annotationPToGenes = new Hashtable<String, Vector<Integer>>();
                        }
                        for (int j = 0; j < selectedGene.annotationsP.size(); j++) {
                            s = ((String[]) selectedGene.annotationsP.elementAt(j))[0];
                            annotation = ((String[]) selectedGene.annotationsP
                                    .elementAt(j))[1];
                            if (!s.equals(Domain.ANNOTATION_P)) {
                                Matrix.perror("Domain.class :: annotation != P in P annotations");
                            }
                            if (!annotationPToGenes.containsKey(annotation)) {
                                selectedGenesAllAnnotationsP.add(annotation);
                                annotationPToGenes.put(annotation,
                                        new Vector<Integer>());
                            }
                            (annotationPToGenes.get(annotation))
                                    .add(selectedGene
                                            .getIdNumber());
                        }
                    }
                    cc.increment();
                }
                cc.end();

                cc = new AGCTCounter("Filling Annotations [3/3]",
                        selectedGenes.size());
                for (Gene selectedGene : selectedGenes) {
                    if (selectedGene.annotationsC != null) {
                        if (annotationCToGenes == null) {
                            annotationCToGenes = new Hashtable<String, Vector<Integer>>();
                        }
                        for (int j = 0; j < selectedGene.annotationsC.size(); j++) {
                            s = ((String[]) selectedGene.annotationsC.elementAt(j))[0];
                            annotation = ((String[]) selectedGene.annotationsC
                                    .elementAt(j))[1];
                            if (!s.equals(Domain.ANNOTATION_C)) {
                                Matrix.perror("Domain.class :: annotation != C in C annotations");
                            }
                            if (!annotationCToGenes.containsKey(annotation)) {
                                selectedGenesAllAnnotationsC.add(annotation);
                                annotationCToGenes.put(annotation,
                                        new Vector<Integer>());
                            }
                            (annotationCToGenes.get(annotation))
                                    .add(selectedGene
                                            .getIdNumber());
                        }
                    }
                    cc.increment();
                }
                cc.end();
                if (selectedGenesAllAnnotationsF.size() > 0) {
                    QuickSort.quicksortString(selectedGenesAllAnnotationsF);
                    someAnnotationsF = true;
                } else {
                    selectedGenesAllAnnotationsF = null;
                }
                if (selectedGenesAllAnnotationsP.size() > 0) {
                    QuickSort.quicksortString(selectedGenesAllAnnotationsP);
                    someAnnotationsP = true;
                } else {
                    selectedGenesAllAnnotationsP = null;
                }
                if (selectedGenesAllAnnotationsC.size() > 0) {
                    QuickSort.quicksortString(selectedGenesAllAnnotationsC);
                    someAnnotationsC = true;
                } else {
                    selectedGenesAllAnnotationsC = null;
                }

                // End: filling Hashtables

                if (someAnnotationsF) {
                    accessToSelectedGenesAllAnnotationsF = new Hashtable<String, Integer>();
                    for (int i = 0; i < selectedGenesAllAnnotationsF.size(); i++) {
                        accessToSelectedGenesAllAnnotationsF.put(
                                selectedGenesAllAnnotationsF
                                        .elementAt(i), i);
                    }
                }

                if (someAnnotationsP) {
                    accessToSelectedGenesAllAnnotationsP = new Hashtable<String, Integer>();
                    for (int i = 0; i < selectedGenesAllAnnotationsP.size(); i++) {
                        accessToSelectedGenesAllAnnotationsP.put(
                                selectedGenesAllAnnotationsP
                                        .elementAt(i), i);
                    }
                }

                if (someAnnotationsC) {
                    accessToSelectedGenesAllAnnotationsC = new Hashtable<String, Integer>();
                    for (int i = 0; i < selectedGenesAllAnnotationsC.size(); i++) {
                        accessToSelectedGenesAllAnnotationsC.put(
                                selectedGenesAllAnnotationsC
                                        .elementAt(i), i);
                    }
                }

                ControlProcess.put("selectedGenesAllAnnotationsOK", true);
            }
        };
        t.start();
    }

    public void initializationFilters() {
        for (Gene gene : getGenes()) {
            gene.setSelected(true);
            gene.setEnabled(true);
        }

        for (Ligand ligand : getLigands()) {
            ligand.setChecked(true);
            ligand.setEnabled(true);
        }

        for (LigandGroup group : getLigands().getGroups()) {
            group.setSelected(true);
        }

    }

    public void updateEnabledLigand() {
        int k;
        boolean found;
        boolean keep;

        for (int i = 0; i < _ligands.size(); i++) {
            LigandGroup group = _ligands.get(i).getGroup();
            keep = true;
            k = 0;
            found = false;
            do {
                LigandGroup gs = getLigands().getGroups()[k];
                if (gs.equals(group)) {
                    found = true;
                } else {
                    k++;
                }
            } while (!found && (k < getLigands().getGroups().length));
            if (!found) {
                Matrix.perror("Could not find ligand name " + group
                        + " in groups\n");
            }
            // if (!selectedGroup[0][k])
            if (!getLigands().getGroups()[k].isSelected()) {
                keep = false;
            }
            getLigands().get(i).setEnabled(keep);
        }
    }

    public void lastUpdateFilterData() {

        int i, index;
        Integer I;

        if (!Prototype.Prototypes_Selected) {
            Matrix.perror("Domain.class :: prototypes not selected");
        }

        numberSelectedGenes = 0;

        if (Prototype.No_Reduction) {
            for (Gene gene : getGenes()) {
                if (gene.isSelected()) {
                    numberSelectedGenes++;
                }
            }
            selectedGeneNumberToGeneNumber = new int[numberSelectedGenes];

            index = 0;
            for (i = 0; i < getGenes().size(); i++) {
                if (getGenes().get(i).isSelected()) {
                    selectedGeneNumberToGeneNumber[index] = i;
                    index++;
                }
            }
        } else {
            numberSelectedGenes = Prototype.Closest_Center_Ordered_List.size();
            selectedGeneNumberToGeneNumber = new int[numberSelectedGenes];

            for (index = 0; index < Prototype.Closest_Center_Ordered_List
                    .size(); index++) {
                I = (Integer) Prototype.Closest_Center_Ordered_List
                        .elementAt(index);
                selectedGeneNumberToGeneNumber[index] = I;
            }
        }

        ControlProcess.put("genesFilteredFor"
                + AGCT.Feature_Method[AGCT.Method_F], true);
    }

    /**
     * ligand , gene の checkedを更新、selectedGenesの作成aa kanno's Way などにより、profile
     * Normalizationも行っている。
     * <p/>
     * <p/>
     * <p/>
     * <p/>
     * 12/28/12 ここを調査しています。Tk こっちが本丸なのかな。。。
     */
    public void updateFilterData() {// used in JAGCTSelectionPane.finally.
        for (Ligand ligand : getLigands()) {
            ligand.setChecked((ligand.isChecked() && ligand.isEnabled()));
        }
        int numberSelectedLigands = getLigands().numberOfSelectedLigands();
        if (numberSelectedLigands == 0) {
            Matrix.perror("Too many constraints on ligands: 0 ligands selected");
        }
        for (Gene gene : getGenes()) {
            if (gene.isSelected() && gene.isEnabled()) {
                gene.setSelected(true);
            } else {
                gene.setSelected(false);
            }
        }

        // Vars to update later after prototype selection

        numberSelectedGenes = -1;
        selectedGeneNumberToGeneNumber = null;

        // TODO -> getSelectedGenes;

        // if (AGCT.MYDEBUG)
        // AGCT.debug(getGenes());

        selectedGenes = new GeneList();

        for (Gene gene : getGenes()) {
            if (gene.isSelected()) {
                selectedGenes.add(gene);
            }
        }

        if (AGCT.getInstance().jNormalizationMenu.jProfileNormalizationMenu.subtractMedian
                .isSelected()) {
            if (AGCT.MYDEBUG) {
                AGCT.debug("SUBT");
            }
            selectedGenes.subtractMedian();
        } else if (AGCT.getInstance().jNormalizationMenu.jProfileNormalizationMenu.kanoWay
                .isSelected()) {
            selectedGenes.kanoWay();
        }

        for (Gene gene : selectedGenes) {
            gene.finishVariables(AGCT.getInstance(), this);
            gene.recomputeSelectedFeatures();
            gene.calcTotalDistance();
        }
        String valret;
        valret = "#Ligands selected: " + numberSelectedLigands + "\n";
        int ns = 0;
        for (Ligand ligand : getLigands()) {
            if (ligand.isChecked()) {
                if (ns > 0) {
                    valret += ", ";
                }
                valret += ligand.getName();
                ns++;
            }
        }
        valret += "\n";
        valret += "#Genes selected before eventual prototype selection: "
                + selectedGenes.size() + "\n";
        for (int i = 0; i < selectedGenes.size(); i++) {
            if (i < AGCT.Max_Number_Of_Genes_Kept) {
                valret += selectedGenes.get(i).getName();
                if (i < selectedGenes.size() - 1) {
                    valret += ", ";
                }
            }
        }
        valret += "\n";
        AGCT.getInstance().setInfoText(valret);

        ControlProcess.put("filtersProcessed", true);

        // if (AGCT.MYDEBUG) {
        // AGCT.debug("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        // AGCT.debug(getGenes().size());
        // AGCT.debug(getLigands().size());
        // AGCT.debug(getSelectedGenes().size());
        // AGCT.debug(numberSelectedLigands);
        // }
    }

    public Gene getSelectedGene(int index) {
        // Integer I = new Integer(index);
        // return (Gene) accessToGenesObject.get((Integer)
        // accessToSelectedGeneNumberToGeneNumber.get(I));
        return selectedGenes.get(index);
    }

    /**
     * *************************************************************************************
     * Feature stuff
     * ***
     */
    private void computeGeneFeaturesWavelets() {
        if (AGCT.MYDEBUG) {
            AGCT.debug("Domain.computeGeneFeatureWavelets()");
        }

        ControlProcess.assertTrue("filtersProcessed");

        // Beware: stamps is used to compute haar basis, to save space

        int numberSelectedLigands = getLigands().numberOfSelectedLigands();
        double[][] stamps = new double[numberSelectedLigands][];
        double[] delta = new double[numberSelectedLigands];
        double[] begintime = new double[numberSelectedLigands];
        for (int i = 0, index = 0; i < _ligands.size(); i++) {
            if (getLigands().get(i).isChecked()) {
                if (AGCT.Number_Of_Wavelet_Stamps == -1) {
                    stamps[index] = Util_Feature.waveletsStampsAutomatic(times);
                } else {
                    stamps[index] = Util_Feature.waveletsStampsUserFixed(times);
                }
                // divide
                // to
                // number
                // of
                // AGCT.Number_Of_Wavelet_Stamps.

                delta[index] = stamps[index][1] - stamps[index][0];
                begintime[index] = stamps[index][0];
                index++;
            }
        }

        AGCTCounter cc = new AGCTCounter("Feature construction with Wavelets",
                selectedGenes.size());
        for (int i = 0; i < selectedGenes.size(); i++) {
            int index = 0;
            for (int j = 0; j < _ligands.size(); j++) {
                if (_ligands.get(j).isChecked()) {
                    double[] cvect = new double[stamps[index].length];
                    Util_Feature.fillTimeSerie(times, getGenes(), cvect,
                            begintime[index], delta, stamps[index].length,
                            // selectedGeneNumberToGeneNumberBefore[i]
                            selectedGenes.get(i).getIdNumber(), index, j);

                    if (AGCT.Method_F == 1) {
                        // ((Gene) genes.get(ig)).finalCoordinates[index] =
                        // Util_Feature.getHaar(cvect);
                        Debug.fdebug("tmp/Haar/cvect.txt", cvect);
                        selectedGenes.get(i).finalCoordinates[index] = Util_Feature
                                .getHaar(cvect);
                        Debug.fdebug("tmp/Haar/haar.txt",
                                selectedGenes.get(i).finalCoordinates[index]);
                    } else if ((AGCT.Method_F >= 2) && (AGCT.Method_F <= 10)) {
                        selectedGenes.get(i).finalCoordinates[index] = Util_Feature
                                .getDn(AGCT.Method_F - 2, cvect);
                    } else {
                        Matrix.perror("Unknown wavelet method");
                    }

                    index++;
                }
            }
            cc.increment();
        }
        cc.end();
    }

    /**
     * slopesに基づく
     */
    public void computeGeneFeaturesSlopes() {
        if (AGCT.MYDEBUG) {
            AGCT.debug("Domain.computeGeneFeatureSlopes()");
        }

        ControlProcess.assertTrue("filtersProcessed");
        // informationFrame = AGCT.getInstance().myInformationFrame;
        AGCTCounter cc = new AGCTCounter("Feature construction with Slopes",
                selectedGenes.size());
        for (int i = 0; i < selectedGenes.size(); i++) {
            // int ig = selectedGeneNumberToGeneNumberBefore[i];
            for (int j = 0, index = 0; j < _ligands.size(); j++) {
                if (getLigands().get(j).isChecked()) {
                    double tab = Util_Feature.slope(times, selectedGenes.get(i)
                            .getRawCoordinate(index), selectedGenes.get(i)
                            .getUndeterminedCoordinate(index));
                    selectedGenes.get(i).finalCoordinates[index] = new double[]{tab};
                    index++;
                }
            }
            cc.increment();
        }
        cc.end();
    }

    /**
     * Raw Coordinatesに基づく
     */
    public void computeGeneFeaturesRawCoordinates() {
        if (AGCT.MYDEBUG) {
            AGCT.debug("Domain.computeGeneFeaturesRawCoordinates()");
        }

        ControlProcess.assertTrue("filtersProcessed");
        
        for (int i = 0; i < selectedGenes.size(); i++) {
            for (int j = 0, index = 0; j < _ligands.size(); j++) {
                if (getLigands().get(j).isChecked()) {
                	// the activation values are actually the coordinates.
                    selectedGenes.get(i).finalCoordinates[index] = selectedGenes.get(i).getRawCoordinate(index).clone();
                    index++;
                }
            }
        }
    }

    public void finishComputeGeneFeatures() {
        // while (!Feature.featuresSelected);

        AGCT.getInstance().data.setDimFeatures(Util_Feature
                .getNumberOfFeaturesAfterSelection());

        AGCT.getInstance().myTabbedPane.mySelectionPane.finishValidate();
        ControlProcess.put("featuresSelected", true);
    }

    public void toTabbedPane() {
        AGCT.getInstance().myTabbedPane.toManifold();
        AGCT.getInstance().myTabbedPane.toPCA();
    }

    public void computeFinalCoordinates_1D() {
        int i;
        AGCT.getInstance().timerTickOn();
        AGCTCounter cc = new AGCTCounter("Computing Gene Coordinates",
                selectedGenes.size());
        for (i = 0; i < selectedGenes.size(); i++) {
            selectedGenes.get(i).computeFinalCoordinates_1D();
            // DoDebug.debug("GeneFoinalCoordinates_1D", i,
            // selectedGenes.get(i).finalCoordinates_1D);
            cc.increment();
        }
        cc.end();
        AGCT.getInstance().timerTickOff();
    }

    /**
     * AGCT.Method_F により、分岐処理を行っている。
     */
    private void featureConstruction() {
        AGCT.getInstance().timerTickOn();
        if (AGCT.Method_F == 0) {
            computeGeneFeaturesSlopes();
        } else if ((AGCT.Method_F >= 1) && (AGCT.Method_F <= 10)) {
            computeGeneFeaturesWavelets();
        } else if (AGCT.Method_F == 11) {
        	computeGeneFeaturesRawCoordinates();
        } else {
            Matrix.perror("Method selected to compute features does not exist !");
        }

        Gene.computeFinalIndex(selectedGenes.get(0).finalCoordinates);
        AGCT.getInstance().timerTickOff();
    }

    /**
     * AGCT.Method_FS により分岐処理を行っている。
     */
    public void featureSelection() {
        if ((AGCT.Method_FS > 0)
                && (AGCT.Max_Number_Of_Features < Util_Feature
                .getNumberOfFeaturesBeforeSelection(getLigands(), times))) {
        	System.err.println("--- Running feature selection ---\n");
            if (AGCT.Max_Number_Of_Features < 0) {
                if (AGCT.DoDebug) {
                    System.out.println("Fixing max number of features to "
                            + AGCT.DEFAULT_Max_Number_Of_Features);
                }
                AGCT.Max_Number_Of_Features = AGCT.DEFAULT_Max_Number_Of_Features;
            }
            Util_Feature.featureSelection(this);
        } else {
        	System.err.println("--- Not running feature selection ---\n");
            Util_Feature.noFeatureSelection(this);
        }
    }

    /**
     * DataForm -> Matrix に関わるいろいろなことをやっている。重要。
     */
    void featureAndPrototypeSelection() {
        System.out.println("Domain.featureAndPrototypeSelection()");
        final Domain me = this;
        Thread t = new Thread() {

            public void run() {
                if (AGCT.getInstance().myTabbedPane.mySelectionPane.keepPrototypes
                        .isSelected()) {
                    if (!getGenes().allNamesAreDistinct()) {
                        Matrix.perror("Domain.class :: some genes have the same name: Stopping.");
                    }
                }

                featureConstruction();

                featureSelection();
                computeFinalCoordinates_1D();

                if (Prototype.Loading_From_Scenario_Begin) {
                    while (!Prototype.Loading_From_Scenario_End)
                        ;
                } else {
                    Prototype.prototypeSelection(me);
                    ControlProcess.put("prototypeSelected", true);
                }

                lastUpdateFilterData();
                if ((someAnnotationsF) || (someAnnotationsP)) {
                    AGCT.getInstance().myAnnotationFrame = new JAnnotationFrame(
                            "AGCT --- Annotation Frame", AGCT.getInstance(),
                            AGCT.getInstance().data.getMyDomain());
                    AGCT.getInstance().myAnnotationFrame.setSize(
                            AGCT.WindowWidth, 600);
                    AGCT.getInstance().myFrame.toFront();
                    AGCT.getInstance().myAnnotationFrame
                            .changeLF_Prototype_Selection();
                    if (!Prototype.No_Reduction) {
                        AGCT.getInstance().myAnnotationFrame.onlyPrototypes
                                .setEnabled(true);
                    }
                }
                toTabbedPane();
                ControlProcess.put("featuresAndPrototypesSelected", true);
                int dim = Util_Feature.getNumberOfFeaturesAfterSelection();
                System.err.println("Got " + dim + " features after selection");
                
                for (int i = 0; i < selectedGenes.size(); i++) {
                	System.out.println("\nDumping gene " + selectedGenes.get(i).myID);
                    for (int j = 0, index = 0; j < _ligands.size(); j++) {
                        if (getLigands().get(j).isChecked()) {
                        	System.out.print("Ligand " + j + " has coordinates:");
                            for (double x : selectedGenes.get(i).finalCoordinates[index]) {
                            	System.out.print(" " + x);
                            }
                            index++;
                        	System.out.println("");
                        } else {
                        	System.out.println("Ligand " + j + " was disabled");
                        }
                        
                    }
                }
                int sign = JOptionPane.showConfirmDialog(null, "Do you have a priori genelist?", "priori genelist", JOptionPane.YES_NO_OPTION);

                if (sign == 0) {

                    new Thread() {

                        public synchronized void run() {
                            AGCT.getInstance().loadResponsiveFile.doClick();

                            Domain.getInstance().saveStat();

                            JLeftIndicator.getInstance().doit();

                            GeneList filtered = new GeneList();

                            for (Gene gene : Domain.getInstance().genes) {
                                if (gene.isVisible()) {
                                    filtered.add(gene);
                                    // DoDebug.debug("Domain : eliminated",gene);
                                }
                            }

                            gene.GeneListSaver saver = new GeneListSaver(
                                    filtered);
                            saver.setFile(new File(getFileName() + "rest.txt"));
                            saver.execute();
                        }
                    }.start();
                }
                // ここでvalidate終了
                if (AGCT.getInstance().AUTO) {

                    AGCT.getInstance().myTabbedPane.myManifoldPane.manifoldEigenButton
                            .doClick();
                } else {
                }
            }
        };
        t.start();
    }

    public void saveStat() {
        GeneList eliminated = new GeneList();

        for (Gene gene : Domain.getInstance().genes) {
            if (!Domain.getInstance().getSelectedGenes().contains(gene)) {
                eliminated.add(gene);
                // DoDebug.debug("Domain : eliminated",gene);
            }
        }

        gene.GeneListSaver saver = new GeneListSaver(eliminated);
        saver.setFile(new File(getFileName() + "_noise.txt"));
        saver.execute();

        GeneList eliminated_r = new GeneList();

        for (Gene gene : Domain.getInstance().genes) {
            if (gene.anyResponsive
                    && !Domain.getInstance().selectedGenes.contains(gene)) {
                eliminated_r.add(gene);
            }
        }

        saver = new GeneListSaver(eliminated_r);
        saver.setFile(new File(getFileName() + "_noise_r.txt"));
        saver.execute();
    }

    public void computeMinMax_ManifoldPnt3D() {
        System.out.println("Domain.computeMinMax_ManifoldPnt3D()");
        // Makes the assumption that gene components in manifold_point3D HAVE
        // BEEN computed OR correctly updated
        Gene gg;
        int i, j, id;
        for (i = 0; i < numberSelectedGenes; i++) {
            id = selectedGeneNumberToGeneNumber[i];
            gg = (Gene) getGenes().get(id);

            for (j = 0; j < 3; j++) {
                if ((i == 0)
                        || (gg.manifold_point3D.coordinates[j] < min_ManifoldPnt3D.coordinates[j])) {
                    min_ManifoldPnt3D.coordinates[j] = gg.manifold_point3D.coordinates[j];
                }
                if ((i == 0)
                        || (gg.manifold_point3D.coordinates[j] > max_ManifoldPnt3D.coordinates[j])) {
                    max_ManifoldPnt3D.coordinates[j] = gg.manifold_point3D.coordinates[j];
                }
            }
        }
        System.out.println("END_Domain.computeMinMax_ManifoldPnt3D()");
    }

    public void computeMinMax_PCAPnt3D() {
        // Makes the assumption that gene components in pca_point3D HAVE BEEN
        // computed OR correctly updated
        Gene gg;
        int i, j, id;
        for (i = 0; i < numberSelectedGenes; i++) {
            id = selectedGeneNumberToGeneNumber[i];
            gg = (Gene) getGenes().get(id);

            for (j = 0; j < 3; j++) {
                if ((i == 0)
                        || (gg.pca_point3D.coordinates[j] < min_PCAPnt3D.coordinates[j])) {
                    min_PCAPnt3D.coordinates[j] = gg.pca_point3D.coordinates[j];
                }
                if ((i == 0)
                        || (gg.pca_point3D.coordinates[j] > max_PCAPnt3D.coordinates[j])) {
                    max_PCAPnt3D.coordinates[j] = gg.pca_point3D.coordinates[j];
                }
            }
        }
    }

    public void fillPCAComponents(Matrix E, double[] ave, double[] sig) {
        AGCTCounter cc = new AGCTCounter("Filling pca comp.",
                numberSelectedGenes);
        int i, id;
        Gene gg;

        for (i = 0; i < numberSelectedGenes; i++) {
            id = selectedGeneNumberToGeneNumber[i];
            gg = (Gene) getGenes().get(id);

            gg.fillPCAComponents(E, ave, sig);
            cc.increment();
        }

        computeMinMax_PCAPnt3D();

        cc.end();
    }

    public void fillManifoldComponents(Matrix M) {
        System.out.println("Domain.fillManifoldComponents()");
        AGCTCounter cc = new AGCTCounter("Filling manifold components",
                numberSelectedGenes);

        int i, j, id;

        double dum;

        Gene gg;

        for (i = 0; i < numberSelectedGenes; i++) {
            cc.increment();

            id = selectedGeneNumberToGeneNumber[i];
            gg = (Gene) getGenes().get(id);

            gg.fillManifoldComponents(M, i);
        }

        computeMinMax_ManifoldPnt3D();

        cc.end();
        System.out.println("END_Domain.fillManifoldComponents()");
    }

    public void addIfNotInside(int[] couple) {
        int i;
        Integer[] c;
        boolean trouve = false;
        Gene gg;

        if (myEdges == null) {
            myEdges = new Vector<Integer[]>();
        } else {
            i = 0;
            do {
                c = myEdges.elementAt(i);
                if (((c[0].intValue() == couple[0]) && (c[1].intValue() == couple[1]))
                        || ((c[0].intValue() == couple[1]) && (c[1].intValue() == couple[0]))) {
                    trouve = true;
                } else {
                    i++;
                }
            } while ((trouve == false) && (i < myEdges.size()));
        }

        if (trouve == false) {
            c = new Integer[2];
            c[0] = new Integer(couple[0]);
            c[1] = new Integer(couple[1]);

            myEdges.addElement(c);

            gg = getGenes().get(
                    selectedGeneNumberToGeneNumber[c[0].intValue()]);
            gg.putNeighbor(c[1]);

            gg = getGenes().get(
                    selectedGeneNumberToGeneNumber[c[1].intValue()]);
            gg.putNeighbor(c[0]);
        }
    }

    public void toEdges() {
        System.out.println("toEdges :: return");
        if (myDelaunayTriangulation != null) {
            AGCTCounter cc = new AGCTCounter("Creating edges",
                    myDelaunayTriangulation.getNeighbors().size());

            JInformationFrame.getInstance()
                    .setTextProgressBar("Creating edges");
            Simplex s;
            Iterator it;
            int i, ip1, ip2, iu1, iu2, nmade = 0, percent;
            int[] couple = new int[2];
            Vector v;

            for (it = myDelaunayTriangulation.getNeighbors().keySet()
                    .iterator(); it.hasNext(); ) {
                cc.increment();

                s = (Simplex) it.next();
                v = new Vector(s.getVertices());
                for (i = 0; i < v.size(); i++) {
                    ip1 = i;
                    ip2 = i + 1;
                    if (i == v.size() - 1) {
                        ip2 = 0;
                    }

                    iu1 = ((Pnt) v.elementAt(ip1)).getGeneIndex();
                    iu2 = ((Pnt) v.elementAt(ip2)).getGeneIndex();

                    if ((iu1 >= 0) && (iu2 >= 0)) {
                        couple[0] = iu1;
                        couple[1] = iu2;
                        addIfNotInside(couple);
                    }
                    if (myEdges != null) {
                        JInformationFrame.getInstance().setTextProgressBar(
                                "Creating edges :: " + myEdges.size()
                                        + " current edges.");
                    }
                }
            }
            cc.end();
            AGCT.getInstance().myTabbedPane.myManifoldPane.visualizationPanel
                    .repaint();
        }
    }

    private class RawFileReader {

        private File rawFile, file;
        int currentLineNumber = -1;
        String currentLine = "";
        Scanner scanner;

        //        BufferedReader reader;
        String readLine() {
            if (!scanner.hasNextLine()) {
                if (currentLine != null) {
                    currentLineNumber++;
                    currentLine = null;
                }
                return null;
            }
            currentLineNumber++;
            currentLine = scanner.nextLine();
            return currentLine;
        }

        public RawFileReader(File rawFile) {
            this.rawFile = rawFile;
        }


        void readDirectory(File dir) throws FileParseException, FileNotFoundException {
            File[] files = dir.listFiles();
            if (files != null) for (File file : files) readFile(file);
        }

        void readFileWithException() throws FileParseException, FileNotFoundException {
            readFile(rawFile);
        }

        /**
         * ReadFile ここで変な変更があった可能性がある。このとつうからGarudaが合流するよ。オーバーライドして
         * ほとんど同じだが別のGaruda用ｗ作ってもいいかもしれない。
         *
         * @param file
         * @return true if the file is successfully loaded.
         */
        private void readFile(File file) throws FileParseException, FileNotFoundException {
            if (AGCT.MYDEBUG) {
                AGCT.debug("Domain.readNewFile(file,agct)");
            }
            if (file.isDirectory()) {
                readDirectory(file);
                return;
            }
            this.file = file;
            currentLineNumber = 0;
            if (file.getName().startsWith(".")) {
                return;
            }
            try {
                scanner = new Scanner(file);
            } catch (FileNotFoundException ex) {
                JInformationFrame.getInstance().setText(
                        "The file you try to load does not exist");
                AGCT.getInstance().data.setRawDomainExists(false);
                throw new FileNotFoundException("The file " + file.getAbsolutePath() + " does not exist.");
            }

            if (AGCT.MYDEBUG) {
                AGCT.debug(file.getName());
            }
            fileNames.add(file.getName().split("\\.")[0]);

            // If multiple files are loaded, ligands in the current file is appended to the ligands that was loaded before.
            // To specify where the beginning of the current file is, current ligand size is needed.
            final int firstLigandSize = getLigands() == null ? 0 : getLigands().size();
            try {
                String line;
                while ((line = readLine()) != null) {
                    if (AGCT.MYDEBUG) {
                        AGCT.debug(line);
                    }
                    parseLine(firstLigandSize, line);
                }
                Debug.debug("Species = " + Domain.this.geneAnnotationsName);
                Debug.debug("#entities = " + _ligands.size());

                scanner.close();
            } catch (FileParseException ex) {
                ex.printStackTrace();
                JInformationFrame.getInstance().setText(
                        "IOError when loading domain");
                AGCT.getInstance().data.setRawDomainExists(false);
                throw ex;
            }

            if (_ligands == null || _ligands.size() == 0) {
                throw FileParseException.make(file, "", currentLineNumber,
                        "No ligand counted. Perhaps a coding problem for your file ?\n"
                                + "Annotation and input file(s) (or directory) must be same folder.");
            }
        }

        private void parseLine(int firstLigandSize, String line) throws FileParseException {
            StringTokenizer tokenizer = new StringTokenizer(line, "\t");
            if (tokenizer.countTokens() > 0) {
                String firstToken = tokenizer.nextToken();
                if ((firstToken.length() >= 1)
                        && (!firstToken.startsWith(AGCTFileWriter.DATA_Comments))) {
                    if (firstToken.startsWith("EG")) { // entity group
                        if (AGCT.MYDEBUG) {
                            AGCT.debug("EG");
                        }
                        readGroups(tokenizer, fileNames.get(fileNames.size() - 1),
                                firstLigandSize);
                    } else if (firstToken.startsWith("E")) { // entity
                        if (AGCT.MYDEBUG) {
                            AGCT.debug("E");
                        }
                        readEntities(tokenizer,
                                fileNames.get(fileNames.size() - 1));
                    } else if (firstToken.startsWith("T")) { // time
                        if (AGCT.MYDEBUG) {
                            AGCT.debug("T");
                        }
                        readTimes_new(tokenizer);
                    } else if (!firstToken.startsWith("Responsive")
                            && firstToken.startsWith("R")) { // Replicate
                        if (AGCT.MYDEBUG) {
                            AGCT.debug("R");
                        }
                        readReplication(tokenizer);
                    } else if (firstToken.startsWith("S")) {  // Species
                        if (AGCT.MYDEBUG) {
                            AGCT.debug("S");
                        }
                        // If no species is specified, it means there is no annotation.
                        if (tokenizer.hasMoreTokens())
                            geneAnnotationsName = tokenizer.nextToken();
                    } else if (firstToken.startsWith("D")) {     // Data
                        if (AGCT.MYDEBUG) {
                            AGCT.debug("D");
                        }
                        readGenes(firstLigandSize);
                    } else {
                        throw parseError("Unknown keyword `" + firstToken +
                                "' was specified.");

                    }
                }
            }
        }

        /**
         * @param st
         * @param fileName
         * @param firstLigandSize
         */
        private void readGroups(StringTokenizer st, String fileName,
                                int firstLigandSize) throws FileParseException {
            if (_ligands == null) {
                System.err.println("put entity BEFORE group!");
                throw new RuntimeException("put entity BEFORE group!");
            }
            for (int id = firstLigandSize; id < _ligands.size(); id++) {
                if (!st.hasMoreTokens()) {
                    throw parseError(_ligands.size() + " entity groups is expected, but there were less tokens.");
                }
                String groupName = fileName + "." + st.nextToken();
                _ligands.get(id).setGroup(LigandGroup.get(groupName));
            }
        }

        private void readReplication(StringTokenizer st) throws FileParseException {
            Replication = readInt(st);
        }

        /**
         * @param st
         */
        private void readTimes_new(StringTokenizer st) throws FileParseException {
        	/* TODO: If the user only wants to use raw coordinates, we do not need time stamps. We abuse
        	 * the number of time stamps (numTime) to mean the length of the raw coordinate vector.
        	 * the actual t_i are irrelevant, so maybe we should not throw an error if they are missing.
        	 */
            if (!st.hasMoreTokens()) throw parseError("Number of time stamps is not specified.");
            String s = st.nextToken();
            int numTime;
            try {
                numTime = Integer.parseInt(s);
            } catch (NumberFormatException e) {
                throw parseError("Number of time stamp is expected being an integer, but was " + s + ".");
            }
            times = new double[numTime];
            for (int i = 0; i < times.length; i++) {
                if (!st.hasMoreTokens()) {
                    throw FileParseException.make(file, currentLine, currentLineNumber,
                            "Number of time stamp is specified as " + numTime + ", but there are less tokens.");
                }
                s = st.nextToken();
                try {
                    times[i] = Double.valueOf(s);
                } catch (NumberFormatException e) {
                    throw parseError("Time stamp is expected being a number, but was " + s + ".");
                }
            }
        }

        private FileParseException parseError(String message) {
            return FileParseException.make(file, currentLine, currentLineNumber, message);
        }

        /**
         * @param st
         * @param fileName
         */
        private void readEntities(StringTokenizer st, String fileName) throws FileParseException {
            int numEntities = readInt(st);
            if (_ligands == null) {
                _ligands = new LigandList();
            }
            for (int i = 0; i < numEntities; i++) {
                if (!st.hasMoreTokens()) {
                    throw parseError("Specified number of entities is " + numEntities + ", but there are less tokens.");
                }
                String name = fileName + "." + st.nextToken();
                Ligand ligand = new Ligand();
                ligand.setName(name);
                ligand.setFileName(fileName);
                _ligands.add(ligand);
            }
        }

        private int readInt(StringTokenizer st) throws FileParseException {
            if (!st.hasMoreTokens())
                throw parseError("An integer was expected, but there was no token.");
            String s = st.nextToken();
            try {
                return Integer.parseInt(s);
            } catch (NumberFormatException ignored) {
                throw parseError("An integer was expected, but was " + s);
            }
        }

        /**
         * read Gene Expression data
         *
         * @param firstLigandSize
         */
        private void readGenes(int firstLigandSize) throws FileParseException {
            Gene.setMyDomain(Domain.this);
            if (genes == null) {
                genes = new GeneList();
                JLeftIndicator.getInstance().setGeneList(genes);
                Gene.setGeneList(genes);
            }
            PrintStream tmp = System.out;
            int idNumber = 0;
            String line;
            while ((line = readLine()) != null) {// while EOF is not reached.
                StringTokenizer st = new StringTokenizer(
                        line.replace('\t', ' '), " ");
                if (st.countTokens() > 0) {
                    String geneId = st.nextToken();
                    Gene gene = genes.createAndAddGene(geneId);
                    gene.setLigands(getLigands());
                    gene.setIdNumber(idNumber++);
                    boolean alreadyInvisible = false;
//                    gene.setVisible(true);
//                    gene.setVisibleForOr(true);
                    for (int i = firstLigandSize; i < _ligands.size(); i++) {
                        double[][] input = new double[times.length][Replication];
                        double[] add = new double[times.length];
                        boolean[] addB = new boolean[times.length];
                        for (int j = 0; j < times.length; j++) {
                            double[] activations = new double[Replication];
                            for (int k = 0; k < activations.length; k++) {
                                input[j][k] = activations[k] = readDouble(st);
                            }
                            add[j] = Util_Math.calcMean(activations, AGCT
                                    .getInstance().getNormalizationData()
                                    .getNormalizeMode());
                            if (Util_Math.calcStandardDivation(activations) >= AGCT
                                    .getInstance().getNormalizationData()
                                    .getSDThreshold()) {
                                addB[j] = true;
                                if(!alreadyInvisible){
//                                    gene.setVisible(false);
//                                    gene.setVisibleForOr(false);
//                                    alreadyInvisible = true;
                                }
                            } else {
                                addB[j] = false;
                            }
                        }
                        gene.addInput(input);
                        gene.addRawCoordinate(add);
                        gene.addInitialCoordinate(add.clone());
                        gene.addUndeterminedCoordinate(addB);
                    }
                    gene.makeFeatures();
                    if (st.hasMoreTokens()) {
                        if (AGCT.MYDEBUG) {
                            AGCT.debug("hasMoreTokensERROR!!");
                        }
                        System.exit(0);
                    }
                    gene.setMyAGCT(AGCT.getInstance());
                }
            }
            System.setOut(tmp);
        }

        private double readDouble(StringTokenizer st) throws FileParseException {
            if (!st.hasMoreTokens()) {
                throw parseError("A real number was expected, but there was no token.");
            }
            String s = st.nextToken();
            try {
                return Double.parseDouble(s);
            } catch (NumberFormatException ignore) {
                throw parseError("A real number was expected, but was " + s + ".");
            }
        }
    }
}
