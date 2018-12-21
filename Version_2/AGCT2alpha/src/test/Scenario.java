package test;

import forDebug.Debug;

import java.io.*;
import java.util.StringTokenizer;
import java.util.Vector;

public class Scenario {
    public static String allKeywords[] = {
            "Load",                         // 0
            "JAGCTSelectionPane_Method_F",  // 1
            "JAGCTSelectionPane_Validate",  // 2
            "JAGCTCheckBox_StateChanged",   // 3
            "AGCT_Modification_Method_F",   // 4
            "AGCT_Modification_Method_W",   // 5
            "AGCT_Modification_Method_S",   // 6
            "AGCT_Modification_Method_N",   // 7
            "AGCT_Modification_Number_Of_Neighbors", // 8
            "AGCT_Modification_Number_Of_Manifold_Components", // 9
            "AGCT_Modification_Number_Of_Triangulation_Dimensions", // 10
            "AGCT_Modification_T_Heat_Kernel",                     // 11
            "AGCT_Modification_T_Heat_Kernel_Similarity_Neighbors",  // 12
            "AGCT_Modification_Sparsify_Statistic",                  // 13
            "JDSScaling_Modification_Bregman_Divergence",            // 14
            "Statistics_Modification_Limit_P", // 15
            "AGCT_Modification_Use_Debug",
            "AGCT_Modification_Warning",
            "AGCT_Modification_Perspective",
            "AGCT_Modification_Sort_Depth",
            "AGCT_Modification_Use_Shadow", // 20
            "AGCT_Manifold_Processed",
            "AGCT_PCA_Processed",
            "AGCT_Manifold_Triangulated",
            "JAGCTClusteringPane_Run_Clustering",
            "AGCT_Modification_Number_Of_Wavelet_Stamps",// 25
            "JAGCTSelectionPane_Max_Number_Of_Features",
            "JAGCTSelectionPane_Feature_Selection_Method",
            "JAGCTSelectionPane_Max_Number_Of_Prototypes",
            "JAGCTSelectionPane_Prototype_Selection_Method",
            "JAGCTSelectionPane_Keep_Prototypes", // 30
            "JAGCTSelectionPane_Prototypes_Begin",
            "JAGCTSelectionPane_Prototypes_End",
            "JAGCTManifoldPane_Keep_Manifold",
            "JAGCTManifoldPane_Manifold_Begin",
            "JAGCTManifoldPane_Manifold_End",// 35
            "JAGCTClusteringPane_Keep_Clustering",
            "JAGCTClusteringPane_Clustering_Begin",
            "JAGCTClusteringPane_Clustering_End",
            "JAGCTClickLigandCheckBox",
            "AGCT_Modification_SD_Threshold"// 40
    };

    private Scenario() {
    }

    /**
     * ****************************************************************************************
     * Class that records and plays scenarios
     * ***
     */

    // Manifold = ORDERED list of gene names, W (last step), M (ready for
    // components)
    public static boolean runningForPrototypes, bClosest_Center_Ordered_List, bClosest_Center_To_Cluster_Number, bCluster_Number_To_Closest_Center,
            bClosest_Center_To_Cluster_Points, bClosest_Center_To_Normalized_Distortions;

    /**
     * manifoldを読んでいる最中である
     */
    public static boolean runningForManifold;
    public static boolean bOrdered_List_Names;
    public static boolean bMatrix_W;
    public static boolean bManifold_Eigenvalues;
    public static boolean bMatrix_M;

    public static String myToken = "|";

    public static String commentString = "//";

    // Load | datafile = loads a datafile
    // JAGCTSelectionPane_Method_F | Integer = codes the Method_F variable
    // JAGCTSelectionPane_Validate | "" = validates on JAGCTSelectionPane
    // JAGCTCheckBox_StateChanged | Integer Type = gives a checkbox whose
    // checked status has changed
    public static Vector<String> allStrings = new Vector<String>();

    public static boolean scenarioTick = false, isRunning = false;

    public static void block(int i) {
        System.out.println("Scenario.block(" + i + ")");
        if (((i >= 1) && (i <= 3)) || ((i >= 21) && (i <= 29))) {
            while (!ControlProcess.hasTrue("displaySelectionPane")) ;
            // while(!ControlProcess.hasTrue("dataLoaded"));
        }
        if (i == 21)
            while (!ControlProcess.hasTrue("featuresSelected")) ;
        if (i == 21)
            while (!ControlProcess.hasTrue("prototypeSelected")) ;
        if (i == 21)
            while (!ControlProcess.hasTrue("featuresAndPrototypesSelected")) ;

        if ((i == 34) || (i == 35))
            while (!ControlProcess.hasTrue("featuresSelected")) ;
        if ((i == 34) || (i == 35))
            while (!ControlProcess.hasTrue("prototypeSelected")) ;
        if ((i == 34) || (i == 35))
            while (!ControlProcess.hasTrue("featuresAndPrototypesSelected")) ;
        if ((i == 34) || (i == 35))
            while (!ControlProcess.hasTrue("annotationsLoaded")) ;

        // if (i==31)
        // while(!ControlProcess.hasTrue("selectedGenesAllAnnotations"));
        // while(!ControlProcess.hasTrue("annotationsLoaded"));

        if ((i == 22) || (i == 23))
            while (!ControlProcess.hasTrue("manifoldProcessed")) ;
        if (i == 24)
            while (!ControlProcess.hasTrue("clusteringAvailable")) ;
    }

    public static void flushAll() {
        runningForPrototypes = runningForManifold = false;
        allStrings = null;
        scenarioTick = false;
        isRunning = false;
    }

    public static void tickOn() {
        scenarioTick = true;
    }

    public static void tickOff() {
        scenarioTick = false;
    }

    public static void tickSwitch() {
        if (AGCT.MYDEBUG)
            AGCT.debug("Scenario.tickSwitch");
        scenarioTick = !scenarioTick;
    }

    public static boolean isTicking() {
        return scenarioTick;
    }

    public static boolean isRunning() {
        return isRunning;
    }

    public static void newScenario() {
        allStrings = new Vector();
    }

    /**
     * 出力に t+\n を追加する
     */
    public static void addThis(String t) {
        Debug.debug("addThis", t);
        //		if(t.contains("ClusteringPane") || t.contains("Manifold")){
        //			try{
        //				throw new RuntimeException();
        //			}catch (Exception e) {
        //				e.printStackTrace();
        //			}
        //		}
        if (allStrings == null)
            newScenario();
        allStrings.addElement(t);
    }

    public static boolean containsKwd(String t) {
        if (t == null)
            return false;

        int i;
        for (i = 0; i < allKeywords.length; i++)
            if (allKeywords[i].equals(t))
                return true;
        return false;
    }

    public static int kwdIndex(String t) {
        if (t == null)
            return -1;

        int i;
        for (i = 0; i < allKeywords.length; i++)
            if (allKeywords[i].equals(t))
                return i;
        return -1;
    }

    public static void add(String kwd) {
        add(kwd, "");
    }

    public static void add(String kwd, int i) {
        add(kwd, i + "");
    }

    public static void add(String kwd, double i) {
        add(kwd, i + "");
    }

    public static void add(String kwd, boolean i) {
        add(kwd, boolean2int(i) + "");
    }

    public static void add(String kwd, String arg) {
        if (scenarioTick) {
            if (!containsKwd(kwd))
                Matrix.perror("Scenario.class :: Keyword " + kwd + " not registered");
            addThis(kwd + myToken + arg);
        }
    }

    public static boolean isComment(String s) {
        if (s.length() < Scenario.commentString.length())
            return false;
        String t = s.substring(0, 2);
        if (t.equals(Scenario.commentString))
            return true;
        return false;
    }

    public static boolean isKeepPrototypes(String s) {
        int ref = allKeywords[30].length();

        if (s.length() < ref)
            return false;
        String t = s.substring(0, ref);
        if (t.equals(allKeywords[30]))
            return true;
        return false;
    }

    public static boolean isKeepManifold(String s) {
        int ref = allKeywords[33].length();

        if (s.length() < ref)
            return false;
        String t = s.substring(0, ref);
        if (t.equals(allKeywords[33]))
            return true;
        return false;
    }

    public static boolean isKeepClustering(String s) {
        int ref = allKeywords[36].length();
        if (s.length() < ref)
            return false;
        String t = s.substring(0, ref);
        if (t.equals(allKeywords[36]))
            return true;
        return false;
    }

    public static void loadExecute(File file, AGCT agct) {
        if (AGCT.MYDEBUG)
            AGCT.debug("Scenario.loadExecute(file,agct)");
        FileReader e;
        BufferedReader br;
        String s;

        try {
            e = new FileReader(file);
        } catch (FileNotFoundException ex) {
            JInformationFrame.getInstance().setText("The scenario you try to load does not exist");
            return;
        }

        newScenario();
        br = new BufferedReader(e);
        try {
            while ((s = br.readLine()) != null) {
                if (!isComment(s))
                    addThis(s);
            }
            e.close();
        } catch (IOException ex) {
            JInformationFrame.getInstance().setText("IOError when loading scenario");
            agct.data.setRawDomainExists(false);
            return;
        }
        e = null;

        if ((allStrings != null) && (allStrings.size() > 0))
            executeScenario(agct);
    }

    /**
     * Scenarioを実行する
     * この前の段階で、allStringsにファイルの、コメント以外のすべての行が入っている
     *
     * @param ap
     */
    private static void executeScenario(final AGCT ap) {
        AGCT.AUTO = false;
        if (AGCT.MYDEBUG)
            AGCT.debug("Scenario.executeScenario");
        Thread t = new Thread() {
            public void run() {
                int i;
                isRunning = true;
                System.out.print("Executing scenario (" + allStrings.size() + " lines)... please wait... ");

                ap.myTabbedPane.setAllEnabled(false);

                AGCTCounter cc = new AGCTCounter("Scenario...", allStrings.size());
                for (i = 0; i < allStrings.size(); i++) {
                    executeLine(allStrings.get(i), ap);
                    cc.setText("Scenario...");
                    cc.increment();
                }
                cc.end();

                flushAll(); // isRunning = false; (saves memory)

                ap.myTabbedPane.setAllEnabled(true);

                ap.scenarioLoadButton.setIcon(ap.scenarioLoadIdleIcon);
            }
        };
        t.start();
    }

    public static int index(String s) {
        int i = 0;
        while (i < allKeywords.length) {
            if (s.equals(allKeywords[i]))
                return i;
            else
                i++;
        }
        return -1;
    }

    public static int indexFirstToken(String ex) {
        StringTokenizer t = new StringTokenizer(ex, Scenario.myToken);
        if (!t.hasMoreTokens())
            return -1;
        return index(t.nextToken());
    }

    public static boolean int2boolean(int dum) {
        if (dum == 0)
            return false;
        else if (dum == 1)
            return true;
        else
            Matrix.perror("Scenario :: not a boolean");
        return false;
    }

    public static int boolean2int(boolean dum) {
        if (dum)
            return 1;
        return 0;
    }

    public static void displayInformationFrame(AGCT ap) {
        String s = "";

        if ((allStrings == null) || (allStrings.size() == 0))
            s = "No scenario to display";
        else {
            s = "Currently recorded scenario:\n";
            int i;
            for (i = 0; i < allStrings.size(); i++) {
                s += allStrings.elementAt(i) + "\n";
            }
        }
        JInformationFrame.getInstance().setText(s);
    }

    private static void executeLine(String ex, AGCT ap) {
        //		System.out.println("Scenario.executeLine() "+ex);
        if (runningForPrototypes)
            executeLinePrototypes(ex, ap);
        else if ((runningForManifold) && ((indexFirstToken(ex) == 34) || (indexFirstToken(ex) == -1)))
            executeLineManifold(ex, ap);
        else
            executeLineNormal(ex, ap);
    }

    public static boolean flagChange_Prototypes(String ex) {
        int i;
        i = Prototype.indexInTokens_Prototypes(ex);
        if (i == -1)
            return false;
        else {
            if (i == 0)
                bClosest_Center_To_Cluster_Number = true;
            else if (i == 1)
                bClosest_Center_To_Cluster_Number = false;
            else if (i == 2)
                bClosest_Center_To_Cluster_Points = true;
            else if (i == 3)
                bClosest_Center_To_Cluster_Points = false;
            else if (i == 4) {
                bClosest_Center_To_Normalized_Distortions = true;
                /*
                 * Enumeration extensions =
				 * Prototype.Closest_Center_To_Cluster_Number.keys(); Integer R,
				 * vp; while (extensions.hasMoreElements()) { R = (Integer)
				 * extensions.nextElement(); vp = (Integer)
				 * Prototype.Closest_Center_To_Cluster_Number.get(R);
				 * System.out.print("(" + R + ":: " + vp + ")"); }
				 */
            } else if (i == 5)
                bClosest_Center_To_Normalized_Distortions = false;
            else if (i == 6)
                bCluster_Number_To_Closest_Center = true;
            else if (i == 7)
                bCluster_Number_To_Closest_Center = false;
            else if (i == 8) {
                runningForPrototypes = false;
                Prototype.Loading_From_Scenario_End = true;
                Prototype.Prototypes_Selected = true;
                ControlProcess.put("prototypeSelected", true);

				/*
                 * Enumeration extensions =
				 * Prototype.Closest_Center_To_Cluster_Number.keys(); Integer R,
				 * vp; while (extensions.hasMoreElements()) { R = (Integer)
				 * extensions.nextElement(); vp = (Integer)
				 * Prototype.Closest_Center_To_Cluster_Number.get(R);
				 * System.out.print("(" + R + ": " + vp + ")"); }
				 * System.exit(0);
				 */
            } else if (i == 9)
                bClosest_Center_Ordered_List = true;
            else if (i == 10)
                bClosest_Center_Ordered_List = false;

            // System.out.println("Bluk proto = " + i);

            return true;
        }
    }

    public static boolean flagChange_Manifold(String ex, AGCT ap) {
        if (ex.startsWith(AGCT.Token_Ordered_List_Names_Begin))
            bOrdered_List_Names = true;
        else if (ex.startsWith(AGCT.Token_Ordered_List_Names_End))
            bOrdered_List_Names = false;
        else if (ex.startsWith(AGCT.Token_Matrix_W_Begin))
            bMatrix_W = true;
        else if (ex.startsWith(AGCT.Token_Matrix_W_End))
            bMatrix_W = false;
        else if (ex.startsWith(AGCT.Token_Manifold_Eigenvalues_Begin))
            bManifold_Eigenvalues = true;
        else if (ex.startsWith(AGCT.Token_Manifold_Eigenvalues_End))
            bManifold_Eigenvalues = false;
        else if (ex.startsWith(AGCT.Token_Matrix_M_Begin))
            bMatrix_M = true;
        else if (ex.startsWith(AGCT.Token_Matrix_M_End))
            bMatrix_M = false;
        else if (ex.startsWith(allKeywords[35])) {
            runningForManifold = false;
            AGCT.Loading_From_Scenario_Manifold_End = true;
            // ControlProcess.put("manifoldProcessed",true);
            ap.justAfterManifold();

            // System.out.println("Bluk mani = " + i);
        } else
            return false;
        return true;
        //		}
    }

    public static void executeLinePrototypes(String ex, AGCT ap) {
        System.out.println("Scenario.executeLinePrototypes()");
        if (!flagChange_Prototypes(ex)) {
            Prototype.executeScenario(ex, ap);
        }
    }

    /**
     * 一行exを読み、その命令を実行する
     *
     * @param ex
     * @param ap
     */
    public static void executeLineManifold(String ex, AGCT ap) {
        if (!flagChange_Manifold(ex, ap)) {
            AGCT.executeManifold(ex, ap);
        }
    }

    private static void executeLineNormal(String ex, AGCT ap) {
        Debug.debug("Scenario # executeLineNormal", ex);
        StringTokenizer t;
        String s, nam;
        int dum, mI, mJ, sel, index;
        double ddum;
        boolean bdum, isSel;
        t = new StringTokenizer(ex, Scenario.myToken);
        s = t.nextToken();
        index = index(s);

        // blocks until conditions are met to run the process
        block(index);

        if (index == 0) {
            ap.goLoadDomain(t.nextToken());
        } else if (index == 1) {
            dum = Integer.parseInt(t.nextToken());
            ap.myTabbedPane.mySelectionPane.setMethod_F(dum);
            ap.myTabbedPane.mySelectionPane.featureMethod.setSelectedIndex(dum);
        } else if (index == 2) {
            ap.myTabbedPane.mySelectionPane.goValidate();
        } else if (index == 3) {
            nam = t.nextToken();
            mI = Integer.parseInt(t.nextToken());
            mJ = Integer.parseInt(t.nextToken());
            sel = Integer.parseInt(t.nextToken());

            isSel = false;
            if (sel == 1)
                isSel = true;
            else if (sel == 0)
                isSel = false;
            else
                Matrix.perror("Scenario :: non boolean selection");

            if (nam.equals(JAGCTSelectionPane.refGene))
                ap.myTabbedPane.mySelectionPane.checkBoxGene.get(mI).setSelected(isSel);
            else if (nam.equals(JAGCTSelectionPane.refLigand))
                ap.myTabbedPane.mySelectionPane.checkBoxLigand.get(mI).setSelected(isSel);
            else if (nam.equals(JAGCTSelectionPane.refGroup))
                ap.myTabbedPane.mySelectionPane.checkBoxGroup.get(mJ).setSelected(isSel);
        } else if (index == 4) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationMethod_F(dum);
        } else if (index == 5) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationMethod_W(dum);
        } else if (index == 6) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationMethod_S(dum);
        } else if (index == 7) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationMethod_N(dum);
        } else if (index == 8) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationNumber_Of_Neighbors(dum);
        } else if (index == 9) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationNumber_Of_Manifold_Components(dum);
        } else if (index == 10) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationNumber_Of_Triangulation_Dimensions(dum);
        } else if (index == 11) {
            ddum = Double.parseDouble(t.nextToken());
            ap.requestModificationT_Heat_Kernel(ddum);
        } else if (index == 12) {
            ddum = Double.parseDouble(t.nextToken());
            ap.requestModificationT_Heat_Kernel_Similarity_Neighbors(ddum);
        } else if (index == 13) {
            ddum = Double.parseDouble(t.nextToken());
            ap.requestModificationSparsify_Statistic(ddum);
        } else if (index == 14) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationBregDiv(dum);
        } else if (index == 15) {
            ddum = Double.parseDouble(t.nextToken());
            ap.requestModificationLimit_P_Chi2(ddum);
        } else if (index == 16) {
            dum = Integer.parseInt(t.nextToken());
            bdum = int2boolean(dum);
            ap.requestModificationUseDebug(bdum);
        } else if (index == 17) {
            dum = Integer.parseInt(t.nextToken());
            bdum = int2boolean(dum);
            ap.requestModificationUseWarning(bdum);
        } else if (index == 18) {
            dum = Integer.parseInt(t.nextToken());
            bdum = int2boolean(dum);
            ap.requestModificationPerspective(bdum);
        } else if (index == 19) {
            dum = Integer.parseInt(t.nextToken());
            bdum = int2boolean(dum);
            ap.requestModificationSort_Depth(bdum);
        } else if (index == 20) {
            dum = Integer.parseInt(t.nextToken());
            bdum = int2boolean(dum);
            ap.requestModificationUse_Shadow(bdum);
        } else if (index == 21) {
            ap.myTabbedPane.myManifoldPane.goManifold();
        } else if (index == 22) {
            ap.myTabbedPane.myPCAPane.goPCA();
        } else if (index == 23) {
            Debug.debug("Scenario : triangulation");
            ap.myTabbedPane.myManifoldPane.goTriangulation();
        } else if (index == 24) {// clustering
            try {
                AGCTFileWriter.loadClustering(new File(t.nextToken()), AGCT.getInstance());
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            //			ap.nonThreadClustering(t.nextToken());
        } else if (index == 25) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationNumber_Of_Wavelet_Stamps(dum);
            ap.myTabbedPane.mySelectionPane.moreNumber_Of_Wavelet_Stamps();
        } else if (index == 26) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationMax_Number_Of_Features(dum);
            ap.myTabbedPane.mySelectionPane.moreMax_Number_Of_Features();
        } else if (index == 27) {
            dum = Integer.parseInt(t.nextToken());
            ap.myTabbedPane.mySelectionPane.playMethod_FS(dum);
        } else if (index == 28) {
            dum = Integer.parseInt(t.nextToken());
            ap.requestModificationMax_Number_Of_Prototypes(dum);
            ap.myTabbedPane.mySelectionPane.moreMax_Number_Of_Prototypes();
        } else if (index == 29) {
            dum = Integer.parseInt(t.nextToken());
            ap.myTabbedPane.mySelectionPane.playMethod_PS(dum);
        } else if (index == 31) {
            runningForPrototypes = true;
            Prototype.flushAll();
            Prototype.Loading_From_Scenario_Begin = true;
            bClosest_Center_Ordered_List = bClosest_Center_To_Cluster_Number = bCluster_Number_To_Closest_Center = bClosest_Center_To_Cluster_Points = bClosest_Center_To_Normalized_Distortions = false;
        } else if (index == 32) {
            runningForPrototypes = false;
        } else if (index == 34) {
            runningForManifold = true;
            AGCT.Loading_From_Scenario_Manifold_Begin = true;
            bOrdered_List_Names = bMatrix_W = bManifold_Eigenvalues = bMatrix_M = false;
            ap.justBeforeManifold();
        } else if (index == 35) {
            ap.justAfterManifold();
            runningForManifold = false;
        } else if (index == 36) {
            //			ap.
        } else if (index == 39) {// "JAGCTClickLigandCheckBox"
            int i = Integer.parseInt(t.nextToken());
            for (int j = 0; JAGCTSelectionPane.getLigandSelectionBoxes().size() <= i; j++)
                if (j % 1000000000 == 0) Debug.debug("JAGCTClickLigandCheckBox", j);
            JAGCTSelectionPane.getLigandSelectionBoxes().get(i).doClick();
        } else if (index == 40) {
            double sd = Double.parseDouble(t.nextToken());
            ap.setSDThreshold(sd);
        }
    }
}