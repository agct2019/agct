package test;

import clustering.AGCTClustering_KM;
import gene.Gene;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

public class Prototype implements Debuggable {
    private Prototype() {
    }

    public static Vector Closest_Center_Ordered_List;
    // Center names are ordered, to be used for matching with M, ...

    public static final String Token_Closest_Center_Ordered_List_Begin = "@Closest_Center_Ordered_List_Begin";
    public static final String Token_Closest_Center_Ordered_List_End = "@Closest_Center_Ordered_List_end";

    public static Hashtable Closest_Center_To_Cluster_Number;
    // Entries = selected genes index (Integer), output = Integers = cluster
    // number

    public static final String Token_Closest_Center_To_Cluster_Number_Begin = "@Closest_Center_To_Cluster_Number_Begin";
    public static final String Token_Closest_Center_To_Cluster_Number_End = "@Closest_Center_To_Cluster_Number_End";

    public static Hashtable<Integer, Vector> Closest_Center_To_Cluster_Points;
    // Entries = selected genes index (Integer), output = Vector = {selected
    // genes belonging to the cluster}

    public static String Token_Closest_Center_To_Cluster_Points_Begin = "@Closest_Center_To_Cluster_Points_Begin";
    public static String Token_Closest_Center_To_Cluster_Points_End = "@Closest_Center_To_Cluster_Points_End";

    public static Hashtable Closest_Center_To_Normalized_Distortions;
    // Entries = selected genes index (Integer), output = Vector = {distortion
    // of selected genes belonging to the cluster wrt to representative,
    // normalized in [0,1] = [min, max]}

    public static String Token_Closest_Center_To_Normalized_Distortions_Begin = "@Closest_Center_To_Normalized_Distortions_Begin";
    public static String Token_Closest_Center_To_Normalized_Distortions_End = "@Closest_Center_To_Normalized_Distortions_End";

    public static Hashtable Cluster_Number_To_Closest_Center;
    // keys = cluster number (Integer), returns the index (Integer) of the
    // closest center

    public static String Token_Cluster_Number_To_Closest_Center_Begin = "@Cluster_Number_To_Closest_Center_Begin";
    public static String Token_Cluster_Number_To_Closest_Center_End = "@Cluster_Number_To_Closest_Center_End";

    // //////////////////////////////////////////////////////////////////////////////////////////////////////
    // WARNING : gene indexes are directly into domainGenes (useless to use
    // selectedGeneNumberToGeneNumber*

    public static final String allTokens_Prototype[] = {Token_Closest_Center_To_Cluster_Number_Begin, Token_Closest_Center_To_Cluster_Number_End, Token_Closest_Center_To_Cluster_Points_Begin,
            Token_Closest_Center_To_Cluster_Points_End, Token_Closest_Center_To_Normalized_Distortions_Begin, Token_Closest_Center_To_Normalized_Distortions_End,
            Token_Cluster_Number_To_Closest_Center_Begin, Token_Cluster_Number_To_Closest_Center_End, Scenario.allKeywords[32], Token_Closest_Center_Ordered_List_Begin,
            Token_Closest_Center_Ordered_List_End};

    public static boolean Prototypes_Selected;
    public static boolean No_Reduction;
    public static boolean Loading_From_Scenario_Begin, Loading_From_Scenario_End;

    public static int indexInTokens_Prototypes(String s) {
        String ref, t;
        int i;
        for (i = 0; i < allTokens_Prototype.length; i++) {
            ref = allTokens_Prototype[i];
            if (s.length() >= ref.length()) {
                t = s.substring(0, ref.length());
                if (t.equals(ref))
                    return i;
            }
        }
        return -1;
    }

    public static void assertComplete() {
        if ((Closest_Center_To_Cluster_Number == null) || (Closest_Center_To_Cluster_Number.size() < 1))
            Matrix.perror("Prototype.class :: Closest_Center_To_Cluster_Number is empty");

        if ((Closest_Center_To_Cluster_Points == null) || (Closest_Center_To_Cluster_Points.size() < 1))
            Matrix.perror("Prototype.class :: Closest_Center_To_Cluster_Points is empty");

        if ((Closest_Center_To_Normalized_Distortions == null) || (Closest_Center_To_Normalized_Distortions.size() < 1))
            Matrix.perror("Prototype.class :: Closest_Center_To_Normalized_Distortions is empty");

        if ((Cluster_Number_To_Closest_Center == null) || (Cluster_Number_To_Closest_Center.size() < 1))
            Matrix.perror("Prototype.class :: Cluster_Number_To_Closest_Center is empty");

        if ((Closest_Center_Ordered_List == null) || (Closest_Center_Ordered_List.size() < 1))
            Matrix.perror("Prototype.class :: Closest_Center_Ordered_List is empty");
    }

    public static void flushAll() {
        if (AGCT.MYDEBUG)
            AGCT.debug("Prototype.flushAll()");
        Prototypes_Selected = false;
        No_Reduction = false;
        Closest_Center_To_Cluster_Number = null;
        Cluster_Number_To_Closest_Center = null;
        Closest_Center_To_Cluster_Points = null;
        Closest_Center_To_Normalized_Distortions = null;
        Closest_Center_Ordered_List = null;
        Loading_From_Scenario_Begin = Loading_From_Scenario_End = false;
    }

    public static void executeScenario(String ex, AGCT ap) {
        StringTokenizer t;
        String NP, s;
        Vector v;
        int nc;

        if (Scenario.bClosest_Center_Ordered_List) {
            Closest_Center_Ordered_List = new Vector();
            t = new StringTokenizer(ex, Scenario.myToken);
            while (t.hasMoreTokens()) {
                Closest_Center_Ordered_List.addElement(new Integer(ap.data.getMyDomain().getGenes().get(t.nextToken()).getIdNumber()));
//						getGeneIndex(t.nextToken())));
            }
        } else if (Scenario.bClosest_Center_To_Cluster_Number) {
            if (Closest_Center_To_Cluster_Number == null)
                Closest_Center_To_Cluster_Number = new Hashtable();
            t = new StringTokenizer(ex, Scenario.myToken);
            NP = t.nextToken();
            nc = Integer.parseInt(t.nextToken());
            Closest_Center_To_Cluster_Number.put(new Integer(ap.data.getMyDomain().getGenes().get(NP).getIdNumber())
//					getGeneIndex(NP))
                    , new Integer(nc));

            // System.out.print("(" + NP + ", " + ap.myDomain.getGeneIndex(NP) +
            // ") ");
        } else if (Scenario.bCluster_Number_To_Closest_Center) {
            if (Cluster_Number_To_Closest_Center == null)
                Cluster_Number_To_Closest_Center = new Hashtable();
            t = new StringTokenizer(ex, Scenario.myToken);
            nc = Integer.parseInt(t.nextToken());
            NP = t.nextToken();
            Cluster_Number_To_Closest_Center.put(new Integer(nc), new Integer(ap.data.getMyDomain().getGenes().get(NP).getIdNumber()));
        } else if (Scenario.bClosest_Center_To_Cluster_Points) {
            if (Closest_Center_To_Cluster_Points == null)
                Closest_Center_To_Cluster_Points = new Hashtable();
            t = new StringTokenizer(ex, Scenario.myToken);
            v = new Vector();
            NP = t.nextToken();
            while (t.hasMoreTokens()) {
                s = t.nextToken();
                v.addElement(new Integer(ap.data.getMyDomain().getGenes().get(s).getIdNumber()));
            }
            Closest_Center_To_Cluster_Points.put(new Integer(ap.data.getMyDomain().getGenes().get(NP).getIdNumber()), v);
        } else if (Scenario.bClosest_Center_To_Normalized_Distortions) {
            if (Closest_Center_To_Normalized_Distortions == null)
                Closest_Center_To_Normalized_Distortions = new Hashtable();
            t = new StringTokenizer(ex, Scenario.myToken);
            v = new Vector();
            NP = t.nextToken();
            while (t.hasMoreTokens()) {
                s = t.nextToken();
                v.addElement(new Double(Double.parseDouble(s)));
            }
            Closest_Center_To_Normalized_Distortions.put(new Integer(ap.data.getMyDomain().getGenes().get(NP).getIdNumber()), v);
        } else
            Matrix.perror("Prototype.class :: not in a prototype scenario");
    }

    public static void prototypeSelection(Domain d) {
        if (AGCT.MYDEBUG)
            AGCT.debug("Prototype.prototypeSelection(domain)");

        int numberOfPrototypesBeforeSelection = getNumberOfPrototypesBeforeSelection(AGCT.getInstance());
        if (AGCT.MYDEBUG)
            AGCT.debug(AGCT.Method_PS, AGCT.Max_Number_Of_Prototypes, numberOfPrototypesBeforeSelection);

        if ((AGCT.Method_PS > 0) && (AGCT.Max_Number_Of_Prototypes < numberOfPrototypesBeforeSelection))
            prototypeAggregation_KM(d);
        else
            No_Reduction = true;

        Prototypes_Selected = true;
    }

    /**
     * Cluster_Center_To_Ordered_Listを計算している。
     * 他にもやっているはず。TODO 追記
     *
     * @param d
     */
    public static void prototypeAggregation_KM(Domain d) {
        Clustering clustering = new Clustering(AGCT.getInstance(), -1);
        clustering.indexReferenceStringCenters = 3;

        if (AGCT.Max_Number_Of_Prototypes < 0)
            Matrix.perror("Negative number of prototypes");

        Closest_Center_To_Cluster_Number = new Hashtable();
        Closest_Center_To_Cluster_Points = new Hashtable();
        Closest_Center_To_Normalized_Distortions = new Hashtable();
        Cluster_Number_To_Closest_Center = new Hashtable();
        Closest_Center_Ordered_List = new Vector();

        AGCT.getInstance().timerTickOn();

        AGCTClustering_KM c = new AGCTClustering_KM(clustering, AGCT.getInstance(), AGCT.Max_Number_Of_Prototypes, 1, false, 0, -1.0, true, 30, 0);
        c.hard_memberships_centers = new Vector();

        c.toClustering();

        c.toClusterNumberToClosestCenter(Cluster_Number_To_Closest_Center, Closest_Center_Ordered_List);
        c.toClosestCenterToClusterNumber(Cluster_Number_To_Closest_Center, Closest_Center_To_Cluster_Number);
        c.toClosestCenterToClusterPoints(Closest_Center_To_Cluster_Number, Cluster_Number_To_Closest_Center, Closest_Center_To_Cluster_Points, Closest_Center_To_Normalized_Distortions);

        orderAndNormalize(d);

        AGCT.Max_Number_Of_Prototypes = c.getNumberOfClusters();
        AGCT.getInstance().myTabbedPane.mySelectionPane.computeLabel();
        // AGCT.getInstance().myAnnotationFrame.onlyPrototypes.setEnabled(true);
        AGCT.getInstance().timerTickOff();
    }

    public static void affiche() {
        Enumeration extensions = Cluster_Number_To_Closest_Center.keys();
        Integer R, S;
        Vector vp, vd;
        int i;
        if (extensions != null) {
            System.out.println("Ordered list of centers ::");
            for (i = 0; i < Closest_Center_Ordered_List.size(); i++)
                System.out.println(Closest_Center_Ordered_List.get(i));
            System.out.println("\n");

            while (extensions.hasMoreElements()) {
                R = (Integer) extensions.nextElement();
                S = (Integer) Cluster_Number_To_Closest_Center.get(R);

                System.out.print("Cluster " + R + " --- center = " + S + " :: ");
                vp = (Vector) Closest_Center_To_Cluster_Points.get(S);
                vd = (Vector) Closest_Center_To_Normalized_Distortions.get(S);

                if (vp != null)
                    for (i = 0; i < vp.size(); i++)
                        System.out.print(" [" + (Integer) vp.get(i) + ", " + (Double) vd.get(i) + "] ");
                System.out.println();
            }
        }
    }

    public static void orderAndNormalize(Domain d) {
        Enumeration extensions = Closest_Center_To_Cluster_Points.keys();
        Integer R;
        Integer TI;
        Double TD;
        Vector vp, vd;
        int i, j;
        double di, dj, max;
        if (extensions != null) {
            AGCTCounter cc = new AGCTCounter("Working on prototypes", Closest_Center_To_Cluster_Points.size());
            while (extensions.hasMoreElements()) {
                R = (Integer) extensions.nextElement();
                vp = (Vector) Closest_Center_To_Cluster_Points.get(R);
                vd = (Vector) Closest_Center_To_Normalized_Distortions.get(R);

                if (vp.size() != vd.size())
                    Matrix.perror("Prototype.class :: size mismatch");

                for (i = 0; i < vp.size() - 1; i++)
                    for (j = i + 1; j < vp.size(); j++) {
                        di = ((Double) vd.get(i)).doubleValue();
                        dj = ((Double) vd.get(j)).doubleValue();
                        if (dj < di) {
                            TD = new Double(((Double) vd.get(i)).doubleValue());
                            vd.setElementAt((Double) vd.get(j), i);
                            vd.setElementAt(TD, j);

                            TI = new Integer(((Integer) vp.get(i)).intValue());
                            vp.setElementAt((Integer) vp.get(j), i);
                            vp.setElementAt(TI, j);
                        }
                    }

                if (vd.size() > 1) {
                    max = ((Double) vd.get(vd.size() - 1)).doubleValue();
                    if (max > 0.0)
                        for (i = 0; i < vd.size(); i++) {
                            di = (((Double) vd.get(i)).doubleValue()) / max;
                            TD = new Double(di);
                            vd.setElementAt(TD, i);
                        }
                }

                cc.increment();
            }
            cc.end();
        }
    }

    public static boolean currentCentersContain(int id) {
        if (Cluster_Number_To_Closest_Center.size() == 0)
            return false;

        Enumeration extensions = Cluster_Number_To_Closest_Center.keys();
        Integer R, S;
        if (extensions != null) {
            while (extensions.hasMoreElements()) {
                R = (Integer) extensions.nextElement();
                S = (Integer) Cluster_Number_To_Closest_Center.get(R);

                if (id == S.intValue())
                    return true;
            }
        }

        return false;
    }

    public static int getNumberOfPrototypesBeforeSelection(AGCT a) {
        return a.myTabbedPane.mySelectionPane.numberCurrentSelectableGenes();
    }

    public static boolean geneOrNeighborsIsReferencedPrototype(AGCT a, int gid) {
        boolean ref = false;
        int i;
        Vector l;
        Gene nn, gg = (Gene) a.data.getMyDomain().getGenes().get(a.data.getMyDomain().selectedGeneNumberToGeneNumber[gid]);
        if (!No_Reduction) {
            l = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(a.data.getMyDomain().selectedGeneNumberToGeneNumber[gid]));
            for (i = 0; i < l.size(); i++) {
                nn = (Gene) a.data.getMyDomain().getGenes().get(((Integer) l.get(i)).intValue());
                if (nn.isReferenced())
                    return true;
            }
            return false;
        } else
            return gg.isReferenced();
    }
}