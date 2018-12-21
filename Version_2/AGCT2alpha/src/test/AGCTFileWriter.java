package test;

import clustering.AGCTClustering_KM;
import forDebug.Debug;
import gene.Gene;
import ligand.Ligand;

import java.io.*;
import java.util.*;

public class AGCTFileWriter implements Debuggable {
    public AGCTFileWriter() {
    }

    public static String DATA_Ligand = new String("@LIGAND"), DATA_Domain = "@DOMAIN_NAME", DATA_Groups = new String("@GROUPS"),
            DATA_Data = new String("@DATA"), DATA_Annotation = new String("@ANNOTATION"), DATA_GeneAnnotations = new String("@GENE_ANNOTATIONS"), DATA_Selection = new String(
            "@SELECTION_REDUCTION_FILE"), DATA_Comments = new String("//"), DATA_QuestionMark = new String("?");

    public static String ANNO_Data = new String("@DATA"), ANNO_Format = new String("@FORMAT"), ANNO_Variable_Definition = new String("@VARIABLES"), ANNO_Empty = new String(
            "@EMPTY");

    public static String ANNO_Format_AGCT = new String("AGCT"), ANNO_Format_GEO = new String("GEO");

    public static String Clustering_Born_String = new String("@Reference_String");

    public static String ANNO_Variable_Definition_Dum = new String("DUM"), ANNO_Variable_Definition_ID = new String("ID"),
            ANNO_Variable_Definition_Ascii_Names = new String("NAMES"), ANNO_Variable_Definition_Ontology_F = new String("F"),
            ANNO_Variable_Definition_Ontology_P = new String("P"), ANNO_Variable_Definition_Ontology_C = new String("C"), ANNO_Variable_Definition_Referenced = new String(
            "REF");

    public static String[][] REPLACE_Old_New = {{"///", "BIG_SEP"}, {"//", "SMALL_SEP"}};

    private AGCT myAGCT;

    AGCTFileWriter(AGCT agct) {
        if (AGCT.MYDEBUG)
            AGCT.debug("AGCTFileWriter.AGCTFileWriter(agct)");
        myAGCT = agct;
    }

    public void toSavingScenario(File rf) {
        FileWriter f = null;
        int i;
        String s;
        if (Scenario.allStrings != null)
            try {
                f = new FileWriter(rf);
                Debug.debug("toSavingScenario : allStrings", Scenario.allStrings);
                for (i = 0; i < Scenario.allStrings.size(); i++) {
                    s = Scenario.allStrings.get(i);
                    if (Scenario.isKeepPrototypes(s))
                        morePrototypes(f);
                    else if (Scenario.isKeepManifold(s))
                        moreManifold(f);
                    else
                        f.write(s + "\n");
                }

                toSavingClustering(f, 0);/*０のみ*/

                f.close();
            } catch (IOException e) {
                Matrix.perror("AGCTFileWriter.class :: Writing scenario error");
            }
    }

    public void toSavingManifold(File rf) {
        FileWriter f = null;
        int i;
        String s;
        try {
            f = new FileWriter(rf);

            ManifoldSave3DPts(f);

            f.close();
        } catch (IOException e) {
            Matrix.perror("AGCTFileWriter.class :: Writing manifold error");
        }
    }

    // frank
    public void ManifoldSave3DPts(FileWriter rf) throws IOException {
        int i, j, dX, dY;
        double d;
        Gene g;

        for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
            g = ((Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]));

            Pnt pt = g.manifold_point3D;
            String str = "";

            for (j = 0; j < pt.coordinates.length - 1; j++) {
                str += pt.coordinates[j] + "\t";
            }

            rf.write(g.myID + "\t" + str + System.getProperty("line.separator"));

        }

    }

    public void moreManifold(FileWriter rf) throws IOException {
        Debug.debug("AGCTFileWriter # moreManifold", rf);
        int i, j, dX, dY;
        double d;
        Gene g;

        if (ControlProcess.hasTrue("manifoldProcessed")) {
            rf.write(Scenario.allKeywords[34] + Scenario.myToken + "\n");

            // Ordered list name
            rf.write(AGCT.Token_Ordered_List_Names_Begin + "\n");
            for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                g = ((Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]));
                rf.write(g.myID);
                if (i < myAGCT.data.getMyDomain().numberSelectedGenes - 1)
                    rf.write(Scenario.myToken);
                else
                    rf.write("\n");
            }
            rf.write(AGCT.Token_Ordered_List_Names_End + "\n");

            // W
            rf.write(AGCT.Token_Matrix_W_Begin + "\n");

            // 0でない部分に関して、id|val という形式で書かれている
            dX = myAGCT.data.getW().rowCount();
            dY = myAGCT.data.getW().colCount();
            for (i = 0; i < dX; i++) {
                boolean first = true;
                for (j = 0; j < dY; j++) {
                    d = myAGCT.data.getW().get(i, j);
                    if (d != 0) {
                        if (!first)
                            rf.write(Scenario.myToken);
                        first = false;
                        rf.write(j + Scenario.myToken + d);
                    }
                }
                rf.write("\n");
            }
            rf.write(AGCT.Token_Matrix_W_End + "\n");

            // Manifold Eigenvalues
            rf.write(AGCT.Token_Manifold_Eigenvalues_Begin + "\n");
            dX = myAGCT.data.getManifold_eigenvalues().length;
            for (i = 0; i < dX; i++) {
                d = myAGCT.data.getManifold_eigenvalues()[i];
                rf.write(d + "");
                if (i < dX - 1)
                    rf.write(Scenario.myToken);
                else
                    rf.write("\n");
            }
            rf.write(AGCT.Token_Manifold_Eigenvalues_End + "\n");

            // M
            rf.write(AGCT.Token_Matrix_M_Begin + "\n");
            // 0でない部分に関して、id|val という形式で書かれている
            dX = myAGCT.data.getM().rowCount();
            dY = myAGCT.data.getM().colCount();
            for (i = 0; i < dX; i++) {
                boolean first = true;
                for (j = 0; j < dY; j++) {
                    d = myAGCT.data.getM().get(i, j);
                    if (d != 0) {
                        rf.write((first ? "" : Scenario.myToken) + j + Scenario.myToken + d);
                        first = false;
                    }
                }
                rf.write("\n");
            }
            rf.write(AGCT.Token_Matrix_M_End + "\n");

            rf.write(Scenario.allKeywords[35] + Scenario.myToken + "\n");
        }
    }

    public void moreClustering(FileWriter rf) throws IOException {
        if (ControlProcess.hasTrue("clusteringProcessed")) {
            rf.write(Scenario.allKeywords[37] + Scenario.myToken + "\n");

            rf.write(Scenario.allKeywords[38] + Scenario.myToken + "\n");
        }
    }

    public void morePrototypes(FileWriter rf) throws IOException {
        int i, j, iter;
        double d;
        Integer I;
        Enumeration extensions;
        Hashtable t;
        Gene gg;
        Vector v;

        if (!Prototype.No_Reduction) {
            Prototype.assertComplete();

            rf.write(Scenario.allKeywords[31] + Scenario.myToken + "\n");
            // Closest_Center_Ordered_List

            rf.write(Prototype.Token_Closest_Center_Ordered_List_Begin + "\n");
            v = Prototype.Closest_Center_Ordered_List;
            for (j = 0; j < v.size(); j++) {
                I = (Integer) v.get(j);
                i = I.intValue();
                gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(i);
                rf.write(gg.myID + "");
                if (j < v.size() - 1)
                    rf.write(Scenario.myToken);
            }
            rf.write("\n");
            rf.write(Prototype.Token_Closest_Center_Ordered_List_End + "\n");

            // Closest_Center_To_Cluster_Number
            rf.write(Prototype.Token_Closest_Center_To_Cluster_Number_Begin + "\n");
            t = Prototype.Closest_Center_To_Cluster_Number;
            extensions = t.keys();
            while (extensions.hasMoreElements()) {
                I = (Integer) extensions.nextElement();
                i = I.intValue();
                gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(i);
                rf.write(gg.myID + Scenario.myToken + t.get(I) + "\n");
            }
            rf.write(Prototype.Token_Closest_Center_To_Cluster_Number_End + "\n");

            // Cluster_Number_To_Closest_Center
            rf.write(Prototype.Token_Cluster_Number_To_Closest_Center_Begin + "\n");
            t = Prototype.Cluster_Number_To_Closest_Center;
            extensions = t.keys();
            while (extensions.hasMoreElements()) {
                I = (Integer) extensions.nextElement();
                i = ((Integer) t.get(I)).intValue();
                gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(i);
                rf.write(I + Scenario.myToken + gg.myID + "\n");
            }
            rf.write(Prototype.Token_Cluster_Number_To_Closest_Center_End + "\n");

            // Closest_Center_To_Cluster_Points
            rf.write(Prototype.Token_Closest_Center_To_Cluster_Points_Begin + "\n");
            t = Prototype.Closest_Center_To_Cluster_Points;
            extensions = t.keys();
            while (extensions.hasMoreElements()) {
                I = (Integer) extensions.nextElement();
                i = I.intValue();
                gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(i);
                rf.write(gg.myID + Scenario.myToken);
                v = (Vector) t.get(I);
                for (i = 0; i < v.size(); i++) {
                    j = ((Integer) v.get(i)).intValue();
                    gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(j);
                    rf.write(gg.myID);
                    if (i < v.size() - 1)
                        rf.write(Scenario.myToken);
                }
                rf.write("\n");
            }
            rf.write(Prototype.Token_Closest_Center_To_Cluster_Points_End + "\n");

            // Closest_Center_To_Normalized_Distortions
            rf.write(Prototype.Token_Closest_Center_To_Normalized_Distortions_Begin + "\n");
            t = Prototype.Closest_Center_To_Normalized_Distortions;
            extensions = t.keys();
            while (extensions.hasMoreElements()) {
                I = (Integer) extensions.nextElement();
                i = I.intValue();
                gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(i);
                rf.write(gg.myID + Scenario.myToken);
                v = (Vector) t.get(I);
                for (i = 0; i < v.size(); i++) {
                    d = ((Double) v.get(i)).doubleValue();
                    rf.write(d + "");
                    if (i < v.size() - 1)
                        rf.write(Scenario.myToken);
                }
                rf.write("\n");
            }
            rf.write(Prototype.Token_Closest_Center_To_Normalized_Distortions_End + "\n");

            rf.write(Scenario.allKeywords[32] + Scenario.myToken + "\n");
        }
    }

    public void toSaving(File rf) {
        FileWriter f = null;
        try {
            int i, j, k, nn = 0, gid, nid, iii, totGenes;
            boolean better;
            Gene gg, hh;
            String s, dums;
            Vector vid, did, littleVector;
            f = new FileWriter(rf);

            f.write(separator());
            f.write("Results for domain " + myAGCT.data.getMyDomain().domainName);
            f.write(separator());

            f.write("\n\n");

            f.write(separator());
            f.write("Processing parameters:\n");
            f.write(myAGCT.getParameterString());
            f.write(separator());

            if ((myAGCT.data.getMyDomain().getGenes() != null) && (myAGCT.data.getMyDomain().getLigands().numberOfSelectedLigands() > 0)) {

                f.write(separator());
                f.write("** Ligands\n\nselected (" + myAGCT.data.getMyDomain().getLigands().numberOfSelectedLigands() + ") :\n");
                for (Ligand ligand : myAGCT.data.getMyDomain().getLigands()) {
                    if (ligand.isChecked()) {
                        s = ligand.getName();
                        f.write(s + " ");
                    }
                }
                f.write(separator());
                f.write("\n\n");
            }

            if ((myAGCT.data.getMyDomain().getGenes() != null) && (myAGCT.data.getMyDomain().numberSelectedGenes > 0)) {
                f.write(separator());
                if (Prototype.No_Reduction) {
                    f.write("** Genes\n\nSelected (" + myAGCT.data.getMyDomain().numberSelectedGenes + "):\n");
                } else {
                    f.write("** Prototypes\n\nSelected (" + myAGCT.data.getMyDomain().numberSelectedGenes + ") + their clusters:\n");
                }

                for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                    gid = myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i];
                    gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(gid);
                    f.write(gg.myID + ((gg.asciiName != null) ? "_" + gg.asciiName + "  " : "  "));
                    if (Prototype.No_Reduction)
                        f.write(" ");
                    else {
                        f.write(" [ ");
                        vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
                        did = (Vector) Prototype.Closest_Center_To_Normalized_Distortions.get(new Integer(gid));
                        for (iii = 0; iii < vid.size(); iii++) {
                            nid = ((Integer) vid.get(iii)).intValue();
                            hh = myAGCT.getDirectGene(nid);
                            f.write("(" + hh.myID + ((hh.asciiName != null) ? "_" + hh.asciiName + "  " : "  ") + "  :  " + did.get(iii) + ")");
                        }
                        f.write(" ]\n");
                    }
                }
                f.write(separator());
                f.write("\n\n");
            }

            if (ControlProcess.hasTrue("manifoldProcessed")) {
                f.write(separator());
                f.write("** Manifold\n\n");

                f.write("Eigenvalues (sorted): ");

                for (i = 0; i < myAGCT.data.getManifold_eigenvalues().length; i++)
                    f.write(DF.format(myAGCT.data.getManifold_eigenvalues()[i]) + " ");

                if (myAGCT.data.isAnnotationsExist()) {
                    f.write("\n\nNatural Neighbors for annotations\n\n");
                    for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                        gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
                        nn = 0;
                        if (((gg.annotationsF != null) || (gg.annotationsP != null) || (gg.annotationsC != null)) && (gg.neighbors != null)) {
                            for (j = 0; j < gg.neighbors.size(); j++) {
                                hh = (Gene) myAGCT.data.getMyDomain().getGenes().get(
                                        myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[((Integer) gg.neighbors.get(j)).intValue()]);
                                if ((hh.annotationsF != null) || (hh.annotationsP != null) || (hh.annotationsC != null)) {
                                    if (nn == 0) {
                                        f.write("  ");
                                        if (gg.annotationsF != null) {
                                            f.write(JAnnotationFrame.TAG_F + " ");
                                            for (k = 0; k < gg.annotationsF.size(); k++) {
                                                if (!((String[]) gg.annotationsF.get(k))[0].equals(Domain.ANNOTATION_F))
                                                    Matrix.perror("AGCTFileWriter.class :: annotation != F in F annotations");
                                                f.write(((String[]) gg.annotationsF.get(k))[1]);
                                                if (k < gg.annotationsF.size() - 1)
                                                    f.write(", ");
                                            }
                                            f.write(";\n");
                                        }
                                        if (gg.annotationsP != null) {
                                            f.write(JAnnotationFrame.TAG_P + " ");
                                            for (k = 0; k < gg.annotationsP.size(); k++) {
                                                if (!((String[]) gg.annotationsP.get(k))[0].equals(Domain.ANNOTATION_P))
                                                    Matrix.perror("AGCTFileWriter.class :: annotation != P in P annotations");
                                                f.write(((String[]) gg.annotationsP.get(k))[1]);
                                                if (k < gg.annotationsP.size() - 1)
                                                    f.write(", ");
                                            }
                                            f.write(";\n");
                                        }
                                        if (gg.annotationsC != null) {
                                            f.write(JAnnotationFrame.TAG_C + " ");
                                            for (k = 0; k < gg.annotationsC.size(); k++) {
                                                if (!((String[]) gg.annotationsC.get(k))[0].equals(Domain.ANNOTATION_C))
                                                    Matrix.perror("AGCTFileWriter.class :: annotation != C in C annotations");
                                                f.write(((String[]) gg.annotationsC.get(k))[1]);
                                                if (k < gg.annotationsC.size() - 1)
                                                    f.write(", ");
                                            }
                                            f.write(";\n");
                                        }
                                        f.write("have natural neighbors :\n");
                                    }
                                    f.write("   " + (j + 1) + " -- ");
                                    if (hh.annotationsF != null) {
                                        f.write(JAnnotationFrame.TAG_F + " ");
                                        for (k = 0; k < hh.annotationsF.size(); k++) {
                                            if (!((String[]) hh.annotationsF.get(k))[0].equals(Domain.ANNOTATION_F))
                                                Matrix.perror("AGCTFileWriter.class :: annotation != F in F annotations");
                                            f.write(((String[]) hh.annotationsF.get(k))[1]);
                                            if (k < hh.annotationsF.size() - 1)
                                                f.write(", ");
                                        }
                                        f.write(";\n");
                                    }
                                    if (hh.annotationsP != null) {
                                        f.write(JAnnotationFrame.TAG_P + " ");
                                        for (k = 0; k < hh.annotationsP.size(); k++) {
                                            if (!((String[]) hh.annotationsP.get(k))[0].equals(Domain.ANNOTATION_P))
                                                Matrix.perror("AGCTFileWriter.class :: annotation != P in P annotations");
                                            f.write(((String[]) hh.annotationsP.get(k))[1]);
                                            if (k < hh.annotationsP.size() - 1)
                                                f.write(", ");
                                        }
                                        f.write(";\n");
                                    }
                                    if (hh.annotationsC != null) {
                                        f.write(JAnnotationFrame.TAG_C + " ");
                                        for (k = 0; k < hh.annotationsC.size(); k++) {
                                            if (!((String[]) hh.annotationsC.get(k))[0].equals(Domain.ANNOTATION_C))
                                                Matrix.perror("AGCTFileWriter.class :: annotation != C in C annotations");
                                            f.write(((String[]) hh.annotationsC.get(k))[1]);
                                            if (k < hh.annotationsC.size() - 1)
                                                f.write(", ");
                                        }
                                        f.write(";\n");
                                    }
                                    f.write("\n");
                                    nn++;
                                }
                            }
                            f.write("\n");
                        }
                    }
                }

                f.write(separator());
                f.write("\n\n");
            }

            if (ControlProcess.hasTrue("pcaProcessed")) {
                f.write(separator());
                f.write("** PCA\n\n");

                f.write("Eigenvalues (sorted): ");

                for (i = 0; i < myAGCT.data.getPca_eigenvalues().length; i++)
                    f.write(DF.format(myAGCT.data.getPca_eigenvalues()[i]) + " ");

                f.write(separator());
                f.write("\n\n");
            }

            if ((ControlProcess.hasTrue("softClusteringProcessed")) || (ControlProcess.hasTrue("hardClusteringProcessed"))) {
                int jj;
                for (jj = 0; jj < myAGCT.data.getNbClustering(); jj++) {
                    f.write(separator());
                    f.write("** Clustering #" + jj + "\n\n");

                    f.write("Method: " + myAGCT.getClustering(jj).myClusteringAlgorithm + "\n\n");

                    f.write("Additional data on method: " + myAGCT.getClustering(jj).myClusteringAlgorithm.getAdditionalSavingData() + "\n\n");

                    if ((AGCT.Referenced_Available) && (myAGCT.data.getNbClustering() == 1)) {
                        f.write("Summary of references by clusters:\n");
                        for (i = 0; i < myAGCT.getClustering(jj).myClusteringAlgorithm.getNumberOfClusters(); i++) {
                            f.write("Cluster " + i + ": ");
                            for (j = 0; j < myAGCT.getClustering(jj).myClusteringAlgorithm.getReferencesByCluster().length; j++)
                                f.write(myAGCT.getClustering(jj).myClusteringAlgorithm.getReferencesByCluster()[j][i] + "\t");
                            f.write("\n");
                        }
                        f.write("\n");
                    }

                    f.write("Memberships:\n");
                    for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                        gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
                        f.write(gg.myID + ": ");
                        for (j = 0; j < myAGCT.getClustering(jj).myClusteringAlgorithm.getNumberOfClusters(); j++)
                            f.write(DF.format(gg.getClusterMemberships(jj, j)) + ", ");
                        f.write("\n");
                    }
                    f.write("\n\n");
                    f.write("Composition of clusters (majority rule for memberships):");
                    for (j = 0; j < myAGCT.getClustering(jj).myClusteringAlgorithm.getNumberOfClusters(); j++) {
                        f.write("\n\nCluster " + j);
                        if (Prototype.No_Reduction)
                            f.write(": ");
                        else
                            f.write("(prototypes + their clusters): ");

                        for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                            gid = myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i];
                            gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(gid);

                            better = true;
                            for (k = 0; k < myAGCT.getClustering(jj).myClusteringAlgorithm.getNumberOfClusters(); k++) {
                                if (gg.getClusterMemberships(jj, j) < gg.getClusterMemberships(jj, k))
                                    better = false;
                            }
                            if (better) {
                                f.write(gg.myID + ((gg.asciiName != null) ? "_" + gg.asciiName + "  " : "  "));

                                if (Prototype.No_Reduction)
                                    f.write(" ");
                                else {
                                    f.write(" [ ");
                                    vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
                                    did = (Vector) Prototype.Closest_Center_To_Normalized_Distortions.get(new Integer(gid));
                                    for (iii = 0; iii < vid.size(); iii++) {
                                        nid = ((Integer) vid.get(iii)).intValue();
                                        hh = myAGCT.getDirectGene(nid);
                                        f.write("(" + hh.myID + ((hh.asciiName != null) ? "_" + hh.asciiName + "  " : "  ") + "  :  " + did.get(iii) + ")");
                                    }
                                    f.write(" ]\n");
                                }
                            }
                        }

                        totGenes = 0;
                        f.write("\n\nCluster names " + j + ", preformatted for highlight:\n");
                        for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                            gid = myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i];
                            gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(gid);
                            better = true;
                            for (k = 0; k < myAGCT.getClustering(jj).myClusteringAlgorithm.getNumberOfClusters(); k++) {
                                if (gg.getClusterMemberships(jj, j) < gg.getClusterMemberships(jj, k))
                                    better = false;
                            }
                            if (better) {
                                // if (totGenes > 0)
                                // f.write(", ");
                                dums = gg.myID;
                                f.write(dums.replace(',', '_') + ",C" + j + "\n");
                                totGenes++;
                                /*
                                 * if ( (gg.asciiName != null) &&
								 * (!gg.asciiName.equals("")) ){ if (totGenes >
								 * 0) f.write(", "); f.write (gg.asciiName +
								 * ", "); totGenes++; }
								 */
                                if (!Prototype.No_Reduction) {
                                    vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
                                    for (iii = 1; iii < vid.size(); iii++) {
                                        nid = ((Integer) vid.get(iii)).intValue();
                                        hh = myAGCT.getDirectGene(nid);
                                        // if (totGenes > 0)
                                        // f.write(", ");
                                        dums = hh.myID;
                                        f.write(dums.replace(',', '_') + ",C" + j + "\n");
                                        totGenes++;
										/*
										 * if ( (hh.asciiName != null) &&
										 * (!hh.asciiName.equals("")) ){ if
										 * (totGenes > 0) f.write(", "); f.write
										 * (hh.asciiName + ", "); totGenes++; }
										 */
                                    }
                                }
                            }
                        }

                        f.write("\n\nCluster annotations " + j + ": ");
                        for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                            gid = myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i];
                            gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(gid);

                            better = true;
                            for (k = 0; k < myAGCT.getClustering(jj).myClusteringAlgorithm.getNumberOfClusters(); k++) {
                                if (gg.getClusterMemberships(jj, j) < gg.getClusterMemberships(jj, k))
                                    better = false;
                            }
                            if (better) {
                                if (gg.annotationsF != null)
                                    for (k = 0; k < gg.annotationsF.size(); k++)
                                        f.write("[" + ((String[]) gg.annotationsF.get(k))[0] + "  " + ((String[]) gg.annotationsF.get(k))[1] + "] --- ");
                                if (gg.annotationsP != null)
                                    for (k = 0; k < gg.annotationsP.size(); k++)
                                        f.write("[" + ((String[]) gg.annotationsP.get(k))[0] + "  " + ((String[]) gg.annotationsP.get(k))[1] + "] --- ");
                                if (gg.annotationsC != null)
                                    for (k = 0; k < gg.annotationsC.size(); k++)
                                        f.write("[" + ((String[]) gg.annotationsC.get(k))[0] + "  " + ((String[]) gg.annotationsC.get(k))[1] + "] --- ");

                                if (Prototype.No_Reduction)
                                    f.write(" ");
                                else {
                                    f.write(" Neighbor annotations :: [ ");
                                    vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
                                    did = (Vector) Prototype.Closest_Center_To_Normalized_Distortions.get(new Integer(gid));
                                    for (iii = 0; iii < vid.size(); iii++) {
                                        nid = ((Integer) vid.get(iii)).intValue();
                                        hh = myAGCT.getDirectGene(nid);
                                        if (hh.annotationsF != null)
                                            for (k = 0; k < hh.annotationsF.size(); k++)
                                                f.write("[" + ((String[]) hh.annotationsF.get(k))[0] + "  " + ((String[]) hh.annotationsF.get(k))[1] + "] --- ");
                                        if (hh.annotationsP != null)
                                            for (k = 0; k < hh.annotationsP.size(); k++)
                                                f.write("[" + ((String[]) hh.annotationsP.get(k))[0] + "  " + ((String[]) hh.annotationsP.get(k))[1] + "] --- ");
                                        if (hh.annotationsC != null)
                                            for (k = 0; k < hh.annotationsC.size(); k++)
                                                f.write("[" + ((String[]) hh.annotationsC.get(k))[0] + "  " + ((String[]) hh.annotationsC.get(k))[1] + "] --- ");
                                    }
                                    f.write(" ]\n");
                                }
                            }
                        }
                    }

                    if (jj < myAGCT.data.getNbClustering() - 1)
                        f.write("\n\n");

                    f.write(separator());
                    f.write("\n\n");
                }
            }

            if ((Statistics.allTests != null) && (Statistics.allTests != "")) {
                int jj;
                f.write(separator());
                f.write("** Chi-square tests performed on clusterings\n\n");
                for (jj = 0; jj < Statistics.allTests.length(); jj++)
                    f.write(Statistics.allTests.charAt(jj));

                f.write("Complete and ordered list of all annotation tags, given as: Type Tag Clustering Cluster P-value P*-value\n");

                if (Statistics.annotationsBH != null)
                    for (i = 0; i < Statistics.annotationsBH.size(); i++) {
                        littleVector = (Vector) Statistics.annotationsBH.get(i);
                        f.write(littleVector.get(0) + "\t" + littleVector.get(1) + "\t" + littleVector.get(2) + "\t" + littleVector.get(3) + "\t" + littleVector.get(4) + "\t"
                                + littleVector.get(5) + "\n");
                    }

                f.write(separator());
                f.write("\n\n");
            }

            f.write(separator());
            f.write("** Processed genes list\n\n");

            for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
                f.write(gg.toString());
            }

            f.write(separator());
            f.write("\n\n");

            f.close();
        } catch (IOException e) {
            Matrix.perror("AGCTFileWriter.class :: Writing error");
        }
    }

    public void toSavingChi2(final File rf) {
        Thread t = new Thread() {
            public void run() {
                FileWriter f = null;
                try {
                    int i, j, k, nc, ncur;
                    Vector littleVector;
                    Vector cumulative = new Vector(), element, vTable;
                    double currentValue, bestValue = -1, pcur, ad, bc;
                    String refann, curann, signe;
                    f = new FileWriter(rf);

                    f.write(separator());
                    f.write("Saving CHI2 Results for domain " + myAGCT.data.getMyDomain().domainName);
                    f.write(separator());

                    f.write("\n\n");

                    f.write(separator());
                    f.write("Processing parameters:\n");
                    f.write(myAGCT.getParameterString());
                    f.write(separator());

                    f.write("** Chi-square tests performed on clusterings\n\n");
                    for (i = 0; i < Statistics.allTests.length(); i++)
                        f.write(Statistics.allTests.charAt(i));

                    f.write("Complete and ordered list of all annotation tags, given as: Type Tag Clustering Cluster P-value P*-value [Chi2 table]\n");

                    if (Statistics.annotationsBH != null) {
                        AGCTCounter cc = new AGCTCounter("Saving Chi2 Stats", 2 * Statistics.annotationsBH.size());
                        for (i = 0; i < Statistics.annotationsBH.size(); i++) {
                            littleVector = (Vector) Statistics.annotationsBH.get(i);
                            f.write(littleVector.get(0) + "\t" + littleVector.get(1) + "\t" + littleVector.get(2) + "\t" + littleVector.get(3) + "\t" + littleVector.get(4) + "\t"
                                    + littleVector.get(6) + "\t");
                            vTable = (Vector) littleVector.get(5);
                            f.write("[" + vTable.get(0) + ", " + vTable.get(1) + ", " + vTable.get(2) + ", " + vTable.get(3) + "]\n");
                            cc.increment();
                        }
                        f.write("\n\n\nSet of tags for each cluster, having P*-Values <= Satistics.LIMIT_P (" + Statistics.LIMIT_P_CHI2 + ")\n");

                        nc = AGCT.getInstance().getClustering(Statistics.refNumClustering).myClusteringAlgorithm.getNumberOfClusters();
                        for (k = 0; k < nc; k++) {
                            f.write("Cluster " + k + ":\n\n");
                            for (j = 0; j < 3; j++) {
                                refann = null;
                                if (j == 0)
                                    refann = Domain.ANNOTATION_F;
                                else if (j == 1)
                                    refann = Domain.ANNOTATION_P;
                                else if (j == 2)
                                    refann = Domain.ANNOTATION_C;
                                f.write("Annotation " + refann + "\n");
                                for (i = 0; i < Statistics.annotationsBH.size(); i++) {
                                    littleVector = (Vector) Statistics.annotationsBH.get(i);
                                    ncur = ((Integer) littleVector.get(3)).intValue();
                                    curann = (String) littleVector.get(0);
                                    pcur = ((Double) littleVector.get(6)).doubleValue();
                                    vTable = (Vector) littleVector.get(5);
                                    ad = ((Double) vTable.get(0)).doubleValue() * ((Double) vTable.get(3)).doubleValue();
                                    bc = ((Double) vTable.get(1)).doubleValue() * ((Double) vTable.get(2)).doubleValue();

                                    if ((ncur == k) && (curann.equals(refann)) && (pcur < Statistics.LIMIT_P_CHI2)) {
                                        if (ad > bc)
                                            f.write("+ ");
                                        else if (ad < bc)
                                            f.write("- ");
                                        else
                                            f.write("/ ");

                                        f.write((String) littleVector.get(1) + "\t" + littleVector.get(6) + "\n");
                                    }
                                }
                                f.write("\n");
                            }
                            f.write("\n");
                        }

                        f.write("Cumulative plot (x,y), in which x=P*-value and y=%tags having P*-value <= x\n");

                        for (i = 0; i < Statistics.annotationsBH.size(); i++) {
                            littleVector = (Vector) Statistics.annotationsBH.get(i);
                            currentValue = ((Double) littleVector.get(6)).doubleValue();
                            if (i == 0) {
                                bestValue = currentValue;
                                element = new Vector();
                                element.addElement(new Double(currentValue));
                                element.addElement(new Double(i + 1));
                                cumulative.addElement(element);
                            } else {
                                if (currentValue == bestValue) {
                                    element = (Vector) cumulative.get(cumulative.size() - 1);
                                    element.setElementAt(new Double(i + 1), 1);
                                } else {
                                    if (currentValue < bestValue)
                                        Matrix.perror("AGCTFileWriter.class :: bad ordering of values");
                                    bestValue = currentValue;
                                    element = new Vector();
                                    element.addElement(new Double(currentValue));
                                    element.addElement(new Double(i + 1));
                                    cumulative.addElement(element);
                                }
                            }
                            cc.increment();
                        }
                        cc.end();

                        currentValue = ((Double) ((Vector) cumulative.get(cumulative.size() - 1)).get(1)).doubleValue();

                        for (i = 0; i < cumulative.size(); i++) {
                            littleVector = (Vector) cumulative.get(i);
                            f.write(littleVector.get(0) + "\t" + (((Double) littleVector.get(1)).doubleValue() / currentValue) + "\n");
                        }
                    }
                    f.write(separator());

                    //

                    f.close();
                } catch (IOException e) {
                    Matrix.perror("AGCTFileWriter.class :: Writing error for Chi2 statistics");
                }
            }
        };
        t.start();
    }

    public void toSavingClustering(final File rf, final int nc) {
        Debug.debug("toSavingClustering");
        FileWriter f;
        Clustering cc = myAGCT.getClustering(nc);
        int i;
        Gene gg;
        try {
            f = new FileWriter(rf);
            f.write(AGCTFileWriter.Clustering_Born_String + Scenario.myToken + cc.bornString + "\n");
            f.write(AGCTFileWriter.DATA_Data + Scenario.myToken + "\n");

            AGCTCounter ccc = new AGCTCounter("Saving clustering", myAGCT.data.getMyDomain().numberSelectedGenes);
            Debug.debug("numberSelectedGenes", myAGCT.data.getMyDomain().numberSelectedGenes);
            for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                gg = myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
                f.write(i + Scenario.myToken +
//						gg.myID+Scenario.myToken+cc.myClusteringAlgorithm.majorityCluster(i)
                        gg.getClusterNumber()
                        + "\n");
                ccc.increment();
            }
            ccc.end();

            f.close();
        } catch (IOException e) {
            Matrix.perror("AGCTFileWriter.class :: Writing error for saving Clustering");
        }
    }

    public void toSavingClustering(final FileWriter f, final int nc) {
        Debug.debug("toSavingClustering");
        //		Thread t=new Thread(){
        //			public void run(){
//		FileWriter f=null;
        Clustering cc = myAGCT.getClustering(nc);
        int i;
        Gene gg;
        try {
//			f=new FileWriter(rf);
            f.write(AGCTFileWriter.Clustering_Born_String + Scenario.myToken + cc.bornString + "\n");
            f.write(AGCTFileWriter.DATA_Data + Scenario.myToken + "\n");

            AGCTCounter ccc = new AGCTCounter("Saving clustering", myAGCT.data.getMyDomain().numberSelectedGenes);
            Debug.debug("numberSelectedGenes", myAGCT.data.getMyDomain().numberSelectedGenes);
            for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
                f.write(i + Scenario.myToken +
//						gg.myID+Scenario.myToken+cc.myClusteringAlgorithm.majorityCluster(i)
                        gg.getClusterNumber()
                        + "\n");
                ccc.increment();
            }
            ccc.end();

            f.close();
        } catch (IOException e) {
            Matrix.perror("AGCTFileWriter.class :: Writing error for saving Clustering");
        }
        //			}
        //		};
        //		t.start();
    }

    public static void loadClustering(File rf, AGCT ap) {
        FileReader e;
        BufferedReader br;
        String st, sc, type = "", name;
        int index, i;
        StringTokenizer t;
        Vector allPoints = new Vector(), element;
        boolean collectingPoints = false;
        // each element = Vector (index, gene name, cluster)
        boolean okToRegister = true;
        Gene gg;

        try {
            e = new FileReader(rf);
        } catch (FileNotFoundException ex) {
            JInformationFrame.getInstance().setText("The clustering you try to load does not exist");
            return;
        }

        br = new BufferedReader(e);
        try {
            while ((st = br.readLine()) != null) {
                if ((st.length() > 1) && (!st.substring(0, AGCTFileWriter.DATA_Comments.length()).equals(AGCTFileWriter.DATA_Comments))) {
                    t = new StringTokenizer(st, Scenario.myToken);
                    if (!collectingPoints) {
                        sc = t.nextToken();
                        if (sc.equals(AGCTFileWriter.DATA_Data))
                            collectingPoints = true;
                        else if (sc.equals(AGCTFileWriter.Clustering_Born_String))
                            type = t.nextToken();
                    } else {
                        if (type.equals(new String("")))
                            okToRegister = false;
                        else if (t.countTokens() == 3) {
                            element = new Vector();
                            element.add(new Integer(Integer.parseInt(t.nextToken())));
                            element.add(new String(t.nextToken()));
                            element.add(new Integer(Integer.parseInt(t.nextToken())));
                            allPoints.addElement(element);
                        } else if (t.countTokens() == 2) {
                            element = new Vector();
                            int id = Integer.parseInt(t.nextToken());
                            element.add(new Integer(id));
                            element.add(Domain.getInstance().getSelectedGenes().get(id).getAnnotation(Gene.ID));
                            element.add(new Integer(Integer.parseInt(t.nextToken())));
                            allPoints.addElement(element);
                        }
                    }
                }
            }
            e.close();
        } catch (IOException ex) {
            JInformationFrame.getInstance().setText("IOError when loading clustering");
            return;
        }
        type = JAGCTClusteringPane.defaultStrings[JAGCTClusteringPane.KM];
        int max = 0;
        for (i = 0; i < allPoints.size(); i++) {
            max = Math.max(max, (Integer) (((Vector) (allPoints.get(i))).get(2)));
        }
        System.err.println("ncluster : " + (max + 1));
        type += (max + 1);

        if (!okToRegister)
            System.out.println("This clustering is not applicable (not found type and number of clusters)");

        if (Clustering.getIndex(type.substring(0, 2)) != 1) {
            System.out.println("Not a K-means algorithm");
//			okToRegister=false;
        }

        if ((okToRegister) && (allPoints.size() != ap.data.getMyDomain().numberSelectedGenes)) {
            System.out.println("This clustering does not contain the same number of points (" + ap.data.getMyDomain().numberSelectedGenes + " != " + allPoints.size());
            okToRegister = false;
        }
        if (okToRegister)
            for (i = 0; i < allPoints.size(); i++) {
                element = (Vector) allPoints.get(i);
                name = (String) element.get(1);
                index = ((Integer) element.get(0)).intValue();
                gg = (Gene) ap.data.getMyDomain().getGenes().get(ap.data.getMyDomain().selectedGeneNumberToGeneNumber[index]);
                if (!name.equals(gg.myID)) {
                    System.out.println("Misatch in gene #" + index + " name for clustering: " + name + " != " + gg.myID);
                    okToRegister = false;
                }
            }
        if (okToRegister) {
            ap.data.setNbClustering(ap.data.getNbClustering() + 1);
            Clustering newClustering;
            newClustering = new Clustering(ap, ap.data.getNbClustering() - 1);
            newClustering.getOptions(type, true);
            ((AGCTClustering_KM) newClustering.myClusteringAlgorithm).toClusteringLite(allPoints);

            if (ap.data.getNbClustering() == 1) {
                ap.data.setAllClusterings(new Vector());
                ap.myTabbedPane.myPCAPane.currentClustering = 0;
                ap.myTabbedPane.myManifoldPane.currentClustering = 0;
                ap.myTabbedPane.myClusteringPane.activateSearch();
            }

            ap.data.getAllClusterings().add(newClustering);
            ap.myClusteringProfileFrame.addClusteringLF(newClustering, ap.data.getNbClustering() - 1);

            ap.updateClusteringLF(JAGCTVisualizationPane.stringClustering(ap.data.getNbClustering() - 1));
            ap.myTabbedPane.myClusteringPane.listClusteringPane.setText(ap.myTabbedPane.myClusteringPane.clusteringList());

            newClustering.generate_VV();
            ap.myTabbedPane.myManifoldPane.setMembershipsButtons(true);
            ap.myTabbedPane.myPCAPane.setMembershipsButtons(true);
            ap.myTabbedPane.myCorrelationPane.setMembershipsButtons(true);

            newClustering.myClusteringAlgorithm.fillGeneMemberships();
        }

        e = null;
    }

    public void toSavingDelaunay(final File rf) {
        Thread t = new Thread() {
            public void run() {

                FileWriter f = null;
                try {
                    int i, j, k, nid;
                    double nangle;
                    Gene gg, ggn;
                    boolean okNeighbor, okRef;
                    Vector nn = null, ll, element;
                    Character pos = new Character('+'), neg = new Character('-'), dumchar = null;

                    f = new FileWriter(rf);

                    f.write(separator());
                    f.write("Saving the *Filtered* Delaunay Triangulation (P-Del = " + Statistics.getLIMIT_P_DELAUNAY() + ") for domain "
                            + myAGCT.data.getMyDomain().domainName);
                    f.write(separator());

                    f.write("\n\n");

                    f.write(separator());
                    f.write("Processing parameters:\n");
                    f.write(myAGCT.getParameterString());
                    f.write(separator());

                    if (myAGCT.data.isTriangulated()) {
                        f.write("[G01] Unordered (Complete) graph, given as: Gene-Name_ID List-of-Delaunay-Neighbors (Gene-Name_IDS)\n");

                        AGCTCounter cc = new AGCTCounter("Saving Complete Filtered Del. Tr.", myAGCT.data.getMyDomain().numberSelectedGenes);
                        for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                            gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);

                            if ((gg.neighbors != null) && (gg.neighbors.size() > 0)) {
                                nn = new Vector();
                                for (j = 0; j < gg.neighbors.size(); j++) {
                                    okNeighbor = false;
                                    nid = ((Integer) gg.neighbors.get(j)).intValue();
                                    ggn = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[nid]);

                                    if (Statistics.isPositiveDelauney(gg, j) || Statistics.isNegativeDelauney(gg, j))
                                        okNeighbor = true;

                                    if (Statistics.isPositiveDelauney(gg, j))
                                        dumchar = pos;
                                    else if (Statistics.isNegativeDelauney(gg, j))
                                        dumchar = neg;

                                    if (okNeighbor) {
                                        element = new Vector();
                                        element.addElement(geneStringDelaunay(ggn, true) + "_(*)");
                                        element.addElement(dumchar);

                                        nn.addElement(element);
                                        if (!Prototype.No_Reduction) {
                                            ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[nid]));
                                            for (k = 1; k < ll.size(); k++) {
                                                Gene gene = (Gene) myAGCT.data.getMyDomain().getGenes().get(((Integer) ll.get(k)).intValue());
                                                element = new Vector();
                                                element.addElement(geneStringDelaunay(gene, false));
                                                element.addElement(dumchar);

                                                nn.addElement(element);
                                            }
                                        }
                                    }
                                }

                                if (nn.size() > 0) {
                                    f.write(geneStringDelaunay(gg, true) + "\t");
                                    for (k = 0; k < nn.size(); k++) {
                                        element = (Vector) nn.get(k);
                                        f.write((String) element.get(0) + " (" + (Character) element.get(1) + ")");
                                        if (k < nn.size() - 1)
                                            f.write("\t");
                                        else
                                            f.write("\n");
                                    }

                                    if (!Prototype.No_Reduction) {
                                        ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]));
                                        for (j = 1; j < ll.size(); j++) {
                                            ggn = (Gene) myAGCT.data.getMyDomain().getGenes().get(((Integer) ll.get(j)).intValue());
                                            f.write(geneStringDelaunay(ggn, false) + "\t");
                                            for (k = 0; k < nn.size(); k++) {
                                                element = (Vector) nn.get(k);
                                                f.write((String) element.get(0) + " (" + (Character) element.get(1) + ")");
                                                if (k < nn.size() - 1)
                                                    f.write("\t");
                                                else
                                                    f.write("\n");
                                            }
                                        }
                                    }
                                }
                            }
                            cc.increment();
                        }
                        cc.end();
                        f.write(separator());
                        f.write("\n\n[G02] Unordered (Prototype + Highlighted) graph, given as: Gene-Name_ID List-of-Delaunay-Neighbors (Gene-Name_IDS)\n");

                        cc = new AGCTCounter("Saving P + H Filtered Del. Tr.", myAGCT.data.getMyDomain().numberSelectedGenes);
                        for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                            gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);

                            if ((gg.neighbors != null) && (gg.neighbors.size() > 0)) {
                                nn = new Vector();
                                for (j = 0; j < gg.neighbors.size(); j++) {
                                    okNeighbor = false;
                                    nid = ((Integer) gg.neighbors.get(j)).intValue();
                                    ggn = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[nid]);
                                    nangle = ((Double) gg.neighbor_angles.get(j)).doubleValue();

                                    if (Statistics.isPositiveDelauney(gg, j) || Statistics.isNegativeDelauney(gg, j))
                                        okNeighbor = true;

                                    if (Statistics.isPositiveDelauney(gg, j))
                                        dumchar = pos;
                                    else if (Statistics.isNegativeDelauney(gg, j))
                                        dumchar = neg;

                                    if (okNeighbor) {
                                        element = new Vector();
                                        element.addElement(geneStringDelaunay(ggn, true));
                                        element.addElement(dumchar);
                                        nn.addElement(element);
                                        if (!Prototype.No_Reduction) {
                                            ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[nid]));
                                            for (k = 1; k < ll.size(); k++) {
                                                Gene gene = (Gene) myAGCT.data.getMyDomain().getGenes().get(((Integer) ll.get(k)).intValue());
                                                if ((AGCT.Referenced_Available) && (gene.isReferenced())) {
                                                    element = new Vector();
                                                    element.addElement(geneStringDelaunay(gene, false));
                                                    element.addElement(dumchar);
                                                    nn.addElement(element);
                                                }
                                            }
                                        }
                                    }
                                }

                                if (nn.size() > 0) {
                                    f.write(geneStringDelaunay(gg, true) + "\t");
                                    for (k = 0; k < nn.size(); k++) {
                                        element = (Vector) nn.get(k);
                                        f.write((String) element.get(0) + " (" + (Character) element.get(1) + ")");
                                        if (k < nn.size() - 1)
                                            f.write("\t");
                                        else
                                            f.write("\n");
                                    }

                                    if (!Prototype.No_Reduction) {
                                        ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]));
                                        for (j = 1; j < ll.size(); j++) {
                                            ggn = (Gene) myAGCT.data.getMyDomain().getGenes().get(((Integer) ll.get(j)).intValue());
                                            if ((AGCT.Referenced_Available) && (ggn.isReferenced())) {
                                                f.write(geneStringDelaunay(ggn, false) + "\t");
                                                for (k = 0; k < nn.size(); k++) {
                                                    element = (Vector) nn.get(k);
                                                    f.write((String) element.get(0) + " (" + (Character) element.get(1) + ")");
                                                    if (k < nn.size() - 1)
                                                        f.write("\t");
                                                    else
                                                        f.write("\n");
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            cc.increment();
                        }
                        cc.end();
                        f.write(separator());
                        if (AGCT.Referenced_Available) {
                            f.write("\n\n[G03] Unordered (Highlighted) graph, given as: Gene-Name_ID List-of-Delaunay-Neighbors (Gene-Name_IDS)\n");

                            cc = new AGCTCounter("Saving P + H Filtered Del. Tr.", myAGCT.data.getMyDomain().numberSelectedGenes);
                            for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                                gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);

                                if ((gg.neighbors != null) && (gg.neighbors.size() > 0)) {
                                    nn = new Vector();
                                    for (j = 0; j < gg.neighbors.size(); j++) {
                                        okNeighbor = false;
                                        nid = ((Integer) gg.neighbors.get(j)).intValue();
                                        ggn = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[nid]);
                                        if (Statistics.isPositiveDelauney(gg, j) || Statistics.isNegativeDelauney(gg, j))
                                            okNeighbor = true;

                                        if (Statistics.isPositiveDelauney(gg, j))
                                            dumchar = pos;
                                        else if (Statistics.isNegativeDelauney(gg, j))
                                            dumchar = neg;

                                        if (okNeighbor) {
                                            element = new Vector();
                                            if (ggn.isReferenced()) {
                                                element.addElement(geneStringDelaunay(ggn, true));
                                                element.addElement(dumchar);
                                                nn.addElement(element);
                                            }
                                            if (!Prototype.No_Reduction) {
                                                ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[nid]));
                                                for (k = 1; k < ll.size(); k++) {
                                                    Gene gene = (Gene) myAGCT.data.getMyDomain().getGenes().get(((Integer) ll.get(k)).intValue());
                                                    if ((AGCT.Referenced_Available) && (gene.isReferenced())) {
                                                        element = new Vector();
                                                        element.addElement(geneStringDelaunay(gene, false));
                                                        element.addElement(dumchar);
                                                        nn.addElement(element);
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    if (nn.size() > 0) {
                                        if (gg.isReferenced()) {
                                            f.write(geneStringDelaunay(gg, true) + "\t");
                                            for (k = 0; k < nn.size(); k++) {
                                                element = (Vector) nn.get(k);
                                                f.write((String) element.get(0) + " (" + (Character) element.get(1) + ")");
                                                if (k < nn.size() - 1)
                                                    f.write("\t");
                                                else
                                                    f.write("\n");
                                            }
                                        }

                                        if (!Prototype.No_Reduction) {
                                            ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]));
                                            for (j = 1; j < ll.size(); j++) {
                                                ggn = (Gene) myAGCT.data.getMyDomain().getGenes().get(((Integer) ll.get(j)).intValue());
                                                if ((AGCT.Referenced_Available) && (ggn.isReferenced())) {
                                                    f.write(geneStringDelaunay(ggn, false) + "\t");
                                                    for (k = 0; k < nn.size(); k++) {
                                                        element = (Vector) nn.get(k);
                                                        f.write((String) element.get(0) + " (" + (Character) element.get(1) + ")");
                                                        if (k < nn.size() - 1)
                                                            f.write("\t");
                                                        else
                                                            f.write("\n");
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                cc.increment();
                            }
                            cc.end();
                            f.write(separator());
                        }
                    }
                    f.close();
                } catch (IOException e) {
                    Matrix.perror("AGCTFileWriter.class :: Writing error for Filtered Delaunay Triangulation");
                }
            }
        };
        t.start();
    }

    private void write_params(FileWriter f) {

    }

    /**
     * clustering とtriangulationが既にすんでいるときに呼ぶ。 triangulation
     * threshold,totaldistance thresholdをともに くぐり抜けた奴らを、属するクラスターごとに表示する。
     *
     * @param rf
     */
    public void toSavingCluseringAndTriangulation(final File rf) {
        ArrayList<Gene> canUse = new ArrayList<Gene>();
        HashSet<Gene> canUseSet = new HashSet<Gene>();
        // double TDthreshold = 0.2;
        ArrayList<Double> tmp = new ArrayList<Double>();
        for (int j = 0; j < myAGCT.data.getMyDomain().numberSelectedGenes; j++) {
            Gene gene = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[j]);
            tmp.add(gene.totalDistance);
        }
        Collections.sort(tmp);
        double thr = tmp.get((int) Math.round((tmp.size() - 1) * (1 - Statistics.LIMIT_P_MAGNITUDE)));
        for (int j = 0; j < myAGCT.data.getMyDomain().numberSelectedGenes; j++) {
            Gene gene = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[j]);
            if (gene.totalDistance >= thr) {
                boolean ok = false;
                for (int k = 0; k < gene.neighbors.size(); k++) {
                    if (Statistics.isPositiveDelauney(gene, k) || Statistics.isNegativeDelauney(gene, k))
                        ok = true;
                }
                if (ok) {
                    canUse.add(gene);
                    canUseSet.add(gene);
                }
            }
        }

        FileWriter f;
        try {
            f = new FileWriter(rf);
            f.write(String.format("LIMIT_P_DELAUNAY = %.2f\t", Statistics.getLIMIT_P_DELAUNAY()));
            f.write(String.format("LIMIT_P_MAGNITUDE = %.2f\n", Statistics.LIMIT_P_MAGNITUDE));
            int clusterCount = myAGCT.getClustering(0).myClusteringAlgorithm.getNumberOfClusters();
            for (int i = 0; i < clusterCount; i++) {
                f.write(String.format("Cluster%d", i));
                ArrayList<Gene> list = new ArrayList<Gene>();
                for (Gene gene : canUse) {
                    if (gene.getClusterMemberships(0, i) == 1) {
                        list.add(gene);
                    }
                }
                f.write("(" + list.size() + ")\n");
                for (Gene gene : list) {
                    // gene の情報を表示
                    f.write(gene.myID + '\t');
                    f.write(gene.asciiName + "\t");
                    f.write(gene.geneTitle + "\t");
                    f.write(gene.refSeq + "\t");
                    if (gene.annotationsP == null || gene.annotationsP.isEmpty())
                        f.write("null" + "\t");
                    else
                        f.write(((String[]) gene.annotationsP.get(0))[1] + "\t");
                    if (gene.annotationsC == null || gene.annotationsC.isEmpty())
                        f.write("null" + "\t");
                    else
                        f.write(((String[]) gene.annotationsC.get(0))[1] + "\t");
                    if (gene.annotationsF == null || gene.annotationsF.isEmpty())
                        f.write("null" + "\n");
                    else
                        f.write(((String[]) gene.annotationsF.get(0))[1] + "\n");
                }
            }
            f.write("binary relations (positive 0 <-> PI negative)\n");//Above
            for (Gene gene : canUse) {
                int cl = -1;
                for (int i = 0; i < clusterCount; i++) {
                    if (gene.getClusterMemberships(0, i) == 1)
                        cl = i;
                }
                // geneのneighborsを表示
                ArrayList<Gene> neis = new ArrayList<Gene>();
                // ArrayList<String> dums = new ArrayList<String>();
                ArrayList<Double> angles = new ArrayList<Double>();
                ArrayList<Integer> cls = new ArrayList<Integer>();
                for (int j = 0; j < gene.neighbors.size(); j++) {
                    int id = (Integer) gene.neighbors.get(j);
                    //					if(id<=gene.getIdNumber())continue;
                    Gene nei = myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[id]);
                    if (!canUseSet.contains(nei) || nei.getIdNumber() < gene.getIdNumber())
                        continue;
                    double nangle = (Double) gene.neighbor_angles.get(j);
                    int c = -1;
                    for (int k = 0; k < clusterCount; k++) {
                        if (nei.getClusterMemberships(0, k) == 1)
                            c = k;
                    }
                    // if (c == i)
                    // c = -1;
                    neis.add(nei);
                    // dums.add(dumchar);
                    angles.add(nangle);
                    cls.add(c);
                    // }
                }

                if (neis.size() > 0) {
                    for (int k = 0; k < neis.size(); k++) {
                        f.write(geneStringDelaunay(gene, true) + " (c" + cl + ")\t");
                        Gene nei = neis.get(k);
                        // String dum = dums.get(k);
                        int c = cls.get(k);
                        f.write(geneStringDelaunay(nei, false) + " (c" + c + ") (" + String.format("%.2f", angles.get(k)) + ")");

                        if (Statistics.isPositiveDelauney(gene, k) || Statistics.isNegativeDelauney(gene, k))
                            f.write("\tabove threshold");
                        else
                            f.write("\tunder threshold");
                        f.write("\n");

                    }
                }
            }
            f.close();
        } catch (IOException e) {
            Matrix.perror("toSavingClusterAndTriangulation : writeError");
        }
    }

    public void toSavingWeb(final File rf) {
        Thread t = new Thread() {
            public void run() {

                FileWriter f = null;
                try {
                    int i, j, d;
                    Gene gg, ggn;
                    Vector ll;
                    double[] compo;

                    f = new FileWriter(rf);
                    AGCTCounter cc = new AGCTCounter("Saving to Webfile", myAGCT.data.getMyDomain().numberSelectedGenes);
                    for (i = 0; i < myAGCT.data.getMyDomain().numberSelectedGenes; i++) {
                        gg = (Gene) myAGCT.data.getMyDomain().getGenes().get(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
                        if (myAGCT.myTabbedPane.getSelectedIndex() == 1)
                            compo = gg.manifold_Pnt3D.coordinates;
                        else
                            compo = gg.pca_Pnt3D.coordinates;

                        f.write(gg.myID + "\t");
                        f.write(gg.asciiName + "\t");
                        for (d = 0; d < 3; d++) {
                            f.write(compo[d] + "");
                            if (d < 2)
                                f.write("\t");
                        }
                        f.write("\n");

                        if (!Prototype.No_Reduction) {
                            ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.data.getMyDomain().selectedGeneNumberToGeneNumber[i]));
                            for (j = 1; j < ll.size(); j++) {
                                ggn = (Gene) myAGCT.data.getMyDomain().getGenes().get(((Integer) ll.get(j)).intValue());
                                f.write(ggn.myID + "\t");
                                f.write(ggn.asciiName + "\t");
                                for (d = 0; d < 3; d++) {
                                    f.write(compo[d] + "");
                                    if (d < 2)
                                        f.write("\t");
                                }
                                f.write("\n");
                            }
                        }
                        cc.increment();
                    }
                    cc.end();
                    f.close();
                } catch (IOException e) {
                    Matrix.perror("AGCTFileWriter.class :: Writing error for Web file");
                }
            }
        };
        t.start();
    }

    public static String geneStringDelaunay(Gene ggn, boolean proto) {
        String dumString = "";
        String refGene = "_(R)";
        if (ggn.asciiName != null)
            dumString += ggn.asciiName + "_";
        dumString += ggn.myID;
        if (proto)
            dumString += "_(*)";
        if ((AGCT.Referenced_Available) && (ggn.isReferenced()))
            dumString += refGene;
        return dumString;
    }

    public static String separator() {
        return "\n----------------------------------------------------------------------------\n";
    }
}