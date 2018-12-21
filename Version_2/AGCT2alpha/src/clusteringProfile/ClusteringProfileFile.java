package clusteringProfile;

import forDebug.Debug;
import gene.Gene;
import gene.GeneList;
import ligand.Ligand;
import test.AGCT;
import test.AGCTFileWriter;
import test.Matrix;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

public class ClusteringProfileFile extends AGCTFileWriter {
    private ClusteringProfileFile() {
    }

    public static void toSavingClusterProfile(final File rf) {
        Debug.debug("ClusteringProfileFile . toSavingClusterProfile", rf);
        Thread t = new Thread() {
            public void run() {
                FileWriter f;
                try {
                    clusteringProfile.JClusteringProfileGridPlot graphPanel = AGCT.getInstance().myClusteringProfileFrame.graphPanel;

                    f = new FileWriter(rf);

                    f.write(separator());
                    f.write("Saving Cluster Profiles for domain "
                            + AGCT.getInstance().data.getMyDomain().domainName);
                    f.write(separator());

                    f.write("\n\n");

                    f.write(separator());
                    f.write("Processing parameters:\n");
                    f.write(AGCT.getInstance().getParameterString());
                    f.write(separator());

                    f
                            .write("*** All average time series for each cluster, given as Time_Stamp Average Q25 Q50 Q75\n\n\n");

                    f.write(separator());
                    ClusteringProfile clusteringProfile = null;
                    for (int i = 0; i < graphPanel.getNumberOfClusters(); i++) {
                        f.write("** Cluster" + i + ":\n\n");
                        int index = 0;
                        for (Ligand ligand : AGCT.getInstance().data.getMyDomain()
                                .getLigands()) {
                            if (ligand.isChecked()) {
                                clusteringProfile = graphPanel.all_stamps[i][index]
                                        .getClusteringProfile();
                                f.write(ligand.getName() + ":\n");
                                for (int timeId = 0; timeId < clusteringProfile.average_profiles.length; timeId++) {
                                    if (!clusteringProfile.undetermined_profiles[timeId])
                                        f.write(clusteringProfile.time(timeId) + "\t"
                                                + clusteringProfile.average_profiles[timeId] + "\t"
                                                + clusteringProfile.q25_profiles[timeId] + "\t"
                                                + clusteringProfile.q50_profiles[timeId] + "\t"
                                                + clusteringProfile.q75_profiles[timeId] + "\n");
                                }
                                index++;
                                f.write("\n");
                            }
                        }
                        f.write("\n");
                    }

                    f.write(separator());

                    if (clusteringProfile != null) {

                        int experimentId = clusteringProfile.experimentId;

                        class Pair implements Comparable<Pair> {
                            Gene gene;
                            double val;

                            Pair(Gene gene, double val) {
                                this.gene = gene;
                                this.val = val;
                            }

                            public int compareTo(Pair o) {
                                return ((Double) val).compareTo(o.val);
                            }
                        }

                        int CHOOSE = 100;
                        f.write("first " + CHOOSE + " genes near average of cluster\n\n");
                        GeneList[] genes = new GeneList[graphPanel.getNumberOfClusters()];
                        for (int i = 0; i < genes.length; i++) {
                            genes[i] = new GeneList();
                        }
                        for (Gene gene : clusteringProfile.getAllGenes()) {
                            int clusterId = gene.getClusterNumber(experimentId);
                            genes[clusterId].add(gene);
                        }
                        for (int i = 0; i < graphPanel.getNumberOfClusters(); i++) {
                            f.write("** Cluster" + i + ":\n\n");
                            int index = 0;
                            for (int ligandId = 0; ligandId < AGCT.getInstance().data.getMyDomain().getLigands().size(); ligandId++) {
                                Ligand ligand = AGCT.getInstance().data.getMyDomain().getLigands().get(ligandId);
                                if (ligand.isChecked()) {
                                    clusteringProfile = graphPanel.all_stamps[i][index]
                                            .getClusteringProfile();
                                    f.write(ligand.getName() + ":\n");   // __Test_570.Cell1:

                                    ArrayList<Pair> list = new ArrayList<Pair>();
                                    for (Gene gene : genes[i]) {
                                        double val = 0;
                                        for (int timeId = 0; timeId < clusteringProfile.average_profiles.length; timeId++) {
                                            if (!clusteringProfile.undetermined_profiles[timeId])
                                                val += Math.pow(clusteringProfile.average_profiles[timeId] - gene.getRawCoordinate(ligandId)[timeId], 2);
                                        }
                                        list.add(new Pair(gene, val));
                                    }
                                    Collections.sort(list);

                                    for (int j = 0; j < Math.min(CHOOSE, list.size()); j++) {
                                        f.write(list.get(j).gene + "\n");
                                    }
                                    f.write("\n");
                                }
                            }
                            f.write("\n");
                        }
                    }
                    f.flush();
                    f.close();
                } catch (IOException e) {
                    Matrix
                            .perror("AGCTFileWriter.class :: Writing error for Cluster profiles");
                }
            }
        };
        t.start();
    }
}
