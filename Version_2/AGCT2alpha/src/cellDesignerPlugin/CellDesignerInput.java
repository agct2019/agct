package cellDesignerPlugin;

import gene.Gene;
import test.AGCT;
import test.Matrix;
import test.Statistics;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

/*
 * Annotationがなければ、Probe名を入れる
 */
public class CellDesignerInput {
    public static void save(final File rf) {
        ArrayList<Gene> canUse = new ArrayList<Gene>();
        HashSet<Gene> canUseSet = new HashSet<Gene>();
        ArrayList<Double> tmp = new ArrayList<Double>();
        for (int j = 0; j < AGCT.getInstance().data.getMyDomain().numberSelectedGenes; j++) {
            Gene gene = (Gene) AGCT.getInstance().data
                    .getMyDomain()
                    .getGenes()
                    .get(
                            AGCT.getInstance().data.getMyDomain().selectedGeneNumberToGeneNumber[j]);
            tmp.add(gene.totalDistance);
        }
        Collections.sort(tmp);
        double thr = tmp.get((int) Math.round((tmp.size() - 1)
                * (1 - Statistics.LIMIT_P_MAGNITUDE)));
        for (int j = 0; j < AGCT.getInstance().data.getMyDomain().numberSelectedGenes; j++) {
            Gene gene = AGCT.getInstance().data
                    .getMyDomain()
                    .getGenes()
                    .get(
                            AGCT.getInstance().data.getMyDomain().selectedGeneNumberToGeneNumber[j]);
            if (gene.isVisible()) {
                boolean ok = false;
                for (int k = 0; k < gene.neighbors.size(); k++) {
                    if (Statistics.isNegativeDelauney(gene, k) || Statistics.isPositiveDelauney(gene, k))
                        ok = true;
                }
                if (ok && gene.isVisible()) {
                    canUse.add(gene);
                    canUseSet.add(gene);
                }
            }
        }

        FileWriter f;
        try {
            f = new FileWriter(rf);
            f.write("type");
            f.write(SEP + "id");
            f.write(SEP + "name");
            f.write(SEP + "reversible");
            f.write(SEP + "fast");
            f.write(SEP + "reactants");
            f.write(SEP + "products");
            f.write(SEP + "modifiers");
            f.write(SEP + "math");
            f.write(SEP + "notes");
            f.write(SEP + "clusterR");
            f.write(SEP + "clusterP");
            f.write(SEP + "width");
            f.write(SEP + "responsive");
            f.write(SEP + "positive/negative");
            f.write(SEP + "reactant is responsive");
            f.write(SEP + "product is responsive");
            f.write("\n");
            int clusterCount = AGCT.getInstance().getClustering(0).myClusteringAlgorithm
                    .getNumberOfClusters();
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
                ArrayList<Integer> clRs = new ArrayList<Integer>();
                ArrayList<Integer> clPs = new ArrayList<Integer>();
                for (int j = 0; j < gene.neighbors.size(); j++) {
                    int id = (Integer) gene.neighbors.get(j);
                    // if(id<=gene.getIdNumber())continue;
                    Gene nei = AGCT.getInstance().data
                            .getMyDomain()
                            .getGenes()
                            .get(AGCT.getInstance().data.getMyDomain().selectedGeneNumberToGeneNumber[id]);
                    //					if (!canUseSet.contains(nei) || nei.getIdNumber() < gene.getIdNumber())	continue;
                    double nangle = (Double) gene.neighbor_angles.get(j);
                    int clR = -1;
                    for (int k = 0; k < clusterCount; k++) {
                        if (nei.getClusterMemberships(0, k) == 1)
                            clR = k;
                    }

                    // if (c == i)
                    // c = -1;
                    if (gene.isVisible() && nei.isVisible()) {
                        neis.add(nei);
                        // dums.add(dumchar);
                        angles.add(nangle);
                        clRs.add(clR);
                        clPs.add(cl);
                    }
                    // }
                }

                if (neis.size() > 0) {
                    for (int k = 0; k < neis.size(); k++) {
                        // showing cluster of gene
                        //						f.write(AGCTFileWriter.geneStringDelaunay(gene, true)
                        //								+ " (c" + cl + ")\t");
                        Gene nei = neis.get(k);
                        // String dum = dums.get(k);
                        int clR = clRs.get(k);
                        int clP = clPs.get(k);
                        //						f.write(AGCTFileWriter.geneStringDelaunay(nei, false)
                        //								+ " (c" + c + ") ("
                        //								+ String.format("%.2f", angles.get(k)) + ")");

                        writeRelationData(f, gene, nei, clR, clP, angles.get(k));

                    }
                }
            }
            f.close();
        } catch (IOException e) {
            Matrix.perror("toSavingClusterAndTriangulation : writeError");
        }
    }

    private static final String SEP = ",";

    static private void writeRelationData(FileWriter f, Gene reactant,
                                          Gene product, int clusterR, int clusterP, double angle) throws IOException {
        if (!reactant.isVisible())
            return;
        if (!product.isVisible())
            return;

        double width = abs(PI / 2 - angle);
        width *= 10 / PI;
        boolean positive = angle < PI / 2;

        f.write("STATE_TRANSITION");
        f.write(SEP + "");
        f.write(SEP + "");
        f.write(SEP + "TRUE");
        f.write(SEP + "FALSE");
        f.write(SEP + (reactant.getAnnotation(Gene.Gene$Symbol).equals(Gene.Not_defined) ? reactant.getAnnotation(Gene.ID) : reactant.getAnnotation(Gene.Gene$Symbol)));
        f.write(SEP + (product.getAnnotation(Gene.Gene$Symbol).equals(Gene.Not_defined) ? product.getAnnotation(Gene.ID) : product.getAnnotation(Gene.Gene$Symbol)));
        f.write(SEP + "");
        f.write(SEP + "");
        f.write(SEP + "");
        f.write(SEP + clusterR);
        f.write(SEP + clusterP);
        f.write(SEP + width);
        f.write(SEP + reactant.isResponsive());
        f.write(SEP + (positive ? "+" : "-"));
        f.write(SEP + reactant.isResponsive());
        f.write(SEP + product.isResponsive());
        f.write("\n");

    }
}
