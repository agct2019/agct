package clusteringProfile;

import gene.Gene;
import gene.GeneList;

import java.util.ArrayList;

public class ClusteringProfile {

    final public int myClusterId;
    final public int myLigandId;
    final public int experimentId;
    final public double[] average_profiles;
    final double[] standard_deviations;
    final public boolean[] undetermined_profiles;
    final public double[] q25_profiles;
    final public double[] q50_profiles;
    final public double[] q75_profiles;
    final public ArrayList<ArrayList<ArrayList<double[]>>> allbelonginggenesplot;
    private final GeneList allGenes;
    final public int numGenes;

    public ClusteringProfile(int myClusterId, int myLigandId, int experimentId,
                             double[] averageProfiles, double[] standard_deviations,
                             boolean[] undeterminedProfiles, double[] q25Profiles,
                             double[] q50Profiles, double[] q75Profiles, ArrayList<ArrayList<ArrayList<double[]>>> allplot, GeneList allGenes, int numGenes) {
        this.myClusterId = myClusterId;
        this.myLigandId = myLigandId;
        this.experimentId = experimentId;
        average_profiles = averageProfiles;
        this.standard_deviations = standard_deviations;
        undetermined_profiles = undeterminedProfiles;
        q25_profiles = q25Profiles;
        q50_profiles = q50Profiles;
        q75_profiles = q75Profiles;
        this.allGenes = allGenes;
        this.allbelonginggenesplot = allplot;
        this.numGenes = numGenes;
    }

    /**
     * @param u
     * @return
     */
    public double time(int u) {
        return ((Double) JClusteringProfileFrame.Time_Stamps_Summaries[myLigandId].get(u)).doubleValue();
    }

    /**
     * @return
     */
    public GeneList getAllGenes() {
        return allGenes;
    }

    public ArrayList<double[]> allgenesexpression() {
        ArrayList<double[]> allExpressionList = new ArrayList<double[]>(allGenes.size());
        for (Gene g : allGenes) {
            allExpressionList.add(g.getRawCoordinate(myLigandId));
        }

        return allExpressionList;
    }

    public ArrayList<double[]> all_genes_plot() {
        return allbelonginggenesplot.get(myClusterId).get(myLigandId);
    }

    public int getNumberOfGenes() {
        return numGenes;
    }
}
