package middleMan;

import clustering.AGCTClustering_KM;
import filtering.JLeftIndicator;
import forDebug.Debug;
import test.*;

public class Clustering_painter {
    public static void plotWithAnnotations(Clustering clustering, JAGCTVisualizationPane v) {
//		if (AGCT.MYDEBUG)
//			AGCT.debug("Clustering_painter.plotWithAnnotations(clustering,jAGCTVisualizationPane)");
        if (clustering.myClusteringAlgorithm != null)
            if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM))
                ((AGCTClustering_KM) clustering.myClusteringAlgorithm).plotClusterCentersWithAnnotations(v);
            else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM))
                ((AGCTClustering_EM) clustering.myClusteringAlgorithm).plotClusterCentersWithAnnotations(v);
            else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP))
                ((AGCTClustering_CP) clustering.myClusteringAlgorithm).plotClusterCentersWithAnnotations(v);
            else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)) {
                AGCTClustering_HC r = ((AGCTClustering_HC) clustering.myClusteringAlgorithm);
                r.plotClusterCentersWithAnnotations(v);
                r.plotTree(v);
            } else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP)) {
                AGCTClustering_AP r = ((AGCTClustering_AP) clustering.myClusteringAlgorithm);
                r.plotClusterCentersWithAnnotations(v);
                r.plotAffinities(v);
            } else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM))
                ((AGCTClustering_NM) clustering.myClusteringAlgorithm).plotClusterCentersWithAnnotations(v);
    }

    public static void plot(Clustering clustering, JAGCTVisualizationPane v) {
//		if (AGCT.MYDEBUG)
//			AGCT.debug("Clustering_painter.plot(clustering,jAGCTVisualizationPane)");
        if (clustering.myClusteringAlgorithm != null)
            if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM)) {
//				if (AGCT.MYDEBUG)
//					AGCT.debug("AGCTClustering_KM.plot(JAGCTVisualizationPane)");
                plotClusterCenters((AGCTClustering_KM) clustering.myClusteringAlgorithm, v);
            } else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM))
                plotClusterCenters((AGCTClustering_EM) clustering.myClusteringAlgorithm, v);
            else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP))
                plotClusterCenters((AGCTClustering_CP) clustering.myClusteringAlgorithm, v);
            else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)) {
                AGCTClustering_HC r = ((AGCTClustering_HC) clustering.myClusteringAlgorithm);
                plotClusterCenters(r, v);
                r.plotTree(v);
            } else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP)) {
                AGCTClustering_AP r = ((AGCTClustering_AP) clustering.myClusteringAlgorithm);
                plotClusterCenters(r, v);
                r.plotAffinities(v);
            } else if (clustering.myClusteringAlgorithm.getMyReferenceName().equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM))
                plotClusterCenters((AGCTClustering_NM) clustering.myClusteringAlgorithm, v);
    }

    static public void plotClusterCenters(AGCTClustering_Algorithm clustering, JAGCTVisualizationPane vv) {
        Debug.debug("Clustering_painter # plotClusterCenters");
        int i;
        View3D v = vv.visualizationPanel.myView3D;
        String sz;
        if (clustering.plotCentersAvailable) {

            if (vv.myReferenceName.equals(JAGCTVisualizationPane.M_P)) {
                clustering.softMembershipsToHardCenters_Manifold_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
                for (i = 0; i < clustering.numberOfClusters; i++)
                    if (vv.structPlottedCluster(i)) {
                        v.g.setColor(clustering.myClustering.colorCluster[i]);
                        if (Prototype.No_Reduction)
                            sz = "" + JLeftIndicator.getInstance().clusterMemberCount(i);
                        else
                            sz = "" + JLeftIndicator.getInstance().clusterMemberCount(i) + "|" + clustering.clusterSizesWhole[i];
                        v.drawCenterVsRef((Pnt3D) clustering.hard_memberships_centers_manifold_Pnt3D.elementAt(i), i, sz);
                    }
            } else if (vv.myReferenceName.equals(JAGCTVisualizationPane.P_P)) {
                clustering.softMembershipsToHardCenters_Pca_Pnt3D(vv.xAxis, vv.yAxis, vv.zAxis);
                for (i = 0; i < clustering.numberOfClusters; i++) {
                    if (vv.structPlottedCluster(i)) {
                        v.g.setColor(clustering.myClustering.colorCluster[i]);
                        if (Prototype.No_Reduction)
                            sz = "" + JLeftIndicator.getInstance().clusterMemberCount(i);
                        else
                            sz = "" + JLeftIndicator.getInstance().clusterMemberCount(i) + "|" + clustering.clusterSizesWhole[i];
                        v.drawCenterVsRef((Pnt3D) clustering.hard_memberships_centers_pca_Pnt3D.elementAt(i), i, sz);
                    }
                }
            }
        }
    }
}
