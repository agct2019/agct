import java.io.*;
import java.util.*;

class AGCTFileWriter implements Debuggable{

    public static boolean LIGHT_SAVING = true;

    public static String DATA_Ligand = new String("@TAG"), 
	    DATA_Domain = new String("@DOMAIN_NAME"), 
	    DATA_Groups = new String("@GROUPS"), 
	    DATA_Data = new String("@DATA"), 
	    DATA_Annotation = new String("@ANNOTATION"), 
	    DATA_GeneAnnotations = new String("@GENE_ANNOTATIONS"),
	    DATA_Selection = new String ("@SELECTION_REDUCTION_FILE"),
	    DATA_Comments = new String ("//"),
	    DATA_QuestionMark = new String("?");

    public static String ANNO_Data = new String("@DATA"),
	ANNO_Format = new String("@FORMAT"),
	ANNO_Variable_Definition = new String("@VARIABLES"),
	ANNO_Empty = new String("@EMPTY");

    public static String ANNO_Format_AGCT = new String("AGCT"),
	ANNO_Format_GEO = new String("GEO");

    public static String Clustering_Born_String = new String("@Reference_String");

    public static String ANNO_Variable_Definition_Dum = new String("DUM"),
	ANNO_Variable_Definition_ID = new String("ID"),
	ANNO_Variable_Definition_Ascii_Names = new String("NAMES"),
	ANNO_Variable_Definition_Ontology_F = new String("F"),
	ANNO_Variable_Definition_Ontology_P = new String("P"),
	ANNO_Variable_Definition_Ontology_C = new String("C"),
	ANNO_Variable_Definition_Referenced = new String("REF");

    public static String [][] REPLACE_Old_New = {{"///", "BIG_SEP"}, {"//", "SMALL_SEP"}};

    AGCT myAGCT;

    AGCTFileWriter(AGCT a){
	myAGCT = a;
    }

    public void toSavingScenario(File rf){
	FileWriter f = null;
	int i;
	String s;
	if (Scenario.allStrings != null)
	    try{
		f = new FileWriter(rf);
		for (i=0;i<Scenario.allStrings.size();i++){
		    s = (String) Scenario.allStrings.elementAt(i);
		    if (Scenario.isKeepPrototypes(s))
			morePrototypes(f);
		    else if (Scenario.isKeepManifold(s))
			moreManifold(f);
		    else
			f.write(s + "\n");
		}
		f.close();
	    }catch(IOException e){
		Matrix.perror("AGCTFileWriter.class :: Writing scenario error");
	    }
    }

    public void moreManifold(FileWriter rf) throws IOException {
	int i,j,dX,dY;
	double d;
	Gene g;

	if (ControlProcess.hasTrue("manifoldProcessed")){
	    rf.write(Scenario.allKeywords[34] + Scenario.myToken + "\n");

	    //Ordered list name
	    rf.write(AGCT.Token_Ordered_List_Names_Begin + "\n");
	    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
		g = ( (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]));
		rf.write(g.name);
		if (i<myAGCT.myDomain.numberSelectedGenes-1)
		    rf.write(Scenario.myToken);
		else
		    rf.write("\n");
	    }
	    rf.write(AGCT.Token_Ordered_List_Names_End + "\n");

	    //W
	    rf.write(AGCT.Token_Matrix_W_Begin + "\n");
	    dX = myAGCT.W.dimX;
	    dY = myAGCT.W.dimY;
	    for (i=0;i<dX;i++){
		for (j=0;j<dY;j++){
		    d = myAGCT.W.coordinates[i][j];
		    rf.write(d + "");
		    if (j<dY-1)
			rf.write(Scenario.myToken);
		    else
			rf.write("\n");
		}
	    }
	    rf.write(AGCT.Token_Matrix_W_End + "\n");

	    //Manifold Eigenvalues
	    rf.write(AGCT.Token_Manifold_Eigenvalues_Begin + "\n");
	    dX = myAGCT.manifold_eigenvalues.length;
	    for (i=0;i<dX;i++){
		d = myAGCT.manifold_eigenvalues[i];
		rf.write(d + "");
		if (i<dX-1)
		    rf.write(Scenario.myToken);
		else
		    rf.write("\n");
	    }
	    rf.write(AGCT.Token_Manifold_Eigenvalues_End + "\n");

	    //M
	    rf.write(AGCT.Token_Matrix_M_Begin + "\n");
	    dX = myAGCT.M.dimX;
	    dY = myAGCT.M.dimY;
	    for (i=0;i<dX;i++){
		for (j=0;j<dY;j++){
		    d = myAGCT.M.coordinates[i][j];
		    rf.write(d + "");
		    if (j<dY-1)
			rf.write(Scenario.myToken);
		    else
			rf.write("\n");
		}
	    }
	    rf.write(AGCT.Token_Matrix_M_End + "\n");

	    rf.write(Scenario.allKeywords[35] + Scenario.myToken + "\n");
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

	if (!Prototype.No_Reduction){
	    Prototype.assertComplete();
		    
	    rf.write(Scenario.allKeywords[31] + Scenario.myToken + "\n");
	    //Closest_Center_Ordered_List

	    rf.write(Prototype.Token_Closest_Center_Ordered_List_Begin + "\n");
	    v = Prototype.Closest_Center_Ordered_List;
	    for (j = 0 ; j < v.size(); j++){
		I = (Integer) v.elementAt(j);
		i = I.intValue();
		gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(i);
		rf.write(gg.name + "");
		if (j<v.size()-1)
		    rf.write(Scenario.myToken);
	    }
	    rf.write("\n");
	    rf.write(Prototype.Token_Closest_Center_Ordered_List_End + "\n");

	    //Closest_Center_To_Cluster_Number
	    rf.write(Prototype.Token_Closest_Center_To_Cluster_Number_Begin + "\n");
	    t = Prototype.Closest_Center_To_Cluster_Number;
	    extensions = t.keys();
	    while (extensions.hasMoreElements()) {
		I = (Integer) extensions.nextElement();
		i = I.intValue();
		gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(i);
		rf.write(gg.name + Scenario.myToken + t.get(I) + "\n");
	    }
	    rf.write(Prototype.Token_Closest_Center_To_Cluster_Number_End + "\n");
	    
	    //Cluster_Number_To_Closest_Center
	    rf.write(Prototype.Token_Cluster_Number_To_Closest_Center_Begin + "\n");
	    t = Prototype.Cluster_Number_To_Closest_Center;
	    extensions = t.keys();
	    while (extensions.hasMoreElements()) {
		I = (Integer) extensions.nextElement();
		i = ( (Integer) t.get(I)).intValue();
		gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(i);
		rf.write(I + Scenario.myToken + gg.name + "\n");
	    }
	    rf.write(Prototype.Token_Cluster_Number_To_Closest_Center_End + "\n");
	    
	    //Closest_Center_To_Cluster_Points
	    rf.write(Prototype.Token_Closest_Center_To_Cluster_Points_Begin + "\n");
	    t = Prototype.Closest_Center_To_Cluster_Points;
	    extensions = t.keys();
	    while (extensions.hasMoreElements()) {
		I = (Integer) extensions.nextElement();
		i = I.intValue();
		gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(i);
		rf.write(gg.name + Scenario.myToken);
		v = (Vector) t.get(I);
		for (i=0;i<v.size();i++){
		    j = ( (Integer) v.elementAt(i) ).intValue();
		    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(j);
		    rf.write(gg.name) ;
		    if (i<v.size()-1)
			rf.write(Scenario.myToken);
		}
		rf.write("\n");
	    }
	    rf.write(Prototype.Token_Closest_Center_To_Cluster_Points_End + "\n");	
	    
	    //Closest_Center_To_Normalized_Distortions
	    rf.write(Prototype.Token_Closest_Center_To_Normalized_Distortions_Begin + "\n");
	    t = Prototype.Closest_Center_To_Normalized_Distortions;
	    extensions = t.keys();
	    while (extensions.hasMoreElements()) {
		I = (Integer) extensions.nextElement();
		i = I.intValue();
		gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(i);
		rf.write(gg.name + Scenario.myToken);
		v = (Vector) t.get(I);
		for (i=0;i<v.size();i++){
		    d = ( (Double) v.elementAt(i) ).doubleValue();
		    rf.write(d + ""); 
		    if (i<v.size()-1)
			rf.write(Scenario.myToken);
		}
		rf.write("\n");
	    }
	    rf.write(Prototype.Token_Closest_Center_To_Normalized_Distortions_End + "\n");

	    rf.write(Scenario.allKeywords[32] + Scenario.myToken + "\n");
	}
    }

    public void toSaving(File rf){
	FileWriter f = null;
	try{
	    int i, j, k, nn = 0, gid, nid, iii, totGenes;
	    boolean better;
	    Gene gg, hh;
	    String s, dums;
	    Vector vid, did, littleVector;
	    f = new FileWriter(rf);
	    
	    f.write(separator());
	    f.write("Results for domain " + myAGCT.myDomain.domainName);
	    f.write(separator());

	    f.write ("\n\n");

	    f.write(separator());
	    f.write("Processing parameters:\n");
	    f.write(myAGCT.getParameterString());
	    f.write(separator()); 

	    if ( (myAGCT.myDomain.domainGenes != null) && (myAGCT.myDomain.numberSelectedLigands > 0) ){

		f.write(separator());
		f.write("** Ligands\n\nselected (" + myAGCT.myDomain.numberSelectedLigands + ") :\n");
		for (i=0;i<myAGCT.myDomain.numberLigands;i++){
		    if (myAGCT.myDomain.selectedLigand[i]){
			s = (String) ((Vector) myAGCT.myDomain.domainLigands.elementAt(i)).elementAt(0);

			f.write(s + " ");
		    }
		}
		f.write(separator());
		f.write("\n\n");
	    }

	    if ( (myAGCT.myDomain.domainGenes != null) && (myAGCT.myDomain.numberSelectedGenes > 0) ){
		f.write(separator());
		if (Prototype.No_Reduction){
		    f.write("** Genes\n\nSelected (" + myAGCT.myDomain.numberSelectedGenes + "):\n");
		}else{
		    f.write("** Prototypes\n\nSelected (" + myAGCT.myDomain.numberSelectedGenes + ") + their clusters:\n");
		}

		for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
		    gid = myAGCT.myDomain.selectedGeneNumberToGeneNumber[i];
		    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(gid);
		    f.write(gg.name + ((gg.asciiName != null) ? "_" + gg.asciiName + "  " : "  "));
		    if (Prototype.No_Reduction)
			f.write(" ");
		    else{
			f.write(" [ ");
			vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
			did = (Vector) Prototype.Closest_Center_To_Normalized_Distortions.get(new Integer(gid));
			for (iii=0;iii<vid.size();iii++){
			    nid = ( (Integer) vid.elementAt(iii) ).intValue();
			    hh = myAGCT.getDirectGene(nid);
			    f.write("(" + hh.name + ((hh.asciiName != null) ? "_" + hh.asciiName + "  "  :  "  ") + "  :  " + did.elementAt(iii) + ")");
			}
			f.write(" ]\n");
		    }
		}
		f.write(separator());
		f.write("\n\n");

		if (AGCT.LOCAL_DIMENSION_COMPUTED){
		    f.write(separator());
		    if (Prototype.No_Reduction){
			f.write("** Local dimension\n\nBy Gene (" + myAGCT.myDomain.numberSelectedGenes + "), as (average) + list of dimensions as f(k):\n");
		    }else{
			f.write("** Local dimension\n\nBy Prototype (" + myAGCT.myDomain.numberSelectedGenes + "), as (average) + list of dimensions as f(k):\n");
		    }

		    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
			gid = myAGCT.myDomain.selectedGeneNumberToGeneNumber[i];
			gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(gid);
			f.write(gg.name + ((gg.asciiName != null) ? "_" + gg.asciiName + "" : "") + "\t(" + gg.averageDimension + ")\t");

			for (iii=AGCT.Min_K_Dimension_Estimation;iii<AGCT.Max_K_Dimension_Estimation;iii++){
			    f.write(( (int) gg.local_dimension[iii]) + "\t");
			}
			f.write("\n");
		    }

		    f.write(separator());
		    f.write("\n\n");   
		}
	    }
	    
	    if (ControlProcess.hasTrue("manifoldProcessed")){
		f.write(separator());
		f.write("** Manifold\n\n");
		
		f.write("Eigenvalues (sorted): ");

		for (i=0;i<myAGCT.manifold_eigenvalues.length;i++)
		    f.write(DF.format(myAGCT.manifold_eigenvalues[i]) + " ");

		if (myAGCT.annotationsExist){
		    f.write("\n\nNatural Neighbors for annotations\n\n");
		    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
			gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
			nn = 0;
			if ( ( (gg.annotationsF != null) || (gg.annotationsP != null) || (gg.annotationsC != null) ) && (gg.neighbors != null) ){
			    for (j=0;j<gg.neighbors.size();j++){
				hh = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[( (Integer) gg.neighbors.elementAt(j)).intValue()]);
				if ( (hh.annotationsF != null) || (hh.annotationsP != null) || (hh.annotationsC != null) ){
				    if (nn == 0){
					f.write("  ");
					if (gg.annotationsF != null){
					    f.write(JAnnotationFrame.TAG_F + " ");
					    for (k=0;k<gg.annotationsF.size();k++){
						if (!((String[])gg.annotationsF.elementAt(k))[0].equals(Domain.ANNOTATION_F))
						    Matrix.perror("AGCTFileWriter.class :: annotation != F in F annotations");
						f.write(((String[])gg.annotationsF.elementAt(k))[1]);
						if (k<gg.annotationsF.size()-1)
						    f.write(", ");
					    }
					    f.write(";\n");
					}
					if (gg.annotationsP != null){
					    f.write(JAnnotationFrame.TAG_P + " ");
					    for (k=0;k<gg.annotationsP.size();k++){
						if (!((String[])gg.annotationsP.elementAt(k))[0].equals(Domain.ANNOTATION_P))
						    Matrix.perror("AGCTFileWriter.class :: annotation != P in P annotations");
						f.write(((String[])gg.annotationsP.elementAt(k))[1]);
						if (k<gg.annotationsP.size()-1)
						    f.write(", ");
					    }
					    f.write(";\n");
					}
					if (gg.annotationsC != null){
					    f.write(JAnnotationFrame.TAG_C + " ");
					    for (k=0;k<gg.annotationsC.size();k++){
						if (!((String[])gg.annotationsC.elementAt(k))[0].equals(Domain.ANNOTATION_C))
						    Matrix.perror("AGCTFileWriter.class :: annotation != C in C annotations");
						f.write(((String[])gg.annotationsC.elementAt(k))[1]);
						if (k<gg.annotationsC.size()-1)
						    f.write(", ");
					    }
					    f.write(";\n");
					}
					f.write("have natural neighbors :\n");
				    }
				    f.write("   " + (j+1) + " -- ");
				    if (hh.annotationsF != null){
					f.write(JAnnotationFrame.TAG_F + " ");
					for (k=0;k<hh.annotationsF.size();k++){
					    if (!((String[])hh.annotationsF.elementAt(k))[0].equals(Domain.ANNOTATION_F))
						Matrix.perror("AGCTFileWriter.class :: annotation != F in F annotations");
					    f.write(((String[])hh.annotationsF.elementAt(k))[1]);
					    if (k<hh.annotationsF.size()-1)
						f.write(", ");
					}
					f.write(";\n");
				    }
				    if (hh.annotationsP != null){
					f.write(JAnnotationFrame.TAG_P + " ");
					for (k=0;k<hh.annotationsP.size();k++){
					    if (!((String[])hh.annotationsP.elementAt(k))[0].equals(Domain.ANNOTATION_P))
						Matrix.perror("AGCTFileWriter.class :: annotation != P in P annotations");
					    f.write(((String[])hh.annotationsP.elementAt(k))[1]);
					    if (k<hh.annotationsP.size()-1)
						f.write(", ");
					}
					f.write(";\n");
				    }
				    if (hh.annotationsC != null){
					f.write(JAnnotationFrame.TAG_C + " ");
					for (k=0;k<hh.annotationsC.size();k++){
					    if (!((String[])hh.annotationsC.elementAt(k))[0].equals(Domain.ANNOTATION_C))
						Matrix.perror("AGCTFileWriter.class :: annotation != C in C annotations");
					    f.write(((String[])hh.annotationsC.elementAt(k))[1]);
					    if (k<hh.annotationsC.size()-1)
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
		    f.write(separator());

		    //Saves the average dimension per highlight performed

		    if ( (AGCT.LOCAL_DIMENSION_COMPUTED) && (AGCT.Referenced_Available) ){
			f.write(separator()); 
			f.write("** Average local dimensions for each highlight tag (last highlights loaded)\n\n");

			for (k = 0; k < myAGCT.myDomain.orderedReferences.size();k++){
			    f.write(myAGCT.myDomain.orderedReferences.elementAt(k) + "\t" + myAGCT.myDomain.highlightAverageDimension.elementAt(k) + "\n");
			}		

			f.write(separator()); 
		    }

		    //Saves the average dimension per tags
		    Vector ddumv;
		    String sstr, sstag;
		    if ( (Statistics.HISTOGRAM_TAG_AVERAGE_DIMENSION != null) && (Statistics.HISTOGRAM_TAG_AVERAGE_DIMENSION.size() > 0) ){
			f.write(separator());
			f.write("** Average local dimensions for each gene tag (P/F/C)\n\n");

			if ( (myAGCT.myTabbedPane.myManifoldPane.statP != null)  ){
			    sstag = Domain.ANNOTATION_P;
			    f.write(" * Tag " + sstag + "\n");
			    for (k = Statistics.HISTOGRAM_TAG_AVERAGE_DIMENSION.size()-1; k>=0; k--){
				ddumv = (Vector) Statistics.HISTOGRAM_TAG_AVERAGE_DIMENSION.elementAt(k);
				sstr = (String) ddumv.elementAt(0);
				if (sstr.equals(sstag)){
				    f.write((String) ddumv.elementAt(1) + " : " + (Double) ddumv.elementAt(2) + "\n");
				}
			    }
			    f.write("\n");
			}
			if ( (myAGCT.myTabbedPane.myManifoldPane.statF != null)  ){
			    sstag = Domain.ANNOTATION_F;
			    f.write(" *** Tag " + sstag + "\n\n");
			    for (k = Statistics.HISTOGRAM_TAG_AVERAGE_DIMENSION.size()-1; k>=0; k--){
				ddumv = (Vector) Statistics.HISTOGRAM_TAG_AVERAGE_DIMENSION.elementAt(k);
				sstr = (String) ddumv.elementAt(0);
				if (sstr.equals(sstag)){
				    f.write((String) ddumv.elementAt(1) + " : " + (Double) ddumv.elementAt(2) + "\n");
				}
			    }
			    f.write("\n");
			}
			if ( (myAGCT.myTabbedPane.myManifoldPane.statC != null)  ){
			    sstag = Domain.ANNOTATION_C;
			    f.write(" * Tag " + sstag + "\n");
			    for (k = Statistics.HISTOGRAM_TAG_AVERAGE_DIMENSION.size()-1; k>=0; k--){
				ddumv = (Vector) Statistics.HISTOGRAM_TAG_AVERAGE_DIMENSION.elementAt(k);
				sstr = (String) ddumv.elementAt(0);
				if (sstr.equals(sstag)){
				    f.write((String) ddumv.elementAt(1) + " : " + (Double) ddumv.elementAt(2) + "\n");
				}
			    }
			    f.write("\n");
			}
			f.write(separator());
		    }
		}
		    
		f.write(separator());
		f.write("\n\n");
	    }

	    if (ControlProcess.hasTrue("pcaProcessed")){
		f.write(separator());
		f.write("** PCA\n\n");
		
		f.write("Eigenvalues (sorted): ");

		for (i=0;i<myAGCT.pca_eigenvalues.length;i++)
		    f.write(DF.format(myAGCT.pca_eigenvalues[i]) + " ");

		f.write(separator());
		f.write("\n\n");
	    }

	    if ( (ControlProcess.hasTrue("softClusteringProcessed")) || (ControlProcess.hasTrue("hardClusteringProcessed")) ){
		int jj;
		for (jj=0;jj<myAGCT.nbClustering;jj++){
		    f.write(separator());
		    f.write("** Clustering #" + jj + "\n\n");
		    
		    f.write("Method: " + myAGCT.getClustering(jj).myClusteringAlgorithm + "\n\n");
		    f.write("Final value for objective function: " + myAGCT.getClustering(jj).myClusteringAlgorithm.finalClusteringObjectiveFunction + "\n\n");

		    if ( ( (AGCT.Filter_Similarities_With_Distances) || (AGCT.Filter_Triangulation_With_Distances) ) && (myAGCT.getClustering(jj).myClusteringAlgorithm.hard_average_distance != null) ){
			f.write("* Statistics on within-cluster distances:\n");
			f.write(myAGCT.getClustering(jj).myClusteringAlgorithm.getDistanceInformation() + "\n\n");
		    }
		    


		    if (AGCT.LOCAL_DIMENSION_COMPUTED){
			f.write(" * average cluster dimension + average dimension for each K\n");

			for (i=0;i<myAGCT.getClustering(jj).myClusteringAlgorithm.nclusters;i++){
			    f.write("C#" + i + " --- Avg = " + myAGCT.getClustering(jj).myClusteringAlgorithm.average_dimension_clusters[i] + " --- for each K : ");
			    for (iii=AGCT.Min_K_Dimension_Estimation;iii<AGCT.Max_K_Dimension_Estimation;iii++){
				f.write(( (int) myAGCT.getClustering(jj).myClusteringAlgorithm.dimension[i][iii]) + " ");
			    }
			    f.write("\n");
			}

			if ( ((AGCTClustering_Algorithm) myAGCT.getClustering(jj).myClusteringAlgorithm).hard_clustering_computed ){
			    f.write("\n * histogram of average dimensions for each cluster\n");
			    Vector vv = (Vector) ((AGCTClustering_Algorithm) myAGCT.getClustering(jj).myClusteringAlgorithm).average_dimension_histogram;
			    Vector vvv, vvvv;
			    for (i=0;i<myAGCT.getClustering(jj).myClusteringAlgorithm.nclusters;i++){
				f.write(" ** cluster# " + i + "\n");
				vvv = (Vector) vv.elementAt(i);
				for (iii=0;iii<vvv.size();iii++){
				    vvvv = (Vector) vvv.elementAt(iii);
				    f.write(vvvv.elementAt(0) + "\t" + vvvv.elementAt(1) + "\n");
				}
				f.write("\n");
			    }
			}
		    }

		    f.write("\n\n");   

		    f.write("Additional data on method: " + myAGCT.getClustering(jj).myClusteringAlgorithm.getAdditionalSavingData() + "\n\n");
		    
		    
		    if ( (AGCT.Referenced_Available) && (myAGCT.nbClustering == 1) ){
			f.write("Summary of references by clusters:\n");
			for (i=0;i<myAGCT.getClustering(jj).myClusteringAlgorithm.nclusters;i++){
			    f.write("Cluster "+i+": ");
			    for (j=0;j<myAGCT.getClustering(jj).myClusteringAlgorithm.referencesByCluster.length;j++)
				f.write(myAGCT.getClustering(jj).myClusteringAlgorithm.referencesByCluster[j][i] + "\t");
			    f.write("\n");
			}
			f.write("\n");
		    }

		    f.write("Memberships:\n");
		    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
			gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
			f.write(gg.name + ": ");
			for (j=0;j<myAGCT.getClustering(jj).myClusteringAlgorithm.nclusters;j++)
			    f.write(DF.format(gg.getClusterMemberships(jj,j)) + ", ");
			f.write("\n");
		    }
		    f.write("\n\n");
		    f.write("Composition of clusters (majority rule for memberships):");
		    for (j=0;j<myAGCT.getClustering(jj).myClusteringAlgorithm.nclusters;j++){
			f.write("\n\nCluster " + j);
			if (Prototype.No_Reduction)
			    f.write(": ");
			else
			    f.write("(prototypes + their clusters): ");

			for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
			    gid = myAGCT.myDomain.selectedGeneNumberToGeneNumber[i];
			    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(gid);
			    
			    better = true;
			    for (k=0;k<myAGCT.getClustering(jj).myClusteringAlgorithm.nclusters;k++){
				if (gg.getClusterMemberships(jj,j) < gg.getClusterMemberships(jj,k))
				    better = false;
			    }
			    if (better){
				f.write(gg.name + ((gg.asciiName != null) ? "_" + gg.asciiName + "  " : "  "));

				if (Prototype.No_Reduction)
				    f.write(" ");
				else{
				    f.write(" [ ");
				    vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
				    did = (Vector) Prototype.Closest_Center_To_Normalized_Distortions.get(new Integer(gid));
				    for (iii=0;iii<vid.size();iii++){
					nid = ( (Integer) vid.elementAt(iii) ).intValue();
					hh = myAGCT.getDirectGene(nid);
					f.write("(" + hh.name + ((hh.asciiName != null) ? "_" + hh.asciiName + "  "  :  "  ") + "  :  " + did.elementAt(iii) + ")");
				    }
				    f.write(" ]\n");
				}
			    }
			}

			totGenes = 0;
			f.write("\n\nCluster names " + j + ", preformatted for highlight:\n");
			for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
			    gid = myAGCT.myDomain.selectedGeneNumberToGeneNumber[i];
			    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(gid);
			    better = true;
			    for (k=0;k<myAGCT.getClustering(jj).myClusteringAlgorithm.nclusters;k++){
				if (gg.getClusterMemberships(jj,j) < gg.getClusterMemberships(jj,k))
				    better = false;
			    }
			    if (better){
				//if (totGenes > 0)
				//  f.write(", ");
				dums = gg.name;
				f.write(dums.replace(',','_') + ",C" + j + "\n");
				totGenes++;
				/*if ( (gg.asciiName != null) && (!gg.asciiName.equals("")) ){
				    if (totGenes > 0)
					f.write(", ");
				    f.write (gg.asciiName + ", ");
				    totGenes++;
				    }*/
				if (!Prototype.No_Reduction){
				    vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
				    for (iii=1;iii<vid.size();iii++){
					nid = ( (Integer) vid.elementAt(iii) ).intValue();
					hh = myAGCT.getDirectGene(nid);
					//if (totGenes > 0)
					//   f.write(", ");
					dums = hh.name;
					f.write(dums.replace(',','_') + ",C" + j + "\n");
					totGenes++;
					/*if ( (hh.asciiName != null) && (!hh.asciiName.equals("")) ){
					    if (totGenes > 0)
						f.write(", ");
					    f.write (hh.asciiName + ", ");
					    totGenes++;
					    }*/
				    }
				}
			    }
			}

			f.write("\n\nCluster annotations " + j + ": ");
			for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
			    gid = myAGCT.myDomain.selectedGeneNumberToGeneNumber[i];
			    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(gid);
			    
			    better = true;
			    for (k=0;k<myAGCT.getClustering(jj).myClusteringAlgorithm.nclusters;k++){
				if (gg.getClusterMemberships(jj,j) < gg.getClusterMemberships(jj,k))
				    better = false;
			    }
			    if (better){
				if (gg.annotationsF != null)
				    for (k=0;k<gg.annotationsF.size();k++)
					f.write("[" + ((String[])gg.annotationsF.elementAt(k))[0] + "  " + ((String[])gg.annotationsF.elementAt(k))[1] + "] --- ");
				if (gg.annotationsP != null)
				    for (k=0;k<gg.annotationsP.size();k++)
					f.write("[" + ((String[])gg.annotationsP.elementAt(k))[0] + "  " + ((String[])gg.annotationsP.elementAt(k))[1] + "] --- ");
				if (gg.annotationsC != null)
				    for (k=0;k<gg.annotationsC.size();k++)
					f.write("[" + ((String[])gg.annotationsC.elementAt(k))[0] + "  " + ((String[])gg.annotationsC.elementAt(k))[1] + "] --- ");
				

				if (Prototype.No_Reduction)
				    f.write(" ");
				else{
				    f.write(" Neighbor annotations :: [ ");
				    vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
				    did = (Vector) Prototype.Closest_Center_To_Normalized_Distortions.get(new Integer(gid));
				    for (iii=0;iii<vid.size();iii++){
					nid = ( (Integer) vid.elementAt(iii) ).intValue();
					hh = myAGCT.getDirectGene(nid);
					if (hh.annotationsF != null)
					    for (k=0;k<hh.annotationsF.size();k++)
						f.write("[" + ((String[])hh.annotationsF.elementAt(k))[0] + "  " + ((String[])hh.annotationsF.elementAt(k))[1] + "] --- ");
					if (hh.annotationsP != null)
					    for (k=0;k<hh.annotationsP.size();k++)
						f.write("[" + ((String[])hh.annotationsP.elementAt(k))[0] + "  " + ((String[])hh.annotationsP.elementAt(k))[1] + "] --- ");
					if (hh.annotationsC != null)
					    for (k=0;k<hh.annotationsC.size();k++)
						f.write("[" + ((String[])hh.annotationsC.elementAt(k))[0] + "  " + ((String[])hh.annotationsC.elementAt(k))[1] + "] --- ");
				    }
				    f.write(" ]\n");
				}
			    }
			}
		    }
		    
		    if (jj<myAGCT.nbClustering-1)
			f.write("\n\n");
		    
		    f.write(separator());
		    f.write("\n\n");
		}
	    }

	    if ( (Statistics.allTests != null) && (Statistics.allTests != "") ){
		int jj;
		f.write(separator());
		f.write("** Chi-square tests performed on clusterings\n\n");
		for (jj=0;jj<Statistics.allTests.length();jj++)
		    f.write(Statistics.allTests.charAt(jj));

		f.write("Complete and ordered list of all annotation tags, given as: Type Tag Clustering Cluster P-value P*-value\n");
		    
		if (Statistics.annotationsBH != null)
		    for (i=0;i<Statistics.annotationsBH.size();i++){
			littleVector = (Vector) Statistics.annotationsBH.elementAt(i);
			f.write(littleVector.elementAt(0) + "\t" + littleVector.elementAt(1) + "\t" + littleVector.elementAt(2) + "\t" + littleVector.elementAt(3) + "\t" + littleVector.elementAt(4) + "\t" + littleVector.elementAt(5) + "\n");
		    }

		f.write(separator());
		f.write("\n\n");
	    }

	    f.write(separator());
	    f.write("** Processed genes list\n\n");

	    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
		gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
		f.write(gg.toString());
	    }

	    f.write(separator());
	    f.write("\n\n");

	    f.close();
	}catch(IOException e){
	    Matrix.perror("AGCTFileWriter.class :: Writing error");
	}
    }

    public void toSavingChi2(final File rf, final AGCT ag){
	Thread t = new Thread(){
		public void run(){
		    FileWriter f = null;
		    try{
			int i, j, k, nc, ncur;
			Vector littleVector;
			Vector cumulative = new Vector(), element, vTable;
			double currentValue, bestValue = -1, pcur, ad, bc;
			String refann, curann, signe;
			f = new FileWriter(rf);
			
			f.write(separator());
			f.write("Saving CHI2 Results for domain " + myAGCT.myDomain.domainName);
			f.write(separator());
			
			f.write ("\n\n");
			
			f.write(separator());
			f.write("Processing parameters:\n");
			f.write(myAGCT.getParameterString());
			f.write(separator()); 
			
			f.write("** Chi-square tests performed on clusterings\n\n");
			for (i=0;i<Statistics.allTests.length();i++)
			    f.write(Statistics.allTests.charAt(i));
			
			f.write("Complete and ordered list of all annotation tags, given as: Type Tag Clustering Cluster P-value P*-value [Chi2 table]\n");
			
			if (Statistics.annotationsBH != null){
			    AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Saving Chi2 Stats", 2*Statistics.annotationsBH.size());
			    for (i=0;i<Statistics.annotationsBH.size();i++){
				littleVector = (Vector) Statistics.annotationsBH.elementAt(i);
				f.write(littleVector.elementAt(0) + "\t" + littleVector.elementAt(1) + "\t" + littleVector.elementAt(2) + "\t" + littleVector.elementAt(3) + "\t" + littleVector.elementAt(4) + "\t" + littleVector.elementAt(6) + "\t");
				vTable = (Vector) littleVector.elementAt(5);
				f.write("[" + vTable.elementAt(0) + ", " + vTable.elementAt(1) + ", " + vTable.elementAt(2) + ", " + vTable.elementAt(3) + "]\n");
				cc.increment();
			    }
			    f.write("\n\n\nSet of tags for each cluster, having P*-Values <= Satistics.LIMIT_P (" + Statistics.LIMIT_P_CHI2 +")\n");
			    

			    nc = ag.getClustering(Statistics.refNumClustering).myClusteringAlgorithm.nclusters;
			    for (k=0;k<nc;k++){
				f.write("Cluster " + k + ":\n\n");
				for (j=0;j<3;j++){
				    refann = null;
				    if (j==0)
					refann = Domain.ANNOTATION_F;
				    else if (j==1)
					refann = Domain.ANNOTATION_P;
				    else if (j==2)
					refann = Domain.ANNOTATION_C;
				    f.write("Annotation " + refann + "\n");
				    for (i=0;i<Statistics.annotationsBH.size();i++){
					littleVector = (Vector) Statistics.annotationsBH.elementAt(i);
					ncur = ((Integer) littleVector.elementAt(3)).intValue();
					curann = (String) littleVector.elementAt(0);
					pcur = ((Double) littleVector.elementAt(6)).doubleValue();
					vTable = (Vector) littleVector.elementAt(5);
					ad = ((Double) vTable.elementAt(0)).doubleValue() * ((Double) vTable.elementAt(3)).doubleValue();
					bc = ((Double) vTable.elementAt(1)).doubleValue() * ((Double) vTable.elementAt(2)).doubleValue();

					if ( (ncur == k) && (curann.equals(refann)) && (pcur < Statistics.LIMIT_P_CHI2) ){
					    if (ad > bc)
						f.write("+ ");
					    else if (ad < bc)
						f.write("- ");
					    else
						f.write("/ ");

					    f.write((String) littleVector.elementAt(1) + "\t" + littleVector.elementAt(6) + "\n");
					}
				    }
				    f.write("\n");
				}
				f.write("\n");
			    }


			    f.write("Cumulative plot (x,y), in which x=P*-value and y=%tags having P*-value <= x\n");
			    
			    for (i=0;i<Statistics.annotationsBH.size();i++){
				littleVector = (Vector) Statistics.annotationsBH.elementAt(i);
				currentValue = ((Double) littleVector.elementAt(6)).doubleValue();
				if (i==0){
				    bestValue = currentValue;
				    element = new Vector();
				    element.addElement(new Double(currentValue));
				    element.addElement(new Double(i+1));
				    cumulative.addElement(element);
				}else{
				    if (currentValue == bestValue){
					element = (Vector) cumulative.get(cumulative.size()-1);
					element.setElementAt(new Double(i+1), 1);
				    }else{
					if (currentValue < bestValue)
					    Matrix.perror("AGCTFileWriter.class :: bad ordering of values");
					bestValue = currentValue;
					element = new Vector();
					element.addElement(new Double(currentValue));
					element.addElement(new Double(i+1));
					cumulative.addElement(element);
				    }
				}
				cc.increment();
			    }
			    cc.end();
			    
			    currentValue = ((Double) ((Vector) cumulative.elementAt(cumulative.size() - 1)).elementAt(1)).doubleValue();
			    
			    for (i=0;i<cumulative.size();i++){
				littleVector = (Vector) cumulative.elementAt(i);
				f.write(littleVector.elementAt(0) + "\t" + ( ( (Double) littleVector.elementAt(1) ).doubleValue() / currentValue) + "\n");
			    }
			}
			f.write(separator());
			
			
			f.close();
		    }catch(IOException e){
			Matrix.perror("AGCTFileWriter.class :: Writing error for Chi2 statistics");
		    }
		}
	    };
	t.start();
    }

    public void toSavingWC(File rf){
	FileWriter f = null;
	int i,j,k,max;
	double [] curv;
	Gene gg;
	try{
	    f = new FileWriter(rf);
	    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
		gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
		f.write(i + ";" + gg.name);
		for (j=0;j<gg.finalCoordinates.length;j++){
		    curv = gg.finalCoordinates[j];
		    if ( (curv != null) && (curv.length > 0) ){
			f.write(";" + (String) ((Vector) myAGCT.myDomain.domainLigands.elementAt(j)).elementAt(0) + ";");
			if (curv.length < AGCT.MAX_WAVELET_COEFFICIENTS_SAVE)
			    max = curv.length;
			else
			    max = AGCT.MAX_WAVELET_COEFFICIENTS_SAVE;
			for (k=0;k<max;k++){
			    if (k>0)
				f.write(";");
			    f.write(curv[k] + "");
			}
			if ( (AGCT.MAX_WAVELET_COEFFICIENTS_SAVE > max) && (AGCT.COMPLETE_WITH_EMPTY_DATA) )
			    for (k=0;k<AGCT.MAX_WAVELET_COEFFICIENTS_SAVE - max;k++)
				f.write(";");
		    }
		}
		f.write("\n");
	    }
	    f.close();
	}catch(IOException e){
	    Matrix.perror("AGCTFileWriter.class :: Writing scenario error");
	}
    }

    public void toSavingClustering(final File rf, final int nc){
	Thread t = new Thread(){
		public void run(){
		    FileWriter f = null;
		    Clustering cc = myAGCT.getClustering(nc);
		    int i;
		    Gene gg;
		    try{
			f = new FileWriter(rf);
			f.write(AGCTFileWriter.Clustering_Born_String + Scenario.myToken + cc.bornString + "\n");
			f.write(AGCTFileWriter.DATA_Data + Scenario.myToken + "\n");

			AGCTCounter ccc = new AGCTCounter(myAGCT.myInformationFrame, "Saving clustering", myAGCT.myDomain.numberSelectedGenes);
			for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
				gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
				f.write (i + Scenario.myToken + gg.name + Scenario.myToken + cc.myClusteringAlgorithm.majorityCluster(i) + "\n");
				ccc.increment();
			}
			ccc.end();

			f.close();
		    }catch(IOException e){
			Matrix.perror("AGCTFileWriter.class :: Writing error for saving Clustering");
		    }
		}
	    };
	t.start();
    }

    public static void loadClustering(File rf, AGCT ap){
	FileReader e;
	BufferedReader br;
	String st,sc,type="",name;
	int index,i;
	StringTokenizer t;
	Vector allPoints = new Vector(), element;
	boolean collectingPoints = false;
	//each element = Vector (index, gene name, cluster)
	boolean okToRegister = true;
	Gene gg;

	try{
	    e = new FileReader(rf);
	}catch(FileNotFoundException ex){
	    ap.myInformationFrame.setText("The clustering you try to load does not exist");
	    return;
	}


	br = new BufferedReader(e);
	try{
	    while ( (st=br.readLine()) != null){
		if ( (st.length()>1) && (!st.substring(0,AGCTFileWriter.DATA_Comments.length()).equals(AGCTFileWriter.DATA_Comments)) ){
		    t = new StringTokenizer(st, Scenario.myToken);
		    if (!collectingPoints){
			sc = t.nextToken();
			if (sc.equals(AGCTFileWriter.DATA_Data))
			    collectingPoints = true;
			else if (sc.equals(AGCTFileWriter.Clustering_Born_String))
			    type = t.nextToken();
		    }else{
			if (type.equals(new String("")))
			    okToRegister = false;
			else{
			    element = new Vector();
			    element.add(new Integer(Integer.parseInt(t.nextToken())));
			    element.add(new String(t.nextToken()));
			    element.add(new Integer(Integer.parseInt(t.nextToken())));
			    allPoints.addElement(element);
			}
		    }
		}
	    }
	    e.close();
	}catch(IOException ex){
	    ap.myInformationFrame.setText("IOError when loading clustering");
	    return;
	}

	if (!okToRegister)
	    System.out.println("This clustering is not applicable (not found type and number of clusters)");

	if (Clustering.getIndex(type.substring(0,2)) != 1){
	    System.out.println("Not a K-means algorithm");
	    okToRegister = false;
	}

	if ( (okToRegister) && (allPoints.size() != ap.myDomain.numberSelectedGenes) ){
	    System.out.println("This clustering does not contain the same number of points (" + ap.myDomain.numberSelectedGenes + " != " + allPoints.size());
	    okToRegister = false;
	}
	if (okToRegister)
	    for (i=0;i<allPoints.size();i++){
		element = (Vector) allPoints.elementAt(i);
		name = (String) element.elementAt(1);
		index = ((Integer) element.elementAt(0)).intValue();
		gg = (Gene) ap.myDomain.domainGenes.elementAt(ap.myDomain.selectedGeneNumberToGeneNumber[index]);
		if (!name.equals(gg.name)){
		    System.out.println("Misatch in gene #" + index + " name for clustering: " + name + " != " + gg.name);
		    okToRegister = false;
		}
	    }
	if (okToRegister){
	    ap.nbClustering++;
	    Clustering newClustering;
	    newClustering = new Clustering(ap,ap.nbClustering-1);
	    newClustering.getOptions(type);
	    ((AGCTClustering_KM) newClustering.myClusteringAlgorithm).toClusteringLite(allPoints);	    

	    if (ap.nbClustering == 1){
		ap.allClusterings = new Vector();
		ap.myTabbedPane.myPCAPane.currentClustering = 0;
		ap.myTabbedPane.myManifoldPane.currentClustering = 0;
		ap.myTabbedPane.myClusteringPane.activateSearch();
	    }
	    
	    ap.allClusterings.add(newClustering);
	    ap.myClusteringProfileFrame.addClusteringLF(newClustering, ap.nbClustering-1);
	
	    ap.updateClusteringLF(JAGCTVisualizationPane.stringClustering(ap.nbClustering-1));
	    ap.myTabbedPane.myClusteringPane.listClusteringPane.setText(ap.myTabbedPane.myClusteringPane.clusteringList());

	    newClustering.generate_VV();
	    ap.myTabbedPane.myManifoldPane.setMembershipsButtons(true);
	    ap.myTabbedPane.myPCAPane.setMembershipsButtons(true);
	    ap.myTabbedPane.myCorrelationPane.setMembershipsButtons(true);
	    
	    newClustering.myClusteringAlgorithm.fillGeneMemberships();
	}

	e = null;
    }


    public void toSavingClusterProfile(final File rf){
	Thread t = new Thread(){
		public void run(){
		    FileWriter f = null;
		    try{
			int i, j, k, index = 0;
			JClusteringProfileGridPlot gp = myAGCT.myClusteringProfileFrame.graphPanel;
			JClusteringProfileStamp js;

			f = new FileWriter(rf);
			
			f.write(separator());
			f.write("Saving Cluster Profiles for domain " + myAGCT.myDomain.domainName);
			f.write(separator());
			
			f.write ("\n\n");
			
			f.write(separator());
			f.write("Processing parameters:\n");
			f.write(myAGCT.getParameterString());
			f.write(separator()); 
			
			f.write("*** All average time series for each cluster, given as Time_Stamp Average Q25 Q50 Q75\n\n\n");
			
			f.write(separator());
			
			for (i=0;i<gp.nclusters;i++){
			    f.write("** Cluster" + i + ":\n\n");
			    index = 0;
			    for (j=0;j<myAGCT.myDomain.numberLigands;j++){
				if (myAGCT.myDomain.selectedLigand[j]){
				    js = gp.all_stamps[i][index];
				    f.write( ((String) ((Vector) myAGCT.myDomain.domainLigands.elementAt(j)).elementAt(0)) + ":\n");
				    for (k=0;k<js.average_profiles.length;k++){
					if (!js.undetermined_profiles[k])
					    f.write(js.time(k) + "\t" + js.average_profiles[k] + "\t" + js.q25_profiles[k] + "\t" + js.q50_profiles[k] + "\t" + js.q75_profiles[k] + "\n");
				    }
				    index++;
				    f.write("\n");
				}
			    }
			    f.write("\n");
			}
			f.close();
		    }catch(IOException e){
			Matrix.perror("AGCTFileWriter.class :: Writing error for Cluster profiles");
		    }
		}
	    };
	t.start();
    }


    public String geneStringDelaunay(Gene ggn, boolean proto){
	String dumString = "";
	String refGene = "_(R)";
	if (ggn.asciiName != null)
	    dumString += ggn.asciiName + "_";
	dumString += ggn.name;
	if (proto)
	    dumString += "_(*)";
	if ( (AGCT.Referenced_Available) && (ggn.referenced) )
	    dumString += refGene;
	return dumString;
    }

    public String keyForNames(String n1, String n2){
	String ret;
	if (n1.compareTo(n2) < 0)
	    ret = n1 + "_" + n2;
	else
	    ret = n2 + "_" + n1;
	return ret;
    }

    public void increment(Vector v, int coord){
	//coord in 0, 1; 0 iff P-significant
	//Vector a 2-dim Integer vector
	addTo(v,coord,1);
    }
    
    public void addTo(Vector v, int coord, int value){
	//coord in 0, 1; 0 iff P-significant
	//Vector a 2-dim Integer vector
	Integer i = (Integer) v.elementAt(coord);
	v.setElementAt(new Integer((i.intValue())+value),coord);
    }
    
    public Vector orderedCorrelationAsciiNames(Hashtable ht){
	//Stores ht in decreasing values of the ratio #sig / #total
	Enumeration extensions = ht.keys();
	int num, den, index = 0, i;
	String R;
	Vector S, ret1 = new Vector(), ret2 = null;
	double [] ratios = new double [ht.size()];
	int [] indexes = new int [ht.size()];
	if(extensions != null) {
	    while (extensions.hasMoreElements()) {
		R = (String) extensions.nextElement();
		S = (Vector) ht.get(R);

		num = ((Integer) S.elementAt(2)).intValue();
		den = ((Integer) S.elementAt(3)).intValue();

		if (den == 0)
		    ratios[index] = 1.0;
		else
		    ratios[index] = 1.0 - ((double) num) / ((double) den);
		indexes[index] = index;

		index++;
		ret1.addElement(S);
	    }
	    QuickSort.quicksort(ratios, indexes);
	    ret2 = new Vector();
	    for (i=0;i<ret1.size();i++)
		ret2.addElement((Vector) ret1.elementAt(indexes[i]));
	}
	return ret2;
    }

    /*public Vector orderedEntropyCorrelationAsciiNamesAggregated(Vector all){
	// input = output of orderedCorrelationAsciiNamesAggregated
	// output = set of ordered entropy, increasing, giving for each element, a Vector
	// v = Species, entropy of distribution of significant correlations

	Vector ret1 = new Vector(), ret2, v, w;
	Hashtable ht;
	int tot, test, cur, i;
	Enumeration extensions;
	String R;
	double contrib;

	double [] entropies = new double [all.size()];
	int [] indexes = new int [all.size()];

	for (i=0;i<all.size();i++){
	    v = (Vector) all.elementAt(i);
	    tot = ((Integer) v.elementAt(1)).intValue();
	    test = 0;
	    ht = (Hashtable) v.elementAt(3);
	    extensions = ht.keys();
	    contrib = 0.0;
	    cur = 0;

	    if(extensions != null) {
		while (extensions.hasMoreElements()) {
		    R = (String) extensions.nextElement();
		    cur = ((Integer) ht.get(R)).intValue();
		    
		    if (cur > 0){
			contrib += -( (double) cur / (double) tot) * Math.log((double) cur / (double) tot);
		    }

		    test += cur;
		}
	    }else
		Matrix.perror("extentions == null for " + (String) v.elementAt(0));
	    if (test != tot)
		Matrix.perror(" cur = " + cur + "; tot = " + tot + " . should be identical");
	    w = new Vector();
	    w.addElement((String) v.elementAt(0));
	    w.addElement(new Double(contrib));
	    ret1.addElement(w);
	    
	    indexes[i] = i;
	    entropies[i] = contrib;
	}

	QuickSort.quicksort(entropies, indexes);
	ret2 = new Vector();
	for (i=0;i<ret1.size();i++)
	    ret2.addElement((Vector) ret1.elementAt(indexes[i]));

	return ret2;
	}*/ // removed to be sure not to use and yet keep code

    public int sizeNotZero(Hashtable ht){
	int ret = 0, cur;
	Enumeration extensions = ht.keys();
	String R;

	if(extensions != null) {
	    while (extensions.hasMoreElements()) {
		R = (String) extensions.nextElement();
		cur = ((Integer) ht.get(R)).intValue();
		if (cur != 0)
		    ret++;
	    }
	}

	return ret;
    }

    public Vector orderedKullbackLeiblerCorrelationAsciiNamesAggregated(Vector all){
	// input = output of orderedCorrelationAsciiNamesAggregated
	// output = set of ordered KL divergences wrt uniform, increasing, giving for each element, a Vector
	// v = Species, KL div of distribution of significant correlations

	Vector ret1 = new Vector(), ret2, v, w;
	Hashtable ht;
	int tot, test, cur, i, count;
	Enumeration extensions;
	String R;
	double contrib;

	double [] kldiv = new double [all.size()];
	int [] indexes = new int [all.size()];

	for (i=0;i<all.size();i++){
	    v = (Vector) all.elementAt(i);
	    tot = ((Integer) v.elementAt(1)).intValue();
	    test = 0;
	    ht = (Hashtable) v.elementAt(3);
	    extensions = ht.keys();
	    contrib = 0.0;
	    cur = 0;
	    count = sizeNotZero(ht);

	    if(extensions != null) {
		while (extensions.hasMoreElements()) {
		    R = (String) extensions.nextElement();
		    cur = ((Integer) ht.get(R)).intValue();
		    
		    if (cur > 0){
			contrib += ( (double) cur / (double) tot) * Math.log( (double) count * (double) cur / (double) tot);
		    }

		    test += cur;
		}
	    }else
		Matrix.perror("extentions == null for " + (String) v.elementAt(0));
	    if (test != tot)
		Matrix.perror(" cur = " + cur + "; tot = " + tot + " . should be identical");
	    w = new Vector();
	    w.addElement((String) v.elementAt(0));
	    w.addElement(new Double(contrib));
	    ret1.addElement(w);
	    
	    indexes[i] = i;
	    kldiv[i] = -contrib;
	}

	QuickSort.quicksort(kldiv, indexes);
	ret2 = new Vector();
	for (i=0;i<ret1.size();i++)
	    ret2.addElement((Vector) ret1.elementAt(indexes[i]));

	return ret2;
    }

    public Vector orderedCorrelationAsciiNamesAggregated(Vector all){
	//Aggregates correlations by summarizing for each asciiNames the TOTAL #sig, and TOTAL #total over all asciiNames
	// returns a vector of Vectors, each element = asciiName, #Tsig, #Ttotal, Ht
	// ordered according to ratio Tsig/Ttotal
	// Ht is a Hashtable giving for each other species the number of sig (thus, the sum over Ht integers = #Tsig)
	Vector allNames = new Vector(), v, v2, ret = new Vector(), ret2;
	int i, j, sig, tot, num, den;
	double [] ratios;
	int [] indexes;
	String n1, n2;
	Hashtable ht;
	for (i=0;i<all.size();i++){
	    v = (Vector) all.elementAt(i);
	    n1 = (String) v.elementAt(0);
	    n2 = (String) v.elementAt(1);

	    if (!(allNames.contains(n1)))
		allNames.addElement(n1);
	    if (!(allNames.contains(n2)))
		allNames.addElement(n2);
	}
	for (i=0;i<allNames.size();i++){
	    v = new Vector();
	    v.addElement(allNames.elementAt(i));
	    v.addElement(new Integer(0));
	    v.addElement(new Integer(0));
	    v.addElement(new Hashtable());
	    ret.addElement(v);
	}

	ratios = new double [ret.size()];
	indexes = new int [ret.size()];

	for (i=0;i<all.size();i++){
	    v = (Vector) all.elementAt(i);
	    n1 = (String) v.elementAt(0);
	    n2 = (String) v.elementAt(1);
	    sig = ((Integer) v.elementAt(2)).intValue();
	    tot = ((Integer) v.elementAt(3)).intValue();
	    
	    j = allNames.indexOf(n1);
	    if (j==-1)
		Matrix.perror("AGCTFileWriter -- " + n1 + " not in asciiNames");
	    v2 = (Vector) ret.elementAt(j);
	    addTo(v2,1,sig);
	    addTo(v2,2,tot);

	    ht = (Hashtable) v2.elementAt(3);
	    if (ht.containsKey(n2))
		Matrix.perror("Key " + n2 + " already present in hashtable of " + n1);
	    ht.put(n2, new Integer(sig));

	    if (!n2.equals(n1)){
		j = allNames.indexOf(n2);
		if (j==-1)
		    Matrix.perror("AGCTFileWriter -- " + n2 + " not in asciiNames");
		v2 = (Vector) ret.elementAt(j);
		addTo(v2,1,sig);
		addTo(v2,2,tot);

		ht = (Hashtable) v2.elementAt(3);
		if (ht.containsKey(n1))
		    Matrix.perror("Key " + n1 + " already present in hashtable of " + n2);
		ht.put(n1, new Integer(sig));
	    }
	}

	for (i=0;i<ret.size();i++){
	    v2 = (Vector) ret.elementAt(i);
	    num = ((Integer) v2.elementAt(1)).intValue();
	    den = ((Integer) v2.elementAt(2)).intValue();
	    if (den == 0)
		ratios[i] = 1.0;
	    else
		ratios[i] = 1.0 - ((double) num) / ((double) den);
	    indexes[i] = i;
	}

	QuickSort.quicksort(ratios, indexes);
	ret2 = new Vector();
	for (i=0;i<ret.size();i++)
	    ret2.addElement((Vector) ret.elementAt(indexes[i]));

	return ret2;
    }


    public Hashtable correlationAsciiNames(boolean usePos, boolean useNeg){
	// returns statistics over asciiNames for couple of genes: given two same asciiNames, returns the number of P-significant Delaunay edges over the number of Delaunau edges, between two genes of these asciiNames
	// More precisely, each element is a Vector, with (0,1) the asciiNames in lexicographic order, (2) the number of sig and (3) the total number
	// "sig" is computed using usePos and useNeg <0 and >0 correlations

	Hashtable ret = new Hashtable();

	int i, j, nid;
	Gene gg, ggn;
	boolean useNeighbor;
	String key, n1, n2;
	Vector v;
	double nangle;

	for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
	    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
	    n1 = gg.asciiName;
				
	    if ( (n1 != null) && (n1 != "") && (gg.neighbors != null) && (gg.neighbors.size()>0) ){
		for (j=0;j<gg.neighbors.size();j++){
		    nid = ((Integer) gg.neighbors.elementAt(j)).intValue();
		    ggn = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[nid]);
		    n2 = ggn.asciiName;

		    if ( (n2 != null) && (n2 != "") ){
			nangle = ((Double) gg.neighbor_angles.elementAt(j)).doubleValue();
			if ( ( (usePos) && (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0) ) || ( (useNeg) && (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0))) ) )
			    useNeighbor = true;
			else
			    useNeighbor = false;
			key = keyForNames(n1, n2);
			if (ret.containsKey(key))
			    v = (Vector) ret.get(key);
			else{
			    v = new Vector();
			    if (n1.compareTo(n2) < 0){
				v.addElement(n1);
				v.addElement(n2);
			    }else{
				v.addElement(n2);
				v.addElement(n1);
			    }
			    v.addElement(new Integer(0));
			    v.addElement(new Integer(0));
			}
			increment(v,3);
			if (useNeighbor)
			    increment(v,2);

			ret.put(key,v);
		    }
		}
	    }
	}


	return ret;
    }

    public void toSavingDelaunay(final File rf){
	Thread t = new Thread(){
		public void run(){

		    FileWriter f = null;
		    try{
			int i, j, k, nid;
			double nangle;
			Gene gg, ggn, ggnn;
			boolean okNeighbor, okRef;
			Vector nn = null, ll, element;
			Character pos = new Character('+'), neg = new Character('-'), dumchar = null;
			
			f = new FileWriter(rf);
			
			f.write(separator());
			f.write("Saving the *Filtered* Delaunay Triangulation (P-Del = " + Statistics.LIMIT_P_DELAUNAY + ") for domain " + myAGCT.myDomain.domainName);
			f.write(separator());
			
			f.write ("\n\n");
			
			f.write(separator());
			f.write("Processing parameters:\n");
			f.write(myAGCT.getParameterString());
			f.write(separator()); 
			
			if (myAGCT.isTriangulated){
			    Vector ht, htel, htag, kldiv;
			    int num, den;
			    double ratio;

			    ht = orderedCorrelationAsciiNames(correlationAsciiNames(true, false));
			    if (ht != null){
				htag = orderedCorrelationAsciiNamesAggregated(ht);
				kldiv = orderedKullbackLeiblerCorrelationAsciiNamesAggregated(htag);

				f.write("[G000+] Kullback-Leibler divergences list, decreasing values, for (>0) correlations, given as: Name -- Entropy\n\n");
				for (i=0;i<kldiv.size();i++){
				    htel = (Vector) kldiv.elementAt(i);
				    
				    f.write((String) htel.elementAt(0) + " -- " + (Double) htel.elementAt(1) + "\n");
				}
				f.write("\n\n");

				f.write("[G001+] aggregated asciiNames Correlation list for (>0) correlations, given as: Name -- RATIO -- #TOTALsignificant_edges_for_P (#TOTALtotal_edges)\n\n");
				for (i=0;i<htag.size();i++){
				    htel = (Vector) htag.elementAt(i);
				    num = ((Integer) htel.elementAt(1)).intValue();
				    den = ((Integer) htel.elementAt(2)).intValue();
				    if (den > 0)
					ratio = ((double) num)/((double) den);
				    else
					ratio = 0.0;

				    f.write((String) htel.elementAt(0) + " -- " + ratio + " -- " + num + " (" + den + ")\n");
				}
				f.write("\n\n");

				f.write("[G002+] asciiNames x asciiNames Correlation list for (>0) correlations, given as: Name1_Name2 -- RATIO -- #significant_edges_for_P (#total_edges)\n\n");
				for (i=0;i<ht.size();i++){
				    htel = (Vector) ht.elementAt(i);
				    num = ((Integer) htel.elementAt(2)).intValue();
				    den = ((Integer) htel.elementAt(3)).intValue();
				    if (den > 0)
					ratio = ((double) num)/((double) den);
				    else
					ratio = 0.0;

				    f.write((String) htel.elementAt(0) + " " + (String) htel.elementAt(1) + " -- " + ratio + " -- " + num + " (" + den + ")\n");
				}
				f.write("\n\n");
			    }

			    ht = orderedCorrelationAsciiNames(correlationAsciiNames(false, true));
			    if (ht != null){
				htag = orderedCorrelationAsciiNamesAggregated(ht);
				kldiv = orderedKullbackLeiblerCorrelationAsciiNamesAggregated(htag);

				f.write("[G000-] Kullback-Leibler divergences list, decreasing values, for (<0) correlations, given as: Name -- Entropy\n\n");
				for (i=0;i<kldiv.size();i++){
				    htel = (Vector) kldiv.elementAt(i);
				    
				    f.write((String) htel.elementAt(0) + " -- " + (Double) htel.elementAt(1) + "\n");
				}
				f.write("\n\n");

				f.write("[G001-] aggregated asciiNames Correlation list for (<0) correlations, given as: Name -- RATIO -- #TOTALsignificant_edges_for_P (#TOTALtotal_edges)\n\n");
				for (i=0;i<htag.size();i++){
				    htel = (Vector) htag.elementAt(i);
				    num = ((Integer) htel.elementAt(1)).intValue();
				    den = ((Integer) htel.elementAt(2)).intValue();
				    if (den > 0)
					ratio = ((double) num)/((double) den);
				    else
					ratio = 0.0;

				    f.write((String) htel.elementAt(0) + " -- " + ratio + " -- " + num + " (" + den + ")\n");
				}
				f.write("\n\n");

				f.write("[G002-] asciiNames x asciiNames Correlation list for (<0) correlations, given as: Name1_Name2 -- RATIO -- #significant_edges_for_P (#total_edges)\n\n");
				for (i=0;i<ht.size();i++){
				    htel = (Vector) ht.elementAt(i);
				    num = ((Integer) htel.elementAt(2)).intValue();
				    den = ((Integer) htel.elementAt(3)).intValue();
				    if (den > 0)
					ratio = ((double) num)/((double) den);
				    else
					ratio = 0.0;

				    f.write((String) htel.elementAt(0) + " " + (String) htel.elementAt(1) + " -- " + ratio + " -- " + num + " (" + den + ")\n");
				}
				f.write("\n\n");
			    }

			    ht = orderedCorrelationAsciiNames(correlationAsciiNames(true, true));
			    if (ht != null){
				htag = orderedCorrelationAsciiNamesAggregated(ht);
				kldiv = orderedKullbackLeiblerCorrelationAsciiNamesAggregated(htag);

				f.write("[G000+-] Kullback-Leibler divergences list, decreasing values, for (>0 & <0) correlations, given as: Name -- Entropy\n\n");
				for (i=0;i<kldiv.size();i++){
				    htel = (Vector) kldiv.elementAt(i);
				    
				    f.write((String) htel.elementAt(0) + " -- " + (Double) htel.elementAt(1) + "\n");
				}
				f.write("\n\n");


				f.write("[G001+-] aggregated asciiNames Correlation list for (>0 & <0) correlations, given as: Name -- RATIO -- #TOTALsignificant_edges_for_P (#TOTALtotal_edges)\n\n");
				for (i=0;i<htag.size();i++){
				    htel = (Vector) htag.elementAt(i);
				    num = ((Integer) htel.elementAt(1)).intValue();
				    den = ((Integer) htel.elementAt(2)).intValue();
				    if (den > 0)
					ratio = ((double) num)/((double) den);
				    else
					ratio = 0.0;

				    f.write((String) htel.elementAt(0) + " -- " + ratio + " -- " + num + " (" + den + ")\n");
				}
				f.write("\n\n");

				f.write("[G002+-] asciiNames x asciiNames Correlation list for (>0 & <0) correlations, given as: Name1_Name2 -- RATIO -- #significant_edges_for_P (#total_edges)\n\n");
				for (i=0;i<ht.size();i++){
				    htel = (Vector) ht.elementAt(i);
				    num = ((Integer) htel.elementAt(2)).intValue();
				    den = ((Integer) htel.elementAt(3)).intValue();
				    if (den > 0)
					ratio = ((double) num)/((double) den);
				    else
					ratio = 0.0;

				    f.write((String) htel.elementAt(0) + " " + (String) htel.elementAt(1) + " -- " + ratio + " -- " + num + " (" + den + ")\n");
				}
				f.write("\n\n");
			    }

			    f.write("[G01] Unordered (Complete) graph, given as: Gene-Name_ID List-of-Delaunay-Neighbors (Gene-Name_IDS)\n");
			
			    AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Saving Complete Filtered Del. Tr.", myAGCT.myDomain.numberSelectedGenes);
			    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
				gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
				
				if ( (gg.neighbors != null) && (gg.neighbors.size()>0) ){
				    nn = new Vector();
				    for (j=0;j<gg.neighbors.size();j++){
					okNeighbor = false;
					nid = ((Integer) gg.neighbors.elementAt(j)).intValue();
					ggn = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[nid]);
					nangle = ((Double) gg.neighbor_angles.elementAt(j)).doubleValue();

					if ( (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0) || (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0))) )
					    okNeighbor = true;

					if (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0)
					    dumchar = pos;
					else if (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0)))
					    dumchar = neg;

					if (okNeighbor){
					    element = new Vector();
					    element.addElement(geneStringDelaunay(ggn, true) + "_(*)");
					    element.addElement(dumchar);

					    nn.addElement(element);
					    if (!Prototype.No_Reduction){
						ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.myDomain.selectedGeneNumberToGeneNumber[nid]));
						for (k=1;k<ll.size();k++){
						    ggnn = (Gene) myAGCT.myDomain.domainGenes.elementAt(((Integer) ll.elementAt(k)).intValue());
						    element = new Vector();
						    element.addElement(geneStringDelaunay(ggnn, false));
						    element.addElement(dumchar);

						    nn.addElement(element);
						}
					    }
					}
				    }
				    
				    if (nn.size()>0){
					f.write(geneStringDelaunay(gg, true) + "\t");
					for (k=0;k<nn.size();k++){
					    element = (Vector) nn.elementAt(k);
					    f.write((String) element.elementAt(0) + " (" + (Character) element.elementAt(1) + ")");
					    if (k<nn.size()-1)
						f.write("\t");
					    else
						f.write("\n");
					}
					
					if (!Prototype.No_Reduction){
					    ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]));
					    for (j=1;j<ll.size();j++){
						ggn = (Gene) myAGCT.myDomain.domainGenes.elementAt(((Integer) ll.elementAt(j)).intValue());
						f.write(geneStringDelaunay(ggn, false) + "\t");
						for (k=0;k<nn.size();k++){
						    element = (Vector) nn.elementAt(k);
						    f.write((String) element.elementAt(0) + " (" + (Character) element.elementAt(1) + ")");
						    if (k<nn.size()-1)
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
			
			    cc = new AGCTCounter(myAGCT.myInformationFrame, "Saving P + H Filtered Del. Tr.", myAGCT.myDomain.numberSelectedGenes);
			    for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
				gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
				
				if ( (gg.neighbors != null) && (gg.neighbors.size()>0) ){
				    nn = new Vector();
				    for (j=0;j<gg.neighbors.size();j++){
					okNeighbor = false;
					nid = ((Integer) gg.neighbors.elementAt(j)).intValue();
					ggn = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[nid]);
					nangle = ((Double) gg.neighbor_angles.elementAt(j)).doubleValue();

					if ( (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0) || (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0))) )
					    okNeighbor = true;

					if (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0)
					    dumchar = pos;
					else if (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0)))
					    dumchar = neg;

					if (okNeighbor){
					    element = new Vector();
					    element.addElement(geneStringDelaunay(ggn, true));
					    element.addElement(dumchar);
					    nn.addElement(element);
					    if (!Prototype.No_Reduction){
						ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.myDomain.selectedGeneNumberToGeneNumber[nid]));
						for (k=1;k<ll.size();k++){
						    ggnn = (Gene) myAGCT.myDomain.domainGenes.elementAt(((Integer) ll.elementAt(k)).intValue());
						    if ( (AGCT.Referenced_Available) && (ggnn.referenced) ){
							element = new Vector();
							element.addElement(geneStringDelaunay(ggnn, false));
							element.addElement(dumchar);
							nn.addElement(element);
						    }
						}
					    }
					}
				    }
				    
				    if (nn.size()>0){
					f.write(geneStringDelaunay(gg, true) + "\t");
					for (k=0;k<nn.size();k++){
					    element = (Vector) nn.elementAt(k);
					    f.write((String) element.elementAt(0) + " (" + (Character) element.elementAt(1) + ")");
					    if (k<nn.size()-1)
						f.write("\t");
					    else
						f.write("\n");
					}
					
					if (!Prototype.No_Reduction){
					    ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]));
					    for (j=1;j<ll.size();j++){
						ggn = (Gene) myAGCT.myDomain.domainGenes.elementAt(((Integer) ll.elementAt(j)).intValue());
						if ( (AGCT.Referenced_Available) && (ggn.referenced) ){
						    f.write(geneStringDelaunay(ggn, false) + "\t");
						    for (k=0;k<nn.size();k++){
							element = (Vector) nn.elementAt(k);
							f.write((String) element.elementAt(0) + " (" + (Character) element.elementAt(1) + ")");
							if (k<nn.size()-1)
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
			    if (AGCT.Referenced_Available){
				f.write("\n\n[G03] Unordered (Highlighted) graph, given as: Gene-Name_ID List-of-Delaunay-Neighbors (Gene-Name_IDS)\n");
				
				cc = new AGCTCounter(myAGCT.myInformationFrame, "Saving P + H Filtered Del. Tr.", myAGCT.myDomain.numberSelectedGenes);
				for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
				    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
				    
				    if ( (gg.neighbors != null) && (gg.neighbors.size()>0) ){
					nn = new Vector();
					for (j=0;j<gg.neighbors.size();j++){
					    okNeighbor = false;
					    nid = ((Integer) gg.neighbors.elementAt(j)).intValue();
					    ggn = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[nid]);
					    nangle = ((Double) gg.neighbor_angles.elementAt(j)).doubleValue();
					    
					    if ( (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0) || (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0))) )
						okNeighbor = true;
					    
					    if (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0)
						dumchar = pos;
					    else if (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0)))
						dumchar = neg;
					    
					    if (okNeighbor){
						element = new Vector();
						if (ggn.referenced){
						    element.addElement(geneStringDelaunay(ggn, true));
						    element.addElement(dumchar);
						    nn.addElement(element);
						}
						if (!Prototype.No_Reduction){
						    ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.myDomain.selectedGeneNumberToGeneNumber[nid]));
						    for (k=1;k<ll.size();k++){
							ggnn = (Gene) myAGCT.myDomain.domainGenes.elementAt(((Integer) ll.elementAt(k)).intValue());
							if ( (AGCT.Referenced_Available) && (ggnn.referenced) ){
							    element = new Vector();
							    element.addElement(geneStringDelaunay(ggnn, false));
							    element.addElement(dumchar);
							    nn.addElement(element);
							}
						    }
						}
					    }
					}
					
					if (nn.size()>0){
					    if (gg.referenced){
						f.write(geneStringDelaunay(gg, true) + "\t");
						for (k=0;k<nn.size();k++){
						    element = (Vector) nn.elementAt(k);
						    f.write((String) element.elementAt(0) + " (" + (Character) element.elementAt(1) + ")");
						    if (k<nn.size()-1)
							f.write("\t");
						    else
							f.write("\n");
						}
					    }
					    
					    if (!Prototype.No_Reduction){
						ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]));
						for (j=1;j<ll.size();j++){
						    ggn = (Gene) myAGCT.myDomain.domainGenes.elementAt(((Integer) ll.elementAt(j)).intValue());
						    if ( (AGCT.Referenced_Available) && (ggn.referenced) ){
							f.write(geneStringDelaunay(ggn, false) + "\t");
							for (k=0;k<nn.size();k++){
							    element = (Vector) nn.elementAt(k);
							    f.write((String) element.elementAt(0) + " (" + (Character) element.elementAt(1) + ")");
							    if (k<nn.size()-1)
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
		    }catch(IOException e){
			Matrix.perror("AGCTFileWriter.class :: Writing error for Filtered Delaunay Triangulation");
		    }
		}
	    };
	t.start();
    }

    public void toSavingWeb(final File rf){
	Thread t = new Thread(){
		public void run(){

		    FileWriter f = null;
		    try{
			int i, j, d;
			Gene gg, ggn;
			Vector ll;
			double [] compo;
			
			f = new FileWriter(rf);
			AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Saving to Webfile", myAGCT.myDomain.numberSelectedGenes);
			for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++){
			    gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
			    if (myAGCT.myTabbedPane.getSelectedIndex() == 1)
				compo = gg.manifold_Pnt3D.coordinates;
			    else
				compo = gg.pca_Pnt3D.coordinates;
			    
			    f.write(gg.name + "\t");
			    f.write(gg.asciiName + "\t");
			    for (d=0;d<3;d++){
				f.write(compo[d] + "");
				if (d<2)
				    f.write("\t");
			    }
			    f.write("\n");

			    if (!Prototype.No_Reduction){
				ll = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]));
				for (j=1;j<ll.size();j++){
				    ggn = (Gene) myAGCT.myDomain.domainGenes.elementAt(((Integer) ll.elementAt(j)).intValue());
				    f.write(ggn.name + "\t");
				    f.write(ggn.asciiName + "\t");
				    for (d=0;d<3;d++){
					f.write(compo[d] + "");
					if (d<2)
					    f.write("\t");
				    }
				    f.write("\n");
				}
			    }
			    cc.increment();
			}
			cc.end();
			f.close();
		    }catch(IOException e){
			Matrix.perror("AGCTFileWriter.class :: Writing error for Web file");
		    }
		}
	    };
	t.start();
    }

    
    public String separator(){
	return "\n----------------------------------------------------------------------------\n";
    }
}
