import java.util.*;
import java.lang.*;
import java.io.*;

class Domain implements Debuggable{
    static final String ANNOTATION_F = new String("F:"), 
	ANNOTATION_P = new String("P:"), 
	ANNOTATION_C = new String("C:"),
	APPEND_NAME = "__";

    AGCT myAGCT;

    String domainName;

    String geneAnnotationsName;
    int maxHierarchy;

    /***** Selection of prototypes, ligands, features, etc *****/
    //Selection and reduction file;
    String selectionFile;
    boolean bSelectionFile;

    /***** LIGANDS *****/

    Vector domainLigands;
    //each entry = ( (String = ligand name), (X Strings = ligand groups) )
    Hashtable accessToLigands;
    //keys are ligand names (String), return the index as Integer in domainLigand for the ligand data
    int numberLigands;

    /***** GROUPS *****/

    Vector domainGroups;
    //each entry = the name of a group as String, in the order as they are built for each ligand
    Vector [] domainGroupElements;
    //each entry = a Vector giving the set of different String for group members
    // ex: for AMPC-CA2+, gives "YN", "NN", "NY", ...
    Hashtable accessToGroups;
    //keys are group names (String), return the index as Integer in domainGroups and domainGroupElements
    int numberGroups;

    /***** TIMES *****/

    Vector domainTimes;
    //each entry=( (String = ligand name), (X Double = time stamps for the ligand) )
    //Order of ligands MUST follow that in domainLigands (to use accurately accessToLigands)
    Vector numberTimes;
    //each entry = an Integer giving the # of time steps for ligand at the same index in domainLigand
    Vector numberOfTimeStamps;
    //each entry = an Integer giving the number of stamps that follow for the features in ligand index

    /***** FEATURES *****/

    Vector domainFeatures;
    //entries are feature names (typically, concatenation of a ligand + a time)
    Hashtable accessToFeatures;
    //keys are the feature names (String), return the index as Integer in featureNames for the feature
    int numberFeatures;

    /***** GENES *****/

    //genes variables
    int numberInitialGenes;
    Vector domainGenes;
    //each entry = an instance from class Gene

    Hashtable accessToGenesIndexBefore;
    //keys are Gene names (String), returns the index of Gene as Integer in domainGenes
    Hashtable accessToGenesObjectBefore;
    //keys are Gene indexes (Integer), returns the Gene in domainGenes

    Hashtable accessToGenesObject;
    //keys are Gene indexes (Integer), returns the Gene in domainGenes
    //this is after prototype selection

    /***** FILTERING VARIABLES *****/

    //Before means before Prototype Selection

    //Booleans for the filtering step
    boolean [] selectedGeneBefore;
    //is gene selected ?
    boolean [] enabledGeneBefore;
    //is gene available for selection ?

    boolean [][] selectedGroup;
    //is (group / ligand) selected ?

    //Following : Selected Genes BEFORE the selection of Prototypes

    int numberSelectedGenesBefore;
    int [] selectedGeneNumberToGeneNumberBefore;
    //maps the index of the selected to the index of the corresponding gene in domainGenes
    Hashtable accessToSelectedGeneNumberToGeneNumberBefore;
    //does the same as selectedGeneNumberToGeneNumberBefore, but with a Hashtable

    //Following : Prototype defined variables (Gene name kept for simplicity)

    int numberSelectedGenes;
    int [] selectedGeneNumberToGeneNumber;
    //maps the index of the selected to the index of the corresponding gene in domainGenes
    Hashtable accessToSelectedGeneNumberToGeneNumber;
    //does the same as selectedGeneNumberToGeneNumber, but with a Hashtable

    // A CHANGER : : avec Prototypes

    Vector selectedGenesAllAnnotationsF, selectedGenesAllAnnotationsP, selectedGenesAllAnnotationsC;
    Hashtable accessToSelectedGenesAllAnnotationsF, accessToSelectedGenesAllAnnotationsP, accessToSelectedGenesAllAnnotationsC;
    //Keys are annotations strings, returns the index as Integer for the feature in the list of checkboxes of the JAnnotationFrame

    Hashtable annotationFToGenes, annotationPToGenes, annotationCToGenes;
    //Keys are annotations, returns an unordered list (Vector) of gene indexes having the annotation

    boolean someAnnotationsF, someAnnotationsP, someAnnotationsC;

    boolean [] selectedLigand;
    //is ligand (in domainLigand) selected ?
    boolean [] enabledLigand;
    //is ligand (in domainLigand) available for selection ?

    int numberSelectedLigands;

    Pnt3D min_ManifoldPnt3D, min_PCAPnt3D;
    Pnt3D max_ManifoldPnt3D, max_PCAPnt3D;

    /***** TRIANGULATION EDGES *****/

    //List of couples of genes
    DelaunayTriangulation myDelaunayTriangulation;
    Vector myEdges;

    /***** Highlighted genes *****/
    Hashtable referencesToInteger;
    Vector orderedReferences;

    public Domain(File rf, AGCT ap){
	AGCT.Referenced_Available = false;
	referencesToInteger = null;
	orderedReferences = null;

	selectedGenesAllAnnotationsF = selectedGenesAllAnnotationsP = selectedGenesAllAnnotationsC = null;
	someAnnotationsF = someAnnotationsP = someAnnotationsC = false;
	accessToSelectedGenesAllAnnotationsF = accessToSelectedGenesAllAnnotationsP = accessToSelectedGenesAllAnnotationsC = null;

	geneAnnotationsName = "";
	selectionFile = "";
	bSelectionFile = false;
	maxHierarchy = -1;

	myAGCT = ap;

	min_ManifoldPnt3D = new Pnt3D();
	max_ManifoldPnt3D = new Pnt3D();

	min_PCAPnt3D = new Pnt3D();
	max_PCAPnt3D = new Pnt3D();

	FileReader e;
	StringTokenizer t;
	String s, currentStep = "", valret;
	BufferedReader br;
	int i, j, k, kk, m, overallTimes, total = 0;
	Vector v, v2, vt;
	boolean collectData = false;
	Gene g;
	Vector distinctNames_list = new Vector();
	Vector distinctNames_counter = new Vector();

	valret = "File opened: " + rf.getName() + " --- statistics:\n";

	//Counting # ligands, groups
	numberLigands = 0;
	numberGroups = 0;
	numberTimes = new Vector();
	overallTimes = 0;
	accessToLigands = new Hashtable();
	domainLigands = new Vector();
	domainGroups = new Vector();
	accessToGroups = new Hashtable();
	domainTimes = new Vector();
	try{
	    e = new FileReader(rf);
	}catch(FileNotFoundException ex){
	    ap.myInformationFrame.setText("The file you try to load does not exist");
	    ap.rawDomainExists = false;
	    return;
	}
	br = new BufferedReader(e);
	try{
	    while ( (s=br.readLine()) != null){
		t = new StringTokenizer(s.replace('\t',' '), " ");
		if (t.countTokens()>0) {
		    s = t.nextToken();
		    if ( (s.length()>1) && (!s.substring(0,AGCTFileWriter.DATA_Comments.length()).equals(AGCTFileWriter.DATA_Comments)) ){
			if (s.equals(AGCTFileWriter.DATA_Ligand)){
			    if (numberGroups == 0)
				Matrix.perror("Groups not computed before the first @LIGAND token");
			    
			    s = t.nextToken();
			    domainLigands.addElement(new Vector());
			    domainTimes.addElement(new Vector());
			    accessToLigands.put(s,new Integer(numberLigands));
			    i = getLigandIndex(s);
			    v = (Vector) domainLigands.elementAt(i);
			    vt = (Vector) domainTimes.elementAt(i);
			    
			    v.addElement(s);
			    vt.addElement(s);
			    for (j=0;j<numberGroups;j++){
				s = t.nextToken();
				v.addElement(s);
				v2 = (Vector) domainGroupElements[j];
				if (!v2.contains(s))
				    v2.addElement(s);
			    }
			    j=0;
			    if (!t.hasMoreTokens())
				Matrix.perror("Ligand " + (String) ((Vector) domainLigands.elementAt(i)).elementAt(0) + " contains no time stamp");
			    while (t.hasMoreTokens()){
				s = t.nextToken();
				vt.addElement(new Double(Double.parseDouble(s)));
				j++;
				overallTimes ++;
			    }
			    numberTimes.addElement(new Integer(j));
			    numberLigands++; 
			}else if (s.equals(AGCTFileWriter.DATA_Groups)){
			    while (t.hasMoreTokens()){
				s = t.nextToken();
				domainGroups.addElement(s);
				accessToGroups.put(s,new Integer(numberGroups));
				numberGroups++;
			    }
			    domainGroupElements = new Vector[numberGroups];
			    for (i=0;i<numberGroups;i++)
				domainGroupElements[i] = new Vector();
			}else if (s.equals(AGCTFileWriter.DATA_Domain)){
			    if (t.hasMoreTokens()){
				s = t.nextToken();
				domainName = s;
			    }else
				domainName = "[No domaine name]";
			}else if (s.equals(AGCTFileWriter.DATA_GeneAnnotations))
			    geneAnnotationsName = t.nextToken();
			else if (s.equals(AGCTFileWriter.DATA_Selection)){
			    setSelectionFile(t.nextToken());
			}
		    }
		}
	    }
	    e.close();
	}catch(IOException ex){
	    ap.myInformationFrame.setText("IOError when loading domain");
	    ap.rawDomainExists = false;
	    return;
	}

	e = null;

	if (numberLigands == 0)
	    Matrix.perror("No ligand counted. Perhaps a coding problem for your file ?");

	if (overallTimes == 0)
	    Matrix.perror("No time stamps counted. Perhaps a coding problem for your file ?");

	//compute feature names
	accessToFeatures = new Hashtable();
	domainFeatures = new Vector();
	numberFeatures=0;
	for (i=0;i<numberLigands;i++)
	    for (j=0;j<( (Integer) numberTimes.elementAt(i)).intValue();j++){
		s = getLigandName(i) + (Double) ((Vector) domainTimes.elementAt(i)).elementAt(j+1);
		accessToFeatures.put(s,new Integer(numberFeatures));
		domainFeatures.addElement(s);
		numberFeatures++;
	    }

	//computes the genes raw components: reopens the file
	numberInitialGenes = 0;
	domainGenes = new Vector();
	accessToGenesIndexBefore = new Hashtable();

	boolean oneUndetermined = false;
	try{
	    e = new FileReader(rf);
	}catch(FileNotFoundException ex){
	    ap.myInformationFrame.setText("The file you try to load does not exist");
	    ap.rawDomainExists = false;
	    return;
	}
	br = new BufferedReader(e);

	double [] vx;
	double [] vy;
	double [] vhaty;
	Vector tstamps;

	try{
	    while ( (s=br.readLine()) != null){
		total++;

		//System.out.println(s);
		t = new StringTokenizer(s.replace('\t',' '), " ");
		if (t.countTokens()>0) {
		    s = t.nextToken();
		    if (s.equals(AGCTFileWriter.DATA_Data))
			collectData = true;
		    else{
			if (collectData == true){

			    if (AGCT.MULTIPLE_ROWS_FOR_A_SAME_GENE_CONSIDERED_DIFFERENT_GENES)
				s = whichName(s, distinctNames_list, distinctNames_counter);
			    if ( (accessToGenesIndexBefore.containsKey(s)) && (!AGCT.MULTIPLE_ROWS_FOR_A_SAME_GENE_DELETED_EXCEPT_FIRST_ONE) )
				Matrix.perror("Domain.class :: Same gene name (" + s + ") found in two distinct data rows");

			    if (accessToGenesIndexBefore.containsKey(s))
				System.out.println("removing copy of gene" + s);

			    if (!accessToGenesIndexBefore.containsKey(s)){
				g = new Gene(s,numberLigands);
				accessToGenesIndexBefore.put(s,new Integer(numberInitialGenes));
				numberInitialGenes++;
			

				for (i=0;i<numberLigands;i++){
				    g.rawCoordinates[i] = new double [( (Integer) numberTimes.elementAt(i)).intValue()];
				    g.undeterminedCoordinates[i] = new boolean [( (Integer) numberTimes.elementAt(i)).intValue()];
				    for (j=0;j<( (Integer) numberTimes.elementAt(i)).intValue();j++){
					s = t.nextToken();
					
					if (s.equals(AGCTFileWriter.DATA_QuestionMark)){
					    g.rawCoordinates[i][j] = -1;
					    g.undeterminedCoordinates[i][j] = true;
					    oneUndetermined = true;
					}else{
					    g.rawCoordinates[i][j] = Double.parseDouble(s);
					    g.undeterminedCoordinates[i][j] = false;
					}
				    }

				    if (AGCT.Correction_Profiles_By_Average){
					//corrects with expected growth LinReg

					int nnn = ( (Integer) numberTimes.elementAt(i)).intValue();

					vx = new double [nnn];
					vy = new double [nnn];

					tstamps = (Vector) domainTimes.elementAt(i);
					for (j=0;j<nnn;j++){
					    vx[j] = ((Double) tstamps.elementAt(j+1)).doubleValue();
					    vy[j] = g.rawCoordinates[i][j];
					}

					vhaty = Statistics.LinReg(vx, vy);

					for (j=0;j<nnn;j++){
					    g.rawCoordinates[i][j] -= vhaty[j];
					}
					

					/****
					//XXX
					//for (j=0;j<nnn;j++)
					//   System.out.print(g.rawCoordinates[i][j] + " ");
					//System.out.println("\n");

					double average = 0.0;
					if (oneUndetermined)
					    Matrix.perror("Cannot correct with distances :: one value not specified in series");
					average = (g.rawCoordinates[i][nnn-1] - g.rawCoordinates[i][0]) / ( (double) nnn - 1 );

					for (j=1;j<nnn;j++){
					    g.rawCoordinates[i][j] -= ( (average) * ( (double) j ) );
					}
					
					//XXX
					//for (j=0;j<nnn;j++)
					//    System.out.print(g.rawCoordinates[i][j] + " ");
					//System.out.println("\n");
					****/
					
				    }
				}
				if (t.hasMoreTokens()){
				    s = t.nextToken();
				    if (!s.equals(AGCTFileWriter.DATA_Annotation))
					Matrix.perror("Annotation for gene " + g.name + " not specified as it should be.");
				    else{
					while (t.hasMoreTokens()){
					    s = t.nextToken();
					    g.addLightAnnotation(s + " ");
					}
				    }
				}
				domainGenes.addElement(g);
			    }
			}
		    }
		}
	    }
	    e.close();
	}catch(IOException ex){
	    ap.myInformationFrame.setText("IOError when loading domain");
	    ap.rawDomainExists = false;
	    return;
	}

	valret += "\nDomain " + domainName + ": \n";

	if (numberGroups > 0){
	    valret += "#Groups: " + numberGroups + " --- Group names: ";
	    for (i=0;i<numberGroups;i++){
		valret += domainGroups.elementAt(i);
		if (i<numberGroups - 1)
		    valret += ", ";
	    }
	    valret += "\n";

	    valret += "Groups composition:";
	    for (i=0;i<numberGroups;i++){
		valret += "Group " + domainGroups.elementAt(i) + " [" + domainGroupElements[i].size() + "] : ";
		for (j=0;j<domainGroupElements[i].size();j++){
		    valret += domainGroupElements[i].elementAt(j);
		    if (j<domainGroupElements[i].size()-1)
			valret += ", ";
		}

		if (i<numberGroups-1)
		    valret += " | ";
		else
		    valret += "\n";
	    }
	}

	valret += "#Ligands: " + numberLigands + " --- Ligand names [Groups: Times]: ";
	for (i=0;i<numberLigands;i++){
	    valret += getLigandName(i) + " ";
	    for (j=0;j<numberGroups;j++){
		if (j==0)
		    valret += "[";
		valret += getLigandGroupName(i,j);
		if (j < numberGroups - 1)
		    valret += ", ";
	    }
	    valret += ": ";
	    for (j=0;j<( (Integer) numberTimes.elementAt(i)).intValue();j++){
		valret += (Double) ((Vector) domainTimes.elementAt(i)).elementAt(j+1);
		if (j < ( (Integer) numberTimes.elementAt(i)).intValue() - 1)
		    valret += ", ";
	    }
	    valret += "] ";
	    if (i<numberLigands-1)
		valret += ", ";
	    else
		valret += "\n";
	}

	valret += "#Features: " + numberFeatures + " --- Feature names: ";
	for (i=0;i<numberFeatures;i++){
	    valret += domainFeatures.elementAt(i);
	    if (i<numberFeatures - 1)
		valret += ", ";
	}
	valret += "\n";

	valret += "#Genes: " + numberInitialGenes + " --- First genes: ";
	m = numberInitialGenes;
	if (m>3)
	    m = 3;
	kk=0;
	for (i=0;i<m;i++){
	    g = (Gene) domainGenes.elementAt(i);
	    valret += g.toString();
	}

	myAGCT.setInfoText(valret);
	
	if (geneAnnotationsName != "")
	    try{
		loadAnnotations(rf.getCanonicalPath().replaceAll(rf.getName(),""));
	    }catch(IOException ee){
		System.out.println("No canonical path for the file !");
	    }
	else
	    ap.annotationsExist = ap.annotationsExistP = ap.annotationsExistF = ap.annotationsExistC = false;

	ControlProcess.put("annotationsLoaded",true);


	if (!AGCT.Referenced_Available){
	    myAGCT.onlyReferenced.setEnabled(false);
	    myAGCT.onlyReferencedEdges.setEnabled(false);
	}else{
	    myAGCT.onlyReferenced.setEnabled(true);
	    myAGCT.onlyReferencedEdges.setEnabled(true);
	}
	myAGCT.loadHighlightFile.setEnabled(true);

	distinctNames_list = null;
	distinctNames_counter = null;
	//Prepare statistics
	Statistics.flushAll(ap);
    }

    public String whichName(String n, Vector dnl, Vector dnc){
	//Makes sure that the name does not exist; otherwise, append a special tag to distinguish it
	int ind, ck;
	String ret = n;

	if (dnl.contains(n)){
	    ind = dnl.indexOf(n);
	    ck = ( (Integer) dnc.elementAt(ind)).intValue();
	    ret += Domain.APPEND_NAME + ck;
	    dnc.setElementAt(new Integer(ck + 1), ind);
	}else{
	    dnl.addElement(ret);
	    dnc.addElement(new Integer(1));
	}
	return ret;
    }

    public void checkNames(){
	//checks that each gene name is different from the others
	boolean okProceed = true;
	String allNames = "";

	if (myAGCT.myTabbedPane.mySelectionPane.keepPrototypes.isSelected()){
	    int i, j;
	    String si, sj;		    
	    //AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Checking genes", numberInitialGenes * (numberInitialGenes - 1) / 2);
	    for (i=0;i<numberInitialGenes-1;i++){
		//cc.setText("Checking genes");
		si = ( (Gene) domainGenes.elementAt(i)).name;
		for (j=i+1;j<numberInitialGenes;j++){
		    sj = ( (Gene) domainGenes.elementAt(j)).name;
		    if (si.equals(sj)){
			okProceed = false;
			allNames += si + "\n";
		    }
		    //cc.increment();
		}
	    }
	    //cc.end();
	}
	if (!okProceed)
	    Matrix.perror("Domain.class :: some genes have the same name: " + allNames + "Stopping.");
    }

    public void setSelectionFile(String t){
	if (!t.equals(new String(""))){
	    selectionFile = t;
	    bSelectionFile = true;
	}
    }

    public String replace(String r){
	int i, j;
	char c;
	String ret = "";
	boolean found;
	i = 0;
	do{
	    j=0;
	    found = false;
	    do{
		if ((r.length()-i >= AGCTFileWriter.REPLACE_Old_New[j][0].length()) && (r.substring(i,i+AGCTFileWriter.REPLACE_Old_New[j][0].length()).equals(AGCTFileWriter.REPLACE_Old_New[j][0])))
		    found = true;
		else
		    j++;
	    }while ((!found) && (j<AGCTFileWriter.REPLACE_Old_New.length));
	    if (!found){
		ret += r.substring(i,i+1);
		i ++;
	    }else{
		ret += AGCTFileWriter.REPLACE_Old_New[j][1];
		i += AGCTFileWriter.REPLACE_Old_New[j][0].length();
	    }
	}while(i<r.length());
	return ret;
    }

    public void loadGeneAnnotationGEO_Format(String sRef, Vector varTypes, String eString){
	StringTokenizer t = new StringTokenizer(sRef, "\t"), dumT;
	String gID = "", gName = "", gF = "", gP = "", gC = "", gRef = "", dumS, refS;
	boolean readingF = false, readingP = false, readingC = false, readRef, bailout = false;
	Gene gg;
	int i, j, refe;
	String [] totalS1;
	String [] totalS2;

	if (t.countTokens() != varTypes.size())
	    Matrix.perror("Domain.class :: " + t.countTokens() + " tokens, but " + varTypes.size() + " variables on line \n" + sRef + "\n(perhaps empty lines in .agct file ?)");
	i = 0;
	while(t.hasMoreTokens()){
	    dumS = t.nextToken();
	    refS = (String) varTypes.elementAt(i);
	    if (refS.equals(AGCTFileWriter.ANNO_Variable_Definition_ID))
		gID = dumS;
	    else if (refS.equals(AGCTFileWriter.ANNO_Variable_Definition_Ascii_Names))
		gName = dumS;
	    else if (refS.equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_F))
		gF = dumS;
	    else if (refS.equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_P))
		gP = dumS;
	    else if (refS.equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_C))
		gC = dumS;
	    else if (refS.equals(AGCTFileWriter.ANNO_Variable_Definition_Referenced))
		gRef = dumS;
	    i++;
	}

	if ( ( (gID.equals(new String(""))) || (gID.equals(eString)) ) && (AGCT.VERBOSE) ){
	    //Matrix.perror("Domain.class :: gene ID, " + gID + " (name = " + gName + ", file string = " + sRef + "), is empty or dummy");
	    //System.out.println("Domain.class :: gene ID, " + gID + " (name = " + gName + ", file string = " + sRef + "), is empty or dummy");
	    System.out.println("Found empty or dummy ID in GEO format");
	    bailout = true;
	}
	if ( (!bailout) && (accessToGenesIndexBefore.containsKey(gID)) ){
	    gg = (Gene) domainGenes.elementAt(((Integer) accessToGenesIndexBefore.get(gID)).intValue());

	    if (!gRef.equals(new String(""))){
		refe = Integer.parseInt(gRef);
		if (refe > -1){
		    gg.putReferenced();
		    gg.typeReferenced(refe);
		    if (refe > JAGCTGraphicsPanel.referencedColors.length)
			Matrix.perror("Domain.class :: referenced number " + refe + " > " + JAGCTGraphicsPanel.referencedColors.length);
		    AGCT.Referenced_Available = true;
		    if (refe > AGCT.Max_Reference)
			AGCT.Max_Reference = refe;
		}
	    }

	    if (!gg.annotationsFilled){
		// F
		if (!gF.equals(eString)){
		    totalS1 = gF.split(AGCTFileWriter.REPLACE_Old_New[0][0]);
		    for (i=0;i<totalS1.length;i++){
			totalS2 = totalS1[i].split(AGCTFileWriter.REPLACE_Old_New[1][0]);
			if (totalS2.length != 3)
			    Matrix.perror("Domain.class :: bad number of tokens for F in line\n" + sRef);
			dumS = totalS2[1]; // take the central token
			readingF = true;
			readingP = readingC = false;
			gg.addAnnotation(readingF, readingP, readingC, dumS.replace(' ','_'));
			myAGCT.annotationsExistF = true;
		    }
		}
		
		// P
		if (!gP.equals(eString)){
		    totalS1 = gP.split(AGCTFileWriter.REPLACE_Old_New[0][0]);
		    for (i=0;i<totalS1.length;i++){
			totalS2 = totalS1[i].split(AGCTFileWriter.REPLACE_Old_New[1][0]);
			if (totalS2.length != 3)
			    Matrix.perror("Domain.class :: bad number of tokens for P in line\n" + sRef);
			dumS = totalS2[1]; // take the central token
			readingP = true;
			readingF = readingC = false;
			gg.addAnnotation(readingF, readingP, readingC, dumS.replace(' ','_'));
			myAGCT.annotationsExistP = true;
		    }
		}
		
		// C
		if (!gC.equals(eString)){
		    totalS1 = gC.split(AGCTFileWriter.REPLACE_Old_New[0][0]);
		    for (i=0;i<totalS1.length;i++){
			totalS2 = totalS1[i].split(AGCTFileWriter.REPLACE_Old_New[1][0]);
			if (totalS2.length != 3)
			    Matrix.perror("Domain.class :: bad number of tokens for C in line\n" + sRef);
			dumS = totalS2[1]; // take the central token
			readingC = true;
			readingF = readingP = false;
			gg.addAnnotation(readingF, readingP, readingC, dumS.replace(' ','_'));
			myAGCT.annotationsExistC = true;
		    }
		}
		
		//Name
		if (!gName.equals(eString))
		    gg.asciiName = gName;
		
		if (gg.annotationsF != null)
		    removeDuplicateStringTables(gg.annotationsF);
		if (gg.annotationsP != null)
		    removeDuplicateStringTables(gg.annotationsP);
		if (gg.annotationsC != null)
		    removeDuplicateStringTables(gg.annotationsC);
		
		if ( (gg.annotationsF != null) && (gg.annotationsF.size() > maxHierarchy) )
		    maxHierarchy = gg.annotationsF.size();
		if ( (gg.annotationsP != null) && (gg.annotationsP.size() > maxHierarchy) )
		    maxHierarchy = gg.annotationsP.size();
		if ( (gg.annotationsC != null) && (gg.annotationsC.size() > maxHierarchy) )
		    maxHierarchy = gg.annotationsC.size();

		gg.annotationsFilled = true;
	    }
	}
    }

    public void removeDuplicateStringTables(Vector v){
	int i = 0, j;
	String ti, tj;
	while(i<v.size()-1){
	    j=i+1;
	    ti = ((String[])v.elementAt(i))[1];
	    do{
		tj = ((String[])v.elementAt(j))[1];
		if (ti.equals(tj)){
		    v.removeElementAt(j);
		}else
		    j++;
	    }while(j<v.size());
	    i++;
	}
    }

    public void loadGeneAnnotationAGCT_Format(String sRef){
	StringTokenizer t = new StringTokenizer(sRef.replace('\t',' '), " ");
	String s;
	boolean readingF = false, readingP = false, readingC = false;
	Gene gg;

	if (t.countTokens()>0) {
	    s = t.nextToken();
	    if ( (s.length() < AGCTFileWriter.DATA_Comments.length()) || (!s.substring(0,AGCTFileWriter.DATA_Comments.length()).equals(AGCTFileWriter.DATA_Comments)) ){
		if (accessToGenesIndexBefore.containsKey(s)){
		    readingF = readingP = readingC = false;
		    gg = (Gene) domainGenes.elementAt(((Integer) accessToGenesIndexBefore.get(s)).intValue());
		    while (t.hasMoreTokens()){
			s = t.nextToken();
			if (s.equals(Domain.ANNOTATION_F))
			    readingF = true;
			else if (s.equals(Domain.ANNOTATION_P))
			    readingP = true;
			else if (s.equals(Domain.ANNOTATION_C))
			    readingC = true;
			else{
			    if (s != ""){
				if (readingF || readingP || readingC){
				    if ( (readingF && readingP)
					 || (readingF && readingC)
					 || (readingP && readingC) )
					System.out.print(" Gene " + gg.name);
				    
				    gg.addAnnotation(readingF, readingP, readingC, s);
				    if (readingF)
					myAGCT.annotationsExistF = true;
				    else if (readingP)
					myAGCT.annotationsExistP = true;
				    else if (readingC)
					myAGCT.annotationsExistC = true;
				    
				    readingF = readingP = readingC = false;
				    gg.annotationsFilled = true;
				}else{
				    if (gg.asciiName == null)
					gg.asciiName = s;
				}
			    }
			}
			
			if ( (gg.annotationsF != null) && (gg.annotationsF.size() > maxHierarchy) )
			    maxHierarchy = gg.annotationsF.size();
			if ( (gg.annotationsP != null) && (gg.annotationsP.size() > maxHierarchy) )
			    maxHierarchy = gg.annotationsP.size();
			if ( (gg.annotationsC != null) && (gg.annotationsC.size() > maxHierarchy) )
			    maxHierarchy = gg.annotationsC.size();

		    }
		}
	    }
	}
    }

    public Vector highlightAverageDimension;
    // average dimension per highlight tag (if avg dimension computed)

    public void loadHighlights(File rf){
	boolean someHighlights = false;
	FileReader e;
	int i, nfound = 0;
	Gene gg;
	BufferedReader br;
	String s, sl, sr, fs;
	boolean trouve;
	StringTokenizer t;
	Integer il, ir, is = null;
	Enumeration extensions;
	double ddim, cdim;

	referencesToInteger = null;
	orderedReferences = null;

	//erase
	AGCT.Referenced_Available = false;
	AGCT.Max_Reference = -1;
	myAGCT.onlyReferenced.setEnabled(false);
	myAGCT.onlyReferencedEdges.setEnabled(false);
	for (i=0;i<numberInitialGenes;i++)
	    if (selectedGeneBefore[i]){
		gg = (Gene) domainGenes.elementAt(i);
		gg.referenced = false;
		gg.typeReferenced = -1;
	    }

	try{
	    e = new FileReader(rf);
	}catch(FileNotFoundException ex){
	    myAGCT.myInformationFrame.setText("The file you try to load does not exist");
	    return;
	}

	Vector highlightAverageDimension_card = null;
	if (AGCT.LOCAL_DIMENSION_COMPUTED){
	    highlightAverageDimension = new Vector();
	    highlightAverageDimension_card = new Vector();
	}

	br = new BufferedReader(e);
	try{
	    while ( (s=br.readLine()) != null){
		if ( (s.length()>1) && (!s.substring(0,AGCTFileWriter.DATA_Comments.length()).equals(AGCTFileWriter.DATA_Comments)) ){
		    t = new StringTokenizer(s, ",");
		    if (t.countTokens()>0) {
			if (referencesToInteger == null)
			    referencesToInteger = new Hashtable();
			sl = t.nextToken();
			sr = t.nextToken();
			trouve = false;
			i=0;
			if (orderedReferences == null)
			    orderedReferences = new Vector();

			if (orderedReferences.size()>0){
			    do{
				if (sr.equals((String) orderedReferences.elementAt(i)))
				    trouve = true;
				else
				    i++;
			    }while( (!trouve) && (i<orderedReferences.size()) );
			}
			if (!trouve)
			    orderedReferences.addElement(sr);

			referencesToInteger.put(sl, new Integer(i));
		    }
		}
	    }
	    e.close();
	}catch(IOException ex){
	    myAGCT.myInformationFrame.setText("IOError when loading highlights");
	    return;
	}
	e = null;

	if (referencesToInteger != null){
	    if (AGCT.LOCAL_DIMENSION_COMPUTED){
		for (i=0;i<orderedReferences.size();i++){
		    highlightAverageDimension.addElement(new Double(0.0));
		    highlightAverageDimension_card.addElement(new Double(0.0));
		}
	    }

	    for (i=0;i<numberInitialGenes;i++)
		if (selectedGeneBefore[i]){
		    gg = (Gene) domainGenes.elementAt(i);
		    sl = gg.name;
		    sr = gg.asciiName;

		    il = (Integer) referencesToInteger.get(sl);		    
		    ir = null;
		    if (sr != null)
			ir = (Integer) referencesToInteger.get(sr);

		    //if ( (il == null) && (ir == null) )
		    //	Matrix.perror("Domain.class :: both il and ir are null");
		    if ( (il != null) || (ir != null) ){
			if (ir != null)
			    is = ir;
			else
			    is = il;
			gg.referenced = true;
			gg.typeReferenced = is.intValue();
			someHighlights = true;

			if (AGCT.LOCAL_DIMENSION_COMPUTED){
			    if (gg.averageDimension > 0.0){
				ddim = ((Double) highlightAverageDimension.elementAt(is.intValue())).doubleValue();
				cdim = ((Double) highlightAverageDimension_card.elementAt(is.intValue())).doubleValue();
				
				cdim += 1.0;
				ddim += (1.0 / gg.averageDimension);
				
				//System.out.println("(" + orderedReferences.elementAt(is.intValue()) + ") -- (" + sl + " / " + sr + ") " + gg.averageDimension + " -> ddim = " + ddim + ", cdim = " + cdim);
				
				highlightAverageDimension.setElementAt(new Double(ddim),is.intValue());
				highlightAverageDimension_card.setElementAt(new Double(cdim),is.intValue());
			    }
			}
		    }
		}
	}

	if (someHighlights){

	    if (AGCT.LOCAL_DIMENSION_COMPUTED){
		for (i=0;i<orderedReferences.size();i++){
		    ddim = ((Double) highlightAverageDimension.elementAt(i)).doubleValue();
		    cdim = ((Double) highlightAverageDimension_card.elementAt(i)).doubleValue();

		    if (cdim > 0.0){
			ddim = 1.0 / ddim;
			ddim *= cdim;
		    }else{
			ddim = -1.0;
		    }
		    
		    highlightAverageDimension.setElementAt(new Double(ddim),i);
		    
		    System.out.println(orderedReferences.elementAt(i) + " --- Expect Dim = " + highlightAverageDimension.elementAt(i));
		}
	    }

	    AGCT.Referenced_Available = true;
	    AGCT.Max_Reference = orderedReferences.size()-1;
	    myAGCT.onlyReferenced.setEnabled(true);
	    myAGCT.onlyReferencedEdges.setEnabled(true);
	}
    }

    public void loadAnnotations(String path){
	FileReader e;
	StringTokenizer t;
	String s, copys;
	BufferedReader br;
	String specialComments = new String ("//"), emptyString = "";
	String fullName = path + geneAnnotationsName;
	boolean collectData = false, foundFormat = false;
	Vector variableTypes = null;
	int typeFormat = -1;

	myAGCT.setInfoText("found annotations file : " + geneAnnotationsName + " ... loading " + fullName);
	
	if (AGCT.Debug)
	    System.out.println("found annotations file : " + geneAnnotationsName + " ... loading " + fullName);

	try{
	    e = new FileReader(new File(fullName));
	}catch(FileNotFoundException ex){
	    myAGCT.myInformationFrame.setText("The annotation file you try to load does not exist");
	    myAGCT.annotationsExist = myAGCT.annotationsExistP = myAGCT.annotationsExistF = myAGCT.annotationsExistC = false;
	    return;
	}
	br = new BufferedReader(e);
	try{
	    while ( (s=br.readLine()) != null){
		if (!collectData){
		    t = new StringTokenizer(s, "\t");

		    if (t.countTokens()>0) {
			s = t.nextToken();

			if (s.equals(AGCTFileWriter.ANNO_Data))
			    collectData = true;
			else{
			    if (s.equals(AGCTFileWriter.ANNO_Format)){
				s = t.nextToken();
				if (s.equals(AGCTFileWriter.ANNO_Format_AGCT)){
				    typeFormat = 0;
				    foundFormat = true;
				}else if (s.equals(AGCTFileWriter.ANNO_Format_GEO)){
				    typeFormat = 1;
				    foundFormat = true;
				}
			    }else if (s.equals(AGCTFileWriter.ANNO_Empty)){
				emptyString = t.nextToken();
			    }else if (s.equals(AGCTFileWriter.ANNO_Variable_Definition)){
				variableTypes = new Vector();
				while(t.hasMoreTokens()){
				    s = t.nextToken();
				    if ( (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Dum))
					 && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_ID))
					 && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Ascii_Names))
					 && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_F))
					 && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_P))
					 && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Ontology_C))
					 && (!s.equals(AGCTFileWriter.ANNO_Variable_Definition_Referenced)) )
					Matrix.perror("Domain.class :: unknown variable tag while scanning file " + path);
				    variableTypes.addElement(s);
				}
			    }
			}
		    }
		}else{
		    if (!foundFormat)
			Matrix.perror("Domain.class :: no file format found inside Annotation file " + path);
			
		    if (typeFormat == 0)
			loadGeneAnnotationAGCT_Format(s);
		    else if (typeFormat == 1){
			if (variableTypes == null)
			    Matrix.perror("Domain.class :: no variable format found inside Annotation file " + path);
			if (emptyString.equals(new String("")))
			    Matrix.perror("Domain.class :: no 'Not_Defined' String defined inside Annotation file " + path);
			loadGeneAnnotationGEO_Format(s, variableTypes, emptyString);
		    }
		    else
			Matrix.perror("Domain.class :: bad format value inside Annotation file " + path);
		}
	    }
	    e.close();
	}catch(IOException ex){
	    myAGCT.myInformationFrame.setText("IOError when loading annotation file");
	    myAGCT.annotationsExist = myAGCT.annotationsExistP = myAGCT.annotationsExistF = myAGCT.annotationsExistC = false;
	    return;
	}
	e = null;

	myAGCT.annotationsExist = true;
	myAGCT.setInfoText("found annotations file : " + geneAnnotationsName + " ... loading... ok.");
   }

    /*****
    public void loadAnnotations(String path){
	FileReader e;
	StringTokenizer t;
	String s;
	BufferedReader br;
	Gene gg;
	boolean readingF = false, readingP = false, readingC = false;
	String specialComments = new String ("//");
	String fullName = path + geneAnnotationsName;

	myAGCT.setInfoText("found annotations file : " + geneAnnotationsName + " ... loading " + fullName);
	
	if (AGCT.Debug)
	    System.out.println("found annotations file : " + geneAnnotationsName + " ... loading " + fullName);

	try{
	    e = new FileReader(new File(fullName));
	}catch(FileNotFoundException ex){
	    myAGCT.myInformationFrame.setText("The annotation file you try to load does not exist");
	    myAGCT.annotationsExist = myAGCT.annotationsExistP = myAGCT.annotationsExistF = myAGCT.annotationsExistC = false;
	    return;
	}
	br = new BufferedReader(e);
	try{
	    while ( (s=br.readLine()) != null){
		t = new StringTokenizer(s.replace('\t',' '), " ");
		if (t.countTokens()>0) {
		    s = t.nextToken();
		    if (!s.substring(0,specialComments.length()).equals(specialComments)){
			if (accessToGenesIndexBefore.containsKey(s)){
			    readingF = readingP = readingC = false;
			    gg = (Gene) domainGenes.elementAt(((Integer) accessToGenesIndexBefore.get(s)).intValue());
			    while (t.hasMoreTokens()){
				s = t.nextToken();
				if (s.equals(Domain.ANNOTATION_F))
				    readingF = true;
				else if (s.equals(Domain.ANNOTATION_P))
				    readingP = true;
				else if (s.equals(Domain.ANNOTATION_C))
				    readingC = true;
				else{
				    if (s != ""){
					if (readingF || readingP || readingC){
					    if ( (readingF && readingP)
						 || (readingF && readingC)
						 || (readingP && readingC) )
						System.out.print(" Gene " + gg.name);

					    gg.addAnnotation(readingF, readingP, readingC, s);
					    if (readingF)
						myAGCT.annotationsExistF = true;
					    else if (readingP)
						myAGCT.annotationsExistP = true;
					    else if (readingC)
						myAGCT.annotationsExistC = true;

					    readingF = readingP = readingC = false;
					}else{
					    if (gg.asciiName == null)
						gg.asciiName = s;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	    e.close();
	}catch(IOException ex){
	    myAGCT.myInformationFrame.setText("IOError when loading annotation file");
	    myAGCT.annotationsExist = myAGCT.annotationsExistP = myAGCT.annotationsExistF = myAGCT.annotationsExistC = false;
	    return;
	}
	e = null;

	myAGCT.annotationsExist = true;
	myAGCT.setInfoText("found annotations file : " + geneAnnotationsName + " ... loading... ok.");
	}*****/


    public void fillSelectedGenesAllAnnotations(){
	Thread t = new Thread(){
		public void run(){
		    
		    selectedGenesAllAnnotationsF = new Vector();
		    selectedGenesAllAnnotationsP = new Vector();
		    selectedGenesAllAnnotationsC = new Vector();

		    annotationFToGenes = annotationPToGenes = annotationCToGenes = null;

		    int i,j;
		    Gene gg;
		    String s, t, sano;
		    Vector dumvect;

		    AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Filling Annotations [1/3]", numberSelectedGenesBefore);
		    for (i=0;i<numberSelectedGenesBefore;i++){
			gg = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumberBefore[i]);
			if (gg.annotationsF != null){
			    if (annotationFToGenes == null)
				annotationFToGenes = new Hashtable();
			    for (j=0;j<gg.annotationsF.size();j++){
				s = ((String[]) gg.annotationsF.elementAt(j))[0];
				sano = ((String[]) gg.annotationsF.elementAt(j))[1];
				if (!s.equals(Domain.ANNOTATION_F))
				    Matrix.perror("Domain.class :: annotation != F in F annotations");
				if (!annotationFToGenes.containsKey(sano)){
				    selectedGenesAllAnnotationsF.add(new String(sano));
				    annotationFToGenes.put(new String(sano), new Vector());
				}
				dumvect = (Vector) annotationFToGenes.get(sano);
				dumvect.add(new Integer(selectedGeneNumberToGeneNumberBefore[i]));
			    }
			}
			cc.increment();
		    }
		    cc.end();

		    cc = new AGCTCounter(myAGCT.myInformationFrame, "Filling Annotations [2/3]", numberSelectedGenesBefore);
		    for (i=0;i<numberSelectedGenesBefore;i++){
			gg = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumberBefore[i]);
			if (gg.annotationsP != null){
			    if (annotationPToGenes == null)
				annotationPToGenes = new Hashtable();
			    for (j=0;j<gg.annotationsP.size();j++){
				s = ((String[]) gg.annotationsP.elementAt(j))[0];
				sano = ((String[]) gg.annotationsP.elementAt(j))[1];
				if (!s.equals(Domain.ANNOTATION_P))
				    Matrix.perror("Domain.class :: annotation != P in P annotations");
				if (!annotationPToGenes.containsKey(sano)){
				    selectedGenesAllAnnotationsP.add(new String(sano));
				    annotationPToGenes.put(new String(sano), new Vector());
				}
				dumvect = (Vector) annotationPToGenes.get(sano);
				dumvect.add(new Integer(selectedGeneNumberToGeneNumberBefore[i]));
			    }
			}
			cc.increment();
		    }
		    cc.end();

		    cc = new AGCTCounter(myAGCT.myInformationFrame, "Filling Annotations [3/3]", numberSelectedGenesBefore);
		    for (i=0;i<numberSelectedGenesBefore;i++){
			gg = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumberBefore[i]);
			if (gg.annotationsC != null){
			    if (annotationCToGenes == null)
				annotationCToGenes = new Hashtable();
			    for (j=0;j<gg.annotationsC.size();j++){
				s = ((String[]) gg.annotationsC.elementAt(j))[0];
				sano = ((String[]) gg.annotationsC.elementAt(j))[1];
				if (!s.equals(Domain.ANNOTATION_C))
				    Matrix.perror("Domain.class :: annotation != C in C annotations");
				if (!annotationCToGenes.containsKey(sano)){
				    selectedGenesAllAnnotationsC.add(new String(sano));
				    annotationCToGenes.put(new String(sano), new Vector());
				}
				dumvect = (Vector) annotationCToGenes.get(sano);
				dumvect.add(new Integer(selectedGeneNumberToGeneNumberBefore[i]));
			    }
			}
			cc.increment();
		    }
		    cc.end();
		    
		    /*cc = new AGCTCounter(myAGCT.myInformationFrame, "Filling Annotations [2/4]", selectedGenesAllAnnotationsF.size());
		    for (i=0;i<selectedGenesAllAnnotationsF.size();i++){
			s = new String((String) selectedGenesAllAnnotationsF.elementAt(i));
			while(selectedGenesAllAnnotationsF.contains(s))
			    selectedGenesAllAnnotationsF.removeElement(s);
			selectedGenesAllAnnotationsF.insertElementAt(s,i);
			cc.increment();
		    }
		    cc.end();*/
		    
		    if (selectedGenesAllAnnotationsF.size() > 0){
			QuickSort.quicksortString(selectedGenesAllAnnotationsF);
			someAnnotationsF = true;
		    }else
			selectedGenesAllAnnotationsF = null;
		    
		    /*cc = new AGCTCounter(myAGCT.myInformationFrame, "Filling Annotations [3/4]", selectedGenesAllAnnotationsP.size());
		    for (i=0;i<selectedGenesAllAnnotationsP.size();i++){
			s = new String((String) selectedGenesAllAnnotationsP.elementAt(i));
			while(selectedGenesAllAnnotationsP.contains(s))
			    selectedGenesAllAnnotationsP.removeElement(s);
			selectedGenesAllAnnotationsP.insertElementAt(s,i);
			cc.increment();
		    }
		    cc.end();*/

		    if (selectedGenesAllAnnotationsP.size() > 0){
			QuickSort.quicksortString(selectedGenesAllAnnotationsP);
			someAnnotationsP = true;
		    }else
			selectedGenesAllAnnotationsP = null;

		    /*cc = new AGCTCounter(myAGCT.myInformationFrame, "Filling Annotations [4/4]", selectedGenesAllAnnotationsC.size());
		    for (i=0;i<selectedGenesAllAnnotationsC.size();i++){
			s = new String((String) selectedGenesAllAnnotationsC.elementAt(i));
			while(selectedGenesAllAnnotationsC.contains(s))
			    selectedGenesAllAnnotationsC.removeElement(s);
			selectedGenesAllAnnotationsC.insertElementAt(s,i);
			cc.increment();
		    }
		    cc.end();*/
		    
		    if (selectedGenesAllAnnotationsC.size() > 0){
			QuickSort.quicksortString(selectedGenesAllAnnotationsC);
			someAnnotationsC = true;
		    }else
			selectedGenesAllAnnotationsC = null;
		    
		    // End: filling Hashtables

		    if (someAnnotationsF){
			accessToSelectedGenesAllAnnotationsF = new Hashtable();
			for (i=0;i<selectedGenesAllAnnotationsF.size();i++){
			    accessToSelectedGenesAllAnnotationsF.put((String) selectedGenesAllAnnotationsF.elementAt(i), new Integer(i));
			}
		    }
		    
		    if (someAnnotationsP){
			accessToSelectedGenesAllAnnotationsP = new Hashtable();
			for (i=0;i<selectedGenesAllAnnotationsP.size();i++){
			    accessToSelectedGenesAllAnnotationsP.put((String) selectedGenesAllAnnotationsP.elementAt(i), new Integer(i));
			}
		    }

		    if (someAnnotationsC){
			accessToSelectedGenesAllAnnotationsC = new Hashtable();
			for (i=0;i<selectedGenesAllAnnotationsC.size();i++){
			    accessToSelectedGenesAllAnnotationsC.put((String) selectedGenesAllAnnotationsC.elementAt(i), new Integer(i));
			}
		    }

		    ControlProcess.put("selectedGenesAllAnnotationsOK",true);
		}
	    };
	t.start();
    }

    public boolean geneExpressionContainsEmptyData(Gene g, int n){
	//returns true iff at least b values NON (?) for EACH selected LIGAND
	//if (ControlProcess.assertTrue("selectedLigands"))
	//Matrix.perror ("You need to validate ligands before processing genes");
	
	int i=0, j, nb, s;
	for (i=0;i<numberLigands;i++)
	    if ( (selectedLigand[i]) && (enabledLigand[i]) ){
		nb = 0;
		for (j=0;j<( (Integer) numberTimes.elementAt(i)).intValue();j++){
		    if (g.undeterminedCoordinates[i][j] == false)
			nb++;
		}
		if (nb < n)
		    return false;
	    }
	return true;
    }

    public boolean geneExpressionIsEligible(Gene g){
	// returns true iff
	if (AGCT.Method_F == 0){
	    //slopes
	    return geneExpressionContainsEmptyData(g, 2);
	}else if ( (AGCT.Method_F >= 1) && (AGCT.Method_F <= 10 ) ){
	    //Haar wavelets MINUS the constant coefficient
	    //(does not represent a variation of the curve)
	    //OR Daubechies interval DX(i), with all coefficients

	    return geneExpressionContainsEmptyData(g, 2);
	}
	Matrix.perror("Feature method invalid");
	return true;
    }

    public int getLigandIndex(String s){
	return ((Integer) accessToLigands.get(s)).intValue();
    }
    public String getLigandName(int i){
	return (String) ((Vector) domainLigands.elementAt(i)).elementAt(0);
    }
    public String getLigandGroupName(int l, int g){
	return (String) ((Vector) domainLigands.elementAt(l)).elementAt(g + 1);
    }
    public String getLigandGroupName(String l, int g){
	return (String) ((Vector) domainLigands.elementAt(getLigandIndex(l))).elementAt(g + 1);
    }
    public int getFeatureIndex(String s){
	return ((Integer) accessToFeatures.get(s)).intValue();
    }

    public int getGeneIndex(String s){
	int i;
	for (i=0;i<domainGenes.size();i++)
	    if ( s.equals( ( (Gene) domainGenes.elementAt(i) ).name ) )
		return i;
	Matrix.perror("Domain.class :: no such gene (" + s + ")");
	return -1;
    }

    public void initializationFilters(){
	int i, j;

	selectedGeneBefore = new boolean [numberInitialGenes];
	enabledGeneBefore = new boolean [numberInitialGenes];

	for (i=0;i<numberInitialGenes;i++){
	    selectedGeneBefore[i] = true;
	    enabledGeneBefore[i] = true;
	}

	selectedLigand = new boolean [numberLigands];
	enabledLigand = new boolean [numberLigands];

	for (i=0;i<numberLigands;i++){
	    selectedLigand[i] = true;
	    enabledLigand[i] = true;
	}

	selectedGroup = new boolean [numberGroups][];
	for (i=0;i<numberGroups;i++){
	    selectedGroup[i] = new boolean [domainGroupElements[i].size()];
	    for (j=0;j<domainGroupElements[i].size();j++)
		selectedGroup[i][j] = true;
	}

    }

    public void updateEnabledGene(){
	int i;
	for (i=0;i<numberInitialGenes;i++)
	    enabledGeneBefore[i] = geneExpressionIsEligible( (Gene) domainGenes.elementAt(i) );
    }


    public void updateEnabledLigand(){
	//updates enabledLigand
	int i, j, k;
	boolean found = false;
	boolean keep = true;
	String ls, gs;

	for (i=0;i<domainLigands.size();i++){
	    j=0;
	    keep = true;
	    while( (j<numberGroups) && (keep == true) ){
		ls = (String) ((Vector) domainLigands.elementAt(i)).elementAt(j+1);
		k = 0;
		found = false;
		do{
		    gs = (String) domainGroupElements[j].elementAt(k);
		    if (gs.equals(ls)){
			found = true;
		    }else
			k++;
		}while( (found == false) && (k < domainGroupElements[j].size()) );
		if (found == false)
		    Matrix.perror("Could not find " + AGCTFileWriter.DATA_Ligand + " name " + ls + " in groups\n");
		if (selectedGroup[j][k] == false)
		    keep = false;
		else
		    j++;
	    }
	    enabledLigand[i] = keep;
	}
    }

    public void lastUpdateFilterData(){
	//Vars updated after prototype selection
	//numberSelectedGenes = -1;
	//selectedGeneNumberToGeneNumber = null;
	//accessToSelectedGeneNumberToGeneNumber = null;
	//accessToGenesObject = null;
	////////////////////////////////////////////////

	int i, index;
	Integer I;

	if (!Prototype.Prototypes_Selected)
	    Matrix.perror("Domain.class :: prototypes not selected");

	numberSelectedGenes = 0;
	accessToSelectedGeneNumberToGeneNumber = new Hashtable();
	accessToGenesObject = new Hashtable();
	
	if (Prototype.No_Reduction){
	    for (i=0;i<numberInitialGenes;i++)
		if (selectedGeneBefore[i])
		    numberSelectedGenes++;
	    selectedGeneNumberToGeneNumber = new int [numberSelectedGenes];
	    
	    index = 0;
	    for (i=0;i<numberInitialGenes;i++)
		if (selectedGeneBefore[i]){
		    selectedGeneNumberToGeneNumber[index] = i;
		    accessToGenesObject.put(new Integer(i), (Gene) domainGenes.elementAt(i));
		    accessToSelectedGeneNumberToGeneNumber.put(new Integer(index), new Integer(i));
		    
		    index ++;
		}
	}else{
	    numberSelectedGenes = Prototype.Closest_Center_Ordered_List.size();
	    selectedGeneNumberToGeneNumber = new int [numberSelectedGenes];
	    
	    for (index=0;index<Prototype.Closest_Center_Ordered_List.size();index++){
		I = (Integer) Prototype.Closest_Center_Ordered_List.elementAt(index);
		i = I.intValue();
		selectedGeneNumberToGeneNumber[index] = i;
		accessToGenesObject.put(new Integer(i), (Gene) domainGenes.elementAt(i));
		accessToSelectedGeneNumberToGeneNumber.put(new Integer(index), new Integer(i));
	    }

	    /*Enumeration extensions = Prototype.Closest_Center_To_Cluster_Number.keys();
	    index = 0;
	    if(extensions != null) {
		while (extensions.hasMoreElements()) {
		    I = (Integer) extensions.nextElement();
		    i = I.intValue();
		    selectedGeneNumberToGeneNumber[index] = i;
		    accessToGenesObject.put(new Integer(i), (Gene) domainGenes.elementAt(i));
		    accessToSelectedGeneNumberToGeneNumber.put(new Integer(index), new Integer(i));

		    index ++;
		}
		}*/
	}

	ControlProcess.put("genesFilteredFor" + AGCT.Feature_Method[AGCT.Method_F],true);
    }

    public void updateFilterData(){
	String valret = "";
	int i, index, ns;

	//Number of Ligands : integrate enabled and selected ligands

	for (i=0;i<numberLigands;i++){
	    if ( (selectedLigand[i]) && (enabledLigand[i]) )
		selectedLigand[i] = true;
	    else
		selectedLigand[i] = false;
	}

	numberSelectedLigands = 0;
	for (i=0;i<numberLigands;i++){
	    if (selectedLigand[i] == true)
		numberSelectedLigands ++;
	}

	if (numberSelectedLigands == 0) 
	    Matrix.perror("Too many constraints on " +  AGCTFileWriter.DATA_Ligand + " : 0 " +  AGCTFileWriter.DATA_Ligand + " selected"); 

	//Number of Genes : integrates enabled and selected genes

	for (i=0;i<numberInitialGenes;i++)
	    if ( (selectedGeneBefore[i]) && (enabledGeneBefore[i]) )
		selectedGeneBefore[i] = true;
	    else
		selectedGeneBefore[i] = false;

	//Vars to update later after prototype selection

	numberSelectedGenes = -1;
	selectedGeneNumberToGeneNumber = null;
	accessToSelectedGeneNumberToGeneNumber = null;
	accessToGenesObject = null;

	////////////////////////////////////////////////

	numberSelectedGenesBefore = 0;
	for (i=0;i<numberInitialGenes;i++)
	    if (selectedGeneBefore[i])
		numberSelectedGenesBefore++;
	selectedGeneNumberToGeneNumberBefore = new int [numberSelectedGenesBefore];
	accessToSelectedGeneNumberToGeneNumberBefore = new Hashtable();
	accessToGenesObjectBefore = new Hashtable();

	index = 0;
	for (i=0;i<numberInitialGenes;i++)
	    if (selectedGeneBefore[i]){
		selectedGeneNumberToGeneNumberBefore[index] = i;
		accessToGenesObjectBefore.put(new Integer(i), (Gene) domainGenes.elementAt(i));
		accessToSelectedGeneNumberToGeneNumberBefore.put(new Integer(index), new Integer(i));

		index ++;
	    }

	//ControlProcess.put("genesFilteredFor" + AGCT.Feature_Method[AGCT.Method_F],true);

	for (i=0;i<numberSelectedGenesBefore;i++){
	    ( (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumberBefore[i])).finishVariables(myAGCT, this);
	    ( (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumberBefore[i])).recomputeSelectedFeatures();
	}

	valret = "#Ligands selected: " + numberSelectedLigands + "\n";
	ns = 0;
	for (i=0;i<numberLigands;i++){
	    if (selectedLigand[i]){
		if (ns > 0)
		    valret += ", ";
		valret += (String) ((Vector) domainLigands.elementAt(i)).elementAt(0);
		ns ++;
	    }
	}
	valret += "\n";
	valret += "#Genes selected before eventual prototype selection: " + numberSelectedGenesBefore + "\n";
	for (i=0;i<numberSelectedGenesBefore;i++){
	    if (i<AGCT.Max_Number_Of_Genes_Kept){
		valret += ( (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumberBefore[i])).name;
		if (i < numberSelectedGenesBefore - 1)
		    valret += ", ";
	    }
	}
	valret += "\n";

	myAGCT.setInfoText(valret);
	ControlProcess.put("filtersProcessed",true);
    }
    
    public Gene getSelectedGene(int index){
	Integer I = new Integer(index);
	Gene g = (Gene) accessToGenesObject.get((Integer) accessToSelectedGeneNumberToGeneNumber.get(I));
	I = null;
	return g;
    }


    /****************************************************************************************
     * Feature stuff
     *****/

    public void computeGeneFeaturesWavelets(){
	int i, ig, index, j, k, nstamps, l, m;

	ControlProcess.assertTrue("filtersProcessed");

	//Beware: stamps is used to compute haar basis, to save space

	double [][] stamps = new double [numberSelectedLigands][];
	double [] delta = new double [numberSelectedLigands];
	double [] begintime = new double [numberSelectedLigands];
	double [] cvect;

	index = 0;
	for (i=0;i<numberLigands;i++){
	    if (selectedLigand[i]){
		if (AGCT.Number_Of_Wavelet_Stamps == -1)
		    stamps[index] = Feature.waveletsStampsAutomatic(myAGCT, domainTimes, i);
		else
		    stamps[index] = Feature.waveletsStampsUserFixed(myAGCT, domainTimes, i);

		delta[index] = stamps[index][1] - stamps[index][0];
		begintime[index] = stamps[index][0];
		index++;
	    }
	}
	
	AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Feature construction with Wavelets", numberSelectedGenesBefore);
	for (i=0;i<numberSelectedGenesBefore;i++){
	    ig = selectedGeneNumberToGeneNumberBefore[i];
	    index = 0;
	    for (j=0;j<numberLigands;j++){
		if (selectedLigand[j]){
		    nstamps = stamps[index].length;
		    cvect = new double[nstamps];
		    Feature.fillTimeSerie(domainGenes, domainTimes, cvect, begintime, delta, nstamps, ig, index, j);

		    if (AGCT.Method_F == 1)
			((Gene) domainGenes.elementAt(ig)).finalCoordinates[index] = Feature.getHaar((Gene) domainGenes.elementAt(ig), cvect);
		    else if ( (AGCT.Method_F >= 2) && (AGCT.Method_F <= 10) )
			((Gene) domainGenes.elementAt(ig)).finalCoordinates[index] = Feature.getDn(AGCT.Method_F - 2,cvect);
		    else
			Matrix.perror("Unknown wavelet method");

		    index++;
		}
	    }
	    cc.increment();
	}
	cc.end();
    }
    
    public void computeGeneFeaturesSlopes(){
	int i, ig, index, j;
	double [] tab;

	ControlProcess.assertTrue("filtersProcessed");

	AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Feature construction with Slopes", numberSelectedGenesBefore);
	for (i=0;i<numberSelectedGenesBefore;i++){
	    ig = selectedGeneNumberToGeneNumberBefore[i];
	    index = 0;
	    for (j=0;j<numberLigands;j++){
		if (selectedLigand[j]){
		    tab = new double[1];
		    tab[0] = Feature.slope(domainTimes, j,((Gene) domainGenes.elementAt(ig)).rawCoordinates[index], ((Gene) domainGenes.elementAt(ig)).undeterminedCoordinates[index]);

		    ((Gene) domainGenes.elementAt(ig)).finalCoordinates[index] = tab;

		    index++;
		}
	    }
	    cc.increment();
	}
	cc.end();
    }

    public void finishComputeGeneFeatures(){
	//while (!Feature.featuresSelected);

	myAGCT.dimFeatures = Feature.getNumberOfFeaturesAfterSelection();

	myAGCT.myTabbedPane.mySelectionPane.finishValidate();
	ControlProcess.put("featuresSelected",true);
    }

    public void toTabbedPane(){

	myAGCT.myTabbedPane.toManifold();
	myAGCT.myTabbedPane.toPCA();

    }

    public void prototypeSelection(){
	Prototype.prototypeSelection(this);
	ControlProcess.put("prototypeSelected",true);
    }

    public void computeFinalCoordinates_1D(){
	int i;
	myAGCT.timerTickOn();	
	AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Computing Gene Coordinates", numberSelectedGenesBefore);
	for (i=0;i<numberSelectedGenesBefore;i++){
	    ((Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumberBefore[i])).computeFinalCoordinates_1D();
	    cc.increment();
	}
	cc.end();
	myAGCT.timerTickOff();
    }

    public void featureConstruction(){
	myAGCT.timerTickOn();

	if (AGCT.Method_F == 0)
	    computeGeneFeaturesSlopes();
	else if ( (AGCT.Method_F >= 1) && (AGCT.Method_F <= 10) )
	    computeGeneFeaturesWavelets();
	else 
	    Matrix.perror("Method selected to compute features does not exist !");

	Gene.computeFinalIndex(((Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumberBefore[0])).finalCoordinates);
	myAGCT.timerTickOff();
    }

    public void featureSelection(){
	if ( (AGCT.Method_FS > 0) && (AGCT.Max_Number_Of_Features < Feature.getNumberOfFeaturesBeforeSelection(myAGCT)) ){
	    if (AGCT.Max_Number_Of_Features < 0){
		if (AGCT.Debug)
		    System.out.println("Fixing max number of features to " + AGCT.DEFAULT_Max_Number_Of_Features);
		AGCT.Max_Number_Of_Features = AGCT.DEFAULT_Max_Number_Of_Features;
	    }

	    Feature.featureSelection(this);
	}else
	    Feature.noFeatureSelection(this);
    }

    public static double [] MAX_LOCAL_DIMENSION;

    public static void FILL_MAX_LOCAL_DIMENSION(AGCT a){
	MAX_LOCAL_DIMENSION = new double [AGCT.Max_K_Dimension_Estimation];
	int i, j;
	double maxd = 0.0;
	Gene g;
	for (i=0;i<AGCT.Max_K_Dimension_Estimation;i++){
	    for (j=0;j<a.dimElements;j++){
		g = (Gene) a.myDomain.domainGenes.elementAt(a.myDomain.selectedGeneNumberToGeneNumber[j]);
		if ( (j == 0) || (g.local_dimension[i] > maxd) )
		    maxd = g.local_dimension[i];
	    }
	    MAX_LOCAL_DIMENSION[i] = maxd;
	}
    }

    public void toNeighborsDimension(){
	AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "computing " + AGCT.Max_K_Dimension_Estimation + " NNs for dimension estimation", myAGCT.dimElements * (myAGCT.dimElements - 1));
	
	int i, j, ind;
	int [] indexes = new int [myAGCT.dimElements-1];
	double [] distances = new double [myAGCT.dimElements-1];
	double distance;
	Gene gi, gj;

	for (i=0;i<myAGCT.dimElements;i++){
	    gi = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumber[i]);
	    ind = 0;
	    for (j=0;j<myAGCT.dimElements;j++){
		if (j!=i){
		    gj = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumber[j]);
		    distance = Math.sqrt(Gene.squaredMagnitudeDifference(gi, gj));
		    indexes[ind] = j;
		    distances[ind] = distance;
		    ind++;
		    
		    cc.increment();
		}
	    }
	    QuickSort.quicksort(distances,indexes);
	    
	    gi.fillNeighborsDimension(distances,indexes);
	}
	for (i=0;i<myAGCT.dimElements;i++){
	    gi = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumber[i]);
	    gi.computeAverageDimension();
	}

	Domain.FILL_MAX_LOCAL_DIMENSION(myAGCT);
	myAGCT.kDimensionEstimation.setEnabled(true);
	myAGCT.displayDimension.setEnabled(true);
	AGCT.LOCAL_DIMENSION_COMPUTED = true;
	myAGCT.clustersToDimension();
	myAGCT.repaint();

	Statistics.fillHistogramTagAverageDimension(myAGCT.myTabbedPane.myManifoldPane.statF,myAGCT.myTabbedPane.myManifoldPane.statP,myAGCT.myTabbedPane.myManifoldPane.statC, myAGCT);

	cc.end();
    }

    public void featureAndPrototypeSelection(){
	Thread t = new Thread(){
		public void run(){
		    checkNames();
		    featureConstruction();
		    featureSelection();
		    computeFinalCoordinates_1D();

		    if (Prototype.Loading_From_Scenario_Begin)
			while (!Prototype.Loading_From_Scenario_End);
		    else
			prototypeSelection();

		    lastUpdateFilterData();
		    if ( (someAnnotationsF) || (someAnnotationsP) ){
			myAGCT.myAnnotationFrame = new JAnnotationFrame("AGCT --- Annotation Frame", myAGCT, myAGCT.myDomain);
			myAGCT.myAnnotationFrame.setSize(AGCT.WindowWidth,600);
			myAGCT.myFrame.toFront();
			myAGCT.myAnnotationFrame.changeLF_Prototype_Selection();
			if (!Prototype.No_Reduction)
			    myAGCT.myAnnotationFrame.onlyPrototypes.setEnabled(true);
		    }

		    toTabbedPane();
		    ControlProcess.put("featuresAndPrototypesSelected", true);
		    myAGCT.dimElements = numberSelectedGenes;
		}
	    };
	t.start();
    }

    public void computeMinMax_ManifoldPnt3D(){
	//Makes the assumption that gene components in manifold_point3D HAVE BEEN computed OR correctly updated
	Gene gg;
	int i, j, id;
	for (i=0;i<numberSelectedGenes;i++){
	    id = selectedGeneNumberToGeneNumber[i];
	    gg = (Gene) domainGenes.elementAt(id);
	    
	    for (j=0;j<3;j++){
		if ( (i==0) || (gg.manifold_point3D.coordinates[j] < min_ManifoldPnt3D.coordinates[j]) )
		    min_ManifoldPnt3D.coordinates[j] = gg.manifold_point3D.coordinates[j];
		if ( (i==0) || (gg.manifold_point3D.coordinates[j] > max_ManifoldPnt3D.coordinates[j]) )
		    max_ManifoldPnt3D.coordinates[j] = gg.manifold_point3D.coordinates[j];
	    }
	}
    }

    public void computeMinMax_PCAPnt3D(){
	//Makes the assumption that gene components in pca_point3D HAVE BEEN computed OR correctly updated
	Gene gg;
	int i, j, id;
	for (i=0;i<numberSelectedGenes;i++){
	    id = selectedGeneNumberToGeneNumber[i];
	    gg = (Gene) domainGenes.elementAt(id);
	    
	    for (j=0;j<3;j++){
		if ( (i==0) || (gg.pca_point3D.coordinates[j] < min_PCAPnt3D.coordinates[j]) )
		    min_PCAPnt3D.coordinates[j] = gg.pca_point3D.coordinates[j];
		if ( (i==0) || (gg.pca_point3D.coordinates[j] > max_PCAPnt3D.coordinates[j]) )
		    max_PCAPnt3D.coordinates[j] = gg.pca_point3D.coordinates[j];
	    }
	}
    }

    public void renormalizePCAComponentsToUnit(int x, int y, int z){
	AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Normalizing pca comp.", 2 * numberSelectedGenes);
	int i, id;
	double norm, maxnorm = 0.0;
	Gene gg;

	for (i=0;i<numberSelectedGenes;i++){
	    id = selectedGeneNumberToGeneNumber[i];
	    gg = (Gene) domainGenes.elementAt(id);

	    norm = (gg.pca_components.coordinates[x]*gg.pca_components.coordinates[x]);
	    norm += (gg.pca_components.coordinates[y]*gg.pca_components.coordinates[y]);
	    norm += (gg.pca_components.coordinates[z]*gg.pca_components.coordinates[z]);
	    norm = Math.sqrt(norm);

	    if ( (i==0) || (norm > maxnorm) )
		maxnorm = norm;

	    cc.increment();
	}
	if (maxnorm == 0.0)
	    Matrix.perror("No PCA renormalization possible: all PCA coordinates = 0 for (x,y,z) = (" + x + "," + y + "," + z + ") for all genes");

	for (i=0;i<numberSelectedGenes;i++){
	    id = selectedGeneNumberToGeneNumber[i];
	    gg = (Gene) domainGenes.elementAt(id);
	    
	    gg.pca_Pnt3D.coordinates[0] = gg.pca_components.coordinates[x] / maxnorm;
	    gg.pca_Pnt3D.coordinates[1] = gg.pca_components.coordinates[y] / maxnorm;
	    gg.pca_Pnt3D.coordinates[2] = gg.pca_components.coordinates[z] / maxnorm;

	    cc.increment();
	}
	cc.end();
    }

    public void fillPCAComponents(Matrix E, double [] ave, double [] sig){
	AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Filling pca comp.", numberSelectedGenes);
	int i, j, id;
	Gene gg;

	double norm, maxnorm;

	for (i=0;i<numberSelectedGenes;i++){
	    id = selectedGeneNumberToGeneNumber[i];
	    gg = (Gene) domainGenes.elementAt(id);

	    gg.fillPCAComponents(E, ave, sig);
	    cc.increment();
	}
	cc.end();

	if (AGCT.NORMALIZE_UNIT_PCA)
	    renormalizePCAComponentsToUnit(0, 1, 2);

	computeMinMax_PCAPnt3D();

    }

    public void fillManifoldComponents(Matrix M){
	AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Filling manifold components", numberSelectedGenes);

	int i, j, id;

	double dum;

	Gene gg;

	for (i=0;i<numberSelectedGenes;i++){
	    cc.increment();

	    id = selectedGeneNumberToGeneNumber[i];
	    gg = (Gene) domainGenes.elementAt(id);

	    gg.fillManifoldComponents(M, i);
	}

	computeMinMax_ManifoldPnt3D();
	cc.end();
    }

    public void addIfNotInside(int [] couple){
	int i;
	Integer [] c;
	boolean trouve = false;
	Gene gg;

	if (myEdges == null){
	    myEdges = new Vector();
	}else{
	    i = 0;
	    do{
		c = (Integer []) myEdges.elementAt(i);
		if ( ( (c[0].intValue() == couple[0]) && (c[1].intValue() == couple[1]) ) 
		     || ( (c[0].intValue() == couple[1]) && (c[1].intValue() == couple[0]) ) )
		    trouve = true;
		else
		    i++;
	    }while( (trouve == false) && (i<myEdges.size()) );
	}

	if (trouve == false){
	    c = new Integer[2];
	    c[0] = new Integer(couple[0]);
	    c[1] = new Integer(couple[1]);

	    myEdges.addElement(c);

	    gg = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumber[c[0].intValue()]);
	    gg.putNeighbor(c[1]);

	    gg = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumber[c[1].intValue()]);
	    gg.putNeighbor(c[0]);
	}
    }

    public boolean addIfNotInsideFilteringDistances(int [] couple){
	int i;
	Integer [] c;
	boolean trouve = false, ret = false;
	double dist;
	Gene gg;

	if (myEdges == null){
	    myEdges = new Vector();
	}else if (myEdges.size() == 0){
	    trouve = false;
	}else{
	    i = 0;
	    do{
		c = (Integer []) myEdges.elementAt(i);
		if ( ( (c[0].intValue() == couple[0]) && (c[1].intValue() == couple[1]) ) 
		     || ( (c[0].intValue() == couple[1]) && (c[1].intValue() == couple[0]) ) )
		    trouve = true;
		else
		    i++;
	    }while( (trouve == false) && (i<myEdges.size()) );
	}

	if (trouve == false){
	    //geneToCoord[c[0].intValue()]

	    c = new Integer[2];
	    c[0] = new Integer(couple[0]);
	    c[1] = new Integer(couple[1]);


	    dist = myAGCT.W.distance( ((Vector) ((Vector) myAGCT.allGenesCoordinates.elementAt(myAGCT.geneToCoord[c[0].intValue()])).elementAt(1)), ((Vector) ((Vector) myAGCT.allGenesCoordinates.elementAt(myAGCT.geneToCoord[c[1].intValue()])).elementAt(1)));

	    if (dist <= AGCT.MAX_DISTANCE_TRIANGULATION){
		myEdges.addElement(c);
		
		gg = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumber[c[0].intValue()]);
		gg.putNeighbor(c[1]);
		
		gg = (Gene) domainGenes.elementAt(selectedGeneNumberToGeneNumber[c[1].intValue()]);
		gg.putNeighbor(c[0]);

		ret = true;
	    }
	}

	return ret;
    }

    public void toEdges(){
	if (myDelaunayTriangulation != null){
	    AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Creating edges", myDelaunayTriangulation.getNeighbors().size());

	    myAGCT.myInformationFrame.setTextProgressBar("Creating edges");
	    Simplex s;
	    Iterator it;
	    int i, ip1, ip2, iu1, iu2, nmade = 0, percent, total = 0, accept = 0;
	    boolean added;
	    int [] couple = new int [2];
	    Vector v;
	    for (it = myDelaunayTriangulation.getNeighbors().keySet().iterator(); it.hasNext();){
		cc.increment();

		s = (Simplex) it.next();
		v = new Vector(s.getVertices());
		for (i=0;i<v.size();i++){
		    total ++;

		    ip1 = i;
		    ip2 = i + 1;
		    if (i == v.size()-1)
			ip2 = 0;
		    
		    iu1 = ((Pnt) v.elementAt(ip1)).getGeneIndex();
		    iu2 = ((Pnt) v.elementAt(ip2)).getGeneIndex();

		    if ( (iu1 >= 0) && (iu2 >= 0) ){
			couple[0] = iu1;
			couple[1] = iu2;

			if (AGCT.Filter_Triangulation_With_Distances){
			    added = addIfNotInsideFilteringDistances(couple);
			    if (added)
				accept++;
			}else
			    addIfNotInside(couple);
		    }
		    if (myEdges != null){
			if (AGCT.Filter_Triangulation_With_Distances)
			    myAGCT.myInformationFrame.setTextProgressBar("Creating edges :: " + myEdges.size() + " current edges --- (" + DF.format(( ((double) accept) / ((double) total)) * 100.0) + " % kept with Distances)");
			else
			    myAGCT.myInformationFrame.setTextProgressBar("Creating edges :: " + myEdges.size() + " current edges.");
		    }
		}
	    }
	    cc.end();
	    myAGCT.myTabbedPane.myManifoldPane.visualizationPanel.repaint();
	}
    }
}
