import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import java.util.*;
import java.text.DecimalFormat;
import javax.swing.Box.Filler;
import java.awt.Rectangle;

class JAGCTCheckBox extends JCheckBox implements ItemListener {
    JAGCTSelectionPane myJAGCTSelectionPane;
    Domain myDomain;
    String publicName;
    String myName;
    int myI, myJ;
    boolean ready;

    static String refGene = "JAGCTCheckBox_gene";
    static String refGroup = "JAGCTCheckBox_group";
    static String refLigand = "JAGCTCheckBox_ligand";

    JAGCTCheckBox(JAGCTSelectionPane ja, Domain dom, String pn, String mn, int i, int j){
	super(pn);
	myJAGCTSelectionPane = ja;
	myDomain = dom;
	myName = mn;
	myI = i;
	myJ = j;
	if (mn.equals(refGroup))
	    ready = true;
	else
	    ready = false;
    }

    public void getReady(){
	ready = true;
    }

    public void itemStateChanged(ItemEvent e){
	int iSel;

	//System.out.println("Ready = " + ready);

	if (ready){
	    if (myName.equals(JAGCTCheckBox.refGene)){
		myDomain.selectedGeneBefore[myI] = isSelected();
	    }else if (myName.equals(JAGCTCheckBox.refGroup)){
		myDomain.selectedGroup[myI][myJ] = isSelected();
		myJAGCTSelectionPane.myDomain.updateEnabledLigand();
		myJAGCTSelectionPane.myDomain.updateEnabledGene();
		myJAGCTSelectionPane.updateLFLigands();
		myJAGCTSelectionPane.updateLFGenes();
	    }else if (myName.equals(JAGCTCheckBox.refLigand)){
		myDomain.selectedLigand[myI] = isSelected();
		myJAGCTSelectionPane.myDomain.updateEnabledGene();
		myJAGCTSelectionPane.updateLFGenes();
	    }else
		Matrix.perror("JCheckBox's name is neither gene nor group\n");
	}
	myJAGCTSelectionPane.computeLabel();
	
	if (isSelected())
	    iSel = 1;
	else
	    iSel = 0;

	if (ready)
	    Scenario.add("JAGCTCheckBox_StateChanged", myName + Scenario.myToken + myI + Scenario.myToken + myJ + Scenario.myToken + iSel);
    }

}

class JAGCTSelectionPane extends JPanel implements JAGCTAbstractPane, ActionListener, Debuggable {

    static String DEFAULT_Selection_Token = "_Selection.sel";

    public static String Token_Prototype_Begin = new String ("@PROTOTYPE_BEGIN");
    public static String Token_Prototype_End = new String ("@PROTOTYPE_END");    

    JButton nextStep;
    AGCT myAGCT;
    Domain myDomain;
    String myName;
    ImageIcon myImageIcon;
    Box mainBox;
    JLabel selectLabel;
    JComboBox fm, fsm, psm;
    JButton nsf, nwc, nsp, wc_save; //number selected features, number wavelet coefficients, number selected prototypes, save wavelet coefficients

    Vector checkBoxGene, checkBoxLigand, checkBoxGroup;

    JCheckBox keepPrototypes;

    static String actionComboFeature = "JAGCTSelectionPane_combo_select";
    static String actionComboFeatureSelection = "JAGCTSelectionPane_feature_selection";
    static String actionComboPrototypeSelection = "JAGCTSelectionPane_prototype_selection";
    static String actionValidate = "JAGCTSelectionPane_validate";
    static String String_NSF = "Max. #features : ",
	String_NSP = "Max. #prototypes : ",
	String_NWC = "#wavelet coeff. : ";

    boolean selectionAvailable;
    //false when no dataset is loaded (no selection of ligands or genes   is possible)

    int xMargin, yMargin, nfeat;

    JAGCTSelectionPane(AGCT a, Domain d, String name, ImageIcon ic){
	checkBoxGene = new Vector();
	checkBoxLigand = new Vector();
	checkBoxGroup = new Vector();
	selectLabel = new JLabel();

	myName = name;
	mainBox = Box.createHorizontalBox();
	mainBox.setAlignmentX(LEFT_ALIGNMENT);

	selectionAvailable = false;

	myAGCT = a;
	myDomain = d;
	myImageIcon = ic;

	xMargin = 20;
	yMargin = 20;
	nfeat = 0;
    }

    public Rectangle getCaptureRectangle(){
	Rectangle bounds = getBounds();
        bounds.setLocation(getLocationOnScreen());
	return bounds;
    }

    public void toTrash(){
	removeAll();

	checkBoxGene = new Vector();
	checkBoxLigand = new Vector();
	checkBoxGroup = new Vector();
	selectLabel = new JLabel();

	mainBox = Box.createHorizontalBox();
	mainBox.setAlignmentX(LEFT_ALIGNMENT);

	selectionAvailable = false;

	xMargin = 20;
	yMargin = 20;
	nfeat = 0;
    }

    public int numberCurrentSelectableGenes(){
	int i, val = 0;
	for (i=0;i<myDomain.numberInitialGenes;i++)
	    if ( (myDomain.selectedGeneBefore[i]) && (myDomain.enabledGeneBefore[i]) )
		val ++;
	return val;
    }

    public int numberCurrentSelectedGenes(){
	int i, val = 0;
	for (i=0;i<myDomain.numberInitialGenes;i++)
	    if ( (myDomain.selectedGeneBefore[i]) && (myDomain.enabledGeneBefore[i]) )
		val ++;
	if ( (AGCT.Max_Number_Of_Prototypes > -1) && (AGCT.Max_Number_Of_Prototypes < val) )
	    val = AGCT.Max_Number_Of_Prototypes;
	return val;
    }

    public int numberCurrentSelectedFeatures(){
	int val = Feature.getNumberOfFeaturesBeforeSelection(myAGCT);
	if ( (AGCT.Max_Number_Of_Features > -1) && (AGCT.Max_Number_Of_Features < val) )
	    val = AGCT.Max_Number_Of_Features;
	return val;
    }

    public void computeLabel(){
	int i, s = 0;
	String v = "";
	v += "<HTML>Summary : current #genes: ";
	v += numberCurrentSelectedGenes() + " ; #time series: ";
	s = 0;
	for (i=0;i<myDomain.numberLigands;i++){
	    if ( (myDomain.selectedLigand[i]) && (myDomain.enabledLigand[i]) )
		s++;
	    //System.out.println("Ligand " + i + " : " + myDomain.selectedLigand[i] + " :: " + myDomain.enabledLigand[i]);
	}
	v += s + " ; #features: ";
	nfeat = numberCurrentSelectedFeatures();
	if (nfeat < 3)
	    v += "<FONT COLOR=FF0000>";
	v += nfeat;
	if (nfeat < 3)
	    v += "</FONT>";
	v += ".</HTML>";
	selectLabel.setText(v);
    }

    public void getReady(){
	int i;
	for (i=0;i<myDomain.numberLigands;i++)
	    ((JAGCTCheckBox) checkBoxLigand.elementAt(i)).getReady();

	if (!AGCT.No_GeneSelection)
	    for (i=0;i<myDomain.numberInitialGenes;i++)
		((JAGCTCheckBox) checkBoxGene.elementAt(i)).getReady();
    }

    public void updateLFLigands(){
	int i;
	for (i=0;i<myDomain.numberLigands;i++)
	    if (myDomain.enabledLigand[i])
		((JAGCTCheckBox) checkBoxLigand.elementAt(i)).setEnabled(true);
	    else
		((JAGCTCheckBox) checkBoxLigand.elementAt(i)).setEnabled(false);
    }

    public void updateLFGenes(){
	if (!AGCT.No_GeneSelection){
	    int i;
	    for (i=0;i<myDomain.numberInitialGenes;i++)
		if (myDomain.enabledGeneBefore[i])
		    ((JAGCTCheckBox) checkBoxGene.elementAt(i)).setEnabled(true);
		else
		    ((JAGCTCheckBox) checkBoxGene.elementAt(i)).setEnabled(false);
	}
    }

    public void setDomain(Domain d){
	myDomain = d;
    }

    public void paintComponent (Graphics g) {
	String defs = "";
	int yRefM = yMargin;

	if (AGCT.SUBMISSION)
	    defs += AGCT.SUBMISSION_TEXT;

	defs += "\n\nNo selection available: no data loaded";

	if (selectionAvailable == false) {
	    super.paintComponent(g);

	    for (String line : defs.split("\n"))
		g.drawString(line, xMargin, yRefM += g.getFontMetrics().getHeight());
	}
    }
    public void filterData(){
	myDomain.initializationFilters();
	selectionAvailable = true;
    }

    public void setButtonLabel(JButton j, String s, int v){
	j.setText(s + v);
    }

    public void displayComponents(){
	Box generalBox, groupBox, superBox;
	
	Box geneSelectionBox, ligandSelectionBox, superLigandSelectionBox, featureBox, prototypeBox, featureSelectionBox, prototypeSelectionBox, horizontalBoxNSF, horizontalBoxNSP, horizontalBoxNWC, horizontalBoxSelect;
	Box [] groupSelectionBox;
	//each element = a vertical Box containing the group elements selected

	String myBoxName;
	JAGCTCheckBox myBox;

	int id, i, j;
	String ss;
	boolean ok_groups = false;

	fm = new JComboBox(AGCT.Feature_Method);
	fm.setSelectedIndex(0);
	fm.setActionCommand(JAGCTSelectionPane.actionComboFeature);
	fm.addActionListener(this);

	nwc = new JButton("");
	setButtonLabel(nwc, JAGCTSelectionPane.String_NWC, AGCT.Number_Of_Wavelet_Stamps);
	nwc.setActionCommand("number_wavelet_coefficients");
	nwc.addActionListener(this);

	horizontalBoxNWC = Box.createHorizontalBox();
	horizontalBoxNWC.add(fm);
	horizontalBoxNWC.add(nwc);

	checkEnabledNWC();

	featureBox = Box.createVerticalBox();
	featureBox.setBorder(BorderFactory.createTitledBorder("1-Feature method"));
	featureBox.add(horizontalBoxNWC);

	fsm = new JComboBox(AGCT.Feature_Selection_Method);
	fsm.setSelectedIndex(0);
	fsm.setActionCommand(JAGCTSelectionPane.actionComboFeatureSelection);
	fsm.addActionListener(this);

	nsf = new JButton("");
	setButtonLabel(nsf, JAGCTSelectionPane.String_NSF, AGCT.Max_Number_Of_Features);
	nsf.setActionCommand("max_number_of_features");
	nsf.addActionListener(this);
	horizontalBoxNSF = Box.createHorizontalBox();
	horizontalBoxNSF.add(fsm);
	horizontalBoxNSF.add(nsf);

	checkEnabledNSF();

	featureSelectionBox = Box.createVerticalBox();
	featureSelectionBox.setBorder(BorderFactory.createTitledBorder("2-Automatic feature selection method"));
	featureSelectionBox.add(horizontalBoxNSF);

	psm = new JComboBox(AGCT.Prototype_Selection_Method);
	psm.setSelectedIndex(0);
	psm.setActionCommand(JAGCTSelectionPane.actionComboPrototypeSelection);
	psm.addActionListener(this);

	nsp = new JButton("");
	setButtonLabel(nsp, JAGCTSelectionPane.String_NSP, AGCT.Max_Number_Of_Prototypes);
	nsp.setActionCommand("max_number_of_prototypes");
	nsp.addActionListener(this);

	keepPrototypes = new JCheckBox("Keep");
	keepPrototypes.setToolTipText("Keep the prototypes computed in scenario");
	keepPrototypes.setEnabled(true);
	keepPrototypes.setSelected(true);

	horizontalBoxNSP = Box.createHorizontalBox();
	horizontalBoxNSP.add(psm);
	horizontalBoxNSP.add(nsp);
	horizontalBoxNSP.add(keepPrototypes);
	//horizontalBoxNSP.add(loadButton);
	//horizontalBoxNSP.add(saveButton);

	checkEnabledNSP();

	prototypeSelectionBox = Box.createVerticalBox();
	prototypeSelectionBox.setBorder(BorderFactory.createTitledBorder("3-Automatic prototype (gene) construction method"));
	prototypeSelectionBox.add(horizontalBoxNSP);

	nextStep = new JButton("<HTML><i>Validate</i> my selection below</HTML>");
	nextStep.setActionCommand(JAGCTSelectionPane.actionValidate);
	nextStep.setAlignmentX(Component.CENTER_ALIGNMENT);
	nextStep.addActionListener(this);

	wc_save = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/disk-haar-save.gif"))));
	wc_save.setToolTipText("save wavelet coefficients");
        wc_save.setActionCommand("wc_save");
	wc_save.addActionListener(this);
	wc_save.setEnabled(true);

	horizontalBoxSelect = Box.createHorizontalBox();
	horizontalBoxSelect.add(nextStep);
	horizontalBoxSelect.add(wc_save);

	geneSelectionBox = Box.createVerticalBox();
	geneSelectionBox.setBorder(BorderFactory.createTitledBorder("Genes"));
	geneSelectionBox.add(new JLabel("<HTML><i>By name</i></HTML>"));

	ss = "";

	if (!AGCT.No_GeneSelection)
	    for (i=0;i<myDomain.numberInitialGenes;i++){
		myBoxName = "gene " + ( (Gene) myDomain.domainGenes.elementAt(i)).name;
		myBox = new JAGCTCheckBox(this, myDomain, myBoxName, JAGCTCheckBox.refGene, i, -1);
		myBox.addItemListener(myBox);
		myBox.setSelected(myDomain.selectedGeneBefore[i]);
		checkBoxGene.addElement(myBox);
		geneSelectionBox.add(myBox); //Remove to prevent gene error
		if (i<AGCT.Max_Number_Of_Genes_Kept)
		    ss += ( (Gene) myDomain.domainGenes.elementAt(i)).toString();   
	    }
	else
	    geneSelectionBox.add(new JLabel("<HTML><b>Desactivated</b></HTML>"));
	
	myAGCT.myInformationFrame.setText(ss);

	ligandSelectionBox = Box.createVerticalBox();
	ligandSelectionBox.add(new JLabel("<HTML><i>By name</i></HTML>"));
	for (i=0;i<myDomain.numberLigands;i++){
	    myBoxName =  AGCTFileWriter.DATA_Ligand + " " + (String) ( (Vector) myDomain.domainLigands.elementAt(i)).elementAt(0);
	    myBox = new JAGCTCheckBox(this, myDomain, myBoxName, JAGCTCheckBox.refLigand, i, -1);
	    myBox.addItemListener(myBox);
	    myBox.setSelected(myDomain.selectedLigand[i]);
	    checkBoxLigand.addElement(myBox);
	    ligandSelectionBox.add(myBox);
	}

	if (myDomain.numberGroups >= 1)
	    ok_groups = true;

	groupBox = null;
	if (ok_groups){
	    groupBox = Box.createHorizontalBox();
	    groupBox.setBorder(BorderFactory.createTitledBorder("by group (AND)"));
	    groupBox.add(Box.createHorizontalStrut(15));
	    groupSelectionBox = new Box[myDomain.numberGroups];
	    for (i=0;i<myDomain.numberGroups;i++){
		groupSelectionBox[i] = Box.createVerticalBox();
		groupSelectionBox[i].setAlignmentY(TOP_ALIGNMENT);
		groupSelectionBox[i].add(new JLabel("<HTML>" + myDomain.domainGroups.elementAt(i) + "</HTML>"));
		checkBoxGroup.addElement(new Vector());
		for (j=0;j<myDomain.domainGroupElements[i].size();j++){
		    myBoxName = "Label " + myDomain.domainGroupElements[i].elementAt(j);
		    myBox = new JAGCTCheckBox(this, myDomain, myBoxName, JAGCTCheckBox.refGroup, i, j);
		    myBox.addItemListener(myBox);
		    myBox.setSelected(myDomain.selectedGroup[i][j]);
		    
		    ((Vector) checkBoxGroup.elementAt(i)).addElement(myBox);
		    groupSelectionBox[i].add(myBox);
		    if (j<myDomain.domainGroupElements[i].size()-1)
			groupSelectionBox[i].add(new JLabel("OR"));
		}
		groupBox.add(groupSelectionBox[i]);
		groupBox.add(Box.createHorizontalStrut(15));
	    }
	}

	geneSelectionBox.setAlignmentY(TOP_ALIGNMENT);
	ligandSelectionBox.setAlignmentY(TOP_ALIGNMENT);

	if (ok_groups)
	    groupBox.setAlignmentY(TOP_ALIGNMENT);

	superLigandSelectionBox = Box.createHorizontalBox();
	superLigandSelectionBox.setAlignmentX(LEFT_ALIGNMENT);
	superLigandSelectionBox.setBorder(BorderFactory.createTitledBorder("Time series"));
	superLigandSelectionBox.add(ligandSelectionBox);
	if (ok_groups)
	    superLigandSelectionBox.add(groupBox);

	generalBox = Box.createHorizontalBox();
	generalBox.add(geneSelectionBox);
	generalBox.add(superLigandSelectionBox);
	//generalBox.setAlignmentY(TOP_ALIGNMENT);
	//generalBox.setAlignmentX(LEFT_ALIGNMENT);

	if (ok_groups)
	    generalBox.setBorder(BorderFactory.createTitledBorder("4- User-fixed gene / time series / group selection"));
	else
	    generalBox.setBorder(BorderFactory.createTitledBorder("4- User-fixed gene / time series selection"));

	superBox = Box.createVerticalBox();
	//superBox.setAlignmentY(TOP_ALIGNMENT);
	//superBox.setAlignmentX(LEFT_ALIGNMENT);
	selectLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
	selectLabel.setSize(nextStep.getSize());


	superBox.add(horizontalBoxSelect);
	superBox.add(selectLabel);
	superBox.add(featureBox);
	superBox.add(featureSelectionBox);
	superBox.add(prototypeSelectionBox);
	superBox.add(generalBox);

	mainBox.add(superBox);
	mainBox.add(new Box.Filler(new Dimension(0, 0), new Dimension(1000, 0), new Dimension (1000, 0)));

        setLayout(new BorderLayout());
	add(new JScrollPane(mainBox), BorderLayout.CENTER);

	
	id = fm.getSelectedIndex();
	AGCT.Method_F = id;
	myDomain.updateEnabledGene();
	//updateLFGenes();
	computeLabel();

	//getReady() etait en commentaire, peut-etre pour la gestion des scenarii ?
	getReady();
	ControlProcess.put("displaySelectionPane",true);
    }
    
    public void offCheckBoxes(){
	int i, j;
	Vector v;
	if ( (checkBoxGene != null) && (checkBoxGene.size() > 0) ){
	    for (i=0;i<checkBoxGene.size();i++){
		( (JAGCTCheckBox) checkBoxGene.elementAt(i) ).setEnabled(false);
	    }
	}

	if ( (checkBoxLigand != null) && (checkBoxLigand.size() > 0) ){
	    for (i=0;i<checkBoxLigand.size();i++){
		( (JAGCTCheckBox) checkBoxLigand.elementAt(i) ).setEnabled(false);
	    }
	}	

	if ( (checkBoxGroup != null) && (checkBoxGroup.size() > 0) ){
	    for (i=0;i<checkBoxGroup.size();i++){
		v = (Vector) checkBoxGroup.elementAt(i);
		if ( ( v != null) && (v.size() > 0) ){
		    for (j=0;j<v.size();j++){
			( (JAGCTCheckBox) v.elementAt(j) ).setEnabled(false);
		    }
		}
	    }
	}	
    }

    public void checkEnabledNWC(){
	if (fm.getSelectedIndex() == 0)
	    nwc.setEnabled(false);
	else
	    nwc.setEnabled(true);
    }

    public void checkEnabledNSF(){
	if (fsm.getSelectedIndex() == 0)
	    nsf.setEnabled(false);
	else
	    nsf.setEnabled(true);
    }
    
    public void checkEnabledNSP(){
	if (psm.getSelectedIndex() == 0)
	    nsp.setEnabled(false);
	else
	    nsp.setEnabled(true);
    }
    
    public void finishValidate(){
	int i;

	ControlProcess.put("geneFeaturesComputed",true);
	myAGCT.computeFeatureNames();
	
	i=0;
	while(myDomain.selectedGeneBefore[i] == false)
	    i++;
	
	myAGCT.myInformationFrame.setText(((Gene) myDomain.domainGenes.elementAt(i)).coordinatesString());
	
	JClusteringProfileFrame.fillTime_Stamps_Summaries(myDomain);
    }

    public void goValidate(){

	fm.setEnabled(false);
	fsm.setEnabled(false);
	nsf.setEnabled(false);
	nwc.setEnabled(false);
	psm.setEnabled(false);
	nsp.setEnabled(false);
	keepPrototypes.setEnabled(false);

	if (keepPrototypes.isSelected())
	    Scenario.add("JAGCTSelectionPane_Keep_Prototypes");

	offCheckBoxes();

	myDomain.updateFilterData();
	myDomain.fillSelectedGenesAllAnnotations();
	
	Thread t = new Thread(){
		public void run(){
		    nextStep.setText("<HTML>Selection is <i>Validated</i></HTML>");
		    nextStep.setEnabled(false);
		}
	    };
	t.start();

	ControlProcess.put("filtersSelected",true);

	myDomain.featureAndPrototypeSelection();
    }
    
    //methods related to feature selection

    public void goComboFeatureSelection(ActionEvent e){
	int id;
	id = ((JComboBox)e.getSource()).getSelectedIndex();
	setMethod_FS(id);
    }

    public void playMethod_FS(int v){
	fsm.setSelectedIndex(v);
	setMethod_FS(v);
    }

    public void setMethod_FS(int v){
	AGCT.Method_FS = v;

	if (v == 0){
	    AGCT.Max_Number_Of_Features = -1;
	    setButtonLabel(nsf, JAGCTSelectionPane.String_NSF, AGCT.Max_Number_Of_Features);
	}

	checkEnabledNSF();
    }

    public int getFeatureSelectionIndex(){
	return fsm.getSelectedIndex();
    }

    public void moreMax_Number_Of_Features(){
	if (AGCT.Max_Number_Of_Features < 1){
	    AGCT.Max_Number_Of_Features = -1;
	    fsm.setSelectedIndex(0);
	}

	setButtonLabel(nsf, JAGCTSelectionPane.String_NSF, AGCT.Max_Number_Of_Features);
	computeLabel();
    }

    //methods related to prototype selection

    public void goComboPrototypeSelection(ActionEvent e){
	int id;
	id = ((JComboBox)e.getSource()).getSelectedIndex();
	setMethod_PS(id);
    }

    public void playMethod_PS(int v){
	psm.setSelectedIndex(v);
	setMethod_PS(v);
    }

    public void setMethod_PS(int v){
	AGCT.Method_PS = v;

	if (v == 0){
	    AGCT.Max_Number_Of_Prototypes = -1;
	    setButtonLabel(nsp, JAGCTSelectionPane.String_NSP, AGCT.Max_Number_Of_Prototypes);
	}

	checkEnabledNSP();
    }

    public int getPrototypeSelectionIndex(){
	return psm.getSelectedIndex();
    }

    public void moreMax_Number_Of_Prototypes(){
	if (AGCT.Max_Number_Of_Prototypes < 1){
	    AGCT.Max_Number_Of_Prototypes = -1;
	    psm.setSelectedIndex(0);
	}

	setButtonLabel(nsp, JAGCTSelectionPane.String_NSP, AGCT.Max_Number_Of_Prototypes);
	computeLabel();
    }

    // methods related to number of features

    public void moreNumber_Of_Wavelet_Stamps(){
	setButtonLabel(nwc, JAGCTSelectionPane.String_NWC, AGCT.Number_Of_Wavelet_Stamps);
	computeLabel();
    }

    public void moreMethod_F(){
	if ( (!Feature.checkAllWaveletStampsAutomatic(myAGCT, myAGCT.myDomain.domainTimes)) && (AGCT.Number_Of_Wavelet_Stamps == -1) ){
	    myAGCT.defaultNWS();
	    setButtonLabel(nwc, JAGCTSelectionPane.String_NWC, AGCT.Number_Of_Wavelet_Stamps);
	}
	
	checkEnabledNWC();
    }

    public void goComboFeature(ActionEvent e){
	int id;
	id = ((JComboBox)e.getSource()).getSelectedIndex();
	AGCT.Method_F = id;
	if ( (!Feature.checkAllWaveletStampsAutomatic(myAGCT, myAGCT.myDomain.domainTimes)) && (AGCT.Number_Of_Wavelet_Stamps == -1) )
	    Scenario.add("AGCT_Modification_Number_Of_Wavelet_Stamps",AGCT.DEFAULT_Number_Of_Wavelet_Stamps);
	setMethod_F(id);
    }

    public void setMethod_F(int id){
	moreMethod_F();
	myDomain.updateEnabledGene();
	updateLFGenes();
	computeLabel();
    }

    public void actionPerformed(ActionEvent e){
	String command = e.getActionCommand();
	if (command.equals(JAGCTSelectionPane.actionValidate)){
	    if (nfeat >= 3){
		goValidate();
		Scenario.add("JAGCTSelectionPane_Validate","");
		wc_save.setEnabled(true);
	    }
	}else if (command.equals(JAGCTSelectionPane.actionComboFeature)){
	    goComboFeature(e);
	    Scenario.add("JAGCTSelectionPane_Method_F", AGCT.Method_F + "");
	}else if (command.equals(JAGCTSelectionPane.actionComboFeatureSelection)){
	    goComboFeatureSelection(e);
	    Scenario.add("JAGCTSelectionPane_Feature_Selection_Method", AGCT.Method_FS + "");
	    if (AGCT.Method_FS == 0)
		Scenario.add("JAGCTSelectionPane_Max_Number_Of_Features", AGCT.Max_Number_Of_Features + "");
	}else if (command.equals(JAGCTSelectionPane.actionComboPrototypeSelection)){
	    goComboPrototypeSelection(e);
	    Scenario.add("JAGCTSelectionPane_Prototype_Selection_Method", AGCT.Method_PS + "");
	    if (AGCT.Method_PS == 0)
		Scenario.add("JAGCTSelectionPane_Max_Number_Of_Prototypes", AGCT.Max_Number_Of_Prototypes + "");
	}else if (command.equals("max_number_of_features")){
	    myAGCT.requestModificationMax_Number_Of_Features();
	    moreMax_Number_Of_Features();
	}else if (command.equals("max_number_of_prototypes")){
	    myAGCT.requestModificationMax_Number_Of_Prototypes();
	    moreMax_Number_Of_Prototypes();
	}else if (command.equals("number_wavelet_coefficients")){
	    myAGCT.requestModificationNumber_Of_Wavelet_Stamps();
	    moreNumber_Of_Wavelet_Stamps();
	}else if (command.equals("wc_save")){
	    myAGCT.requestSave_WC();
	}
    }
}
