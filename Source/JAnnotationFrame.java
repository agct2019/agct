import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

class JAnnotationCheckBox extends JCheckBox{
    static String leftB = "<b>", rightB = "</b>";

    String type;
    boolean light; //true iff string belongs to Selected prototypes

    JAnnotationCheckBox(String n, String t){
	super(n);
	type = t;
	light = false;
    }

    public void changeLF(boolean b){
	String t = getText();
	light = true;
	if (b)
	    setText("<HTML>" + JAnnotationCheckBox.leftB + t + JAnnotationCheckBox.rightB + "</HTML>");
    }
}

class JAnnotationFrame extends JFrame implements ActionListener, ItemListener{
    
    static String anyOfOrNoneFString = "Any function -OR None-", 
	anyOfOrNonePString = "Any process -OR None-",
	anyOfOrNoneCString = "Any compartment -OR None-",
	anyFString = "Any function",
	anyPString = "Any process",
	anyCString = "Any compartment",
	checkedFString = "Checked functions",
	checkedPString = "Checked processes",
	checkedCString = "Checked compartments",
	notCheckedFString = "Unchecked functions",
	notCheckedPString = "Unchecked processes",
	notCheckedCString = "Unchecked compartments",
	andString,
	orString;

    static String switchPrototypes = "JAnnotationFrame_switchPrototype", showCard = "show Cardinals", rasPString = "Clear P", rasFString = "Clear F", rasCString = "Clear C";

    static String TAG_F = "Biological Function", TAG_P = "Biological Process", TAG_C = "Cellular Compartment";

    JButton andButton, rasP, rasF, rasC;
    JRadioButton anyOfOrNoneF, anyOfOrNoneP, anyOfOrNoneC, 
	anyF, anyP, anyC, 
	checkedF, checkedP, checkedC, 
	notCheckedF, notCheckedP, notCheckedC;
    Vector boxF, boxP, boxC;

    boolean bAnyOfOrNoneF, bAnyOfOrNoneP, bAnyOfOrNoneC,
	bAnyP, bAnyF, bAnyC,
	bCheckedF, bCheckedP, bCheckedC,
	bNotCheckedF, bNotCheckedP, bNotCheckedC,
	bAnd;

    AGCT myAGCT;
    Domain myDomain;

    int nCheckF, nCheckP, nCheckC,
	nNotCheckedF, nNotCheckedP, nNotCheckedC;

    JCheckBox onlyPrototypes, showNumber;

    Box selectionBoxP, selectionBoxF, selectionBoxC;
    Hashtable stringToBoxP, stringToBoxF, stringToBoxC;

    public static String makeAndOrString(String ao, boolean f, boolean p, boolean c){
	String ret = ao + " (logical combinations of ";
	if (f)
	    ret +="functions";
	if (f && (p || c))
	    ret += ", ";
	if (p)
	    ret +="processes";
	if ((f || p) && c)
	    ret += ", ";
	if (c)
	    ret +="compartments";
	ret += ")";
	return ret;
    }

    public static void initAndString(boolean f, boolean p, boolean c){
	andString = makeAndOrString("AND", f, p, c);
    }

    public static void initOrString(boolean f, boolean p, boolean c){
	orString = makeAndOrString("OR", f, p, c);
    }

    JAnnotationFrame(String name, AGCT c, Domain d){
	super(name);

	bAnyOfOrNoneF = bAnyOfOrNoneP = bAnyOfOrNoneC = bAnd = true;
	bAnyP = bAnyF = bAnyC = bCheckedF = bCheckedP = bCheckedC = bNotCheckedF = bNotCheckedP = bNotCheckedC = false;

	myAGCT = c;
	myDomain = d;

	nCheckF = nCheckP = nCheckC = nNotCheckedF = nNotCheckedP = nNotCheckedC = -1;
	stringToBoxP = stringToBoxF = stringToBoxC = null;

	displayInfo();
	setVisible(true);
	setResizable(true);
    }

    public static boolean geneIsDisplayed(AGCT ag, Gene gg){
	Domain myDomain = ag.myDomain;

	boolean bF = false, bP = false, bC = false;
	boolean bNotFound = false, val;
	JAnnotationCheckBox jc;
	String st;
	int i, j;

	if (ag.myAnnotationFrame.bAnyOfOrNoneF)
	    bF = true;
	else{
	    if (ag.myAnnotationFrame.boxF != null){
		bF = false;
		if ( (gg.annotationsF != null) &&
		     ( ( (ag.myAnnotationFrame.nCheckF == 0) && (ag.myAnnotationFrame.bCheckedF) ) 
		       || ( (ag.myAnnotationFrame.nNotCheckedF == 0) && (ag.myAnnotationFrame.bNotCheckedF) ) ) )
		    bF = true;
		else if ( (gg.annotationsF != null) 
		     && (ag.myAnnotationFrame.boxF != null) ){
		    i = 0;
		    do{
			bNotFound = true;
			if (!((String[]) gg.annotationsF.elementAt(i))[0].equals(Domain.ANNOTATION_F))
			    Matrix.perror("JAnnotationFrame.class :: annotation != F in F annotations");
			st = ((String[]) gg.annotationsF.elementAt(i))[1];
			jc = (JAnnotationCheckBox) ag.myAnnotationFrame.boxF.elementAt( ( (Integer) myDomain.accessToSelectedGenesAllAnnotationsF.get(st)).intValue() );
			if ( ( (jc.isSelected()) && (ag.myAnnotationFrame.bCheckedF) )
			     || ( (!jc.isSelected()) && (ag.myAnnotationFrame.bNotCheckedF) ) )
			    bNotFound = false;
			if (bNotFound)
			    i++;
		    }while( (bNotFound) && (i<gg.annotationsF.size()) );

		    if (bNotFound)
			bF = false;
		    else
			bF = true;
		}
	    }
	}

	if (ag.myAnnotationFrame.bAnyOfOrNoneP)
	    bP = true;
	else{
	    if (ag.myAnnotationFrame.boxP != null){
		bP = false;
		if ( (gg.annotationsP != null) &&
		     ( ( (ag.myAnnotationFrame.nCheckP == 0) && (ag.myAnnotationFrame.bCheckedP) ) 
		       || ( (ag.myAnnotationFrame.nNotCheckedP == 0) && (ag.myAnnotationFrame.bNotCheckedP) ) ) )
		    bP = true;
		else if ( (gg.annotationsP != null) 
		     && (ag.myAnnotationFrame.boxP != null) ){
		    i = 0;

		    do{
			bNotFound = true;
			if (!((String[]) gg.annotationsP.elementAt(i))[0].equals(Domain.ANNOTATION_P))
			    Matrix.perror("JAnnotationFrame.class :: annotation != P in P annotations");
			st = ((String[]) gg.annotationsP.elementAt(i))[1];
			jc = (JAnnotationCheckBox) ag.myAnnotationFrame.boxP.elementAt( ( (Integer) myDomain.accessToSelectedGenesAllAnnotationsP.get(st)).intValue() );
			if ( ( (jc.isSelected()) && (ag.myAnnotationFrame.bCheckedP) )
			     || ( (!jc.isSelected()) && (ag.myAnnotationFrame.bNotCheckedP) ) )
			    bNotFound = false;
			if (bNotFound)
			    i++;
		    }while( (bNotFound) && (i<gg.annotationsP.size()) );

		    if (bNotFound)
			bP = false;
		    else
			bP = true;
		}
	    }
	}

	if (ag.myAnnotationFrame.bAnyOfOrNoneC)
	    bC = true;
	else{
	    if (ag.myAnnotationFrame.boxC != null){
		bC = false;
		if ( (gg.annotationsC != null) &&
		     ( ( (ag.myAnnotationFrame.nCheckC == 0) && (ag.myAnnotationFrame.bCheckedC) ) 
		       || ( (ag.myAnnotationFrame.nNotCheckedC == 0) && (ag.myAnnotationFrame.bNotCheckedC) ) ) )
		    bC = true;
		else if ( (gg.annotationsC != null) 
		     && (ag.myAnnotationFrame.boxC != null) ){
		    i = 0;

		    do{
			bNotFound = true;
			if (!((String[]) gg.annotationsC.elementAt(i))[0].equals(Domain.ANNOTATION_C))
			    Matrix.perror("JAnnotationFrame.class :: annotation != C in C annotations");
			st = ((String[]) gg.annotationsC.elementAt(i))[1];
			jc = (JAnnotationCheckBox) ag.myAnnotationFrame.boxC.elementAt( ( (Integer) myDomain.accessToSelectedGenesAllAnnotationsC.get(st)).intValue() );
			if ( ( (jc.isSelected()) && (ag.myAnnotationFrame.bCheckedC) )
			     || ( (!jc.isSelected()) && (ag.myAnnotationFrame.bNotCheckedC) ) )
			    bNotFound = false;
			if (bNotFound)
			    i++;
		    }while( (bNotFound) && (i<gg.annotationsC.size()) );

		    if (bNotFound)
			bC = false;
		    else
			bC = true;
		}
	    }
	}

	if (ag.myAnnotationFrame.bAnd)
	    val = bF && bP && bC;
	else
	    val = bF || bP || bC;

	return val;
    }

    public static int numberOfDisplayedGenes(AGCT ag, int id){
	if (ag.myAnnotationFrame == null)
	    return 1;

	Vector l;
	int k, tot = 0;
	Gene gg = (Gene) ag.myDomain.domainGenes.elementAt(ag.myDomain.selectedGeneNumberToGeneNumber[id]);

	if ( (Prototype.No_Reduction) || (ag.myAnnotationFrame.onlyPrototypes.isSelected()) ){
	    if (geneIsDisplayed(ag, gg))
		tot = 1;
	}else{
	    l = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(ag.myDomain.selectedGeneNumberToGeneNumber[id]));
	    for (k=0;k<l.size();k++){
		gg = (Gene) ag.myDomain.domainGenes.elementAt(((Integer) l.elementAt(k)).intValue());
		if (geneIsDisplayed(ag, gg))
		    tot = tot + 1;
	    }
	}

	return tot;
    }


    public static boolean prototypeIsDisplayed(AGCT ag, int id){
	if (ag.myAnnotationFrame == null)
	    return true;

	Vector l;
	int k;
	Gene gg = (Gene) ag.myDomain.domainGenes.elementAt(ag.myDomain.selectedGeneNumberToGeneNumber[id]);

	if ( (Prototype.No_Reduction) || (ag.myAnnotationFrame.onlyPrototypes.isSelected()) )
	    return geneIsDisplayed(ag, gg);
	else{
	    l = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(ag.myDomain.selectedGeneNumberToGeneNumber[id]));
	    for (k=0;k<l.size();k++){
		gg = (Gene) ag.myDomain.domainGenes.elementAt(((Integer) l.elementAt(k)).intValue());
		if (geneIsDisplayed(ag, gg)){
		    return true;
		}
	    }
	}
	return false;
    }

    public void modifyP(String s){
	if (s.equals(anyPString)){
	    bAnyP = true;
	    bAnyOfOrNoneP = false;
	    bCheckedP = false;
	    bNotCheckedP = false;
	}else if (s.equals(anyOfOrNonePString)){
	    bAnyP = false;
	    bAnyOfOrNoneP = true;
	    bCheckedP = false;
	    bNotCheckedP = false;
	}else if (s.equals(checkedPString)){
	    bAnyP = false;
	    bAnyOfOrNoneP = false;
	    bCheckedP = true;
	    bNotCheckedP = false;
	}else if (s.equals(notCheckedPString)){
	    bAnyP = false;
	    bAnyOfOrNoneP = false;
	    bCheckedP = false;
	    bNotCheckedP = true;
	}
	myAGCT.repaint();
    }

    public void modifyF(String s){
	if (s.equals(anyFString)){
	    bAnyF = true;
	    bAnyOfOrNoneF = false;
	    bCheckedF = false;
	    bNotCheckedF = false;
	}else if (s.equals(anyOfOrNoneFString)){
	    bAnyF = false;
	    bAnyOfOrNoneF = true;
	    bCheckedF = false;
	    bNotCheckedF = false;
	}else if (s.equals(checkedFString)){
	    bAnyF = false;
	    bAnyOfOrNoneF = false;
	    bCheckedF = true;
	    bNotCheckedF = false;
	}else if (s.equals(notCheckedFString)){
	    bAnyF = false;
	    bAnyOfOrNoneF = false;
	    bCheckedF = false;
	    bNotCheckedF = true;
	}
	myAGCT.repaint();
    }

    public void modifyC(String s){
	if (s.equals(anyCString)){
	    bAnyC = true;
	    bAnyOfOrNoneC = false;
	    bCheckedC = false;
	    bNotCheckedC = false;
	}else if (s.equals(anyOfOrNoneCString)){
	    bAnyC = false;
	    bAnyOfOrNoneC = true;
	    bCheckedC = false;
	    bNotCheckedC = false;
	}else if (s.equals(checkedCString)){
	    bAnyC = false;
	    bAnyOfOrNoneC = false;
	    bCheckedC = true;
	    bNotCheckedC = false;
	}else if (s.equals(notCheckedCString)){
	    bAnyC = false;
	    bAnyOfOrNoneC = false;
	    bCheckedC = false;
	    bNotCheckedC = true;
	}
	myAGCT.repaint();
    }



    public void displayInfo(){
	int i;
	JCheckBox jc;

	anyOfOrNoneF = new JRadioButton(anyOfOrNoneFString); 
	anyOfOrNoneP = new JRadioButton(anyOfOrNonePString); 
	anyOfOrNoneC = new JRadioButton(anyOfOrNoneCString); 
	checkedF = new JRadioButton(checkedFString); 
	checkedP = new JRadioButton(checkedPString); 
	checkedC = new JRadioButton(checkedCString); 
	notCheckedF = new JRadioButton(notCheckedFString); 
	notCheckedP = new JRadioButton(notCheckedPString); 
	notCheckedC = new JRadioButton(notCheckedCString); 

	ButtonGroup groupF = new ButtonGroup(), groupP = new ButtonGroup(), groupC = new ButtonGroup();
	
	groupF.add(anyOfOrNoneF);
	groupF.add(checkedF);
	groupF.add(notCheckedF);

	groupP.add(anyOfOrNoneP);
	groupP.add(checkedP);
	groupP.add(notCheckedP);

	groupC.add(anyOfOrNoneC);
	groupC.add(checkedC);
	groupC.add(notCheckedC);

	JAnnotationFrame.initAndString(myDomain.someAnnotationsF, myDomain.someAnnotationsP, myDomain.someAnnotationsC);
	JAnnotationFrame.initOrString(myDomain.someAnnotationsF, myDomain.someAnnotationsP, myDomain.someAnnotationsC);

	andButton = new JButton(JAnnotationFrame.andString);
	andButton.setActionCommand(JAGCTSelectionPane.actionValidate);
        andButton.setAlignmentX(Component.CENTER_ALIGNMENT);
	andButton.addActionListener(this);

	onlyPrototypes = new JCheckBox("Only prototypes");
	onlyPrototypes.setActionCommand(JAnnotationFrame.switchPrototypes);
	onlyPrototypes.setSelected(true);
	onlyPrototypes.setEnabled(false);
	onlyPrototypes.addActionListener(this);

	showNumber = new JCheckBox("Show Card");
	showNumber.setActionCommand(JAnnotationFrame.showCard);
	showNumber.setSelected(false);
	showNumber.setEnabled(true);
	showNumber.addActionListener(this);

	Box generalBox = Box.createHorizontalBox();

	if (myDomain.someAnnotationsP){

	    rasP = new JButton(JAnnotationFrame.rasPString);
	    rasP.setActionCommand(JAnnotationFrame.rasPString);
	    rasP.setAlignmentX(Component.LEFT_ALIGNMENT);
	    rasP.addActionListener(this);

	    boxP = new Vector();
	    stringToBoxP = new Hashtable();
	    for (i=0;i<myDomain.selectedGenesAllAnnotationsP.size();i++){
		jc = new JAnnotationCheckBox((String) myDomain.selectedGenesAllAnnotationsP.elementAt(i), Domain.ANNOTATION_P);
		jc.setFont(JAGCTGraphicsPanel.SMALL_FONT);
		jc.addItemListener(this);
		jc.setSelected(false);
		boxP.add(jc);

		stringToBoxP.put(new String(jc.getText()), jc);
	    }

	    Box buttonBoxP = Box.createVerticalBox();
	    buttonBoxP.add(anyOfOrNoneP);
	    buttonBoxP.add(checkedP);
	    buttonBoxP.add(notCheckedP); 
	    buttonBoxP.add(rasP);

	    anyOfOrNoneP.addActionListener(new ActionListener() {
		    public void actionPerformed( ActionEvent e )
		    {
			modifyP(anyOfOrNoneP.getText());
		    }
		});
	    anyOfOrNoneP.setSelected(true);

	    checkedP.addActionListener(new ActionListener() {
		    public void actionPerformed( ActionEvent e )
		    {
			modifyP(checkedP.getText());
		    }
		});

	    notCheckedP.addActionListener(new ActionListener() {
		    public void actionPerformed( ActionEvent e )
		    {
			modifyP(notCheckedP.getText());
		    }
		});

	    selectionBoxP = Box.createVerticalBox();
	    selectionBoxP.setAlignmentY(TOP_ALIGNMENT);
	    selectionBoxP.setAlignmentX(LEFT_ALIGNMENT);


	    nCheckP = 0;
	    nNotCheckedP = boxP.size();
	    selectionBoxP.setBorder(BorderFactory.createTitledBorder(getTitleP()));
 
	    for (i=0;i<myDomain.selectedGenesAllAnnotationsP.size();i++){
		selectionBoxP.add((JCheckBox) boxP.elementAt(i));
	    }

	    Box sBoxP = Box.createVerticalBox();
	    sBoxP.setAlignmentY(TOP_ALIGNMENT);
	    sBoxP.add(buttonBoxP);
	    sBoxP.add(selectionBoxP);

	    generalBox.add(sBoxP);

	}

	if (myDomain.someAnnotationsF){

	    rasF = new JButton(JAnnotationFrame.rasFString);
	    rasF.setActionCommand(JAnnotationFrame.rasFString);
	    rasF.setAlignmentX(Component.LEFT_ALIGNMENT);
	    rasF.addActionListener(this);

	    boxF = new Vector();
	    stringToBoxF = new Hashtable();

	    for (i=0;i<myDomain.selectedGenesAllAnnotationsF.size();i++){
		jc = new JAnnotationCheckBox((String) myDomain.selectedGenesAllAnnotationsF.elementAt(i), Domain.ANNOTATION_F);
		jc.setFont(JAGCTGraphicsPanel.SMALL_FONT);
		jc.addItemListener(this);
		jc.setSelected(false);
		boxF.add(jc);

		stringToBoxF.put(new String(jc.getText()), jc);
	    }

	    Box buttonBoxF = Box.createVerticalBox();
	    buttonBoxF.add(anyOfOrNoneF);
	    buttonBoxF.add(checkedF);
	    buttonBoxF.add(notCheckedF); 
	    buttonBoxF.add(rasF);

	    anyOfOrNoneF.addActionListener(new ActionListener() {
		    public void actionPerformed( ActionEvent e )
		    {
			modifyF(anyOfOrNoneF.getText());
		    }
		});
	    anyOfOrNoneF.setSelected(true);

	    checkedF.addActionListener(new ActionListener() {
		    public void actionPerformed( ActionEvent e )
		    {
			modifyF(checkedF.getText());
		    }
		});

	    notCheckedF.addActionListener(new ActionListener() {
		    public void actionPerformed( ActionEvent e )
		    {
			modifyF(notCheckedF.getText());
		    }
		});

	    selectionBoxF = Box.createVerticalBox();
	    selectionBoxF.setAlignmentY(TOP_ALIGNMENT);
	    selectionBoxF.setAlignmentX(LEFT_ALIGNMENT);

	    nCheckF = 0;
	    nNotCheckedF = boxF.size();
	    selectionBoxF.setBorder(BorderFactory.createTitledBorder(getTitleF()));
	    selectionBoxF.add(buttonBoxF);

	    for (i=0;i<myDomain.selectedGenesAllAnnotationsF.size();i++){
		selectionBoxF.add((JCheckBox) boxF.elementAt(i));
	    }

	    Box sBoxF = Box.createVerticalBox();
	    sBoxF.setAlignmentY(TOP_ALIGNMENT);
	    sBoxF.add(buttonBoxF);
	    sBoxF.add(selectionBoxF);

	    generalBox.add(sBoxF);

	}

	if (myDomain.someAnnotationsC){

	    rasC = new JButton(JAnnotationFrame.rasCString);
	    rasC.setActionCommand(JAnnotationFrame.rasCString);
	    rasC.setAlignmentX(Component.LEFT_ALIGNMENT);
	    rasC.addActionListener(this);

	    boxC = new Vector();
	    stringToBoxC = new Hashtable();

	    for (i=0;i<myDomain.selectedGenesAllAnnotationsC.size();i++){
		jc = new JAnnotationCheckBox((String) myDomain.selectedGenesAllAnnotationsC.elementAt(i), Domain.ANNOTATION_C);
		jc.setFont(JAGCTGraphicsPanel.SMALL_FONT);
		jc.addItemListener(this);
		jc.setSelected(false);
		boxC.add(jc);

		stringToBoxC.put(new String(jc.getText()), jc);
	    }

	    Box buttonBoxC = Box.createVerticalBox();
	    buttonBoxC.add(anyOfOrNoneC);
	    buttonBoxC.add(checkedC);
	    buttonBoxC.add(notCheckedC); 
	    buttonBoxC.add(rasC);

	    anyOfOrNoneC.addActionListener(new ActionListener() {
		    public void actionPerformed( ActionEvent e )
		    {
			modifyC(anyOfOrNoneC.getText());
		    }
		});
	    anyOfOrNoneC.setSelected(true);

	    checkedC.addActionListener(new ActionListener() {
		    public void actionPerformed( ActionEvent e )
		    {
			modifyC(checkedC.getText());
		    }
		});

	    notCheckedC.addActionListener(new ActionListener() {
		    public void actionPerformed( ActionEvent e )
		    {
			modifyC(notCheckedC.getText());
		    }
		});

	    selectionBoxC = Box.createVerticalBox();
	    selectionBoxC.setAlignmentY(TOP_ALIGNMENT);
	    selectionBoxC.setAlignmentX(LEFT_ALIGNMENT);

	    nCheckC = 0;
	    nNotCheckedC = boxC.size();
	    selectionBoxC.setBorder(BorderFactory.createTitledBorder(getTitleC()));
	    selectionBoxC.add(buttonBoxC);

	    for (i=0;i<myDomain.selectedGenesAllAnnotationsC.size();i++){
		selectionBoxC.add((JCheckBox) boxC.elementAt(i));
	    }

	    Box sBoxC = Box.createVerticalBox();
	    sBoxC.setAlignmentY(TOP_ALIGNMENT);
	    sBoxC.add(buttonBoxC);
	    sBoxC.add(selectionBoxC);

	    generalBox.add(sBoxC);

	}


	Box superBox = Box.createVerticalBox();
	superBox.setBorder(BorderFactory.createTitledBorder("Annotations"));
	superBox.add(generalBox);

	Box aboveBox = Box.createHorizontalBox();
	aboveBox.add(onlyPrototypes);
	aboveBox.add(andButton);
	aboveBox.add(showNumber);

        Container pane = getContentPane();
        pane.setLayout(new BorderLayout());
	pane.add(aboveBox, BorderLayout.NORTH);
	pane.add(new JScrollPane(superBox), BorderLayout.CENTER);

	addWindowListener(new FermetureListener("Closing AGCT's AnnotationFrame\n"));

	Image img = Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/information.png"));
	setIconImage(img);
    }

    public String getTitleP(){
	return "[" + nCheckP + "|" + nNotCheckedP + "]";
    }

    public String getTitleF(){
	return "[" + nCheckF + "|" + nNotCheckedF + "]";
    }

    public String getTitleC(){
	return "[" + nCheckC + "|" + nNotCheckedC + "]";
    }

    public void actionPerformed(ActionEvent e){
	String command = e.getActionCommand();
	if (command.equals(JAGCTSelectionPane.actionValidate)){
	    if (bAnd){
		bAnd = false;
		andButton.setText(JAnnotationFrame.orString);
	    }else{
		bAnd = true;
		andButton.setText(JAnnotationFrame.andString);
	    }
	    myAGCT.repaint();
	}else if (command.equals(JAnnotationFrame.showCard)){
	    myAGCT.repaint();
	}else if (command.equals(JAnnotationFrame.rasPString)){
	    rasP();
	}else if (command.equals(JAnnotationFrame.rasFString)){
	    rasF();
	}else if (command.equals(JAnnotationFrame.rasCString)){
	    rasC();
	}
    }

    public void rasP(){
	JAnnotationCheckBox jc;
	int i;
	for (i=0;i<boxP.size();i++){
	    jc = (JAnnotationCheckBox) boxP.elementAt(i);
	    jc.setSelected(false);
	}
	nCheckP = 0;
    }

    public void rasF(){
	JAnnotationCheckBox jc;
	int i;
	for (i=0;i<boxF.size();i++){
	    jc = (JAnnotationCheckBox) boxF.elementAt(i);
	    jc.setSelected(false);
	}
	nCheckF = 0;
    }

    public void rasC(){
	JAnnotationCheckBox jc;
	int i;
	for (i=0;i<boxC.size();i++){
	    jc = (JAnnotationCheckBox) boxC.elementAt(i);
	    jc.setSelected(false);
	}
	nCheckC = 0;
    }

    public void itemStateChanged(ItemEvent e){
	JAnnotationCheckBox jc = (JAnnotationCheckBox) e.getSource();

	if (jc.isSelected()){
	    if (jc.type.equals(Domain.ANNOTATION_P)){
		nCheckP ++;
		nNotCheckedP --;
		selectionBoxP.setBorder(BorderFactory.createTitledBorder(getTitleP()));
	    }else if (jc.type.equals(Domain.ANNOTATION_F)){
		nCheckF ++;
		nNotCheckedF --;
		selectionBoxF.setBorder(BorderFactory.createTitledBorder(getTitleF()));
	    }else if (jc.type.equals(Domain.ANNOTATION_C)){
		nCheckC ++;
		nNotCheckedC --;
		selectionBoxC.setBorder(BorderFactory.createTitledBorder(getTitleC()));
	    }else
		Matrix.perror("CheckBox type not recognized");

	    if ( ( (nCheckP<0) || (nNotCheckedP<0) ) && (myAGCT.annotationsExistP) )
		Matrix.perror("Negative number of CheckBoxes(P)");
	    if ( ( (nCheckF<0) || (nNotCheckedF<0) ) && (myAGCT.annotationsExistF) )
		Matrix.perror("Negative number of CheckBoxes(F)");
	    if ( ( (nCheckC<0) || (nNotCheckedC<0) ) && (myAGCT.annotationsExistC) )
		Matrix.perror("Negative number of CheckBoxes(C)");

	}else{
	    if (jc.type.equals(Domain.ANNOTATION_P)){
		nCheckP --;
		nNotCheckedP ++;
		selectionBoxP.setBorder(BorderFactory.createTitledBorder(getTitleP()));
	    }else if (jc.type.equals(Domain.ANNOTATION_F)){
		nCheckF --;
		nNotCheckedF ++;
		selectionBoxF.setBorder(BorderFactory.createTitledBorder(getTitleF()));
	    }else if (jc.type.equals(Domain.ANNOTATION_C)){
		nCheckC --;
		nNotCheckedC ++;
		selectionBoxC.setBorder(BorderFactory.createTitledBorder(getTitleC()));
	    }else
		Matrix.perror("CheckBox type not recognized");

	    if ( ( (nCheckP<0) || (nNotCheckedP<0) ) && (myAGCT.annotationsExistP) )
		Matrix.perror("Negative number of CheckBoxes(P)");
	    if ( ( (nCheckF<0) || (nNotCheckedF<0) ) && (myAGCT.annotationsExistF) )
		Matrix.perror("Negative number of CheckBoxes(F)");
	    if ( ( (nCheckC<0) || (nNotCheckedC<0) ) && (myAGCT.annotationsExistC) )
		Matrix.perror("Negative number of CheckBoxes(C)");

	}

	myAGCT.repaint();
    }

    public void changeLF_Prototype_Selection(){
	if (!Prototype.No_Reduction){
	    int i, j;
	    Gene gg;
	    JAnnotationCheckBox ja;
	    String s;
	    if ( (myDomain.someAnnotationsF) || (myDomain.someAnnotationsP) || (myDomain.someAnnotationsC) ){
		AGCTCounter cc = new AGCTCounter(myAGCT.myInformationFrame, "Completing Annotations", myDomain.numberSelectedGenes);
		for (i=0;i<myDomain.numberSelectedGenes;i++){
		    gg = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[i]);
		    if (gg.annotationsP != null)
			for (j=0;j<gg.annotationsP.size();j++){
			    s = ((String[]) gg.annotationsP.elementAt(j))[0];
			    if (!s.equals(Domain.ANNOTATION_P))
				Matrix.perror("JAnnotationFrame.class :: annotation != P in P annotations for prototypes");
			    ja = (JAnnotationCheckBox) stringToBoxP.get(((String[]) gg.annotationsP.elementAt(j))[1]);
				if (!ja.light)
				    ja.changeLF(true);
			}
		    if (gg.annotationsF != null)
			for (j=0;j<gg.annotationsF.size();j++){
			    s = ((String[]) gg.annotationsF.elementAt(j))[0];
			    if (!s.equals(Domain.ANNOTATION_F))
				Matrix.perror("JAnnotationFrame.class :: annotation != F in F annotations for prototypes");
			    ja = (JAnnotationCheckBox) stringToBoxF.get(((String[]) gg.annotationsF.elementAt(j))[1]);
				if (!ja.light)
				    ja.changeLF(true);
			}
		    if (gg.annotationsC != null)
			for (j=0;j<gg.annotationsC.size();j++){
			    s = ((String[]) gg.annotationsC.elementAt(j))[0];
			    if (!s.equals(Domain.ANNOTATION_C))
				Matrix.perror("JAnnotationFrame.class :: annotation != C in C annotations for prototypes");
			    ja = (JAnnotationCheckBox) stringToBoxC.get(((String[]) gg.annotationsC.elementAt(j))[1]);
				if (!ja.light)
				    ja.changeLF(true);
			}
		    cc.increment();
		}
		cc.end();
	    }
	}
    }
}