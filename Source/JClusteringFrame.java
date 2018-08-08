import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;
import javax.imageio.ImageIO;

class KeyChange extends KeyAdapter {
    public static int DELTA = 50;
    public static double FACTOR = 2.0;

    JClusteringFrame myg;

    KeyChange(JClusteringFrame g){
	myg = g;
    }

    public void keyPressed(KeyEvent e) {
	
	if (myg.graphPanel.zoomAuthorized){
	    
	    switch (e.getKeyCode()) {
	    case KeyEvent.VK_ENTER:
		myg.graphPanel.magnification = 1.0;
		myg.graphPanel.Xoffset = 0;
		myg.graphPanel.Yoffset = 0;
		break;
	    case KeyEvent.VK_LEFT:
		myg.graphPanel.Xoffset += KeyChange.DELTA;
		break;
	    case KeyEvent.VK_RIGHT:
		myg.graphPanel.Xoffset -= KeyChange.DELTA;
		break;
	    case KeyEvent.VK_DOWN:
		myg.graphPanel.Yoffset -= KeyChange.DELTA;
		break;
	    case KeyEvent.VK_UP:
		myg.graphPanel.Yoffset += KeyChange.DELTA;
		break;
	    case KeyEvent.VK_P:
		myg.graphPanel.magnification *= KeyChange.FACTOR;
		break;
	    case KeyEvent.VK_M:
		myg.graphPanel.magnification /= KeyChange.FACTOR;
		break;
	    }


	    if (myg.listPlots.getSelectedIndex()>=0){
		myg.onePlot = true;
		myg.graphPanel.update = true;
		myg.graphPanel.repaint();
	    }
	}
    }
}

class JClusteringPlot extends JPanel implements MouseMotionListener, MouseListener, Debuggable{

    public static Color[] colors = {Color.blue,Color.red,Color.pink,Color.cyan,
				Color.black,Color.yellow};


    public static Color bgColor = Color.white;
    public static Color lineColor = Color.red;
    public static Color structColor = Color.green;
    public static Color pointColor = Color.blue;
    public static Color textColor = Color.pink;

    public static Color highlightColor = Color.magenta;

    public static int pointRadius = 3;
    public static int highlightRadius = 10;

    public static double leftMarginPercent = 10.0,
	rightMarginPercent = 10.0,
	upMarginPercent = 10.0,
	downMarginPercent = 10.0;

    JClusteringFrame myFrame;
    AGCT myAGCT;
    boolean update;

    boolean zoomAuthorized;

    double magnification;
    int Xoffset, Yoffset, pMouseX, pMouseY, cMouseX, cMouseY;

    private KeyChange listener;

    JClusteringPlot(JClusteringFrame jf, AGCT ma){
	myFrame = jf;
	myAGCT = ma;
	update = zoomAuthorized = false;
	magnification = 1.0;
	Xoffset = Yoffset = 0;
    }
    

    public void mousePressed (MouseEvent e) {
	pMouseX = e.getX();
	pMouseY = e.getY();
	requestFocus();
    }

    public void mouseEntered (MouseEvent e) {
	pMouseX = e.getX();
	pMouseY = e.getY();
    }
    
    public void mouseExited (MouseEvent e) {}

    public void mouseClicked (MouseEvent e) {
    }

    public void mouseReleased (MouseEvent e) {}


    public void mouseDragged(MouseEvent e){
	cMouseX = e.getX();
	cMouseY = e.getY();

	int deltaX = (cMouseX - pMouseX);
	int deltaY = (cMouseY - pMouseY);

	Xoffset -= deltaX;
	Yoffset -= deltaY;

	pMouseX = cMouseX;
	pMouseY = cMouseY;

	if (myFrame.listPlots.getSelectedIndex()>=0){
	    myFrame.onePlot = true;
	    update = true;
	    repaint();
	}
    }

    public void mouseMoved(MouseEvent e){
    }

    public Rectangle getCaptureRectangle(){
	Rectangle bounds = getBounds();
        bounds.setLocation(getLocationOnScreen());
	return bounds;
    }

    public void paintComponent (Graphics g) {
	if ( (myFrame.onePlot) && (myFrame.plotAuthorized) && (update) ){
	    plotClustering(g);
	}
	revalidate();
    }

    public int refX(AGCTClustering_Algorithm aa, int kORp){
	int ktarg;
	if ( (kORp == 1) && (aa.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP)) )
	    ktarg = ((AGCTClustering_AP) aa).vP;
	else
	    ktarg = aa.nclusters;
	return ktarg;
    }

    public double refY(AGCTClustering_Algorithm aa, int kORpORn, int dORs){
	//kORpORn = -K or -P or Clustering Number
	//dORs = Distortion or Similarity

	double ktarg;
	if ( (kORpORn == 1) && (aa.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP)) )
	    ktarg = (double) aa.nclusters;
	else if ( (kORpORn == 2) || (dORs == 0) )
	    ktarg = aa.returnFinalObjectiveFunction();
	else
	    ktarg = aa.kernelSim();
	return ktarg;
    }

    public Vector getData(String ss, int kORpORn, int dORs){
	/*Returns a vector over all clusterings ss containing in this order:
	  listx = new Integer[indices.size()]; : all X data values (k or p or number)
	  listy = new Double[indices.size()]; : all Y data (averages)
	  lists = new Double[indices.size()]; : all Sigmas over Y
	  listnb = new Integer [indices.size()]; : the # Y values for each X 
	  listminmax = new Vector(4); : xmin, xmax, ymin, ymax;
	  listValues = new Vector(); : all values as double for each X

	Tables ordered following increasing list?x*/
	
	// kORpORn = 0 for -K = X & distort = Y, 1 for -P = X & -K = Y, 2 for Clustering Number = X, distort = Y

	Vector ret = new Vector();
	Vector indices = new Vector();
	Vector minmax = new Vector();
	Vector allValues = null;
	AGCTClustering_Algorithm aa;
	int i, j, kcur, ktarg = -1, ni, refVal;
	boolean add;
	Integer idum;
	double dni, dni2, miny = -1.0, maxy = -1.0;
	
	for (i=0;i<myAGCT.allClusterings.size();i++){
	    aa = ( (Clustering) myAGCT.allClusterings.elementAt(i)).myClusteringAlgorithm;
	    if (aa.myReferenceName.equals(ss)){
		if (kORpORn < 2)
		    ktarg = refX(aa, kORpORn);
		else if (kORpORn == 2)
		    ktarg = i;
		else
		    Matrix.perror("JClusteringFrame.class :: bad value for kORpORn");
		add = true;
		j=0;
		if (indices.size()>0)
		    while ( (add == true) && (j<indices.size()) ){
			kcur = ( (Integer) indices.elementAt(j) ).intValue();
			if (kcur == ktarg)
			    add = false;
			else
			    j++;
		    }
		if (add)
		    indices.add(new Integer(ktarg));
	    }
	}

	Integer listx[] = new Integer[indices.size()];
	Double listy[] = new Double[indices.size()];
	Double lists[] = new Double[indices.size()];
	Integer listnb[] = new Integer[indices.size()];
	Double sumSquared[] = new Double[indices.size()];

	for (i=0;i<indices.size()-1;i++){
	    listx[i] = null;
	    listy[i] = null;
	    lists[i] = null;
	    listnb[i] = null;
	}

	if (kORpORn == 2){
	    allValues = new Vector();
	    for (i=0;i<indices.size();i++){
		aa = ( (Clustering) myAGCT.allClusterings.elementAt( ((Integer) indices.elementAt(i)).intValue() )).myClusteringAlgorithm;
		dni = refY(aa,kORpORn,dORs);
		listx[i] = new Integer(((Integer) indices.elementAt(i)).intValue());
		listy[i] = new Double (dni);
		lists[i] = new Double(0.0);
		listnb[i] = new Integer(1);
		allValues.addElement(new Integer(aa.nclusters));

		if ( (i==0) || (dni < miny) )
		    miny = dni;
		
		if ( (i==0) || (dni > maxy) )
		    maxy = dni;
	    }
	    minmax.add(new Integer(((Integer) indices.elementAt(0)).intValue()));
	    minmax.add(new Integer(((Integer) indices.elementAt(indices.size()-1)).intValue()));
	    minmax.add(new Double(miny));
	    minmax.add(new Double(maxy));
	}else if (kORpORn < 2){
	    for (i=0;i<indices.size()-1;i++)
		for (j=i+1;j<indices.size();j++)
		    if ( ( (Integer) indices.elementAt(i) ).intValue() > ( (Integer) indices.elementAt(j) ).intValue() ){
			idum = (Integer) indices.elementAt(i);
			indices.setElementAt(indices.elementAt(j),i);
			indices.setElementAt(idum,j);
		    }
	    
	    for (i=0;i<indices.size();i++)
		listx[i] = (Integer) indices.elementAt(i);
	    minmax.add(new Integer(listx[0].intValue()));
	    minmax.add(new Integer(listx[indices.size()-1].intValue()));
	    
	    for (j=0;j<indices.size();j++){
		ktarg = ( (Integer) indices.elementAt(j) ).intValue();
		for (i=0;i<myAGCT.allClusterings.size();i++){
		    aa = ( (Clustering) myAGCT.allClusterings.elementAt(i)).myClusteringAlgorithm;
		    refVal = refX(aa, kORpORn);
		    
		    if ( (aa.myReferenceName.equals(ss)) && (refVal == ktarg) ){
			
			if (allValues == null){
			    allValues = new Vector();
			    allValues.addElement(new Vector());
			}else if (allValues.size() == j)
			    allValues.addElement(new Vector());
			
			((Vector) allValues.elementAt(j)).addElement(new Double (refY(aa,kORpORn,dORs)));
			
			if (listnb[j] == null)
			    listnb[j] = new Integer(1);
			else{
			    ni = 1 + listnb[j].intValue();
			    listnb[j] = new Integer(ni);
			}
			
			if (listy[j] == null)
			    listy[j] = new Double (refY(aa,kORpORn,dORs));
			else{
			    dni = refY(aa,kORpORn,dORs) + listy[j].doubleValue();
			    listy[j] = new Double(dni);
			}
			
			if (sumSquared[j] == null)
			    sumSquared[j] = new Double (refY(aa,kORpORn,dORs) * refY(aa,kORpORn,dORs));
			else{
			    dni = (refY(aa,kORpORn,dORs) * refY(aa,kORpORn,dORs))
				+ sumSquared[j].doubleValue();
			    sumSquared[j] = new Double(dni);
			}
		    }
		}
	    }
	    
	    for (j=0;j<indices.size();j++){
		dni = (listy[j].doubleValue() / (double) listnb[j].intValue());
		listy[j] = new Double(dni);
		
		if ( (j==0) || (dni < miny) )
		    miny = dni;
		
		if ( (j==0) || (dni > maxy) )
		    maxy = dni;
		
		dni2 = (sumSquared[j].doubleValue() / (double) listnb[j].intValue()) - (dni * dni);
		dni = Math.sqrt(dni2);
		lists[j] = new Double(dni);
	    }
	    minmax.add(new Double(miny));
	    minmax.add(new Double(maxy));
	}
 
	ret.add(listx);
	ret.add(listy);
	ret.add(lists);
	ret.add(listnb);
	ret.add(minmax);
	ret.add(allValues);
	
	indices = null;
	sumSquared = null;

	return ret;
    }

    public void plotAllClusterings(Graphics g, int dORs){
	myAGCT.myInformationFrame.setText("Processing All Clusterings\nThis may take some time...");

	int [] rX = new int[2];
	double [] rY = new double[2];
	double rat;
	int i, j;
	AGCTClustering_Algorithm aa;
	boolean tested;
	for (i=0;i<myAGCT.allClusterings.size();i++){
	    aa = ( (Clustering) myAGCT.allClusterings.elementAt(i)).myClusteringAlgorithm;
	    if ( (i==0) || (aa.nclusters < rX[0]) )
		rX[0] = aa.nclusters;
	    if ( (i==0) || (aa.nclusters > rX[1]) )
		rX[1] = aa.nclusters;

	    rat = refY(aa,0,dORs);
	    if ( (i==0) || (rat < rY[0]) )
		rY[0] = rat;
	    if ( (i==0) || (rat > rY[1]) )
		rY[1] = rat;
	}

	for (j=0;j<Clustering.clusteringAlgorithms.length;j++){
	    tested = false;
	    for (i=0;i<myAGCT.allClusterings.size();i++){
		aa = ( (Clustering) myAGCT.allClusterings.elementAt(i)).myClusteringAlgorithm;
		if ( (!tested) && (aa.myReferenceName.equals(Clustering.clusteringAlgorithms[j])) ){
		    tested = true;
		    plotClustering(g, aa.myReferenceName, JClusteringPlot.colors[j], false, rX, rY, true, false, dORs);
		}
	    }
	}
	if (dORs == 0)
	    myFrame.updateBorder(JClusteringFrame.String_All_Distortions);
	else
	    myFrame.updateBorder(JClusteringFrame.String_All_Similarities);
	myAGCT.myInformationFrame.appendText(" done.");
    }

    public void initPanel(Graphics g){
	g.setColor(bgColor);
	super.paintComponent(g);
	g.fillRect(0,0,getWidth(),getHeight());
    }

    public void plotClustering(Graphics g, String kfc, Color col, boolean verticalLines, int [] referencesX, double [] referencesY, boolean useReferences, boolean makeStats, int dORs){

	int i, j, k, xx, yy, xx2, yy2, diff, si, sj, xmin, xmax, isiz, kORpORn;
	AGCTClustering_Algorithm aa;
	boolean keepIt;
	double dcur, dtcur, dxx, dyy, dsyy, dymin = -1, dymax = -1, dpt, delta, prob, dii, djj;
	String sX, sY;
	double listy[], lists[], datai[], dataj[];
	int listx[], listnb[];
	Vector data, setOfValues;

	if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM)){
	    diff = 3*pointRadius;
	    myFrame.updateBorder((String) myFrame.listPlots.getSelectedItem());
	    //myFrame.updateBorder(JClusteringFrame.String_KM);
	}else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM)){
	    diff = 3*pointRadius;
	    myFrame.updateBorder(JClusteringFrame.String_EM);
	}else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM)){
	    diff = 3*pointRadius;
	    myFrame.updateBorder(JClusteringFrame.String_NM);
	}else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP) ){
	    diff = 3*pointRadius;
	    myFrame.updateBorder(JClusteringFrame.String_CP);
	}else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)){
	    diff = 3*pointRadius;
	    myFrame.updateBorder(JClusteringFrame.String_HC_HC);
	}else if ( ((String) myFrame.listPlots.getSelectedItem()).equals(JClusteringFrame.String_AP_KP)){
	    diff = 3*pointRadius;
	    myFrame.updateBorder(JClusteringFrame.String_AP_KP);
	}else{
	    diff = 3*pointRadius;
	    myFrame.updateBorder(JClusteringFrame.String_AP);
	}

	kORpORn = 0;
	if ( (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP) )
	     && ( ((String) myFrame.listPlots.getSelectedItem()).equals(JClusteringFrame.String_AP_KP)) )
	    kORpORn = 1;
	else if ( ((String) myFrame.listPlots.getSelectedItem()).equals(JClusteringFrame.String_KM_Number) )
	    kORpORn = 2;

	if (kORpORn == 0){
	    sX = "k = ";
	    sY = "m = ";
	}else if (kORpORn == 1){
	    sX = "p = ";
	    sY = "k = ";
	}else{
	    sX = "#";
	    sY = "d = ";
	}

	data = getData(kfc,kORpORn,dORs);
	isiz = ((Integer [])data.elementAt(0)).length;
	listx = new int[isiz];
	listy = new double[isiz];
	lists = new double[isiz];
	listnb = new int[isiz];
	
	for (i=0;i<isiz;i++){
	    //System.out.print(" " + ((Integer []) data.elementAt(0))[i]);

	    listx[i] = ((Integer []) data.elementAt(0))[i].intValue();
	    listy[i] = ((Double []) data.elementAt(1))[i].doubleValue();
	    lists[i] = ((Double []) data.elementAt(2))[i].doubleValue();
	    listnb[i] = ((Integer []) data.elementAt(3))[i].intValue();
	}
	if (useReferences){
	    xmin = referencesX[0];
	    xmax = referencesX[referencesX.length-1];

	    dymin = referencesY[0];
	    dymax = referencesY[referencesY.length-1];
	}else{
	    xmin = ((Integer) ((Vector) data.elementAt(4)).elementAt(0)).intValue();
	    xmax = ((Integer) ((Vector) data.elementAt(4)).elementAt(1)).intValue();

	    dymin = ((Double) ((Vector) data.elementAt(4)).elementAt(2)).doubleValue();
	    dymax = ((Double) ((Vector) data.elementAt(4)).elementAt(3)).doubleValue();
	}
	setOfValues = (Vector) data.elementAt(5);
	data = null;
	
	//Structure
	
	g.setColor(structColor);
	g.setFont(JAGCTGraphicsPanel.SMALL_FONT);
	
	xx = (int) X(( (double) xmin ), ( (double) xmin ), ( (double) xmax ));
	yy = (int) Y(dymin, dymin, dymax);
	xx2 = (int) X(( (double) xmax ), ( (double) xmin ), ( (double) xmax ));
	yy2 = (int) Y(dymin, dymin, dymax);
	g.drawLine(xx, yy, xx2, yy2);
	
	xx2 = (int) X(( (double) xmin ), ( (double) xmin ), ( (double) xmax ));
	yy2 = (int) Y(dymax, dymin, dymax);
	g.drawLine(xx, yy, xx2, yy2);
	
	for (i=0;i<isiz;i++){
	    xx = (int) X(( (double) listx[i] ), ( (double) xmin ), ( (double) xmax ));
	    yy = (int) Y(dymin, dymin, dymax);

	    g.fillOval(xx - pointRadius, 
		       yy - pointRadius,
		       pointRadius + pointRadius,
		       pointRadius + pointRadius);
	    if ( listx[i] == xmin )
		g.drawString(sX + listx[i], xx - diff, yy + 2 * diff);
	    else
		g.drawString("" + listx[i], xx - diff, yy + 2 * diff);
	    
	    if (verticalLines){
		xx2 = (int) X(( (double) listx[i] ), ( (double) xmin ), ( (double) xmax ));
		yy2 = (int) Y(listy[i], dymin, dymax);
		g.drawLine(xx, yy, xx2, yy2);
	    }
	}
	
	//Data
	
	for (i=0;i<isiz;i++){
	    
	    g.setColor(col);
	    if (lists[i] < 0.0)
		Matrix.perror("Negative std dev.");
	    
	    if (lists[i] > 0.0){
		delta = 2*lists[i];
		
		xx = (int) X(( (double) listx[i] ), ( (double) xmin ), ( (double) xmax ));
		yy = (int) Y(listy[i] - delta, dymin, dymax);
		
		xx2 = (int) X(( (double) listx[i] ), ( (double) xmin ), ( (double) xmax ));
		yy2 = (int) Y(listy[i] + delta, dymin, dymax);
		
		g.drawLine(xx, yy, xx2, yy2);
		g.drawLine(xx - pointRadius, yy, xx + pointRadius, yy);
		g.drawLine(xx2 - pointRadius, yy2, xx2 + pointRadius, yy2);
		
	    }
	    xx = (int) X(( (double) listx[i] ), ( (double) xmin ), ( (double) xmax ));
	    yy = (int) Y(listy[i], dymin, dymax);

	    g.fillOval(xx - pointRadius, 
		       yy - pointRadius,
		       pointRadius + pointRadius,
		       pointRadius + pointRadius);

	    if ( (i==0) && (useReferences) ){
		g.setFont(JAGCTGraphicsPanel.HIGHLIGHT_FONT);
		g.drawString(kfc, xx, yy - diff);
	    }
	    
	    dyy = listy[i];
	    
	    if (!useReferences)
		g.setColor(textColor);
	    g.setFont(JAGCTGraphicsPanel.SMALL_FONT);
	    
	    if ( (kORpORn == 0) || (kORpORn == 2) )
		if (listnb[i]>1)
		    g.drawString(sY + DF.format(listy[i]) + " ; s = " + DF.format(lists[i]) + " (n = " + listnb[i] + ")", xx, yy);
		else if (kORpORn == 0)
		    g.drawString(sY + DF.format(listy[i]) + " (n = " + listnb[i] + ")", xx, yy);
		else
		    g.drawString(sY + DF.format(listy[i]) + " (k = " + setOfValues.elementAt(i) + ")", xx, yy);
	    else
		g.drawString(sY + DF.format(listy[i]), (int) X( (double) xmin, (double) xmin, (double) xmax) + diff, yy);
	    
	    //g.drawString("" + DF.format(dyy) + "(n = " + (int) ((Double) nbValues.elementAt(i)).doubleValue() + ")", xx, yy);
	}
	
	if (!useReferences)
	    g.setColor(lineColor);
	for (i=0;i<isiz-1;i++){
	    xx = (int) X(( (double) listx[i] ), ( (double) xmin ), ( (double) xmax ));
	    yy = (int) Y(listy[i], dymin, dymax);
	    
	    xx2 = (int) X(( (double) listx[i+1] ), ( (double) xmin ), ( (double) xmax ));
	    yy2 = (int) Y(listy[i+1], dymin, dymax);
	    
	    g.drawLine(xx, yy, xx2, yy2);
	}
	
	//statistical tests
	
	if ( (kORpORn == 0) && (makeStats) )
	    for (i=0;i<isiz;i++){
		dii = listy[i];
		keepIt = true;
		j=0;
		
		while ( (keepIt) && (j<i) ){
		    djj = listy[j];
		    if (j!=i){
			si = listnb[i];
			sj = listnb[j];
			
			datai = new double[si];
			dataj = new double[sj];
			
			for (k=0;k<si;k++)
			    datai[k] = ( (Double) ( (Vector) setOfValues.elementAt(i) ).elementAt(k) ).doubleValue();
			for (k=0;k<sj;k++)
			    dataj[k] = ( (Double) ( (Vector) setOfValues.elementAt(j) ).elementAt(k) ).doubleValue();
			prob = Statistics.tutest(datai, dataj);
			
			//System.out.println("i = " + i + ", dii = " + dii + ", j = " + j + ", djj = " + djj + ", prob = " + prob);
			
			if ( ( kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM) ) || ( kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP) ) || ( kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP) ) || ( kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC) ) || ( kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM) ) )
			    if ( (djj < dii) || (prob > Statistics.LIMIT_P_CHI2) )
				keepIt = false;
			
			if ( kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM) )
			    if ( (djj > dii) || (prob > Statistics.LIMIT_P_CHI2) )
				keepIt = false;
		    }
		    j++;
		}
		
		if ( listnb[i] == 1.0 )
		    keepIt = false;
		
		if (keepIt){
		    g.setColor(highlightColor);
		    xx = (int) X(( (double) listx[i] ), ( (double) xmin ), ( (double) xmax ));
		    yy = (int) Y(listy[i], dymin, dymax);
		    
		    g.drawArc(xx - highlightRadius, 
			      yy - highlightRadius,
			      highlightRadius + highlightRadius,
			      highlightRadius + highlightRadius, 0, 360);
		}
	    }
	    
	data = setOfValues = null;
    }

    public void plotClustering(Graphics g){
	initPanel(g);

	int ind;
	String kfc = "", hcs = "";
	boolean ok = true;

	try{
	    ind = myFrame.listPlots.getSelectedIndex();
	    kfc = (String) myFrame.keysForCombo.elementAt(ind);
	    hcs = (String) myFrame.listPlots.getSelectedItem();
	}catch (NullPointerException e){
	    ok = false;
	}

	if (ok){
	    if ( (zoomAuthorized == false) && (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)) )
		zoomAuthorized = true;
	    if ( (zoomAuthorized == true) && (!kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)) )
		zoomAuthorized = false;

	    if ( (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM))
		 || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM))
		 || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP))
		 || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP)) 
		 || (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM)) )
		plotClustering(g, kfc, pointColor, true, null, null, false, true, 0);
	    else if (hcs.equals(JClusteringFrame.String_HC_HC))
		plotClustering(g, kfc, pointColor, true, null, null, false, false, 0);
	    else if (kfc.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC))
		plotClusteringHC(g);
	    else if ( (kfc.equals(AGCTClustering_Algorithm.ALL_CLUSTERS)) && (hcs.equals(JClusteringFrame.String_All_Distortions)) )
		plotAllClusterings(g,0);
	    else if ( (kfc.equals(AGCTClustering_Algorithm.ALL_CLUSTERS)) && (hcs.equals(JClusteringFrame.String_All_Similarities)) )
		plotAllClusterings(g,1);
	    update = false;
	}
    }

    public void plotClusteringHC(Graphics g){
	myFrame.updateBorder(JClusteringFrame.String_HC);

	int ind = myFrame.listPlots.getSelectedIndex(), 
	    indhc = ( (Integer) myFrame.indexHC.elementAt(ind) ).intValue(),
	    i, j, k, g1, g2, tper, cli, clj, minp = -1;

	double xx, yy, oxx, oyy, xx1, oxx1, yy1, oyy1, xx2, oxx2, yy2, oyy2, nxx, onxx,
	    minDisplayX = 0.0, maxDisplayX = (double) getWidth(), minDisplayY = 0.0, maxDisplayY = (double) getHeight();

	double dxmax = ( (AGCTClustering_HC) ( (Clustering) myAGCT.allClusterings.elementAt(indhc)).myClusteringAlgorithm).totdistort;
	double delta = dxmax * 0.10;
	dxmax += delta;

	boolean first = true;

	int nclust = ( (Clustering) myAGCT.allClusterings.elementAt(indhc)).nclusters;

	double dxmin = 0.0;
	double dymax = (double) myAGCT.myDomain.numberSelectedGenes;
	double dymin = 0.0;

	double dxx, dyy, dxx1, dyy1, dxx2, dyy2;

	Vector allX = new Vector(), allY = new Vector();
	Vector dummyHierarchy = ( (AGCTClustering_HC) ( (Clustering) myAGCT.allClusterings.elementAt(indhc)).myClusteringAlgorithm).dummyHierarchy;
	Vector nonDummyHierarchy = ( (AGCTClustering_HC) ( (Clustering) myAGCT.allClusterings.elementAt(indhc)).myClusteringAlgorithm).nonDummyHierarchy;
	Cluster ccur, ci, cj;
	Vector v1, v2;

	int [] place = new int [myAGCT.myDomain.numberSelectedGenes];

	Matrix sf = ( (AGCTClustering_HC) ( (Clustering) myAGCT.allClusterings.elementAt(indhc)).myClusteringAlgorithm).soft_memberships;
	Gene gg;
	boolean isTop;

	for (i=0;i<myAGCT.myDomain.numberSelectedGenes;i++)
	    place[i] = i;

	for (i=dummyHierarchy.size()-1;i>=0;i--){
	    ccur = (Cluster) dummyHierarchy.elementAt(i);
	    if ( (ccur.indexPrevious1 != -1) || (ccur.indexPrevious2 != -1) ){
		v1 = ( (Cluster) dummyHierarchy.elementAt(ccur.indexPrevious1) ).allGenes;
		v2 = ( (Cluster) dummyHierarchy.elementAt(ccur.indexPrevious2) ).allGenes;

		for (j=0;j<v1.size();j++){
		    g1 = ( (Integer) v1.elementAt(j) ).intValue();

		    for (k=0;k<v2.size();k++){
			g2 = ( (Integer) v2.elementAt(k) ).intValue();
			if ( (k==0) || (place[minp] > place[g2]) )
			    minp = g2;
		    }

		    if (place[g1] > place[minp]){
			tper = place[g1];
			place[g1] = place[minp];
			place[minp] = tper;
		    }
		}
	    }
	}

	for (i=0;i<nonDummyHierarchy.size();i++){
	    ccur = (Cluster) nonDummyHierarchy.elementAt(i);
	    isTop = ccur.top;
	    dxx = dxmax - ccur.curdist;
	    
	    if ( (ccur.indexPrevious1 == -1) || (ccur.indexPrevious2 == -1) ){
		dyy = (double) place[i];
	    }else{
		dyy = ( (Double) allY.elementAt(ccur.indexPrevious1) ).doubleValue()
		    + ( (Double) allY.elementAt(ccur.indexPrevious2) ).doubleValue();
		dyy /= 2.0;
		
		g.setColor(structColor);
		
		dxx1 = ( (Double) allX.elementAt(ccur.indexPrevious1) ).doubleValue();
		dyy1 = ( (Double) allY.elementAt(ccur.indexPrevious1) ).doubleValue();
		xx = X(dxx, dxmin, dxmax);
		xx1 = X(dxx1, dxmin, dxmax);
		yy1 = Y(dyy1, dymin, dymax);

		oxx = scaledX(xx,magnification,Xoffset);
		oxx1 = scaledX(xx1,magnification,Xoffset);
		oyy1 = scaledY(yy1,magnification,Yoffset);

		g.setColor(( (Clustering) myAGCT.allClusterings.elementAt(indhc)).colorCluster[ccur.clusterIndex]);

		// horizontal line down
		g.drawLine((int) oxx, (int) oyy1, (int) oxx1, (int) oyy1);
		
		dxx2 = ( (Double) allX.elementAt(ccur.indexPrevious2) ).doubleValue();
		dyy2 = ( (Double) allY.elementAt(ccur.indexPrevious2) ).doubleValue();
		xx2 = X(dxx2, dxmin, dxmax);
		yy2 = Y(dyy2, dymin, dymax);

		oxx2 = scaledX(xx2,magnification,Xoffset);
		oyy2 = scaledY(yy2,magnification,Yoffset);

		// horizontal line up
		g.drawLine((int) oxx, (int) oyy2, (int) oxx2, (int) oyy2);
		
		// vertical line
		g.drawLine((int) oxx, (int) oyy1, (int) oxx, (int) oyy2);

		if ( (first) || ( oxx < minDisplayX ) )
		    minDisplayX = oxx;
		if ( (first) || ( oxx > maxDisplayX ) )
		    maxDisplayX = oxx;

		if ( (first) || ( oxx1 < minDisplayX ) )
		    minDisplayX = oxx1;
		if ( (first) || ( oxx1 > maxDisplayX ) )
		    maxDisplayX = oxx1;

		if ( (first) || ( oxx2 < minDisplayX ) )
		    minDisplayX = oxx2;
		if ( (first) || ( oxx2 > maxDisplayX ) )
		    maxDisplayX = oxx2;

		if ( (first) || ( oyy1 < minDisplayY ) )
		    minDisplayY = oyy1;
		if ( (first) || ( oyy1 > maxDisplayY ) )
		    maxDisplayY = oyy1;

		if ( (first) || ( oyy2 < minDisplayY ) )
		    minDisplayY = oyy2;
		if ( (first) || ( oyy2 > maxDisplayY ) )
		    maxDisplayY = oyy2;
	    }
	    
	    allX.add(new Double(dxx));
	    allY.add(new Double(dyy));
	    
	    xx = X(dxx, dxmin, dxmax);
	    yy = Y(dyy, dymin, dymax);

	    oxx = scaledX(xx,magnification,Xoffset);
	    oyy = scaledY(yy,magnification,Yoffset);

	    if ( (first) || ( oxx < minDisplayX ) )
		minDisplayX = oxx;
	    if ( (first) || ( oxx > maxDisplayX ) )
		maxDisplayX = oxx;

	    if ( (first) || ( oyy < minDisplayY ) )
		minDisplayY = oyy;
	    if ( (first) || ( oyy > maxDisplayY ) )
		maxDisplayY = oyy;

	    if (isTop){
		g.setColor(( (Clustering) myAGCT.allClusterings.elementAt(indhc)).colorCluster[ccur.clusterIndex]);
		g.setFont(JAGCTGraphicsPanel.SMALL_FONT);
		nxx = X(dxmin, dxmin, dxmax);

		onxx = scaledX(nxx,magnification,Xoffset);

		if ( (first) || ( onxx < minDisplayX ) )
		    minDisplayX = onxx;
		if ( (first) || ( onxx > maxDisplayX ) )
		    maxDisplayX = onxx;

		g.drawLine((int) onxx, (int) oyy, (int) oxx, (int) oyy);
		g.drawString(Clustering.Center_Name + ccur.clusterIndex, (int) onxx, (int) oyy - (pointRadius));
	    }

	    g.setColor(( (Clustering) myAGCT.allClusterings.elementAt(indhc)).colorCluster[ccur.clusterIndex]);

	    //g.setColor(pointColor);
	    //g.fillOval(xx - pointRadius, 
	    //	       yy - pointRadius,
	    //       pointRadius + pointRadius,
	    //       pointRadius + pointRadius);
	    
	    if ( (ccur.indexPrevious1 == -1) || (ccur.indexPrevious2 == -1) ){
		gg = (Gene) myAGCT.myDomain.domainGenes.elementAt(myAGCT.myDomain.selectedGeneNumberToGeneNumber[i]);
		
		g.setFont(JAGCTGraphicsPanel.SMALL_FONT);
		g.drawString(gg.name + "(" + gg.asciiName + ")", (int) oxx, (int) oyy + pointRadius);
	    }

	    first = false;
	}

	int dX = getWidth() / 20;
	int dY = getHeight() / 20;
	int lim = 5;
	
	double ratX = (double) getWidth() / ( (double) maxDisplayX - (double) minDisplayX );
	if (ratX > 1.0)
	    ratX = 1.0;
	double ratY = (double) getHeight() / ( (double) maxDisplayY - (double) minDisplayY );
	if (ratY > 1.0)
	    ratY = 1.0;

	ratX *= (double) dX;
	ratY *= (double) dY;

	if (ratX < lim)
	    ratX = lim;

	if (ratY < lim)
	    ratY = lim;

	double off = 15.0;
	
	int pX = (int) ( ( ( - minDisplayX + ((double) getWidth()/2.0)) / ( (double) maxDisplayX - (double) minDisplayX ) ) * (double) dX );
	
	int pY = (int) ( ( ( - minDisplayY + ((double) getHeight()/2.0)) / ( (double) maxDisplayY - (double) minDisplayY ) ) * (double) dY );

	double mX, mY, MX, MY;
	mX = (double) pX - (ratX / 2.0);
	MX = (double) pX + (ratX / 2.0);
	mY = (double) pY - (ratY / 2.0);
	MY = (double) pY + (ratY / 2.0);
	
	if (mX < 0.0) 
	    mX = 0.0;
	if (mX > dX) 
	    mX = dX;
	if (MX < 0.0) 
	    MX = 0.0;
	if (MX > dX) 
	    MX = dX;

	if (mY < 0.0) 
	    mY = 0.0;
	if (mY > dY) 
	    mY = dY;
	if (MY < 0.0) 
	    MY = 0.0;
	if (MY > dY) 
	    MY = dY;

	Graphics2D gg2 = (Graphics2D) g;

	g.setColor(Color.black);

	gg2.fill(new java.awt.geom.Rectangle2D.Double(off + 0.0, off + 0.0,
						      (double) dX,
						      (double) dY));

	g.setColor(Color.white);

	gg2.fill(new java.awt.geom.Rectangle2D.Double(off + mX, off + mY, MX - mX, MY - mY));


	g.setColor(Color.black);

	gg2.draw(new java.awt.geom.Rectangle2D.Double(off + 0.0, off + 0.0,
						      (double) dX,
						      (double) dY));

	if (pX < 0)
	    pX = 0;

	if (pY < 0)
	    pY = 0;

	if (pX > dX)
	    pX = dX;

	if (pY > dY)
	    pY = dY;

	/*g.fillOval((int) off + pX - pointRadius, 
		   (int) off + pY - pointRadius,
	           pointRadius + pointRadius,
	           pointRadius + pointRadius);*/
    }

    public double scaledX(double xx, double mag, int Xoffset) {
	return  (((((double) xx - ((double) getWidth() / 2.0)) * mag) + (double) getWidth() / 2.0) - (double) Xoffset);
    }	

    public double scaledY(double yy, double mag, int Yoffset) {
	return  (((((double) yy - ((double) getHeight() / 2.0)) * mag) + (double) getHeight() / 2.0) - (double) Yoffset);
    }	

    public double X(double xx, double xmin, double xmax){
	double w = ( (double) getWidth() * (100.0 - leftMarginPercent - rightMarginPercent) / 100.0 );
	double delta = ( (double) getWidth() * leftMarginPercent) / 100.0;
	double v = delta + ( ( (xx - xmin) * w ) / (xmax - xmin) );
	return v;
    }

    public double Y(double yy, double ymin, double ymax){
	double h = ( (double) getHeight() * (100.0 - upMarginPercent - downMarginPercent) / 100 );
	double delta = ( (double) getHeight() * upMarginPercent) / 100.0;
	double v = h + delta - ( ( (yy - ymin) * h ) / (ymax - ymin) );
	return v;
    }
}

class JClusteringFrame extends JFrame{
    
    public static String Default_Border = "(no plot selected)";
    public static String Default_Combo[] = {"(no clustering)"};
    public static String String_KM = "All K-Means Distortions = f(K) (regardless of -Point)";
    public static String String_KM_Number = "All K-Means Distortions = f(Clustering Number) (regardless of -Point)";
    public static String String_EM = "All EM Log-likelihoods = f(K) (regardless of -Point)";
    public static String String_CP = "All CP Distortions = f(K) (regardless of -Point)";
    public static String String_AP = "All AP Distortions = f(K) (regardless of -Point)";
    public static String String_NM = "All NM Distortions = f(K) (regardless of -Point)";
    public static String String_HC_HC = "All HC Distortions = f(K) (regardless of -Point)";
    public static String String_HC = "Dendrogram for Hierarchical Clustering ";
    public static String String_AP_KP = "K = f(P) for All AP";
    public static String String_All_Distortions = "All Distortions = f(K) (regardless of -Point, for all algorithms)";
    public static String String_All_Similarities = "Kernel similarities = f(K) (initial space, for all algorithms)";

    String borderString;
    JComboBox listPlots;
    JClusteringPlot graphPanel;
    JButton goButton, cameraButton;
    Box selectionBox;

    Vector keysForCombo;
    //contains the vector of reference clustering class names for the combo indices
    
    Vector indexHC;
    //contains Integers = index of the HC in the Clusterings in AGCT

    AGCT myAGCT;
    boolean onePlot, plotAuthorized, oneAlgo;

    public void flushAll(){
	listPlots.removeAllItems();
	listPlots.addItem(JClusteringFrame.Default_Combo);
	keysForCombo = null;
	indexHC = null;
	listPlots.setEnabled(false);
	onePlot = oneAlgo = false;
	plotAuthorized = false;
    }

    public void captureAndSave(){
	Rectangle rect = graphPanel.getCaptureRectangle();
	BufferedImage screencapture = null, framecapture = null;
	try{
	    framecapture = new Robot().createScreenCapture(rect);
	}catch(AWTException e){}
	
	JFileChooser chooser = new JFileChooser();
	ExampleFileFilter filter = new ExampleFileFilter();
	filter.addExtension("png");
	filter.setDescription(".png Files");
	chooser.setFileFilter(filter);
	chooser.setApproveButtonText("Save");
	int returnVal = chooser.showSaveDialog(this);
	if( (chooser.getSelectedFile() != null) && (returnVal == JFileChooser.APPROVE_OPTION) ){
	    try{
		ImageIO.write(framecapture, "png", chooser.getSelectedFile());
	    }catch(IOException e){}
	    myAGCT.myInformationFrame.setText("Saving clustering pane to file: " +
				       chooser.getSelectedFile().getName());
	    ControlProcess.put("frameCapturedClustering",true);
	}
    }

    public void displayInfo(){

	cameraButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/camera.png"))));
	cameraButton.setToolTipText("capture the visible clustering plot");
        cameraButton.setActionCommand("capture");
	cameraButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){
		    captureAndSave();
		}
	    });
	cameraButton.setEnabled(false);

	goButton = new JButton(new ImageIcon(Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/go.png"))));
	goButton.setToolTipText("run clustering");
        goButton.setActionCommand("clustering");
	goButton.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){
		    if (listPlots.getSelectedIndex()>=0){
			onePlot = true;
			graphPanel.update = true;
			graphPanel.repaint();
		    }
		}
	    });
	goButton.setEnabled(false);

	graphPanel = new JClusteringPlot(this, myAGCT);
	graphPanel.setBorder(BorderFactory.createTitledBorder(JClusteringFrame.Default_Border));

	graphPanel.addMouseListener(graphPanel);
	graphPanel.addMouseMotionListener(graphPanel);


	listPlots = new JComboBox(JClusteringFrame.Default_Combo);
	listPlots.setToolTipText("select which clustering data to plot");
	listPlots.setSelectedIndex(0);
	listPlots.setEnabled(false);

	KeyChange listener = new KeyChange(this);
	cameraButton.addKeyListener(listener);
	goButton.addKeyListener(listener);
	listPlots.addKeyListener(listener);
	graphPanel.addKeyListener(listener);

	selectionBox = Box.createHorizontalBox();
	selectionBox.add(cameraButton);
	selectionBox.add(listPlots);
	selectionBox.add(goButton);

        Container pane = getContentPane();
        pane.setLayout(new BorderLayout());
	pane.add(selectionBox, BorderLayout.NORTH);
	pane.add(graphPanel, BorderLayout.CENTER);

	addWindowListener(new FermetureListener("Closing AGCT's ClusteringFrame\n"));

	Image img = Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/information.png"));
	setIconImage( img );
    }

    public void updateBorder(String t){
	graphPanel.setBorder(BorderFactory.createTitledBorder(t));
    }

    public void updateClusteringLF(){
	int i, nhc = 0;
	AGCTClustering_Algorithm aa;
	flushAll();
	requestFocus();

	if (myAGCT.allClusterings != null){
	    listPlots.removeAllItems();
	    if (!isVisible()){
		setVisible(true);
		myAGCT.myFrame.toFront();
	    }

	    listPlots.setEnabled(true);
	    keysForCombo = new Vector();
	    indexHC = new Vector();
	    cameraButton.setEnabled(true);
	    goButton.setEnabled(true);

	    for (i=0;i<myAGCT.allClusterings.size();i++){
		aa = ( (Clustering) myAGCT.allClusterings.elementAt(i)).myClusteringAlgorithm;
		if (oneAlgo == false){
		    keysForCombo.add(AGCTClustering_Algorithm.ALL_CLUSTERS);
		    listPlots.addItem(JClusteringFrame.String_All_Distortions);
		    keysForCombo.add(AGCTClustering_Algorithm.ALL_CLUSTERS);
		    listPlots.addItem(JClusteringFrame.String_All_Similarities);
		    oneAlgo = true;
		    indexHC.add(new Integer(-1));
		    indexHC.add(new Integer(-1));
		}
		if ( (!aa.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)) && (!keysForCombo.contains(aa.myReferenceName)) ){
		    keysForCombo.add(aa.myReferenceName);
		    indexHC.add(new Integer(-1));
		    if (aa.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_KM)){
			listPlots.addItem(JClusteringFrame.String_KM);
			keysForCombo.add(aa.myReferenceName);
			indexHC.add(new Integer(-1));
			listPlots.addItem(JClusteringFrame.String_KM_Number);
		    }else if (aa.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_EM))
			listPlots.addItem(JClusteringFrame.String_EM);
		    else if (aa.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_CP))
			listPlots.addItem(JClusteringFrame.String_CP);
		    else if (aa.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_NM))
			listPlots.addItem(JClusteringFrame.String_NM);
		    else if (aa.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_AP)){
			listPlots.addItem(JClusteringFrame.String_AP);

			indexHC.add(new Integer(-1));
			keysForCombo.add(aa.myReferenceName);
			listPlots.addItem(JClusteringFrame.String_AP_KP);
		    }
		}
	    }

	    oneAlgo = false;
	    for (i=0;i<myAGCT.allClusterings.size();i++){
		aa = ( (Clustering) myAGCT.allClusterings.elementAt(i)).myClusteringAlgorithm;
		if (aa.myReferenceName.equals(AGCTClustering_Algorithm.REFERENCE_NAME_HC)){
		    if (!oneAlgo){
			keysForCombo.add(aa.myReferenceName);
			listPlots.addItem(JClusteringFrame.String_HC_HC);
			indexHC.add(new Integer(-1));
			oneAlgo = true;
		    }

		    keysForCombo.add(aa.myReferenceName);
		    listPlots.addItem(JClusteringFrame.String_HC + nhc);
		    indexHC.add(new Integer(i));
		    nhc++;
		}
	    }

	    plotAuthorized = true;
	}

    }


    JClusteringFrame(String name, AGCT c){
	super(name);
	setSize(AGCT.WindowWidth,600);
	myAGCT = c;
	displayInfo();
	setVisible(false);
	setResizable(true);
	onePlot = false;

        addComponentListener(new ComponentAdapter() {
		public void componentResized(ComponentEvent evt) {
		    graphPanel.update = true;
		    if (onePlot)
			graphPanel.repaint();
		}

		public void componentMoved(ComponentEvent evt) {
		    graphPanel.update = true;
		    if (onePlot)
			graphPanel.repaint();
		}
	    });
    }
}