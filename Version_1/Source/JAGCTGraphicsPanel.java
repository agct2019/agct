import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import java.util.*;
import java.text.DecimalFormat;

class JAGCTGraphicsPanel extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener, ComponentListener, Debuggable {
    

    public static int DISPLAY_K_DIMENSION_ESTIMATION = 2;
    // 0 - No display of the dimension
    // 1 - Plots with pies, angle proportional to the max local dimension
    // 2 - Displays the rounding of the dimension

    static double RATIO_ZOOM = 1.2;
    static int ALPHA_MIN = 2;
    static int OFFSET_MIN = 10;
    static int SPLIT_POINTS = 20;
    static final int SENSITIVITY = 3;

    static final float MIN_STROKE = 1.0F;
    static final float MAX_STROKE = 1.5F;

    static final int VERY_SMALL_FONT_SIZE = 6;
    static final int SMALL_FONT_SIZE = 9;
    static final int NORMAL_FONT_SIZE = 12;
    static final int HIGHLIGHT_FONT_SIZE = 20;
    static final int SEARCH_FONT_SIZE = 20;

    static final Font VERY_SMALL_FONT = new Font("Times Roman",Font.ITALIC,VERY_SMALL_FONT_SIZE);
    static final Font SMALL_FONT = new Font("Times Roman",Font.PLAIN,SMALL_FONT_SIZE);
    static final Font NORMAL_FONT = new Font("Times Roman",Font.PLAIN,NORMAL_FONT_SIZE);
    static final Font HIGHLIGHT_FONT = new Font("Times Roman",Font.BOLD,HIGHLIGHT_FONT_SIZE);
    static final Font SEARCH_FONT = new Font("Times Roman",Font.ITALIC,SEARCH_FONT_SIZE);

    JAGCTVisualizationPane myVisualizationPane;
    View3D myView3D;

    AGCT myAGCT;
    Domain myDomain;
    String myName;

    int pointRadius = 2, pointSquare = 4, pMouseX, pMouseY, deltaX, deltaY, zoomingFactor;
    int xOffset, yOffset;
    Pnt3D reference_Pnt3D;
    //point around which are made all rotations
    String currentView;

    int hierarchyP, hierarchyF, hierarchyC;

    boolean plotAvailable, showVarIds = true, zooming, dragging;
    //plotAvailable = false when no plot is possible

    public static Color myColor = Color.black;
    public static Color highlightColor = Color.orange; //Color.red;
    public static Color paramColor = Color.blue;
    public static Color structColor = Color.lightGray; //Color.green;
    public static Color [] referencedColors = {Color.green,
					       Color.red,
					       Color.cyan,
					       Color.pink,
					       Color.orange,
					       Color.blue,
					       Color.magenta,
					       Color.darkGray,
					       darkGoldenRod,
					       darkSalmon};/*{aquaMarine, 
					       bisque, 
					       brown, 
					       chartreuse, 
					       coral, 
					       darkGoldenRod, 
					       darkOliveGreen, 
					       darkSalmon, 
					       midnightBlue, 
					       peru};*/

    public static Color bufferedGenesColor = peru;

    public static Color positiveCorrelationColor = Color.green;
    public static Color negativeCorrelationColor = Color.red;

    public static Color minColor = new Color(240, 240, 240);

    int index_data_pointed;

    JAGCTGraphicsPanel(JAGCTVisualizationPane myj, AGCT a, Domain d, String name){
	Pnt3D [] pt;

	myVisualizationPane = myj;

	currentView = JAGCTVisualizationPane.allViews[0]; //Change if not XYZ !

	plotAvailable = false;
	myAGCT = a;
	myDomain = d;
	myName = name;
	
	pt = JAGCTVisualizationPane.getViewPoint(currentView);
	if (pt[1] == null)
	    myView3D = new View3D(pt[0], this);
	else
	    myView3D = new View3D(pt[0], pt[1], this);

	//myView3D.setPerspective(false);
	reference_Pnt3D = Pnt3D.o;

	zooming = false;
	dragging = false;

	index_data_pointed = -1;
	hierarchyP = hierarchyF = hierarchyC = JAGCTVisualizationPane.DEFAULT_HIERARCHY_LEVEL;
    }

    public void modifyHierarchy(int v){
	hierarchyP = hierarchyF = hierarchyC = v;
    }

    public void hiPlus(){
	if ( (hierarchyP != hierarchyF) || (hierarchyF != hierarchyC) )
	    Matrix.perror("JAGCTGraphicsPanel.class :: different hierarchies values");
	if (hierarchyF < myDomain.maxHierarchy-1){
	    hierarchyF++;
	    hierarchyP++;
	    hierarchyC++;
	}
    }

    public void hiMinus(){
	if ( (hierarchyP != hierarchyF) || (hierarchyF != hierarchyC) )
	    Matrix.perror("JAGCTGraphicsPanel.class :: different hierarchies values");
	if (hierarchyF > -1){
	    hierarchyF--;
	    hierarchyP--;
	    hierarchyC--;
	}
    }

    public Rectangle getCaptureRectangle(){
	Rectangle bounds = getBounds();
        bounds.setLocation(getLocationOnScreen());
	return bounds;
    }

    public void initAll(boolean zoom, boolean xO, boolean yO){
	if (plotAvailable == true){

	}
    }

    public void changeView(int i){
	currentView = JAGCTVisualizationPane.allViews[i];
	Pnt3D [] pt = JAGCTVisualizationPane.getViewPoint(currentView);
	if (pt[1] == null)
	    myView3D.set(new View3DInfo(pt[0], false));
	else
	    myView3D.set(new View3DInfo(pt[0], pt[1], false));

	//myView3D.set(new View3DInfo(JAGCTVisualizationPane.getViewPoint(currentView), false));

	repaint();
    }

    public void setDomain(Domain d){
	myDomain = d;
    }

    public void paintComponent (Graphics g) {
	myView3D.setGraphics(g, getWidth(), getHeight());
	if (!plotAvailable){
	    g.drawString("No visualisation possible: no data processed", 20, 20);
	}else{
	    super.paintComponent(g);
	    
	    if (zooming){
		myView3D.zoom(Math.pow(1.01,-zoomingFactor));
		zooming = false;
	    }
	    
	    if (dragging){
		myView3D.pan(deltaX, deltaY);
		myView3D.adjustCamera(JAGCTGraphicsPanel.SENSITIVITY*deltaX, JAGCTGraphicsPanel.SENSITIVITY*deltaY);

		if ( (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) && (AGCT.Sort_Depth) )
		    myVisualizationPane.updateDepth(myView3D);
		dragging = false;
	    }
	    
	    myView3DUpdateDepth();

	    if (!(myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) && (myAGCT.myClusteringProfileFrame.atLeastTwoGenes()) )
		drawListGenes(g);

	    plotPointsAndEdges(g);
	    if (myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
		drawEdgesFeatures(g);

	    drawStructure(g);
	    if ( !(myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) && (myVisualizationPane.clusterStructure.isSelected()) )
		myVisualizationPane.plotClustering(g);
	}
    }

    public void drawListGenes(Graphics g){
	Vector lg = myAGCT.myClusteringProfileFrame.listGenes;
	Domain d = myAGCT.myDomain;
	int i, j, num1, num2;
	Color old;
	Pnt3D p1, p2;
	old = g.getColor();

	Graphics2D g2d = (Graphics2D)g;
	g2d.setStroke(new BasicStroke(2));

	if (!myAGCT.Use_Shadow)
	    g.setColor(JAGCTGraphicsPanel.bufferedGenesColor);
	for (i=0;i<lg.size()-1;i++){
	    num1 = ((Integer) lg.elementAt(i)).intValue();
	    num2 = ((Integer) lg.elementAt(i+1)).intValue();
	    p1 = myVisualizationPane.getRightPoint(num1); 
	    p2 = myVisualizationPane.getRightPoint(num2); 
	    if (AGCT.Use_Shadow)
		myView3D.drawShadowLineVsRef(p1, p2, JAGCTGraphicsPanel.minColor, JAGCTGraphicsPanel.bufferedGenesColor, true);
	    else
		myView3D.drawLineVsRef(p1, p2, -1, false);
	}
	g.setColor(old);

	g2d.setStroke(new BasicStroke());
    }

    public void myView3DUpdateDepth(){
	if (myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)){
	    myView3D.minDepth = -1.0;
	    myView3D.maxDepth = 1.0;
	}else{
	    int i;
	    myView3D.initDepth();
	    for (i=0;i<myVisualizationPane.getTotalFeatures();i++){
		myView3D.updateDepth(myVisualizationPane.getRightPoint(i));
	    }
	}
    }

    public void drawEdgesFeatures(Graphics g){
	if (plotAvailable){
	    int i;
	    Pnt3D pt;
	    for (i=0;i<myAGCT.feature_pca_Pnt3D.size();i++){
		pt = (Pnt3D) myAGCT.feature_pca_Pnt3D.elementAt(i);
		if (i != index_data_pointed)
		    myView3D.drawShadowLineVsRef(Pnt3D.o, pt, JAGCTGraphicsPanel.minColor, Color.black, true);
		else
		    myView3D.drawShadowLineVsRef(Pnt3D.o, pt, JAGCTGraphicsPanel.minColor, highlightColor, true);
	    }

	    if ( (ControlProcess.hasTrue("softClusteringProcessed") || ControlProcess.hasTrue("hardClusteringProcessed")) && (myVisualizationPane.memberships.isSelected()) && (myVisualizationPane.currentClustering > -1) ){
		Pnt3D pcl = new Pnt3D();
		double norm;
		Clustering cc = myAGCT.getClustering(myVisualizationPane.currentClustering);
		
		for (i=0;i<cc.nclusters;i++){
		    pcl.coordinates[0] = cc.VV.coordinates[i][myVisualizationPane.xAxis];
		    pcl.coordinates[1] = cc.VV.coordinates[i][myVisualizationPane.yAxis];
		    pcl.coordinates[2] = cc.VV.coordinates[i][myVisualizationPane.zAxis];
		    
		    norm = Math.sqrt(pcl.dot(pcl));
		    pcl.coordinates[0]/=norm;
		    pcl.coordinates[1]/=norm;
		    pcl.coordinates[2]/=norm;
		    
		    g.setColor(new Color(cc.redCluster[i], cc.greenCluster[i], cc.blueCluster[i]));
		    myView3D.drawLineVsRef(Pnt3D.o, pcl,-1,false);
		    
		    if (showVarIds)
			myView3D.drawStringVsRef("Cluster_" + i,pcl);
		}
	    }
	}
    }
    
    public void plotPointsAndEdges(Graphics g){
	int ii, i, id, j, rad, corel;
	Pnt3D p, p2;
	String nm;

	Gene gg;
	Color c, prevc, nextc, cap = null;
	String name;
	Vector e;
	boolean bb, dum = false;
	double nangle;
	boolean okPlot, vetoRef, vetoTriangulated, vetoEdges;
	Vector element;

	if (plotAvailable == true){
	    for (ii=0;ii<myVisualizationPane.getTotalFeatures();ii++){
		id = myVisualizationPane.getRightIndex(ii);
		if ( (myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) 
		     || ( (JAnnotationFrame.prototypeIsDisplayed(myAGCT, id)) && (myVisualizationPane.genePlottedCluster(id)) ) ){
		    p = myVisualizationPane.getRightPoint(id); 

		    vetoRef = false;
		    vetoTriangulated = false;
		    vetoEdges = false;
		    if (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)){
			gg = ( (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[id]) );
			cap = null;
			if (gg.referenced)
			    cap = JAGCTGraphicsPanel.referencedColors[gg.typeReferenced];

			if ( (!myView3D.onlyReferenced) || (!gg.referenced) ) //( (myView3D.onlyReferenced) && (!gg.referenced) )
			    vetoRef = true;

			vetoEdges = ( (myView3D.onlyReferencedEdges) && (!Prototype.geneOrNeighborsIsReferencedPrototype(myAGCT,id)) );

			if ( (AGCT.Show_Only_Genes_With_Visible_Edge) && (myAGCT.isTriangulated) ){
			    vetoTriangulated = true;
			    if ( (gg.neighbors != null) && (!vetoEdges) ){
				for (j=0;j<gg.neighbors.size();j++){
				    p2 = myVisualizationPane.getRightPoint(((Integer) gg.neighbors.elementAt(j)).intValue());
				    nangle = ((Double) gg.neighbor_angles.elementAt(j)).doubleValue();
				    if ( (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0) && (AGCT.Show_Correlation != 0) ){
					vetoTriangulated = false;
					break;
				    }else if ( (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0))) && (AGCT.Show_Correlation != 1) ){
					vetoTriangulated = false;
					break;
				    }
				}
			    }
			}
		    }else{
			cap = Color.black;
			gg = null;
		    }
	
		    if ( (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) && (JAGCTGraphicsPanel.DISPLAY_K_DIMENSION_ESTIMATION == 1) ){
			gg = ( (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[id]) );
		    	plotDimension(g, gg, p);
		    }

		    if ( (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
			 && (myAGCT.allClusterings != null)
			 && (myAGCT.allClusterings.size() > 0)
			 && (myVisualizationPane.currentClustering > -1)
			 && ( (ControlProcess.hasTrue("softClusteringProcessed")) || (ControlProcess.hasTrue("hardClusteringProcessed")) ) ){
			gg = ( (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[id]) );

			plotMemberships(g, gg, p);
		    }
			
		    setRightFontAndColor(g, id);
		    if ( (!AGCT.Use_Shadow) || (id == index_data_pointed) ){
			if (id == index_data_pointed)
			    bb = true;
			else
			    bb = false;
			myView3D.drawPointVsRef(p, id, bb, vetoRef, vetoTriangulated, cap, gg);
		    }else
			myView3D.drawShadowPntVsRef(p, id, JAGCTGraphicsPanel.minColor, myColor, vetoRef, vetoTriangulated, cap, gg);//myView3D.drawShadowLineVsRef(p, p, JAGCTGraphicsPanel.minColor, myColor, true);
		    if (myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
			dum = true;
		    else if ( (id != index_data_pointed) || (!myVisualizationPane.stars.isSelected()) )
			dum = true;
		    else
			dum = false;
			
		    if ( (dum) && (!vetoTriangulated) && (AGCT.LOCAL_DIMENSION_COMPUTED) && (JAGCTGraphicsPanel.DISPLAY_K_DIMENSION_ESTIMATION == 2) )
			myVisualizationPane.drawDimensionVsRef(myView3D, id, this);

		    if ( (showVarIds) && (dum) && (!vetoTriangulated) )
			myVisualizationPane.drawLabelVsRef(myView3D, id, this, true);

		    if ( (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) && (!myAGCT.triangulationTimerIsTicking) ){
			gg = ( (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[id]) );
			    
			if ( (myDomain.myEdges != null) 
			     && (gg.neighbors != null) 
			     && ( ( (myVisualizationPane.informationManifold.equals(JAGCTVisualizationPane.allManifoldButtonString)) 
				    && (!vetoEdges) )
				  || ( (myVisualizationPane.informationManifold.equals(JAGCTVisualizationPane.onlyPointedButtonString))
				       && (id == index_data_pointed) ) ) ){

			    for (j=0;j<gg.neighbors.size();j++){
				    
				p2 = myVisualizationPane.getRightPoint(((Integer) gg.neighbors.elementAt(j)).intValue());
				if (id == index_data_pointed)
				    if (AGCT.Use_Shadow)
					myView3D.drawShadowLineVsRef(p, p2, JAGCTGraphicsPanel.minColor, highlightColor, true);
				    else{
					prevc = g.getColor();
					g.setColor(highlightColor);
					myView3D.drawLineVsRef(p, p2, -1, true);
					g.setColor(prevc);
				    }
				else{
				    okPlot = false;
					
				    nangle = ((Double) gg.neighbor_angles.elementAt(j)).doubleValue();
				    prevc = g.getColor();
				    nextc = null;
				    corel = 0;
				    if (nangle < Math.PI * Statistics.LIMIT_P_DELAUNAY/2.0){
					nextc = positiveCorrelationColor;
					okPlot = true;
					corel = 1;
				    }else if (nangle > Math.PI * (1.0 - (Statistics.LIMIT_P_DELAUNAY/2.0))){
					nextc = negativeCorrelationColor;
					okPlot = true;
					corel = -1;
				    }

				    
				    if ( (AGCT.Show_Correlation == 1) && (nextc == negativeCorrelationColor) )
					okPlot = false;
				    else if ( (AGCT.Show_Correlation == 0) && (nextc == positiveCorrelationColor) )
					okPlot = false;

				    if (okPlot){
					g.setColor(nextc);
					if (AGCT.Use_Shadow)
					    myView3D.drawShadowLineVsRef(p, p2, JAGCTGraphicsPanel.minColor, nextc, true);
					else
					    myView3D.drawLineVsRef(p, p2, -1, true);
				    }
				    g.setColor(prevc);
				    
				    /*if (AGCT.Use_Shadow)
				      myView3D.drawShadowLineVsRef(p, p2, JAGCTGraphicsPanel.minColor, myColor, true);
				      else
				      myView3D.drawLineVsRef(p, p2, -1, true);*/
				    
				}
			    }
			}
		    }
		}
	    }

	    if (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)){
		id = index_data_pointed;
		if ( (id != -1) && (JAnnotationFrame.prototypeIsDisplayed(myAGCT, id)) ){
		    gg = ( (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[id]) );
		    vetoRef = false;
		    if  ( (!myView3D.onlyReferenced) || (!gg.referenced) )// ( (myView3D.onlyReferenced) && (!gg.referenced) )
			vetoRef = true;

		    if ( (myDomain.myEdges != null) 
			 && (gg.neighbors != null) 
			 && (!myVisualizationPane.informationManifold.equals(JAGCTVisualizationPane.noPointedButtonString)) ){
			p = myVisualizationPane.getRightPoint(id); 
			for (j=0;j<gg.neighbors.size();j++){
			    p2 = myVisualizationPane.getRightPoint(((Integer) gg.neighbors.elementAt(j)).intValue());
			    if (AGCT.Use_Shadow)
				myView3D.drawShadowLineVsRef(p, p2, JAGCTGraphicsPanel.minColor, highlightColor, true);
			    else{
				prevc = g.getColor();
				g.setColor(highlightColor);
				myView3D.drawLineVsRef(p, p2, -1, true);
				g.setColor(prevc);
			    }
			}
		    }
		}
	    }

	    g.setFont(JAGCTGraphicsPanel.NORMAL_FONT);
	    g.setColor(myColor);
	}
    }
    
    public void plotDimension(Graphics g, Gene gg, Pnt3D p){
	double beginAngle, widthAngle, totalProp, curProp, cm;
	int j;

	if (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)){

	    beginAngle = 0.0;
	    widthAngle = ( 360.0 * gg.local_dimension[AGCT.K_Dimension_Estimation-1] / Domain.MAX_LOCAL_DIMENSION[AGCT.K_Dimension_Estimation-1] );
	    
	    g.setColor(Color.gray);

	    myView3D.fillShadowArcVsRef(p, JAGCTGraphicsPanel.minColor, Color.darkGray, (int) (1.5*myVisualizationPane.radiusMemberships), (int) beginAngle, (int) widthAngle);
	}
    }


    public void plotMemberships(Graphics g, Gene gg, Pnt3D p){
	double beginAngle, widthAngle, totalProp, curProp, cm;
	int j;
	Clustering cc = myAGCT.getClustering(myVisualizationPane.currentClustering);

	if ( (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
	     && (myVisualizationPane.memberships.isSelected()) ){
	    
	    if (myVisualizationPane.membershipViewSelect.getSelectedIndex() == 0){
		totalProp = 0.0;
		for (j=0;j<cc.nclusters;j++){
		    cm = gg.getClusterMemberships(myVisualizationPane.currentClustering, j);
		    totalProp += cm;
		}
		
		beginAngle = 0.0;
		for (j=0;j<cc.nclusters;j++){
		    curProp = gg.getClusterMemberships(myVisualizationPane.currentClustering, j); //myAGCT.myClustering.myClusteringAlgorithm.soft_memberships.coordinates[j][id];
		    widthAngle = ( 360.0 * curProp / totalProp );
		    if ( j == cc.nclusters )
			widthAngle = 360.0 - beginAngle;
		    
		    g.setColor(clusterColor(p, new Color(cc.redCluster[j], cc.greenCluster[j], cc.blueCluster[j])));

		    myView3D.fillArcVsRef(p, myVisualizationPane.radiusMemberships, (int) beginAngle, (int) widthAngle);
		    beginAngle = beginAngle + widthAngle;
		}
	    }else if (myVisualizationPane.membershipViewSelect.getSelectedIndex() == 1){
		for (j=0;j<cc.nclusters;j++){
		    g.setColor(clusterColor(p, new Color(cc.redCluster[j], cc.greenCluster[j], cc.blueCluster[j])));
		    cm = gg.getClusterMemberships(myVisualizationPane.currentClustering, j);
		    myView3D.drawCircleVsRef(p, (int) ( ( (double) myVisualizationPane.radiusMemberships ) * cm) );
		}
	    }
	}
    }

    public Color clusterColor(Pnt3D p, Color c){

	if (!AGCT.Use_Shadow){
	    return c;
	}else{
	    double prop = 1.0 - ( (myView3D.clampDepth(p)  - myView3D.minDepth ) / (myView3D.maxDepth - myView3D.minDepth) );
	    if (prop < 0.0)
		prop = 0.0;
	    if (prop > 1.0)
		prop = 1.0;

	    double dr = (255.0 - (double) c.getRed()) * prop;
	    double dg = (255.0 - (double) c.getGreen()) * prop;
	    double db = (255.0 - (double) c.getBlue()) * prop;

	    int vr, vg, vb;

	    if (c.getRed() >= JAGCTVisualizationPane.Max_Shadow_Color)
		vr = c.getRed();
	    else{
		vr = (int) (dr + (double) c.getRed());
		if (vr > JAGCTVisualizationPane.Max_Shadow_Color)
		    vr = JAGCTVisualizationPane.Max_Shadow_Color;
	    }

	    if (c.getGreen() >= JAGCTVisualizationPane.Max_Shadow_Color)
		vg = c.getGreen();
	    else{
		vg = (int) (dg + (double) c.getGreen());
		if (vg > JAGCTVisualizationPane.Max_Shadow_Color)
		    vg = JAGCTVisualizationPane.Max_Shadow_Color;
	    }

	    if (c.getBlue() >= JAGCTVisualizationPane.Max_Shadow_Color)
		vb = c.getBlue();
	    else{
		vb = (int) (db + (double) c.getBlue());
		if (vb > JAGCTVisualizationPane.Max_Shadow_Color)
		    vb = JAGCTVisualizationPane.Max_Shadow_Color;
	    }

	    Color cp = new Color( vr, vg, vb );

	    return cp;
	}
    }

    public void setRightFontAndColor(Graphics g, int id){
	g.setFont(JAGCTGraphicsPanel.NORMAL_FONT);
	if ( (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) && ( ( (myVisualizationPane.geneF.isSelected()) && (myVisualizationPane.geneF.isEnabled()) ) ||  ( (myVisualizationPane.geneP.isSelected()) && (myVisualizationPane.geneP.isEnabled()) ) ) )
	    g.setFont(JAGCTGraphicsPanel.SMALL_FONT);
	g.setColor(myColor);

	if (id == myVisualizationPane.searchIndex){
	    g.setFont(JAGCTGraphicsPanel.SEARCH_FONT);
	    g.setColor(structColor);
	}
	
	if (id == index_data_pointed){
	    g.setFont(JAGCTGraphicsPanel.HIGHLIGHT_FONT);
	    g.setColor(highlightColor);
	}
    }

    public void up(){
	yOffset += OFFSET_MIN;
	repaint();
    }
    
    public void down(){
	yOffset -= OFFSET_MIN;
	repaint();
    }
    
    public void left(){
	xOffset += OFFSET_MIN;
	repaint();
    }
    
    public void right(){
	xOffset -= OFFSET_MIN;
	repaint();
    }
    
    public void drawStructure(Graphics g){
	g.setColor(structColor);

	if (!AGCT.Perspective){
	    myView3D.drawStringVsRef("x",Pnt3D.i);
	    myView3D.drawStringVsRef("y",Pnt3D.j);
	    myView3D.drawStringVsRef("z",Pnt3D.k);
	    
	    if (AGCT.Use_Shadow){
		myView3D.drawShadowLineVsRef(Pnt3D.o,Pnt3D.i, JAGCTGraphicsPanel.minColor, structColor, false);
		myView3D.drawShadowLineVsRef(Pnt3D.o,Pnt3D.j, JAGCTGraphicsPanel.minColor, structColor, false);
		myView3D.drawShadowLineVsRef(Pnt3D.o,Pnt3D.k, JAGCTGraphicsPanel.minColor, structColor, false);
	    }else{
		myView3D.drawLineVsRef(Pnt3D.o,Pnt3D.i, -1, true);
		myView3D.drawLineVsRef(Pnt3D.o,Pnt3D.j, -1, true);
		myView3D.drawLineVsRef(Pnt3D.o,Pnt3D.k, -1, true);
	    }
	    
	    g.setColor(structColor);
	    if (myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P))
		myView3D.drawCircleVsRef(Pnt3D.o, myView3D.v);
	}
    }

    public void mousePressed (MouseEvent e) {
	pMouseX = e.getX();
	pMouseY = e.getY();
	requestFocus();
    }

    public void mouseEntered (MouseEvent e) {
	pMouseX = e.getX();
	pMouseY = e.getY();
	dragging = true;
    }
    
    public void mouseExited (MouseEvent e) {}

    public void mouseClicked (MouseEvent e) {
	if (index_data_pointed >= 0){
	    xOffset = yOffset = 0;
	    reference_Pnt3D = myVisualizationPane.getRightPoint(index_data_pointed);
	    myVisualizationPane.searchIndex = index_data_pointed;
	    repaint();
	}
    }

    public void mouseReleased (MouseEvent e) {}

    public void keyPressed(KeyEvent e){
	if (e.getKeyCode() == KeyEvent.VK_T){
	    AGCT.Show_Only_Genes_With_Visible_Edge = !AGCT.Show_Only_Genes_With_Visible_Edge;
	}else
	    keyReleased(e);
    }

    public void keyReleased(KeyEvent e){
    	if (e.getComponent() != this) return;

	if (e.getKeyCode() == KeyEvent.VK_UP){
	    up();
	}else if (e.getKeyCode() == KeyEvent.VK_DOWN){
	    down();
	}else if (e.getKeyCode() == KeyEvent.VK_LEFT){
	    left();
	}else if (e.getKeyCode() == KeyEvent.VK_RIGHT){
	    right();
	}else if (e.getKeyCode() == KeyEvent.VK_ENTER){
	    xOffset = yOffset = 0;
	    reference_Pnt3D = Pnt3D.o;
	    myVisualizationPane.searchIndex = -1;
	}else if (e.getKeyCode() == KeyEvent.VK_F1){
	    adjustScales();
	}else if (e.getKeyCode() == KeyEvent.VK_Z){
	    zoomIn();
	}else if (e.getKeyCode() == KeyEvent.VK_A){
	    zoomOut();
	}else if (e.getKeyCode() == KeyEvent.VK_C){
	    AGCT.Show_Correlation ++;
	    if (AGCT.Show_Correlation >= 2)
		AGCT.Show_Correlation = -1;
	}else if (e.getKeyCode() == KeyEvent.VK_G){
	    if ( (!myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) && (index_data_pointed != -1) ){
		myAGCT.myClusteringProfileFrame.addGene(index_data_pointed);
		myAGCT.myClusteringProfileFrame.goPlotGenes();
	    }
	}else if (e.getKeyCode() == KeyEvent.VK_F){
	    myAGCT.myClusteringProfileFrame.flushGenes();
	    myAGCT.myClusteringProfileFrame.repaint();
	}
	repaint();
    }

    public void zoomIn(){
	zooming = true;
	zoomingFactor = SENSITIVITY;
	repaint();
    }

    public void zoomOut(){
	zooming = true;
	zoomingFactor = -SENSITIVITY;
	repaint();
    }

    public void adjustScales(){
	if ( (xOffset < - getWidth() / 2) || (xOffset > getWidth() / 2) )
	    xOffset = 0;
	if ( (yOffset < - getHeight() / 2) || (yOffset > getHeight() / 2) )
	    yOffset = 0;
	
	Thread t = new Thread(){
		public void run(){
		    boolean inside, zoomIn = true, test = true, stop = false;
		    Point pnt;
		    int i;
		    do{
			inside = true;
			i = 0;
			do{
			    pnt = myView3D.toPointVsRef(myVisualizationPane.getRightPoint(i));
			    if ( (pnt.x < 0) || (pnt.x > getWidth() ) || (pnt.y < 0) || (pnt.y > getHeight() ) )
				inside = false;
			    else
				i++;
			}while ( (inside == true) && ( i < myVisualizationPane.getTotalFeatures() ) );
			if (test == true){
			    if (inside)
				zoomIn = true;
			    else
				zoomIn = false;
			    test = false;
			}else{
			    if ( (zoomIn) && (inside) ){
				zooming = true;
				zoomingFactor = SENSITIVITY;
				repaint();
			    }else if ( (!zoomIn) && (!inside) ){
				zooming = true;
				zoomingFactor = - SENSITIVITY;
				repaint();
			    }else
				stop = true;
			}

			try { Thread.sleep(10); } catch (Exception e) {}
		    }while(!stop);
		}
	    };
	t.start();
    }

    public void keyTyped(KeyEvent e){}

    public void mouseDragged(MouseEvent e){
	myView3D.initDepth();
	

	int cMouseX, cMouseY;

	cMouseX = e.getX();
	cMouseY = e.getY();

	deltaX = (cMouseX - pMouseX);
	deltaY = (cMouseY - pMouseY);

	pMouseX = cMouseX;
	pMouseY = cMouseY;

	dragging = true;

	repaint();
    }

    public void mouseMoved(MouseEvent e){
	if (plotAvailable){
	    pMouseX = e.getX();
	    pMouseY = e.getY();

	    Point q = new Point(pMouseX, pMouseY);
	    int dcur, dmin = 0;
	    Point pnt;
	    int oi = -1, i, j, id, step = 0;
	    Gene gg;
	    String s1 = "";

	    for (i=0;i<myVisualizationPane.getTotalFeatures();i++){
		id = myVisualizationPane.getRightIndex(i);
		if ( (myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) 
		     || ( (JAnnotationFrame.prototypeIsDisplayed(myAGCT, id)) && (myVisualizationPane.genePlottedCluster(id)) ) ){

		    pnt = myView3D.toPointVsRef(myVisualizationPane.getRightPoint(id));
		    dcur = ( (pnt.x - (q.x + xOffset) ) * (pnt.x - (q.x + xOffset) ) ) + ( (pnt.y - (q.y + yOffset) ) * (pnt.y - (q.y + yOffset) ) );

		    if ( (step==0) || (dcur < dmin) ){
			oi = id;
			dmin = dcur;
		    }
		    step ++;
		}
	    }

	    if ( (oi != -1) && (index_data_pointed != oi) ){
		index_data_pointed = oi;
		s1 = "";
		
		if (myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)){
		    s1 = "Variable pointed: " + myVisualizationPane.getFeatureName(index_data_pointed);
		    myAGCT.myInformationFrame.setText(s1);
		}else if ( (myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.P_P)) || (myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.M_P)) ){
		    gg = (Gene) myDomain.domainGenes.elementAt(myDomain.selectedGeneNumberToGeneNumber[oi]);

		    myAGCT.displayHTMLGene(gg);
		}
		repaint();
	    }
	}
    }

    public void mouseWheelMoved(MouseWheelEvent e){
	zooming = true;
	int r = e.getWheelRotation();
	if (r < 0)
	    zoomingFactor = -SENSITIVITY;
	else if (r > 0)
	    zoomingFactor = SENSITIVITY;
	repaint();
    }

    public void componentHidden(ComponentEvent e){}

    public void componentMoved(ComponentEvent e){}

    public void componentResized(ComponentEvent e){
	initAll(true, true, true);
	myView3D.initDepth();
    }
    
    public void componentShown(ComponentEvent e){}
}
