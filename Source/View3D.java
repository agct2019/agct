/** This class implements a 3D view.
  It provides a mapping from points in XYZ space to pixels on the screen.
@author Richard Nock & Tim Lambert
 */

import java.awt.*;
import java.util.Vector;
import java.awt.font.*;
import javax.swing.JCheckBox;

 class View3D extends View3DInfo {
  public Graphics g;  //The graphics object to draw in
  private int width,height; //Width and height of g
  private double adjwx,adjwy;  //Adjusted window origin
  public double xscale,yscale; // scale factors in window -> viewport mapping
  public static Color bgcolor = Color.white; // background colour to use in clear() method
  protected Color ambient = new Color(50,50,50);
  protected Color lightColor = new Color(205,205,205);
  public Color[] colors = {Color.blue,Color.green,Color.red,Color.yellow,
			      Color.black,Color.black,Color.cyan,Color.magenta,
			      Color.cyan,Color.magenta};
  protected int defaultColor = 0;
  protected Pnt3D lightDirection = null; //direction of light source

    public JAGCTGraphicsPanel myPanel;
    double minDepth, maxDepth;

    public static Rectangle SCREEN_RECTANGLE = new Rectangle(Toolkit.getDefaultToolkit().getScreenSize());

  /** Creates a 3D view given a graphics object to draw in, a viewport
    within that graphics object, and a background colour
    */

    public View3D(Pnt3D dirn, JAGCTGraphicsPanel j){
    super(dirn);
    myPanel = j;
    minDepth = maxDepth = 0.0;
    setLightDirection(w);
    onlyReferenced = false;
    onlyReferencedEdges = false;
  }

    public View3D(Pnt3D dirn, Pnt3D up, JAGCTGraphicsPanel j){
	super(dirn, up);
    myPanel = j;
    minDepth = maxDepth = 0.0;
    setLightDirection(w);
    onlyReferenced = false;
    onlyReferencedEdges = false;
  }

  public View3D(Pnt3D dirn){
    super(dirn);
    setLightDirection(w);
    onlyReferenced = false;
    onlyReferencedEdges = false;
  }

  public void set(View3DInfo v){
    super.set(v);
    setLightDirection(w);
    onlyReferenced = false;
    onlyReferencedEdges = false;
  }

    public void initDepth(){
	minDepth = maxDepth = 0.0;
    }

  public void setGraphics(Graphics g,int width, int height){
    this.g = g;
    this.width = width;
    this.height = height;
    setWindow(wx,wy,wwidth,wheight);
    //g.clipRect(0,0,width,height);
  }    

  public void setWindow(double wx,double wy,double wwidth,double wheight){
    super.setWindow(wx,wy,wwidth,wheight);
    if (width==0) {
      return;
    }
    xscale = width/wwidth;     
    yscale = height/wheight;
    if (yscale > xscale) {
      yscale = xscale;     
      double newwheight = height/yscale;
      adjwy += (wheight - newwheight)/2.0;
      adjwx = wx;
    } else {
      xscale = yscale;     
      double newwidth = width/xscale;
      adjwx += (wwidth - newwidth)/2.0;
      adjwy = wy;
    }
  }
  
  public void adjustCamera(double deltax,double deltay){
    adjustCameraDI(-deltax/width,deltay/height);
    setLightDirection(w);
  }
 
  
  public void pan(int deltax, int deltay){
    panDI(-(double)deltax/(xscale*wwidth),(double)deltay/(yscale*wheight));
  }
  
  /** set direction of light source */
  public void setLightDirection(Pnt3D dirn){
    lightDirection = dirn.normalize();
  }

  
  /** set a colour for objects */
  public void setColor(int i,Color c){
    colors[i] = c;
  }

  /**get a colour index - translating -1 to default colour */
  public int getColorIndex(int i){
    return i==-1?defaultColor:i;
  }

  /** get a colour for objects */
  public Color getColor(int i){
    return colors[getColorIndex(i)];
  }


  /** set a default colour for objects */
  public int setDefaultColor(int i){
    int save = defaultColor;
    defaultColor = i;
    return save;
  }
  	
  /** clamp a colour component to 0..255 range */
  private static int clamp(double c){
    if (c >= 255) {
      return 255;
    } else if (c<=0) {
      return 0;
    } else {
      return (int) c;
    }
  }
      

  /** compute shading colour */
  public Color shade(int c1, int c2, Pnt3D n, Pnt3D p) {
    if (colors[getColorIndex(c1)] == null){
      return null;
    }
    double k;
    if (dinverse == 0) {
      k = n.dot(w);
    } else {
      k = n.dot(w.scale(1/dinverse).subtract(p).normalize());
    }
    int i;
    if (k>0) {
      i = c1;
    } else {
      i = c2;
      k = -k;
    }
    if (i<0) {
      i = defaultColor;
    }
    Color c = colors[i];
    Color lc = lightColor;
    Color a = ambient;
    return new Color(
		     clamp(a.getRed()+k*c.getRed()*lc.getRed()/255.0),
		     clamp(a.getGreen()+k*c.getGreen()*lc.getGreen()/255.0),
		     clamp(a.getBlue()+k*c.getBlue()*lc.getBlue()/255.0));
  }
  

    public void updateDepth(Pnt3D p){
	double d = depth(p);
	
	if (d < minDepth)
	    minDepth = d;
	
	if (d > maxDepth)
	    maxDepth = d;
    }

    public Point toPoint(Pnt3D p){

	double divisor = 1-dinverse*w.dot(p);
	return new Point((int)((u.dot(p)/divisor-wx)*xscale),
			 height-(int)((v.dot(p)/divisor-wy)*yscale));
    }

  /* Clear the drawing area
   */
  public void clear(){
    g.setColor(bgcolor);
    g.fillRect(0,0,width,height);
  }


    public int [] toInt(Color c){
	int [] v = new int [3];
	v[0] = c.getRed();
	v[1] = c.getGreen();
	v[2] = c.getBlue();
	return v;
    }
    
    public double clampDepth(Pnt3D p){
	double d = depth(p);

	if (d < minDepth)
	    d = minDepth;
	if (d > maxDepth)
	    d = maxDepth;

	return d;
    }

    /************************************************************************************
     * These methods DO NOT use reference points OR panel's offsets
     *****/


    //Inside the viewing rectangle ?

    public boolean inside(Point p){
	return ( (p.x >= 0) && (p.x <= width) && (p.y >= 0) && (p.y <= height) );
    }

    public void drawLine(Pnt3D p1, Pnt3D p2){
	Point s1 = toPoint(p1);
	Point s2 = toPoint(p2);
	if ( (inside(s1) || inside(s2)) )
	    g.drawLine(s1.x,s1.y,s2.x,s2.y);
    }
    
    public void drawString(String s, Pnt3D p1){
	Point s1 = toPointVsRef(p1);
	g.drawString(s,s1.x,s1.y);
    }
    
    public void drawStringBelow(String s, Pnt3D p1){
	Point s1 = toPoint(p1);
	FontMetrics fm = g.getFontMetrics();
	int w = fm.stringWidth(s);
	int h = fm.getAscent();
	g.drawString(s,s1.x-w/2,s1.y+h+2);
    }


    /************************************************************************************
     * These methods use reference points AND panel's offsets
     *****/


    public boolean insideVsRef(Point p){
	return ( (p.x - myPanel.xOffset>= 0) && (p.x - myPanel.xOffset <= width) && (p.y - myPanel.yOffset >= 0) && (p.y - myPanel.yOffset <= height) );
    }

    public void drawCircleVsRef(Pnt3D center, Pnt3D border){
	Point c = toPointVsRef(center);
	Point b = toPointVsRef(border);

	int dx = ( (c.x - b.x) * (c.x - b.x) );
	int dy = ( (c.y - b.y) * (c.y - b.y) );
	
	int rad = (int) Math.sqrt(dx + dy);
	
	int lx = c.x - rad;
	int ly = c.y - rad;

	g.drawArc( c.x - rad - myPanel.xOffset, c.y - rad - myPanel.yOffset, rad + rad, rad + rad, 0, 360);
    }

    public void drawCircleVsRef(Pnt3D center, int radius){
	Point c = toPointVsRef(center);
	if (insideVsRef(c)){
	    g.drawArc( c.x - radius - myPanel.xOffset, c.y - radius - myPanel.yOffset, radius + radius, radius + radius, 0, 360);
	}
    }
	
    public void fillArcVsRef(Pnt3D center, int radius, int beginAngle, int endAngle){
	Point c = toPointVsRef(center);
	if (insideVsRef(c)){
	    g.fillArc( c.x - radius - myPanel.xOffset, c.y - radius - myPanel.yOffset, radius + radius, radius + radius, beginAngle, endAngle);
	}
    }

    public void drawPointVsRef(Pnt3D pnt, int id, boolean bb, boolean vetoRef, boolean vetoTriangulated, Color cap, Gene gg){
	Point s = toPointVsRef(pnt);
	int rad = myPanel.pointRadius, tot = 0;
	float x, w;
	Graphics2D g2 = (Graphics2D) g;
	FontRenderContext frc = g2.getFontRenderContext();
	Font font = g2.getFont();
	String st;

	if ( (gg != null) && (myPanel.myAGCT.myAnnotationFrame != null) && (myPanel.myAGCT.myAnnotationFrame.showNumber.isSelected()) )
	    tot = JAnnotationFrame.numberOfDisplayedGenes(myPanel.myAGCT, id);

	if ( (bb) || (!vetoTriangulated) ){
	    if (insideVsRef(s)){	    

		if ( (tot>0) && (myPanel.myAGCT.myAnnotationFrame != null) && (!myPanel.myAGCT.myAnnotationFrame.onlyPrototypes.isSelected()) ){
		    st = "" + tot;
		    w = (float) font.getStringBounds(st, frc).getWidth();
		    x = s.x - myPanel.xOffset;
		    x -= w;
		    g2.drawString(st, x, s.y - myPanel.yOffset);
		}

		g.fillOval( s.x-rad - myPanel.xOffset, s.y-rad - myPanel.yOffset, rad + rad, rad + rad);
		if (!vetoRef)
		    referencedPlot(s, cap);
	    }
	    
	    if ( (Prototype.Prototypes_Selected) && (!Prototype.No_Reduction) && (!myPanel.myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P))  && (myPanel.myVisualizationPane.stars.isSelected()) ){
		funkyStar(s, id, bb);
	    }
	}
    }

    public void referencedPlot(Point s, Color cap){
	//draws a funny figure to mean that point is among the referenced
        int rad2 = 3 * myPanel.pointRadius;
	Color old = g.getColor();
	g.setColor(cap);
	g.fillOval(s.x-rad2 - myPanel.xOffset, s.y-rad2 - myPanel.yOffset, rad2 + rad2, rad2 + rad2);
	g.setColor(old);
	//g.draw3DRect(s.x-rad2 - myPanel.xOffset, s.y-rad2 - myPanel.yOffset, rad2 + rad2, rad2 + rad2, true);
    }

    public Color getShadowColor(Pnt3D pnt, Color miColor, Color maColor){
	Color c;
	int [] minColor = toInt(miColor), maxColor = toInt(maColor);
	double [] ratColor = new double[3];
	int k;
	double dd, prop;
	dd = clampDepth(pnt);
	prop = ( ( ( dd / 2.0 ) - minDepth ) / (maxDepth - minDepth) );
	
	for (k=0;k<3;k++){
	    if (maxDepth == minDepth)
		ratColor[k] = 0.5;
	    else
		ratColor[k] = prop;
	}
	c = new Color(minColor[0] + (int) ( (double) (maxColor[0] - minColor[0]) * ratColor[0] ),
		      minColor[1] + (int) ( (double) (maxColor[1] - minColor[1]) * ratColor[1] ),
		      minColor[2] + (int) ( (double) (maxColor[2] - minColor[2]) * ratColor[2] ));
	return c;
    }


    public void fillShadowArcVsRef(Pnt3D center, Color miColor, Color maColor, int radius, int beginAngle, int endAngle){
	Point c = toPointVsRef(center);
	Color ref = g.getColor();
	if (insideVsRef(c)){
	    g.setColor(getShadowColor(center, miColor, maColor));
	    g.fillArc( c.x - radius - myPanel.xOffset, c.y - radius - myPanel.yOffset, radius + radius, radius + radius, beginAngle, endAngle);
	    g.setColor(ref);
	}
    }

    public void drawShadowPntVsRef(Pnt3D pnt, int id, Color miColor, Color maColor, boolean vetoRef, boolean vetoTriangulated, Color cap, Gene gg){
	Point s = toPointVsRef(pnt);
	int rad = myPanel.pointRadius, tot = 0;
	float w, x;
	Graphics2D g2 = (Graphics2D) g;
	FontRenderContext frc = g2.getFontRenderContext();
	Font font = g2.getFont();
	String st;

	if ( (gg != null) && (myPanel.myAGCT.myAnnotationFrame != null) && (myPanel.myAGCT.myAnnotationFrame.showNumber.isSelected()) )
	    tot = JAnnotationFrame.numberOfDisplayedGenes(myPanel.myAGCT, id);

	if (!vetoTriangulated){
	    g.setColor(getShadowColor(pnt, miColor, maColor));
	    if (insideVsRef(s)){

		if ( (tot>0) && (myPanel.myAGCT.myAnnotationFrame != null) && (!myPanel.myAGCT.myAnnotationFrame.onlyPrototypes.isSelected()) ){
		    st = "" + tot;
		    w = (float) font.getStringBounds(st, frc).getWidth();
		    x = s.x - myPanel.xOffset;
		    x -= w;
		    g2.drawString(st, x, s.y - myPanel.yOffset);
		}

		g.fillOval( s.x-rad - myPanel.xOffset, s.y-rad - myPanel.yOffset, rad + rad, rad + rad);
		if (!vetoRef)
		    referencedPlot(s, cap);
	    }
	    
	    if ( (Prototype.Prototypes_Selected) && (!Prototype.No_Reduction) && (!myPanel.myVisualizationPane.myReferenceName.equals(JAGCTVisualizationPane.C_C_P)) && (myPanel.myVisualizationPane.stars.isSelected()) ){
		funkyStar(s, id, false);
	    }
	}
    }


    public void funkyStar(Point s, int id, boolean bb){
	int gid = myPanel.myAGCT.myDomain.selectedGeneNumberToGeneNumber[id], nid;
	int i, rad = myPanel.pointRadius / 2;
	double alpha = 0.0, ref, inc, factor, lmin = 0.0, lmax = 0.0, length = ((double) myPanel.myVisualizationPane.radiusMemberships) * JAGCTVisualizationPane.factorRadiusPrototypeCluster;
	Vector vid = (Vector) Prototype.Closest_Center_To_Cluster_Points.get(new Integer(gid));
	Vector did = (Vector) Prototype.Closest_Center_To_Normalized_Distortions.get(new Integer(gid));
	Color cur = g.getColor();
	Point p = new Point();
	boolean left = false, highl = false;;
	Color cap = null;

	ref = ( (double) myPanel.myVisualizationPane.sliderAlphaStars.getValue() / (double) JAGCTVisualizationPane.numberTicksAlphaStar ) * 2.0 * Math.PI;
	alpha = ref;

	if (bb){
	    if (myPanel.getWidth() < myPanel.getHeight())
		lmin = (double) myPanel.getWidth() / 8.0;
	    else
		lmin = (double) myPanel.getHeight() / 8.0;
	    lmax = lmin * 2.0;
	}

	if (rad == 0)
	    rad = 1;

	if (vid.size() >= 1){
	    g.setFont(JAGCTGraphicsPanel.SMALL_FONT);
	    inc = 2*Math.PI / (double) vid.size();
	    for (i=0;i<vid.size();i++){
		if (bb){
		    factor = ( (Double) did.elementAt(i) ).doubleValue();
		    length = lmin + ( factor * (lmax - lmin) );
		}
		nid = ( (Integer) vid.elementAt(i) ).intValue();
		if ( (!onlyReferenced) || (myPanel.myVisualizationPane.getDirectGene(nid).referenced) ){
		    p.x = s.x + (int) (length*Math.cos(alpha));
		    p.y = s.y + (int) (length*Math.sin(alpha));
		    
		    if ( (bb) ){ //(!myPanel.myVisualizationPane.memberships.isSelected()) &&
			g.setColor( myPanel.highlightColor );
			g.drawLine(s.x - myPanel.xOffset, s.y - myPanel.yOffset, p.x - myPanel.xOffset, p.y - myPanel.yOffset);
			g.setColor(cur);
		    }
		    
		    g.fillOval( p.x-rad - myPanel.xOffset, p.y-rad - myPanel.yOffset, rad + rad, rad + rad);
		    
		    if ( (onlyReferenced) && (myPanel.myVisualizationPane.getDirectGene(nid).referenced) && (i>0) ){
			cap = JAGCTGraphicsPanel.referencedColors[myPanel.myVisualizationPane.getDirectGene(nid).typeReferenced];
			referencedPlot(p, cap);
		    }
		    
		    if ( ( (myPanel.showVarIds) && (bb) ) || ( (onlyReferenced) && (myPanel.myVisualizationPane.getDirectGene(nid).referenced) && (i>0) ) ){
			if ( ( (alpha > Math.PI / 2.0) && (alpha < 3.0 * Math.PI / 2.0) )
			     || ( (alpha > (2.0 * Math.PI) + (Math.PI / 2.0) ) && (alpha < (2.0 * Math.PI) + (3.0 * Math.PI / 2.0) ) ) )
			    left = true;
			else
			    left = false;
			if (i==0)
			    highl = false;
			else
			    highl = true;
			drawGeneVsRef(myPanel.myVisualizationPane.getDirectGene(nid), p, myPanel, left, highl);
		    }
		}
		alpha += inc;
	    }
	}
    }

    public void drawCenterVsRef(Pnt3D pnt, int nclust, String sz){
	Point s = toPointVsRef(pnt);
	int rad = myPanel.pointSquare;

	if (insideVsRef(s)){	    
	    g.fillRect( s.x-rad - myPanel.xOffset, s.y-rad - myPanel.yOffset, rad + rad, rad + rad);
	    g.drawString(Clustering.Center_Name + nclust + "(" + sz + ")",s.x - myPanel.xOffset,s.y - myPanel.yOffset);
	}
    }

    public void drawCenterVsRefWithAnnotations(Pnt3D pnt, int nclust, String sz){
	Point s = toPointVsRef(pnt);
	int rad = myPanel.pointSquare, delta = 0;

	if (insideVsRef(s)){	    
	    g.fillRect( s.x-rad - myPanel.xOffset, s.y-rad - myPanel.yOffset, rad + rad, rad + rad);
	    g.drawString(Clustering.Center_Name + nclust + "(" + sz +")",s.x - myPanel.xOffset,s.y - myPanel.yOffset);
	    if (Statistics.allChi2Tests){
		delta += JAGCTGraphicsPanel.SMALL_FONT_SIZE;
		pileDisplayCenter(s, delta, nclust);
	    }
	}
    }

    public void pileDisplayCenter(Point s1, int deltaRef, int nc){
	int i, delta = deltaRef, t;
	float x, w;
	String refStringAnno = "", chiString = "";
	Graphics2D g2 = (Graphics2D) g;
	FontRenderContext frc = g2.getFontRenderContext();
	Font font = g2.getFont();
	Vector annoStatDum = null;
	
	g.setFont(JAGCTGraphicsPanel.SMALL_FONT);
	for (t=0;t<3;t++){
	    if ( (t==0) && ( (Statistics.annotationsF == null) || (Statistics.annotationsF.size() < 1) ) )
		continue;
	    if ( (t==1) && ( (Statistics.annotationsP == null) || (Statistics.annotationsP.size() < 1) ) )
		continue;
	    if ( (t==2) && ( (Statistics.annotationsC == null) || (Statistics.annotationsC.size() < 1) ) )
		continue;

	    if ( (t==0) && (myPanel.myDomain.someAnnotationsF) ){
		annoStatDum = (Vector) Statistics.annotationsF.elementAt(nc);
		refStringAnno = Domain.ANNOTATION_F;
	    }else if ( (t == 1) && (myPanel.myDomain.someAnnotationsP) ){
		annoStatDum = (Vector) Statistics.annotationsP.elementAt(nc);
		refStringAnno = Domain.ANNOTATION_P;
	    }else if ( (t == 2) && (myPanel.myDomain.someAnnotationsC) ){
		annoStatDum = (Vector) Statistics.annotationsC.elementAt(nc);
		refStringAnno = Domain.ANNOTATION_C;
	    }

	    if ( ( (t==0) && (myPanel.myDomain.someAnnotationsF) )
		 || ( (t == 1) && (myPanel.myDomain.someAnnotationsP) )
		 || ( (t == 2) && (myPanel.myDomain.someAnnotationsC) ) )
		continue;
	    else{
		chiString = refStringAnno;
		for (i=0;i<annoStatDum.size();i++)
		    chiString += (String) annoStatDum.elementAt(i) + "; ";
		
		w = (float) font.getStringBounds(chiString, frc).getWidth();
		x = s1.x - myPanel.xOffset;
		g2.drawString(chiString, x, s1.y - myPanel.yOffset + delta);	 
		delta += JAGCTGraphicsPanel.SMALL_FONT_SIZE;
	    }
	}
	g.setFont(JAGCTGraphicsPanel.NORMAL_FONT);
    }

    public void drawShadowCenterVsRef(Pnt3D pnt, int nclust, Color miColor, Color maColor){
	Point s = toPointVsRef(pnt);
	int rad = myPanel.pointSquare, k;
	double prop, dc = depth(pnt);
	Color c;
	int [] minColor = toInt(miColor), maxColor = toInt(maColor);
	double [] ratColor = new double[3];

	if (insideVsRef(s)){	   

	    prop = ( (dc - minDepth ) / (maxDepth - minDepth) );
	    if (prop < 0.0)
		prop = 0.0;
	    if (prop > 1.0)
		prop = 1.0;
	    for (k=0;k<3;k++){
		if (maxDepth == minDepth)
		    ratColor[k] = 0.5;
		else
		    ratColor[k] = prop;
	    }

	    c = new Color(minColor[0] + (int) ( (double) (maxColor[0] - minColor[0]) * ratColor[0] ),
			  minColor[1] + (int) ( (double) (maxColor[1] - minColor[1]) * ratColor[1] ),
			  minColor[2] + (int) ( (double) (maxColor[2] - minColor[2]) * ratColor[2] ));
 
	    g.setColor(c);

	    g.fillRect( s.x-rad - myPanel.xOffset, s.y-rad - myPanel.yOffset, rad + rad, rad + rad);
	    g.drawString(Clustering.Center_Name + nclust,s.x - myPanel.xOffset,s.y - myPanel.yOffset);
	}
    }


    public void drawLineVsRef(Pnt3D p1, Pnt3D p2, float stroke, boolean test){
	Point s1 = toPointVsRef(p1);
	Point s2 = toPointVsRef(p2);

	
	Graphics2D g2d = (Graphics2D)g;
	
	if (stroke != -1)
	    g2d.setStroke(new BasicStroke(stroke));
	

	if ( (!test) || (insideVsRef(s1)) || (insideVsRef(s2)) )
	    g.drawLine(s1.x - myPanel.xOffset, s1.y - myPanel.yOffset, s2.x - myPanel.xOffset, s2.y - myPanel.yOffset);

	if (stroke != -1)
	    g2d.setStroke(new BasicStroke());
    }

    public void drawStringVsRef(String s, Pnt3D p1){
	Point s1 = toPointVsRef(p1);
	
	g.drawString(s,s1.x - myPanel.xOffset,s1.y - myPanel.yOffset);
    }

    public void drawDimensionVsRef(Gene gg, Pnt3D p1, JAGCTGraphicsPanel pane){
	Point s1 = toPointVsRef(p1);
	
	int i, delta = 0;
	
	Graphics2D g2 = (Graphics2D) g;
	FontRenderContext frc = g2.getFontRenderContext();
	Font font = g2.getFont();
	LineMetrics metrics;
	String s, sdim;
	float x, y, w;
	Color old = null, cap = null;

	if ( (pane.myVisualizationPane.geneID.isEnabled()) && (gg.name != null) ){
		old = g2.getColor();

		g2.setColor(Color.red);
		sdim = "" + ((int) gg.local_dimension[AGCT.K_Dimension_Estimation-1]);
		g2.setFont(JAGCTGraphicsPanel.VERY_SMALL_FONT);
		g2.drawString(sdim, s1.x, s1.y - myPanel.yOffset - JAGCTGraphicsPanel.VERY_SMALL_FONT_SIZE);

		g2.setColor(old);
	}
	g2.setFont(font);
    }

    public void drawGeneVsRef(Gene gg, Point s1, JAGCTGraphicsPanel pane, boolean left, boolean highlightRef){
	// left = false -> plot at right
	// left = true -> plot at left

	int i, delta = 0;
	
	Graphics2D g2 = (Graphics2D) g;
	FontRenderContext frc = g2.getFontRenderContext();
	Font font = g2.getFont();
	LineMetrics metrics;
	String s, sdim;
	float x, y, w;
	Color old = null, cap = null;

	if ( (pane.myVisualizationPane.geneID.isSelected()) && (pane.myVisualizationPane.geneID.isEnabled()) && (gg.name != null) ){
	    s = gg.name;
	    w = (float) font.getStringBounds(s, frc).getWidth();
	    x = s1.x - myPanel.xOffset;
	    if (left == true)
		x -= w;

	    g2.drawString(s, x, s1.y - myPanel.yOffset);
	    delta += JAGCTGraphicsPanel.SMALL_FONT_SIZE;

	    // Displaying local dimensions if computed
	}	    
	
	g.setFont(JAGCTGraphicsPanel.SMALL_FONT);
	if ( ( ( (gg.referenced) && (onlyReferenced) ) || ( (pane.myVisualizationPane.geneTrueName.isSelected()) && (pane.myVisualizationPane.geneTrueName.isEnabled()) ) ) && (gg.asciiName != null) ){
	    s = gg.asciiName;
	    w = (float) font.getStringBounds(s, frc).getWidth();
	    x = s1.x - myPanel.xOffset;
	    if (left == true)
		x -= w;

	    if ( (gg.referenced) && (highlightRef) ){
		old = g2.getColor();
		cap = JAGCTGraphicsPanel.referencedColors[gg.typeReferenced];
		g2.setColor(cap);
	    }
	    g2.drawString(s, x, s1.y - myPanel.yOffset + delta);
	    if ( (gg.referenced) && (highlightRef) ){
		g2.setColor(old);
	    }
	    delta += JAGCTGraphicsPanel.SMALL_FONT_SIZE;
	}
	
	pileDisplayGeneAnnotations(gg, s1, g2, delta, pane, left);

	/*****
	      if (gg.annotations != null)
	    for (i=0;i<gg.annotations.size();i++){
		if ( (pane.myVisualizationPane.geneF.isSelected()) && (pane.myVisualizationPane.geneF.isEnabled()) && (((String[])gg.annotations.elementAt(i))[0].equals(Domain.ANNOTATION_F)) ){
		    s = ((String[])gg.annotations.elementAt(i))[0] + "  " + ((String[])gg.annotations.elementAt(i))[1];
		    w = (float) font.getStringBounds(s, frc).getWidth();
		    x = s1.x - myPanel.xOffset;
		    if (left == true)
			x -= w;

		    g2.drawString(s, x, s1.y - myPanel.yOffset + delta);	 
		    delta += JAGCTGraphicsPanel.SMALL_FONT_SIZE;
		}else if ( (pane.myVisualizationPane.geneP.isSelected()) && (pane.myVisualizationPane.geneP.isEnabled()) && (((String[])gg.annotations.elementAt(i))[0].equals(Domain.ANNOTATION_P)) ){
		    s = ((String[])gg.annotations.elementAt(i))[0] + "  " + ((String[])gg.annotations.elementAt(i))[1];
		    w = (float) font.getStringBounds(s, frc).getWidth();
		    x = s1.x - myPanel.xOffset;
		    if (left == true)
			x -= w;

		    g2.drawString(s, x, s1.y - myPanel.yOffset + delta);	 
		    delta += JAGCTGraphicsPanel.SMALL_FONT_SIZE;
		}
	    }
	*****/
    }

    public void pileDisplayGeneAnnotations(Gene gg, Point s1, Graphics2D g2, int deltaRef, JAGCTGraphicsPanel pane, boolean left){
	//Remark: we do not plot anymore the significant CHI2 annotations (too heavy)

	int i, delta = deltaRef, t;
	float x, w;
	String s = "", refStringAnno = "", chiString = "", dumString, dumdum;
	FontRenderContext frc = g2.getFontRenderContext();
	Font font = g2.getFont();
	Vector annoDum = null, annoStatDum = null;
	JCheckBox refBoxAnno = null;
	Color cref = g2.getColor();
	
	for (t=0;t<3;t++){
	    if (t==0){
		annoDum = gg.annotationsF;
		annoStatDum = Statistics.annotationsF;
		refStringAnno = Domain.ANNOTATION_F;
		refBoxAnno = pane.myVisualizationPane.geneF;
	    }else if (t == 1){
		annoDum = gg.annotationsP;
		annoStatDum = Statistics.annotationsP;
		refStringAnno = Domain.ANNOTATION_P;
		refBoxAnno = pane.myVisualizationPane.geneP;
	    }else if (t == 2){
		annoDum = gg.annotationsC;
		annoStatDum = Statistics.annotationsC;
		refStringAnno = Domain.ANNOTATION_C;
		refBoxAnno = pane.myVisualizationPane.geneC;
	    }


	    if ( (annoDum != null) && (annoDum.size() > 0) && (refBoxAnno.isSelected()) && (refBoxAnno.isEnabled()) ){
		if (!((String[])annoDum.elementAt(0))[0].equals(refStringAnno))
		    Matrix.perror("View3D.class :: annotation != that for " + refStringAnno + "  annotations");
		s = ((String[])annoDum.elementAt(0))[0] + " ";
		if (pane.hierarchyF >= 0){
		    if (annoDum.size() > pane.hierarchyF)
			i = pane.hierarchyF;
		    else
			i = annoDum.size() - 1;
		    if (!((String[])annoDum.elementAt(i))[0].equals(refStringAnno))
			Matrix.perror("View3D.class :: annotation != F in F annotations");
		    s += ((String[])annoDum.elementAt(i))[1];
		}else{
		    for (i=0;i<annoDum.size();i++){
			dumdum = ((String[])annoDum.elementAt(i))[1];
			if (Statistics.refNumClustering == pane.myVisualizationPane.currentClustering){
			    if (Statistics.containsInBH(dumdum, refStringAnno))
				s+=dumdum + " ";
			}else{
			    s += dumdum;
			    if (i<annoDum.size()-1)
				s += ", ";
			}
		    }
		}

		w = (float) font.getStringBounds(s, frc).getWidth();
		x = s1.x - myPanel.xOffset;
		if (left == true)
		    x -= w;
		
		g2.drawString(s, x, s1.y - myPanel.yOffset + delta);	 
		delta += JAGCTGraphicsPanel.SMALL_FONT_SIZE;
	    }

	    /*****if ( (annoDum != null) && (annoDum.size() > 0) && (Statistics.allChi2Tests) && (pane.myVisualizationPane.geneStat.isSelected()) && (pane.myVisualizationPane.geneStat.isEnabled()) ) {
		for (i=0;i<annoDum.size();i++){
		    dumString = ((String[])annoDum.elementAt(i))[1];
		    if (annoStatDum.contains(dumString)){
			chiString += refStringAnno;
			chiString += " ";
			chiString += dumString;
			chiString += "; ";
		    }
		}
		}*****/
	}

	/*****if ( (Statistics.allChi2Tests) && (pane.myVisualizationPane.geneStat.isSelected()) && (pane.myVisualizationPane.geneStat.isEnabled()) ){
	    g2.setColor(Statistics.CHI2_COLOR);
	    w = (float) font.getStringBounds(chiString, frc).getWidth();
	    x = s1.x - myPanel.xOffset;
	    if (left == true)
		x -= w;
	    
	    g2.drawString(chiString, x, s1.y - myPanel.yOffset + delta);	 
	    delta += JAGCTGraphicsPanel.SMALL_FONT_SIZE;
	    g2.setColor(cref);
	    }*****/
    }

    public void drawGeneVsRef(Gene gg, Pnt3D p1, JAGCTGraphicsPanel pane, boolean left, boolean highlightRef){
	Point s1 = toPointVsRef(p1);
	drawGeneVsRef(gg, s1, pane, left, highlightRef);
    }

    public Point toPointVsRef(Pnt3D p){
	updateDepth(p);

	Point dp = toPoint(p), dref = toPoint(myPanel.reference_Pnt3D);
	return new Point(dp.x - (dref.x - (width/2)), dp.y - (dref.y - (height/2)));
    }


    public void drawShadowLineVsRef(Pnt3D pnt1, Pnt3D pnt2, Color miColor, Color maColor, boolean clamp){
	Point s1 = toPointVsRef(pnt1);
	Point s2 = toPointVsRef(pnt2);

	if ( (insideVsRef(s1) || insideVsRef(s2)) ){
	    Pnt3D curPnt, nextPnt;
	    double [] ratColor = new double[3];
	    double sX, sY, sZ, ratX, ratY, ratZ, ratD, d1, d2, dc, ratS, prop, div;
	    Gene gg1, gg2;
	    int j, k;
	    float stroke;
	    int [] minColor = toInt(miColor), maxColor = toInt(maColor);
	    Color c;
	    
	    nextPnt = new Pnt3D(); 
	    curPnt = new Pnt3D(); 
	    curPnt.coordinates[0] = pnt1.coordinates[0];
	    curPnt.coordinates[1] = pnt1.coordinates[1];
	    curPnt.coordinates[2] = pnt1.coordinates[2];
	    
	    ratX = (pnt2.coordinates[0] - pnt1.coordinates[0]) / ( (double) JAGCTGraphicsPanel.SPLIT_POINTS );
	    ratY = (pnt2.coordinates[1] - pnt1.coordinates[1]) / ( (double) JAGCTGraphicsPanel.SPLIT_POINTS );
	    ratZ = (pnt2.coordinates[2] - pnt1.coordinates[2]) / ( (double) JAGCTGraphicsPanel.SPLIT_POINTS );
	    
	    if (clamp){
		d1 = clampDepth(pnt1);
		d2 = clampDepth(pnt2);
	    }else{
		d1 = depth(pnt1);
		d2 = depth(pnt2);
	    }

	    ratD = (d2 - d1) / (double) JAGCTGraphicsPanel.SPLIT_POINTS;
	    
	    dc = d1;

	    if (pnt1.equals(pnt2))
		div = 2.0;
	    else
		div = 1.0;
	    
	    for (j=0;j<JAGCTGraphicsPanel.SPLIT_POINTS;j++){
		prop = ( ( ( (dc + ratD) / 2.0 ) - minDepth ) / (maxDepth - minDepth) );
		
		for (k=0;k<3;k++){
		    if (maxDepth == minDepth)
			ratColor[k] = 0.5;
		    else
			ratColor[k] = prop;
		}
		
		ratS = (prop) * ((double) JAGCTGraphicsPanel.MAX_STROKE - (double) JAGCTGraphicsPanel.MIN_STROKE);
		stroke = JAGCTGraphicsPanel.MIN_STROKE + (float) ratS;
		
		dc += ratD;
		
		c = new Color(minColor[0] + (int) ( (double) (maxColor[0] - minColor[0]) * ratColor[0] ),
			      minColor[1] + (int) ( (double) (maxColor[1] - minColor[1]) * ratColor[1] ),
			      minColor[2] + (int) ( (double) (maxColor[2] - minColor[2]) * ratColor[2] ));
		
		nextPnt.coordinates[0] = curPnt.coordinates[0] + ratX;
		nextPnt.coordinates[1] = curPnt.coordinates[1] + ratY;
		nextPnt.coordinates[2] = curPnt.coordinates[2] + ratZ;
		
		g.setColor(c);
		
		//drawLineVsRef(curPnt, nextPnt, stroke); TOO SLOW !
		drawLineVsRef(curPnt, nextPnt, -1, false);
		
		curPnt.coordinates[0] = nextPnt.coordinates[0];
		curPnt.coordinates[1] = nextPnt.coordinates[1];
		curPnt.coordinates[2] = nextPnt.coordinates[2];
	    }
	}
    }
}
