/** This class represents the device independant part of  a 3D view.
@author Tim Lambert & Richard Nock
 */
import java.util.*;

 class View3DInfo{
  protected double wx,wy,wwidth,wheight; //the window on UV plane
  protected double dinverse;  // 1/distance of centre of projection from UV plane
  public Pnt3D u; // These three vectors are the 
  public Pnt3D v; // basis of the UVW coordinate
  public Pnt3D w; // system
    public boolean onlyReferenced, onlyReferencedEdges;
 
  public View3DInfo() {
    this(true);
  }

  public View3DInfo(boolean persp) {
      setWindow(-1,-1,2,2);  //setWindow(-1,-1,2,2);
    setPerspective(persp);
  }
	
  public View3DInfo(Pnt3D dirn, Pnt3D up) {
    this();
    setCamera(dirn,up);
  }

  public View3DInfo(Pnt3D dirn, Pnt3D up, double dinverse) {
    this(dirn,up);
    this.dinverse = dinverse;
  }
	
  public View3DInfo(Pnt3D dirn, Pnt3D up, boolean persp) {
    this(persp);
    setCamera(dirn,up);
  }

  public View3DInfo(Pnt3D dirn) {
    this(dirn,true);
  }
	
  public View3DInfo(Pnt3D dirn, boolean persp) {
    this(persp);
    setCamera(dirn);
  }
	
  public View3DInfo(Pnt3D dirn, Pnt3D up,double dinverse,
		    double wx,double wy,double wwidth,double wheight) {
    setWindow(wx,wy,wwidth,wheight);
    setCamera(dirn,up);
    this.dinverse = dinverse;
  }

    public void setOnlyReferenced(boolean b){
	onlyReferenced = b;
    }

    public void setOnlyReferencedEdges(boolean b){
	onlyReferencedEdges = b;
    }
    
  public void set(View3DInfo vi) {
    setWindow(vi.wx,vi.wy,vi.wwidth,vi.wheight);
    w = vi.w;
    v = vi.v;
    u = vi.u;
    dinverse = vi.dinverse;
  }

  /* is this projection orthographic*/
  public boolean orthographic(){
    return dinverse==0;
  }

  /* set projection to perspective (if arg is true) else orthographic */
  public void setPerspective(boolean persp){
    if (persp) {
      dinverse = 1/wwidth;
    } else {
      dinverse = 0;
    }
  }


  /** Calculate UVW coordinate system, given view direction, and up vector
   */
  public void setCamera(Pnt3D dirn, Pnt3D up){
    w = dirn.normalize();
    v = (up.subtract(w.scale(w.dot(up)))).normalize();
    u = v.cross(w);
  }
  
  /** Calculate UVW coordinate system, given view direction, using default up
    vector */
  public void setCamera(Pnt3D dirn){
    setCamera(dirn,new Pnt3D(0,1,0));
  }

  /** Adjust the position of the camera, given a displacement of the view
    position relative to U and V vectors adjustment is given in Device Independant coords
    */

  public void adjustCameraDI(double dx,double dy){
    setCamera(w.add(u.scale(dx).add(v.scale(dy))),v);
  }

  public void panDI(double dx,double dy){
    wx += dx*wwidth;
    wy += dy*wheight;
  }
    
  
  public void zoom(double scale){
    setWindow(wx-(scale-1)*wwidth/2,
	      wy-(scale-1)*wheight/2,
	      scale*wwidth,scale*wheight);
  }

  public void dolly(double scale){
    dinverse /= scale;
  }

  /** Set the window on the UV plane
   */
  public void setWindow(double wx,double wy,double wwidth,double wheight){
    this.wx = wx;
    this.wy = wy;
    this.wwidth = wwidth;
    this.wheight = wheight;
    if (!orthographic()){
      dinverse = 1/wwidth;
    }
  }

  /* Depth of a point in UVW system.
     Useful for doing a depth sort
     */
  public double depth(Pnt3D x){
    if (dinverse == 0){
      return w.dot(x);
    } else {
      return -w.scale(1/dinverse).subtract(x).length();
    }
  }

  public boolean frontFace(Pnt3D normal, Pnt3D p) {
    if (dinverse == 0) {
      return normal.dot(w)>0;
    } else {
      return normal.dot(w.scale(1/dinverse).subtract(p))>0;
    }
  }

  public String toString(){
    return(w+" "+v+" "+(float)dinverse)+" "+(float)wx+" "+(float)wy+" "+
      (float)wwidth+" "+(float)wheight;
  }
}
