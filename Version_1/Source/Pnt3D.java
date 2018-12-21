import java.util.*;
 class Pnt3D extends Pnt{
  public static Pnt3D o = new Pnt3D(0,0,0);
  public static Pnt3D i = new Pnt3D(1,0,0);
  public static Pnt3D j = new Pnt3D(0,1,0);
  public static Pnt3D k = new Pnt3D(0,0,1);
  public static Pnt3D ijk = new Pnt3D(1,1,1);

  public Pnt3D(){
      super(3);
  }

  public Pnt3D(double x, double y, double z){
    this();
    this.coordinates[0] = x;
    this.coordinates[1] = y;
    this.coordinates[2] = z;
  }

  public static Pnt3D fromSpherical(double r,double theta,double phi){
    return new Pnt3D(r*Math.cos(theta)*Math.cos(phi),
		       r*Math.sin(theta)*Math.cos(phi),
		       r*Math.sin(phi));
  }
  
  public static Pnt3D fromCylindrical(double r,double theta,double y){
    return new Pnt3D(r*Math.cos(theta),
		       y,
		       r*Math.sin(theta));
  }
  
  public double x(){
    return coordinates[0];
  }

  public double y(){
    return coordinates[1];
  }

  public double z(){
    return coordinates[2];
  }

  public double theta(){
    return Math.atan2(coordinates[0],coordinates[2]);
  }

  public double r(){
    return Math.sqrt(coordinates[0]*coordinates[0]+coordinates[2]*coordinates[2]);
  }


  public String toString(){
    String s = "coordinates["+(float)coordinates[0]+","+(float)coordinates[1]+","+(float)coordinates[2]+"]";
    return s;
  }

  public static Pnt3D fromString(String s) throws NumberFormatException{
    StringTokenizer st = new StringTokenizer(s,"[,]");
    try{
      st.nextToken(); //get rid of leading coordinates
      double x = Double.valueOf(st.nextToken()).doubleValue();
      double y = Double.valueOf(st.nextToken()).doubleValue();
      double z = Double.valueOf(st.nextToken()).doubleValue();
      return new Pnt3D(x,y,z);
    } catch (NoSuchElementException e){
      throw new NumberFormatException();
    }
  }

  public Pnt3D add(Pnt3D x){
    Pnt3D a = new Pnt3D();
    for (int i=0;i<3;i++){
      a.coordinates[i] = coordinates[i] + x.coordinates[i];
    }
    return a;
  }

  public Pnt3D subtract(Pnt3D x){
    Pnt3D a = new Pnt3D();
    for (int i=0;i<3;i++){
      a.coordinates[i] = coordinates[i] - x.coordinates[i];
    }
    return a;
  }

  public Pnt3D scale(double x){
    Pnt3D a = new Pnt3D();
    for (int i=0;i<3;i++){
      a.coordinates[i] = x*coordinates[i];
    }
    return a;
  }

  public Pnt3D scale(double x, double y, double z){
    return new Pnt3D(x*coordinates[0],y*coordinates[1],z*coordinates[2]);
  }

  public Pnt3D normalize(){
    return scale(1/length());
  }

  public double length(){
    return Math.sqrt(dot(this));
  }

  public Pnt3D cross(Pnt3D x){
    return new Pnt3D(coordinates[1]*x.coordinates[2]-x.coordinates[1]*coordinates[2],
		       coordinates[2]*x.coordinates[0]-x.coordinates[2]*coordinates[0],
		       coordinates[0]*x.coordinates[1]-x.coordinates[0]*coordinates[1]);
  }
}
