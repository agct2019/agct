import java.lang.*;
import java.util.*;

class ControlProcess {

    /*******************************************************************************************
     * Class that can be helpful to control the current state of the program
     *****/

    static Vector controls = new Vector();
    //Each element = a Boolean giving a state for the execution of AGCT
    static Hashtable accessToBooleans = new Hashtable();
    //Each key = a String, returns an Integer corresponding to the Boolean of this String control

    public static void removeAll(){
	controls = new Vector();
	accessToBooleans = new Hashtable();
    }

    public static void put(String s, boolean b){
	int i;
	if (accessToBooleans.containsKey(s)){
	    i = ( (Integer) accessToBooleans.get(s) ).intValue();
	    controls.setElementAt(new Boolean(b),i);
	}else{
	    controls.addElement(new Boolean(b));
	    accessToBooleans.put(s,new Integer(controls.size()-1));
	}
    }

    public static boolean get(String s){
	int i = ( (Integer) accessToBooleans.get(s)).intValue();
	return ( (Boolean) controls.elementAt(i)).booleanValue();
    }

    public static void assertContains(String s){
	if (!accessToBooleans.containsKey(s))
	    Matrix.perror("ControlProcess: " + s + " is not in the keys");
    }

    public static void assertTrue(String s){
	assertContains(s);
	if (!get(s))
	    Matrix.perror("ControlProcess: " + s + " is FALSE in the keys");
    }

    public static boolean hasTrue(String s){
	return ( (accessToBooleans.containsKey(s)) && (get(s)) );
    }
}