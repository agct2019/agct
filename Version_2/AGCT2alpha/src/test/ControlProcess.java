package test;

import forDebug.Debug;

import java.util.Hashtable;
import java.util.Vector;

public class ControlProcess {

    /**
     * ****************************************************************************************
     * Class that can be helpful to control the current state of the program
     * ***
     */

    static Vector<Boolean> controls = new Vector<Boolean>();
    // Each element = a Boolean giving a state for the execution of AGCT
    static Hashtable<String, Integer> accessToBooleans = new Hashtable<String, Integer>();

    public static void removeAll() {
        controls = new Vector<Boolean>();
        accessToBooleans = new Hashtable<String, Integer>();
    }

    public static void put(String s, boolean b) {
        Debug.debug("ControlProcess : put", s, b);
        int i;
        if (accessToBooleans.containsKey(s)) {
            i = accessToBooleans.get(s);
            controls.setElementAt(b, i);
        } else {
            controls.addElement(b);
            accessToBooleans.put(s, controls.size() - 1);
        }
    }

    public static boolean get(String s) {
        int i = accessToBooleans.get(s);
        return controls.elementAt(i);
    }

    public static void assertContains(String s) {
        if (!accessToBooleans.containsKey(s))
            Matrix.perror("ControlProcess: " + s + " is not in the keys");
    }

    public static void assertTrue(String s) {
        assertContains(s);
        if (!get(s))
            Matrix.perror("ControlProcess: " + s + " is FALSE in the keys");
    }

    public static boolean hasTrue(String s) {
        return ((accessToBooleans.containsKey(s)) && (get(s)));
    }
}