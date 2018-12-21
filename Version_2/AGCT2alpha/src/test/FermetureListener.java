package test;
//Contains all miscellaneous classes for AGCT

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class FermetureListener extends WindowAdapter {
    String myText;

    public FermetureListener(String s) {
        super();
        myText = s;
    }

    public void windowClosing(WindowEvent e) {
        System.out.println(myText);
        System.exit(0);
    }
}

