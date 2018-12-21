package test;

import javax.swing.*;

class JBugFreeFrame extends JFrame {
    // This class solves Bug ID 6264013
    public JBugFreeFrame() {
        super();
        if (AGCT.MYDEBUG)
            AGCT.debug("JBugFreeFrame.JBugFreeFrame()");
        pack();
    }

    public JBugFreeFrame(String t) {
        super(t);
        if (AGCT.MYDEBUG)
            AGCT.debug("JBugFreeFrame.JBugFreeFrame(st)");
        pack();
    }

    public void dispose() {
        super.dispose();
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                LookAndFeel laf = UIManager.getLookAndFeel();
                if (laf != null) {
                    laf.uninitialize();
                }
            }
        });
    }
}
