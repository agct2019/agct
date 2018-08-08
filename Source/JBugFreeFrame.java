import javax.swing.*;

 class JBugFreeFrame extends JFrame {
    //This class solves Bug ID 6264013
    public JBugFreeFrame() {
	super();
	pack();
    }
    public JBugFreeFrame(String t) {
	super(t);
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

