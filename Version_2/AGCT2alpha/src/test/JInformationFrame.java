package test;

import forDebug.Debug;

import javax.swing.*;
import java.awt.*;

public class JInformationFrame extends JFrame {

    static String defaultProgressString = "(no computation)";

    JEditorPane informationPane;
    // Misc infos.

    MemoryMonitorBar mem;
    JProgressBar progressBar;
    JLabel labelProgressBar;
    String textLabelProgressBar;

    String myText;

    AGCT myAGCT;

    public void setTitle(String s) {
        informationPane.setBorder(BorderFactory.createTitledBorder(s));
    }

    public void setText(String s) {
        if (s == null) s = "";
        myText = s;
        try {
            informationPane.setText(myText);
        } catch (NullPointerException e) {
            Debug.debug("NullPointException");
            Debug.debug(e.getStackTrace());
        }
    }

    public void appendText(String s) {
        myText += s;
        informationPane.setText(myText);
    }

    public void setTextProgressBar(String s) {
        textLabelProgressBar = s;
        labelProgressBar.setText(textLabelProgressBar + ": ");
        repaint();
    }

    public void setProgressString(boolean b) {
        progressBar.setStringPainted(b);
    }

    public void setTextProgressString(String s) {
        textLabelProgressBar = s;
        progressBar.setString(textLabelProgressBar);
        repaint();
    }

    public void resetTextProgressString() {
        progressBar.setString(null);
        repaint();
    }

    public void setValueProgressBar(int v) {
        progressBar.setValue(v);
        repaint();
    }

    public void displayInfo() {
        progressBar = new JProgressBar(0, 100);
        progressBar.setValue(0);
        progressBar.setStringPainted(true);

        labelProgressBar = new JLabel(textLabelProgressBar);

        // Frank
        informationPane = new JEditorPane("text/html", "");// ();
        informationPane.setEditable(false);
//        informationPane.setBorder(BorderFactory.createTitledBorder("(misc. informations || modified Frank, October 2009)"));

        mem = new MemoryMonitorBar();

        JPanel progressPanel = new JPanel();
        progressPanel.setLayout(new BorderLayout());
        progressPanel.add(labelProgressBar, BorderLayout.WEST);
        progressPanel.add(progressBar, BorderLayout.CENTER);

        JPanel infraPane = new JPanel();
        infraPane.setLayout(new BorderLayout());
        infraPane.add(new JScrollPane(informationPane), BorderLayout.CENTER);
        infraPane.add(progressPanel, BorderLayout.SOUTH);

        Container pane = getContentPane();
        pane.setLayout(new BorderLayout());
        pane.add(mem, BorderLayout.EAST);
        pane.add(infraPane, BorderLayout.CENTER);
        // pane.add(new JScrollPane(informationPane), BorderLayout.CENTER);

//		addWindowListener(new FermetureListener("Closing AGCT's InformationFrame\n"));

        Image img = Toolkit.getDefaultToolkit().getImage(java.net.URLClassLoader.getSystemResource("Images/help.png"));
        setIconImage(img);
    }

    static public JInformationFrame getInstance() {
        return singleton;
    }

    static private JInformationFrame singleton = new JInformationFrame();

    private JInformationFrame() {
        super("AGCT --- Information Frame");
        setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        myText = "";
        textLabelProgressBar = JInformationFrame.defaultProgressString;

        displayInfo();
        setVisible(true);
        setResizable(true);
    }
}