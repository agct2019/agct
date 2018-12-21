package gui;

import javax.swing.*;

public class JProfileNormalizationMenu extends JMenu {
    public JRadioButtonMenuItem subtractMedian = new JRadioButtonMenuItem("subtract median");
    public JRadioButtonMenuItem kanoWay = new JRadioButtonMenuItem("Subtract circadian");

    public JProfileNormalizationMenu() {
        super("profile normalization");
        // TODO Auto-generated constructor stub
        // ここ、もしかして未完成？tk
        // add(subtractMedian);
        add(kanoWay);
//		kanoWay.doClick();
    }
}
