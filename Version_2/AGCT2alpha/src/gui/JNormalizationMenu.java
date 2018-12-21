package gui;

import javax.swing.*;
import java.awt.event.ActionListener;

public class JNormalizationMenu extends JMenu {
    public JReplicateNormalizationMenu jReplicateNormalizationMenu = new JReplicateNormalizationMenu();
    public JProfileNormalizationMenu jProfileNormalizationMenu = new JProfileNormalizationMenu();

    public JNormalizationMenu() {
        super("Normalization");
        add(jReplicateNormalizationMenu);
        add(jProfileNormalizationMenu);
    }

    public void addSDThresholdChangeListener(ActionListener listener) {
        jReplicateNormalizationMenu.addSDChangeListener(listener);
    }
}
