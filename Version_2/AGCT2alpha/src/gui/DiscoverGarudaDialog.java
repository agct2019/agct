/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class DiscoverGarudaDialog extends JDialog {
    private Button yesButton, noButton;

    // コンストラクタ
    public DiscoverGarudaDialog(Frame host) {
        super(host, "Will you discover more software in garuda platform?", true);

        setBounds(128, 256, 300, 60);

        Container cp = getContentPane();
        cp.setLayout(new BoxLayout(cp, BoxLayout.LINE_AXIS));

//		Font f = new Font("ＭＳ 明朝", Font.PLAIN, 10);

        yesButton = new Button("Yes");
//		yesButton.setFont(f);
        cp.add(yesButton);

        noButton = new Button("No");
//		noButton.setFont(f);
        cp.add(noButton);

        DiscoverGarudaAction action = new DiscoverGarudaAction(this);
        yesButton.addActionListener(action);
        yesButton.setActionCommand("GARUDACALL");
        noButton.addActionListener(action);
        noButton.setActionCommand("EXIT");

        setVisible(true);

        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                setVisible(false);
            }
        });
    }

}
