package gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class ClosingDialog extends JDialog {
    private Button yesButton, noButton;

    // コンストラクタ
    public ClosingDialog(Frame host) {
        super(host, "Will you stop running AGCT?", true);

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

        ClosingDialogAction action = new ClosingDialogAction(this);
        yesButton.addActionListener(action);
        yesButton.setActionCommand("EXIT");
        noButton.addActionListener(action);
        noButton.setActionCommand("NOT_EXIT");

        setVisible(true);

        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                setVisible(false);
            }
        });
    }

}