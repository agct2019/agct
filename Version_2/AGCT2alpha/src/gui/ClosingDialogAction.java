package gui;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * $Id: Dialog.html,v 1.1 2006/05/17 04:17:11 kishi Exp kishi $
 *
 * @author KISHI Yasuhiro
 */

public class ClosingDialogAction extends AbstractAction {

    private JDialog dialog;

    public ClosingDialogAction(JDialog dialog) {
        super();

        this.dialog = dialog;
    }

    public void actionPerformed(ActionEvent event) {

//        Component component = ( Component ) event.getSource();
        String command = event.getActionCommand();

//        System.out.println( command );

        if ("EXIT".equals(command)) {
            System.exit(0);

        } else {
            dialog.dispose();
        }

    }

}