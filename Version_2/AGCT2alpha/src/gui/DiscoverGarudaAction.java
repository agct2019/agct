/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gui;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * $Id: Dialog.html,v 1.1 2006/05/17 04:17:11 kishi Exp kishi $
 *
 * @author KISHI Yasuhiro
 */

public class DiscoverGarudaAction extends AbstractAction {

    private JDialog dialog;

    public DiscoverGarudaAction(JDialog dialog) {
        super();

        this.dialog = dialog;
    }

    public void actionPerformed(ActionEvent event) {

//        Component component = ( Component ) event.getSource();
        String command = event.getActionCommand();

//        System.out.println( command );

        if ("GARUDACALL".equals(command)) {

//            System.exit( 0 );

        } else {
            dialog.dispose();
        }

    }

}