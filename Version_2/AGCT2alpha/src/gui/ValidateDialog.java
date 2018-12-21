/**
 *
 */
package gui;

import javax.swing.*;
import java.awt.*;

/**
 * @author takashi
 */
public class ValidateDialog extends JFrame {

    /**
     * @throws HeadlessException
     */
    public ValidateDialog() throws HeadlessException {
        getContentPane().setLayout(new FlowLayout());

        JLabel label = new JLabel("Validation Finished");
        getContentPane().add(label);

        JButton button = new JButton("OK") {

        };


        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setTitle("Validated");
        setSize(200, 100);
        setVisible(true);
    }

    /**
     * @param gc
     */
    public ValidateDialog(GraphicsConfiguration gc) {
        super(gc);
    }

    /**
     * @param title
     * @throws HeadlessException
     */
    public ValidateDialog(String title) throws HeadlessException {
        super(title);
    }

    /**
     * @param title
     * @param gc
     */
    public ValidateDialog(String title, GraphicsConfiguration gc) {
        super(title, gc);
    }

}
