package test;

import javax.swing.*;
import java.awt.*;

class MemoryMonitorBar extends JPanel {
    public MemoryMonitorBar() {
        setLayout(new FlowLayout(FlowLayout.CENTER, 5, 0));

        this.label = new JLabel();
        add(label);

        Thread thread =
                new Thread() {
                    public void run() {
                        while (true) {
                            repaint();
                            try {
                                Thread.sleep(2000);
                            } catch (Exception e) {
                            }
                        }
                    }
                };
        thread.start();
    }

    public void paintComponent(Graphics g) {
        int rr, gg, pmemu, pmemm;

        super.paintComponent(g);

        long freeMemory = Runtime.getRuntime().freeMemory();
        long totalMemory = Runtime.getRuntime().totalMemory();
        long maxMemory = Runtime.getRuntime().maxMemory();

        long memoryUsed = totalMemory - freeMemory;

        pmemu = (int) (memoryUsed / (1024 * 1024));
        pmemm = (int) (maxMemory / (1024 * 1024));

        label.setText("<html><font size=\"1\">Mb:<br>" + pmemu + "<br>/<br>" + pmemm + "</P>");

        rr = (int) ((memoryUsed * 255) / maxMemory);
        gg = 255 - rr;

        Insets insets = getInsets();

        int width = getWidth() - insets.left - insets.right;
        int height = getHeight() - insets.top - insets.bottom;

        int x0 = insets.left;
        int x1 = insets.left + (int) (width * memoryUsed / maxMemory);
        int x2 = insets.left + (int) (width * totalMemory / maxMemory);
        int x3 = insets.left + width;

        int y0 = insets.top;
        int y1 = insets.top + height;

        int nx0 = insets.left;
        int nx1 = insets.left + width;

        int ny0 = insets.top + height;
        int ny3 = insets.top;
        int sz1 = (int) (height * memoryUsed / maxMemory);
        int ny1 = insets.top + height - sz1;
        int sz2 = (int) (height * totalMemory / maxMemory);
        int ny2 = insets.top + height - sz2;

        Graphics2D g2d = (Graphics2D) g;
        Paint oldPaint = g2d.getPaint();

        g2d.setColor(new Color(rr, gg, 0));
        g2d.fillRect(nx0, ny1, nx1, sz1);

        g2d.setColor(new Color(0x999999));
        g2d.fillRect(nx0, ny2, nx1, (sz2 - sz1));

        g2d.setColor(getBackground());
        g2d.fillRect(nx0, ny3, nx1, height - sz2);

    }

    private JLabel label;

}
