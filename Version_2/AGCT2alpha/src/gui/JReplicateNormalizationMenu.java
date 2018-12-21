package gui;

import normalizationData.AbstractNormalizationData;
import normalizationData.NormalizationData;
import normalizationData.NormalizationData.NormalizeMode;
import test.Scenario;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;

public class JReplicateNormalizationMenu extends JMenu implements NormalizationDataObserver {
    static private String title = "replicate nomalization";
    private JMeanModeMenu jMeanMode = new JReplicateNormalizationMenu.JMeanModeMenu();
    private JMenuItem jSDthreshold = new JMenuItem("SDthreshold");
    private NormalizationData normalizationData;
    private List<ActionListener> sDChangeListeners = new ArrayList<ActionListener>();

    public void addSDChangeListener(ActionListener listener) {
        sDChangeListeners.add(listener);
    }

    public JReplicateNormalizationMenu() {
        super(title);
        add(jMeanMode);
        add(jSDthreshold);
        jSDthreshold.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent arg0) {
                requestModificationSDThreshold();
                for(ActionListener listener : sDChangeListeners) listener.actionPerformed(null);
            }
        });
    }

    public void setNormalizationData(NormalizationData normalizationData) {
        this.normalizationData = normalizationData;
    }

    void requestModificationSDThreshold() {
        String ret = JOptionPane.showInputDialog(this, "SD Threshold:", new Double(normalizationData.getSDThreshold()));
        if (ret != null) {
            double val = -1;
            try {
                val = Double.parseDouble(ret);//ここで変なものが入っているようです。
            } catch (NumberFormatException e) {
                System.err.println("Not an double!");
            }
            if (val >= 0) {
                normalizationData.setSDThreshold(val);
                Scenario.add("AGCT_Modification_SD_Threshold", val);
            }
        }
    }

    public void update(AbstractNormalizationData normalizationData) {
        // TODO Auto-generated method stub

    }

    class JMeanModeMenu extends JMenu {
        static final private String title = "Mean Mode";
        JRadioButtonMenuItem meanMode_ArithmeticMean = new JRadioButtonMenuItem("Arithmetic Mean");
        JRadioButtonMenuItem meanMode_GeometricMean = new JRadioButtonMenuItem("Geometric Mean");

        public JMeanModeMenu() {
            super(title);
            ButtonGroup groupMeanMode = new ButtonGroup();
            groupMeanMode.add(meanMode_ArithmeticMean);
            groupMeanMode.add(meanMode_GeometricMean);

            add(meanMode_ArithmeticMean);
            add(meanMode_GeometricMean);
            meanMode_ArithmeticMean.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent arg0) {
                    // TODO Auto-generated method stub
                    normalizationData.setNormalizeMode(NormalizeMode.Arithmetic_mean);
                }
            });
            meanMode_GeometricMean.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent arg0) {
                    // TODO Auto-generated method stub
                    normalizationData.setNormalizeMode(NormalizeMode.Geometric_mean);
                }
            });
            meanMode_GeometricMean.setSelected(true);
        }
    }
}
