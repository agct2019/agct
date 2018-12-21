package gene;

import javax.swing.*;
import javax.swing.border.LineBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;

public class JGeneSaveFrame extends JFrame {

    private final GeneListSaver setting;

    private class ExplainText extends JLabel {

        private static final String txt =
                "<html>"
                        + "<font size=4><b>"
                        + "&nbsp;Export Entity list<br>"
                        + "</b></font>"
                        + "<font size=3>"
                        + "&nbsp;&nbsp;Export the entity list with associated data and annotations as a plain text&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<br>"
                        + "&nbsp;&nbsp;file. Choose an interpretation to export normalized and/or raw signals<br>" + "&nbsp;&nbsp;from the experiment.<br>" + "</font>" + "</html>";

        private ExplainText() {
            super(txt);
            setBorder(new LineBorder(Color.BLACK, 1));
            setBackground(Color.WHITE);
            setOpaque(true);
        }
    }

    private class DataColumnSelectionPane extends JPanel {

        public DataColumnSelectionPane() {
            JCheckBox checkBox = new JCheckBox("Raw signal values");
            checkBox.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    setting.setRawSignalValueColumn(((JCheckBox) e.getSource()).isSelected());
                }
            });
            this.add(checkBox);
            checkBox = new JCheckBox("Normalized signal values");
            checkBox.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    setting.setNormalizedSignalValuesColumn(((JCheckBox) e.getSource()).isSelected());
                }
            });
            this.add(checkBox);
            checkBox = new JCheckBox("Flags");
            checkBox.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    setting.setFlagsColumn(((JCheckBox) e.getSource()).isSelected());
                }
            });
            this.add(checkBox);
            checkBox = new JCheckBox("Entitylist data");
            checkBox.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    setting.setEntitylistDataColumn(((JCheckBox) e.getSource()).isSelected());
                }
            });
            this.add(checkBox);

            setTitle("Select Data Columns");
            setLayout(new GridLayout(4, 1));
            setBorder(new TitledBorder(new LineBorder(Color.GRAY, 1), getTitle(), TitledBorder.LEFT, TitledBorder.TOP));
            setPreferredSize(new Dimension(480, 120));
        }
    }

    private class SelectAnnationColumns extends JPanel {

        private class Items extends JPanel implements GeneListSaverObserver {

            private final DefaultListModel listModel = new DefaultListModel();
            private final JList list = new JList(listModel);
            private final boolean condition;

            private Items(String name, boolean condition) {
                setLayout(new BorderLayout());
                add(new JLabel(name), BorderLayout.NORTH);
                this.condition = condition;
                final JScrollPane scr = new JScrollPane(list);
                add(scr);
                list.setVisibleRowCount(8);
                list.setFixedCellWidth(130);
                update();
                setting.addObserver(this);
            }

            public void update() {
                for (final String prop : Gene.ANNO_LIST) {
                    if (condition == setting.isSelected(prop)) {
                        if (!listModel.contains(prop)) {
                            listModel.addElement(prop);
                        }
                    } else {
                        if (listModel.contains(prop)) {
                            listModel.removeElement(prop);
                        }
                    }
                }
            }

            private String[] getSelectedStrings() {
                final Object[] os = list.getSelectedValues();
                final int n = os.length;
                final String[] res = new String[n];
                for (int i = 0; i < n; i++) {
                    res[i] = (String) os[i];
                }
                return res;
            }

            private void addListSelectionListener(ListSelectionListener listener) {
                list.addListSelectionListener(listener);
            }

            /**
             * selectedItemをup=trueなら上の、up=falseなら下のitemとswapする。
             * swapできなければ何もしない
             *
             * @param up
             */
            private void swap(boolean up) {
                assert list.getSelectedIndices().length == 1;
                final int i = list.getSelectedIndex();
                final int j = up ? i - 1 : i + 1;
                if (0 <= j && j < listModel.size()) {
                    final Object o = listModel.get(i);
                    listModel.set(i, listModel.get(j));
                    listModel.set(j, o);
                    list.setSelectedIndex(j);
                }
            }
        }

        private SelectAnnationColumns() {
            setTitle("Select Annotation Columns");
            final Items availableItems = new Items("Available items", false);
            final Items selectedItems = new Items("Selected items", true);
            final JButton right = new JButton("→");
            right.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent arg0) {
                    final String[] selected = availableItems.getSelectedStrings();
                    for (final String s : selected) {
                        setting.setSelected(s, true);
                    }
                }
            });
            final JButton left = new JButton("←");
            left.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent arg0) {
                    final String[] selected = selectedItems.getSelectedStrings();
                    for (final String s : selected) {
                        setting.setSelected(s, false);
                    }
                }
            });
            final JButton up = new JButton("↑");
            final JButton down = new JButton("↓");
            up.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    assert selectedItems.getSelectedStrings().length <= 1;
                    final String[] selected = selectedItems.getSelectedStrings();
                    if (selected.length != 1) {
                        return;
                    }
                    selectedItems.swap(true);
                }
            });
            down.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    assert selectedItems.getSelectedStrings().length <= 1;
                    final String[] selected = selectedItems.getSelectedStrings();
                    if (selected.length != 1) {
                        return;
                    }
                    selectedItems.swap(false);
                }
            });

            selectedItems.addListSelectionListener(new ListSelectionListener() {

                @Override
                public void valueChanged(ListSelectionEvent e) {
                    final JList jlist = (JList) e.getSource();
                    if (jlist.getSelectedValues().length != 1) {
                        up.setEnabled(false);
                        down.setEnabled(false);
                    } else {
                        up.setEnabled(true);
                        down.setEnabled(true);
                    }
                }
            });
            final JPanel rlPanel = new JPanel();
            rlPanel.setLayout(new BorderLayout());
            rlPanel.add(right, BorderLayout.NORTH);
            rlPanel.add(left, BorderLayout.SOUTH);
            final JPanel udPanel = new JPanel();
            udPanel.setLayout(new BorderLayout());
            udPanel.add(up, BorderLayout.NORTH);
            udPanel.add(down, BorderLayout.SOUTH);
            add(availableItems);
            add(rlPanel);
            add(selectedItems);
            add(udPanel);
            setPreferredSize(new Dimension(480, 200));
            setLayout(new FlowLayout());
            setBorder(new TitledBorder(getTitle()));
        }
    }

    private static int defaultCloseOperation = WindowConstants.DISPOSE_ON_CLOSE;

    public JGeneSaveFrame(GeneList geneList) {
        super("Gene Export");
        this.setting = new GeneListSaver(geneList);

        final JPanel mainPanel = new JPanel();
        final JPanel pane = (JPanel) getContentPane();
        mainPanel.add(new ExplainText());
        mainPanel.add(new DataColumnSelectionPane());
        mainPanel.add(new SelectAnnationColumns());
        mainPanel.setLayout(new FlowLayout());
        mainPanel.setPreferredSize(new Dimension(480, 600));

        final JPanel buttonPanel = new JPanel();
        final JButton go = new JButton("OK");
        go.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                final JFileChooser fileChooser = new JFileChooser(".");
                final int ret = fileChooser.showSaveDialog((Component) e.getSource());
                if (ret != JFileChooser.APPROVE_OPTION) {
                    return;
                }
                final File file = fileChooser.getSelectedFile();
                setting.setFile(file);
                final JProgressBar progressBar = new JProgressBar(0, 100);
                setting.addPropertyChangeListener(new PropertyChangeListener() {

                    @Override
                    public void propertyChange(PropertyChangeEvent evt) {
                        progressBar.setValue(setting.getProgress());
                    }
                });
//				try{
                setting.execute();
//				}
//				catch(final IOException e2){
//					JOptionPane.showMessageDialog((Component)e.getSource(),e2.toString());
//				}catch (Exception e2) {
//					e2.printStackTrace();
//				}
            }
        });

        final JButton watch = new JButton("Watch");
        watch.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {

            }
        });


        buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
        buttonPanel.add(go);


        pane.setLayout(new BorderLayout());
        pane.add(mainPanel, BorderLayout.NORTH);
        pane.add(buttonPanel, BorderLayout.SOUTH);

        // 何か追加したければここに入れればよい。centerが余っている。


        pane.setPreferredSize(new Dimension(480, 720));
        this.setDefaultCloseOperation(defaultCloseOperation);
        this.setSize(480, 720);
        //		this.setLayout(new BorderLayout());
        this.setResizable(false);
        this.setVisible(true);
    }

    public static void main(String[] args) {
        GeneList list = new GeneList();
        Gene piyo = new Gene("piyo");
        list.add(piyo);

        new JGeneSaveFrame(list);
    }
}