package filtering;

import forDebug.Debug;
import gene.Gene;
import gene.GeneList;
import test.AGCT;
import test.Domain;
import test.JAGCTGraphicsPanel;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;


public class JLeftIndicator extends JToolBar {
    public ArrayList<JCheckBox> getCheckBoxes() {
        System.out.println("JLeftIndicator.getCheckBoxes()");
        Debug.debug("components", getComponents());
        ArrayList<JCheckBox> res = new ArrayList<JCheckBox>();
        for (Component comp : this.getComponents()) {
            if (comp instanceof JCheckBox) {
                res.add((JCheckBox) comp);
            } else if (comp instanceof Box) {
                for (Component comp2 : ((Box) comp).getComponents()) {
                    if (comp2 instanceof JCheckBox) {
                        res.add((JCheckBox) comp2);
                    }
                }
            }
        }
        return res;
    }

    private GeneList geneList;
    //	private GeneList selectedGeneList

    private Box selectionBox;

    public static JLeftIndicator getInstance() {
        return singleton;
    }

    private boolean loaded = false;

    public boolean getLoaded() {
        return loaded;
    }

    private static JLeftIndicator singleton = new JLeftIndicator().init();

    private JLeftIndicator() {
        super(JToolBar.VERTICAL);
    }

    private JRadioButton jand, jor;

    private JLeftIndicator init() {
        selectionBox = Box.createVerticalBox();
        selectionBox.setVisible(true);
        selectionBox.setName("Filters");
        selectionBox.add(new JLabel("Filters"));
        add(selectionBox);
        selectionBox.add(Box.createHorizontalGlue());
        setFloatable(true);
        setVisible(true);

        jand = new JRadioButton("and");
        jor = new JRadioButton("or");
        ButtonGroup group = new ButtonGroup();
        group.add(jand);
        group.add(jor);
        jand.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent arg0) {
                Gene.setAndMode(true);
                rep();
            }
        });
        jor.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                Gene.setAndMode(false);
                rep();
            }
        });
        jand.setSelected(true);
        selectionBox.add(jand);
        selectionBox.add(jor);

        return this;
    }

    public void setGeneList(GeneList geneList) {
        this.geneList = geneList;
    }

    //	public void setSelectedGeneList(GeneList selectedGeneList) {
    //		this.selectedGeneList = selectedGeneList;
    //	}

    public void readResponsiveGeneFile(File[] files) {
        for (File file : files) {
            if (file.isDirectory()) {
                readResponsiveGeneFile(file.listFiles());
            } else {
                try {
                    BufferedReader br = new BufferedReader(new FileReader(file));
                    ArrayList<Gene> list = new ArrayList<Gene>();
                    HashSet<String> set = new HashSet<String>();
                    while (true) {
                        String s = br.readLine();
                        if (s == null)
                            break;
                        set.add(s.split("\t")[0]);
                        if (AGCT.MYDEBUG)
                            AGCT.debug(s.split("\t")[0]);
                    }
                    for (Gene gene : geneList) {
                        if (set.contains(gene.getName())) {
                            list.add(gene);
                        }
                    }
                    for (Gene gene : Domain.getInstance().getGenes()) {
                        if (set.contains(gene.getName())) {
                            gene.setAnyResponsive(true);
                        }
                    }
                    ResponsiveCheckBox checkBox = new ResponsiveCheckBox(file.getName().split("\\.")[0]);
                    checkBox.setResponsiveGenes(list);
                    checkBox.setVisible(true);
                    checkBox.doClick();
                    selectionBox.add(checkBox);
                    loaded = true;
                    geneList.get(0).myAGCT.myTabbedPane.stateChanged((ChangeEvent) null);
                } catch (IOException e) {
                    // TODO: handle exception
                    if (AGCT.MYDEBUG)
                        AGCT.debug("JResponsiveGeneSetter.readResponsiveGeneFile");
                    e.printStackTrace();
                    System.exit(0);
                }
            }
        }
    }

    private static int nowId = 0;

    private static String colorToString(Color color) {
        if (AGCT.MYDEBUG)
            AGCT.debug(color);
        return String.format("%02x%02x%02x", color.getRed(), color.getGreen(), color.getBlue());
    }

    private class ResponsiveCheckBox extends JCheckBox implements ActionListener {
        private ArrayList<Gene> responsiveGenes;
        private HashSet<Gene> set;
        private final int id;

        void setResponsiveGenes(ArrayList<Gene> responsiveGenes) {
            this.responsiveGenes = responsiveGenes;
            set = new HashSet<Gene>();
            for (Gene gene : this.responsiveGenes) {
                set.add(gene);
                gene.setResponsive(id, isSelected());
            }
            rep();
        }

        ResponsiveCheckBox(String caption) {
            super("<HTML><FONT SIZE=4 COLOR=" + colorToString(JAGCTGraphicsPanel.responsiveGeneColors[nowId]) + ">\u25CB</FONT>" + caption + "</HTML>");
            this.id = nowId++;
            addActionListener(this);
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            for (Gene gene : responsiveGenes) {
                gene.setResponsive(id, true);
            }
            for (Gene gene : Gene.getGeneList()) {
                if (!set.contains(gene)) {
                    if (this.isSelected()) {
                        gene.setVisible(false);
                    } else {
                        gene.setVisible(true);
                    }
                } else {
                    if (this.isSelected()) {
                        gene.setVisibleForOr(true);
                    } else {
                        gene.setVisibleForOr(false);
                    }
                }
            }
            rep();
        }
    }

    @Override
    public JLeftIndicator clone() {
        // TODO Auto-generated method stub
        // JResponsiveGeneSetter res = new JResponsiveGeneSetter();
        // res.selectionBox = this.selectionBox;
        // res.add(selectionBox);
        // res.setFloatable(true);
        // res.setVisible(true);
        // return res;
        throw new UnsupportedOperationException();
        // return singleton;
    }

    /**
     * このペーンが更新されたら呼ぶこと。
     */
    public void rep() {
        geneList.get(0).myAGCT.repaint();
    }

    private Box highlightIndicator = null;

    public void setHighlightIndicator(double[] times) {
        if (highlightIndicator == null) {
            highlightIndicator = Box.createVerticalBox();
            highlightIndicator.add(new JLabel("Peak time"));
            for (int i = 0; i < times.length; i++) {
                JHighlightCheckBox checkBox = new JHighlightCheckBox("<HTML><FONT SIZE=4 COLOR=" + colorToString(JAGCTGraphicsPanel.referencedColors[i])
                        + ">\u25CF</FONT> : time " + times[i] + "</HTML>");
                highlightIndicator.add(checkBox);
                checkBox.doClick();
            }
            JHighlightOtherCheckBox checkBox = new JHighlightOtherCheckBox("\u25CF : others");
            checkBox.doClick();
            highlightIndicator.add(checkBox);
            highlightIndicator.setVisible(true);
            this.add(highlightIndicator);
        }
    }

    static int nowId2 = 0;

    class JHighlightOtherCheckBox extends JCheckBox implements ActionListener {

        JHighlightOtherCheckBox(String caption) {
            super(caption);
            addActionListener(this);
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            // TODO Auto-generated method stub
            AGCT.highlightOther = isSelected();
            if (AGCT.MYDEBUG)
                AGCT.debug(AGCT.highlightOther);
            rep();
        }
    }

    class JHighlightCheckBox extends JCheckBox implements ActionListener {
        final int id;

        JHighlightCheckBox(String caption) {
            super(caption);
            id = nowId2++;
            addActionListener(this);
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            // TODO Auto-generated method stub
            if (isSelected()) {
                AGCT.HighlightMask |= 1 << id;
            } else {
                AGCT.HighlightMask &= ~(1 << id);
            }
            if (AGCT.MYDEBUG)
                AGCT.debug(Integer.toString(AGCT.HighlightMask, 2));
            rep();
        }
    }

    class JClusterCheckBox extends JCheckBox implements ActionListener {
        int[] clusterIds;

        JClusterCheckBox(int[] clusterIds) {
            super("cluster " + Arrays.toString(clusterIds));
            this.clusterIds = clusterIds;
            addActionListener(this);
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            int experimentId = AGCT.getInstance().getCurrentExperimentId();
            for (Gene gene : geneList) {
                int i = gene.getClusterNumber(experimentId);
                boolean found = false;
                for (int j : clusterIds)
                    if (i == j)// 表示すべきクラスタに含まれていたら
                        found = true;
                if (!found) {
                    gene.setVisible(!isSelected());
                }
            }
            rep();
        }
    }

    private Box clusterFilterIndicator = null;

    public void setCluster(int[] clusterIds) {
        System.out.println("JLeftIndicator.setCluster()");
        if (clusterFilterIndicator == null) {
            clusterFilterIndicator = Box.createVerticalBox();
            clusterFilterIndicator.add(new JLabel("clusterid"));// TODO clusteridが残っている
            clusterFilterIndicator.setVisible(true);
            this.add(clusterFilterIndicator);
        }

        JCheckBox checkBox = new JClusterCheckBox(clusterIds);
        checkBox.setVisible(true);
        checkBox.doClick();
        clusterFilterIndicator.add(checkBox);
        loaded = true;
        geneList.get(0).myAGCT.myTabbedPane.stateChanged((ChangeEvent) null);
        rep();
    }

    class JRangeCheckBox extends JCheckBox implements ActionListener {
        final int doseId;
        final double from, to;

        JRangeCheckBox(int doseId, double from, double to) {
            super(Domain.getInstance().getLigands().get(doseId).getSimpleName() + " range(" + from + "-" + to + ")");
            this.doseId = doseId;
            this.from = from;
            this.to = to;

            addActionListener(this);
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            if (isSelected()) {
                for (Gene gene : geneList) {
                    double peak = gene.peakInitialValue(doseId);
                    if (peak < from || to < peak) {
                        gene.setVisible(false);
                    }
                }
            } else {
                for (Gene gene : geneList) {
                    double peak = gene.peakInitialValue(doseId);
                    if (peak < from || to < peak) {
                        gene.setVisible(true);
                    }
                }
            }
            rep();
        }
    }

    private Box rangeFilterIndicator = null;

    public void setRange(int doseId, double from, double to) {
        System.out.println("JLeftIndicator.setRange()");
        if (rangeFilterIndicator == null) {
            rangeFilterIndicator = Box.createVerticalBox();
            rangeFilterIndicator.add(new JLabel("range of value"));
            rangeFilterIndicator.setVisible(true);
            this.add(rangeFilterIndicator);
        }

        JCheckBox checkBox = new JRangeCheckBox(doseId, from, to);
        checkBox.setVisible(true);
        checkBox.doClick();
        rangeFilterIndicator.add(checkBox);
        loaded = true;
        geneList.get(0).myAGCT.myTabbedPane.stateChanged((ChangeEvent) null);
        rep();
    }

    public int clusterMemberCount(int cid) {
        Debug.debug("BEGIN: JLeftIndicator # clusterMemberCount", cid);
        int res = 0;
        int experimentId = AGCT.getInstance().getCurrentExperimentId();
        int inCluster = 0;
        for (Gene gene : geneList) {
            if (gene.isVisible() && gene.isSelected()) {
                if (gene.getClusterNumber(experimentId) == cid)
                    res++;
            }
            if(gene.getClusterNumber(experimentId) == cid) inCluster++;
        }
        Debug.debug("END: JLeftIndicator # clusterMemberCount", cid, res, inCluster);
        return res;
    }

    class JPeakTimeCheckBox extends JCheckBox implements ActionListener {
        final int doseId;
        final int peakTimeId;

        public JPeakTimeCheckBox(int doseId, int peakTimeId) {
            super(Domain.getInstance().getLigands().get(doseId).getSimpleName() + " peak time(" + Domain.getInstance().getTimes()[peakTimeId] + ")");
            this.doseId = doseId;
            this.peakTimeId = peakTimeId;
            Debug.debug("JPeakTimeCheckBox", doseId, peakTimeId);
            addActionListener(this);
        }

        public void actionPerformed(ActionEvent e) {
            for (Gene gene : geneList) {
                int peak = gene.peakTimeId(doseId);
                Debug.debug(peak);
                if (peak != peakTimeId) {
                    if (this.isSelected()) {
                        gene.setVisible(false);
                    } else {
                        gene.setVisible(true);
                    }
                } else {
                    if (this.isSelected()) {
                        gene.setVisibleForOr(true);
                    } else {
                        gene.setVisibleForOr(false);
                    }
                }
            }
            rep();
        }
    }

    private Box peakTimeFilterIndicator = null;

    public void setPeakTime(int doseid, int peakTimeId) {
        System.out.println("JLeftIndicator.setPeakTime()");
        if (peakTimeFilterIndicator == null) {
            Box box = Box.createVerticalBox();
            box.add(new JLabel("peak time"));
            box.setVisible(true);
            this.peakTimeFilterIndicator = box;
            this.add(box);
        }

        JCheckBox checkBox = new JPeakTimeCheckBox(doseid, peakTimeId);
        checkBox.setVisible(true);
        checkBox.doClick();
        peakTimeFilterIndicator.add(checkBox);
        loaded = true;
        AGCT.getInstance().myTabbedPane.stateChanged((ChangeEvent) null);
        rep();
    }

    class JRangeOfCopyPerCellCheckBox extends JCheckBox implements ActionListener {
        final double minCopy;
        final double maxCopy;

        public JRangeOfCopyPerCellCheckBox(double minCopy, double maxCopy) {
            super(String.format("range(%.2f,%.2f)", minCopy, maxCopy));
            this.minCopy = minCopy;
            this.maxCopy = maxCopy;
            Debug.debug("JRangeOfCopyPerCellCheckBox", minCopy, maxCopy);
            addActionListener(this);
        }

        public void actionPerformed(ActionEvent e) {
            for (Gene gene : geneList) {
                double val = gene.calcAverageCopyPerCell();
                if (Double.isNaN(val))
                    continue;
                if (val < minCopy || maxCopy < val) {
                    if (this.isSelected()) {
                        gene.setVisible(false);
                    } else {
                        gene.setVisible(true);
                    }
                } else {
                    if (this.isSelected()) {
                        gene.setVisibleForOr(true);
                    } else {
                        gene.setVisibleForOr(false);
                    }
                }
            }
            rep();
        }
    }

    private Box rangeOfCopyPerCellFilterIndicator = null;

    public JCheckBox setRangeOfCopyPerCell(double mincopy, double maxcopy) {
        System.out.println("JLeftIndicator.setRangeOfCopyPerCell()");
        if (rangeOfCopyPerCellFilterIndicator == null) {
            Box box = Box.createVerticalBox();
            box.add(new JLabel("range of copy per cell"));
            box.setVisible(true);
            this.rangeOfCopyPerCellFilterIndicator = box;
            this.add(box);
        }
        JCheckBox checkBox = new JRangeOfCopyPerCellCheckBox(mincopy, maxcopy);
        checkBox.setVisible(true);
        checkBox.doClick();
        rangeOfCopyPerCellFilterIndicator.add(checkBox);
        loaded = true;
        AGCT.getInstance().myTabbedPane.stateChanged((ChangeEvent) null);
        rep();
        return checkBox;
    }

    class JRangeOfMagnitudeCheckBox extends JCheckBox implements ActionListener {
        final double minMagnitude;
        final double maxMagnitude;

        public JRangeOfMagnitudeCheckBox(double minMagnitude, double maxMagnitude) {
            super(String.format("magnitude(%.2f,%.2f)", minMagnitude, maxMagnitude));
            this.minMagnitude = minMagnitude;
            this.maxMagnitude = maxMagnitude;
            Debug.debug("JRangeOfMagnitudeCheckBox", minMagnitude, maxMagnitude);
            addActionListener(this);
        }

        public void actionPerformed(ActionEvent e) {
            for (Gene gene : geneList) {
                double magnitude = gene.calcTotalDistance();
                if (magnitude < minMagnitude || maxMagnitude < magnitude) {
                    if (this.isSelected()) {
                        gene.setVisible(false);
                    } else {
                        gene.setVisible(true);
                    }
                } else {
                    if (this.isSelected()) {
                        gene.setVisibleForOr(true);
                    } else {
                        gene.setVisibleForOr(false);
                    }
                }
            }
            rep();
        }
    }

    class JElimSurcCheckBox extends JCheckBox implements ActionListener {
        public JElimSurcCheckBox() {
            super("No series of copy per cell under 1");
            Debug.debug("JElimSurcCheckBox");
            addActionListener(this);
        }

        public void actionPerformed(ActionEvent e) {
            for (Gene gene : geneList) {
                //				double magnitude=gene.calcTotalDistance();
                if (!gene.doContainsSeriesOfAbove1()) {
                    if (this.isSelected()) {
                        gene.setVisible(false);
                    } else {
                        gene.setVisible(true);
                    }
                } else {
                    if (this.isSelected()) {
                        gene.setVisibleForOr(true);
                    } else {
                        gene.setVisibleForOr(false);
                    }
                }
            }
            rep();
        }
    }

    private Box rangeOfMagnitudeFilterIndicator = null;

    public JCheckBox setRangeOfMagnitude(double minMagnitude, double maxMagnitude) {
        System.out.println("JLeftIndicator.setRangeOfMagnitude()");
        if (rangeOfMagnitudeFilterIndicator == null) {
            Box box = Box.createVerticalBox();
            box.add(new JLabel("range of magnitude"));
            box.setVisible(true);
            this.rangeOfMagnitudeFilterIndicator = box;
            this.add(box);
        }
        JCheckBox checkBox = new JRangeOfMagnitudeCheckBox(minMagnitude, maxMagnitude);
        checkBox.setVisible(true);
        checkBox.doClick();
        rangeOfMagnitudeFilterIndicator.add(checkBox);
        loaded = true;
        AGCT.getInstance().myTabbedPane.stateChanged((ChangeEvent) null);
        rep();
        return checkBox;
    }

    class JRangeOfDelaunayCheckBox extends JCheckBox implements ActionListener {
        final double min;
        final double max;

        public JRangeOfDelaunayCheckBox(double min, double max) {
            super(String.format("delaunay(%.2f,%.2f)", min, max));
            this.min = min;
            this.max = max;
            Debug.debug("JRangeOfDelaunayCheckBox", min, max);
            addActionListener(this);
        }

        public void actionPerformed(ActionEvent e) {
            for (Gene gene : geneList) {
                double delaunay = gene.getMaxDelaunay();
                if (delaunay < min || max < delaunay) {
                    if (this.isSelected()) {
                        gene.setVisible(false);
                    } else {
                        gene.setVisible(true);
                    }
                } else {
                    if (this.isSelected()) {
                        gene.setVisibleForOr(true);
                    } else {
                        gene.setVisibleForOr(false);
                    }
                }
            }
            rep();
        }
    }

    private Box rangeOfDelaunayFilterIndicator = null;

    public JCheckBox setRangeOfDelaunay(double min, double max) {
        System.out.println("JLeftIndicator.setRangeOfDelaunay()");
        if (rangeOfDelaunayFilterIndicator == null) {
            Box box = Box.createVerticalBox();
            box.add(new JLabel("range of delaunay"));
            box.setVisible(true);
            this.rangeOfDelaunayFilterIndicator = box;
            this.add(box);
        }
        JCheckBox checkBox = new JRangeOfDelaunayCheckBox(min, max);
        checkBox.setVisible(true);
        checkBox.doClick();
        rangeOfDelaunayFilterIndicator.add(checkBox);
        loaded = true;
        AGCT.getInstance().myTabbedPane.stateChanged((ChangeEvent) null);
        rep();
        return checkBox;
    }

    private void setDelRange() {
        boolean ok = false;
        double min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;

        //		GeneList visible = new GeneList();

        for (Gene gene : Gene.getGeneList()) {
            if (gene.isVisible()) {
                ok = true;
                double val = gene.getMaxDelaunay();
                if (!Double.isNaN(val)) {
                    min = Math.min(min, val);
                    max = Math.max(max, val);
                }
            }
        }

        if (!ok) {
            return;
        }
        if (min == Integer.MAX_VALUE) {
            return;
        }
        JLeftIndicator.getInstance().setRangeOfDelaunay(min - 1E-3, max + 1E-3).doClick();
    }

    private Box elimSurcIndicator = null;

    private void setElimSurc() {
//        System.out.println("JLeftIndicator.elimSurc()");
//        if (elimSurcIndicator == null) {
//            Box box = Box.createVerticalBox();
//            box.add(new JLabel("elimSurc"));
//            box.setVisible(true);
//            this.elimSurcIndicator = box;
//            this.add(box);
//        }
//        JCheckBox checkBox = new JElimSurcCheckBox();
//        checkBox.setVisible(true);
//        checkBox.doClick();
//        elimSurcIndicator.add(checkBox);
//        loaded = true;
//        AGCT.getInstance().myTabbedPane.stateChanged((ChangeEvent) null);
//        rep();
//        checkBox.doClick();
//        return;
    }

    private void setMugRange() {
        boolean ok = false;
        double min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;
        for (Gene gene : Gene.getGeneList()) {
            if (gene.isVisible()) {
                ok = true;
                double val = gene.calcTotalDistance();
                min = Math.min(min, val);
                max = Math.max(max, val);
            }
        }
        if (!ok) {
            return;
        }
        JLeftIndicator.getInstance().setRangeOfMagnitude(min, max).doClick();
    }

    private void setCopyPerCellRange() {
        boolean ok = false;
        double min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;
        for (Gene gene : Gene.getGeneList()) {
            if (gene.isVisible()) {
                double val = gene.calcAverageCopyPerCell();
                Debug.debug("JLeftInducator : CopyPerCell val", val);
                if (Double.isNaN(val))
                    continue;
                ok = true;
                min = Math.min(min, val);
                max = Math.max(max, val);
            }
        }
        min = Math.max(min, 1);
        if (!ok) {
            //				JOptionPane.showMessageDialog(null,"No gene selected.");
            return;
        }
        JLeftIndicator.getInstance().setRangeOfCopyPerCell(min, max).doClick();
    }

    /**
     * とりあえずの実装。
     * すでにanchorのみが入っていることを前提とし、
     * それぞれのanchorを使用して、filter by range, filter by delauneyを施し、
     */
    public void doit() {
        System.out.println("JLeftIndicator.doit()");
        ArrayList<JCheckBox> boxes = getCheckBoxes();
        Debug.debug("boxes.size()", boxes.size());
        for (JCheckBox tmp : boxes) {
            if (!tmp.isSelected()) {
                tmp.doClick();
            }
        }
        // all selected

        jor.doClick();
        //		jand.doClick();
        //		for(JCheckBox box:boxes){
        //			box.doClick();
        setDelRange();
        setMugRange();
        setCopyPerCellRange();
        //			box.doClick();
        //		}
        setElimSurc();
        jand.doClick();
        //		jand.doClick();

        for (JCheckBox tmp : boxes) {
            if (tmp.isSelected()) {
                tmp.doClick();
            }
        }

        ArrayList<JCheckBox> after = getCheckBoxes();
        for (JCheckBox tmp : after) {
            if (!boxes.contains(tmp)) {
                tmp.doClick();
            }
        }
    }
}
