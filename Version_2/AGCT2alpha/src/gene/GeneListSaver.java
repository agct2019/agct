package gene;

import javax.swing.*;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

// GUIとは無関係なクラス
public class GeneListSaver extends SwingWorker<Void, Gene> {
    private final ArrayList<GeneListSaverObserver> observers = new ArrayList<GeneListSaverObserver>();

    public void addObserver(GeneListSaverObserver observer) {
        observers.add(observer);
    }

    private void update() {
        for (GeneListSaverObserver observer : observers) {
            observer.update();
        }
    }

    private final GeneList geneList;
    private final Set<String> selectedAnnotations = new HashSet<String>();

    public GeneListSaver(GeneList geneList) {
        this.geneList = geneList;
        for (String s : Gene.ANNO_LIST) {
            this.selectedAnnotations.add(s);
        }
    }

    private boolean rawSignalValueColumn;
    private boolean normalizedSignalValuesColumn;
    private boolean flagsColumn;
    private boolean entitylistDataColumn;

    public boolean isRawSignalValueColumn() {
        return rawSignalValueColumn;
    }

    public void setRawSignalValueColumn(boolean rawSignalValueColumn) {
        System.out.println("GeneListFile.setRawSignalValueColumn()");
        System.out.println(rawSignalValueColumn);
        this.rawSignalValueColumn = rawSignalValueColumn;
    }

    public boolean isNormalizedSignalValuesColumn() {
        return normalizedSignalValuesColumn;
    }

    public void setNormalizedSignalValuesColumn(boolean normalizedSignalValuesColumn) {
        this.normalizedSignalValuesColumn = normalizedSignalValuesColumn;
    }

    public boolean isFlagsColumn() {
        return flagsColumn;
    }

    public void setFlagsColumn(boolean flagsColumn) {
        this.flagsColumn = flagsColumn;
    }

    public boolean isEntitylistDataColumn() {
        return entitylistDataColumn;
    }

    public void setEntitylistDataColumn(boolean entitylistDataColumn) {
        this.entitylistDataColumn = entitylistDataColumn;
    }

    boolean isSelected(String s) {
        assert Gene.ANNO_LIST.contains(s);
        if (!Gene.ANNO_LIST.contains(s)) throw new IllegalArgumentException();
        return selectedAnnotations.contains(s);
    }

    void setSelected(String s, boolean value) {
        assert Gene.ANNO_LIST.contains(s);
        if (!Gene.ANNO_LIST.contains(s)) throw new IllegalArgumentException();
        if (value) {
            assert !selectedAnnotations.contains(s);
            selectedAnnotations.add(s);
        } else {
            assert selectedAnnotations.contains(s);
            selectedAnnotations.remove(s);
        }
        update();
    }

    private File file;

    public void setFile(File file) {
        this.file = file;
    }

    /**
     * save all genes in geneList.
     */
    @Override
    protected Void doInBackground() throws Exception {// save処理
        assert file != null;
        final PrintWriter pw = new PrintWriter(file);
        pw.println();
        final StringBuilder line = new StringBuilder();
        for (final String prop : Gene.ANNO_LIST) {
            if (isSelected(prop)) {
                if (line.length() > 0) {
                    line.append('\t');
                }
                line.append(prop);
            }
        }
        pw.println(line);
        pw.flush();
        int cnt = 0;
        for (Gene gene : geneList) {
            if (isCancelled()) {
                break;
            }

            line.setLength(0);
            for (String prop : Gene.ANNO_LIST) {
                if (isSelected(prop)) {
                    if (line.length() > 0) {
                        line.append('\t');
                    }
                    line.append(gene.getAnnotation(prop));
                }
            }
            pw.println(line);
            setProgress(++cnt * 100 / geneList.size());
        }
        if (!isCancelled()) {
            setProgress(100);
        }
        pw.flush();
        pw.close();
        return null;
    }
}
