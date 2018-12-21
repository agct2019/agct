package combinationRun;

import gene.Gene;
import gene.GeneList;
import test.AGCT;
import test.Domain;
import test.Matrix;

import javax.swing.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class IterativeRun {
    public static enum SortMode {
        CommonGenes, KBCoverage;
    }

    public SortMode sortMode = SortMode.KBCoverage;

    public GeneList initialGenes;
    private Domain domain;

    public void setDomain(Domain domain) {
        this.domain = domain;
    }

    double[] percentiles = {0.25, 0.50, 0.75};
    ArrayList<State> list = new ArrayList<State>();

    public void iterativeRun() {
        HashSet<Gene> responsiveGeneSet = initialGenes.getResponsiveGeneSet();
        if (AGCT.MYDEBUG)
            AGCT.debug("responsiveGeneSet.size() : " + responsiveGeneSet.size());
        // if (AGCT.MYDEBUG) {
        // AGCT.debug("responsiveGeneSet");
        // for (Gene gene : responsiveGeneSet) {
        // if (AGCT.MYDEBUG)
        // AGCT.debug(gene.name);
        // }
        // AGCT.debug();
        // }

        int fileNumber = domain.getFileCount();
        int n = initialGenes.get(0).getNumberOfLigand() / fileNumber;//ひとつのfileにおけるligandの数
        for (int i = 0; i < 1 << n; i++)
            if (Integer.bitCount(i) >= 2) {
                int mask = 0;
                for (int k = 0; k < fileNumber; k++) {
                    mask <<= n;
                    mask |= i;
                }
                BitSet bitMask = getBitSet(mask);
                for (double percentile : percentiles) {
                    GeneList samePeakGenes = initialGenes.getSamePeakGenes(bitMask, percentile);
                    int common = 0;
                    for (Gene gene : samePeakGenes) {
                        if (responsiveGeneSet.contains(gene)) {
                            common++;

                            if (AGCT.MYDEBUG)
                                AGCT.debug(gene.getName());
                        }
                    }
                    double val = 0;
                    samePeakGenes.commonPercentile = (double) samePeakGenes.size() / initialGenes.size();
                    samePeakGenes.coverage = (double) common / (responsiveGeneSet.size());
                    list.add(new State(bitMask, percentile, samePeakGenes));
                    if (sortMode == SortMode.CommonGenes) {
                        val = samePeakGenes.commonPercentile;
                    } else if (sortMode == SortMode.KBCoverage) {
                        val = samePeakGenes.coverage;
                    } else {
                        if (AGCT.MYDEBUG)
                            AGCT.debug("IterativeRun.iterativeRun");
                        System.exit(0);
                    }
                    samePeakGenes.evaluation = val;
                }
            }
        Collections.sort(list, new Comparator<State>() {
            @Override
            public int compare(State o1, State o2) {
                // TODO Auto-generated method stub
                if (o1.percentile != o2.percentile)
                    return (int) Math.signum(o1.percentile - o2.percentile);
                return -(int) Math.signum(o1.geneList.evaluation - o2.geneList.evaluation);
            }
        });

        saveFile("IterativeTest_" + domain.getFileName() + ".txt");
        String statisticsOfOtherResponsiveGenes = "StatisticsOfResponsive_" + domain.getFileName() + ".txt";
        State max = list.get(0);
        initialGenes.saveAnnotationFile(new File(max.bitSet + "_" + max.percentile * 100 + ".txt"), max.bitSet, max.percentile);

        domain.loadHighlights(new File(max.bitSet + "_" + max.percentile * 100 + ".txt"));

        saveStatisticsOfResponsiveGenes(statisticsOfOtherResponsiveGenes, responsiveGeneSet, max.bitSet);

        String statisticsOfHighlightedGenes = "StatisticsOfHighlight_" + domain.getFileName() + ".txt";
        saveStatisticsOfHighlightedGenes(statisticsOfHighlightedGenes, max);
        final StringBuilder message = new StringBuilder();
        message.append("max bitSet: " + max.bitSet + "\n");
        message.append("max percentile: " + max.percentile * 100 + "\n");
        message.append("percentage of chosen genes: " + String.format("%.3f%%\n", max.geneList.commonPercentile * 100));
        message.append("KB Coverage: " + String.format("%.3f%%\n", max.geneList.coverage * 100));
        message.append("loaded highlight File: " + max.bitSet + "_" + max.percentile * 100 + ".txt\n");
        message.append("saved test result on " + "\"IterativeTest_" + domain.getFileName() + ".txt\"\n");
        message.append("saved statistics of responsive genes on \"" + statisticsOfOtherResponsiveGenes + "\"\n");
        message.append("saved statistics of higihlighted genes on \"" + statisticsOfHighlightedGenes + "\"\n");
        new Thread() {
            public void run() {
                JOptionPane.showMessageDialog(null, message.toString());
            }
        }.start();
    }

    boolean ok = false;

    int toi(BitSet bitset) {
        int res = 0;
        int i = 0;
        while (i < bitset.length()) {
            if (bitset.get(i)) {
                res |= 1 << i;
                // bitset.clear(i);
            }
            i++;
        }
        return res;
    }

    private void saveStatisticsOfHighlightedGenes(String fileName, State max) {
        HashSet<Gene> responsiveGeneSet = initialGenes.getResponsiveGeneSet();
        GeneList list = initialGenes.getSamePeakGenes(max.bitSet, max.percentile);
        File file = new File(fileName);
        int numberOfLigands = initialGenes.getNumberOfLigand();
        int numberOfTimes = initialGenes.getNumberOfTimes();
        // ArrayList<Gene> list = new ArrayList<Gene>();
//		for (Gene gene : responsiveGeneSet)
//			list.add(gene);
        Collections.sort(list, new Comparator<Gene>() {
            @Override
            public int compare(Gene o1, Gene o2) {
                // TODO 自動生成されたメソッド・スタブ
                return o1.getName().compareTo(o2.getName());
                // return 0;
            }
        });
        try {
            FileWriter fr = new FileWriter(file);
            fr.write("Genename\tGene title\tGene symbol\tresponsive Gene");
            for (int i = 0; i < 32; i++)
                if (max.bitSet.get(i)) {
                    // for (int j = 0; j < numberOfTimes; j++) {
                    // fr.write("\tligand" + i + "_time" + j);
                    // }
                    fr.write("\tligand" + i + "_TotalVariation");
                }
            for (Gene gene : list) {
                fr.write("\n" + gene.getName());
                fr.write("\t" + gene.geneTitle);
                fr.write("\t" + gene.getAnnotation(Gene.Gene$Symbol));
                if (gene.isResponsive()) {
                    fr.write("\tYES");
                    // if (gene.hasSamePeak(max.bitSet)) {
                    // fr.write("\t-");
                    // } else {
                    // Matrix.perror("IterativeRun.saveStaticticsOfResponsiveGenes");
                    // }
                } else {
                    fr.write("\tNO");
                    // fr.write("\t");
                    // fr.write(gene.hasSamePeak(max.bitSet) ? "YES" : "NO");
                }
                for (int i = 0; i < 32; i++)
                    if (max.bitSet.get(i)) {
                        // for (int j = 0; j < numberOfTimes; j++) {
                        // fr.write("\t" + gene.getRawCoordinate(i)[j]);
                        // }
                        fr.write("\t" + gene.features.get(i).calcTotalVariation());
                    }
            }
            fr.flush();
        } catch (IOException e) {
            // TODO: handle exception
            AGCT.debug("failed in IterativeRun.saveStatisticsOfResponsive");
        }
    }

    private void saveStatisticsOfResponsiveGenes(String fileName, HashSet<Gene> responsiveGeneSet, BitSet ligandBitmask) {
        File file = new File(fileName);
        int numberOfLigands = initialGenes.getNumberOfLigand();
        int numberOfTimes = initialGenes.getNumberOfTimes();
        ArrayList<Gene> list = new ArrayList<Gene>();
        for (Gene gene : responsiveGeneSet)
            list.add(gene);
        Collections.sort(list, new Comparator<Gene>() {
            @Override
            public int compare(Gene o1, Gene o2) {
                // TODO 自動生成されたメソッド・スタブ
                if (o1.isReferenced() != o2.isReferenced()) {
                    return o1.isReferenced() ? 1 : -1;
                }
                return o1.getName().compareTo(o2.getName());
                // return 0;
            }
        });
        try {
            // System.setOut(new PrintStream(file));
            FileWriter fr = new FileWriter(file);
            fr.write("Genename\tGene title\tGene symbol\tcontained\tsame peak");
            for (int i = 0; i < 32; i++)
                if (ligandBitmask.get(i)) {
                    for (int j = 0; j < numberOfTimes; j++) {
                        fr.write("\tligand" + i + "_time" + j);
                    }
                    fr.write("\tligand" + i + "_TotalVariation");
                }
            for (Gene gene : list) {
                fr.write("\n" + gene.getName());
                fr.write("\t" + gene.geneTitle);
                fr.write("\t" + gene.getAnnotation(Gene.Gene$Symbol));
                if (gene.isReferenced()) {
                    fr.write("\tYES");
                    if (gene.hasSamePeak(ligandBitmask)) {
                        fr.write("\t-");
                    } else {
                        Matrix.perror("IterativeRun.saveStaticticsOfResponsiveGenes");
                    }
                } else {
                    fr.write("\tNO");
                    fr.write("\t");
                    fr.write(gene.hasSamePeak(ligandBitmask) ? "YES" : "NO");
                }
                for (int i = 0; i < 32; i++)
                    if (ligandBitmask.get(i)) {
                        for (int j = 0; j < numberOfTimes; j++) {
                            fr.write("\t" + gene.getRawCoordinate(i)[j]);
                        }
                        fr.write("\t" + gene.features.get(i).calcTotalVariation());
                    }
            }
            fr.flush();
        } catch (IOException e) {
            // TODO: handle exception
            AGCT.debug("failed in IterativeRun.saveStatisticsOfResponsive");
        }
    }

    private void saveFile(String fileName) {
        File file = new File(fileName);
        ArrayList<State> tmp = new ArrayList<State>();
        for (int i = 0; i < list.size(); i++) {
            tmp.add(list.get(i));
        }
        Collections.sort(tmp, new Comparator<State>() {
            @Override
            public int compare(State o1, State o2) {
                // TODO Auto-generated method stub
                if (o1.percentile != o2.percentile)
                    return (int) Math.signum(o1.percentile - o2.percentile);
                if (o1.geneList.coverage != o2.geneList.coverage)
                    return (int) Math.signum(o1.geneList.coverage - o2.geneList.coverage);
                return toi(o1.bitSet) - toi(o2.bitSet);

                // return toi(o1.bitSet) - toi(o2.bitSet);
            }
        });
        try {
            FileWriter fr = new FileWriter(file);
            fr.write("bitset\tKB coverage\tcommon\tpercentile\n");
            for (State state : tmp) {

                fr.write(state.bitSet + "\t" + state.geneList.coverage * 100 + "\t" + state.geneList.commonPercentile * 100 + "\t" + state.percentile * 100 + "\n");
                if (AGCT.MYDEBUG)
                    AGCT.debug(state.bitSet + "\t" + state.geneList.coverage * 100 + "\t" + state.geneList.commonPercentile * 100 + "\t" + state.percentile * 100 + "\n");
            }
            fr.flush();
        } catch (IOException e) {
            // TODO: handle exception
            e.printStackTrace();
        }
    }

    private class State {
        BitSet bitSet;
        double percentile;
        GeneList geneList;

        State(BitSet bitSet, double percentile, GeneList geneList) {
            this.bitSet = bitSet;
            this.percentile = percentile;
            this.geneList = geneList;
        }
    }

    private BitSet getBitSet(int mask) {
        BitSet res = new BitSet();
        for (int i = 0; i < 31; i++) {
            if (((mask >> i) & 1) == 1)
                res.set(i);
        }
        return res;
    }
}
