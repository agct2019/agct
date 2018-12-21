package gene;

import test.AGCT;
import test.Util_Feature;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author Takashi
 */
public class GeneList extends ArrayList<Gene> implements Set<Gene> {

    private static final long serialVersionUID = 6342589729443428143L;
    public double coverage;
    public double evaluation;

    public int getNumberOfTimes() {
        return get(0).getRawCoordinate(0).length;
    }

    public int getNumberOfLigand() {
        return get(0).getNumberOfLigand();
    }

    public int getNumberOfDoses() {
        return get(0).getNumberOfDoses();
    }

    HashMap<String, Gene> map = new HashMap<String, Gene>();

    public Gene get(String name) {
        return map.get(name);
    }

    public Gene createAndAddGene(String name) {
        if (map.containsKey(name)) {
            return map.get(name);
        }
        final Gene gene = new Gene(name);
        map.put(name, gene);
        add(gene);
        return gene;
    }

    public boolean contains(String name) {
        return map.containsKey(name);
    }

    public HashSet<Gene> getResponsiveGeneSet() {
        final HashSet<Gene> res = new HashSet<Gene>();
        for (final Gene gene : this) {
            if (gene.isResponsive()) {
                res.add(gene);
            }
        }
        return res;
    }

    public double commonPercentile;

    public GeneList getSamePeakGenes(BitSet bitMask, double X) {
        final HashSet<Feature> set = Util_Feature.getXofUpperTotalVariationFeatures(this, X);
        final GeneList res = new GeneList();
        for (final Gene gene : this) {
            if (gene.hasSamePeak(bitMask)) {
                boolean contain = false;
                for (final Feature feature : gene.features) {
                    if (set.contains(feature)) {
                        contain = true;
                        break;
                    }
                }
                if (contain) {
                    res.add(gene);
                }
            }
        }
        return res;
    }

    // TODO refactoring
    public void saveAnnotationFile(File _file, BitSet bitMask, double X) {
        if (AGCT.MYDEBUG) {
            AGCT.debug(bitMask, X);
        }
        try {
            final FileWriter file = new FileWriter(_file);
            final GeneList samePeakGenes = getSamePeakGenes(bitMask, X);
            if (AGCT.MYDEBUG) {
                AGCT.debug(samePeakGenes.size());
            }
            for (final Gene gene : samePeakGenes) {
                if (AGCT.MYDEBUG) {
                    AGCT.debug(gene.getName() + "," + gene.features.get(bitMask.nextSetBit(0)).calcPeakTimeId() + "," + gene.geneTitle + "," + gene.getAnnotation(Gene.Gene$Symbol) + "," + gene.refSeq + "\n");
                }
                file.write(gene.getName() + "," + gene.features.get(bitMask.nextSetBit(0)).calcPeakTimeId() + "," + gene.geneTitle + "," + gene.getAnnotation(Gene.Gene$Symbol) + "," + gene.refSeq + "\n");
            }
            file.close();
        } catch (final IOException e) {
            if (AGCT.MYDEBUG) {
                AGCT.debug("IOERROR!");
            }
            System.exit(0);
            // TODO: handle exception
        }
    }

    private double percentageOfCoveringResponsiveGenes = -1;

    public void updateEnabledGene() {
        for (final Gene gene : this) {
            gene.setEnabled(gene.isEligible());
        }
    }

    public boolean allNamesAreDistinct() {
        final HashSet<String> names = new HashSet<String>();
        for (final Gene gene : this) {
            if (names.contains(gene.getName())) {
                return false;
            } else {
                names.add(gene.getName());
            }
        }
        return true;
    }

    public double getPercentageOfCoveringResponsiveGenes() {// [0,100]
        if (percentageOfCoveringResponsiveGenes != -1) {
            return percentageOfCoveringResponsiveGenes;
        }
        final HashSet<Gene> responsiveGenes = getResponsiveGenes();
        final HashSet<Gene> commonGenes = new HashSet<Gene>();
        for (final Gene gene : this) {
            if (responsiveGenes.contains(gene)) {
                commonGenes.add(gene);
            }
        }
        percentageOfCoveringResponsiveGenes = (double) commonGenes.size() / (double) (this.size() + responsiveGenes.size() - commonGenes.size()) * 100;
        return percentageOfCoveringResponsiveGenes;
    }

    HashSet<Gene> responsiveGenes = null;

    private HashSet<Gene> getResponsiveGenes() {
        if (responsiveGenes != null) {
            return responsiveGenes;
        }
        responsiveGenes = new HashSet<Gene>();
        for (final Gene gene : this) {
            if (gene.isResponsive()) {
                responsiveGenes.add(gene);
            }
        }
        return responsiveGenes;
    }

    /**
     * profile normalizationのメソッドの実装っぽい tk, 12/28
     * <p/>
     * circadianの（これで文法的にあってるのかな）メディアンをとってそれを全体から引いているということなんじゃないか
     */
    public void kanoWay() {
        // time 0が含まれていない場合は何もしない

        // time 0が存在するかの問い合わせ関数: Domainに存在。どうにかしてこっちにもってこよう。この関数はkanoの内包しよう。

        // 除外中
        boolean jogai = false;
        if (jogai) {
            final int numberOfLigand = get(0).getNumberOfLigand();
            final int numberOfTimes = get(0).getRawCoordinate(0).length;// TODO これらが入力ファイル中のどこなのかを確かめる。
            final int repelcate = get(0).getInput(0)[0].length;//
            for (int k = 0; k < size(); k++) {// number of genes
                // 下の４行のコメントを外してみます。-> privateを参照している。だめだった。仕様変更があったっぽい (12/28/12, tk)
//			get(k).initialCoordinates = new ArrayList<double[]>();
//			for(int i=0;i<get(k).rawCoordinates.size();i++){
//				get(k).initialCoordinates.add(get(k).rawCoordinates.get(i).clone());
//			}

                final double[] doze0 = new double[numberOfTimes * repelcate];
                for (int i = 0; i < numberOfTimes; i++) {
                    final double[] reps = new double[repelcate];
                    for (int j = 0; j < repelcate; j++) {
                        // 0を前提にしているんじゃないか
                        // get(k).getInput(0)[i][j]が何かを調べる。
                        reps[j] = doze0[i * repelcate + j] = get(k).getInput(0)[i][j];

                    }
                    Arrays.sort(reps);
                    final double sub = reps[repelcate / 2];// 注目しています。tk, 核心部
                    for (int j = 0; j < numberOfLigand; j++) {
                        get(k).getRawCoordinate(j)[i] -= sub;// ここで処理
                    }
                }
                Arrays.sort(doze0);
                final double absoluteLevel = doze0[numberOfTimes * repelcate / 2];
                for (int i = 0; i < numberOfTimes; i++) {
                    for (int j = 0; j < repelcate; j++) {
                        get(k).getRawCoordinate(j)[i] += absoluteLevel;
                    }
                }
            }
        }
    }

    public void subtractMedian() {
        final int numberOfLigand = get(0).getNumberOfLigand();
        final int numberOfTimes = get(0).getRawCoordinate(0).length;
        final double[][][] vals = new double[numberOfLigand][numberOfTimes][size()];
        for (int i = 0; i < numberOfLigand; i++) {
            for (int j = 0; j < numberOfTimes; j++) {
                for (int k = 0; k < size(); k++) {
                    vals[i][j][k] = get(k).getRawCoordinate(i)[j];
                }
                Arrays.sort(vals[i][j]);
                final double sub = vals[i][j][size() / 2];
                for (int k = 0; k < size(); k++) {
                    get(k).getRawCoordinate(i)[j] -= sub;// 上との違いを考察している、tk
                }
            }
        }
    }
}
