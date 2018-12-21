///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package parser;
//
////import listener.ProgressListener;
//
//import javax.swing.*;
//import java.io.File;
//import java.io.FileInputStream;
//import java.io.FileNotFoundException;
//import java.io.InputStream;
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Scanner;
//
///**
// * - GeneProfile fileについて
// * - タブ区切り
// * - "hoge" は,""なしで,hogeと書くことを意味する
// * - ()でくくられている文字が先頭にあれば,それと認識される.
// * - 順番はこのとおりでないといけない
// * - <hoge> は,hogeが意味するものを記述することを意味する
// * - 行の先頭にあるのがキーワードで,キーワード以降に,データを記述する.
// * - データ終了後にキーワードがこなかった場合,それはコメントと見なされ無視される.
// * - activation_i_j_k は,dose i,time j,replicate k におけるactivationを意味する.
// * goo
// * "(E)ntity"   <number of doses> <dose1>  <dose2>  ...
// * "(T)ime" <number of times> <time1>  <time2>  ...
// * "(R)eplicate"   <number of replicate>
// * "(D)ata"
// * <gene name1>   <activation_0_0_0>   <activation_0_0_1>   ...   <activation_0_1_0>   ...   <activation_1_0_0>   ...
// * <gene name2>   ...
// * ...
// * EOF
// *
// * @return
// * @throws java.io.IOException,java.util.IllegalFormatException
// */
//public class GeneProfileParser implements Parser<Void> {
//    private static final String WHITE = "\\p{javaWhitespace}+";
//    private File file;
////    private GeneProfile geneProfile;
//    private MyScanner sc;
////    private ProgressListener progressListener;
//
////    private GeneProfileParser(ProgressListener progressListener) {
////        this.progressListener = progressListener;
////    }
//
////    public GeneProfileParser() {
////        this(new ProgressListener() {
////
////            public void setProgress(int progress) {
////                do nothing
////            }
////
////            public void setMessage(String message) {
////                do nothing
////            }
////        });
////    }
//
////    public void setProgressListener(ProgressListener listener) {
////        this.progressListener = listener;
////    }
//
//    public void parse() throws FileParseException, FileNotFoundException {
//        if (file == null) throw new FileNotFoundException("file is null.");
//        if (!file.exists()) throw new FileNotFoundException("file does not exist.");
//        InputStream in = new FileInputStream(file);
//        sc = new MyScanner(in);
//        FileParseException info = parse0();
//        if (info != null) throw info;
//    }
//
//    private FileParseException make(String line, String message) {
//        return FileParseException.make(file, line, sc.getLineIndex(), message);
//    }
//
//    private void setProgress(String message, int progress) {
////        progressListener.setMessage(message);
////        progressListener.setProgress(progress);
//    }
//
////    public GeneProfile getResult() {
////        return geneProfile;
////    }
//    public Void getResult() {
//        return null;
//    }
//
//    public void setFile(File file) {
//        this.file = file;
//    }
//
//    private static class MyScanner {
//
//        private Scanner sc;
//        private int lineIndex = 0;
//        private boolean eof = false;
//
//        MyScanner(InputStream in) {
//            this.sc = new Scanner(in);
//        }
//
//        String nextLine() {
//            if (eof) return null;
//            lineIndex++;
//            if (!sc.hasNextLine()) {
//                eof = true;
//                return null;
//            }
//            return sc.nextLine();
//        }
//
//        int getLineIndex() {
//            return lineIndex;
//        }
//    }
//
//    private FileParseException parse0() {
//        // Entity
//        setProgress("parse entity", 0);
//        String line = sc.nextLine();
//        if (line == null) return make(line, "First row must start with \"Entity\".");
//        String[] elements = line.split(WHITE);
//        if (elements.length <= 0 || !line.startsWith("E")) return make(line, "First row must start with \"Entity\".");
//        if (elements.length <= 1) return make(line, "Entity is empty.");
//        String[] entites = new String[elements.length - 1];
//        for (int i = 0; i < entites.length; i++) entites[i] = elements[i + 1];
//        // Time
//        setProgress("parse time", 10);
//        line = sc.nextLine();
//        if (line == null) return make(line, "Second row must start with \"Time\".");
//        elements = line.split(WHITE);
//        if (elements.length <= 0 || !elements[0].startsWith("T"))
//            return make(line, "Second row must start with \"Time\".");
//        if (elements.length <= 1) return make(line, "Time is empty.");
//        String[] timeStrings = new String[elements.length - 1];
//        double[] times = new double[elements.length - 1];
//        for (int i = 0; i < timeStrings.length; i++) {
//            timeStrings[i] = elements[i + 1];
//            times[i] = parseTime(timeStrings[i]);
//            if (Double.isNaN(times[i])) {
//                return make(line, "Time format is wrong : ex. 2d, 10h, 30m.");
//            } else if (times[i] < 0) {
//                return make(line, "Time is negative");
//            }
//        }
//        // Replicate
//        setProgress("parse replicate", 20);
//        line = sc.nextLine();
//        if (line == null) return make(line, "Third row must start with \"Replicate\".");
//        String[] reps = line.split(WHITE);
//        if (reps.length <= 0 || !reps[0].startsWith("R")) return make(line, "Third row must start with \"Replicate\".");
//        if (reps.length <= 1) return make(line, "Relicate is empty");
//        int replicate = -1;
//        try {
//            replicate = Integer.valueOf(reps[1]);
//        } catch (NumberFormatException ignore) {
//        }
//        if (replicate <= 0) return make(line, "Replicate must be a positive integer");
//        // Gene
//        setProgress("parse data", 30);
//        line = sc.nextLine();
//        if (line == null) return make(line, "Data is empty");
//        List<Gene> dataList = new ArrayList<Gene>();
//        List<double[][][]> valueList = new ArrayList<double[][][]>();
//        int curProgress = 30;
//        do {
//            elements = line.split(WHITE);
//            Gene dat = new Gene();
//            // text[p] is already trimed, and text[p].length > 0. so ss.length > 0.
//            dat.setName(elements[0]);
//            if (elements.length - 1 < entites.length * times.length * replicate) {
//                return make(line, String.format(
//                        "Expected data value count is %d * %d * %d = %d, but observed is %d. Check your input information.",
//                        entites.length, times.length, replicate, entites.length * times.length * replicate, elements.length - 1));
//            }
//            double[][][] value = new double[entites.length][times.length][replicate];
//            int q = 1;
//            for (int j = 0; j < entites.length; j++) {
//                for (int k = 0; k < times.length; k++) {
//                    for (int l = 0; l < replicate; l++) {
//                        try {
//                            value[j][k][l] = Double.valueOf(elements[q++]);
//                        } catch (NumberFormatException ignore) {
//                            return make(line, "data value must be a number");
//                        }
//                    }
//                }
//            }
//            dataList.add(dat);
//            valueList.add(value);
//            // 30 + log_6(len) * 10
//            int nextProgress = (int) (30 + 10 * Math.log(dataList.size()) / Math.log(6.0));
//            if (nextProgress > curProgress) {
//                progressListener.setProgress(nextProgress);
//                curProgress = nextProgress;
//            }
//        } while ((line = sc.nextLine()) != null);
//        Gene[] data = dataList.toArray(new Gene[0]);
//        double[][][][] values = valueList.toArray(new double[0][][][]);
//        geneProfile = new GeneProfile(null, new ExperimentResult(entites, times, data, replicate, values));
//        setProgress("finished", 100);
//        return null;
//    }
//
//    /**
//     * return the corresponding time in minute.
//     * <p/>
//     * ex. "1m" -> 1, "3.14h" -> 3.14 * 60, "1.2d" -> 1.2 * 24 * 60,"0.0" -> 0.0
//     *
//     * @param time
//     * @return
//     */
//    // invalid -> NaN
//    private double parseTime(String time) {
//        try {
//            double t = Double.valueOf(time);
//            if (t == 0) return 0;
//        } catch (NumberFormatException ignore) {
//        }
//        if (time.length() < 2) return Double.NaN;
//        long mul;
//        switch (time.charAt(time.length() - 1)) {
//            case 'd':
//            case 'D':
//                mul = 60 * 24;
//                break;
//            case 'h':
//            case 'H':
//                mul = 60;
//                break;
//            case 'm':
//            case 'M':
//                mul = 1;
//                break;
//            default:
//                return Double.NaN;
//        }
//        try {
//            return Double.valueOf(time.substring(0, time.length() - 1)) * mul;
//        } catch (NumberFormatException ignore) {
//            return Double.NaN;
//        }
//    }
//
//    void error() {
//        try {
//
//        } catch (FileParseException ex) {
//            StringBuilder sb = new StringBuilder();
//            File file = ex.getFile();
//            if (file == null) throw new AssertionError();
//            String fileName = file.getName();
//            int lineIndex = ex.getLineIndex();
//            String line = ex.getLine();
//            String message = ex.getMessage();
//            sb.append("Parse error on line " + lineIndex + " of the file " + fileName + ".\n");
//            sb.append("line:\n");
//            sb.append(line + "\n");
//            sb.append("reason:\n");
//            sb.append(message);
//            JOptionPane.showMessageDialog(null, sb.toString().replace('\t', ' '));
//        } catch (FileNotFoundException ex) {
//            JOptionPane.showMessageDialog(null, ex.getMessage());
//        }
//    }
//}
//}
