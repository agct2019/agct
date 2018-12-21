package parser;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.IllegalFormatException;
import java.util.Scanner;

public class RawDataParser {
    public void check(Scanner sc) throws IOException, IllegalFormatException, ParseException {

    }

    private void check2(Scanner sc) throws IOException, IllegalFormatException, ParseException {
        while (!sc.next().startsWith("S")) ;
        // Species
        String annotationFileName = sc.next();
//        String realAnnotationFilePath = curDir.getAbsolutePath()+"/"+annotationFileName;
//	realAnnotationFilePath = GOFileName;

        while (!sc.next().startsWith("E")) ;
        // Entities
        int numberOfDoses = sc.nextInt();
        String[] doses = new String[numberOfDoses];
        for (int i = 0; i < numberOfDoses; i++)
            doses[i] = sc.next();

        while (!sc.next().startsWith("T")) ;
        // Times
        int numberOfTimes = sc.nextInt();
        double[] times = new double[numberOfTimes];
        for (int i = 0; i < numberOfTimes; i++)
            times[i] = sc.nextDouble();

        while (!sc.next().startsWith("R")) ;
        // Replicate
        int replicate = sc.nextInt();

        while (!sc.nextLine().startsWith("D")) ;

//        int numGene = sc.nextInt();
        // Datas,Genes
        ArrayList<double[][][]> activationList = new ArrayList<double[][][]>();
        ArrayList<String> geneIdList = new ArrayList<String>();
        HashMap<String, Integer> dejaGeneSet = new HashMap<String, Integer>();
        while (sc.hasNext()) {
            String geneId = sc.next();
            if (dejaGeneSet.containsKey(geneId)) {
                int cnt = dejaGeneSet.get(geneId);
                cnt++;
                dejaGeneSet.put(geneId, cnt);
                geneId += "(" + cnt + ")";
            } else {
                dejaGeneSet.put(geneId, 1);
            }
            geneIdList.add(geneId);
            double[][][] activation = new double[numberOfDoses][numberOfTimes][replicate];
            for (int i = 0; i < numberOfDoses; i++)
                for (int j = 0; j < numberOfTimes; j++)
                    for (int k = 0; k < replicate; k++)
                        activation[i][j][k] = sc.nextDouble();
            activationList.add(activation);
        }
//        assert activationList.size()==geneIdList.size();
//        int numberOfGenes = geneIdList.size();
//        Gene[] genes = Gene.read(geneIdList.toArray(new String[0]), realAnnotationFilePath);
//        double[][][][] activations = new double[numberOfGenes][numberOfDoses][numberOfTimes][replicate];
//        for (int i = 0; i<numberOfGenes; i++)
//            activations[i] = activationList.get(i);
//        return new GeneProfile(annotationFileName, doses, times, genes, replicate, activations);
    }
}
