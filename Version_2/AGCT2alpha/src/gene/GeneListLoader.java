package gene;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class GeneListLoader {
    private static final String BEGIN = "!platform_table_begin";
    private final GeneList geneList;

    GeneListLoader(GeneList geneList) {
        this.geneList = geneList;
    }

    public void load(String path) throws Exception {
        File file = new File(path);
        BufferedReader br = new BufferedReader(new FileReader(file));
        String[] annos = null;
        for (; ; ) {

            String[] ss = br.readLine().split("\t");
            if (annos == null) {
                if (ss[0].equals(BEGIN)) {
                    annos = br.readLine().split("\t");
                    for (String anno : annos) {
                        if (!Gene.ANNO_LIST.contains(anno)) {
                            throw new Exception(String.format("%s is not a valid annotation name.", anno));
                        }
                    }
                }
                continue;
            }
            Gene gene = null;
            for (int i = 0; i < annos.length; i++) {
                if (annos[i].equals(Gene.ID)) {
                    gene = geneList.get(ss[i]);
                }
            }
            for (int i = 0; i < annos.length; i++) {
                gene.setAnnotation(annos[i], ss[i]);
            }
        }
    }
}
