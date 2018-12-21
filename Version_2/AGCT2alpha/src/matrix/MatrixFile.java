package matrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class MatrixFile {
    private MatrixFile() {
    }

    public static boolean writeTo(AbstMatrix mat, String fileName) {
        try {
            FileWriter fw = new FileWriter(fileName);
            int m = mat.rowCount(), n = mat.colCount();
            fw.write(m + " " + n + "\n");
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    double val = mat.get(i, j);
                    if (val != 0) {
                        fw.write(i + " " + j + " " + val + "\n");
                    }
                }
            }
            fw.flush();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            return false;
        }
        return true;
    }

    /**
     * テスト用入力をfilenameからとってくる 形式: m n x1 y1 val ...
     *
     * @param sc
     */
    public static AbstMatrix getFrom(String filename) {
        //		Scanner sc = null;
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String[] ss = br.readLine().split(" ");
//			Debug.debug(ss);
            int m = Integer.valueOf(ss[0]), n = Integer.valueOf(ss[1]);
            AbstMatrix res = new GeneralMatrix(m, n);
            while (true) {
                String s = null;
                try {
                    s = br.readLine();
                } catch (IOException e) {
                }
                if (s == null)
                    break;
                ss = s.split(" ");
//				Debug.debug(ss);
                res.set(Integer.valueOf(ss[0]), Integer.valueOf(ss[1]), Double
                        .valueOf(ss[2]));
            }
            return res;
        } catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
    }
}
