package test;

import junit.framework.TestCase;

public class MatrixTest extends TestCase {
    public void testSpasify() {
        double val = 0.02;
        AGCT.setSparsify_Statistic(val);
        Matrix m = new Matrix(10, 10);
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                m.set(i, j, i * 10 + j + 1);
            }
        }
        m.sparsifyMatrix();
        int cnt = 0;
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                if (m.get(i, j) == 0) {
                    cnt++;
                }
            }
        }
        System.out.println(m);
        assertEquals(val, (double) cnt / 100, 0.01);
    }
}
