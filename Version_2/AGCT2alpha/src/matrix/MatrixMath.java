package matrix;

import forDebug.Debug;
import test.AGCT;
import test.AGCTCounter;
import test.JInformationFrame;
import test.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class MatrixMath {

    private MatrixMath() {
    }

    /**
     * @param A     三重対角化したい行列. 実行後に PtAP = T なる P に書き換えられている.
     * @param alpha Tの対角成分
     * @param beta  Tの対角の隣
     */
    public static void lanczos(AbstMatrix A, double[] alpha, double[] beta) {
        int dimX = A.rowCount();
        int dimY = A.colCount();

        assert dimX == dimY;
        JInformationFrame.getInstance().setTextProgressBar("Lanczos");
        assert A.isSymmetry();
        AbstMatrix B = A.myClone();
        double[] random = new double[dimX];
        Random rand = new Random(0);
        double norm = 0;
        for (int i = 0; i < dimX; i++) {
            random[i] = rand.nextDouble();
            norm += random[i] * random[i];
        }
        norm = Math.sqrt(norm);
        Vector v0 = new SparseVector(dimX), v1 = new SparseVector(dimX);
        for (int i = 0; i < dimX; i++) {
            v1.set(i, random[i] / norm);
        }
        Vector w;
        beta[0] = 0;
        for (int i = 0; i < dimX; i++) {
            for (int j = 0; j < dimY; j++) {
                A.set(i, j, 0);
            }
        }
        for (int j = 0; j < dimY; j++) {
            assert eq(v1.norm(), 1);
            //			DoDebug.debug("lanczos", i);
            JInformationFrame.getInstance().setValueProgressBar(
                    (int) (100. / j * dimY));
            for (int i = 0; i < dimX; i++) {
                A.set(i, j, v1.get(i));
            }
            w = B.mul(v1).sub(v0.mul(beta[j]));
            alpha[j] = w.dot(v1);
            if (j == dimY - 1) {
                break;
            }
            w = w.sub(v1.mul(alpha[j]));
            beta[j + 1] = w.norm();
            if (beta[j + 1] == 0) {
                break;
            }
            v0 = v1;
            v1 = w.div(beta[j + 1]);
        }
    }

    private final static boolean lanczos = true;

    private static boolean useLanczos(int dim) {
        return lanczos && dim > MatrixConfig.MAX_COL_OF_EIGENVECTOR;
    }

    public static void Tridiagonalize(AbstMatrix mat, double d[], double e[]) {
        if (useLanczos(d.length)) {
            lanczos(mat, d, e);
            return;
        }

        int dimX = mat.rowCount();
        int dimY = mat.colCount();
        AGCTCounter cc = new AGCTCounter("tred2", 2 * dimX - 2);

        if (AGCT.DoDebug) {
            System.out.print("Householder tridiagonalization on Matrix " + " ("
                    + dimX + " x " + dimY + ") (start)... ");
        }

        assert mat.isSymmetry();

        int l, k, j, i, n = dimX;
        double scale, hh, h, g, f;

        for (i = n; i >= 2; i--) {
            cc.increment();

            l = i - 1;
            h = scale = 0.0;
            if (l > 1) {
                for (k = 1; k <= l; k++) {
                    scale += Math.abs(mat.get(i - 1, k - 1));
                }
                if (scale == 0.0) {
                    e[i - 1] = mat.get(i - 1, l - 1);
                } else {
                    for (k = 1; k <= l; k++) {
                        mat.set(i - 1, k - 1, mat.get(i - 1, k - 1) / scale); // coordinates[i
                        // -
                        // 1][k
                        // - 1]
                        // /=
                        // scale;
                        h += mat.get(i - 1, k - 1) * mat.get(i - 1, k - 1);
                    }
                    f = mat.get(i - 1, l - 1);
                    g = (f >= 0.0 ? -Math.sqrt(h) : Math.sqrt(h));
                    e[i - 1] = scale * g;
                    h -= f * g;
                    mat.set(i - 1, l - 1, f - g);
                    f = 0.0;
                    for (j = 1; j <= l; j++) {
                        mat.set(j - 1, i - 1, mat.get(i - 1, j - 1) / h);
                        g = 0.0;
                        for (k = 1; k <= j; k++) {
                            g += mat.get(j - 1, k - 1) * mat.get(i - 1, k - 1);
                        }
                        for (k = j + 1; k <= l; k++) {
                            g += mat.get(k - 1, j - 1) * mat.get(i - 1, k - 1);
                        }
                        e[j - 1] = g / h;
                        f += e[j - 1] * mat.get(i - 1, j - 1);
                    }
                    hh = f / (h + h);
                    for (j = 1; j <= l; j++) {
                        f = mat.get(i - 1, j - 1);
                        e[j - 1] = g = e[j - 1] - hh * f;
                        for (k = 1; k <= j; k++) {
                            mat.set(j - 1, k - 1, mat.get(j - 1, k - 1)
                                    - (f * e[k - 1] + g * mat.get(i - 1, k - 1)));
                        }
                        // coordinates[j - 1][k - 1] -= (f * e[k - 1] + g *
                        // get(i - 1,k - 1));
                    }
                }
            } else {
                e[i - 1] = mat.get(i - 1, l - 1);
            }
            d[i - 1] = h;
        }
        d[0] = 0.0;
        e[0] = 0.0;
        for (i = 1; i <= n; i++) {
            cc.increment();

            l = i - 1;
            if (d[i - 1] != 0.0) {
                for (j = 1; j <= l; j++) {
                    g = 0.0;
                    for (k = 1; k <= l; k++) {
                        g += mat.get(i - 1, k - 1) * mat.get(k - 1, j - 1);
                    }
                    for (k = 1; k <= l; k++) // coordinates[k - 1][j - 1] -= g * get(k - 1,i - 1);
                    {
                        mat.set(k - 1, j - 1, mat.get(k - 1, j - 1) - g
                                * mat.get(k - 1, i - 1));
                    }
                }
            }
            d[i - 1] = mat.get(i - 1, i - 1);
            mat.set(i - 1, i - 1, 1.0);
            for (j = 1; j <= l; j++) {
                mat.set(j - 1, i - 1, 0);
                mat.set(i - 1, j - 1, 0.0);
            }
        }

        if (AGCT.DoDebug) {
            System.out.print("ok.\n");
        }
        cc.end();
    }

    static int Numerical_Iterations_Max = 1000;

    public static int getNumericalIterationsMax() {
        return Numerical_Iterations_Max;
    }

    public static void setNumerical_Iterations_Max(int num) {
        Numerical_Iterations_Max = num;
    }

    private static double sqr(double a) {
        return a * a;
    }

    private static double sign(double a, double b) {
        return ((b) >= 0.0 ? Math.abs(a) : -Math.abs(a));
    }

    private static double pythag(double a, double b) {
        double absa, absb;
        absa = Math.abs(a);
        absb = Math.abs(b);
        if (absa > absb) {
            return absa * Math.sqrt(1.0 + sqr(absb / absa));
        } else {
            return (absb == 0.0 ? 0.0 : absb * Math.sqrt(1.0 + sqr(absa / absb)));
        }
    }

    //	private static boolean QR = false;
//
//	private static void QRMethod(AbstMatrix mat, double[] alpha, double[] beta) {
//
//	}
    public static void calcEigenVectorOfTridiagonalizedMatrix(AbstMatrix mat,
                                                              double[] alpha, double[] beta) {
        int dimX = mat.rowCount();
        AGCTCounter cc = new AGCTCounter("QL-Implicit Shifts", dimX);

        assert mat.rowCount() == mat.colCount();

        int m, L, iter, i, k, n = dimX;
        double s, r, p, g, f, dd, c, b, sf;
        for (i = 2; i <= n; i++) {
            beta[i - 2] = beta[i - 1];
        }
        beta[n - 1] = 0.0;

        for (L = 1; L <= n; L++) {
            cc.increment();

            iter = 0;
            do {
                for (m = L; m <= n - 1; m++) {
                    dd = Math.abs(alpha[m - 1]) + Math.abs(alpha[m]);
                    sf = (double) (Math.abs(beta[m - 1]) + dd);
                    if (sf == dd) {
                        break;
                    }
                }
                if (m != L) {
                    if (iter++ == Numerical_Iterations_Max) {
                        Matrix.perror("Too many iterations in tqli");
                    }
                    g = (alpha[L] - alpha[L - 1]) / (2.0 * beta[L - 1]);
                    r = pythag(g, 1.0);
                    g = alpha[m - 1] - alpha[L - 1] + beta[L - 1] / (g + sign(r, g));
                    s = c = 1.0;
                    p = 0.0;
                    for (i = m - 1; i >= L; i--) {
                        f = s * beta[i - 1];
                        b = c * beta[i - 1];
                        beta[i] = (r = pythag(f, g));
                        if (r == 0.0) {
                            alpha[i] -= p;
                            beta[m - 1] = 0.0;
                            break;
                        }
                        s = f / r;
                        c = g / r;
                        g = alpha[i] - p;
                        r = (alpha[i - 1] - g) * s + 2.0 * c * b;
                        alpha[i] = g + (p = s * r);
                        g = c * r - b;
                        for (k = 1; k <= n; k++) {
                            f = mat.get(k - 1, i);
                            mat.set(k - 1, i, s * mat.get(k - 1, i - 1) + c * f);
                            mat.set(k - 1, i - 1, c * mat.get(k - 1, i - 1) - s * f);
                        }
                    }
                    if (r == 0.0 && i >= L) {
                        continue;
                    }
                    alpha[L - 1] -= p;
                    beta[L - 1] = g;
                    beta[m - 1] = 0.0;
                }
            } while (m != L);
        }
        if (AGCT.DoDebug) {
            System.out.print("ok.\n");
        }
        cc.end();
    }
//

    private static void lnczosEigenSort(AbstMatrix mat, double[] eigen) {
        int n = mat.rowCount();
        assert mat.rowCount() == mat.colCount();
        Entry[] es = new Entry[n];
        for (int i = 0; i < n; i++) {
            es[i] = new Entry(i, eigen[i]);
        }
        Arrays.sort(es);
        ArrayList<Integer> same = new ArrayList<Integer>();
        ArrayList<Integer> diff = new ArrayList<Integer>();
        for (int i = 0; i < n; ) {
            diff.add(es[i].i);
            int j = i + 1;
            while (j < n && eq(eigen[es[i].i], eigen[es[j].i])) {
                same.add(es[j].i);
                j++;
            }
            i = j;
        }
        AbstMatrix dupM = mat.myClone();
        double[] dupE = eigen.clone();
        for (int j = 0; j < diff.size(); j++) {
            for (int i = 0; i < n; i++) {
                mat.set(i, j, dupM.get(i, diff.get(j)));
            }
            eigen[j] = dupE[diff.get(j)];
        }
        for (int j = 0; j < same.size(); j++) {
            for (int i = 0; i < n; i++) {
                mat.set(i, j + diff.size(), dupM.get(i, same.get(j)));
            }
            eigen[j + diff.size()] = dupE[same.get(j)];
        }
    }

    static class Entry implements Comparable<Entry> {// 降順にソート

        int i;
        double d;

        Entry(int i, double d) {
            this.d = d;
            this.i = i;
        }

        @Override
        public int compareTo(Entry o) {
            return (int) Math.signum(o.d - d);
        }
    }

    public static void eigenSort(AbstMatrix mat, double[] eigen) {
        if (useLanczos(eigen.length)) {
            lnczosEigenSort(mat, eigen);
            return;
        }

        if (AGCT.DoDebug) {
            System.out.print("Sorting eigenvalues" + " by descending order... ");
        }

        assert mat.rowCount() == mat.colCount();

        int k, j, i, n = mat.rowCount();
        double p;
        for (i = 0; i <= n - 1; i++) {
            p = eigen[k = i];
            for (j = i + 1; j <= n - 1; j++) {
                if (eigen[j] >= p) {
                    p = eigen[k = j];
                }
            }
            if (k != i) {
                eigen[k] = eigen[i];
                eigen[i] = p;
                for (j = 0; j <= n - 1; j++) {
                    p = mat.get(j, i);
                    mat.set(j, i, mat.get(j, k));
                    mat.set(j, k, p);
                }
            }
        }
        if (AGCT.DoDebug) {
            System.out.print("ok.\n");
        }
    }

    static void transpose(AbstMatrix mat) {
        assert mat.rowCount() == mat.colCount();
        int n = mat.rowCount();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double d = mat.get(i, j);
                mat.set(i, j, mat.get(j, i));
                mat.set(j, i, d);
            }
        }
    }

    static boolean newEigenSystem = true;

    /**
     * eigenValuesは基本的に降順に並ぶ。重複は後の方に置かれる。mat.col(i)が対応するeigenVectorになる
     *
     * @param mat
     * @param eigenValues
     */
    public static void getEigenSystem(AbstMatrix mat, double[] eigenValues) {// mat.row(i)が、対応するeigenVectorになる
        if (newEigenSystem) {
            NewEigenSystem.getEigenSystem(mat, eigenValues);
            return;
        }

        assert mat.isSymmetry();
        mat.isSymmetry();
        int n = mat.colCount();
        double[] d = new double[n];
        Tridiagonalize(mat, eigenValues, d);
        calcEigenVectorOfTridiagonalizedMatrix(mat, eigenValues, d);
        eigenSort(mat, eigenValues);
        transpose(mat);
    }
    /*
    public static void main(String[] args) {
    int n = 50;
    String fileName = "test_Matrix" + n + "by" + n + ".txt";
    //		lanczos = false;
    //		String fileName = "M10.txt";
    final AbstMatrix N = MatrixFile.getFrom(fileName);
    AbstMatrix E = N.myClone();
    double[] eigenValues = new double[E.rowCount()];
    getEigenSystem(E, eigenValues);
    DoDebug.debug("tred2 eigenValues", eigenValues);
    Vector[] eigenVectors = new Vector[E.rowCount()];
    for (int i = 0; i < E.rowCount(); i++) {
    eigenVectors[i] = E.getRow(i);
    }
    //		System.out.println(testEigenVector(N, eigenVectors, eigenValues));

    //		String fileName = "M34.txt";
    //		final AbstMatrix N = MatrixFile.getFrom(fileName);
    E = N.myClone();
    //		lanczos = true;
    double[] eigenValues2 = new double[E.rowCount()];
    getEigenSystem(E, eigenValues2);
    DoDebug.debug("lanczos eigenValues", eigenValues2);
    eigenVectors = new Vector[E.rowCount()];
    for (int i = 0; i < E.rowCount(); i++) {
    eigenVectors[i] = E.getRow(i);
    }
    DoDebug.debug("tred2 diff");
    showDiff(eigenValues, eigenValues2);
    DoDebug.debug("lanczos diff");
    showDiff(eigenValues2, eigenValues);

    //		System.out.println(testEigenVector(N, eigenVectors, eigenValues));

    System.exit(0);
    }
     */

    private static void showDiff(double[] a, double[] b) {
        int n = a.length;
        for (int i = 0; i < n; i++) {
            boolean bo = true;
            for (int j = 0; j < n; j++) {
                if (eq(a[i], b[j])) {
                    bo = false;
                }
            }
            for (int j = 0; j < i; j++) {
                if (eq(a[i], b[j])) {
                    bo = false;
                }
            }
            if (bo) {
                Debug.debug(a[i]);
            }
        }
    }

    private static double EPS = 1e-5;

    private static boolean testEigenVector(AbstMatrix mat, Vector[] eigenVectors,
                                           double[] eigenValues) {
        boolean res = true;
        for (int i = 0; i < eigenValues.length; i++) {
            Debug.debug("eigenValue", eigenValues[i]);
            Debug.debug("eigenVector", eigenVectors[i]);
            if (!(Math.abs(1 - eigenVectors[i].norm()) < EPS)) {
                Debug.debug(eigenVectors[i].norm());
            }
            assert eq(eigenVectors[i].norm(), 1);
            Vector v1 = mat.mul(eigenVectors[i]);
            Vector v2 = eigenVectors[i].mul(eigenValues[i]);
            if (Math.abs(v1.sub(v2).norm()) > EPS) {
                Debug.debug("v1", v1);
                Debug.debug("v2", v2);
                res = false;
            }
        }
        return res;
    }

    private static boolean eq(double a, double b) {
        return Math.abs(a - b) < EPS;
    }
}
