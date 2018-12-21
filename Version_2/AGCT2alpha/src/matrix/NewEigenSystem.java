package matrix;

import forDebug.Debug;
import matrix.MatrixMath.Entry;
import test.AGCT;
import test.AGCTCounter;
import test.JInformationFrame;
import test.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class NewEigenSystem {

    /**
     * [v0,Av0,A^2v0,... ] が、T = PtAP なるP.
     *
     * @param A     三重対角化したい行列. 変更しない
     * @param alpha Tの対角成分
     * @param beta  Tの対角の隣
     * @return v0
     */
    private static Vector lanczos(AbstMatrix A, double[] alpha, double[] beta,
                                  int iter) {
        assert A.isSymmetry();
        int n = A.rowCount();
        AbstMatrix B = A.myClone();
        JInformationFrame.getInstance().setTextProgressBar("Lanczos");
        double[] random = new double[n];
        Random rand = new Random(14214987);
        double norm = 0;
        for (int i = 0; i < n; i++) {
            random[i] = rand.nextDouble();
            norm += random[i] * random[i];
        }
        norm = Math.sqrt(norm);
        Vector v0 = new SparseVector(n), v1 = new SparseVector(n);
        Vector res = new DenseVector(n);
        for (int i = 0; i < n; i++) {
            v1.set(i, random[i] / norm);
            res.set(i, random[i] / norm);
        }
        Vector w;
        beta[0] = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A.set(i, j, 0);
            }
        }
        for (int j = 0; j < iter; j++) {
            assert eq(v1.norm(), 1);
            JInformationFrame.getInstance().setValueProgressBar(
                    (int) (100. / j * iter));
            for (int i = 0; i < n; i++) {
                A.set(i, j, v1.get(i));
            }
            w = B.mul(v1).sub(v0.mul(beta[j]));
            alpha[j] = w.dot(v1);
            if (j == n - 1)
                break;
            w = w.sub(v1.mul(alpha[j]));
            beta[j + 1] = w.norm();
            if (beta[j + 1] == 0)
                break;
            v0 = v1;
            v1 = w.div(beta[j + 1]);
        }
        return res;
    }

//    private static int Numerical_Iterations_Max = 1000;

    private static double sqr(double a) {
        return ((Sqrarg = (a)) == 0.0 ? 0.0 : Sqrarg * Sqrarg);
    }

    static double Sqrarg;

    private static double sign(double a, double b) {
        return ((b) >= 0.0 ? Math.abs(a) : -Math.abs(a));
    }

    private static double pythag(double a, double b) {
        double absa, absb;
        absa = Math.abs(a);
        absb = Math.abs(b);
        if (absa > absb)
            return absa * Math.sqrt(1.0 + sqr(absb / absa));
        else
            return (absb == 0.0 ? 0.0 : absb * Math.sqrt(1.0 + sqr(absa / absb)));
    }

    private static boolean QR = false;

    private static void QRMethod(AbstMatrix mat, Vector v1, double[] alpha,
                                 double[] beta, int row) {

    }

    /**
     * A.row(i)が対応するeigenVectorになる。最大でも大きい方からMatrixConfig.MAX_COL_OF_EIGENVECTOR個のeigenVectorしか計算しない
     *
     * @param A
     * @param v1
     * @param alpha
     * @param beta
     * @param row
     */
    private static void calcEigenVectorOfTridiagonalizedMatrix(AbstMatrix A,
                                                               double[] alpha, double[] beta, int row) {
        //		if (QR) {
        //			QRMethod(A, v1, alpha, beta, row);
        //			return;
        //		}

        //		row = A.rowCount();

        int n = A.rowCount();
        AGCTCounter cc = new AGCTCounter("QL-Implicit Shifts", row);

        assert A.rowCount() == A.colCount();

        int m, L, iter, i, k;

        double s, r, p, g, f, dd, c, b, sf;
        for (i = 2; i <= n; i++)
            beta[i - 2] = beta[i - 1];
        beta[n - 1] = 0.0;

        for (L = 1; L <= row; L++) {
            cc.increment();

            iter = 0;
            do {
                for (m = L; m <= row - 1 - 1; m++) {
                    dd = Math.abs(alpha[m - 1]) + Math.abs(alpha[m]);
                    sf = (double) (Math.abs(beta[m - 1]) + dd);
                    if (sf == dd)
                        break;
                }
                if (m != L) {
                    if (iter++ == MatrixMath.Numerical_Iterations_Max)
                        Matrix.perror("Too many iterations in tqli");
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
                        for (k = 1; k <= n; k++) {// dimXまでが必須
                            f = A.get(k - 1, i);
                            A.set(k - 1, i, s * A.get(k - 1, i - 1) + c * f);
                            A.set(k - 1, i - 1, c * A.get(k - 1, i - 1) - s * f);
                        }
                    }
                    if (r == 0.0 && i >= L)
                        continue;
                    alpha[L - 1] -= p;
                    beta[L - 1] = g;
                    beta[m - 1] = 0.0;
                }
            } while (m != L);
        }
        if (AGCT.DoDebug)
            System.out.print("ok.\n");
        cc.end();
    }

    private static void lanczosEigenSort(AbstMatrix mat, double[] eigen, int iter) {//先頭からiterのみ計算する
        int n = mat.rowCount();
        assert mat.rowCount() == mat.colCount();
        Entry[] es = new Entry[iter];
        for (int i = 0; i < iter; i++) {
            es[i] = new Entry(i, eigen[i]);
        }
        Arrays.sort(es);
        ArrayList<Integer> same = new ArrayList<Integer>();
        ArrayList<Integer> diff = new ArrayList<Integer>();
        for (int i = 0; i < iter; ) {
            diff.add(es[i].i);
            int j = i + 1;
            while (j < iter && eq(eigen[es[i].i], eigen[es[j].i])) {
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

    /**
     * eigenValuesは基本的に降順に並ぶ。重複は後の方に置かれる。mat.col(i)が対応するeigenVectorになる
     *
     * @param mat
     * @param eigenValues
     */
    public static void getEigenSystem(AbstMatrix mat, double[] eigenValues) {// mat.row(i)が、対応するeigenVectorになる
        assert mat.isSymmetry();
        mat.isSymmetry();
        int n = mat.colCount();
        double[] beta = new double[n];
        //		AbstMatrix A = mat.myClone();
        int iter = Math.min(MatrixConfig.MAX_COL_OF_EIGENVECTOR, mat.rowCount());
        Debug.debug("lanczos...");
        lanczos(mat, eigenValues, beta, iter);
        System.gc();
        Debug.debug("calcEigenVectorOfTridiagonalizedMatirx...");
        calcEigenVectorOfTridiagonalizedMatrix(mat, eigenValues, beta, iter);
        System.gc();
        Debug.debug("eigenSort...");
        lanczosEigenSort(mat, eigenValues, iter);
        System.gc();
        Debug.debug("transpose...");
        transpose(mat);
        System.gc();
        Debug.debug("END getEigenSystem.");
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

    private static boolean eq(double a, double b) {
        //		if (Math.abs(a - b) < EPS)
        //			return true;
        if (Math.abs(1 - a / b) < EPS)
            return true;
        if (Math.abs(1 - b / a) < EPS)
            return true;
        return false;
    }

    private static final double EPS = 1e-2;

    public static void main(String[] args) {
        MatrixMath.newEigenSystem = false;
        MatrixConfig.MAX_COL_OF_EIGENVECTOR = 100;
        int n = 1000;
        //		String fileName = "test_Matrix" + n + "by" + n + ".txt";
        String fileName = "M38686.txt";
        final AbstMatrix N = MatrixFile.getFrom(fileName);
        n = N.colCount();
        AbstMatrix E = N.myClone();
        double[] eigenValues = new double[E.rowCount()];
        getEigenSystem(E, eigenValues);
        Debug.debug("my eigenValues", eigenValues);
        Vector[] eigenVectors = new Vector[E.rowCount()];
        for (int i = 0; i < E.rowCount(); i++) {
            eigenVectors[i] = E.getRow(i);
        }
        Debug.debug("eigenVector");
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < Math.min(100, n); j++) {
                System.err.print(eigenVectors[i].get(j) + ", ");
            }
            System.err.println();
            //			DoDebug.debug(eigenVectors2[i]);
        }
        //		System.out.println(testEigenVector(N, eigenVectors, eigenValues));

        //		String fileName = "M34.txt";
        //		final AbstMatrix N = MatrixFile.getFrom(fileName);
        E = N.myClone();
        double[] eigenValues2 = new double[E.rowCount()];
        MatrixMath.getEigenSystem(E, eigenValues2);
        Debug.debug("math eigenValues", eigenValues2);
        Vector[] eigenVectors2 = new Vector[E.rowCount()];
        for (int i = 0; i < E.rowCount(); i++) {
            eigenVectors2[i] = E.getRow(i);
        }
        Debug.debug("eigenVector2");
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < Math.min(100, n); j++) {
                System.err.print(eigenVectors2[i].get(j) + ", ");
            }
            System.err.println();
        }

        System.exit(0);
    }

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
                if (eq(a[i], b[j]))
                    bo = false;
            }
            if (bo) {
                Debug.debug(a[i]);
            }
        }
    }
}
