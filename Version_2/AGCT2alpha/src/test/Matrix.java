package test;

import forDebug.Debug;
import gene.Gene;
import matrix.*;
import matrix.Vector;

import javax.swing.*;
import java.io.FileWriter;
import java.util.*;

/**
 * class matrix
 *
 * @author Takashi(comment only)
 */
public class Matrix implements Debuggable, AbstMatrix {

    @Override
    public int colCount() {
        return mat.colCount();
    }

    @Override
    public Vector getCol(int y) {
        return mat.getCol(y);
    }

    @Override
    public Vector getRow(int x) {
        return mat.getRow(x);
    }

    public boolean isSymmetry() {
        return mat.isSymmetry();
    }

    public AbstMatrix mul(AbstMatrix mat) {
        return mat.mul(mat);
    }

    public Vector mul(Vector v) {
        return mat.mul(v);
    }

    public AbstMatrix myClone() {
        return mat.myClone();
    }

    public int nonZeroCount() {
        return mat.nonZeroCount();
    }

    public int rowCount() {
        return mat.rowCount();
    }

    public AbstMatrix turn() {
        return mat.turn();
    }

    String name;

    private AbstMatrix mat;

    public double get(int x, int y) {
        return mat.get(x, y);
    }

    public void set(int x, int y, double d) {
        mat.set(x, y, d);
    }

    Matrix(int dimX, int dimY) {
        this.dimX = dimX;
        this.dimY = dimY;
        mat = new GeneralMatrix(dimX, dimY);
    }

    Matrix choleski_decomposition;

    private final int dimX;// number of rows
    private final int dimY;// number of columns

    boolean isVector; // true iff dimX = 1 (vectors are ROW)
    boolean hasDefaultContent; // true iff the content has defaultValue
    boolean isSquared; // true iff dimX = dimY
    boolean isSymmetric; // true iff isSquared && symmetric

    static double Tiny = 1.0E-30;

    /**
     * *************************************************************************************************
     * Class Methods
     * ***
     */

    public static int Hamming(int[] x, int[] y) {
        assert x.length == y.length;

        int val = 0;
        int i;
        for (i = 0; i < x.length; i++)
            if (x[i] != y[i])
                val++;

        return val;
    }

    public void dot(Matrix w, Matrix h) {
        // makes the product inside this
        if (dimX != w.dimX)
            Matrix.perror("Matrix.class :: Dimension mismatch (this,W)");
        if (dimY != h.dimY)
            Matrix.perror("Matrix.class :: Dimension mismatch (this,H)");
        if (w.dimY != h.dimX)
            Matrix.perror("Matrix.class :: Dimension mismatch (W,H)");
        int i, j, k;
        double sum;

        for (i = 0; i < w.dimX; i++) {
            for (j = 0; j < h.dimY; j++) {
                sum = 0.0;
                for (k = 0; k < w.dimY; k++)
                    sum += (w.get(i, k) * h.get(k, j));
                set(i, j, sum);
            }
        }

        isVector = w.isVector;
        hasDefaultContent = false;
        if (dimX == dimY)
            isSquared = true;
        else
            isSquared = false;
        isSymmetric = false;
    }

    public static void assertSymmetric(Matrix m) {
        assertSquare(m);

        if (m.isSymmetric() == false)
            Matrix.perror("Matrix " + m.name + " is not symmetric...");
    }

    /**
     * 縦横の長さが同じことをチェック
     *
     * @param m
     */
    public static void assertSquare(Matrix m) {
        if (m.dimX != m.dimY)
            Matrix.perror("Matrix " + m.name + " is not square...");
    }

    public static void assertColumnStochastic(Matrix m) {
        int col = m.dimY;
        int lig = m.dimX;
        int i, j;
        double sum;
        boolean ok = true;
        String out = "";

        for (i = 0; i < col; i++) {
            sum = 0.0;
            for (j = 0; j < lig; j++) {
                if (m.get(j, i) < 0.0) {
                    ok = false;
                    out += "coordinate[" + i + "][" + j + "] <0. ";
                }
                sum += m.get(j, i);
            }
            if ((sum < 1.0 - MatrixConfig.Precision_For_Stochasticity) || (sum > 1.0 + MatrixConfig.Precision_For_Stochasticity)) {
                ok = false;
                out += "sum for col " + i + " = " + sum + " != 1. ";
            }
        }

        if (ok == false)
            Matrix.perror("Matrix " + m.name + " is not column stochastic: " + out + " (required precision: " + MatrixConfig.Precision_For_Stochasticity + " )");
    }

    public static void assertRowStochastic(Matrix m) {
        int col = m.dimY;
        int lig = m.dimX;
        int i, j;
        double sum;
        boolean ok = true;
        String out = "";

        for (i = 0; i < lig; i++) {
            sum = 0.0;
            for (j = 0; j < col; j++) {
                if (m.get(i, j) < 0.0) {
                    ok = false;
                    out += "coordinate[" + i + "][" + j + "] <0. ";
                }
                sum += m.get(i, j);
            }
            if ((sum < 1.0 - MatrixConfig.Precision_For_Stochasticity) || (sum > 1.0 + MatrixConfig.Precision_For_Stochasticity)) {
                ok = false;
                out += "sum for row " + i + " = " + sum + " != 1. ";
            }
        }

        if (ok == false)
            Matrix.perror("Matrix " + m.name + " is not row stochastic: " + out + " (required precision: " + MatrixConfig.Precision_For_Stochasticity + " )");
    }

    public static void assertDoublyStochastic(Matrix m) {
        assertRowStochastic(m);
        assertColumnStochastic(m);
    }

    public static void perror(String error_text) {
        JOptionPane.showMessageDialog(null, error_text, "Error", JOptionPane.WARNING_MESSAGE);
        System.out.println(error_text);
        return;
        //		System.out.println("...now exiting to system...\n");
        //		System.exit(1);
    }

    /**
     * クラスタリングに関係しているらしい
     *
     * @param W
     * @param H
     * @param r
     * @param c
     * @return
     */
    public static double multiplyComponent(Matrix W, Matrix H, int r, int c) {
        // returns (WH)_{rc}
        double v = 0.0;
        int i;
        for (i = 0; i < W.dimY; i++)
            v += (W.get(r, i) * H.get(i, c));
        return v;
    }

    /**
     * クラスタリングに関係しているらしい
     *
     * @param V
     * @param W
     * @param H
     * @param param
     * @param alpha
     * @return
     */
    public static double BregmanError(Matrix V, Matrix W, Matrix H, int param, int alpha) {
        // returns Bregman(V||WH)
        int i, j;
        Pnt p = new Pnt(1);
        Pnt q = new Pnt(1);
        double val = 0.0, dum;

        if (V.dimX != W.dimX)
            perror("Matrix.class :: Dimension mismatch (V,W)");
        if (V.dimY != H.dimY)
            perror("Matrix.class :: Dimension mismatch (V,H)");
        if (W.dimY != H.dimX)
            perror("Matrix.class :: Dimension mismatch (W,H)");

        for (i = 0; i < V.dimX; i++)
            for (j = 0; j < V.dimY; j++) {
                p.coordinates[0] = V.get(i, j);
                q.coordinates[0] = multiplyComponent(W, H, i, j);
                dum = Distortion.Bregman(p, q, param, alpha);
                val += dum;
            }

        return val;
    }

    public Matrix(String nm, int dX, int dY) {

        if ((dX == 0) || (dY == 0)) {
            Matrix.perror("Error(s) in Matrix Constructor");
        }

        name = nm;

        if (dX == 1)
            isVector = true;
        else
            isVector = false;

        hasDefaultContent = true;

        isSymmetric = true;

        if (dX == dY)
            isSquared = true;
        else {
            isSquared = false;
            isSymmetric = false;
        }

        dimX = dX;
        dimY = dY;
        mat = new GeneralMatrix(dimX, dimY);
        int i, j;

        // coordinates = new double[dimX][];
        // for (i = 0; i < dimX; i++)
        // coordinates[i] = new double[dimY];
        // rows

        // for (i = 0; i < dimX; i++)
        // for (j = 0; j < dimY; j++)
        // set(i, j, DefaultValue);
    }

    /**
     * クラスタリングに関係
     * <p/>
     * 実対称行列 A を受け取り,
     * A = LL^t と分解(コレスキー分解) して,
     * Aが正定値行列(固有値がすべて正)ならば trueを返す.
     *
     * @param a
     * @return
     */
    public static boolean test_choltmp(Matrix a) {
        if (a.dimX != a.dimY)
            return false;
        int i, j, k, n = a.dimX;
        double sum;
        Matrix test_choleski = new Matrix("Choleski decomposition of " + a.name, n, n);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                test_choleski.set(i, j, a.get(i, j));
        // aのコピー.

        for (i = 0; i < n; i++) {
            for (j = i; j < n; j++) {
                for (sum = test_choleski.get(i, j), k = i - 1; k >= 0; k--)
                    sum -= test_choleski.get(i, k) * test_choleski.get(j, k);
                if (i == j) {
                    if (sum <= 0.0)
                        return false;
                    test_choleski.set(i, i, Math.sqrt(sum));
                } else
                    test_choleski.set(j, i, sum / test_choleski.get(i, i));
            }
        }
        for (i = 0; i < n; i++)
            for (j = 0; j < i; j++)
                test_choleski.set(j, i, 0.0);
        test_choleski = null;
        return true;
    }

    /**
     * a.choleski_decomposition に,
     * aのコレスキー分解 L (a = LL^t)を入れる.
     * aは正定値行列 であること.
     *
     * @param a
     */
    public static void choltmp(Matrix a) {
        if (a.dimX != a.dimY)
            perror("Dimension mismatch");
        int i, j, k, n = a.dimX;
        double sum;
        a.choleski_decomposition = new Matrix("Choleski decomposition of " + a.name, n, n);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                a.choleski_decomposition.set(i, j, a.get(i, j));

        for (i = 0; i < n; i++) {
            for (j = i; j < n; j++) {
                for (sum = a.choleski_decomposition.get(i, j), k = i - 1; k >= 0; k--)
                    sum -= a.choleski_decomposition.get(i, k) * a.choleski_decomposition.get(j, k);
                if (i == j) {
                    if (sum <= 0.0)
                        perror(a.choleski_decomposition.name + " is not positive definite");
                    a.choleski_decomposition.set(i, i, Math.sqrt(sum));
                } else
                    a.choleski_decomposition.set(j, i, sum / a.choleski_decomposition.get(i, i));
            }
        }
        for (i = 0; i < n; i++)
            for (j = 0; j < i; j++)
                a.choleski_decomposition.set(j, i, 0.0);
    }

    /**
     * ******************************************************************************************************
     * Instance methods
     * ***
     */

    public String affiche(int n) {
        // prints a n x n submatrix
        int i, j, v;
        String val = "";

        if (n != -1)
            v = Math.min(n, Math.min(dimX, dimY));
        else
            v = Math.min(dimX, dimY);

        for (i = 0; i < v; i++) {
            for (j = 0; j < v; j++)
                val += Matrix.DF.format(get(i, j)) + " ";
            val += "\n";
        }
        return val;
    }

    public String toString() {
        int i, j;
        StringBuilder sb = new StringBuilder();
        sb.append("Matrix").append(name).append("\n");
        // String val = "Matrix " + name + "\n";

        for (i = 0; i < dimX; i++) {
            for (j = 0; j < dimY; j++) {
                if (j > 0) sb.append(' ');
                sb.append(String.format("%.03f", get(i, j)));
//		sb.append(Matrix.DF.format(get(i,j)));
            }
            // val += Matrix.DF.format(get(i, j)) + " ";
            sb.append('\n');
            // val += "\n";
        }
        return sb.toString();
    }

    /**
     * この行列 <- m
     *
     * @param m
     */
    public void copyOnThis(Matrix m) {
        if ((dimX != m.dimX) || (dimY != m.dimY))
            Matrix.perror("Dimension mismatch to copy matrices");

        int i, j;
        for (i = 0; i < dimX; i++)
            for (j = 0; j < dimY; j++)
                set(i, j, m.get(i, j));

        hasDefaultContent = m.hasDefaultContent;
    }

    public boolean isSymmetric() {
        return mat.isSymmetry();
    }

    public void computeSymmetric() {
        int i, j;
        double m;
        boolean found = false;
        if (dimX != dimY)
            isSymmetric = false;

        for (i = 0; i < dimX - 1; i++)
            for (j = i + 1; j < dimY; j++)
                if (Math.abs(get(i, j) - get(j, i)) > MatrixConfig.Precision_For_Eigensystems) {
                    isSymmetric = false;
                    found = true;
                    ;
                } else if (get(i, j) != get(j, i)) {
                    m = (get(i, j) + get(j, i)) / 2.0;
                    set(i, j, m);
                    set(j, i, m);
                }
        if (found == false)
            isSymmetric = true;
    }

    public synchronized void toD(Domain d, Matrix W) {// 対角行列　Dii = sum Wij.
        int dim;
        double sum;
        int i, j, nmade = 0, percent;

        if (dimX != dimY)
            Matrix.perror("Cannot generate D on a NON SQUARE matrix");

        dim = dimX;

        if ((dimX != W.dimX) || (dimY != W.dimY))
            Matrix.perror("Dimension mismatch between W and the number of selected Genes");

        if (AGCT.DoDebug)
            System.out.print("Generating matrix D (" + dim + " x " + dim + ")... ");
        JInformationFrame.getInstance().setTextProgressBar("computing D");

        for (i = 0; i < dim; i++) {
            sum = 0.0;

            //			sum = W.getRow(i).sum();
            for (j = 0; j < dim; j++) {
                percent = (int) (100.0 * ((double) nmade / (double) (dim * dim)));
                JInformationFrame.getInstance().setValueProgressBar(percent);

                sum += W.get(i, j);
                if (i != j)
                    set(i, j, 0.0);
                nmade++;
            }
            set(i, i, sum);
        }

        Matrix.assertSymmetric(this);

        if (AGCT.DoDebug)
            System.out.println("ok.");
        JInformationFrame.getInstance().setValueProgressBar(0);
        JInformationFrame.getInstance().setTextProgressBar(JInformationFrame.defaultProgressString);
    }

    public synchronized void toW(AGCT ag) {// 近さに基づきマトリックスを作る
        int dim = dimX, nmade = 0, percent;
        double kappa;
        if (dimX != dimY)
            Matrix.perror("Cannot generate W on a NON SQUARE matrix");

        if (dimX != ag.data.getMyDomain().numberSelectedGenes)
            Matrix.perror("Dimension mismatch between W and the number of selected Genes");

        if (AGCT.DoDebug)
            System.out.print("Generating matrix W (" + dim + " x " + dim + ") ... ");
        JInformationFrame.getInstance().setTextProgressBar("computing W");

        Gene[] genes = new Gene[dim];
        for (int i = 0; i < dim; i++) {
            genes[i] = (Gene) ag.data.getMyDomain().getGenes().get(ag.data.getMyDomain().selectedGeneNumberToGeneNumber[i]);
        }

        if (AGCT.Method_W == 0) {
            for (int i = 0; i < dim; i++) {
                percent = (int) (100.0 * ((double) nmade / (double) (dim * dim)));
                JInformationFrame.getInstance().setValueProgressBar(percent);

                for (int j = 0; j < dim; j++) {
                    set(i, j, Gene.getSimilarity(genes[i], genes[j]));
                    nmade++;
                }
            }
        } else if (AGCT.Method_W == 1) {//default
            AGCTCounter cc = new AGCTCounter("computing W", dim);

            for (int i = 0; i < dim; i++) {
                class A implements Comparable<A> {
                    double d;
                    int i;

                    A(double d, int i) {
                        this.d = d;
                        this.i = i;
                    }

                    public int compareTo(A o) {
                        return -(int) Math.signum(d - o.d);
                    }

                    ;

                    @Override
                    public String toString() {
                        return "A{" + "d=" + d + ", i=" + i + '}';
                    }
                }
                A[] as = new A[dim];
                for (int j = 0; j < dim; j++) {
                    as[j] = new A(Gene.getSimilarity(genes[i], genes[j]), j);
                    if (i == j) as[j] = new A(Double.POSITIVE_INFINITY, j);
                }
                Arrays.sort(as);

//		for(int j=0;j<dim/2;j++){
//		    double tmp=as[j].d;
//		    as[j].d=as[dim-1-j].d;
//		    as[dim-1-j].d=tmp;
//
//		    int tmp2=as[j].i;
//		    as[j].i=as[dim-1-j].i;
//		    as[dim-1-j].i=tmp2;
//		}

//		Debug.debug("Matrix#toW$as",as);

                for (int j = 0; j < AGCT.Number_Of_Neighbors; j++)
                    if (j + 1 < dim && as[j + 1].i < dimY) {
                        // avoid to Wii > 0
                        set(i, as[j + 1].i, 1.0);
                        set(as[j + 1].i, i, 1.0);
                    } else if (AGCT.MYDEBUG)
                        AGCT.debug("!?");
                cc.increment();
            }

            cc.end();
        }

        if (AGCT.Method_N == 3) {
            ag.myScaling.memorize(this);
            kappa = ag.myScaling.kappa();
            getScaling(kappa);
        }
        // getScaling();
        // last checks

        hasDefaultContent = false;
        computeSymmetric();
        assert isSymmetric();

        if (AGCT.DoDebug)
            System.out.println("ok.");
        JInformationFrame.getInstance().setValueProgressBar(0);
        JInformationFrame.getInstance().setTextProgressBar(JInformationFrame.defaultProgressString);
    }

    public synchronized void toC(Domain domain, boolean beforeSelection) {
        if (AGCT.MYDEBUG)
            AGCT.debug("Matrix.toC(domain,beforeSelection)");

        int dim = dimX, x, y, z;
        if (AGCT.MYDEBUG)
            AGCT.debug(toString());
        double slx, sly, avex, avey;
        Gene gg;
        if (dimX != dimY)
            Matrix.perror("Cannot generate C on a NON SQUARE matrix");

        if (dimX != AGCT.getInstance().data.getDimFeatures())
            Matrix.perror("Dimension mismatch");

        if (domain.getGenes() != null) {

            AGCTCounter cc = new AGCTCounter("computing " + name, (2 * dimX * domain.numberSelectedGenes) + (dimX * (dimX - 1) * domain.numberSelectedGenes / 2));

            if (AGCT.DoDebug)
                System.out.print("Generating matrix C (" + dim + " x " + dim + ") ... ");

            AGCT.getInstance().data.setAverage_features(new double[dimX]);
            AGCT.getInstance().data.setSigma_features(new double[dimX]);

            for (x = 0; x < dimX; x++) {
                //				AGCT.getInstance().data.getAverage_features()[x]=0.0;
                if (AGCT.MYDEBUG)
                    AGCT.debug("numberSelectedGenes:" + domain.numberSelectedGenes);
                for (y = 0; y < domain.numberSelectedGenes; y++) {
                    gg = domain.getSelectedGene(y); // (Gene)
                    // d.domainGenes.get(d.selectedGeneNumberToGeneNumber[y]);
                    AGCT.getInstance().data.getAverage_features()[x] += gg.getFinalCoordinates(x, beforeSelection);

                    cc.increment();
                }
            }

            for (x = 0; x < dimX; x++)
                AGCT.getInstance().data.getAverage_features()[x] /= (double) domain.numberSelectedGenes;

            for (x = 0; x < dimX; x++) {
                AGCT.getInstance().data.getSigma_features()[x] = 0.0;
                for (y = 0; y < domain.numberSelectedGenes; y++) {
                    gg = domain.getSelectedGene(y); // (Gene)
                    // d.domainGenes.get(d.selectedGeneNumberToGeneNumber[y]);
                    AGCT.getInstance().data.getSigma_features()[x] += ((gg.getFinalCoordinates(x, beforeSelection) - AGCT.getInstance().data.getAverage_features()[x]) * (gg
                            .getFinalCoordinates(x, beforeSelection) - AGCT.getInstance().data.getAverage_features()[x]));

                    cc.increment();
                }
            }

            for (x = 0; x < dimX; x++)
                AGCT.getInstance().data.getSigma_features()[x] = Math.sqrt(AGCT.getInstance().data.getSigma_features()[x]);

            for (x = 0; x < dimX; x++) {
                for (y = x; y < dimX; y++) {
                    set(x, y, 0);
                    // coordinates[x][y] = 0.0;
                    if (y != x)
                        set(y, x, 0);
                    // coordinates[y][x] = 0.0;
                    for (z = 0; z < domain.numberSelectedGenes; z++) {
                        gg = domain.getSelectedGene(z); // (Gene)
                        // d.domainGenes.get(d.selectedGeneNumberToGeneNumber[z]);
                        slx = gg.getFinalCoordinates(x, beforeSelection);
                        sly = gg.getFinalCoordinates(y, beforeSelection);
                        avex = AGCT.getInstance().data.getAverage_features()[x];
                        avey = AGCT.getInstance().data.getAverage_features()[y];

                        set(x, y, get(x, y) + ((slx - avex) * (sly - avey)));
                        // coordinates[x][y] += ((slx - avex) * (sly - avey));

                        cc.increment();
                    }
                    set(x, y, get(x, y) / (AGCT.getInstance().data.getSigma_features()[x] * AGCT.getInstance().data.getSigma_features()[y]));
                    // coordinates[x][y] /=
                    // (AGCT.getInstance().data.getSigma_features()[x] *
                    // AGCT.getInstance().data.getSigma_features()[y]);
                    if ((AGCT.getInstance().data.getSigma_features()[x] == 0.0) || (AGCT.getInstance().data.getSigma_features()[y] == 0.0)) {
                        if (x != y)
                            set(x, y, 0);
                            // coordinates[x][y] = 0.0;
                        else
                            set(x, y, domain.numberSelectedGenes);
                        // coordinates[x][y] = ((double)
                        // domain.numberSelectedGenes);
                    }
                    if (y != x)
                        set(y, x, get(x, y));
                    // coordinates[y][x] = coordinates[x][y];
                }
            }

            Matrix.assertSymmetric(this);
            cc.end();
        }

        if (AGCT.DoDebug)
            System.out.println("ok.");
        if (AGCT.MYDEBUG)
            AGCT.debug(toString());
    }

    public synchronized void toE(Domain d, Matrix C) {
        if (dimX != dimY)
            Matrix.perror("Cannot generate E on a NON SQUARE matrix");

        if (dimX != AGCT.getInstance().data.getDimFeatures())
            Matrix.perror("Dimension mismatch");

        AGCT.getInstance().data.setPca_eigenvalues(new double[dimX]);
        double[] e = new double[dimX];

        copyOnThis(C);
        Debug.debug("symm", isSymmetric());
        tred2(d, AGCT.getInstance().data.getPca_eigenvalues(), e);
        tqli(d, AGCT.getInstance().data.getPca_eigenvalues(), e);

        eigenSort(AGCT.getInstance().data.getPca_eigenvalues());
        Debug.debug(AGCT.getInstance().data.getPca_eigenvalues());
        transposeSquare();

        // checkOrthonormality(d,-1);
        // Matrix.checkEigensystem(d,cop,this,AGCT.getInstance().pca_eigenvalues,-1);
    }

    public synchronized void toV(Domain d) {
        double no, nol, ave, sig, ps, slo, cos, pco;
        int i, x, y;
        Gene gg;

        if (dimX != dimY)
            Matrix.perror("Cannot generate V on a NON SQUARE matrix");

        if (dimX != AGCT.getInstance().data.getDimFeatures())
            Matrix.perror("Dimension mismatch");

        if (d.getGenes() != null) {

            AGCTCounter cc = new AGCTCounter("computing " + name, dimX);

            for (i = 0; i < dimX; i++) {
                no = 0.0;
                for (x = 0; x < d.numberSelectedGenes; x++) {
                    gg = (Gene) d.getGenes().get(d.selectedGeneNumberToGeneNumber[x]);
                    pco = gg.pca_components.coordinates[i];
                    pco *= pco;
                    no += pco;
                }

                no = Math.sqrt(no);

                for (y = 0; y < dimX; y++) {
                    ave = AGCT.getInstance().data.getAverage_features()[y];
                    sig = AGCT.getInstance().data.getSigma_features()[y];

                    nol = 0.0;
                    ps = 0.0;

                    for (x = 0; x < d.numberSelectedGenes; x++) {
                        gg = (Gene) d.getGenes().get(d.selectedGeneNumberToGeneNumber[x]);
                        slo = gg.getFinalCoordinates(y, true);

                        if (sig > 0.0) {
                            nol += (((slo - ave) / sig) * ((slo - ave) / sig));
                            ps += ((slo - ave) / sig) * gg.pca_components.coordinates[i];
                        }
                    }

                    if (sig > 0.0) {
                        nol = Math.sqrt(nol);
                        cos = ps / (no * nol);
                    } else {
                        cos = 0.0;
                    }
                    set(y, i, cos);
                    // coordinates[y][i] = cos;
                }
                cc.increment();
            }
            cc.end();
        }
    }

    public synchronized void toVV(Domain d, Clustering cl) {
        if (d.getGenes() == null)
            Matrix.perror("No Genes !");

        AGCTCounter cc = new AGCTCounter("computing " + name, dimX * dimY * d.numberSelectedGenes);
        if (AGCT.DoDebug)
            System.out.print("Generating matrix VV (" + dimX + " x " + dimY + ") ... ");

        int i, j, k;
        Gene gg;

        double K = (double) dimX;
        double N = (double) d.numberSelectedGenes;
        Matrix soft = cl.myClusteringAlgorithm.getSoft_memberships();

        for (i = 0; i < dimX; i++)
            for (j = 0; j < dimY; j++)
                set(i, j, 0.0);

        for (i = 0; i < dimX; i++) {
            for (j = 0; j < dimY; j++) {
                for (k = 0; k < d.numberSelectedGenes; k++) {
                    gg = (Gene) d.getGenes().get(d.selectedGeneNumberToGeneNumber[k]);

                    set(i, j, get(i, j) + soft.get(i, k) * gg.pca_components.coordinates[j] * K / N);

                    // coordinates[i][j] += soft.get(i,k) *
                    // gg.pca_components.coordinates[j] * K / N;

                    cc.increment();
                }
            }
        }
        cc.end();
        if (AGCT.DoDebug)
            System.out.println("ok.");
    }

    public synchronized void toN(AGCT a, Matrix D, Matrix W) {
        int dim = dimX, nmade = 0, percent;

        if (dimX != dimY)
            Matrix.perror("Cannot generate N on a NON SQUARE matrix");

        if (dimX != a.data.getMyDomain().numberSelectedGenes)
            Matrix.perror("Dimension mismatch between N and the number of selected Genes");

        if (dimX != D.dimX)
            Matrix.perror("Dimension mismatch between N and D");

        if (dimX != W.dimX)
            Matrix.perror("Dimension mismatch between N and W");

        if (AGCT.DoDebug)
            System.out.print("Generating matrix N (" + dim + " x " + dim + ")... by Method " + AGCT.Method_N);
        JInformationFrame.getInstance().setTextProgressBar("computing N");
        double[] d = new double[D.rowCount()];
        for (int i = 0; i < d.length; i++) {
            d[i] = D.get(i, i);
        }

        int x, y, max = (dim * dim);
        if (AGCT.Method_N == 1)// default
            max *= 2;

        if (AGCT.Method_N == 1)
            a.data.setDWD(new Matrix("DWD", dim, dim));

        if ((AGCT.Method_N == 0) || (AGCT.Method_N == 1)) {
            for (x = 0; x < dim; x++) {
                for (y = 0; y < dim; y++) {
                    if (AGCT.Method_N == 0) {
                        set(x, y, W.get(x, y) / (Math.sqrt(d[x]) * Math.sqrt(d[y])));
                    } else if (AGCT.Method_N == 1) {
                        a.data.getDWD().set(x, y, W.get(x, y) / (Math.sqrt(d[x]) * Math.sqrt(d[y])));
                    }

                    percent = (int) (100.0 * ((double) nmade / (double) (max)));
                    JInformationFrame.getInstance().setValueProgressBar(percent);
                    nmade++;
                }
            }

            if (AGCT.Method_N == 1) {
                for (x = 0; x < dim; x++) {
                    for (y = 0; y < dim; y++) {
                        set(x, y, W.get(x, y) / d[x]);
                        percent = (int) (100.0 * ((double) nmade / (double) (max)));
                        JInformationFrame.getInstance().setValueProgressBar(percent);
                        nmade++;
                    }
                }
            }
        } else if ((AGCT.Method_N == 2) || (AGCT.Method_N == 3)) {
            for (x = 0; x < dim; x++) {
                for (y = 0; y < dim; y++) {
                    // coordinates[x][y] = W.coordinates[x][y];
                    set(x, y, W.get(x, y));
                    percent = (int) (100.0 * ((double) nmade / (double) (max)));
                    JInformationFrame.getInstance().setValueProgressBar(percent);
                    nmade++;
                }
            }
        }

        JInformationFrame.getInstance().setValueProgressBar(0);
        JInformationFrame.getInstance().setTextProgressBar(JInformationFrame.defaultProgressString);

        if ((AGCT.Method_N == 2) || (AGCT.Method_N == 3))
            fastDoublyStochasticApproximation();

        if (AGCT.Method_N == 0) {
            Matrix.assertSymmetric(this);
        } else if (AGCT.Method_N == 1) {
            Matrix.assertRowStochastic(this);
            Matrix.assertSymmetric(a.data.getDWD());
        } else if ((AGCT.Method_N == 2) || (AGCT.Method_N == 3)) {
            Matrix.assertDoublyStochastic(this);
        }
        Debug.debug("N", affiche(10));
        if (a.data.getDWD() == null) {
            Debug.debug("Matrix # toN, a.data.getDWD() = " + a.data.getDWD());
        } else {
            Debug.debug("Matrix # toN, a.data.getDWD() = ", a.data.getDWD().affiche(10));
        }

        if (AGCT.DoDebug)
            System.out.println("ok.");
    }

    /**
     * @param i   rowId
     * @param num - 大きい方からnum個をのこす
     */
    private void sparsifyRow(int i, int num) {
        double[] vs = new double[dimY - 1];
        for (int j = 0, k = 0; j < dimY; j++)
            if (i != j) {
                vs[k++] = get(i, j);
            }
        Arrays.sort(vs);
        num = Math.min(num, dimY - 1);
        double th = vs[dimY - 1 - num];//これ以上は残す
        for (int j = 0; j < dimY; j++)
            if (i != j) {
                if (get(i, j) < th) {
                    set(i, j, 0);
                }
            }
    }

    /**
     * @param d
     * @param clu
     * @param alpha
     * @param valSKK ii要素
     */
    public void toP(Domain d, Clustering clu, double alpha, double valSKK) {
        if (d.getGenes() == null)
            Matrix.perror("No Genes !");

        AGCTCounter cc = new AGCTCounter("computing " + name, dimX * dimY);
        if (AGCT.DoDebug)
            System.out.print("Generating matrix P (" + dimX + " x " + dimY + ") ... ");

        int n = d.numberSelectedGenes;
        if ((dimX != n) || (dimY != n))
            Matrix.perror("Dimension mismatch");

        for (int i = 0; i < dimX; i++)
            for (int j = 0; j < dimY; j++)
                set(i, j, 0);

        Pnt[] points = new Pnt[dimX];
        for (int i = 0; i < dimX; i++) {
            points[i] = clu.getRightPoint(i, clu.myClusteringAlgorithm.getBefore());
        }

        double[] phis = new double[n];
        Pnt[] grads = new Pnt[n];
        for (int i = 0; i < n; i++) {
            phis[i] = Distortion.phi_pnorm(points[i], alpha);
            grads[i] = Distortion.gradient_pnorm(points[i], alpha);
        }

        long st = System.currentTimeMillis();

        class E implements Comparable<E> {
            int i;
            double d;

            E(int i, double d) {
                this.i = i;
                this.d = d;
            }

            public int compareTo(E o) {
                return (int) Math.signum(d - o.d);
            }

            ;
        }
        ArrayList<E> list = new ArrayList<E>();
        for (int i = 0; i < dimX; i++) {
            for (int j = 0; j < dimY; j++) {
                if (i != j) {
                    //										set(i,j,-Distortion.Bregman(points[i],points[j],4,alpha));// 距離のマイナス→近いほど値が大きい,最も遠い→-inf
                    //下は、上のコメントアウトしたコードと同値
                    //					set(i,j,-(phis[i]-phis[j]-points[i].subtract(points[j]).dot(grads[j])));
                    list.add(new E(j, -(phis[i] - phis[j] - points[i].subtract(points[j]).dot(grads[j]))));
                } else {
                }
                cc.increment();
            }
            set(i, i, valSKK);
            //							sparsifyRow(i,MatrixConfig.MAX_KEEP_NUMBER_OF_AP);// TODO 独自に設定できるようにする
            Collections.sort(list);
            int keep = Math.min(list.size(), MatrixConfig.MAX_KEEP_NUMBER_OF_AP);
            // 大きい方からkeep個を残す
            int m = list.size();
            for (int j = 0; j < keep; j++) {
                set(i, list.get(m - 1 - j).i, list.get(m - 1 - j).d);
            }
            list.clear();
        }
        Debug.debug(System.currentTimeMillis() - st);

        Debug.debug("density of P", density());
        cc.end();
        if (AGCT.DoDebug)
            System.out.println("ok.");
    }

    public synchronized void getScaling(double kappa) {
        // Doubly stochastic scaling
        if (AGCT.DoDebug)
            System.out.println("Doubly Stochastic Scaling for Matrix " + name + " with value kappa = " + DF.format(kappa) + "...");
        int i, j;
        for (i = 0; i < dimX; i++)
            for (j = 0; j < dimY; j++)
                // coordinates[i][j] *= kappa;
                set(i, j, get(i, j) * kappa);
    }

    /**
     * matrixの非ゼロ要素のうち小さい　(sparsify_statistic * 100) % を除去する
     */
    public void sparsifyMatrix() {
        if (AGCT.DoDebug)
            System.out.print("Sparsifying matrix (" + (AGCT.Sparsify_Statistic * 100) + " % removed)... ");

        PriorityQueue<Double> que = new PriorityQueue<Double>(100, new Comparator<Double>() {
            @Override
            public int compare(Double arg0, Double arg1) {
                return arg1.compareTo(arg0);
            }

            ;
        });

        int ndiff = 0;
        for (int i = 0; i < dimX; i++) {
            for (int j = 0; j < dimY; j++) {
                if (get(i, j) != 0)
                    ndiff++;
            }
        }
        double nremove = AGCT.Sparsify_Statistic * (double) ndiff;
        for (int i = 0; i < dimX; i++)
            for (int j = 0; j < dimY; j++) {
                if (get(i, j) != 0)
                    que.offer(get(i, j));
                while (que.size() > nremove)
                    que.poll();
            }

        if (que.isEmpty())
            return;
        double threshold = que.peek();

        for (int i = 0; i < dimX; i++)
            for (int j = 0; j < dimY; j++) {
                if (get(i, j) < threshold)
                    set(i, j, 0.0);
            }
        if (AGCT.DoDebug)
            System.out.print("ok.\n");
    }

    /**
     * この行列自体を、この行列に「近い」二重確率行列に書き換える。
     * (二重確率行列：各行、列に関して、その列、行の数字を全て足すと１になるような行列。)
     * FIXME 対称行列は対称行列に変換されるはずだが、そうなっていない。
     */
    public synchronized void fastDoublyStochasticApproximation() {
        JInformationFrame.getInstance().setTextProgressBar("Fast DSA of " + name);

        if (AGCT.DoDebug)
            System.out.print("Computing the fast doubly stochastic approximation for system " + name + " (" + dimX + " x " + dimY + ")... ");
        double lastfrob = 0.0, lastmin = 0.0, allSum;
        int niter = 0, niterid = 0, i, j, n = dimX;

        //		int lastp=100+Inper;

        Matrix.assertSquare(this);

        boolean[] f_Stop = new boolean[1];
        double[] f_Sumrow = new double[n];
        double[] f_Sumcol = new double[n];
        double[] f_Sumgen = new double[1];
        double[] f_Frob = new double[1];
        double[] f_Min = new double[1];
        allSum = 0.0;

        f_Stop[0] = false;

        for (i = 0; i < n; i++) {
            f_Sumrow[i] = 0.0;
            f_Sumcol[i] = 0.0;
        }

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                f_Sumrow[i] += get(i, j);
                f_Sumcol[i] += get(j, i);
                allSum += get(i, j);
            }
        }

        f_Sumgen[0] = allSum;

        do {
            fastDsaP1P2(f_Stop, f_Sumrow, f_Sumcol, f_Sumgen, f_Frob, f_Min);

            // if ( (AGCT.DoDebug) && (niter % 10 == 0) )
            // System.out.print("\n. Frob = " + f_Frob[0] + " Min = " + f_Min[0]
            // + " (prec. = " + Precision_For_Doubly_Stochastic_Approximation +
            // ", niterid = " + niterid + " ");

            niter++;

            // if (niter%100 == 0)
            // System.out.println("Current niter = " + niter);

            if (niter == 1) {
                lastfrob = f_Frob[0];
                lastmin = f_Min[0];
                niterid = 0;
            } else {
                if ((lastfrob > f_Frob[0] - MatrixConfig.Precision_For_Doubly_Stochastic_Approximation)
                        && (lastfrob < f_Frob[0] + MatrixConfig.Precision_For_Doubly_Stochastic_Approximation)
                        && (lastmin > f_Min[0] - MatrixConfig.Precision_For_Doubly_Stochastic_Approximation)
                        && (lastmin < f_Min[0] + MatrixConfig.Precision_For_Doubly_Stochastic_Approximation))
                    niterid++;
                else {
                    if (lastfrob != f_Frob[0]) {
                        lastfrob = f_Frob[0];
                        niterid = 0;
                    }
                    if (lastmin != f_Min[0]) {
                        lastmin = f_Min[0];
                        niterid = 0;
                    }
                }
            }

            // System.out.println("Niterid = " + niterid + "; DS_Max = " +
            // Doubly_Stochastic_Iterations_Id_Max + ".");

            if (niterid == Doubly_Stochastic_Iterations_Id_Max)
                f_Stop[0] = true;

            JInformationFrame.getInstance().setValueProgressBar((int) (100.0 * (double) niterid / (double) Doubly_Stochastic_Iterations_Id_Max));
            JInformationFrame.getInstance().setTextProgressString(niterid + " out of " + Doubly_Stochastic_Iterations_Id_Max + " in plateau");

        } while (f_Stop[0] == false);

        JInformationFrame.getInstance().resetTextProgressString();
        JInformationFrame.getInstance().setValueProgressBar(0);
        JInformationFrame.getInstance().setTextProgressBar(JInformationFrame.defaultProgressString);
        if (AGCT.DoDebug)
            System.out.print(" [#Iterations : " + niter + "] ok.\n");
    }

    public synchronized void fastDsaP1P2(boolean[] f_Stop, double[] f_Sumrow, double[] f_Sumcol, double[] f_Sumgen, double[] f_Frob, double[] f_Min) {
        int n = f_Sumrow.length;
        assert f_Stop.length == 1;
        assert f_Sumrow.length == n;
        assert f_Sumcol.length == n;
        assert n == dimX;
        assert f_Sumgen.length == 1;
        assert f_Frob.length == 1;
        assert f_Min.length == 1;

        int i, j;
        double nextsumgen = 0.0;
        double v, num, den;

        f_Frob[0] = 0.0;
        f_Min[0] = 0.0;

        boolean allpos = true;

        double[] nextsumrow = new double[n];
        double[] nextsumcol = new double[n];

        for (i = 0; i < n; i++) {
            nextsumrow[i] = f_Sumrow[i];
            nextsumcol[i] = f_Sumcol[i];
        }
        nextsumgen = f_Sumgen[0];

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                num = (get(i, j) * n * n) + n + f_Sumgen[0] - ((f_Sumrow[i] + f_Sumcol[j]) * n);
                den = n * n;

                v = (num / den);

                f_Frob[0] += ((v - get(i, j)) * (v - get(i, j)));

                nextsumgen -= get(i, j);
                nextsumrow[i] -= get(i, j);
                nextsumcol[j] -= get(i, j);

                if (v < 0.0) {
                    allpos = false;
                    if (v < f_Min[0])
                        f_Min[0] = v;
                    set(i, j, 0.0);
                } else {
                    set(i, j, v);

                    nextsumgen += v;
                    nextsumrow[i] += v;
                    nextsumcol[j] += v;
                }
            }
        }

        if (allpos == true)
            f_Stop[0] = true;

        for (i = 0; i < n; i++) {
            f_Sumrow[i] = nextsumrow[i];
            f_Sumcol[i] = nextsumcol[i];
        }
        f_Sumgen[0] = nextsumgen;
    }

    public synchronized void tred2(Domain dd, double d[], double e[]) {
        MatrixMath.Tridiagonalize(this, d, e);
    }

    /**
     * 大事っぽい
     *
     * @param ddd
     * @param d
     * @param e
     */
    public synchronized void tqli(Domain ddd, double d[], double e[]) {
        MatrixMath.calcEigenVectorOfTridiagonalizedMatrix(this, d, e);
    }

    // O(n^2)
    public void eigenSort(double d[]) {
        MatrixMath.eigenSort(this, d);
    }

    public void transpose(Matrix w) {
        // puts the transpose in this
        // does not check for all other parameters
        if ((dimX != w.dimY) || (dimY != w.dimX))
            Matrix.perror("Matrix.class :: Dimension mismatch to transpose matrice");

        int i, j;
        for (i = 0; i < dimX; i++)
            for (j = 0; j < dimY; j++)
                set(i, j, w.get(j, i));

        if (dimX == 1)
            isVector = true;
        else
            isVector = false;
        hasDefaultContent = w.hasDefaultContent;
        if (dimX == dimY)
            isSquared = true;
        else
            isSquared = false;
        isSymmetric = w.isSymmetric;
    }

    public void transposeSquare() {
        // if (AGCT.DoDebug) System.out.print("Transposing square matrix " + name
        // + "... ");

        Matrix.assertSquare(this);

        int x, y, n = dimX;
        double val;
        for (x = 0; x < n; x++) {
            for (y = x + 1; y < n; y++) {
                val = get(x, y);
                set(x, y, get(y, x));
                set(y, x, val);
            }
        }
        // if (AGCT.DoDebug) System.out.print("ok.\n");
    }

    public synchronized void ludcmp(int[] indx, double[] d, boolean[] singular) {
        if ((dimX != dimY) || (dimX != indx.length))
            Matrix.perror("Dimension mismatch");
        if ((d.length != 1) || (singular.length != 1))
            Matrix.perror("Dimension of d is not 1");
        int i, imax = -1, j, k, n = dimX;
        double big, dum, sum, temp;
        double[] vv = new double[dimX];
        d[0] = 1.0;

        singular[0] = false;

        for (i = 1; i <= n; i++) {
            big = 0.0;
            for (j = 1; j <= n; j++)
                if ((temp = Math.abs(get(i - 1, j - 1))) > big)
                    big = temp;
            if (big == 0.0) {
                singular[0] = true;
                return;
            }
            vv[i - 1] = 1.0 / big;
        }
        for (j = 1; j <= n; j++) {
            for (i = 1; i < j; i++) {
                sum = get(i - 1, j - 1);
                for (k = 1; k < i; k++)
                    sum -= get(i - 1, k - 1) * get(k - 1, j - 1);
                set(i - 1, j - 1, sum);
            }
            big = 0.0;
            for (i = j; i <= n; i++) {
                sum = get(i - 1, j - 1);
                for (k = 1; k < j; k++)
                    sum -= get(i - 1, k - 1) * get(k - 1, j - 1);
                set(i - 1, j - 1, sum);
                if ((dum = vv[i - 1] * Math.abs(sum)) >= big) {
                    big = dum;
                    imax = i;
                }
            }
            if (j != imax) {
                for (k = 1; k <= n; k++) {
                    dum = get(imax - 1, k - 1);
                    set(imax - 1, k - 1, get(j - 1, k - 1));
                    set(j - 1, k - 1, dum);
                }
                d[0] = -d[0];
                vv[imax - 1] = vv[j - 1];
            }
            indx[j - 1] = imax;
            if (get(j - 1, j - 1) == 0.0)
                set(j - 1, j - 1, Matrix.Tiny);
            if (j != n) {
                dum = 1.0 / get(j - 1, j - 1);
                for (i = j + 1; i <= n; i++)
                    set(i - 1, j - 1, get(i - 1, j - 1) * dum);
                // coordinates[i - 1][j - 1] *= dum;
            }
        }
    }

    /**
     * Reference:
     * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.78.2759&rep=rep1&type=pdf
     * A Unifying Approach to Hard and Probabilistic Clustering
     * <p/>
     * CP factorization partを行っている.
     * <p/>
     * <p/>
     * look for CP factorization minimize || F - GGt ||^2.
     * ( F is DSSimMat, G is this(soft_memberships))
     * where G>=0,G is n*k matrix
     * where n is number of genes,k is number of clusters
     * <p/>
     * 2011/01/19 現在,Zass_Shashuaからしか呼ばれない.
     * 自身は,soft_memberships.
     * soft_memberships.get(clusterId,geneId) = probability of the gene is belong to the cluster.
     *
     * @param myDomain
     * @param F        doubly stochastic similarity matrix.
     */
    public void completePositiveFactorization(Domain myDomain, Matrix F) {
        int k = dimX, n = dimY, niter = 0;
        // height : number of clusters
        // width : number of genes
        double lastObjectiveFunction = 0.0;
        boolean stop = false;

        Matrix.assertSquare(F);

        AGCTCounter cc = new AGCTCounter("Pos. fact. of " + name, 100);

        if (F.dimX != n)
            Matrix.perror("Dimension mismatch for CP factorization.");

        if (AGCT.DoDebug)
            System.out.print("Computing the complete positive factorization of system " + name + " (" + dimX + " x " + dimY + ") ... ");

        do {
            for (int r = 0; r < k; r++) {
                for (int s = 0; s < n; s++) {
                    // clusterId,geneId はこの時点で固定.
                    double num = 0.0;
                    for (int i = 0; i < n; i++) {
                        if (i != s) {
                            num += get(r, i) * F.get(s, i);
                        }
                    }
                    num *= get(r, s);
                    // num = sum( soft_memberships.get(clusterId,geneId) * softmemberships.get(clusterId,geneId2) * sim(geneId1,geneId2) )

                    double den = 0.0;
                    for (int j = 0; j < k; j++) {
                        for (int i = 0; i < n; i++) {
                            if (i != s)
                                den += (get(j, s) * get(j, i) * get(r, i));
                        }
                    }
                    // den = sum( soft_memberships.get(clusterId2,geneId) * soft_memberships.get(clusterId2,geneId2) * soft_memberships.get(clusterId,geneId2).

                    if (den != 0.0) {
                        // 途中で書き換えているが良いのか？
                        set(r, s, num / den);
                    } else {
                        if (num == 0.0)
                            set(r, s, 0.0);
                        else {
                            Matrix.perror("CP-factorization: num = " + num + ", den = " + den);
                        }
                    }
                }
            }

            double objectiveFunction = 0.0;
            for (int geneId = 0; geneId < n; geneId++) {
                for (int geneId2 = 0; geneId2 < n; geneId2++) {
                    if (geneId2 != geneId) {
                        double sum = 0.0;
                        for (int clusterId = 0; clusterId < k; clusterId++) {
                            sum += (get(clusterId, geneId) * get(clusterId, geneId2));
                        }
                        objectiveFunction += ((F.get(geneId, geneId2) - sum) * (F.get(geneId, geneId2) - sum));
                    }
                }
            }

            if (niter == 0)
                lastObjectiveFunction = objectiveFunction;
            else {
                double delta = Math.abs(lastObjectiveFunction - objectiveFunction);

                int iq = (int) (100.0 * MatrixConfig.Precision_For_Doubly_Stochastic_Approximation / delta);
                if (iq > 100)
                    iq = 100;
                cc.setCounter(iq);

                if (delta < MatrixConfig.Precision_For_Doubly_Stochastic_Approximation)
                    stop = true;

                lastObjectiveFunction = objectiveFunction;
            }

            niter++;
        } while (stop == false);

        cc.end();
    }

    public void normalize() {
        int i, j, nlig = dimX, ncol = dimY;
        double sum = 0.0;

        for (j = 0; j < ncol; j++) {
            sum = 0.0;
            for (i = 0; i < nlig; i++) {
                if (get(i, j) < MatrixConfig.Precision_For_Eigensystems)
                    set(i, j, 0.0);
                sum += get(i, j);
            }

            if (sum == 0.0)
                Matrix.perror("The sum, " + sum + ", is not strictly positive");

            for (i = 0; i < nlig; i++)
                // coordinates[i][j] /= sum;
                set(i, j, get(i, j) / sum);
        }
    }

    /**
     * thisはcholtmpで,コレスキー分解されていること.
     * 固有値のlogをとったものを返す.
     *
     * @return
     */
    public double logdet() {
        double sum = 0.0;
        if (dimX != dimY)
            Matrix.perror("Dimension mismatch");
        for (int i = 0; i < dimX; i++)
            sum += Math.log(choleski_decomposition.get(i, i));
        return 2.0 * sum;
    }

    public void elsolve(double[] b, double[] y) {
        int i, j, n = dimX;
        double sum;
        if ((b.length != dimX) || (y.length != dimX) || (dimX != dimY))
            Matrix.perror("Dimension mismatch");
        for (i = 0; i < n; i++) {
            for (sum = b[i], j = 0; j < i; j++)
                sum -= choleski_decomposition.get(i, j) * y[j];
            y[i] = sum / choleski_decomposition.get(i, i);
        }
    }

    public boolean writeTo(String fileName) {
        try {
            FileWriter fw = new FileWriter(fileName);
            int m = dimX, n = dimY;
            fw.write(m + " " + n + "\n");
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    double val = get(i, j);
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

    @Override
    public double density() {
        int a = 0;
        for (int i = 0; i < dimX; i++) {
            for (int j = 0; j < dimY; j++) {
                if (get(i, j) != 0)
                    a++;
            }
        }
        return (double) a / (dimX * dimY);
    }

    public void show(int h, int w) {
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                System.out.print(get(i, j));
                if (j < w - 1) System.out.print(" ");
            }
            System.out.println();
        }
    }
}
