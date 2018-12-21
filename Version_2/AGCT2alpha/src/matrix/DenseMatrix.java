package matrix;

/**
 * 隣接行列表現で行列を格納する。
 *
 * @author okakeigo
 */
class DenseMatrix implements AbstMatrix {
    private double EPS = 1e-7;// problem specific.
    private int m, n;
    private double[][] mat;
    private int nonZero;

    DenseMatrix(int m, int n) {
        this.m = m;
        this.n = n;
        this.nonZero = 0;
        mat = new double[m][n];
    }

    @Override
    public Vector getCol(int y) {
        DenseVector res = new DenseVector(m);
        for (int i = 0; i < m; i++) {
            res.set(i, get(i, y));
        }
        return res;
    }

    @Override
    public Vector getRow(int x) {
        DenseVector res = new DenseVector(n);
        for (int i = 0; i < n; i++) {
            res.set(i, get(x, i));
        }
        return res;
    }

    @Override
    public Vector mul(Vector v) {
        Vector res = new DenseVector(m);
        for (int i = 0; i < m; i++) {
            res.set(i, v.dot(getRow(i)));
        }
        return res;
    }

    @Override
    public AbstMatrix myClone() {
        AbstMatrix res = new DenseMatrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res.set(i, j, get(i, j));
            }
        }
        return res;
    }

    @Override
    public int nonZeroCount() {
        return nonZero;
    }

    @Override
    public int rowCount() {
        return m;
    }

    @Override
    public int colCount() {
        return n;
    }

    @Override
    public double get(int x, int y) {
        return mat[x][y];
    }

    private boolean zero(double d) {
        return -EPS < d && d < EPS;
    }

    @Override
    public void set(int x, int y, double value) {
        if (zero(mat[x][y]) && !zero(value))
            nonZero++;
        else if (!zero(mat[x][y]) && zero(value))
            nonZero--;
        mat[x][y] = value;
    }

    @Override
    public AbstMatrix mul(AbstMatrix mat) {
        assert colCount() == mat.rowCount();
        int m = rowCount(), n = colCount(), l = mat.colCount();
        AbstMatrix res = new SparseMatrix(m, l);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < l; j++) {
                double value = 0;
                for (int k = 0; k < n; k++) {
                    value += get(i, k) * mat.get(k, j);
                }
                res.set(i, j, value);
            }
        }
        return res;
    }

    @Override
    public AbstMatrix turn() {
        DenseMatrix res = new DenseMatrix(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res.set(j, i, get(i, j));
            }
        }
        return res;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                sb.append(String.format("%7.3f", get(i, j)));
                if (j < n - 1)
                    sb.append(' ');
            }
            sb.append('\n');
        }
        return sb.toString();
    }

    public boolean isSymmetry() {
        int i, j;
        double d;
        if (m != n)
            return false;

        for (i = 0; i < m - 1; i++)
            for (j = i + 1; j < n; j++)
                if (Math.abs(get(i, j) - get(j, i)) > MatrixConfig.Precision_For_Eigensystems) {
                    return false;
                } else if (get(i, j) != get(j, i)) {
                    d = (get(i, j) + get(j, i)) / 2.0;
                    set(i, j, d);
                    set(j, i, d);
                }
        return true;
    }

    ;

    @Override
    public double density() {
        return 1;
    }

}
