package matrix;

//import utilities.Debug;

public class GeneralMatrix implements AbstMatrix {
    private AbstMatrix mat;
    private boolean sparse;
    private final int m, n;

    public GeneralMatrix(int m, int n) {
        this.m = m;
        this.n = n;
        mat = new SparseMatrix(m, n);
        sparse = true;
    }

    public int colCount() {
        return mat.colCount();
    }

    public double get(int x, int y) {
        return mat.get(x, y);
    }

    public Vector getCol(int y) {
        return mat.getCol(y);
    }

    public Vector getRow(int x) {
        return mat.getRow(x);
    }

    public AbstMatrix mul(AbstMatrix b) {
        AbstMatrix res = new GeneralMatrix(m, n);
        AbstMatrix a = mat.mul(b);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res.set(i, j, a.get(i, j));
            }
        }
        return res;
    }

    public Vector mul(Vector v) {
        return mat.mul(v);
    }

    public AbstMatrix myClone() {
        AbstMatrix res = new GeneralMatrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res.set(i, j, get(i, j));
            }
        }
        return res;
    }

    public int nonZeroCount() {
        return mat.nonZeroCount();
    }

    public int rowCount() {
        return mat.rowCount();
    }

    public void set(int x, int y, double value) {
        mat.set(x, y, value);
        if (sparse && mat.nonZeroCount() * 16 > m * n) {
            AbstMatrix nmat = new DenseMatrix(m, n);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    nmat.set(i, j, mat.get(i, j));
                }
            }
            mat = nmat;
            sparse = false;
        } else if (!sparse && mat.nonZeroCount() * 16 < m * n) {
            AbstMatrix nmat = new SparseMatrix(m, n);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    nmat.set(i, j, mat.get(i, j));
                }
            }
            mat = nmat;
            sparse = true;
        }
    }

    public AbstMatrix turn() {
        AbstMatrix res = new GeneralMatrix(m, n);
        AbstMatrix a = mat.turn();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res.set(i, j, a.get(i, j));
            }
        }
        return res;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                sb.append(String.format("%.3f ", get(i, j)));
            }
            sb.append('\n');
        }
        return sb.toString();
    }

    @Override
    public boolean isSymmetry() {
        return mat.isSymmetry();
    }

    @Override
    public double density() {
        return mat.density();
    }

}
