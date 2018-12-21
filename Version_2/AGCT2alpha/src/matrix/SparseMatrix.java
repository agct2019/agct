package matrix;

import java.util.HashMap;
import java.util.Map.Entry;

import static java.lang.Math.abs;

/**
 * 空間計算量 O(要素数)のMatrix.
 * 一次元的に情報を格納している。
 * setされていない場合0と考える
 *
 * @author okakeigo
 */
class SparseMatrix implements AbstMatrix {
    private final static double EPS = 1e-7; // problem specific.

    private static boolean eq(double a, double b) {
        return abs(a - b) < EPS;
    }

    private final int m;
    private final int n;

    private long toi(int x, int y) {
        return (((long) x) << 32) | y;
    }

    private int getX(long id) {
        return (int) (id >> 32);
    }

    private int getY(long id) {
        return (int) id;
    }

    private HashMap<Long, Double> map;
    private Vector[] rows;

    public SparseMatrix(int m, int n) {
        this.m = m;
        this.n = n;
        map = new HashMap<Long, Double>(m);
        rows = new SparseVector[m];
        for (int i = 0; i < m; i++) {
            rows[i] = new SparseVector(n);
        }
    }

    private boolean in(int x, int y) {
        return 0 <= x && x < m && 0 <= y && y < n;
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
    public int nonZeroCount() {
        return map.size();
    }

    public void set(int x, int y, double value) {
        assert in(x, y);
        rows[x].set(y, value);
        long id = toi(x, y);
        if (eq(value, 0)) {
            if (map.containsKey(id)) {
                map.remove(id);
            }
        } else {
            map.put(id, value);
        }
    }

    public double get(int x, int y) {
        assert in(x, y);
        long id = toi(x, y);
        if (map.containsKey(id))
            return map.get(id);
        else
            return 0;
    }

    /**
     * 時間計算量はO(各ベクトルの内積の計算量の和).
     */
    @Override
    public AbstMatrix mul(AbstMatrix mat) {
        assert colCount() == mat.rowCount();
        int m = rowCount(), n = mat.colCount();
        SparseMatrix res = new SparseMatrix(m, n);
        Vector[] vs = new Vector[m];
        Vector[] us = new Vector[n];
        for (int i = 0; i < m; i++) {
            vs[i] = getRow(i);
        }
        for (int i = 0; i < n; i++) {
            us[i] = mat.getCol(i);
        }

        for (int i = 0; i < m; i++) {
            //			Debug.debug(i);
            for (int j = 0; j < n; j++) {
                double value = vs[i].dot(us[j]);
                res.set(i, j, value);
            }
        }
        return res;
    }

    //TODO sparse性を利用するように変更！！
    @Override
    public Vector mul(Vector v) {
        //		Vector res = null;
        //		try {
        //			Class<?>[] types = { int.class };
        //			Constructor<? extends Vector> constructor = v.getClass().getConstructor(
        //					types);
        //			res = constructor.newInstance(n);
        //		} catch (Exception e) {
        //			throw new RuntimeException(e);
        //		}
        //		for (int i = 0; i < m; i++) {
        //			double val = 0;
        //			for (int j = 0; j < n; j++) {
        //				val += get(i, j) * v.get(j);
        //			}
        //			res.set(i, val);
        //		}
        //		return res;
        Vector res = new DenseVector(m);
        for (int i = 0; i < m; i++) {
            res.set(i, rows[i].dot(v));
        }
        return res;
    }

    /*
     * SparseVectorを返す。O(n)
     */
    @Override
    public Vector getRow(int x) {
        Vector res = new SparseVector(n);
        for (int y = 0; y < n; y++) {
            if (map.containsKey(toi(x, y)))
                res.set(y, get(x, y));
        }
        return res;
    }

    /**
     * SparseVectorを返す。O(m)
     */
    @Override
    public Vector getCol(int y) {
        Vector res = new SparseVector(m);
        for (int x = 0; x < m; x++) {
            if (map.containsKey(toi(x, y)))
                res.set(x, get(x, y));
        }
        return res;
    }

    @Override
    public AbstMatrix turn() {
        AbstMatrix res = new SparseMatrix(n, m);
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

    @Override
    public AbstMatrix myClone() {
        SparseMatrix res = new SparseMatrix(m, n);
        for (Entry<Long, Double> e : map.entrySet()) {
            res.set(getX(e.getKey()), getY(e.getKey()), e.getValue());
        }
        return res;
    }

    @Override
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
        return (double) map.size() / m / n;
    }
}
