import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

/**
 * 疎なベクトルに対応するデータ構造
 * 　hashmapを用いている。mapに要素が存在しない場合、０と考える。
 *
 * @author okakeigo
 */
 class SparseVector implements MyVector {
    private double EPS = 1e-7; // TODO 適切に設定
    private HashMap<Integer, Double> map;
    private int n;

    public SparseVector(int n) {
        this.n = n;
        map = new HashMap<Integer, Double>();
    }

    @Override
    public double get(int i) {
        if (map.containsKey(i))
            return map.get(i);
        return 0;
    }

    private boolean in(int i) {
        return 0 <= i && i < n;
    }

    @Override
    public void set(int i, double value) {
        assert in(i);
        if (Math.abs(value) < EPS) {
            if (map.containsKey(i))
                map.remove(i);
        } else
            map.put(i, value);
    }

    @Override
    public MyVector sub(MyVector v) {
        MyVector res = new SparseVector(n);
        for (int i = 0; i < n; i++) {
            res.set(i, get(i) - v.get(i));
        }
        return res;
    }

    /**
     * O(要素数)で計算できる。
     */
    @Override
    public double dot(MyVector v) {
        assert size() == v.size();
        double res = 0;
        for (Entry<Integer, Double> entry : map.entrySet()) {
            res += entry.getValue() * v.get(entry.getKey());
        }
        return res;
    }

    public MyVector mul(double d) {
        MyVector res = new SparseVector(n);
        for (int i = 0; i < n; i++) {
            res.set(i, get(i) * d);
        }
        return res;
    }

    @Override
    public MyVector div(double d) {
        MyVector res = new SparseVector(n);
        for (int i = 0; i < n; i++) {
            res.set(i, get(i) / d);
        }
        return res;
    }

    @Override
    public double norm() {
        return Math.sqrt(dot(this));
    }

    @Override
    public int size() {
        return n;
    }

    @Override
    public String toString() {
        double[] ds = new double[n];
        for (int i = 0; i < n; i++) {
            ds[i] = get(i);
        }
        return Arrays.toString(ds);
    }

    public static void main(String[] args) {
        MyVector v = new SparseVector(400000);
        for (int i = 0; i < 400000; i++) {
            if (i % 1 == 0)
                v.set(i, i);
        }
        System.out.println(v.dot(v));
    }
}
