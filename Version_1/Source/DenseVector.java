import java.util.Arrays;

 class DenseVector implements MyVector {
    private final int n;
    private double[] ds;

    public DenseVector(int n) {
        this.n = n;
        ds = new double[n];
    }

    private boolean in(int i) {
        return 0 <= i && i < n;
    }

    @Override
    public double get(int i) {
        assert in(i);
        return ds[i];
    }

    @Override
    public void set(int i, double value) {
        assert in(i);
        ds[i] = value;
    }

    @Override
    public int size() {
        return n;
    }

    @Override
    public double dot(MyVector v) {
        double res = 0;
        for (int i = 0; i < n; i++) {
            res += get(i) * v.get(i);
        }
        return res;
    }

    @Override
    public MyVector mul(double d) {
        MyVector res = new DenseVector(n);
        for (int i = 0; i < n; i++) {
            res.set(i, get(i) * d);
        }
        return res;
    }

    @Override
    public MyVector div(double d) {
        MyVector res = new DenseVector(n);
        for (int i = 0; i < n; i++) {
            res.set(i, get(i) / d);
        }
        return res;
    }

    @Override
    public double norm() {
        double res = 0;
        for (int i = 0; i < n; i++) {
            res += get(i) * get(i);
        }
        return Math.sqrt(res);
    }

    @Override
    public MyVector sub(MyVector v) {
        MyVector res = new DenseVector(n);
        for (int i = 0; i < n; i++) {
            res.set(i, get(i) - v.get(i));
        }
        return res;
    }

    @Override
    public String toString() {
        double[] ds = new double[n];
        for (int i = 0; i < n; i++) {
            ds[i] = get(i);
        }
        return Arrays.toString(ds);
    }


}
