public interface MyVector {
    double get(int i);

    void set(int i, double value);

    int size();

    double dot(MyVector v);

    double norm();

    MyVector sub(MyVector v);

    MyVector mul(double d);

    MyVector div(double d);
}
