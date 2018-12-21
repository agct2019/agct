package matrix;

public interface Vector {
    double get(int i);

    void set(int i, double value);

    int size();

    double dot(Vector v);

    double norm();

    Vector sub(Vector v);

    Vector mul(double d);

    Vector div(double d);
}
