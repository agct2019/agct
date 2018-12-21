package matrix;

public interface AbstMatrix {
    int rowCount();

    int colCount();

    int nonZeroCount();

    double get(int x, int y);

    void set(int x, int y, double value);

    AbstMatrix mul(AbstMatrix mat);

    Vector mul(Vector v);

    Vector getRow(int x);

    Vector getCol(int y);

    /**
     * 非破壊的
     *
     * @return
     */
    AbstMatrix turn();

    AbstMatrix myClone();

    boolean isSymmetry();

    double density();
}
