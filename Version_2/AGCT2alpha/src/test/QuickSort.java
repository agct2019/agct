package test;

import java.util.Vector;

public class QuickSort {
    private static long comparisons = 0;
    private static long exchanges = 0;

    /**
     * ********************************************************************
     * Quicksort code from Sedgewick 7.1, 7.2. (Ascending order)
     * *********************************************************************
     */
    public static void quicksort(double[] a) {
        shuffle(a); // to guard against worst-case
        quicksort(a, 0, a.length - 1);
    }

    public static void quicksort(double[] a, int[] ind) {
        shuffle(a, ind);
        quicksort(a, 0, a.length - 1, ind);
    }

    public static void quicksortString(Vector v) {
        shuffleString(v);
        quicksortString(v, 0, v.size() - 1);
    }

    public static void quicksortVectorDouble(Vector v, int i) {
        shuffleVector(v);
        quicksortVectorDouble(v, 0, v.size() - 1, i);
    }

    public static void quicksortVectorInteger(Vector v, int i) {
        shuffleVector(v);
        quicksortVectorInteger(v, 0, v.size() - 1, i);
    }

    public static void quicksort(double[] a, int left, int right) {
        if (right <= left)
            return;
        int i = partition(a, left, right);
        quicksort(a, left, i - 1);
        quicksort(a, i + 1, right);
    }

    public static void quicksortString(Vector v, int left, int right) {
        if (right <= left)
            return;
        int i = partitionString(v, left, right);
        quicksortString(v, left, i - 1);
        quicksortString(v, i + 1, right);
    }

    public static void quicksortVectorDouble(Vector v, int left, int right,
                                             int index) {
        if (right <= left)
            return;
        int i = partitionVectorDouble(v, left, right, index);
        quicksortVectorDouble(v, left, i - 1, index);
        quicksortVectorDouble(v, i + 1, right, index);
    }

    public static void quicksortVectorInteger(Vector v, int left, int right,
                                              int index) {
        if (right <= left)
            return;
        int i = partitionVectorInteger(v, left, right, index);
        quicksortVectorInteger(v, left, i - 1, index);
        quicksortVectorInteger(v, i + 1, right, index);
    }

    public static void quicksort(double[] a, int left, int right, int[] ind) {
        if (right <= left)
            return;
        int i = partition(a, left, right, ind);
        quicksort(a, left, i - 1, ind);
        quicksort(a, i + 1, right, ind);
    }

    private static int partition(double[] a, int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (less(a[++i], a[right]))
                // find item on left to swap
                ; // a[right] acts as sentinel
            while (less(a[right], a[--j]))
                // find item on right to swap
                if (j == left)
                    break; // don't go out-of-bounds
            if (i >= j)
                break; // check if pointers cross
            exch(a, i, j); // swap two elements into place
        }
        exch(a, i, right); // swap with partition element
        return i;
    }

    private static int partitionString(Vector v, int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (lessString((String) v.elementAt(++i), (String) v.elementAt(right)))
                // find item on left to swap
                ; // a[right] acts as sentinel
            while (lessString((String) v.elementAt(right), (String) v.elementAt(--j)))
                // find item on right to swap
                if (j == left)
                    break; // don't go out-of-bounds
            if (i >= j)
                break; // check if pointers cross
            exchString(v, i, j); // swap two elements into place
        }
        exchString(v, i, right); // swap with partition element
        return i;
    }

    private static int partitionVectorDouble(Vector v, int left, int right,
                                             int index) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (lessVectorDouble((Vector) v.elementAt(++i), (Vector) v
                    .elementAt(right), index))
                // find item on left to swap
                ; // a[right] acts as sentinel
            while (lessVectorDouble((Vector) v.elementAt(right), (Vector) v
                    .elementAt(--j), index))
                // find item on right to swap
                if (j == left)
                    break; // don't go out-of-bounds
            if (i >= j)
                break; // check if pointers cross
            exchVector(v, i, j); // swap two elements into place
        }
        exchVector(v, i, right); // swap with partition element
        return i;
    }

    private static int partitionVectorInteger(Vector v, int left, int right,
                                              int index) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (lessVectorInteger((Vector) v.elementAt(++i), (Vector) v
                    .elementAt(right), index))
                // find item on left to swap
                ; // a[right] acts as sentinel
            while (lessVectorInteger((Vector) v.elementAt(right), (Vector) v
                    .elementAt(--j), index))
                // find item on right to swap
                if (j == left)
                    break; // don't go out-of-bounds
            if (i >= j)
                break; // check if pointers cross
            exchVector(v, i, j); // swap two elements into place
        }
        exchVector(v, i, right); // swap with partition element
        return i;
    }

    private static int partition(double[] a, int left, int right, int[] ind) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (less(a[++i], a[right]))
                // find item on left to swap
                ; // a[right] acts as sentinel
            while (less(a[right], a[--j]))
                // find item on right to swap
                if (j == left)
                    break; // don't go out-of-bounds
            if (i >= j)
                break; // check if pointers cross
            exch(a, i, j, ind); // swap two elements into place
        }
        exch(a, i, right, ind); // swap with partition element
        return i;
    }

    // is x < y ?
    private static boolean less(double x, double y) {
        comparisons++;
        return (x < y);
    }

    // is x < y ?
    private static boolean lessString(String x, String y) {
        comparisons++;
        return ((x.compareTo(y)) < 0);
    }

    // is x < y ?
    private static boolean lessVectorDouble(Vector x, Vector y, int index) {
        comparisons++;
        double dx = ((Double) x.elementAt(index)).doubleValue();
        double dy = ((Double) y.elementAt(index)).doubleValue();
        return (dx < dy);
    }

    private static boolean lessVectorInteger(Vector x, Vector y, int index) {
        comparisons++;
        int dx = ((Integer) x.elementAt(index)).intValue();
        int dy = ((Integer) y.elementAt(index)).intValue();
        return (dx < dy);
    }

    private static void reverse(double[] a) {
        double t;
        int i;
        for (i = 0; i < (a.length / 2); i++) {
            t = a[i];
            a[i] = a[a.length - i];
            a[a.length - i] = t;
        }
    }

    // exchange a[i] and a[j]
    private static void exch(double[] a, int i, int j) {
        exchanges++;
        double swap = a[i];
        a[i] = a[j];
        a[j] = swap;
    }

    private static void exchString(Vector v, int i, int j) {
        exchanges++;
        String swap = (String) v.elementAt(i);
        v.setElementAt((String) v.elementAt(j), i);
        v.setElementAt(swap, j);
    }

    private static void exchVector(Vector v, int i, int j) {
        exchanges++;
        Vector swap = (Vector) v.elementAt(i);
        v.setElementAt((Vector) v.elementAt(j), i);
        v.setElementAt(swap, j);
    }

    private static void exch(double[] a, int i, int j, int[] ind) {
        exchanges++;
        double swap = a[i];
        int t = ind[i];
        a[i] = a[j];
        ind[i] = ind[j];
        a[j] = swap;
        ind[j] = t;
    }

    // shuffle the array a
    private static void shuffle(double[] a) {
        int N = a.length;
        for (int i = 0; i < N; i++) {
            int r = i + (int) (Math.random() * (N - i)); // between i and N-1
            exch(a, i, r);
        }
    }

    private static void shuffleString(Vector v) {
        int N = v.size();
        for (int i = 0; i < N; i++) {
            int r = i + (int) (Math.random() * (N - i)); // between i and N-1
            exchString(v, i, r);
        }
    }

    private static void shuffleVector(Vector v) {
        int N = v.size();
        for (int i = 0; i < N; i++) {
            int r = i + (int) (Math.random() * (N - i)); // between i and N-1
            exchVector(v, i, r);
        }
    }

    private static void shuffle(double[] a, int[] ind) {
        int N = a.length;
        for (int i = 0; i < N; i++) {
            int r = i + (int) (Math.random() * (N - i)); // between i and N-1
            exch(a, i, r, ind);
        }
    }

    // test client
    public static void main(String[] args) {
        int N = Integer.parseInt(args[0]);

        // generate N random real numbers between 0 and 1
        long start = System.currentTimeMillis();
        double[] a = new double[N];
        for (int i = 0; i < N; i++)
            a[i] = Math.random();
        long stop = System.currentTimeMillis();
        double elapsed = (stop - start) / 1000.0;
        System.out.println("Generating input:  " + elapsed + " seconds");

        // sort them
        start = System.currentTimeMillis();
        quicksort(a);
        stop = System.currentTimeMillis();
        elapsed = (stop - start) / 1000.0;
        System.out.println("Quicksort:   " + elapsed + " seconds");

        // print statistics
        System.out.println("Comparisons: " + comparisons);
        System.out.println("Exchanges:   " + exchanges);
    }
}
