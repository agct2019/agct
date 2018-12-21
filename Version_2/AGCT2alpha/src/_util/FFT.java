package _util;
//INF555
//Frank Nielsen, 2D FFT demo program
//
//Fourier transform: Glassman's algorithm

public class FFT {
    /**
     * 順方向のFFTを行い,同じ要素数の周波数成分を返す.
     *
     * @param Re -
     * @param Im -
     */
    static public void forwardFFT(double[] Re, double[] Im) {
        if (Re.length != Im.length)
            return;

        int n = Re.length;

        int a = 1;
        int b = n;
        int c = 1;

        boolean inu = true;

        double[] tempRe = new double[n];
        double[] tempIm = new double[n];

        while (b > 1) {
            a *= c;
            c = 2;
            while (b % c != 0)
                c++;
            b /= c;
            // call Glassman's Fast Fourier Transform routine
            if (inu) {
                primeFactorTransform(a, b, c, Re, Im, tempRe, tempIm);
            } else {
                primeFactorTransform(a, b, c, tempRe, tempIm, Re, Im);
            }
            inu = !inu;
        }

        // if odd number of calls to primeFactorTransform copy temp
        if (!inu) {
            for (int i = 0; i < n; i++) {
                Re[i] = tempRe[i];
                Im[i] = tempIm[i];
            }
        }
    }

    public static void forwardFFT(double[][] Re, double[][] Im) {
        // 2D FFT using 1D FFT

        int nCols = Re.length;
        int nRows = Re[0].length;

        // transform each of the columns
        for (int i = 0; i < nCols; i++)
            forwardFFT(Re[i], Im[i]);

        // transform each of the rows
        double[] rowRe = new double[nCols];
        double[] rowIm = new double[nCols];
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                rowRe[j] = Re[j][i];
                rowIm[j] = Im[j][i];
            }
            forwardFFT(rowRe, rowIm);
            for (int j = 0; j < nCols; j++) {
                Re[j][i] = rowRe[j];
                Im[j][i] = rowIm[j];
            }
        }
    }

    public static void inverseFFT(double[] Re, double[] Im) {
        int n = Re.length;
        for (int i = 0; i < n; i++)
            Im[i] = -Im[i];
        forwardFFT(Re, Im);
        for (int i = 0; i < Re.length; i++) {
            Re[i] = Re[i] / n;
            Im[i] = -Im[i] / n;
        }
    }

    public static void inverseFFT(double[][] Re, double[][] Im) {
        int nRows = Re.length;
        int nCols = Re[0].length;
        for (int i = 0; i < nRows; i++)
            inverseFFT(Re[i], Im[i]);

        double[] colRe = new double[nRows];
        double[] colIm = new double[nRows];
        for (int i = 0; i < nCols; i++) {
            for (int j = 0; j < nRows; j++) {
                colRe[j] = Re[j][i];
                colIm[j] = Im[j][i];
            }
            inverseFFT(colRe, colIm);
            for (int j = 0; j < nRows; j++) {
                Re[j][i] = colRe[j];
                Im[j][i] = colIm[j];
            }
        }
    }

    private static void primeFactorTransform(int a, int b, int c,
                                             double[] inRe, double[] inIm, double[] outRe, double[] outIm) {

        int inOffset = b * (c + 1) + 1;
        int outOffset = b * (a + 1) + 1;

        double angle = 6.28318530717958 / (double) (a * c);
        double omegaRe = 1.0;
        double omegaIm = 0.0;
        double deltaRe = (double) Math.cos(angle);
        double deltaIm = -(double) Math.sin(angle);
        double sumRe, sumIm;

        for (int ic = 1; ic <= c; ++ic) {
            for (int ia = 1; ia <= a; ++ia) {
                for (int ib = 1; ib <= b; ++ib) {
                    int index = ib + (c + ia * c) * b;
                    sumRe = inRe[index - inOffset];
                    sumIm = inIm[index - inOffset];
                    for (int jcr = 2; jcr <= c; ++jcr) {
                        int jc = c + 1 - jcr;
                        index = ib + (jc + ia * c) * b;
                        double tRe = inRe[index - inOffset] + omegaRe * sumRe - omegaIm * sumIm;
                        sumIm = inIm[index - inOffset] + omegaRe * sumIm + omegaIm * sumRe;
                        sumRe = tRe;
                    }
                    index = ib + (ia + ic * a) * b;
                    outRe[index - outOffset] = sumRe;
                    outIm[index - outOffset] = sumIm;
                }
                double tRe = deltaRe * omegaRe - deltaIm * omegaIm;
                omegaIm = deltaRe * omegaIm + deltaIm * omegaRe;
                omegaRe = tRe;
            }
        }
    }

    static void print(double[][] tab) {
        int i, j;

        for (i = 0; i < tab.length; i++) {

            for (j = 0; j < tab[i].length; j++) {
                System.out.print(tab[i][j] + "\t");
            }
            System.out.println("");
        }

    }

    public static void main(String[] args) {
        System.out.println("INF555 FFT demo program:");
        int d = 4;
        double[][] Re = new double[d][d];
        double[][] Im = new double[d][d];

        int i, j;

        for (i = 0; i < d; i++)
            for (j = 0; j < d; j++) {
                Re[i][j] = Math.random();
                Im[i][j] = 0.0;
            }

        print(Re);
        print(Im);
        System.out.println("------FFT-------");
        forwardFFT(Re, Im);
        print(Re);
        print(Im);

        System.out.println("-----Inverse---------");
        inverseFFT(Re, Im);
        print(Re);
        print(Im);
    }
}