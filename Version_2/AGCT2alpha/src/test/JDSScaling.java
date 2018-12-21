package test;

class JDSScaling implements Debuggable {
    public static double EPSILON_DICHO = 1.0E-10;

    public static String divergences[] = {"L22", "Kullback-Leibler",
            "Itakura-Saito", "Amari (alpha = 0)"};

    public int indexSelected;
    double[] memoryLine;
    double[] memoryColumn;
    double memoryTot, binf, bsup;
    boolean memorized;

    AGCT myAGCT;

    JDSScaling(AGCT a) {
        myAGCT = a;
        indexSelected = 0;
        memorized = false;
        memoryLine = memoryColumn = null;
        memoryTot = 0.0;
    }

    public static String getString() {
        String val = "";
        int i;
        for (i = 0; i < divergences.length; i++)
            val += (i + " : " + divergences[i] + "\n");
        return val;
    }

    public static boolean isValid(int i) {
        if ((i >= 0) && (i < divergences.length))
            return true;
        return false;
    }

    public void flushAll() {
        memoryLine = memoryColumn = null;
        memoryTot = binf = bsup = 0.0;
        indexSelected = 0;
        memorized = false;
    }

    public void memorize(Matrix M) {
        int i, j;
        memoryLine = new double[M.rowCount()];
        memoryColumn = new double[M.colCount()];
        boolean binff = false;
        boolean bsupf = false;

        for (i = 0; i < M.rowCount(); i++)
            memoryLine[i] = 0.0;
        for (i = 0; i < M.colCount(); i++)
            memoryColumn[i] = 0.0;
        memoryTot = 0.0;

        for (i = 0; i < M.rowCount(); i++)
            for (j = 0; j < M.colCount(); j++) {
                memoryLine[i] += M.get(i, j);
                memoryColumn[j] += M.get(i, j);
                memoryTot += M.get(i, j);
            }

        for (i = 0; i < M.rowCount(); i++)
            if ((memoryLine[i] > 0.0)
                    && ((binf > 1.0 / memoryLine[i]) || (!binff))) {
                binf = 1.0 / memoryLine[i];
                binff = true;
            }

        for (j = 0; j < M.colCount(); j++)
            if ((memoryColumn[j] > 0.0)
                    && ((bsup < 1.0 / memoryColumn[j]) || (!bsupf))) {
                bsup = 1.0 / memoryColumn[j];
                bsupf = true;
            }

        // System.out.println("binf = " + binf + ", bsup = " + bsup);

        memorized = true;
    }

    public double lhs(double kappa) {
        if (!memorized)
            Matrix.perror("W not memorized");
        double val = 0.0;
        int i, j;
        for (i = 0; i < memoryLine.length; i++)
            val += (memoryLine[i] * Distortion.gradientX(this,
                    (kappa * memoryLine[i])));

        for (j = 0; j < memoryColumn.length; j++)
            val += (memoryColumn[j] * Distortion.gradientX(this,
                    (kappa * memoryColumn[j])));

        val -= (2.0 * memoryTot * Distortion.gradientX(this, 1.0));
        return val;
    }

    public double kappa() {
        if (!memorized)
            Matrix.perror("W not memorized");
        double bi = binf, bs = bsup, med, val, lastmed = 0.0, curmed = 0.0;
        boolean stop = false;
        int niter = 0;
        do {
            if (niter > 0)
                lastmed = curmed;

            med = (bi + bs) / 2.0;
            curmed = med;

            val = lhs(med);

            // System.out.println("bi = " + bi + ", med = " + med + ", bs = " +
            // bs + ", val = " + val);

            if (Math.abs(val) < EPSILON_DICHO)
                stop = true;
            else {
                if (val < 0.0)
                    bi = med;
                else
                    bs = med;
            }
            niter++;

            if ((niter > 0) && (lastmed == curmed))
                stop = true;

        } while (!stop);
        return med;
    }
}