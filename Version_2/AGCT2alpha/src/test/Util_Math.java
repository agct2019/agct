package test;

import normalizationData.NormalizationData;


/**
 * @author Takashi(Only comment I wrote, the true author is Natalia Polouliakh and Keigo Oka.)
 */
public class Util_Math {

    /**
     * @param vals            Readed data with expression level
     * @param normalizeOption Check T/F in the circadian's option.<- 誤解では？
     * @return
     */
    final public static double calcMean(double[] vals, NormalizationData.NormalizeMode normalizeOption) {
        switch (normalizeOption) {
            case Geometric_mean:
                double res = 1;
                for (double val : vals) {
                    res *= val;
                }
                return Math.pow(res, 1. / vals.length);
            case Arithmetic_mean:
                res = 0;
                for (double val : vals) {
                    res += val;
                }
                return res / vals.length;
            default:
                throw new IllegalArgumentException();
        }
    }

    final public static double calcStandardDivation(double[] vals) {
        double mean = 0;
        for (double val : vals) {
            mean += val;
        }
        mean /= vals.length;
        double res = 0;
        for (double val : vals) {
            res += Math.pow(Math.abs(val - mean), 2);
        }
        return Math.sqrt(res / vals.length);
    }
}
