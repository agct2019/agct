package _util;

/**
 * よく使う関数をstatic methodとしていれておく雑多なクラス.
 * 確実によく使う簡単目な関数だけを入れておくこと.
 * 適宜分類すること.
 *
 * @author okakeigo
 */
public class Util {
    private Util() {
    }

    /**
     * |x| means (arithmetic) mean of xs.<br>
     * n = xs.length<br>
     * standard deviation of xs = sqrt((sum (x-|x|)^2) / n).<br>
     *
     * @param xs - xs.length must be > 0.
     * @return standardDeviation
     */
    public static double standardDeviation(double[] xs) {
        assert xs.length > 0;
        int n = xs.length;
        double sum = 0;
        for (double x : xs) sum += x;
        double mean = sum / n;
        double res = 0;
        for (double x : xs) res += (x - mean) * (x - mean);
        return Math.sqrt(res / n);
    }
}
