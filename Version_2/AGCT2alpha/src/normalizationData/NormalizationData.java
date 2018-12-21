package normalizationData;


public class NormalizationData extends AbstractNormalizationData {
    public static enum NormalizeMode {
        Geometric_mean, Arithmetic_mean;
    }

    private double SDThreshold;
    private NormalizeMode normalizeMode;

    public void setSDThreshold(double SDThreshold) {
        this.SDThreshold = SDThreshold;
        nortifyObservers();
    }

    public double getSDThreshold() {
        return SDThreshold;
    }

    public void setNormalizeMode(NormalizeMode normalizeMode) {
        this.normalizeMode = normalizeMode;
        nortifyObservers();
    }

    public NormalizeMode getNormalizeMode() {
        return normalizeMode;
    }

    public NormalizationData() {
        setSDThreshold(2.0);
        setNormalizeMode(NormalizeMode.Geometric_mean);
    }
}
