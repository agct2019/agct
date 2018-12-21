package gene;

import java.util.Arrays;

public class Feature {
    private final double[] activation;// correspond to each time.
    private final double[] initialActivation;// meanのみをとったactivation.

    public Feature(double[] activation, double[] initialActivation) {
        this.activation = activation;
        this.initialActivation = initialActivation;
    }

    public final double calcTotalVariation() {
        double res = 0;
        for (int i = 0; i < activation.length - 1; i++)
            res += Math.abs(activation[i + 1] - activation[i]);
        return res;
    }

    public final int calcPeakTimeId() {
        int res = 0;
        double val = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < activation.length; i++) {
            if (activation[i] > val) {
                val = activation[i];
                res = i;
            }
        }
        return res;
    }

    //	double getActivation(int timeId){
//		return activation[timeId];
//	}
    double getInitialActivation(int timeId) {
        return initialActivation[timeId];
    }

    private int hash = -1;

    @Override
    public int hashCode() {
        if (hash != -1)
            return hash;
        final int prime = 31;
        hash = 1;
        hash = prime * hash + Arrays.hashCode(activation);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        Feature other = (Feature) obj;
        if (!Arrays.equals(activation, other.activation))
            return false;
        return true;
    }

}
