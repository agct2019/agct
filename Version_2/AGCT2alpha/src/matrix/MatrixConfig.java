package matrix;

public class MatrixConfig {
    private MatrixConfig() {
    }

    public static final double Precision_For_Eigensystems = 10E-9;
    // precision required to compute eigensystems

    public static final double Precision_For_Cosine = 10E-6;
    // precision required to compute cosine

    public static final double Precision_For_Doubly_Stochastic_Approximation = Precision_For_Eigensystems / 100.0;
    // precision required to stop doubly stochastic approximation

    public static final double Precision_For_Stochasticity = Precision_For_Eigensystems * 100.0;

    public static int MAX_COL_OF_EIGENVECTOR = 50;

    public static int MAX_KEEP_NUMBER_OF_AP = 200;
}
