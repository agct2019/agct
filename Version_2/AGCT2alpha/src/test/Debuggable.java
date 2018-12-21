package test;

import java.awt.*;
import java.text.DecimalFormat;

public interface Debuggable {

    static Color aquaMarine = new Color(127, 255, 212),
            bisque = new Color(255, 228, 196),
            brown = new Color(165, 42, 42),
            chartreuse = new Color(127, 255, 0),
            coral = new Color(255, 127, 80),
            darkGoldenRod = new Color(184, 134, 11),
            darkOliveGreen = new Color(85, 107, 47),
            darkSalmon = new Color(233, 150, 122),
            midnightBlue = new Color(25, 25, 112),
            peru = new Color(205, 133, 63);

    static DecimalFormat DF = new DecimalFormat("#0.000");

    ///////////////////////////////////////////////////////////////////////////////////////////

    static int Inper = 5;
    // for Verbose: prints the computations every inper% on numerical functions

    ///////////////////////////////////////////////////////////////////////////////////////////

    // For numerical methods that rely on convergence

    static int Step_Jacobi = 10000;
    // For verbose print on Jacobi rotations

    static int Doubly_Stochastic_Iterations_Id_Max = 100;
    // #times doubly stochastic iteration still runs while the stop bounds have been reached

    static int Affinity_Propagation_Id_Max = 20;
    // #times Affinity Propagation still runs while the exemplars do not change

    ///////////////////////////////////////////////////////////////////////////////////////////

    // precision required to check (row / col)-stochastic matrix entries

    ///////////////////////////////////////////////////////////////////////////////////////////

    static int Zoom_Choice = 1;
    // 0 : the diagonal of drag-drop rectangle is used to define the pos-neg Zoom
    // 1 : the rectangle of drag-drop represents the zoom
}