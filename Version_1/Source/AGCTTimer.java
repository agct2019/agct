import java.awt.event.ActionListener;

class AGCTTimer extends javax.swing.Timer {
    public static String GeneralTimerString = "generalTimer";
    public static String TriangulationTimerString = "triangulationTimer";
    public static String ScenarioTimerString = "scenarioTimer";

    public String name;

    AGCTTimer(int delay, ActionListener listener, String n){
	super(delay, listener);
	name = n;
    }
}