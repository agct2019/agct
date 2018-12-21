package test;

public class AGCTCounter {
    //	private JInformationFrame myJ;
    String myString;
    int counter;
    int max;

    public AGCTCounter(String s, int m) {
//		myJ = j;
        max = m;
        myString = s;
        counter = 0;
        JInformationFrame.getInstance().setTextProgressBar(s);
    }

    public void setText(String s) {
        JInformationFrame.getInstance().setTextProgressBar(s);
    }

    public void setCounter(int v) {
        int percent;
        counter = v;
        percent = (int) (100.0 * ((double) counter / (double) max));
        JInformationFrame.getInstance().setValueProgressBar(percent);
    }

    public void setBound(int niterid) {
        JInformationFrame.getInstance().setTextProgressString(niterid + " out of " + max + " in plateau");
    }

    public void increment() {
        int percent;
        counter++;
        percent = (int) (100.0 * ((double) counter / (double) max));
        JInformationFrame.getInstance().setValueProgressBar(percent);
    }

    public void setPercent(int p) {
        int q = p;
        if (q < 0)
            q = 0;
        if (q > 100)
            q = 100;

        JInformationFrame.getInstance().setValueProgressBar(q);
    }

    public void end() {
        JInformationFrame.getInstance().setValueProgressBar(0);
        JInformationFrame.getInstance().setTextProgressBar(JInformationFrame.defaultProgressString);
        JInformationFrame.getInstance().resetTextProgressString();
    }
}