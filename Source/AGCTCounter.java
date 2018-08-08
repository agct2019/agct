
class AGCTCounter{
    JInformationFrame myJ;
    String myString;
    int counter;
    int max;

    AGCTCounter(JInformationFrame j, String s, int m){
	myJ = j;
	max = m;
	myString = s;
	counter = 0;
	myJ.setTextProgressBar(s);
    }

    public void setText(String s){
	myJ.setTextProgressBar(s);
    }

    public void setCounter(int v){
	int percent;
	counter = v;
	percent = (int) (100.0 * ((double) counter / (double) max));
	myJ.setValueProgressBar(percent);
    }
 
    public void setBound(int niterid){
	myJ.setTextProgressString(niterid + " out of " + max + " in plateau");
    }
   
    public void increment(){
	int percent;
	counter ++;
	percent = (int) (100.0 * ((double) counter / (double) max));
	myJ.setValueProgressBar(percent);
    }

    public void setPercent(int p){
	int q = p;
	if (q < 0)
	    q = 0;
	if (q > 100)
	    q = 100;

	myJ.setValueProgressBar(q);
    }

    public void end(){
	myJ.setValueProgressBar(0);
	myJ.setTextProgressBar(JInformationFrame.defaultProgressString);
	myJ.resetTextProgressString();
    }
}