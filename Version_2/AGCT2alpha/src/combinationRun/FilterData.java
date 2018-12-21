package combinationRun;

public class FilterData {
    private int bitMask;
    private double percentile;

    public int getBitMask() {
        return bitMask;
    }

    public double getPercentile() {
        return percentile;
    }

    public FilterData(int bitMask, double percentile) {
        this.bitMask = bitMask;
        this.percentile = percentile;
        // TODO Auto-generated constructor stub
    }

    @Override
    public int hashCode() {
        // TODO Auto-generated method stub
        return Integer.valueOf(bitMask).hashCode() + Double.valueOf(percentile).hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        // TODO Auto-generated method stub
        FilterData o = (FilterData) obj;
        return bitMask == o.bitMask && percentile == o.percentile;
    }

}
