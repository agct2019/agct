package _util;

import junit.framework.Assert;
import org.junit.Test;


public class UtilTest {
    /*
     * sample of http://en.wikipedia.org/wiki/Standard_deviation
     */
    double[] xs = new double[]{2, 4, 4, 4, 5, 5, 7, 9};

    @Test
    public void testStandardDeviation() {
        Assert.assertEquals(2, Util.standardDeviation(xs), 0);
    }
}
