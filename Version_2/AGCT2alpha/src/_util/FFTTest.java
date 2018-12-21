package _util;

import junit.framework.Assert;
import org.junit.Test;

import static java.lang.Math.abs;
import static java.lang.Math.max;

public class FFTTest {
    private double[] polyMul(double[] as, double[] bs) {
        int n = as.length, m = bs.length;
        double[] res = new double[n + m - 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                res[i + j] += as[i] * bs[j];
            }
        }
        return res;
    }

    private double[] fastPolyMul(double[] as, double[] bs) {
        int N = 1;
        while (N < max(as.length, bs.length))
            N *= 2;
        N *= 2;
        double[] ra = new double[N], ia = new double[N];
        double[] rb = new double[N], ib = new double[N];
        for (int i = 0; i < as.length; i++)
            ra[i] = as[i];
        for (int i = 0; i < bs.length; i++)
            rb[i] = bs[i];
        FFT.forwardFFT(ra, ia);
        FFT.forwardFFT(rb, ib);

        double[] real = new double[N];
        double[] imag = new double[N];
        for (int i = 0; i < N; i++) {
            real[i] = ra[i] * rb[i] - ia[i] * ib[i];
            imag[i] = ra[i] * ib[i] + ia[i] * rb[i];
        }
        FFT.inverseFFT(real, imag);
//    for(int i=0;i<N;i++){
//      real[i]/=N;
//      imag[i]/=N;
//    }
        double[] res = new double[as.length + bs.length - 1];
        for (int i = 0; i < res.length; i++)
            res[i] = real[i];
        return res;
    }

    //  private void fft(int sign,double[] real,double[] imag){
//    int n=real.length,d=Integer.numberOfLeadingZeros(n)+1;
//    double theta=sign*2*PI/n;
//    for(int m=n;m>=2;m>>=1,theta*=2){
//      for(int i=0,mh=m>>1;i<mh;i++){
//        double wr=cos(i*theta),wi=sin(i*theta);
//        for(int j=i;j<n;j+=m){
//          int k=j+mh;
//          double xr=real[j]-real[k],xi=imag[j]-imag[k];
//          real[j]+=real[k];
//          imag[j]+=imag[k];
//          real[k]=wr*xr-wi*xi;
//          imag[k]=wr*xi+wi*xr;
//        }
//      }
//    }
//    for(int i=0;i<n;i++){
//      int j=Integer.reverse(i)>>>d;
//      if(j<i){
//        double tr=real[i];
//        real[i]=real[j];
//        real[j]=tr;
//        double ti=imag[i];
//        imag[i]=imag[j];
//        imag[j]=ti;
//      }
//    }
//  }
    private double EPS = 1e-8;

    private boolean eq(double[] as, double[] bs) {
        if (as.length != bs.length) return false;
        for (int i = 0; i < as.length; i++) {
            if (abs(as[i] - bs[i]) > EPS) return false;
        }
        return true;
    }


    @Test
    public void testFFT() {
        int n = 3000, m = 3000;
        double[] as = new double[n];
        double[] bs = new double[m];
        for (int i = 0; i < n; i++) {
            as[i] = Math.random();
        }
        for (int i = 0; i < m; i++) {
            bs[i] = Math.random();
        }
        Assert.assertTrue(eq(polyMul(as, bs), fastPolyMul(as, bs)));
    }
}
