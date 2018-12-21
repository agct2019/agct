package test;

public class Distortion implements Debuggable {

    public static double MAX_DISTORT = 1E10;
    //For KL, IS

    //order of param =
    // 0 -> L22
    // 1 -> KL
    // 2 -> IS
    // 3 -> Amari(alpha)
    // 4 -> p-norm(alpha)

    public static String DIVERGENCES_NAME[] = {"L22", "Kullback-Leibler", "Itakura-Saito", "Amari (alpha)", "p-norm (p)"};

    public static double phi(Pnt p, int param, double alpha) {
        if (param == 0)
            return phi_l22(p);
        else if (param == 1)
            return phi_KL(p);
        else if (param == 2)
            return phi_IS(p);
        else if (param == 3)
            return phi_Amari(p, alpha);
        else if (param == 4)
            return phi_pnorm(p, alpha);
        else
            Matrix.perror("Distortion.class :: No param choice available for phi");
        return 0.0;
    }

    public static Pnt gradient(Pnt p, double param, double alpha) {
        if (param == 0)
            return gradient_l22(p);
        else if (param == 1)
            return gradient_KL(p);
        else if (param == 2)
            return gradient_IS(p);
        else if (param == 3)
            return gradient_Amari(p, alpha);
        else if (param == 4)
            return gradient_pnorm(p, alpha);
        else
            Matrix.perror("Distortion.class :: No param choice available for gradient");
        return null;
    }

    public static double phi_l22(Pnt p) {
        return p.dot(p);
    }

    public static double phi_KL(Pnt p) {
        int i;
        double val = 0.0;
        for (i = 0; i < p.coordinates.length; i++) {
            if (p.coordinates[i] < 0.0)
                Matrix.perror("Distortion.class :: Negative values in KL divergence");
            else if (p.coordinates[i] != 0.0)
                val += ((p.coordinates[i] * Math.log(p.coordinates[i])) - p.coordinates[i]);
        }
        return val;
    }

    public static double phi_IS(Pnt p) {
        int i;
        double val = 0.0;
        for (i = 0; i < p.coordinates.length; i++) {
            if (p.coordinates[i] < 0.0)
                Matrix.perror("Distortion.class :: Negative values in KL divergence");
            else if (p.coordinates[i] == 0.0) {
                val = MAX_DISTORT;
                break;
            } else
                val += (-Math.log(p.coordinates[i]));
        }
        return val;
    }

    public static double phi_Amari(Pnt p, double alpha) {
        if ((alpha < -1.0) || (alpha > 1.0))
            Matrix.perror("Distortion.class :: invalid values for alpha in Amari alpha divergence");

        int i;
        double val = 0.0, x, dum;
        for (i = 0; i < p.coordinates.length; i++) {
            if (p.coordinates[i] < 0.0)
                Matrix.perror("Distortion.class :: Negative values in KL divergence");

            x = p.coordinates[i];
            if (x != 0.0)
                dum = Math.pow(x, ((1.0 + alpha) / 2.0)); //Math.exp( ((1+alpha)/2.0) * Math.log(x));
            else
                dum = 0.0;

            val += 4.0 * (1 - dum) / (1 - (alpha * alpha));
        }
        return val;
    }

    /*
     * q = alpha/(alpha-1)
     * phi(xs) = ((sum |x_i|^q)^(1/q))^2 / 2
     */
    public static double phi_pnorm(Pnt p, double alpha) {
        assert alpha > 1.0 : "Distortion.class :: bad value for p in p-norm divergence";
        int i;
        double val = 0.0, dum, s = alpha / (alpha - 1.0);
        if (!p.all_null()) {
            for (i = 0; i < p.coordinates.length; i++) {
                dum = Math.abs(p.coordinates[i]);
                if (dum > 0.0)
                    dum = Math.pow(dum, s); //Math.exp(alpha*Math.log(dum));
                val += dum;
            }

            /***** Enleve: on peut la calculer meme lorsque la norme est nulle, mais il faut faire attention !
             if (val == 0.0){
             Matrix.perror("Distortion.class :: phi p-norm trouble");
             }
             *****/

            if (val != 0.0) {
                val = Math.pow(val, (1.0 / s)); //Math.exp( (1.0/alpha) * Math.log(val) );
                val = (val * val) / 2.0;
            }
        }
        return val;
    }

    /**
     * http://ibisforest.org/index.php?Bregman%E3%83%80%E3%82%A4%E3%83%90%E3%83%BC%E3%82%B8%E3%82%A7%E3%83%B3%E3%82%B9
     * <p/>
     * D(p,q) = phi(p) - qhi(q) + (p-q).dot(div(phi(q))).
     *
     * @param p
     * @param r
     * @param param num  param         phi(p)                        dist(p,q)
     *              // 0 -> L22           ||p||^2                       ||p-q||^2
     *              // 1 -> KL            sum((p_i log p_i) - p_i)      sum(p_i log(p_i/q_i))
     *              // 2 -> IS            -sum(log p_i)                 sum(p_i/q_i - log(p_i/q_i) -1)
     *              // 3 -> Amari(alpha)  sum(4(1-x^((1−α)/2))/(1-α^2))
     *              // 4 -> p-norm(alpha)
     * @param alpha
     * @return
     */
    public static double Bregman(Pnt p, Pnt r, int param, double alpha) {
        //		double breg4=-1;
        if (param == 4) {
            return Bregman4(p, r, alpha);
        }
        double phip = phi(p, param, alpha);
        double phir = phi(r, param, alpha);
        //		Debug.debug("phip",phip);
        //		Debug.debug("phir",phir);
        Pnt subs = p.subtract(r);
        //		Debug.debug("subs",subs);
        Pnt grad = gradient(r, param, alpha);
        //		Debug.debug("grad",grad);
        double bregdiv = phip - phir - subs.dot(grad);

        subs = grad = null;

        //		if(breg4!=-1){
        //			if(breg4!=bregdiv){
        //				Debug.debug(breg4,bregdiv);
        //				throw new RuntimeException();
        //			}
        //		}

        return bregdiv;
    }

    private static double Bregman4(Pnt p, Pnt r, double alpha) {
        double phip = phi_pnorm(p, alpha);
        double phir = phi_pnorm(r, alpha);
        //		Debug.debug("phip",phip);
        //		Debug.debug("phir",phir);
        Pnt subs = p.subtract(r);
        //		Debug.debug(subs);
        Pnt grad = gradient_pnorm(r, alpha);
        //		Debug.debug("grad",grad);
        double bregdiv = phip - phir - subs.dot(grad);
        subs = grad = null;
        return bregdiv;
    }

    public static double gradientX(JDSScaling j, double x) {
        double val = -1.0;
        if (j.indexSelected == 0)
            val = grad_l22(x);
        else if (j.indexSelected == 1)
            val = grad_KL(x);
        else if (j.indexSelected == 2)
            val = grad_IS(x);
        else if (j.indexSelected == 3)
            val = grad_Amari(x, 0.0);
        else
            Matrix.perror("Distortion.class :: No gradient applicable");
        return val;
    }

    public static double grad_l22(double x) {
        return 2.0 * x;
    }

    public static double grad_KL(double x) {
        if (x <= 0.0)
            Matrix.perror("Distortion.class :: non positive values for KL grad");
        return Math.log(x);
    }

    public static double grad_IS(double x) {
        if (x <= 0.0)
            Matrix.perror("Distortion.class :: non positive values for IS grad");
        return -1.0 / x;
    }

    public static double grad_Amari(double x, double alpha) {
        double fact = 0.0;
        if (x != 0.0)
            fact = (2.0 / (1.0 - alpha)) * Math.pow(x, ((alpha - 1.0) / 2.0));
        double pen = (4.0 / (1.0 - (alpha * alpha)));

        return pen - fact; //Math.exp( ((alpha-1.0)/2.0) * Math.log(x) );
    }

    public static Pnt gradient_l22(Pnt p) {
        Pnt g = new Pnt(p.coordinates.length);
        int i;
        for (i = 0; i < p.coordinates.length; i++)
            g.coordinates[i] = grad_l22(p.coordinates[i]);
        return g;
    }

    public static Pnt gradient_KL(Pnt p) {
        Pnt g = new Pnt(p.coordinates.length);
        int i;
        for (i = 0; i < p.coordinates.length; i++)
            g.coordinates[i] = grad_KL(p.coordinates[i]);
        return g;
    }

    public static Pnt gradient_IS(Pnt p) {
        Pnt g = new Pnt(p.coordinates.length);
        int i;
        for (i = 0; i < p.coordinates.length; i++)
            g.coordinates[i] = grad_IS(p.coordinates[i]);
        return g;
    }

    public static Pnt gradient_Amari(Pnt p, double alpha) {
        Pnt g = new Pnt(p.coordinates.length);
        int i;
        for (i = 0; i < p.coordinates.length; i++)
            g.coordinates[i] = grad_Amari(p.coordinates[i], alpha);
        return g;
    }

    /**
     * s = alpha/(alpha-1)
     * div phi(p)_i = (sum |x_i|^(2/s-1)) * |x_i|^(s-1) * sign(x_i)
     */
    public static Pnt gradient_pnorm(Pnt p, double alpha) {
        if (alpha <= 1.0)
            Matrix.perror("Distortion.class :: bad value for p in p-norm");

        Pnt g = new Pnt(p.coordinates.length);
        int i;
        double dum, val = 0.0, s = alpha / (alpha - 1.0), sig, x;

        if (!p.all_null()) {
            for (i = 0; i < p.coordinates.length; i++) {
                dum = Math.abs(p.coordinates[i]);
                if (dum > 0.0)
                    dum = Math.pow(dum, s); //Math.exp(q*Math.log(dum));
                val += dum;
            }
            if (val == 0.0)
                for (i = 0; i < p.coordinates.length; i++)
                    g.coordinates[i] = 0.0;
                //Matrix.perror("Distortion.class :: p-norm gradient trouble");
            else {
                // あとで割り算しているのでこれで良い
                val = Math.pow(val, ((s - 2.0) / s)); //Math.exp( ((q-2.0)/q) * Math.log(val) );

                for (i = 0; i < p.coordinates.length; i++) {
                    x = p.coordinates[i];
                    if (x < 0.0)
                        sig = -1.0;
                    else
                        sig = 1.0;
                    g.coordinates[i] = 0.0;
                    if (x != 0.0)
                        g.coordinates[i] = (sig * Math.pow(Math.abs(x), (s - 1.0))) / val; //(sig * Math.exp( (q-1.0) * Math.log(Math.abs(x)))) / val;
                }
            }
        } else {
            for (i = 0; i < p.coordinates.length; i++)
                g.coordinates[i] = 0.0;
        }

        return g;
    }

    public static double distortion_l22(Pnt p, Pnt q) {
        return Bregman(p, q, 0, -1);
    }

}