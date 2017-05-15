/*
* Dev. version
*/
package bequ.stat;

/**
 *
 * @author v.s.arnautov@yandex.ru 
 */
import java.util.Random;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;






public class SampleN {
    
    private final Random             r                  = new Random(); 
    private       TDistribution      tDistribution;
    private       NormalDistribution NormalDistribution = new NormalDistribution (); 
    private       double             SE;
    private       double             critVal;
    private       double             ip;                                        //iteration power
    private       double             sn;                                        //Удачных итераций
    private       double             CV;
    private       int                num;
    
    private       boolean            ready              = false;
    private final int                maxNum             = 120;
   

    //PUBLIC
    
    private       double             power              = 0.8;
    private       double             alpha              = 0.05;
    private       double             diff               = 0.0;
    private       int                iteration          = 1000000;
    private       double             CIh                = 1.25;
    private       double             CIl                = 0.8;
   
    
    public SampleN (double cv) {
        
        CV = cv;
        ready = true;
    }
    
    public void set_cv (double cv) {
        CV = cv;
    }
    

    public void calculate () {
        
        double logCIh;
        double logCIl;
        double TRdif;
        double lCIh = Math.log(CIh);
        double lCIl = Math.log(CIl);

        
        if (ready) {  
            for (int n = 12; n < maxNum; n = n + 2){
                tDistribution       = new TDistribution(n - 2);
                NormalDistribution  = new NormalDistribution ();
                critVal         = tDistribution.inverseCumulativeProbability(1 - alpha);
                sn = 0.0;
                for (int i = 0; i < iteration; i++) {
                    SE = Math.sqrt(2*Math.log(Math.pow(CV, 2)+ 1)/n); 
                    TRdif  = SE*r.nextGaussian() + Math.log(1.0 + diff); 
                    logCIh = TRdif + SE*critVal;
                    logCIl = TRdif - SE*critVal;
                    if (logCIh < lCIh && logCIl > lCIl) sn++;  
                }
                ip = sn / (double) iteration;
                if (ip >= 0.8) {
                    num = n;
                    return;
                }
            }
        } 
    }
    
    public int get_num () {
        return num;
    }
    
    public double get_power () {
        return ip;
    }
    
}
