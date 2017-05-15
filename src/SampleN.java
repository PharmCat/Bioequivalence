/*
* Dev. version
*/
package bequ.stat;

/**
 * Based on Simple Power Calc
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
    

    public void calculate () throws BEException {
        
        double logCIh;
        double logCIl;
        double TRdif;
        double lCIh = Math.log(CIh);
        double lCIl = Math.log(CIl);
        double pc;
        Power powerCalc = new Power ();
        powerCalc.CV        = CV;
        powerCalc.alpha     = 0.05;
        powerCalc.diff      = -0.05;
        powerCalc.iteration = 1000000;
        powerCalc.CIh       = 1.25;
        powerCalc.CIl       = 0.8;
        
        
        if (ready) {  
            for (int n = 12; n < maxNum; n = n + 2){
                powerCalc.num = n;     
                powerCalc.apply_data();
                pc = powerCalc.calculate();
                if ( pc >= 0.8) {
                    num = n;
                    ip = pc;
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
