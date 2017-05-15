/*
* Dev. version
*/
package bequ.stat;

import cern.jet.random.AbstractContinousDistribution;
import java.util.Random;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import cern.jet.random.ChiSquare;
import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

/**
 *
 * @author v.s.arnautov@yandex.ru 
 */
public class Power {
    
    RandomEngine engine = new MersenneTwister();
    
    
    private final Random             r                  = new Random(); 
    private       TDistribution      tDistribution;
    private       NormalDistribution NormalDistribution = new NormalDistribution (); 
    private       double             SD = 123.0;
    private       double             MSE = 123.0;
    private       double             critVal;
    private       boolean            apply              = false;
    private       double             power;
    private       AbstractContinousDistribution    cdist;
    private       AbstractContinousDistribution    ndist;
    //PUBLIC
    public        int                num                = 24;
    public        double             CV                 = 0.30;
    public        double             alpha              = 0.05;
    public        double             diff               = 0.0;
    public        int                iteration          = 1000000;
    public        double             CIh                = 1.25;
    public        double             CIl                = 0.8;
   
  
    
    public boolean apply_data () {
        tDistribution   = new TDistribution(num - 2);
        critVal         = tDistribution.inverseCumulativeProbability(1 - alpha);
        cdist = new ChiSquare(num - 2,engine);
        ndist = new Normal(0.0, 1.0, engine);
        apply = true;
        return true;
    }
    
    public double calculate () throws BEException {
        
        double logCIh;
        double logCIl;
        double TRdif;
        double lCIh = Math.log(CIh);
        double lCIl = Math.log(CIl);
        power = 0.0;
        
        if (apply) {
            MSE = Math.log(CV*CV+ 1);
            double msq = Math.sqrt(0.5*(1/Math.floor((double) num/2) + 1/Math.ceil((double) num/2)));
            SD = Math.sqrt(MSE)*msq; 
            double lnDiff = Math.log(1.0 + diff);
            for (int i = 0; i < iteration; i++) {
                
                double SDVar  = Math.sqrt(MSE*cdist.nextDouble()/(num-2))*msq; 
                TRdif  = SD*ndist.nextDouble() + lnDiff; 
                logCIh = TRdif + SDVar*critVal;
                logCIl = TRdif - SDVar*critVal;
                if (logCIh < lCIh && logCIl > lCIl) power++;  
            }
            return power / (double) iteration;
            
        }
        
        return 0.0;
    }
    
}
