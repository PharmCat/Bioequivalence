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
    private       Double             SD = 123.0;
    private       Double             MSE = 123.0;
    private       Double             critVal;
    private       boolean            apply              = false;
    private       Double             power;
    private       AbstractContinousDistribution    cdist;
    private       AbstractContinousDistribution    ndist;
    //PUBLIC
    public        int                num                = 24;
    public        Double             CV                 = 0.30;
    public        Double             alpha              = 0.05;
    public        Double             diff               = 0.0;
    public        int                iteration          = 1000000;
    public        Double             CIh                = 1.25;
    public        Double             CIl                = 0.8;
   
  
    
    public boolean apply_data () {
        tDistribution   = new TDistribution(num - 2);
        critVal         = tDistribution.inverseCumulativeProbability(1 - alpha);
        cdist = new ChiSquare(num - 2,engine);
        ndist = new Normal(0.0, 1.0, engine);
        apply = true;
        return true;
    }
    
    public double calculate () {
        
        Double logCIh;
        Double logCIl;
        Double TRdif;
        Double lCIh = Math.log(CIh);
        Double lCIl = Math.log(CIl);
        power = 0.0;
        
        if (apply) {
            
            for (int i = 0; i < iteration; i++) {
                MSE = Math.log(Math.pow(CV, 2)+ 1);
                double msq = Math.sqrt(0.5*(1/Math.floor((double) num/2) + 1/Math.ceil((double) num/2)));
                SD = Math.sqrt(MSE)*msq; 
                double SDVar  = Math.sqrt(MSE*cdist.nextDouble()/(num-2))*msq; 
                TRdif  = SD*ndist.nextDouble() + Math.log(1.0 + diff); 
                logCIh = TRdif + SDVar*critVal;
                logCIl = TRdif - SDVar*critVal;
                if (logCIh < lCIh && logCIl > lCIl) power++;  
            }
            return power / (double) iteration;
            
        }
        
        return 0.0;
    }
    
}
