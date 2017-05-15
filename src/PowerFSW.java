/*
* Dev. version
*/
package bequ.stat;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;



/**
 *<p>Power calculation by Full-subject Simulation</p>
 * @author FrozenCat
 * @since 0.1
 */
public class PowerFSW {

    private       double             power;
    //PUBLIC
    public        int                num                = 24;
    public        double             CV                 = 0.30;
    public        double             alpha              = 0.05;
    public        double             diff               = 0.0;
    public        int                iteration          = 1000000;
    public        double             CIh                = 1.25;
    public        double             CIl                = 0.8;
    private       Bioequivalence     dataset; 
    private       DescriptiveStatistics cvstat;
    BEGeneration BEGen                                  = new BEGeneration ();
  
    public void apply_data () {}
    
    public void calculate () throws BEException {
        
        power                                          = 0.0;
        cvstat                                         = new DescriptiveStatistics();
        dataset                                        = new Bioequivalence ();
        for (int i = 0; i < iteration; i++) {
        
            dataset.loadData(BEGen.getData(num, CV, diff));
            dataset.calculate();
            cvstat.addValue(dataset.getIntraCV());
            if (dataset.getPassBE()) power++;
        }
        power = power / (double) iteration;
    }
    
    
    
    public double getPower () {
        return power;
    }
    
    public double getCV () {
        return cvstat.getMean();
    }
}
