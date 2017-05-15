/*
* Dev. version
*/
package bequ.stat;



/**
 *<p>Power calculation by Full-subject Simulation</p>
 * @author FrozenCat
 * @since 0.1
 */
public class PowerFS {

    private       double             power;
    //PUBLIC
    public        int                num                = 24;
    public        double             CV                 = 0.30;
    public        double             alpha              = 0.05;
    public        double             diff               = 0.0;
    public        int                iteration          = 100000;
    public        double             CIh                = 1.25;
    public        double             CIl                = 0.8;
    private       Bioequivalence     dataset; 
    BEGeneration BEGen                                  = new BEGeneration ();
  
    public void apply_data () {}
    
    public double calculate () throws BEException {
        
        power   = 0.0; 
        dataset = new Bioequivalence ();
        for (int i = 0; i < iteration; i++) {
            dataset.loadData(BEGen.getData(num, CV, diff));
            dataset.calculate();
            if (dataset.getPassBE()) power++;
        }
        return power / (double) iteration;
    }
}
