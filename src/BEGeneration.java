/*
* Dev. version
*/
package bequ.stat;

import cern.jet.random.AbstractContinousDistribution;
import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.MersenneTwister64;
import cern.jet.random.engine.RandomEngine;
import java.nio.ByteBuffer;
import java.security.SecureRandom;

/**
 * BE Generation methods
 * @author v.s.arnautov@yandex.ru 
 */
public class BEGeneration {
        
    private         RandomEngine                     engine;
    private         AbstractContinousDistribution    ndist;
    private         int[][]                          sq;  
    
        
    
    public BEGeneration () {
        engine  = new MersenneTwister(new java.util.Date());
        ndist   = new Normal(0.0, 1.0, engine);
        setMatrix ();
    
    }
    
    public void setMatrix (String type) {
        if (type.equals("2X2")){
            sq            = new int[2][2];
            sq[0][0]      = 0; //TR
            sq[0][1]      = 1; //TR
            sq[1][0]      = 1; //RT
            sq[1][1]      = 0; //RT
        } else {
            sq[0][0]      = 0; //TR
            sq[0][1]      = 1; //TR
            sq[1][0]      = 1; //RT
            sq[1][1]      = 0; //RT
        
        }
    }

    public void setMatrix () {
        sq            = new int[2][2];
        sq[0][0]      = 0; //TR
        sq[0][1]      = 1; //TR
        sq[1][0]      = 1; //RT
        sq[1][1]      = 0; //RT
    }
    
    private int getIntSeed() {
        SecureRandom sec = new SecureRandom();
        return ByteBuffer.wrap(sec.generateSeed(4)).getInt();
    }

  /**
  * <p>DataSet Generation</p>
  *
  * @param param subject num
  * @return BESubject[] - DataSet
  * @since 0.1
  * @see GenerationParam
  */
    public BESubject[] getData (GenerationParam param) {
       
    
        BESubject[] subjset           = new BESubject[param.snum];
        int     subjn     = 1;
        int     scnt      = 0;
        double  M         = 2;
        double  interVar  = 0;
        double  intraVar  = 0;
        double  SD        = StrictMath.sqrt(StrictMath.log(StrictMath.pow(param.intraCV, 2) + 1));
        
        // Субъект / Subject
        for (int sbcnt = 0; sbcnt < param.snum; sbcnt++){ 
            subjset[sbcnt]              = new BESubject (2);
            if (param.interSubjectCV) interVar = param.interCV*ndist.nextDouble();   // Интер вариация
            if (scnt == 0) {
                subjset[sbcnt].sequence = scnt;
                scnt = 1;
            } else {
                subjset[sbcnt].sequence = scnt;
                scnt = 0;
            }
            subjset[sbcnt].id = subjn;
            // Период / Period
            for (int pcnt = 0; pcnt < 2; pcnt++) {     
                intraVar = SD*ndist.nextDouble();  
                subjset[sbcnt].param[pcnt] = M + interVar + intraVar + (sq[subjset[sbcnt].sequence][pcnt] == 0 ? param.ratio : 0);   
            }   
            subjn++;
        }        
        return subjset;
    }
    
    
  /**
  * <p>DataSet Generation</p>
  *
  * @param i subject num
  * @param cv intersubject cv
  * @param ratio ln (ratio)
  * @return BESubject[] - DataSet
  * @since 0.1
  * 
  */
    public BESubject[] getData (int i, double cv, double ratio) {
        
        BESubject[] subjset           = new BESubject[i];
        int         subjn             = 1;
        int         scnt              = 0;
        double      M                 = 2;
        double      SD                = StrictMath.sqrt(StrictMath.log(cv*cv + 1));  
        // Субъект / Subject
        for (int sbcnt = 0; sbcnt < i; sbcnt++){ 
            subjset[sbcnt]              = new BESubject (2);  
            if (scnt == 0) {
                subjset[sbcnt].sequence = scnt;
                scnt = 1;
            } else {
                subjset[sbcnt].sequence = scnt;
                scnt = 0;
            }
            subjset[sbcnt].id = subjn;
            // Период / Period
            for (int pcnt = 0; pcnt < 2; pcnt++) {                          
                subjset[sbcnt].param[pcnt] = M + SD*ndist.nextDouble() + (sq[subjset[sbcnt].sequence][pcnt] == 0 ? ratio : 0);        
            }   
            subjn++;
        }        
        return subjset;
    }
    
}
