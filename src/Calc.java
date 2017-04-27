/*
* Dev. version
*/
package bequ.stat;

import java.math.BigDecimal;
import java.math.RoundingMode;

/**
 *
 * @author v.s.arnautov@yandex.ru 
 */
/**
  * <p>Calculations.</p>
  */
public class Calc {
 /**
  * <p>Returns rounded number.</p>
  *
  * @param d  number 
  * @param n  digits 
  * @return String 
  */
    public String StrRound (double d, int n) {
        return String.valueOf (new BigDecimal(d).setScale(n, RoundingMode.HALF_UP));
    }
    
  /**
  * <p>Returns square</p>
  *
  * @param d number (double) 
  * @return double 
  */
    public double sq (double d) {
        return d*d;
    }

    
}
