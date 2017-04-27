/*
* Dev. version
*/
package bequ.stat;

/**
 * BE Exceptions
 * @author v.s.arnautov@yandex.ru 
 */
public class BEException extends Exception {
    
    public String description;
    
    BEException(Throwable e) { 

    }

    BEException() {
        
    }
    
    BEException(String str) {
        
    }


    BEException (int i) {
        if (i == 1) {
            description = "Fast calculation enabled!";
        } else description = "Uncnown error!";
   
    }
    
    BEException(int i, String str) {
        
    }
}
