/*
* Dev. version
*/
package bequ.stat;

/**
 * Public:
 * intraCV        - intra-subject CV
 * ratio          - ln (t/r) ratio
 * snum           - subject num
 * design         - design ("2X2"/""/"")
 * interSubjectCV - inter-subject simulation
 * interCV        - inter-subject CV
 * 
 * @author v.s.arnautov@yandex.ru 
 * @since 0.1
 */
public class GenerationParam {
    
    public double  intraCV;
    public double  ratio = 0.0;
    public int     snum;
    public String  design = "2X2";
    public boolean interSubjectCV = false;
    public double  interCV = 0.0;
    
}
