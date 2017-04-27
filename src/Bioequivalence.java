/*
 * Licenced under Creative Commons «Attribution-NonCommercial-ShareAlike»
 * https://creativecommons.org/licenses/by-nc-sa/4.0/ 
 */
package bequ.stat;
import java.util.Arrays;
import org.apache.commons.math3.distribution.TDistribution;

/**
 *
 * @author v.s.arnautov@yandex.ru 
 */
public class Bioequivalence {
    
    
    // Formulation T - 0; R - 1;

    
    
    private TDistribution tDistribution;                                        // Распределение
    private double critVal;                                                     // Критическое значение
    private Calc          calc        = new Calc ();
    
    private boolean       completed   = false;                                  // Расчет выполнен
    private boolean       ready       = false;                                  // Данные загружены
    private boolean       fastCalc    = true;                                  // Расчет только ключевых параметров
    private boolean       noerror     = true;                                   // Нет ошибок
    private boolean       txtOutEn    = false;                                  // Текстовый вывод   
    private boolean       wsCalc      = false;                                  // Welch - Satterwhite

    //Служебные параметры и переменные
    
    // Базовые параметры
    private BESubject[]    subjset;                                             // Данные субъектов
    private int            snum        = 0;                                     // Количество субъектов
    private int            pnum        = 2;                                     // Количество периодов
    private int            sqnum       = 2;                                     // Количество последовательностей
    private int            fnum        = 2;                                     // Количество препаратов
    private double         alpha       = 0.05;                                  // Альфа     

    private double         lowRange    = 0.8;
    private double         highRange   = 1.25;
    private double         lnLowRange;
    private double         lnHighRange;
    
    
    
  
    private int[][]        sq          = new int[sqnum][pnum];                  // Матрица [Последовательность][Период] = Препарат
    private String[]       forml       = new String[fnum];                      // Название препаратов
    private int[]          periodl     = new int [fnum];                        // Названия периодов
    private int[]          ssqn        = new int [sqnum];                       // Количество субъектов в последовательности
    private double [][]    Mmatrix     = new double[sqnum][pnum];               // Матрица средних [Последовательность][Период]
    private int[][]        MmatrixN    = new int[sqnum][pnum];                  // Количество элементов в матрице средних
   // private double[]       Dij;                                                 // Этап расчета вариации
  //  private double[]       Di          = new double[2];                         // Этап расчета вариации

    private String         log         = "";                                    // Лог
    private String         txtOut      = "";    
    
    // Результаты расчета
    
    private double        SSE         = 0.0;                                    // SS Error
    private double        MSE         = 0.0;                                    // MS Error
    private double        SE          = 0.0;                                    // SE
    private double        intraCV     = 0.0;                                    // Интра вариация
    
    //private double[]      LSM         = new double[pnum];                       // LSM [препарат]
    private double        lnCIl       = 0.0;                                     // Ln Нижняя граница
    private double        lnCIh       = 0.0;                                     // Ln Верхняя граница
    private double        CIl         = 0.0;                                    // CI Нижняя граница
    private double        CIh         = 0.0;                                    // CI Верхняя граница
    private double        TRdif       = 0.0;                                    // TR Разность
    private double        TRratio     = 0.0;                                    // TR Отношение
    private double        SDF         = 0.0;                                    // Satterwhite DF
    private double        S_MSE       = 0.0;
    private double        S_CVintra   = 0.0;
   

    
    
  /**
  * <p>Empty constructor.</p>
  */
    public Bioequivalence () {
        
        lnLowRange    = StrictMath.log(lowRange);
        lnHighRange   = StrictMath.log(highRange);

 
    }
    
    public Bioequivalence (int i) {        
        snum          = i;
        tDistribution = new TDistribution(i - 2);
        critVal       = tDistribution.inverseCumulativeProbability(1 - alpha);
        lnLowRange    = StrictMath.log(lowRange);
        lnHighRange   = StrictMath.log(highRange);
        matrix ();
 
    }
    
    
    
    private void matrix () {
        forml[0]      = "T";
        forml[1]      = "R";

        
        periodl[0]    = 1;
        periodl[1]    = 2;
        
        // [Последовательность][Период] = Препарат
        sq[0][0]      = 0; //TR
        sq[0][1]      = 1;
        sq[1][0]      = 1; //RT
        sq[1][1]      = 0;
        
        for (int scnt = 0; scnt < sqnum; scnt++) {
            for (int pcnt = 0; pcnt < pnum; pcnt++) { 
                Mmatrix[scnt][pcnt]  = 0.0;
                MmatrixN[scnt][pcnt] = 0;
            }
        }
    }

    
  /**
  * <p>Load Data to BE calculation class</p>
  *
  * @param importsubj Subject Data Set Aray
  * @since 0.1
  */
    public void loadData (BESubject[] importsubj) {
        
        
        
         
            

        ssqn              = new int [sqnum];                                     // Количество субъектов в последовательности
        
        if (importsubj.length != snum){
            snum          = importsubj.length;
            subjset       = new BESubject[snum];
            tDistribution = new TDistribution(snum - 2);
            critVal       = tDistribution.inverseCumulativeProbability(1 - alpha);
        }
        
        for (int i = 0; i < snum; i++) {
            subjset[i]          = new BESubject (2);    
            subjset[i].id       = importsubj[i].id;
            subjset[i].sequence = importsubj[i].sequence;
            subjset[i].param    = Arrays.copyOf(importsubj[i].param, 2);
            if (importsubj[i].sequence == 0) ssqn[0]++;
            if (importsubj[i].sequence == 1) ssqn[1]++;
        }
        matrix ();
        ready         = true;
    
    
    
    }
   
  /**
  * <p>Calculate loaded data.</p>
  *
  * @return true if successful
  */
    public boolean calculate () {
        
        if (!ready) return false;
        
        int iVarf1          = 0; 
        int iVarf2          = 1; 
        double svar;
        double sdf1         = 0.0;
        double sdf2         = 0.0;
        //private double[]       Dij;                                                 // Этап расчета вариации
        double[] Dij                 = new double [snum];
        // private double[]       Di          = new double[2];                         // Этап расчета вариации
        double[] Di          = new double[2];  
        double[] Sseq       = new double [sqnum];                                  // Вариация для каждой группы
        double   preSS       = 0.0;
        //private double[]      LSM         = new double[pnum];                       // LSM [препарат]
        double[]      LSM         = new double[pnum];
        LSM[0]        = 0.0;
        LSM[1]        = 0.0;
        for (int scnt = 0; scnt < sqnum; scnt ++) {
            Di[scnt]   = 0.0;
            Sseq[scnt] = 0.0;
        }
        
        
        for (int scnt = 0; scnt < sqnum; scnt++) {
            for (int pcnt = 0; pcnt < pnum; pcnt++) {
                Mmatrix[scnt][pcnt]  = 0.0;
                MmatrixN[scnt][pcnt] = 0;
            }
        }
        
        if(txtOutEn) txtOut += "Subject DELTA calculation<br>";   
        
        for (int sbcnt = 0; sbcnt < snum; sbcnt++){
            
            for (int pcnt = 0; pcnt < 2; pcnt++) {
                Mmatrix[subjset[sbcnt].sequence][pcnt] += subjset[sbcnt].param[pcnt];
                MmatrixN[subjset[sbcnt].sequence][pcnt] ++;            
            }
            
            Dij[sbcnt] = get_param_by_form(subjset[sbcnt], iVarf1) - get_param_by_form(subjset[sbcnt], iVarf2);              // Delta calculate
            Di[subjset[sbcnt].sequence] += Dij[sbcnt]/ssqn[subjset[sbcnt].sequence];                                         // Mean of Delata for group
            
            if(txtOutEn) txtOut += "Subject" + String.valueOf(subjset[sbcnt].id) + " - " + String.valueOf(Dij[sbcnt]) + "<br>";  
        }
        
        for (int sbcnt = 0; sbcnt < snum; sbcnt++){
            svar  = calc.sq(Dij[sbcnt] - Di[subjset[sbcnt].sequence]); 
            preSS += svar;
            if (wsCalc){
                Sseq[subjset[sbcnt].sequence]  += svar/ssqn[subjset[sbcnt].sequence];                                            //Group VAR 
            }
        }
        
        for (int scnt = 0; scnt < sqnum; scnt++) {
            for (int pcnt = 0; pcnt < pnum; pcnt++) {
                Mmatrix[scnt][pcnt] = Mmatrix[scnt][pcnt]/MmatrixN[scnt][pcnt];
            }
            if (wsCalc){
                sdf1 += Sseq[scnt]/ssqn[scnt];
                sdf2 += calc.sq(Sseq[scnt]/ssqn[scnt])/(ssqn[scnt] - 1);
            }
        }
        
        for (int scnt = 0; scnt < sqnum; scnt++) {
            for (int pcnt = 0; pcnt < pnum; pcnt++) { 
                if (sq[scnt][pcnt] == 0) LSM[0]  +=  Mmatrix[scnt][pcnt];   
                if (sq[scnt][pcnt] == 1) LSM[1]  +=  Mmatrix[scnt][pcnt];  
            }
        }
        

        LSM[0]    = LSM[0] / sqnum;
        LSM[1]    = LSM[1] / sqnum; 
        TRdif     = LSM[0] - LSM[1];
        SSE       = preSS/2.0;
        MSE       = SSE/(snum - 2);
        SE        = StrictMath.sqrt(0.5*MSE*(1.0/ssqn[0] + 1.0/ssqn[1]));
        lnCIh    = TRdif + SE*critVal;
        lnCIl    = TRdif - SE*critVal;
        
        intraCV   = StrictMath.sqrt(StrictMath.exp(MSE) - 1);        
        
        
        if (!fastCalc) {
            TRratio   = StrictMath.exp(TRdif);
            CIh       = StrictMath.exp(lnCIh);
            CIl       = StrictMath.exp(lnCIl);
        }
        
        
        // Satterwhite 
        if (wsCalc){
        //SDF       = calc.sq(sdf1)/sdf2; 
        //S_MSE     = SSE/SDF;
        //S_CVintra = StrictMath.sqrt(StrictMath.exp(S_MSE) - 1);          
        }
        
        completed = true;
        
        return true;
    } 
    
    
    private double get_param_by_form(BESubject sub, int form) {
        
        for (int pcnt = 0; pcnt < pnum; pcnt++) {
            if (sq[sub.sequence][pcnt] == form) return sub.param[pcnt];
        }
        return 0.0;
        
    }
    
    // GET
    
  /**
  * <p>Returns completed status</p>
  *
  * @return true if completed
  */
    public boolean getCpl ()  {
        return completed;
    }
    
  /**
  * <p>Returns intra-subject CV.</p>
  *
  * @return true if completed
  * @throws bequ.stat.BEException
  */
    public double getIntraCV ()  throws BEException {
        if (!completed) throw new BEException ();
        return intraCV; 
    }
  /**
  * <p>Returns current subject number.</p>
  *
  * @return int snum 
  * @throws bequ.stat.BEException 
  */  
    public int getSnum () throws BEException {
        if (!completed) throw new BEException ();
        return snum; 
    }
    
  /**
  * <p>Returns SE.</p>
  *
  * @return double SE
  * @throws bequ.stat.BEException
  */     
    public double getSE ()  throws BEException {
        if (!completed) throw new BEException ();
        return SE; 
    }
    
  /**
  * <p>Returns MSE.</p>
  *
  * @return double MSE
  * @throws bequ.stat.BEException
  */
    
    public double getMSE ()  throws BEException {
        if (!completed) throw new BEException ();
        return MSE;
    }
  /**
  * <p>Returns High CI (%).</p>
  *
  * @return double Confidence Interval
  * @throws bequ.stat.BEException
  */    
    public double getHighCI () throws BEException  {
        if (!completed) throw new BEException ();
        if (fastCalc)  throw new BEException (1);
        return CIh*100.0; 
    }
  /**
  * <p>Returns Lower CI (%).</p>
  *
  * @return double Confidence Interval
  * @throws bequ.stat.BEException
  */      
    public double getLowCI () throws BEException  {
        if (!completed) throw new BEException ();
        if (fastCalc)  throw new BEException (1);
        return CIl*100.0; 
    }
    
      /**
  * <p>Returns LN High CI.</p>
  *
  * @return double Confidence Interval
  * @throws bequ.stat.BEException
  */    
    public double getLnHighCI () throws BEException  {
        if (!completed) throw new BEException ();
        return lnCIh; 
    }
  /**
  * <p>Returns Ln Lower CI.</p>
  *
  * @return double Confidence Interval
  * @throws bequ.stat.BEException
  */      
    public double getLnLowCI () throws BEException  {
        if (!completed) throw new BEException ();
        return lnCIl; 
    }
    
    
    
   /**
  * <p>Returns difference.</p>
  *
  * @return double difference
  * @throws bequ.stat.BEException
  */     
    public double getDif () throws BEException  {
        if (!completed) throw new BEException ();
        return TRdif; 
    }
 /**
  * <p>Returns ratio.</p>
  *
  * @return double ratio
  * @throws bequ.stat.BEException
  */     
    public double getRatio () throws BEException  {
        if (!completed) throw new BEException ();
        if (fastCalc)  throw new BEException (1);
        return TRratio;
    }
    
    
 /**
  * <p>Returns true if BE Passed.</p>
  *
  * @return boolean BE result
  * @throws bequ.stat.BEException
  */  
    
  
    public boolean getPassBE () throws BEException  {
        if (!completed) throw new BEException ();
        return lnCIh < lnHighRange && lnCIl > lnLowRange;
    }
   
   /*public boolean getPassBE () {
        
        if (completed){ 
            return CIh < 1.25 && CIl > 0.8;
        } else return false;
    
    }
*/
 /**
  * <p>Returns result BE object.</p>
  * <p>In project...</p>
  *
  * @return float BE result
  * @throws bequ.stat.BEException
  */ 
    public Object getResult () throws BEException {
        if (!completed) throw new BEException ();
        return new BEResult ();
    }
    // SET
    
    
    public void set_Alpha (Double a) {
        alpha         = a;
        critVal       = tDistribution.inverseCumulativeProbability(1 - alpha);
    }
    
    
}
