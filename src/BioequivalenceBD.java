/*
* Dev. version
*/
package bequ.stat;

import java.util.Arrays;
import org.apache.commons.math3.distribution.TDistribution;
import org.nevec.rjm.BigDecimalMath;

/**
 *
 * @author v.s.arnautov@yandex.ru 
 */

public class BioequivalenceBD {
    
    private BigDecimalMath BDM;
    
    // Formulation T - 0; R - 1;
    //private         RandomEngine                     engine  = new MersenneTwister(new java.util.Date());
    //private         AbstractContinousDistribution    ndist   = new Normal(0.0, 1.0, engine);
    //private final Random  r           = new Random();                           // Генератор случайных чисел                      
    private boolean       completed   = false;                                  // Расчет выполнен
    private boolean       ready       = false;                                  // Данные загружены
    private boolean       noerror     = true;                                   // Нет ошибок
    
    //private boolean       simulation  = false;
    //private boolean       calculation = false;   // тип работы
    
    // Базовые параметры
    private BESubject[]    subjset;                                             // Данные субъектов
    private int            snum       = 0;                                      // Количество субъектов
    private int            pnum       = 2;                                      // Количество периодов
    private int            sqnum      = 2;                                      // Количество последовательностей
    private int            fnum       = 2;                                      // Количество препаратов
    private double         alpha      = 0.05;                                   // Альфа     
    private boolean        txtOutEn   = false;                                  // Текстовый вывод
    
  
    private int[][]        sq          = new int[sqnum][pnum];                  // Матрица [Последовательность][Период] = Препарат
    private String[]       forml       = new String[2];                         // Название препаратов
    private int[]          periodl     = new int [fnum];                        // Названия периодов
    private int[]          ssqn        = new int [sqnum];                       // Количество субъектов в последовательности
    private String         log         = "";                                    // Лог
    // Параметры симуляции
    //private double         cvIntra;                                             // Интра вариация
    //private double         cvInter     = 0.3;                                   // Интер вариация
    //private double         Fratio      = 0.0;
    //private double         Sqratio     = 0.01;
    //private double         Pratio      = 0.01;
    //private double         M           = 1.5;
    //private double         SD          = 1.5;
    
    // Результаты расчета
    private double        SSE         = 0.0;                                    // SS Error
    private double        MSE         = 0.0;                                    // MS Error
    private double        SE          = 0.0;                                    // SE
    private double        intraCV     = 0.0;                                    // Интра вариация
    
    private double[]      LSM         = new double[pnum];                       // LSM [препарат]
    private double        logCIl      = 0.0;                                    // Log Нижняя граница
    private double        logCIh      = 0.0;                                    // Log Верхняя граница
    private double        CIl         = 0.0;                                    // CI Нижняя граница
    private double        CIh         = 0.0;                                    // CI Верхняя граница
    private double        TRdif       = 0.0;                                    // TR Разность
    private double        TRratio     = 0.0;                                    // TR Отношение
    private double        SDF         = 0.0;                                    // Satterwhite DF
    private double        S_MSE       = 0.0;
    private double        S_CVintra   = 0.0;
    private String        txtOut      = "";
    
    private double        dCV         = 0.0;                                    //Отклонение от целевого CV
    
    //private int           errn        = 0;
    
    
    private double [][]   Mmatrix     = new double[sqnum][pnum];                // Матрица средних [Последовательность][Период]
    private int[][]       MmatrixN    = new int[sqnum][pnum];                   // Количество элементов в матрице
    
    private double[]       Dij;                                                 // Этап расчета вариации
    private double[]       Di         = new double[2];                          // Этап расчета вариации

    private TDistribution tDistribution;                                        // Распределение
    private double critVal;                                                     // Критическое значение
    
    
  /**
  * <p>Empty constructor.</p>
  */
    public BioequivalenceBD () {

 
    }
    
    public BioequivalenceBD (int i) {        
        snum          = i;
        tDistribution = new TDistribution(i - 2);
        critVal       = tDistribution.inverseCumulativeProbability(1 - alpha);
        matrix ();
 
    }
    
    
    
    private void matrix () {
        forml[0]      = "T";
        forml[1]      = "R";
        
        LSM[0]        = 0.0;
        LSM[1]        = 0.0;
        
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
   /* 
    public boolean refill () {
        if (!simulation) return false;
        completed     = false;
        ready         = false;
        ssqn[0]       = 0;
        ssqn[1]       = 0;
        erase ();
        fill (snum);
        return true;
    }
  */
    
  /**
  * <p>Load Data to BE calculation class</p>
  *
  * @param importsubj Subject Data Set Aray
  * @since 0.1
  */
    public void loadData (BESubject[] importsubj) {
        BDM = new BigDecimalMath ();
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
        
        Dij               = new double [snum];
        double [] Sseq    = new double [sqnum];                                  // Вариация для каждой группы
        double  preSS     = 0.0;
        
        
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
            svar  = (Dij[sbcnt] - Di[subjset[sbcnt].sequence])*(Dij[sbcnt] - Di[subjset[sbcnt].sequence]); 
            preSS += svar;
            Sseq[subjset[sbcnt].sequence]  += svar/ssqn[subjset[sbcnt].sequence];                                            //Group VAR    
        }
        
        for (int scnt = 0; scnt < sqnum; scnt++) {
            for (int pcnt = 0; pcnt < pnum; pcnt++) {
                Mmatrix[scnt][pcnt] = Mmatrix[scnt][pcnt]/MmatrixN[scnt][pcnt];
            }
            sdf1 += Sseq[scnt]/ssqn[scnt];
            sdf2 += Sseq[scnt]/ssqn[scnt]*Sseq[scnt]/ssqn[scnt]/(ssqn[scnt] - 1);
        }
        
        for (int scnt = 0; scnt < sqnum; scnt++) {
            for (int pcnt = 0; pcnt < pnum; pcnt++) { 
                if (sq[scnt][pcnt] == 0) LSM[0]  +=  Mmatrix[scnt][pcnt];   
                if (sq[scnt][pcnt] == 1) LSM[1]  +=  Mmatrix[scnt][pcnt];  
            }
        }
        
        SDF       = sdf1*sdf1/sdf2;
        
        LSM[0]    = LSM[0] / sqnum;
        LSM[1]    = LSM[1] / sqnum;
        
        TRdif     = LSM[0] - LSM[1];
        TRratio   = StrictMath.exp(TRdif);
        
        
        SSE       = preSS/2.0;
        MSE       = SSE/(snum - 2);
   
        S_MSE     = SSE/SDF;
          
        intraCV   = StrictMath.sqrt(StrictMath.exp(MSE) - 1);
        
        S_CVintra = StrictMath.sqrt(StrictMath.exp(S_MSE) - 1);
        
       
        SE        = StrictMath.sqrt(0.5*MSE*(1.0/ssqn[0] + 1.0/ssqn[1]));

        
        logCIh    = TRdif + SE*critVal;
        logCIl    = TRdif - SE*critVal;
        CIh       = StrictMath.exp(logCIh);
        CIl       = StrictMath.exp(logCIl);
        
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
  */
    public double getIntraCV ()  {
        return intraCV; 
    }
    
    public int get_snum ()  {
         return snum; 
    }
    
    public double get_SE ()  {
         return SE; 
    }
    
    public double get_MSE ()  {
        return MSE;
    }
    
    public double get_CIh ()  {
        return CIh; 
    }
    
    public double get_CIl ()  {
        return CIl; 
    }
    
    public double get_TRdif ()  {
        return TRdif; 
    }
    
    public double get_TRratio ()  {
        return TRratio;
    }
    
    public double get_SDF ()  {
        return SDF; 
    }
    
    public double get_S_CVintra ()  {
        return S_CVintra; 
    }
    
    
    // SET
    
    
    public void set_Alpha (Double a) {
        alpha         = a;
        critVal       = tDistribution.inverseCumulativeProbability(1 - alpha);
    }
    
    
    
    public boolean BE_pass () {
        
        if (completed){ 
            if (CIh < 1.25 && CIl > 0.8) {
                return true;
            } else return false;
        } else return false;
    
    
    }
    

    
    private void erase () {
        
            SSE         = 0.0;                                   
            MSE         = 0.0;                                    
            SE          = 0.0;                                   
            intraCV     = 0.0;                                    
            LSM         = new double[pnum];                     
            logCIl      = 0.0;                                    
            logCIh      = 0.0;                                    
            CIl         = 0.0;                                    
            CIh         = 0.0;                                   
            TRdif       = 0.0;                                    
            TRratio     = 0.0;                                    
            Mmatrix     = new double[sqnum][pnum];                
            MmatrixN    = new int[sqnum][pnum];                  
            Dij         = new double [snum];                                                 
            Di          = new double[2];                           
    }
    
    
}
