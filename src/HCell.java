import java.io.*;
import java.applet.Applet;
import java.awt.*;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
//jfree libs
import org.jfree.chart.*;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.general.DefaultPieDataset;
import org.jfree.data.xy.*;
import org.jfree.data.*;
import org.jfree.data.xy.XYDataset;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.jfree.data.xy.XYSeries;

class HCell{
   public static final double E=java.lang.Math.E;
   double time=30,h=0.0001,hh=0.1,t=0.0,tt=0.0;//time=sim time, h=time step, hh= output time step,t= current time,tt=output time
   double[] tData=new double[301];//time array
   double[][] YData=new double[15][301];//state variables array
   double[] Y=new double[15];
   double[] dY=new double[15];
// 1: four_AP_sensitive_currents_q_gate_q (dimensionless)
// 2: four_AP_sensitive_currents_r_gate_r (dimensionless)
// 3: hyperpolarisation_activated_current_y_gate_y (dimensionless)
// 4: L_type_Ca_channel_d_gate_d_L (dimensionless)
// 5: L_type_Ca_channel_f_gate_f_L (dimensionless)
// 6: membrane_V (millivolt)
// 7: rapid_delayed_rectifying_potassium_current_P_af_gate_P_af (dimensionless)
// 8: rapid_delayed_rectifying_potassium_current_P_as_gate_P_as (dimensionless)
// 9: rapid_delayed_rectifying_potassium_current_P_i_gate_P_i (dimensionless)
// 10: slow_delayed_rectifying_potassium_current_xs_gate_xs (dimensionless)
// 11: sodium_current_h_gate_h1 (dimensionless)
// 12: sodium_current_h_gate_h2 (dimensionless)
// 13: sodium_current_m_gate_m (dimensionless)
// 14: T_type_Ca_channel_d_gate_d_T (dimensionless)
// 15: T_type_Ca_channel_f_gate_f_T (dimensionless)

//initiation variables
    double calcium_background_current_g_b_Ca_Centre, calcium_background_current_g_b_Ca_Periphery,  four_AP_sensitive_currents_g_sus_Centre,four_AP_sensitive_currents_g_sus_Periphery;
    double four_AP_sensitive_currents_g_to_Centre,four_AP_sensitive_currents_g_to_Periphery,hyperpolarisation_activated_current_g_f_K_Centre;
    double hyperpolarisation_activated_current_g_f_K_Periphery,hyperpolarisation_activated_current_g_f_Na_Centre,hyperpolarisation_activated_current_g_f_Na_Periphery;
    double ionic_concentrations_Ca_i,ionic_concentrations_Ca_o,ionic_concentrations_K_i,ionic_concentrations_K_o,ionic_concentrations_Na_i,ionic_concentrations_Na_o;
    double L_type_Ca_channel_E_Ca_L,L_type_Ca_channel_g_Ca_L_Centre,L_type_Ca_channel_g_Ca_L_Periphery,membrane_CmCentre,membrane_CmPeriphery,membrane_dCell;
    double membrane_F,membrane_R,membrane_T,persistent_calcium_current_i_Ca_p_max_Centre,persistent_calcium_current_i_Ca_p_max_Periphery;
    double potassium_background_current_g_b_K_Centre,potassium_background_current_g_b_K_Periphery,rapid_delayed_rectifying_potassium_current_g_K_r_Centre;
    double rapid_delayed_rectifying_potassium_current_g_K_r_Periphery,rapid_delayed_rectifying_potassium_current_P_i_gate_tau_P_i;
    double slow_delayed_rectifying_potassium_current_g_K_s_Centre,slow_delayed_rectifying_potassium_current_g_K_s_Periphery,sodium_background_current_g_b_Na_Centre;
    double sodium_background_current_g_b_Na_Periphery,sodium_calcium_exchanger_d_NaCa,sodium_calcium_exchanger_gamma_NaCa,sodium_calcium_exchanger_k_NaCa_Centre;
    double sodium_calcium_exchanger_k_NaCa_Periphery,sodium_current_g_Na_Centre,sodium_current_g_Na_Periphery,sodium_potassium_pump_i_p_max_Centre,sodium_potassium_pump_i_p_max_Periphery;
    double sodium_potassium_pump_K_m_K,sodium_potassium_pump_K_m_Na,T_type_Ca_channel_E_Ca_T,T_type_Ca_channel_g_Ca_T_Centre,T_type_Ca_channel_g_Ca_T_Periphery;
   
 //-------------------------------------------------------------------------------
// Computed variables
//-------------------------------------------------------------------------------

// calcium_background_current_g_b_Ca   (microS)
// calcium_background_current_i_b_Ca   (nanoA)
// four_AP_sensitive_currents_g_sus   (microS)
// four_AP_sensitive_currents_g_to   (microS)
// four_AP_sensitive_currents_i_sus   (nanoA)
// four_AP_sensitive_currents_i_to   (nanoA)
// four_AP_sensitive_currents_q_gate_q_infinity   (dimensionless)
// four_AP_sensitive_currents_q_gate_tau_q   (second)
// four_AP_sensitive_currents_r_gate_r_infinity   (dimensionless)
// four_AP_sensitive_currents_r_gate_tau_r   (second)
// hyperpolarisation_activated_current_g_f_K   (microS)
// hyperpolarisation_activated_current_g_f_Na   (microS)
// hyperpolarisation_activated_current_i_f_K   (nanoA)
// hyperpolarisation_activated_current_i_f_Na   (nanoA)
// hyperpolarisation_activated_current_y_gate_alpha_y   (per_second)
// hyperpolarisation_activated_current_y_gate_beta_y   (per_second)
// L_type_Ca_channel_d_gate_alpha_d_L   (per_second)
// L_type_Ca_channel_d_gate_beta_d_L   (per_second)
// L_type_Ca_channel_d_gate_d_L_infinity   (dimensionless)
// L_type_Ca_channel_d_gate_tau_d_L   (second)
// L_type_Ca_channel_f_gate_alpha_f_L   (per_second)
// L_type_Ca_channel_f_gate_beta_f_L   (per_second)
// L_type_Ca_channel_f_gate_f_L_infinity   (dimensionless)
// L_type_Ca_channel_f_gate_tau_f_L   (second)
// L_type_Ca_channel_g_Ca_L   (microS)
// L_type_Ca_channel_i_Ca_L   (nanoA)
// membrane_Cm   (microF)
// membrane_FCell   (dimensionless)
// persistent_calcium_current_i_Ca_p   (nanoA)
// persistent_calcium_current_i_Ca_p_max   (nanoA)
// potassium_background_current_g_b_K   (microS)
// potassium_background_current_i_b_K   (nanoA)
// rapid_delayed_rectifying_potassium_current_g_K_r   (microS)
// rapid_delayed_rectifying_potassium_current_i_K_r   (nanoA)
// rapid_delayed_rectifying_potassium_current_P_a   (dimensionless)
// rapid_delayed_rectifying_potassium_current_P_af_gate_P_af_infinity   (dimensionless)
// rapid_delayed_rectifying_potassium_current_P_af_gate_tau_P_af   (second)
// rapid_delayed_rectifying_potassium_current_P_as_gate_P_as_infinity   (dimensionless)
// rapid_delayed_rectifying_potassium_current_P_as_gate_tau_P_as   (second)
// rapid_delayed_rectifying_potassium_current_P_i_gate_P_i_infinity   (dimensionless)
// reversal_and_equilibrium_potentials_E_Ca   (millivolt)
// reversal_and_equilibrium_potentials_E_K   (millivolt)
// reversal_and_equilibrium_potentials_E_K_s   (millivolt)
// reversal_and_equilibrium_potentials_E_Na   (millivolt)
// slow_delayed_rectifying_potassium_current_g_K_s   (microS)
// slow_delayed_rectifying_potassium_current_i_K_s   (nanoA)
// slow_delayed_rectifying_potassium_current_xs_gate_alpha_xs   (per_second)
// slow_delayed_rectifying_potassium_current_xs_gate_beta_xs   (per_second)
// sodium_background_current_g_b_Na   (microS)
// sodium_background_current_i_b_Na   (nanoA)
// sodium_calcium_exchanger_i_NaCa   (nanoA)
// sodium_calcium_exchanger_k_NaCa   (nanoA)
// sodium_current_g_Na   (microS)
// sodium_current_h_gate_F_Na   (dimensionless)
// sodium_current_h_gate_h   (dimensionless)
// sodium_current_h_gate_h1_infinity   (dimensionless)
// sodium_current_h_gate_h2_infinity   (dimensionless)
// sodium_current_h_gate_tau_h1   (second)
// sodium_current_h_gate_tau_h2   (second)
// sodium_current_i_Na   (nanoA)
// sodium_current_m_gate_m_infinity   (dimensionless)
// sodium_current_m_gate_tau_m   (second)
// sodium_potassium_pump_i_p   (nanoA)
// sodium_potassium_pump_i_p_max   (nanoA)
// T_type_Ca_channel_d_gate_alpha_d_T   (per_second)
// T_type_Ca_channel_d_gate_beta_d_T   (per_second)
// T_type_Ca_channel_d_gate_d_T_infinity   (dimensionless)
// T_type_Ca_channel_d_gate_tau_d_T   (second)
// T_type_Ca_channel_f_gate_alpha_f_T   (per_second)
// T_type_Ca_channel_f_gate_beta_f_T   (per_second)
// T_type_Ca_channel_f_gate_f_T_infinity   (dimensionless)
// T_type_Ca_channel_f_gate_tau_f_T   (second)
// T_type_Ca_channel_g_Ca_T   (microS)
// T_type_Ca_channel_i_Ca_T   (nanoA)
    
 double membrane_FCell,calcium_background_current_g_b_Ca,reversal_and_equilibrium_potentials_E_Ca,calcium_background_current_i_b_Ca,four_AP_sensitive_currents_g_to,four_AP_sensitive_currents_g_sus,reversal_and_equilibrium_potentials_E_K,four_AP_sensitive_currents_i_to,
        four_AP_sensitive_currents_i_sus,four_AP_sensitive_currents_q_gate_q_infinity,four_AP_sensitive_currents_q_gate_tau_q,four_AP_sensitive_currents_r_gate_r_infinity,
        four_AP_sensitive_currents_r_gate_tau_r,hyperpolarisation_activated_current_g_f_Na,reversal_and_equilibrium_potentials_E_Na,hyperpolarisation_activated_current_i_f_Na,
        hyperpolarisation_activated_current_g_f_K,hyperpolarisation_activated_current_i_f_K,hyperpolarisation_activated_current_y_gate_alpha_y,
        hyperpolarisation_activated_current_y_gate_beta_y,L_type_Ca_channel_g_Ca_L,L_type_Ca_channel_i_Ca_L,L_type_Ca_channel_d_gate_d_L_infinity,
        L_type_Ca_channel_d_gate_alpha_d_L,L_type_Ca_channel_d_gate_beta_d_L,L_type_Ca_channel_d_gate_tau_d_L,L_type_Ca_channel_f_gate_f_L_infinity,
        L_type_Ca_channel_f_gate_alpha_f_L,L_type_Ca_channel_f_gate_beta_f_L,L_type_Ca_channel_f_gate_tau_f_L,membrane_Cm,sodium_current_g_Na,sodium_current_h_gate_F_Na,
        sodium_current_h_gate_h,sodium_current_i_Na,T_type_Ca_channel_g_Ca_T,T_type_Ca_channel_i_Ca_T,rapid_delayed_rectifying_potassium_current_g_K_r,rapid_delayed_rectifying_potassium_current_P_a,rapid_delayed_rectifying_potassium_current_i_K_r,
        slow_delayed_rectifying_potassium_current_g_K_s,reversal_and_equilibrium_potentials_E_K_s,slow_delayed_rectifying_potassium_current_i_K_s,
        sodium_background_current_g_b_Na,sodium_background_current_i_b_Na,potassium_background_current_g_b_K,potassium_background_current_i_b_K,
        sodium_calcium_exchanger_k_NaCa,sodium_calcium_exchanger_i_NaCa,sodium_potassium_pump_i_p_max,sodium_potassium_pump_i_p,persistent_calcium_current_i_Ca_p_max,
        persistent_calcium_current_i_Ca_p,rapid_delayed_rectifying_potassium_current_P_af_gate_P_af_infinity,rapid_delayed_rectifying_potassium_current_P_af_gate_tau_P_af,
        rapid_delayed_rectifying_potassium_current_P_as_gate_P_as_infinity,rapid_delayed_rectifying_potassium_current_P_as_gate_tau_P_as,
        rapid_delayed_rectifying_potassium_current_P_i_gate_P_i_infinity,slow_delayed_rectifying_potassium_current_xs_gate_alpha_xs,slow_delayed_rectifying_potassium_current_xs_gate_beta_xs,
        sodium_current_h_gate_h1_infinity,sodium_current_h_gate_tau_h1,sodium_current_h_gate_h2_infinity,sodium_current_h_gate_tau_h2,
        sodium_current_m_gate_m_infinity,sodium_current_m_gate_tau_m,T_type_Ca_channel_d_gate_d_T_infinity,T_type_Ca_channel_d_gate_alpha_d_T,T_type_Ca_channel_d_gate_beta_d_T,
        T_type_Ca_channel_d_gate_tau_d_T,T_type_Ca_channel_f_gate_f_T_infinity,T_type_Ca_channel_f_gate_alpha_f_T,T_type_Ca_channel_f_gate_beta_f_T,T_type_Ca_channel_f_gate_tau_f_T;
    
    
   public void init(){
       Y[0]=0.29760539675;
       Y[1]=0.064402950262;
       Y[2]=0.03889291759;
       Y[3]=0.04804900895;      
       Y[4]=0.48779845203;      
       Y[5]=-39.013558536;     
       Y[6]= 0.13034201158;    
       Y[7]=0.46960956028;   
       Y[8]=0.87993375273;  
       Y[9]=0.082293827208;   
       Y[10]=0.015905380261;        
       Y[11]=0.01445216109;        
       Y[12]=0.092361701692;      
       Y[13]=0.42074047435;       
       Y[14]= 0.038968420558;
        calcium_background_current_g_b_Ca_Centre = 1.32e-5;   // microS
        calcium_background_current_g_b_Ca_Periphery = 4.3e-5;   // microS
        four_AP_sensitive_currents_g_sus_Centre = 6.65e-5;   // microS
        four_AP_sensitive_currents_g_sus_Periphery = 0.0114;   // microS
        four_AP_sensitive_currents_g_to_Centre = 0.00491;   // microS
        four_AP_sensitive_currents_g_to_Periphery = 0.03649;   // microS
        hyperpolarisation_activated_current_g_f_K_Centre = 0.000548;   // microS
        hyperpolarisation_activated_current_g_f_K_Periphery = 0.0069;   // microS
        hyperpolarisation_activated_current_g_f_Na_Centre = 0.000548;   // microS
        hyperpolarisation_activated_current_g_f_Na_Periphery = 0.0069;   // microS
        ionic_concentrations_Ca_i = 0.0001;   // millimolar
        ionic_concentrations_Ca_o = 2.0;   // millimolar
        ionic_concentrations_K_i = 140.0;   // millimolar
        ionic_concentrations_K_o = 5.4;   // millimolar
        ionic_concentrations_Na_i = 8.0;   // millimolar
        ionic_concentrations_Na_o = 140.0;   // millimolar
        L_type_Ca_channel_E_Ca_L = 46.4;   // millivolt
        L_type_Ca_channel_g_Ca_L_Centre = 0.0058;   // microS
        L_type_Ca_channel_g_Ca_L_Periphery = 0.0659;   // microS
        membrane_CmCentre = 2.0e-5;   // microF
        membrane_CmPeriphery = 6.5e-5;   // microF
        membrane_dCell = 0.0;   // dimensionless
        membrane_F = 96845.0;   // coulomb_per_mole
        membrane_R = 8314.0;   // millijoule_per_mole_kelvin
        membrane_T = 310.0;   // kelvin
        persistent_calcium_current_i_Ca_p_max_Centre = 0.0;   // nanoA
        persistent_calcium_current_i_Ca_p_max_Periphery = 0.0;   // nanoA
        potassium_background_current_g_b_K_Centre = 2.52e-5;   // microS
        potassium_background_current_g_b_K_Periphery = 8.19e-5;   // microS
        rapid_delayed_rectifying_potassium_current_g_K_r_Centre = 0.000797;   // microS
        rapid_delayed_rectifying_potassium_current_g_K_r_Periphery = 0.016;   // microS
        rapid_delayed_rectifying_potassium_current_P_i_gate_tau_P_i = 0.002;   // second
        slow_delayed_rectifying_potassium_current_g_K_s_Centre = 0.000518;   // microS
        slow_delayed_rectifying_potassium_current_g_K_s_Periphery = 0.0104;   // microS
        sodium_background_current_g_b_Na_Centre = 5.8e-5;   // microS
        sodium_background_current_g_b_Na_Periphery = 0.000189;   // microS
        sodium_calcium_exchanger_d_NaCa = 0.0001;   // dimensionless
        sodium_calcium_exchanger_gamma_NaCa = 0.5;   // dimensionless
        sodium_calcium_exchanger_k_NaCa_Centre = 2.7e-6;   // nanoA
        sodium_calcium_exchanger_k_NaCa_Periphery = 8.8e-6;   // nanoA
        sodium_current_g_Na_Centre = 0.0;   // microS
        sodium_current_g_Na_Periphery = 1.2e-6;   // microS
        sodium_potassium_pump_i_p_max_Centre = 0.0478;   // nanoA
        sodium_potassium_pump_i_p_max_Periphery = 0.16;   // nanoA
        sodium_potassium_pump_K_m_K = 0.621;   // millimolar
        sodium_potassium_pump_K_m_Na = 5.64;   // millimolar
        T_type_Ca_channel_E_Ca_T = 45.0;   // millivolt
        T_type_Ca_channel_g_Ca_T_Centre = 0.0043;   // microS
        T_type_Ca_channel_g_Ca_T_Periphery = 0.0139;   // microS
      }
   
   public double log(double a){
       return java.lang.Math.log(a);
   }
   public double pow(double a, double b){
       return java.lang.Math.pow(a, b);
    }
   public double exp(double a){
       return java.lang.Math.pow(E, a);
   }
   public double setY(int i,double a){
       Y[i]=a;
       return Y[i];
       }
    public double setYData(int i,int count,double a){
       YData[i][count]=a;
       return YData[i][count];
       }

public void compute(){
 membrane_FCell = 1.07*(3.0*membrane_dCell-0.1)/(3.0*(1.0+0.7745*pow(E,-(3.0*membrane_dCell-2.05)/0.295)));
 System.out.println(Y[5]);
 calcium_background_current_g_b_Ca = calcium_background_current_g_b_Ca_Centre+membrane_FCell*(calcium_background_current_g_b_Ca_Periphery-calcium_background_current_g_b_Ca_Centre);
 reversal_and_equilibrium_potentials_E_Ca = membrane_R*membrane_T/(2.0*membrane_F)*log(ionic_concentrations_Ca_o/ionic_concentrations_Ca_i);
 calcium_background_current_i_b_Ca = calcium_background_current_g_b_Ca*(Y[5]-reversal_and_equilibrium_potentials_E_Ca);
 four_AP_sensitive_currents_g_to = four_AP_sensitive_currents_g_to_Centre+membrane_FCell*(four_AP_sensitive_currents_g_to_Periphery-four_AP_sensitive_currents_g_to_Centre);
 four_AP_sensitive_currents_g_sus = four_AP_sensitive_currents_g_sus_Centre+membrane_FCell*(four_AP_sensitive_currents_g_sus_Periphery-four_AP_sensitive_currents_g_sus_Centre);
 reversal_and_equilibrium_potentials_E_K = membrane_R*membrane_T/membrane_F*log(ionic_concentrations_K_o/ionic_concentrations_K_i);
 four_AP_sensitive_currents_i_to = four_AP_sensitive_currents_g_to*Y[0]*Y[1]*(Y[5]-reversal_and_equilibrium_potentials_E_K);
 four_AP_sensitive_currents_i_sus = four_AP_sensitive_currents_g_sus*Y[1]*(Y[5]-reversal_and_equilibrium_potentials_E_K);
 four_AP_sensitive_currents_q_gate_q_infinity = 1.0/(1.0+exp((Y[5]+59.37)/13.1));
 four_AP_sensitive_currents_q_gate_tau_q = 0.0101+0.06517/(0.57*exp(-0.08*(Y[5]+49.0)))+2.4e-5*exp(0.1*(Y[5]+50.93));
dY[0] =  (four_AP_sensitive_currents_q_gate_q_infinity-Y[0])/four_AP_sensitive_currents_q_gate_tau_q;
 four_AP_sensitive_currents_r_gate_r_infinity = 1.0/(1.0+exp(-(Y[5]-10.93)/19.7));
 four_AP_sensitive_currents_r_gate_tau_r = 0.001*(2.98+15.59/(1.037*exp(0.09*(Y[5]+30.61))+0.369*exp(-0.12*(Y[5]+23.84))));
dY[1] = (four_AP_sensitive_currents_r_gate_r_infinity-Y[1])/four_AP_sensitive_currents_r_gate_tau_r;
 hyperpolarisation_activated_current_g_f_Na = hyperpolarisation_activated_current_g_f_Na_Centre+membrane_FCell*(hyperpolarisation_activated_current_g_f_Na_Periphery-hyperpolarisation_activated_current_g_f_Na_Centre);
 reversal_and_equilibrium_potentials_E_Na = membrane_R*membrane_T/membrane_F*log(ionic_concentrations_Na_o/ionic_concentrations_Na_i);
 hyperpolarisation_activated_current_i_f_Na = hyperpolarisation_activated_current_g_f_Na*Y[2]*(Y[5]-reversal_and_equilibrium_potentials_E_Na);
 hyperpolarisation_activated_current_g_f_K = hyperpolarisation_activated_current_g_f_K_Centre+membrane_FCell*(hyperpolarisation_activated_current_g_f_K_Periphery-hyperpolarisation_activated_current_g_f_K_Centre);
 hyperpolarisation_activated_current_i_f_K = hyperpolarisation_activated_current_g_f_K*Y[2]*(Y[5]-reversal_and_equilibrium_potentials_E_K);
 hyperpolarisation_activated_current_y_gate_alpha_y = exp(-(Y[5]+78.91)/26.62);
 hyperpolarisation_activated_current_y_gate_beta_y = exp((Y[5]+75.13)/21.25);
dY[2] = hyperpolarisation_activated_current_y_gate_alpha_y*(1.0-Y[2])-hyperpolarisation_activated_current_y_gate_beta_y*Y[2];
 L_type_Ca_channel_g_Ca_L  = L_type_Ca_channel_g_Ca_L_Centre+membrane_FCell*(L_type_Ca_channel_g_Ca_L_Periphery-L_type_Ca_channel_g_Ca_L_Centre);
 L_type_Ca_channel_i_Ca_L = L_type_Ca_channel_g_Ca_L*(Y[4]*Y[3]+0.006/(1.0+exp(-(Y[5]+14.1)/6.0)))*(Y[5]-L_type_Ca_channel_E_Ca_L);
 L_type_Ca_channel_d_gate_d_L_infinity = 1.0/(1.0+exp(-(Y[5]+23.1)/6.0));
 L_type_Ca_channel_d_gate_alpha_d_L = -28.38*(Y[5]+35.0)/(exp(-(Y[5]+35.0)/2.5)-1.0)-84.9*Y[5]/(exp(-0.208*Y[5])-1.0);
 L_type_Ca_channel_d_gate_beta_d_L = 11.42*(Y[5]-5.0)/(exp(0.4*(Y[5]-5.0))-1.0);
 L_type_Ca_channel_d_gate_tau_d_L = 2.0/(L_type_Ca_channel_d_gate_alpha_d_L+L_type_Ca_channel_d_gate_beta_d_L);
dY[3] = (L_type_Ca_channel_d_gate_d_L_infinity-Y[3])/L_type_Ca_channel_d_gate_tau_d_L;
 L_type_Ca_channel_f_gate_f_L_infinity = 1.0/(1.0+exp((Y[5]+45.0)/5.0));
 L_type_Ca_channel_f_gate_alpha_f_L = 3.12*(Y[5]+28.0)/(exp((Y[5]+28.0)/4.0)-1.0);
 L_type_Ca_channel_f_gate_beta_f_L = 25.0/(1.0+exp(-(Y[5]+28.0)/4.0));
 L_type_Ca_channel_f_gate_tau_f_L = 1.0/(L_type_Ca_channel_f_gate_alpha_f_L+L_type_Ca_channel_f_gate_beta_f_L);
dY[4] =  (L_type_Ca_channel_f_gate_f_L_infinity-Y[4])/L_type_Ca_channel_f_gate_tau_f_L;
 membrane_Cm = membrane_CmCentre+membrane_FCell*(membrane_CmPeriphery-membrane_CmCentre);
 sodium_current_g_Na = sodium_current_g_Na_Centre+membrane_FCell*(sodium_current_g_Na_Periphery-sodium_current_g_Na_Centre);
 sodium_current_h_gate_F_Na = 0.0952*exp(-0.063*(Y[5]+34.4))/(1.0+1.66*exp(-0.225*(Y[5]+63.7)))+0.0869;
 sodium_current_h_gate_h = (1.0-sodium_current_h_gate_F_Na)*Y[10]+sodium_current_h_gate_F_Na*Y[11];
 sodium_current_i_Na = sodium_current_g_Na*pow(Y[12],3.0)*sodium_current_h_gate_h*ionic_concentrations_Na_o*pow(membrane_F,2.0)/(membrane_R*membrane_T)*(exp((Y[5]-reversal_and_equilibrium_potentials_E_Na)*membrane_F/(membrane_R*membrane_T))-1.0)/(exp(Y[5]*membrane_F/(membrane_R*membrane_T))-1.0)*Y[5];
 T_type_Ca_channel_g_Ca_T = T_type_Ca_channel_g_Ca_T_Centre+membrane_FCell*(T_type_Ca_channel_g_Ca_T_Periphery-T_type_Ca_channel_g_Ca_T_Centre);
 T_type_Ca_channel_i_Ca_T = T_type_Ca_channel_g_Ca_T*Y[13]*Y[14]*(Y[5]-T_type_Ca_channel_E_Ca_T);
 rapid_delayed_rectifying_potassium_current_g_K_r = rapid_delayed_rectifying_potassium_current_g_K_r_Centre+membrane_FCell*(rapid_delayed_rectifying_potassium_current_g_K_r_Periphery-rapid_delayed_rectifying_potassium_current_g_K_r_Centre);
 rapid_delayed_rectifying_potassium_current_P_a = 0.6*Y[6]+0.4*Y[7];
 rapid_delayed_rectifying_potassium_current_i_K_r = rapid_delayed_rectifying_potassium_current_g_K_r*rapid_delayed_rectifying_potassium_current_P_a*Y[8]*(Y[5]-reversal_and_equilibrium_potentials_E_K);
 slow_delayed_rectifying_potassium_current_g_K_s = slow_delayed_rectifying_potassium_current_g_K_s_Centre+membrane_FCell*(slow_delayed_rectifying_potassium_current_g_K_s_Periphery-slow_delayed_rectifying_potassium_current_g_K_s_Centre);
 reversal_and_equilibrium_potentials_E_K_s = membrane_R*membrane_T/membrane_F*log((ionic_concentrations_K_o+0.12*ionic_concentrations_Na_o)/(ionic_concentrations_K_i+0.12*ionic_concentrations_Na_i));
 slow_delayed_rectifying_potassium_current_i_K_s = slow_delayed_rectifying_potassium_current_g_K_s*pow(Y[9],2.0)*(Y[5]-reversal_and_equilibrium_potentials_E_K_s);
 sodium_background_current_g_b_Na = sodium_background_current_g_b_Na_Centre+membrane_FCell*(sodium_background_current_g_b_Na_Periphery-sodium_background_current_g_b_Na_Centre);
 sodium_background_current_i_b_Na = sodium_background_current_g_b_Na*(Y[5]-reversal_and_equilibrium_potentials_E_Na);
 potassium_background_current_g_b_K = potassium_background_current_g_b_K_Centre+membrane_FCell*(potassium_background_current_g_b_K_Periphery-potassium_background_current_g_b_K_Centre);
 potassium_background_current_i_b_K = potassium_background_current_g_b_K*(Y[5]-reversal_and_equilibrium_potentials_E_K);
 sodium_calcium_exchanger_k_NaCa = sodium_calcium_exchanger_k_NaCa_Centre+membrane_FCell*(sodium_calcium_exchanger_k_NaCa_Periphery-sodium_calcium_exchanger_k_NaCa_Centre);
 sodium_calcium_exchanger_i_NaCa = sodium_calcium_exchanger_k_NaCa*(pow(ionic_concentrations_Na_i,3.0)*ionic_concentrations_Ca_o*exp(0.03743*Y[5]*sodium_calcium_exchanger_gamma_NaCa)-pow(ionic_concentrations_Na_o,3.0)*ionic_concentrations_Ca_i*exp(0.0374*Y[5]*(sodium_calcium_exchanger_gamma_NaCa-1.0)))/(1.0+sodium_calcium_exchanger_d_NaCa*(ionic_concentrations_Ca_i*pow(ionic_concentrations_Na_o,3.0)+ionic_concentrations_Ca_o*pow(ionic_concentrations_Na_i,3.0)));
 sodium_potassium_pump_i_p_max = sodium_potassium_pump_i_p_max_Centre+membrane_FCell*(sodium_potassium_pump_i_p_max_Periphery-sodium_potassium_pump_i_p_max_Centre);
 sodium_potassium_pump_i_p = sodium_potassium_pump_i_p_max*pow((ionic_concentrations_Na_i/(sodium_potassium_pump_K_m_Na+ionic_concentrations_Na_i)),3.0)*pow((ionic_concentrations_K_o/(sodium_potassium_pump_K_m_K+ionic_concentrations_K_o)),2.0)*1.6/(1.5+exp(-(Y[5]+60.0)/40.0));
 persistent_calcium_current_i_Ca_p_max = persistent_calcium_current_i_Ca_p_max_Centre+membrane_FCell*(persistent_calcium_current_i_Ca_p_max_Periphery-persistent_calcium_current_i_Ca_p_max_Centre);
 persistent_calcium_current_i_Ca_p = persistent_calcium_current_i_Ca_p_max*ionic_concentrations_Ca_i/(ionic_concentrations_Ca_i+0.0004);
dY[5] =  -1.0/membrane_Cm*(sodium_current_i_Na+L_type_Ca_channel_i_Ca_L+T_type_Ca_channel_i_Ca_T+four_AP_sensitive_currents_i_to+four_AP_sensitive_currents_i_sus+rapid_delayed_rectifying_potassium_current_i_K_r+slow_delayed_rectifying_potassium_current_i_K_s+hyperpolarisation_activated_current_i_f_Na+hyperpolarisation_activated_current_i_f_K+sodium_background_current_i_b_Na+calcium_background_current_i_b_Ca+potassium_background_current_i_b_K+sodium_calcium_exchanger_i_NaCa+sodium_potassium_pump_i_p+persistent_calcium_current_i_Ca_p);
 rapid_delayed_rectifying_potassium_current_P_af_gate_P_af_infinity = 1.0/(1.0+exp(-(Y[5]+14.2)/10.6));
 rapid_delayed_rectifying_potassium_current_P_af_gate_tau_P_af = 1.0/(37.2*exp((Y[5]-9.0)/15.9)+0.96*exp(-(Y[5]-9.0)/22.5));
dY[7] = (rapid_delayed_rectifying_potassium_current_P_as_gate_P_as_infinity-Y[7])/rapid_delayed_rectifying_potassium_current_P_as_gate_tau_P_as;
 rapid_delayed_rectifying_potassium_current_P_i_gate_P_i_infinity = 1.0/(1.0+exp((Y[5]+18.6)/10.1));
dY[8] = (rapid_delayed_rectifying_potassium_current_P_i_gate_P_i_infinity-Y[8])/rapid_delayed_rectifying_potassium_current_P_i_gate_tau_P_i;
 slow_delayed_rectifying_potassium_current_xs_gate_alpha_xs = 14.0/(1.0+pow(E,-(Y[5]-40.0)/9.0));
 slow_delayed_rectifying_potassium_current_xs_gate_beta_xs = pow(E,-Y[5]/45.0);
dY[9] =  (rapid_delayed_rectifying_potassium_current_P_i_gate_P_i_infinity-Y[8])/rapid_delayed_rectifying_potassium_current_P_i_gate_tau_P_i;
 sodium_current_h_gate_h1_infinity = 1.0/(1.0+exp((Y[5]+66.1)/6.4));
 sodium_current_h_gate_tau_h1 = 3.717e-6*exp(-0.2815*(Y[5]+17.11))/(1.0+0.003732*exp(-0.3426*(Y[5]+37.76)))+0.0005977;
dY[10] =  (sodium_current_h_gate_h1_infinity-Y[10])/sodium_current_h_gate_tau_h1;
 sodium_current_h_gate_h2_infinity = sodium_current_h_gate_h1_infinity;
 sodium_current_h_gate_tau_h2 = 3.186e-8*exp(-0.6219*(Y[5]+18.8))/(1.0+7.189e-5*exp(-0.6683*(Y[5]+34.07)))+0.003556;
dY[11] =  (sodium_current_h_gate_h2_infinity-Y[11])/sodium_current_h_gate_tau_h2;
 sodium_current_m_gate_m_infinity = pow((1.0/(1.0+exp(-Y[5]/5.46))),(1.0/3.0));
 sodium_current_m_gate_tau_m = 0.0006247/(0.832*exp(-0.335*(Y[5]+56.7))+0.627*exp(0.082*(Y[5]+65.01)))+4.0e-5;
dY[12] = (sodium_current_m_gate_m_infinity-Y[12])/sodium_current_m_gate_tau_m;
 T_type_Ca_channel_d_gate_d_T_infinity = 1.0/(1.0+exp(-(Y[5]+37.0)/6.8));
 T_type_Ca_channel_d_gate_alpha_d_T = 1068.0*exp((Y[5]+26.3)/30.0);
 T_type_Ca_channel_d_gate_beta_d_T = 1068.0*exp(-(Y[5]+26.3)/30.0);
 T_type_Ca_channel_d_gate_tau_d_T = 1.0/(T_type_Ca_channel_d_gate_alpha_d_T+T_type_Ca_channel_d_gate_beta_d_T);
dY[13] = (T_type_Ca_channel_d_gate_d_T_infinity-Y[13])/T_type_Ca_channel_d_gate_tau_d_T;
 T_type_Ca_channel_f_gate_f_T_infinity = 1.0/(1.0+exp((Y[5]+71.0)/9.0));
 T_type_Ca_channel_f_gate_alpha_f_T = 15.3*exp(-(Y[5]+71.7)/83.3);
 T_type_Ca_channel_f_gate_beta_f_T = 15.0*exp((Y[5]+71.7)/15.38);
 T_type_Ca_channel_f_gate_tau_f_T = 1.0/(T_type_Ca_channel_f_gate_alpha_f_T+T_type_Ca_channel_f_gate_beta_f_T);
dY[14] = (T_type_Ca_channel_f_gate_f_T_infinity-Y[14])/T_type_Ca_channel_f_gate_tau_f_T;

}

 public void txtout(){
     
    try {
      PrintStream out = new PrintStream(new FileOutputStream("Data.txt"));
      
      
      for(int j=0;j<14;j++){
          out.println("membrane_V"+"\n");
      for(int i=0;i<301;i++){
         out.println(YData[j][i]+" at time  "+tData[i]);
      }
      }
      out.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }
 }
 
 public void visualise(){
     //visualise data
  XYSeries series = new XYSeries("Membrane mVolt");
  XYSeries series2 = new XYSeries("potassium_channel_n_gate_n");
  XYSeries series3 = new XYSeries("sodium_channel_h_gate_h");
  XYSeries series4 = new XYSeries("sodium_channel_m_gate_m");
  //add data to charts
  for(int l=0;l<301;l++){
  series.add(tData[l], YData[0][l]);
  series2.add(tData[l], YData[1][l]);
  series3.add(tData[l], YData[2][l]);
  series4.add(tData[l], YData[3][l]);
  }
  //create lines
  XYDataset xyDataset = new XYSeriesCollection(series);
  final XYSeriesCollection dataset = new XYSeriesCollection();
  dataset.addSeries(series2);
  dataset.addSeries(series3);
  dataset.addSeries(series4);
  //visualize
  JFreeChart chart = ChartFactory.createXYLineChart("Membrane", "time", "data",xyDataset, PlotOrientation.VERTICAL, true, true, false);
  JFreeChart chart5 = ChartFactory.createXYLineChart("Channels", "time", "data",dataset, PlotOrientation.VERTICAL, true, true, false);
  //open window
  ChartFrame frame1=new ChartFrame(" ",chart);
  ChartFrame frame5=new ChartFrame(" ",chart5);
  //make window
  frame1.setVisible(true);
  frame1.setSize(900,900);
  frame5.setVisible(true);
  frame5.setSize(900,900);
 }
 

public static void main(String[] args) {
      HCell a= new HCell();
       a.init();
       int count=0;
       
    while(a.t<=a.time+a.h/2.0){
        if(a.t>=a.tt-a.h/2.0){
          a.tData[count]=a.t;
          a.YData[0][count]=a.setYData(0, count, a.Y[0]);
          a.YData[1][count]=a.setYData(1, count, a.Y[1]);
          a.YData[2][count]=a.setYData(2, count, a.Y[2]);
          a.YData[3][count]=a.setYData(3, count, a.Y[3]);
          a.YData[4][count]=a.setYData(4, count, a.Y[4]);
          a.YData[5][count]=a.setYData(5, count, a.Y[5]);
          a.YData[6][count]=a.setYData(6, count, a.Y[6]);
          a.YData[7][count]=a.setYData(7, count, a.Y[7]);
          a.YData[8][count]=a.setYData(8, count, a.Y[8]);
          a.YData[9][count]=a.setYData(9, count, a.Y[9]);
          a.YData[10][count]=a.setYData(10, count, a.Y[10]);
          a.YData[11][count]=a.setYData(11, count, a.Y[11]);
          a.YData[12][count]=a.setYData(12, count, a.Y[12]);
          a.YData[13][count]=a.setYData(13, count, a.Y[13]);
          a.YData[14][count]=a.setYData(14, count, a.Y[14]);
          a.tt=a.tt+a.hh;
          count++;
         //System.out.println(a.YData[4][count]);
        }
        a.compute();
        
        for(int i=0;i<14;i++){
            a.setY(i,a.Y[i]+(a.h*a.dY[i]));
            
        }
       a.t=a.t+a.h;
      }
 
//a.visualise();
a.txtout();
  
  }
    
}

