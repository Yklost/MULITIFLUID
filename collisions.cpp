#include <iostream> 
#include "collisions.h"
#include "constants.h"
#include "math.h"

void calc_coll_sources(double** x, double** y, double*** u, double*** res, double time)
{ 

    int i,j,k;
    
    double v_x_i, v_x_n, v_y_i, v_y_n, v_z_i, v_z_n;
    
    double v_x_r, v_y_r, v_z_r;
    
    double v_n2, v_i2;
    

    double rho_i_change, m_x_i_change, m_y_i_change, m_z_i_change, ene_i_change;
    double rho_n_change, m_x_n_change, m_y_n_change, m_z_n_change, ene_n_change;
    double ene_e_change;
    
    double rho_rec, rho_ion;
    
    double i_change_coef, n_change_coef;
    
    double v_i_dot_v_r, v_n_dot_v_r;
    
    
//    get_collision_rates(x, y, u);
    
    
    for (i = 2; i < n_x - 2; i++)
    { 
        for (j = 2; j < n_y - 2; j++)
        {
	
	    v_x_i = u[v_x_i_][i][j] ; //u[m_x_i_][i][j]/u[rho_i_][i][j];
	    v_y_i = u[v_y_i_][i][j] ; //u[m_y_i_][i][j]/u[rho_i_][i][j];
	    v_z_i = u[v_z_i_][i][j] ; //u[m_z_i_][i][j]/u[rho_i_][i][j];	    
	    
	    v_x_n = u[v_x_n_][i][j] ; //u[m_x_n_][i][j]/u[rho_n_][i][j];
	    v_y_n = u[v_y_n_][i][j] ; //u[m_y_n_][i][j]/u[rho_n_][i][j];
	    v_z_n = u[v_z_n_][i][j] ; //u[m_z_n_][i][j]/u[rho_n_][i][j];
	    
	    
	    v_x_r = v_x_i - v_x_n;
	    v_y_r = v_y_i - v_y_n;
	    v_z_r = v_z_i - v_z_n;

            v_n2 = v_x_n * v_x_n + v_y_n * v_y_n + v_z_n * v_z_n;
	    v_i2 = v_x_i * v_x_i + v_y_i * v_y_i + v_z_i * v_z_i;

	    
	    //Does not conserve full mass.	    
	    
	    rho_i_change = - u[k_rec_][i][j] * u[rho_n_][i][j] * u[rho_i_][i][j] / (me+mi) + u[k_ion_][i][j] * u[rho_i_][i][j] * u[rho_n_][i][j] / (me+mi);  //recombination - ionization

	    rho_n_change =   u[k_rec_][i][j] * u[rho_n_][i][j] * u[rho_i_][i][j] / mi      - u[k_ion_][i][j] * u[rho_i_][i][j] * u[rho_n_][i][j] / mi;  //recombination - ionization
	    
	    
	    
	    i_change_coef = - u[a_ni_][i][j] - u[rho_i_][i][j] * u[rho_n_][i][j] / (me+mi) * u[k_ion_][i][j];
	    
	    n_change_coef =   u[a_ne_][i][j] + u[a_ni_][i][j] + u[rho_n_][i][j] * u[rho_i_][i][j] / mi * u[k_ion_][i][j];
	    
	    
	    m_x_i_change =  i_change_coef * v_x_r;
	    
	    m_y_i_change =  i_change_coef * v_y_r;
	    
	    m_z_i_change =  i_change_coef * v_z_r;
	    
	    
	    
	    m_x_n_change =  n_change_coef * v_x_r; //note minus
	    
	    m_y_n_change =  n_change_coef * v_y_r;
	    
	    m_z_n_change =  n_change_coef * v_z_r;	    
	    
				
				
	    v_i_dot_v_r = u[v_x_i_][i][j] * v_x_r + u[v_y_i_][i][j] * v_y_r + u[v_z_i_][i][j] * v_z_r;			
				
	    v_n_dot_v_r = u[v_x_n_][i][j] * v_x_r + u[v_y_n_][i][j] * v_y_r + u[v_z_n_][i][j] * v_z_r;				
				
				
					       
						       
	    ene_i_change =  u[QQQ_][i][j] - u[a_ni_][i][j] * v_i_dot_v_r + 0.5 * u[rho_i_][i][j] * v_n2 * u[k_ion_][i][j] * u[rho_n_][i][j] / (me+mi) -
	                                                                   0.5 * u[rho_i_][i][j] * v_i2 * u[k_rec_][i][j] * u[rho_n_][i][j] / (me+mi);				       
						       
						       
            ene_n_change =  u[QQQ_][i][j] + (u[a_ni_][i][j] + u[a_ne_][i][j]) * v_n_dot_v_r - 0.5 * u[rho_n_][i][j] * v_n2 * u[k_ion_][i][j] * u[rho_i_][i][j] / mi +
	                                                                                      0.5 * u[rho_n_][i][j] * v_i2 * u[k_rec_][i][j] * u[rho_i_][i][j] / mi;
	    
	    
	    ene_e_change =  u[QQQ_][i][j] - u[a_ne_][i][j] * v_i_dot_v_r;// - beta_0 * u[rho_i_][i][j] / mi * u[vidte_][i][j];// - 
	                                 //  (me_mi * 0.5 * u[rho_i_][i][j] * v_n2 + u[rho_i_][i][j] / mi * III_erg) * u[k_ion_][i][j] * u[rho_n_][i][j] / (me+mi);







	    	

            res[rho_i_][i][j] = res[rho_i_][i][j] - rho_i_change;  // -(rec - ion)

            res[rho_n_][i][j] = res[rho_n_][i][j] - rho_n_change;  // +(rec - ion)

	    
	    res[m_x_i_][i][j] = res[m_x_i_][i][j] -  m_x_i_change ;
	    
	    res[m_y_i_][i][j] = res[m_y_i_][i][j] -  m_y_i_change ;

	    res[m_z_i_][i][j] = res[m_z_i_][i][j] -  m_z_i_change ;


	    res[m_x_n_][i][j] = res[m_x_n_][i][j] -  m_x_n_change ;
	    
	    res[m_y_n_][i][j] = res[m_y_n_][i][j] -  m_y_n_change ;

	    res[m_z_n_][i][j] = res[m_z_n_][i][j] -  m_z_n_change ;
	    
	    
	    res[ene_i_][i][j] = res[ene_i_][i][j] -  ene_i_change ;
	    
	    res[ene_n_][i][j] = res[ene_n_][i][j] -  ene_n_change ;

	 
	    res[ene_e_][i][j] = res[ene_e_][i][j] -  ene_e_change ;



        } //j
    } //i

}






void get_collision_rates(double** x, double** y, double*** u)
{
int i,j;
double dte_x, dte_y, dte_z;
double dpe_x, dpe_y, dpe_z;
double dpn_x, dpn_y, dpn_z;
double dpi_x, dpi_y, dpi_z;

double dpepi_x, dpepi_y;

double dxx = x[1][0]-x[0][0];
double dyy = y[0][1]-y[0][0];

double rho_t;

double rho_n_rho_t;
double rho_i_rho_t;


double lambda, tau_e, sigma_tr_ion, sigma_coul, freq_ei, a_ei;

double v_x_r, v_y_r, v_z_r;

double velAverTher_e; //average thermal velocity of electron
double sigmaIoniz_e; //characterisitic ionization section


    
    for (i = 2; i < n_x - 2; i++)
    { 
        for (j = 2; j < n_y - 2; j++)
        {

//rho: u[rho_i_][i][j],u[rho_n_][i][j]
//tem: u[tem_i_][i][j],n,e
//mi, me, me_mi

	       
 	       
	       dte_x = 0.5 * (u[tem_e_][i+1][j] - u[tem_e_][i-1][j]) / dxx;
	       dte_y = 0.5 * (u[tem_e_][i][j+1] - u[tem_e_][i][j-1]) / dyy;
	       dte_z = 0.0;	       

	       dpi_x = 0.5 * (u[pre_i_][i+1][j] - u[pre_i_][i-1][j]) / dxx;
	       dpi_y = 0.5 * (u[pre_i_][i][j+1] - u[pre_i_][i][j-1]) / dyy;
	       dpe_z = 0.0;
	       
	       dpn_x = 0.5 * (u[pre_n_][i+1][j] - u[pre_n_][i-1][j]) / dxx;
	       dpn_y = 0.5 * (u[pre_n_][i][j+1] - u[pre_n_][i][j-1]) / dyy;
	       dpn_z = 0.0;
	       
	       dpe_x = 0.5 * (u[pre_e_][i+1][j] - u[pre_e_][i-1][j]) / dxx;
	       dpe_y = 0.5 * (u[pre_e_][i][j+1] - u[pre_e_][i][j-1]) / dyy;
	       dpe_z = 0.0;	
	       
	       
	       dpepi_x = dpe_x + dpi_x;
	       dpepi_y = dpe_y + dpi_y;
	       
	       	       
	       v_x_r = u[v_x_i_][i][j] - u[v_x_n_][i][j];
	       v_y_r = u[v_y_i_][i][j] - u[v_y_n_][i][j];
	       v_z_r = u[v_z_i_][i][j] - u[v_z_n_][i][j];
	       	       
	       
	       rho_t = u[rho_i_][i][j] * me_mi + u[rho_i_][i][j] + u[rho_n_][i][j];  //rho_e + rho_i + rho_n
	       
	       rho_n_rho_t = u[rho_n_][i][j] / rho_t;

	       rho_i_rho_t = u[rho_i_][i][j] / rho_t;	

	       //velAverTher_e = sqrt(8.0*u[tem_e_][i][j]*eVtoerg/PI/me); 

               //sigmaIoniz_e = 2.0e-17*u[tem_e_][i][j]; //here C_i = 2.0e-17 is a coeffitient, units cm^2/eV     
     
	   	       
	       lambda = 17.0 ;//23.4 - 1.15 * log10(u[rho_i_][i][j] / mi) + 3.45 * log10(u[tem_e_][i][j]);
	       
	       tau_e = (3.5e5 / lambda) * pow(u[tem_e_][i][j] / 11600.0, 1.5) / (u[rho_i_][i][j] / mi);  
	       
	       u[k_ion_][i][j] = 2e-17 * sqrt(8.0 * k_b * u[tem_e_][i][j] / PI / me) * (III_erg / k_b / u[tem_e_][i][j] + 2.e0) * exp (-III_erg / k_b / u[tem_e_][i][j]); //tem in K 
	       
	       u[k_rec_][i][j] = 0.0; 
		
	       sigma_tr_ion = 2.0 * sqrt(2) * 0.88e-16 * sqrt(4.5 * e * e * PI / 2.0 / 0.5e-8 / 4.0 / k_b / u[tem_i_][i][j]); //tem in K, sigma tranport
		
	       sigma_coul  = 4.0 * PI * pow(e, 4) * lambda / 9.0 / pow(k_b*u[tem_e_][i][j], 2.0); //tem in K, sigma coulomb
		
               freq_ei = 4.0 * sqrt(3.0) * pow(e,4.0) * u[rho_i_][i][j] * lambda / 9.0 / mi / sqrt(me) / pow(k_b * u[tem_e_][i][j], 1.5) ; //tem in K, coulomb collisions 
	       
	       u[a_ni_][i][j] = mi * u[rho_i_][i][j] * u[rho_n_][i][j] * sigma_tr_ion * v_y_r / mi / mi;  
	       	       
	       u[a_ne_][i][j] = me * u[rho_i_][i][j] / mi * 10e-16 * u[rho_n_][i][j] / mi * 4.2e7 * pow(u[tem_e_][i][j] / 11600.0, 0.5) ; // tem in K 
		       
	       a_ei = me * u[rho_i_][i][j] * freq_ei / mi;     
  
  
               u[vidte_][i][j] = u[v_x_i_][i][j] * dte_x + u[v_y_i_][i][j] * dte_y + u[v_z_i_][i][j] * dte_z;
	       
	       
	       u[QQQ_][i][j] = 1.0 / (u[a_ni_][i][j] + u[a_ne_][i][j]) * (rho_i_rho_t * rho_i_rho_t * (dpn_x   * dpn_x   + dpn_y   * dpn_y  ) -
	                                                                  rho_n_rho_t * rho_n_rho_t * (dpepi_x * dpepi_x + dpepi_y * dpepi_y));
	       
	       

	       
///MOVE q TO FLUXES! and reshuffle.	       
	       
	       u[q_x_][i][j] = (u[rho_i_][i][j]/mi) * u[tem_e_][i][j] / 11600.0 * ( beta_0 * v_x_r - (tau_e / me) * gamma_0 * dte_x );
	       u[q_y_][i][j] = (u[rho_i_][i][j]/mi) * u[tem_e_][i][j] / 11600.0 * ( beta_0 * v_y_r - (tau_e / me) * gamma_0 * dte_y );
	       u[q_z_][i][j] = 0.0;
	       
	       
	       
	                  


        }
    }
}




