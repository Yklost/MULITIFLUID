#include <iostream> 
#include <cmath>
#include "constants.h"
#include "mpi.h"
#include "physics.h"


void get_pflux(int direction, double* u, double* fl)
{

    double v_x_i, v_y_i, v_z_i;
    double v_x_n, v_y_n, v_z_n;
    double pre_i, pre_n, pre_e;

    v_x_i = u[m_x_i_] / u[rho_i_]; 
    v_y_i = u[m_y_i_] / u[rho_i_];
    v_z_i = u[m_z_i_] / u[rho_i_];
    
    v_x_n = u[m_x_n_] / u[rho_n_]; 
    v_y_n = u[m_y_n_] / u[rho_n_];
    v_z_n = u[m_z_n_] / u[rho_n_];    
    
    pre_n =  ( u[ene_n_] - 0.5 * ( u[m_x_n_]*v_x_n + u[m_y_n_]*v_y_n + u[m_z_n_]*v_z_n )) * (gam - 1.) ;
 
    pre_i =  ( u[ene_i_] - 0.5 * ( u[m_x_i_]*v_x_i + u[m_y_i_]*v_y_i + u[m_z_i_]*v_z_i )) * (gam - 1.) ;
  
    pre_e =  ( u[ene_e_] - 0.5 * me_mi * ( u[m_x_i_]*v_x_i + u[m_y_i_]*v_y_i + u[m_z_i_]*v_z_i ) - u[rho_i_] / mi * III_erg ) * (gam - 1.);
 
 
    if (direction == 0) {      
    fl[rho_i_] = u[m_x_i_]; 

    fl[m_x_i_] = u[m_x_i_] * v_x_i + pre_i; 
    fl[m_y_i_] = u[m_y_i_] * v_x_i; 
    fl[m_z_i_] = u[m_z_i_] * v_x_i;
    
    fl[ene_i_] = v_x_i * (u[ene_i_] + pre_i); 


    fl[rho_n_] = u[m_x_n_]; 

    fl[m_x_n_] = u[m_x_n_] * v_x_n + pre_n; 
    fl[m_y_n_] = u[m_y_n_] * v_x_n ; 
    fl[m_z_n_] = u[m_z_n_] * v_x_n ;
    
    fl[ene_n_] = v_x_n * (u[ene_n_] + pre_n) ;   
    
    fl[ene_e_] = v_x_i * (u[ene_e_] + pre_e) + u[q_x_];

    }
    


    if (direction == 1) {
    fl[rho_i_]=u[m_y_i_];

    fl[m_x_i_] = u[m_x_i_] * v_y_i;
    fl[m_y_i_] = u[m_y_i_] * v_y_i + pre_i;
    fl[m_z_i_] = u[m_z_i_] * v_y_i;
    
    fl[ene_i_] = v_y_i * (u[ene_i_] + pre_i);

    fl[rho_n_]=u[m_y_n_];

    fl[m_x_n_] = u[m_x_n_] * v_y_n ;
    fl[m_y_n_] = u[m_y_n_] * v_y_n + pre_n ;
    fl[m_z_n_] = u[m_z_n_] * v_y_n ;
    
    fl[ene_n_] = v_y_n * (u[ene_n_] + pre_n);
    
    fl[ene_e_] = v_y_i * (u[ene_e_] + pre_e) + u[q_y_] ;

    }
    
    
    
    if (direction == 2) {
    fl[rho_i_]=u[m_z_i_];

    fl[m_x_i_] = u[m_x_i_] * v_z_i ;
    fl[m_y_i_] = u[m_y_i_] * v_z_i ;
    fl[m_z_i_] = u[m_z_i_] * v_z_i + pre_i;
    
    fl[ene_i_] = v_z_i * (u[ene_i_] + pre_i);
    
    fl[rho_n_]=u[m_z_n_];

    fl[m_x_n_] = u[m_x_n_] * v_z_n ;
    fl[m_y_n_] = u[m_y_n_] * v_z_n ;
    fl[m_z_n_] = u[m_z_n_] * v_z_n + pre_n ;
    
    fl[ene_n_] = v_z_n * (u[ene_n_] + pre_n);
    
    fl[ene_e_] = v_z_i * (u[ene_e_] + pre_e) +u[q_z_];   

    }   

}


double astar(double* u1, double* u2, int dir)
{
double fs1, fs2, aaa1, aaa2, aaa;
 
    fs1=fast_cfl(u1); 
    fs2=fast_cfl(u2); 

//    aaa=fmax(fabs(u1[m_x_+dir]/u1[rho_]+fs1),fabs(u2[dir]/u2[rho_]+fs2)); 
//    aaa=fmax(aaa, fabs(u1[m_x_+dir]/u1[rho_]-fs1));
//    aaa=fmax(aaa, fabs(u2[m_x_+dir]/u2[rho_]-fs2));
    aaa1=fmax(fs1+fabs(u1[m_x_i_+dir]/u1[rho_i_]),fs2+fabs(u2[m_x_i_+dir]/u2[rho_i_])); //ad hoc?
    aaa2=fmax(fs1+fabs(u1[m_x_n_+dir]/u1[rho_n_]),fs2+fabs(u2[m_x_n_+dir]/u2[rho_n_])); //ad hoc?
    
    aaa = fmax (aaa1, aaa2);

    return aaa;

}



double cs_cfl2(double* u)
{
double pl_1, pl_2, cs2_1, cs2_2, cs2;
    pl_1=(gam-1.0)*(u[ene_i_]-(u[m_x_i_]*u[m_x_i_]+u[m_y_i_]*u[m_y_i_]+u[m_z_i_]*u[m_z_i_])/(2.0*u[rho_i_]));
    cs2_1=gam*pl_1/u[rho_i_]; //speed of sound fluid 1 magnetic
    
    pl_2=(gam-1.0)*(u[ene_n_]-(u[m_x_n_]*u[m_x_n_]+u[m_y_n_]*u[m_y_n_]+u[m_z_n_]*u[m_z_n_])/(2.0*u[rho_n_]));
    cs2_2=gam*pl_2/u[rho_n_]; //speed of sound fluid 2 non-magnetic    
    
    cs2 = fmax(cs2_1, cs2_2);
    
    return cs2;
}



double fast_cfl(double* u)
{
 double cs2;
 cs2 = cs_cfl2(u);
 return sqrt(cs2);
}





double timestep(double*** u)
{
    int i, j;
    int i_var;
    double dt_local, dt, c_max;    
    double v_fast, v_flow, v_flow_1, v_flow_2;
    double u_loc[n_var];
    
    c_max=0.0;
    
for (i = 0; i < n_x; i++)
{
  for (j = 0; j < n_y; j++)
  {
	  for (i_var = 0; i_var < n_consvar; i_var++) u_loc[i_var]=u[i_var][i][j];  
	  
	    v_fast = fast_cfl(u_loc);
            v_flow_1 = sqrt(u[v_x_i_][i][j]*u[v_x_i_][i][j] + u[v_y_i_][i][j]*u[v_y_i_][i][j] + u[v_z_i_][i][j]*u[v_z_i_][i][j]);
            v_flow_2 = sqrt(u[v_x_n_][i][j]*u[v_x_n_][i][j] + u[v_y_n_][i][j]*u[v_y_n_][i][j] + u[v_z_n_][i][j]*u[v_z_n_][i][j]);
	    v_flow = fmax(v_flow_1, v_flow_2);
            c_max=fmax(c_max, v_fast + v_flow); //max speed of fast speed
  }
} 

dt_local=CFL*fmin(dx,dy)/c_max;

MPI_Allreduce(&dt_local, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); 

return dt;

}





void prim_to_con(double*** u)
{
int i, j;

for (i = 0; i < n_x; i++)
{
  for (j = 0; j < n_y; j++)
  {
        u[m_x_i_][i][j] = u[v_x_i_][i][j] * u[rho_i_][i][j];
        u[m_y_i_][i][j] = u[v_y_i_][i][j] * u[rho_i_][i][j];
	u[m_z_i_][i][j] = u[v_z_i_][i][j] * u[rho_i_][i][j];
        u[m_x_n_][i][j] = u[v_x_n_][i][j] * u[rho_n_][i][j];
        u[m_y_n_][i][j] = u[v_y_n_][i][j] * u[rho_n_][i][j];
	u[m_z_n_][i][j] = u[v_z_n_][i][j] * u[rho_n_][i][j];
	
        u[ene_i_][i][j] = u[pre_i_][i][j] / (gam -1.) + 0.5* (u[m_x_i_][i][j] * u[v_x_i_][i][j] + 
	                                                      u[m_y_i_][i][j] * u[v_y_i_][i][j] + 
							      u[m_z_i_][i][j] * u[v_z_i_][i][j]);
							      
        u[ene_n_][i][j] = u[pre_n_][i][j] / (gam -1.) + 0.5* (u[m_x_n_][i][j] * u[v_x_n_][i][j] + 
	                                                      u[m_y_n_][i][j] * u[v_y_n_][i][j] + 
							      u[m_z_n_][i][j] * u[v_z_n_][i][j]);
							      
	u[ene_e_][i][j] = u[pre_e_][i][j] / (gam -1.) + me_mi * 0.5*(u[v_x_i_][i][j] * u[m_x_i_][i][j] + 
	                                                             u[v_y_i_][i][j] * u[m_y_i_][i][j] + 
						                     u[v_z_i_][i][j] * u[m_z_i_][i][j]) +
								     (u[rho_i_][i][j] * III_erg / mi);
								     
	u[tem_i_][i][j] = u[pre_i_][i][j] / u[rho_i_][i][j] * mi / k_b;
	u[tem_n_][i][j] = u[pre_n_][i][j] / u[rho_n_][i][j] * (mi+me) / k_b;
	u[tem_e_][i][j] = u[pre_e_][i][j] / u[rho_i_][i][j] * mi / k_b;


  }
} 
       
}




void con_to_prim(double*** u)
{
int i, j;

for (i = 0; i < n_x; i++)
{
  for (j = 0; j < n_y; j++)
  {

        u[v_x_i_][i][j] = u[m_x_i_][i][j] / u[rho_i_][i][j];
        u[v_y_i_][i][j] = u[m_y_i_][i][j] / u[rho_i_][i][j];
	u[v_z_i_][i][j] = u[m_z_i_][i][j] / u[rho_i_][i][j];
        u[v_x_n_][i][j] = u[m_x_n_][i][j] / u[rho_n_][i][j];
        u[v_y_n_][i][j] = u[m_y_n_][i][j] / u[rho_n_][i][j];
	u[v_z_n_][i][j] = u[m_z_n_][i][j] / u[rho_n_][i][j];	

        u[pre_i_][i][j] = (u[ene_i_][i][j] - 0.5*(u[v_x_i_][i][j] * u[m_x_i_][i][j] + 
	                                          u[v_y_i_][i][j] * u[m_y_i_][i][j] + 
						  u[v_z_i_][i][j] * u[m_z_i_][i][j]))  * (gam -1.);

        u[pre_n_][i][j] = (u[ene_n_][i][j] - 0.5*(u[v_x_n_][i][j] * u[m_x_n_][i][j] + 
	                                          u[v_y_n_][i][j] * u[m_y_n_][i][j] + 
						  u[v_z_n_][i][j] * u[m_z_n_][i][j]))  * (gam -1.);

        u[pre_e_][i][j] = (u[ene_e_][i][j] - me_mi * 0.5*(u[v_x_i_][i][j] * u[m_x_i_][i][j] + 
	                                                  u[v_y_i_][i][j] * u[m_y_i_][i][j] + 
						          u[v_z_i_][i][j] * u[m_z_i_][i][j]) - III_erg * u[rho_i_][i][j] / mi )  * (gam -1.);
	
//	u[tem_i_][i][j] = u[pre_i_][i][j] / u[rho_i_][i][j] * mi / k_b;
	
//	u[tem_n_][i][j] = u[pre_n_][i][j] / u[rho_n_][i][j] * (mi+me) / k_b;
	
//	u[tem_e_][i][j] = u[pre_e_][i][j] / u[rho_i_][i][j] * mi / k_b;						  


        u[tem_i_][i][j] = u[pre_i_][i][j] / u[rho_i_][i][j] * mi / k_b;
	u[tem_n_][i][j] = u[pre_n_][i][j] / u[rho_n_][i][j] * (mi+me) / k_b;
	u[tem_e_][i][j] = u[pre_e_][i][j] / u[rho_i_][i][j] * mi / k_b;	




  }
} 

}




void calc_grav_sources(double** x, double** y, double*** u, double*** res, double time)
{ 

    int i,j,k;
    double radius;
    double dmom;

    for (i = 2; i < n_x - 2; i++)
    { 
        for (j = 2; j < n_y - 2; j++)
        {
               res[m_y_i_][i][j]=res[m_y_i_][i][j]-u[rho_i_][i][j]*ggg;
               res[ene_i_][i][j]=res[ene_i_][i][j]-u[m_y_i_][i][j]*ggg;

               res[m_y_n_][i][j]=res[m_y_n_][i][j]-u[rho_n_][i][j]*ggg;
               res[ene_n_][i][j]=res[ene_n_][i][j]-u[m_y_n_][i][j]*ggg;

        } //j
    } //i

}







/*
void calc_visc_sources(double** x, double** y,  double*** u, double*** res, double time)
{ 

    int i,j,k;
    double radius;
    double dmom;

    for (i = 2; i < n_x - 2; i++)
    { 
        for (j = 2; j < n_y - 2; j++)
        {

	       //	   res[m_x_][i][j] = res[m_x_][i][j] - u[rho_][i][j] * nu_visc * (//( u[v_x_][i+1][j] - 2.0*u[v_x_][i][j] + u[v_x_][i-1][j] ) / ( dx * dx ) 
	       //	                                                                + ( u[v_x_][i][j+1] - 2.0*u[v_x_][i][j] + u[v_x_][i][j-1] ) / ( dy * dy ));

	       //	   res[m_y_][i][j] = res[m_y_][i][j] - u[rho_][i][j] * nu_visc * (( u[v_y_][i+1][j] - 2.0*u[v_y_][i][j] + u[v_y_][i-1][j] ) / ( dx * dx ) );
	                                                                	       //+ ( u[v_y_][i][j+1] - 2.0*u[v_y_][i][j] + u[v_y_][i][j-1] ) / ( dy * dy ));	   

	       // Make is a separate routine     
	        	   radius=sqrt(x[i][j]*x[i][j]+y[i][j]*y[i][j]);	   
			  if (radius < 0.1) {
			      dmom=0.01*sin(10.0*time);
			      u[m_x_][i][j]=u[m_x_][i][j]+dmom;
			      u[ene_][i][j]=u[ene_][i][j]+dmom*dmom/u[rho_][i][j];
        		  }

        } //j
    } //i

}
*/

