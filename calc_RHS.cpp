#include <iostream>
#include <cmath>
#include "constants.h"
#include "memory.h"
#include "physics.h"
#include "numerical_flux.h"
#include "apply_BC.h"
#include "collisions.h"

void calc_residuals(double** x, double** y, double*** u0, double*** res, double**** nflux)
{ 

int i_var;   
int i, j, k;


//double ullx[n_consvar], ulx[n_consvar], uux[n_consvar] ,urx[n_consvar];
//double ully[n_consvar], uly[n_consvar], uuy[n_consvar] ,ury[n_consvar];

double ullx[n_var], ulx[n_var], uux[n_var] ,urx[n_var];
double ully[n_var], uly[n_var], uuy[n_var] ,ury[n_var];

double flux_x[n_consvar], flux_y[n_consvar];


for (i = ghosts; i < n_x-ghosts+1; i++)
    {    
    for (j = ghosts; j < n_y-ghosts+1; j++)
        {
//	     for (i_var = 0; i_var < n_consvar; i_var++)   
	     for (i_var = 0; i_var < n_var; i_var++)   
                 {       
                           ullx[i_var]= u0[i_var][i-2][j]; 
                           ulx[i_var] = u0[i_var][i-1][j];
                           uux[i_var] = u0[i_var][i][j];
                           urx[i_var] = u0[i_var][i+1][j];
  
                           ully[i_var]= u0[i_var][i][j-2]; 
                           uly[i_var] = u0[i_var][i][j-1];
                           uuy[i_var] = u0[i_var][i][j];
                           ury[i_var] = u0[i_var][i][j+1];

                 } 

             numerical_flux(0,ullx,ulx,uux,urx,flux_x);
             numerical_flux(1,ully,uly,uuy,ury,flux_y);

             for (i_var = 0; i_var < n_consvar; i_var++)
                 {                  
                	   nflux[0][i_var][i][j]=flux_x[i_var]; 
                	   nflux[1][i_var][i][j]=flux_y[i_var];
                 } 
		     
        } //j
    }//i

    for (i_var = 0; i_var < n_consvar; i_var++)   
    {
	for (i = ghosts; i < n_x - ghosts; i++)
	{    
	    for (j = ghosts; j < n_y - ghosts; j++)
	    {
		   res[i_var][i][j] = (nflux[0][i_var][i+1][j]-nflux[0][i_var][i][j])/(x[i+1][j]-x[i][j]) + 
        		              (nflux[1][i_var][i][j+1]-nflux[1][i_var][i][j])/(y[i][j+1]-y[i][j]);
	    }
	}
    }
}





void calc_sources(double** x, double** y, double*** u, double*** res, double time)
{

  calc_grav_sources(x, y, u, res, time);
  calc_coll_sources(x, y, u, res, time);
  
  }









void calc_RHS(double** x, double** y, double*** u, double*** u_init, double*** res, double**** flux, double time)
{

//allocate_5D(flux,   n_dim, n_var, n_x, n_y, n_z);

        con_to_prim(u);

        apply_BC(u, u_init, x, y);
		
        get_collision_rates(x, y, u);
	
        calc_residuals(x, y, u, res, flux); 

        calc_sources(x, y, u, res, time);


//deallocate_5D(flux,   n_dim, n_var, n_x, n_y, n_z);

}



void advance_solution(double*** u1, double*** u0, double*** res, double dtime)
{
 int i_var,i,j,k;
  
 for (i_var=0; i_var<n_consvar; i_var++)
  {
    for (i = ghosts; i < n_x-ghosts; i++)
      { 
        for (j = ghosts; j < n_y-ghosts; j++)
        {
            u1[i_var][i][j]= u0[i_var][i][j] - res[i_var][i][j]*dtime;
	}        
      }
   }
}



void add_arrays(double*** u1, double*** u2, double a1, double a2)
{
 int i_var,i,j,k;
  
 for (i_var = 0; i_var < n_consvar; i_var++)
  {
    for (i = ghosts; i < n_x-ghosts; i++)
    { 
        for (j = ghosts; j <= n_y-ghosts; j++)
        {
            u1[i_var][i][j]= a1 * u1[i_var][i][j] + a2 * u2[i_var][i][j];
        }  
    }
  }
}



void copy_arrays(double*** u1, double*** u2)
{
 int i_var,i,j,k;
  
 for (i_var = 0; i_var < n_consvar; i_var++)
  {
    for (i = ghosts; i < n_x-ghosts; i++)
    { 
        for (j = ghosts; j <= n_y-ghosts; j++)
        {
            u1[i_var][i][j]= u2[i_var][i][j];
        }  
    }
  }
}





void RK4step(double** x, double** y, double*** u0, double*** u1, double*** u2, double*** u_init, double*** res0, double**** flux, double time, double dtime)
{

calc_RHS(x, y, u0, u_init, res0, flux, time);   //f0
advance_solution(u0, u0, res0, dtime);

/*
copy_arrays(u2, u0); 
copy_arrays(u1, u0); //u1 and u2 have a copy of u0


calc_RHS(x, y, u0, u_init, res0, flux, time);   //f0


advance_solution(u1, u2, res0, 0.5*dtime);      //u1 = u2 + dt*f0/2.0
advance_solution(u0, u0, res0, dtime/6.0);     //solution

calc_RHS(x, y, u1, u_init, res0, flux, time + 0.5*dtime);  // f1


advance_solution(u1, u2, res0, 0.5*dtime);    //u2
advance_solution(u0, u0, res0, dtime/3.0);   //solution

calc_RHS(x, y, u1, u_init, res0, flux, time + 0.5*dtime);  //f2


advance_solution(u1, u2, res0, dtime);  //u3
advance_solution(u0, u0, res0, dtime/3.0);   //solution

calc_RHS(x, y, u1, u_init, res0, flux, time + dtime); //f3

advance_solution(u0, u0, res0, dtime/6.0);    //solution
*/


}

















