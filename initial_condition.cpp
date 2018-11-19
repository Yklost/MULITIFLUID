#include <iostream> 
#include <cmath>
#include "constants.h"
#include "mpigeom.h"
#include "physics.h"
#include "collisions.h"


void initial_condition(double*** u, double*** u_init, double** x, double** y, int cpuindex, int* i_t_0, double* time_0)
{ 
int i,j,k;
int i_var;
double rad, rad_0, vabs, b0;

int index_x, index_y;
double x_loc_start, y_loc_start, z_loc_start;

*i_t_0 = 0; //start iterations from iteration 0 and time 0 unless initialtype=255, which is specified later
*time_0 = 0.0;

indto3d(cpuindex, &index_x, &index_y);
	
x_loc_start = start_x + (index_x*n_act_x - ghosts)*dx;
y_loc_start = start_y + (index_y*n_act_y - ghosts)*dy;

	
for (i = 0; i < n_x; i++)
{ 
   for (j = 0; j < n_y; j++)
   {

           x[i][j] = x_loc_start + i*dx;
	   y[i][j] = y_loc_start + j*dy;

   }
}
    
    
rad_0 = (end_x - start_x) * 0.25;


//Sod shock
if (initialtype == 1) {

            for (i = 0; i < n_x; i++)
                { 
                for (j = 0; j < n_y; j++)
                    {
		    
		    rad = sqrt((x[i][j] - 0.5 * (end_x - start_x))*(x[i][j] - 0.5 * (end_x - start_x)) + (y[i][j] - 0.5 * (end_y - start_y))*(y[i][j] - 0.5 * (end_y - start_y)));
		    

//                            if (x[i][j] < 0.5 * (end_x - start_x)) //0.25)
                            if (rad < rad_0)
                            {
                            u[rho_i_][i][j]=5.74e11 * mi;
                            u[pre_i_][i][j]=0.167;
			    
                            u[rho_n_][i][j]=1.1e10 * (mi+me); 
                            u[pre_n_][i][j]=0.167;
			    
			    u[pre_e_][i][j]=0.167;
			    		    
			    u[v_x_i_][i][j]=0.0;
			    u[v_y_i_][i][j]=0.0;
			    u[v_z_i_][i][j]=0.0;
			    u[v_x_n_][i][j]=0.0;
			    u[v_y_n_][i][j]=0.0;
			    u[v_z_n_][i][j]=0.0;
			    
                            }
			   
//                            if (x[i][j] >= 0.5 * (end_x - start_x)) //0.25)
                            if (rad >= rad_0)
                            {
                            u[rho_i_][i][j]=5.74e11 * mi;
                            u[pre_i_][i][j]=0.167;
			    
                            u[rho_n_][i][j]=0.5* 1.1e10 * (mi+me);
                            u[pre_n_][i][j]=0.4* 0.167;
			    
			    u[pre_e_][i][j]=0.167;
			    
			    u[v_x_i_][i][j]=0.0;
			    u[v_y_i_][i][j]=0.0;
			    u[v_z_i_][i][j]=0.0;
			    u[v_x_n_][i][j]=0.0;
			    u[v_y_n_][i][j]=0.0;
			    u[v_z_n_][i][j]=0.0;
			    
			    
                            }
		    }
                }
}

  
  


//Kelvin-Helmholtz    
if (initialtype == 4) {
    
            for (i = 0; i < n_x; i++)
                { 
                for (j = 0 ; j< n_y ; j++)
                    {
                            u[rho_i_][i][j]=1.0; 
                            u[pre_i_][i][j]=2.5;
                            u[rho_n_][i][j]=1.0; 
                            u[pre_n_][i][j]=2.5;
			    u[v_x_i_][i][j]=0.5;
			    u[v_y_i_][i][j]=0.0;
			    u[v_z_i_][i][j]=0.0;	
			    u[v_x_n_][i][j]=0.5;
			    u[v_y_n_][i][j]=0.0;
			    u[v_z_n_][i][j]=0.0;

                            if ((y[i][j] >= 0.4+0.01*sin(2.0*PI*x[i][j])) && (y[i][j] <= 0.6+0.01*sin(2.0*PI*x[i][j])))
                           {
                            u[v_x_i_][i][j]=-0.5;
			    u[rho_i_][i][j]=2.0;
                            u[v_x_n_][i][j]=-0.5;
			    u[rho_n_][i][j]=2.0;
                            }
			    
			    u[pre_e_][i][j]=3.0;
		    }
                }
} 







if (initialtype == 100) {

    int n_z;
    double zmin, zmax, dz;
    double dum;
    int i_in;

    double z;
    int i_z;
    double interp_coef;

    FILE *fp = fopen("initmod.dat","r");
    fscanf(fp, "%i\n", &n_z);

    double *zzz = new double[n_z];
    double *ppp_i = new double[n_z];
    double *ppp_n = new double[n_z];
    double *ppp_e = new double[n_z];
    double *rhr_e = new double[n_z];
    double *rhr_i = new double[n_z];
    double *rhr_n = new double[n_z];
    double *v = new double[256]; //testcase

    fscanf(fp, "%lf %lf %lf \n", &zmin, &zmax, &dz);
    

    if (cpuindex == 0) {
    std::cout << "EXTERNAL RP FOR INE TABLE INITIAL CONDITION" << std::endl;
    std::cout << "N_Z=" << n_z << std::endl;
    std::cout << "ZMIN=" << zmin << " ZMAX=" << zmax << " DZ=" << dz << std::endl;
    }

    
    for (i_in = 0; i_in < n_z; i_in++) {
                                       fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf \n",
				              &zzz[i_in], &rhr_e[i_in], &rhr_i[i_in], &rhr_n[i_in], &ppp_e[i_in], &ppp_i[i_in], &ppp_n[i_in]);
	
			               }
    fclose(fp);

    
    for (k = 0; k < 256; k++){v[k] = 255-k;}//testcase

    for (i = 0; i < n_x; i++)
	{ 
	for (j = 0; j < n_y; j++)
            {
		 i_z = y[i][j] / dz;
		 interp_coef = (y[i][j] - i_z * dz) / dz;

		 u[rho_i_][i][j] = exp(log(rhr_i[i_z]) + (log(rhr_i[i_z+1]) - log(rhr_i[i_z])) * interp_coef);
		 u[rho_n_][i][j] = exp(log(rhr_n[i_z]) + (log(rhr_n[i_z+1]) - log(rhr_n[i_z])) * interp_coef);    
		 u[pre_i_][i][j] = exp(log(ppp_i[i_z]) + (log(ppp_i[i_z+1]) - log(ppp_i[i_z])) * interp_coef);
		 u[pre_n_][i][j] = exp(log(ppp_n[i_z]) + (log(ppp_n[i_z+1]) - log(ppp_n[i_z])) * interp_coef);     
		 u[pre_e_][i][j] = exp(log(ppp_e[i_z]) + (log(ppp_e[i_z+1]) - log(ppp_e[i_z])) * interp_coef); 
                 //u[v_y_i_][i][j] = exp(log(v[i_z]) + (log(v[i_z+1])-log(v[i_z]))*interp_coef);   //testcase 
	     }
	}


 
    delete [] zzz;
    delete [] ppp_i;
    delete [] ppp_n;
    delete [] ppp_e;
    delete [] rhr_i;
    delete [] rhr_n;
    delete [] rhr_e;
    delete [] v;  //testcase

}



//INITIAL MODEL FROM FILE
if (initialtype == 255) {



}





//END INITIAL MODELS


prim_to_con(u);
get_collision_rates(x, y, u);


//copy ICs to u_init
for (i_var = 0; i_var < n_var; i_var++)
    {
         for (i = 0; i < n_x; i++)
         { 
             for (j = 0 ; j< n_y ; j++)
                 {
			  u_init[i_var][i][j]=u[i_var][i][j];
		 }
	 }
    }


}
