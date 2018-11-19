#include<iostream>
#include<cmath>
#include "constants.h"
#include "mpi.h"

void allocate_2D(double** &var, int n_x, int n_y)
{
   int i;
   var = new double*[n_x];
   for (i = 0; i < n_x; i++) var[i] = new double[n_y];
}

void allocate_3D(double*** &var, int n_x, int n_y, int n_z)
{
   int i,j;
   var = new double**[n_x];
   for (i = 0; i < n_x; i++) allocate_2D(var[i],n_y,n_z);   
}

void allocate_4D(double**** &var, int n_v, int n_x, int n_y, int n_z)
{
   int i,j,k;
   var = new double***[n_v];
   for (i = 0; i < n_v; i++) allocate_3D(var[i],n_x,n_y,n_z);
}

void allocate_5D(double***** &var, int n_d, int n_v, int n_x, int n_y, int n_z)
{
   int i,j,k,l;
   var = new double****[n_d];
   for (i = 0; i < n_d; i++) allocate_4D(var[i],n_v,n_x,n_y,n_z);  
}

void deallocate_2D(double** var, int n_x, int n_y)
{
   int i;
   for (i = 0; i < n_x; i++) delete [] var[i];
   delete [] var;
}
  
void deallocate_3D(double*** var, int n_x, int n_y, int n_z)
{
   int i,j;
   for (i = 0; i < n_x; i++) deallocate_2D(var[i],n_y,n_z);
   delete [] var;
}
	  
void deallocate_4D(double**** var, int n_v, int n_x, int n_y, int n_z)
{
   int i,j,k;
   for (i = 0; i < n_v; i++) deallocate_3D(var[i],n_x,n_y,n_z);
   delete [] var;   
}

void deallocate_5D(double***** var, int n_d, int n_v, int n_x, int n_y, int n_z)
{
   int i,j,k,l;
   for (i = 0; i < n_d; i++) deallocate_4D(var[i],n_v,n_x,n_y,n_z);
   delete [] var;
}



               //(u0,             u1,            u2,             u_init,             res0,             res1,             flux, x, y, z);
void allocate_memory3D(double*** &u0, double*** &u1, double*** &u2, double*** &u_init, double*** &res0, double*** &res1, double**** &flux, 
                       double** &x, double** &y) 
{

allocate_3D(u_init, n_var, n_x, n_y);
allocate_3D(u0,     n_var, n_x, n_y);
allocate_3D(u1,     n_var, n_x, n_y);
allocate_3D(u2,     n_var, n_x, n_y);
allocate_3D(res0,   n_var, n_x, n_y);
allocate_3D(res1,   n_var, n_x, n_y);

allocate_4D(flux,   n_dim, n_var, n_x, n_y);

allocate_2D(x,      n_x, n_y);
allocate_2D(y,      n_x, n_y);

}



void deallocate_memory3D(double*** u0, double*** u1, double*** u2, double*** u_init, double*** res0, double*** res1, double**** flux, 
                        double** x, double** y) 
{

deallocate_3D(u_init, n_var, n_x, n_y);
deallocate_3D(u0,     n_var, n_x, n_y);
deallocate_3D(u1,     n_var, n_x, n_y);
deallocate_3D(u2,     n_var, n_x, n_y);
deallocate_3D(res0,   n_var, n_x, n_y);
deallocate_3D(res1,   n_var, n_x, n_y);

deallocate_4D(flux,   n_dim, n_var, n_x, n_y);

deallocate_2D(x,      n_x, n_y);
deallocate_2D(y,      n_x, n_y);

}


void ABORT()
{
std::cout << "NEED TO ABORT!" << std::endl;
MPI_Finalize();
abort();
}
