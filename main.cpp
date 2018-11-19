#include <cstdlib>
#include <cmath>
#include <array>
#include <locale>
#include <iostream>

#include <string.h>
#include <fstream>
#include <sstream>

#include "constants.h"
#include "mpi.h"
#include "mpiio.h"
#include "memory.h"
#include "vars.h"
#include "initial_condition.h"
#include "physics.h"
#include "calc_RHS.h"

//FYI:
//double ***x, ***y, ***z, ****u0, ****u1, ****u2, ****res0, ****res1, ****res2, ****u_init, *****flux;

int main(int argc, char** argv) {    


int i,j,k;

double dtime, time, time_0;

int n_cpu, i_cpu;

int i_t, i_t_0;

    
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &n_cpu);
MPI_Comm_rank(MPI_COMM_WORLD, &i_cpu);
if (nd_x*nd_y != n_cpu) {ABORT();}

    allocate_memory3D(u0, u1, u2, u_init, res0, res1, flux, x, y);
    
    if (i_cpu == 0) std::cout << "Memory allocated." << std::endl;
    
    initial_condition(u0, u_init, x, y, i_cpu, &i_t_0, &time_0);

    if (i_cpu == 0) {std::cout << "Box initialised." << std::endl;
                     std::cout << "GLOBAL GRID: " << n_glob_x << " X " << n_glob_y << " X " <<  std::endl;
                     std::cout << "LOCAL GRID: " << n_x << " X " << n_y << " X " << std::endl;
                     std::cout << "Doing " << n_t << " steps." << std::endl;
		     std::cout << "Starting iteration " << i_t_0 << std::endl;
		     std::cout << "Starting time " << time_0 << std::endl;}
    
    time=time_0; 

    for (i_t = i_t_0; i_t < n_t; i_t++)
    {
        
	dtime = timestep(u0);

        if (i_t % save_freq == 0) save_solution(i_cpu, i_t, time, dtime, x, y, u0);

//        con_to_prim(u0);

        if (i_cpu == 0) printf("t=%15.8f     dt=%10.8f     i_t=%8d \n",time,dtime,i_t);

        RK4step(x, y, u0, u1, u2, u_init, res0, flux, time, dtime);
     
        time=time+dtime; //physical time
	
    }
    
    
    
    deallocate_memory3D(u0, u1, u2, u_init, res0, res1, flux, x, y);
    MPI_Finalize();  
    if (i_cpu == 0) std::cout << "Simulation completed." << std::endl;        
    return 0;
}

