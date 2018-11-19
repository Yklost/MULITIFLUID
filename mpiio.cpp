#include <iostream> //needs to be included!
#include <cmath>
#include "constants.h"
#include <string.h>
#include <fstream>
#include <sstream>
#include "mpi.h"
#include "mpigeom.h"
#include "mpiio.h"



void save_solution(int i_cpu, int i_t, double time, double dtime, double** x, double** y, double*** u)
{
     if (i_cpu == 0) save_info(i_t, time, dtime);
     
     save_binary(i_t, x, y, u);
}



void save_info(int i_t, double time, double dtime)
{
//Here we save two files: current and solution.log. Current contains the state of the latest saved snapshot only; solution.log is appended.
FILE *cur_file;
FILE *log_file;

cur_file = fopen("current", "w");
log_file = fopen("solution.log", "a");

fprintf(cur_file,"%8d %15.8f %10.8f \n",i_t, time, dtime);
fprintf(log_file,"i_t=%8d     t=%15.8f     dt=%10.8f \n", i_t, time, dtime);

fclose(cur_file);
fclose(log_file);

}




void load_info(int i_t, double time, double dtime)
{





}




void save_binary(int i_t, double** x, double** y, double*** u)
{  

    std::string str1;
    int len, pos;
    int ioerr;

    MPI_File outfile;
    MPI_Status status;
    MPI_Offset offset;
    int index;
    int index_x, index_y;

    double *obuffer;
    
    int i,j;
    int i_var;
    int bufsize;
    
    int n_save_var;


    int glodim[n_dim];
    int locdim[n_dim];
    int start[n_dim];
    
    MPI_Datatype bufd;
    MPI_Offset filesize;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &index);
    
    glodim[0] = n_glob_x;
    glodim[1] = n_glob_y;

    locdim[0] = n_glob_x/nd_x;
    locdim[1] = n_glob_y/nd_y;

    
    indto3d(index, &index_x, &index_y);
    
    start[0] = locdim[0] * index_x;
    start[1] = locdim[1] * index_y;
    

    bufsize = (n_x - 2 * ghosts) * (n_y - 2 * ghosts);

    
    obuffer = new double[bufsize];
    

    str1 = std::to_string(i_t);

    len = str1.length();

    std::string num = std::string(7, '0');

    std::string filename = "./data_" + num.replace(7 - len, len, str1);



    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Type_create_subarray(n_dim, glodim, locdim, start, MPI_ORDER_FORTRAN, MPI_DOUBLE, &bufd);
    
    MPI_Type_commit(&bufd);

    MPI_File_open(MPI_COMM_WORLD, filename.c_str() , MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfile);

    if (save_grid == 1) {
    //save x coordinates 

	for (i = ghosts; i < n_x - ghosts; i++)
            {
	     for (j = ghosts; j < n_y - ghosts; j++) 
	             obuffer[(i - ghosts) + (n_x - 2*ghosts) * (j - ghosts)] = x[i][j];
	}

	MPI_File_get_size(outfile, &filesize);

	MPI_File_set_view(outfile, filesize, MPI_DOUBLE, bufd, "native", MPI_INFO_NULL);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File_write_all(outfile, obuffer, bufsize, MPI_DOUBLE, &status);

	MPI_Barrier(MPI_COMM_WORLD);


    //save y coordinates

	for (i = ghosts; i < n_x - ghosts; i++)
            {
	     for (j = ghosts; j < n_y - ghosts; j++) 
	             obuffer[(i - ghosts) + (n_x - 2*ghosts) * (j - ghosts)] = y[i][j];
	}

	MPI_File_get_size(outfile, &filesize);

	MPI_File_set_view(outfile, filesize, MPI_DOUBLE, bufd, "native", MPI_INFO_NULL);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File_write_all(outfile, obuffer, bufsize, MPI_DOUBLE, &status);

	MPI_Barrier(MPI_COMM_WORLD);

    }    
    
    n_save_var = save_prim == 1 ? n_var : n_consvar;
        
//save data block    

    for (i_var = 0; i_var < n_save_var; i_var++)
	{
	for (i = ghosts; i < n_x - ghosts; i++)
            {
	     for (j = ghosts; j < n_y - ghosts; j++) 
	             obuffer[(i - ghosts) + (n_x - 2*ghosts) * (j - ghosts)] = u[i_var][i][j];
	}

    MPI_File_get_size(outfile, &filesize);

    MPI_File_set_view(outfile, filesize, MPI_DOUBLE, bufd, "native", MPI_INFO_NULL);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_File_write_all(outfile, obuffer, bufsize, MPI_DOUBLE, &status);

    MPI_Barrier(MPI_COMM_WORLD);

    }

    MPI_File_close(&outfile);
    MPI_Type_free(&bufd);
   
    delete [] obuffer;
    

}





//NOTE HERE: THERE ARE ABSOLUTELY NO MEANS TO CHECK VALIDITY OF THE LOADED DATA
void load_binary(int i_t, double** x, double** y, double*** u)
{  

    std::string str1;
    int len, pos;
    int ioerr;

    MPI_File infile;
    MPI_Status status;
    MPI_Offset offset;
    int index;
    int index_x, index_y;

    double *ibuffer;
    
    int i,j;
    int i_var;
    int bufsize;
    
    int n_load_var;


    int glodim[n_dim];
    int locdim[n_dim];
    int start[n_dim];
    
    MPI_Datatype bufd;
    MPI_Offset filesize;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &index);
    
    glodim[0] = n_glob_x;
    glodim[1] = n_glob_y;

    locdim[0] = n_glob_x/nd_x;
    locdim[1] = n_glob_y/nd_y;

    
    indto3d(index, &index_x, &index_y);
    
    start[0] = locdim[0] * index_x;
    start[1] = locdim[1] * index_y;
    

    bufsize = (n_x - 2 * ghosts) * (n_y - 2 * ghosts);

    
    ibuffer = new double[bufsize];
    

    str1 = std::to_string(i_t);

    len = str1.length();

    std::string num = std::string(7, '0');

    std::string filename = "./data_" + num.replace(7 - len, len, str1);




    MPI_Type_create_subarray(n_dim, glodim, locdim, start, MPI_ORDER_FORTRAN, MPI_DOUBLE, &bufd);
    
    MPI_Type_commit(&bufd);

    MPI_File_open(MPI_COMM_WORLD, filename.c_str() , MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);

    if (load_grid == 1) {
    //load x coordinates 



	MPI_File_get_size(infile, &filesize);

	MPI_File_set_view(infile, filesize, MPI_DOUBLE, bufd, "native", MPI_INFO_NULL);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File_read_all(infile, ibuffer, bufsize, MPI_DOUBLE, &status);

	for (i = ghosts; i < n_x - ghosts; i++)
            {
	     for (j = ghosts; j < n_y - ghosts; j++) 
	             x[i][j] = ibuffer[(i - ghosts) + (n_x - 2*ghosts) * (j - ghosts)];
	}

	MPI_Barrier(MPI_COMM_WORLD);


    //load y coordinates


	MPI_File_get_size(infile, &filesize);

	MPI_File_set_view(infile, filesize, MPI_DOUBLE, bufd, "native", MPI_INFO_NULL);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_File_read_all(infile, ibuffer, bufsize, MPI_DOUBLE, &status);
	
	for (i = ghosts; i < n_x - ghosts; i++)
            {
	     for (j = ghosts; j < n_y - ghosts; j++) 
	             y[i][j] = ibuffer[(i - ghosts) + (n_x - 2*ghosts) * (j - ghosts)];
	}

	MPI_Barrier(MPI_COMM_WORLD);

    }    
        
//load data block    

    MPI_File_get_size(infile, &filesize);

    MPI_File_set_view(infile, filesize, MPI_DOUBLE, bufd, "native", MPI_INFO_NULL);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_File_read_all(infile, ibuffer, bufsize, MPI_DOUBLE, &status);
    
    for (i_var = 0; i_var < n_consvar; i_var++)
    {
    for (i = ghosts; i < n_x - ghosts; i++)
        {
	 for (j = ghosts; j < n_y - ghosts; j++) 
	          u[i_var][i][j] = ibuffer[(i - ghosts) + (n_x - 2*ghosts) * (j - ghosts)];
    }


    MPI_Barrier(MPI_COMM_WORLD);

    }

    MPI_File_close(&infile);
    MPI_Type_free(&bufd);
   
    delete [] ibuffer;
    

}










