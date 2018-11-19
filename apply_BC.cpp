#include <iostream>
#include <cmath>
#include "constants.h"
#include "mpi.h"
#include "mpigeom.h"

void apply_BC(double*** u, double*** u_init, double** x, double** y)
{
    int i,j;
    
    double radius;

    int bctype[n_var][n_dim][2];
       
    int i_var, ib, idim, i_gc;
    
    int nbp1;
    int ind1;
    
    int bufsize;
    
    int i_x_in, i_y_in, i_x_p, i_y_p, i_x, i_y;
    
    int i_cpu, i_neighbour_to, i_neighbour_fr;
    int b_type;  
    int thisbc;  
    
    int tag_send, tag_recv;
    
    double *sendbuf, *recvbuf;
    
    MPI_Request req_send, req_recv;
    MPI_Status status;
 
    int ierr;
 
 
 //move boundary types to somewhere somehow
 for (i_var = 0; i_var < n_var; i_var++) {
    bctype[i_var][0][0]=0;
    bctype[i_var][0][1]=0;
    bctype[i_var][1][0]=1;
    bctype[i_var][1][1]=1;
    }
 
    
    
    MPI_Comm_rank(MPI_COMM_WORLD, &i_cpu);

    for (idim = 0; idim < n_dim; idim++) {  
        for (ib = 0; ib <= 1; ib++) {     

             if (idim == 0) nbp1=n_y;
		       
	     if (idim == 1) nbp1=n_x;

	     
	     bufsize = n_var * nbp1 * ghosts; //change to n_consvar later if possible
		
	     sendbuf = new double [bufsize]; 
	     recvbuf = new double [bufsize]; 

	     for (ind1 = 0; ind1 < nbp1; ind1++)
              {
               for (i_gc = 0; i_gc < ghosts; i_gc++)
        	   {
                    if ((idim == 0) && (ib == 0)) { //left x boundary
			i_x = ghosts + i_gc;
			i_y = ind1;
                	}

		    if ((idim == 0) && (ib == 1)) { //right x boundary
			i_x = n_x - 2*ghosts + i_gc;
			i_y = ind1;
			}

		    if ((idim == 1) && (ib == 0)) { //bottom y boundary
			i_x = ind1;
			i_y = ghosts + i_gc;
			}

		    if ((idim == 1) && (ib == 1)) { //top y boundary
			i_x = ind1;
			i_y = n_y - 2*ghosts + i_gc;
			}
		
        	   for (i_var = 0; i_var < n_var; i_var++) {
                	sendbuf[i_var + n_var*(i_gc + ghosts*ind1)] = u[i_var][i_x][i_y];
		       }   
		   } //i_gc
	       } //ind1


               i_neighbour_to = neighbour(i_cpu, idim, ib);    //which cpu send to
               i_neighbour_fr = neighbour(i_cpu, idim, 1-ib);  //which cpu recv from
	       tag_send = 1;//2*i_neighbour+ib;
	       tag_recv = 1;//2*i_cpu+ib;

//std::cout << "CPU=" << i_cpu << " SENDING DIRECTION=" << idim << " BOUNDARY=" << ib << " FROM CPU=" << i_cpu << " TO CPU=" << i_neighbour_to << " TAGSEND=" << tag_send  << std::endl;
//std::cout << "CPU=" << i_cpu << " RECVING DIRECTION=" << idim << " BOUNDARY=" << ib << " FROM CPU=" << i_neighbour_fr << " TO CPU=" << i_cpu << " TAGRECV=" << tag_recv  << std::endl;

               MPI_Sendrecv(sendbuf, bufsize, MPI_DOUBLE, i_neighbour_to, tag_send, recvbuf, bufsize, MPI_DOUBLE, i_neighbour_fr, tag_recv, 
				                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

 				      
	       for (ind1 = 0; ind1 < nbp1; ind1++)
                {
                 for (i_gc = 0; i_gc < ghosts; i_gc++)
                     {
                      if ((idim == 0) && (ib == 1)) { //left x boundary
			  i_x = 0 + i_gc;
			  i_y = ind1;
                          }

		      if ((idim == 0) && (ib == 0)) { //right x boundary
			  i_x = n_x - ghosts + i_gc;
			  i_y = ind1;
			  }

		      if ((idim == 1) && (ib == 1)) { //bottom y boundary
			  i_x = ind1;
			  i_y = 0 + i_gc;
			  }

		      if ((idim == 1) && (ib == 0)) { //top y boundary
			  i_x = ind1;
			  i_y = n_y - ghosts + i_gc;
			  }


                     for (i_var = 0; i_var < n_var; i_var++) {
                          u[i_var][i_x][i_y]=recvbuf[i_var + n_var*(i_gc + ghosts*ind1)];
		         }   
		     } //i_gc
		 } //ind1	      

		 delete [] sendbuf;
		 delete [] recvbuf;	       
		
                 MPI_Barrier(MPI_COMM_WORLD);	
		
	    
    } //ib

} //idim

//end MPI communication			
			
			
			
			
for (idim = 0; idim < n_dim; idim++) {  
        for (ib = 0; ib <= 1; ib++) { 
  
        if (isinner(i_cpu, idim, ib) == 0) {
	
             if (idim == 0) { nbp1=n_y; }
		       
	     if (idim == 1) { nbp1=n_x; }
		

             for (ind1 = 0; ind1 < nbp1; ind1++)
                 {
                  for (i_gc = 0; i_gc < ghosts; i_gc++)
                      {
                       if ((idim == 0) && (ib == 0)) { //left x boundary
		           i_x_in = ghosts;
			   i_y_in = ind1;
			   i_x = 0 + i_gc;
			   i_y = ind1;
			   i_x_p = n_x - 2*ghosts + i_gc;
			   i_y_p = ind1;
                           }
		       
		       if ((idim == 0) && (ib == 1)) { //right x boundary
		           i_x_in = n_x - ghosts - 1;
			   i_y_in = ind1;
			   i_x = n_x - ghosts + i_gc;
			   i_y = ind1;
			   i_x_p = ghosts + i_gc;
			   i_y_p = ind1;
			   }
			   
		       if ((idim == 1) && (ib == 0)) { //bottom y boundary
		           i_x_in = ind1;
			   i_y_in = ghosts;
			   i_x = ind1;
			   i_y = 0 + i_gc;
			   i_x_p = ind1;
			   i_y_p = n_y - 2*ghosts + i_gc;
			   }
			   
		       if ((idim == 1) && (ib == 1)) { //top y boundary
		           i_x_in = ind1;
			   i_y_in = n_y - ghosts - 1;
			   i_x = ind1;
			   i_y = n_y - ghosts + i_gc;
			   i_x_p = ind1;
			   i_y_p = ghosts + i_gc;
			   }


                      for (i_var = 0; i_var < n_var; i_var++) {

                	       thisbc = bctype[0][idim][ib];

                	       switch (thisbc) {
        		       case(0):    //constant extrapolation from i_x_b, i_y_b
				      u[i_var][i_x][i_y] = u[i_var][i_x_in][i_y_in];
				      break;
        		       case(1):    //fixed, copy initial values
        			      u[i_var][i_x][i_y] = u_init[i_var][i_x][i_y];
				      break;
        		       case(2):    //zeros
        			      u[i_var][i_x][i_y] = 0.0;
				      break;
			       case(3):    //periodic is already done by MPI communication; all other boundaries are local.
        			//      u[i_var][i_x][i_y] = u[i_var][i_x_p][i_y_p];
				      break;
                               case(4):
                                      if (i_var == 0 || i_var == 4 || i_var == 15) {
                                        u[i_var][i_x][i_y] = u[i_var][i_x_in][i_y_in] + (y[0][1]-y[0][0])*1000.0;  
                                      } else {	
                                        u[i_var][i_x][i_y] = u[i_var][i_x_in][i_y_in];
                                      }  
	        		}
                           }	// i_var 	     
                   } //i_gc
		} //ind1
	
      } //isinner			    
			    
    } //ib

} //idim

		        

MPI_Barrier(MPI_COMM_WORLD);
		  

}

