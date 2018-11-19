#include <iostream>
#include <cmath>
#include "constants.h"
#include "mpi.h"

void indto3d(int index, int *index_x, int *index_y)
{
*index_x = index % nd_x;
*index_y = (index % (nd_x * nd_y)) / nd_x;
}

int neighbour(int cpuindex, int direction, int leftorright)
{
int index;

int index_x, index_y;
int index_x_l, index_x_r;
int index_y_l, index_y_r;

indto3d(cpuindex, &index_x, &index_y);

if (leftorright == 1) {   //if right boundary, next neighbour
switch (direction) {
 case(0):  index_x_r = index_x + 1;
           index_x_r = index_x_r > nd_x-1 ? 0 : index_x_r;
	   index = index_x_r + index_y * nd_x;
           break;
 case(1):  index_y_r = index_y + 1;
           index_y_r = index_y_r > nd_y-1 ? 0 : index_y_r;
           index = index_x + index_y_r * nd_x;
	   break;
	   }
}

if (leftorright == 0) {   //if left boundary, previous neighbour
switch (direction) {
 case(0):  index_x_l = index_x - 1;
           index_x_l = index_x_l < 0 ? nd_x-1 : index_x_l;
           index = index_x_l + index_y*nd_x;
	   break;
 case(1):  index_y_l = index_y - 1;
           index_y_l = index_y_l < 0 ? nd_y-1 : index_y_l;
           index = index_x + index_y_l * nd_x;
	   break;
	   }		   
}   
	 	    
return(index);
}
 

int isinner(int cpuindex, int direction, int leftorright)
{
int is_inner;
int cpu_max, cpu_ind;
int index_x, index_y;

indto3d(cpuindex, &index_x, &index_y);

switch (direction) {
 case(0): cpu_max = nd_x;
          cpu_ind = index_x;
          break;
 case(1): cpu_max = nd_y;
          cpu_ind = index_y;
          break;
}

is_inner = 1;

if ((cpu_ind == 0) && (leftorright == 0)) {  is_inner = 0; }
if ((cpu_ind == cpu_max - 1) && (leftorright == 1)) {  is_inner = 0; }

return(is_inner);
}
