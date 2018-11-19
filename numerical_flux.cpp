#include <iostream> 
#include <cmath>
#include "constants.h"
#include "physics.h"

double limiter(double r) 
{
return (1.50*(r*r+r)/(r*r+r+1.0)); //ospre
//return(fmax(0,fmin(1,r))); //minmod
//return(fmax(0,fmax(fmin(2.0*r,1.0),fmin(r,2.0)))); //superbee - unstable for stiff problems
}

void numerical_flux(int dir, double* ull, double* ul, double* uu, double* ur, double* flux)
{ 
    double limiter(double r);
    
    int k;

//    double ru[n_consvar]; 
//    double rl[n_consvar];
//    double ulmh[n_consvar];
//    double urmh[n_consvar];
//    double fll[n_consvar];
//    double flr[n_consvar];
 
 
    double ru[n_var]; 
    double rl[n_var];
    double ulmh[n_var];
    double urmh[n_var];
    double fll[n_var];
    double flr[n_var];
 
   
    double nom1, denom1, nom2, denom2;
    double aaa;
     
//    for (k = 0; k < n_consvar; k++) 
    for (k = 0; k < n_var; k++)
    {
        nom1 = uu[k] - ul[k]; //takes flux from i - (i-1) :
        denom1 = ur[k] - uu[k]; //i+1 - i 

        nom2 = ul[k] - ull[k];
        denom2 = uu[k] - ul[k];

        if (fabs(nom1) < SMALL) { nom1=0.0; denom1=1.0;}
	if (fabs(denom1) < SMALL) {
           if (nom1 > SMALL)  { nom1=BIG;  denom1=1.0; }
           if (nom1 < -SMALL) { nom1=-BIG; denom1=1.0; }
	}
	 
        if (fabs(nom2) < SMALL){ nom2=0.0; denom2=1.0;}
	if (fabs(denom2) < SMALL) {
           if (nom2 > SMALL)  { nom2=BIG;  denom2=1.0; }
           if (nom2 < -SMALL) { nom2=-BIG; denom2=1.0; }
	}
	
	ru[k] = nom1 / denom1; 
        rl[k] = nom2 / denom2; 

	ulmh[k]=ul[k]+0.50*limiter(rl[k])*(uu[k]-ul[k]); 
	urmh[k]=uu[k]-0.50*limiter(ru[k])*(ur[k]-uu[k]);
    }

    aaa = astar(ulmh, urmh, dir);

    get_pflux(dir,urmh,flr);
    get_pflux(dir,ulmh,fll);

    for (k = 0; k < n_consvar; k++) flux[k] = 0.50*(flr[k]+fll[k]-aaa*(urmh[k]-ulmh[k])); 


}

