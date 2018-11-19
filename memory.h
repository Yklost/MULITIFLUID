void allocate_memory3D  (double*** &, double*** &, double*** &, double*** &, double *** &, double *** &, double **** &, double ** &, double ** &);
void deallocate_memory3D(double***,   double***,   double***,   double***,   double ***,   double ***,   double ****,   double **,   double **);
void allocate_2D(double** &, int , int );
void allocate_3D(double*** &, int , int , int );
void allocate_4D(double**** &, int , int , int , int );
void allocate_5D(double***** &, int , int , int , int , int );
void deallocate_2D(double** , int , int );
void deallocate_3D(double*** , int , int , int );
void deallocate_4D(double**** , int , int , int , int );
void deallocate_5D(double***** , int , int , int , int , int );

void ABORT();
