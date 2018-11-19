double fast_cfl(double* );
double astar(double* , double* , int );
void calc_grav_sources(double** , double** , double*** , double*** , double );
void calc_coll_sources(double** , double** , double*** , double*** , double );
void get_collision_rates(double*** );
void con_to_prim(double*** );
void prim_to_con(double*** );
double tem(double, double);
//void calc_visc_sources(double** , double** , double*** , double*** , double );

double timestep(double*** );

void get_pflux(int, double*, double*);
