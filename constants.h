const double PI=3.141592653589793;
const double SMALL = 1.0e-100;
const double BIG   = 1.0/SMALL;

const double me = 9.109e-28 ; //g
const double mi = 1.6726e-24 ; //g
const double e = 4.8032e-10; //statC

const double k_b=1.381e-16 ; //erg/K

const double me_mi = me/mi;

const double III_erg = 10.e0 * 1.6e-12 ; //eV to erg, first ionisation potential
const double III_eV = 10.e0;

const double eVtoerg = 1.6e-12;
const double eVtoerg2 = (1.0/6.2415e11)*(1.0/6.2415e11);


const double alpha_0 = 0.5;
const double beta_0 = 0.7;
const double gamma_0 = 3.2;


const int n_var=37; 
const int n_consvar = 14; 
const int n_dim = 2; //number dimensions

const int n_glob_x = 256;  // global grid x
const int n_glob_y = 256;  // global grid y


const int ghosts = 2; // number of ghost cells

const int nd_x = 2;  // number of subdomains in x
const int nd_y = 2;  // number of subdomains in y


const int n_x = n_glob_x/nd_x + 2*ghosts; //whole subdomain
const int n_y = n_glob_y/nd_y + 2*ghosts; 


const int n_act_x = n_glob_x/nd_x;  //active part of the subdomain
const int n_act_y = n_glob_y/nd_y;  


const int n_t = 1000000;

const int rho_i_ = 0;
const int m_x_i_ = 1;
const int m_y_i_ = 2;
const int m_z_i_ = 3;
const int ene_i_ = 4;

const int e_x_ = 5;
const int e_y_ = 6;
const int e_z_ = 7;

const int rho_n_ = 8;
const int m_x_n_ = 9;
const int m_y_n_ = 10;
const int m_z_n_ = 11;
const int ene_n_ = 12;

const int ene_e_ = 13;

const int v_x_i_ = 14;
const int v_y_i_ = 15;
const int v_z_i_ = 16;
const int pre_i_ = 17;

const int v_x_n_ = 18;
const int v_y_n_ = 19;
const int v_z_n_ = 20;
const int pre_n_ = 21;

const int pre_e_ = 22;

const int k_ion_ = 23;
const int k_rec_ = 24;
const int a_ni_ = 25;
const int a_ne_ = 26;
const int a_ei_ = 27;

const int QQQ_  = 28;
const int Q_N_  = 29;

const int q_x_ = 30;
const int q_y_ = 31;
const int q_z_ = 32;

const int tem_i_ = 33;
const int tem_n_ = 34;
const int tem_e_ = 35;

const int vidte_ = 36;



const double ggg = -1.0;//-1.0; //gravity
const double gam = 5.0/3.0; /// 


const double start_x = 0.0;
const double start_y = 8e7;

const double end_x = 2e8;
const double end_y = 1.5e8;


const double dx = (end_x - start_x) / (1.0*n_glob_x);
const double dy = (end_y - start_y) / (1.0*n_glob_y);


const int initialtype=100;

const int save_freq=1000;

const double CFL=0.00001;

   
const double eta=0.1; //magnetic resistivity

const double nu_visc = 0.01;


//OUTPUT
const int save_grid = 1;
const int save_prim = 1;

//INPUT
const int load_grid = 1;


