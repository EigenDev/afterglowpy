// offaxis.h

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

// some physical and mathematical constants
#define PI          3.14159265358979323846
#define v_light     2.99792458e10   // speed of light in cm / s
#define invv_light  3.335640952e-11 // inverse speed of light s / cm
#define m_e         9.1093897e-28   // electron mass in g
#define m_p         1.6726231e-24   // proton mass in g
#define invm_e      1.097768383e27  // inverse electron mass in 1/g
#define invm_p      5.978633202e23  // inverse proton mass in 1/g
#define h_planck    6.6260755e-27   // Planck's constant in erg / s
#define h_bar       1.05457266e-27  // Planck's constant / 2 PI in erg /s
#define k_B         1.380658e-16    // Boltzmann's constant in erg / K
#define e_e         4.803e-10       // electron charge in Gaussian cgs units
#define sigma_T     6.65e-25        // Thomson cross section free electron cm^2
#define cgs2mJy     1e26            // quantity in cgs to mJy
#define mJy2cgs     1e-26           // quantity in mJy to cgs
#define deg2rad     0.017453292     // quantity in degrees to radians
#define rad2deg     57.29577951     // quantity in radians to degrees
#define sec2day     0.000011574     // quantity in seconds to days
#define day2sec     86400           // quantity in days to seconds
#define parsec      3.0857e18       // quantity in parsec to cm
#define Hz2eV       4.13566553853599E-15
#define eV2Hz       2.417991e+14

#define _tophat -1
#define _Gaussian 0
#define _powerlaw 1
#define _Gaussian_core 2 // has a core as well

struct fluxParams
{
    double theta;
    double phi;
    
    double theta_obs;
    double t_obs;
    double nu_obs;
    double d_L;

    double E_iso;
    double n_0;

    double p;
    double epsilon_E;
    double epsilon_B;
    double ksi_N;

    double theta_h;
    
    double current_theta_cone_hi;
    double current_theta_cone_low;
    double theta_obs_cur;

    double Rt0;
    double Rt1;

    double C_BMsqrd;
    double C_STsqrd;

    double t_NR;

    double *t_table;
    double *R_table;
    double *mu_table;
    double *alpha_table;
    int table_entries;
};


double dmin(const double a, const double b);
double get_lfacbetashocksqrd(double a_t_e, double C_BMsqrd, double C_STsqrd);
double get_lfacbetasqrd(double a_t_e, double C_BMsqrd, double C_STsqrd);
double Rintegrand(double a_t_e, void* params);
void make_R_table(struct fluxParams *pars);
// get R from table
double get_R(double a_t_e, double Rt0, double Rt1, double *R_table, 
                double *t_table, double *alpha_table, int table_entries);
void make_mu_table(struct fluxParams *pars);
double get_t_e(double a_mu, double t_obs, double *mu_table, double *t_table,
                int table_entries);
double theta_integrand(double a_theta, void* params); // inner integral
double phi_integrand(double a_phi, void* params); // outer integral
double flux(struct fluxParams *pars); // determine flux for a given t_obs

double flux_cone(double t_obs, double nu_obs, double E_iso, double theta_h,
                    double theta_cone_low, double theta_cone_hi,
                    struct fluxParams *pars);
void lc_tophat(double *t, double *nu_obs, double *F, int Nt,
                double E_iso, double theta_h, struct fluxParams *pars);
void lc_powerlaw(double *t, double *nu_obs, double *F, int Nt,
                    double E_iso_core, double theta_h_core, 
                    double theta_h_wing, double beta,
                    double *theta_c_arr, double *E_iso_arr,
                    int res_cones, struct fluxParams *pars);
void lc_Gaussian(double *t, double *nu_obs, double *F, int Nt,
                        double E_iso_core, 
                        double theta_h_core, double theta_h_wing,
                        double *theta_c_arr, double *E_iso_arr,
                        int res_cones, struct fluxParams *pars);
void lc_GaussianCore(double *t, double *nu_obs, double *F, int Nt,
                        double E_iso_core,
                        double theta_h_core, double theta_h_wing,
                        double *theta_c_arr, double *E_iso_arr,
                        int res_cones, struct fluxParams *pars);
void calc_flux_density(int jet_type, double *t, double *nu, double *Fnu, int N,
                            double theta_obs, double E_iso_core,
                            double theta_h_core, double theta_h_wing, 
                            double n_0, double p, double epsilon_E,
                            double epsilon_B, double ksi_N, double d_L);

void setup_fluxParams(struct fluxParams *pars,
                    double d_L,
                    double theta_obs,
                    double n_0,
                    double p,
                    double epsilon_E,
                    double epsilon_B, 
                    double ksi_N,
                    double Rt0, double Rt1,
        int table_entries);
void set_jet_params(struct fluxParams *pars, double E_iso, double theta_h);
void set_obs_params(struct fluxParams *pars, double t_obs, double nu_obs,
                        double theta_obs_cur, double current_theta_cone_hi, 
                        double current_theta_cone_low);
void free_fluxParams(struct fluxParams *pars);