#ifndef OFFAXIS_STRUCT_HPP
#define OFFAXIS_STRUCT_HPP

// offaxis.h

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "integrate.hpp"

#define MSG_LEN 4096
#define DUMP_MSG_LEN_MAX 16384 // overkill: 200 lines * 80c per line = 16000

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// some physical and mathematical constants
constexpr double pi = 3.14159265358979323846;
constexpr double v_light = 2.99792458e10;      // speed of light in cm / s
constexpr double invv_light = 3.335640952e-11; // inverse speed of light s / cm
constexpr double m_e = 9.1093897e-28;          // electron mass in g
constexpr double m_p = 1.6726231e-24;          // proton mass in g
constexpr double invm_e = 1.097768383e27;      // inverse electron mass in 1/g
constexpr double invm_p = 5.978633202e23;      // inverse proton mass in 1/g
constexpr double h_planck = 6.6260755e-27;     // Planck's constant in erg / s
constexpr double h_bar = 1.05457266e-27;       // Planck's constant / 2 PI in erg /s
constexpr double k_B = 1.380658e-16;           // Boltzmann's constant in erg / K
constexpr double e_e = 4.803e-10;              // electron charge in Gaussian cgs units
constexpr double sigma_T = 6.65e-25;           // Thomson cross section free electron cm^2
constexpr double cgs2mJy = 1e26;               // quantity in cgs to mJy
constexpr double mJy2cgs = 1e-26;              // quantity in mJy to cgs
constexpr double deg2rad = 0.017453292;        // quantity in degrees to radians
constexpr double rad2deg = 57.29577951;        // quantity in radians to degrees
constexpr double sec2day = 0.000011574;        // quantity in seconds to days
constexpr double day2sec = 86400;              // quantity in days to seconds
constexpr double parsec = 3.0857e18;           // quantity in parsec to cm
constexpr double Hz2eV = 4.13566553853599E-15;
constexpr double eV2Hz = 2.417991e+14;

enum class EVOL_TYPE
{
    CONE,
    TOPHAT,
    GAUSSIAN,
    POWERLAWCORE,
    GAUSSIANCORE,
    SPHERICAL,
    EXPONENTIAL,
    TWOCOMPONENT,
    EXPONENTIAL2,
    RING,
};

#define _cone -2
#define _tophat -1
#define _Gaussian 0
#define _powerlaw_core 1 // has a core as well
#define _Gaussian_core 2 // has a core as well
#define _spherical 3
#define _powerlaw 4
#define _exponential 5
#define _twocomponent 6
#define _exponential2 7
#define _ring 8

#define IC_COOLING_FLAG 1
#define EPS_E_BAR_FLAG 2
#define SSA_SMOOTH_FLAG 4
#define SSA_SHARP_FLAG 8

enum class INTEGRAL_TYPE
{
    INT_TRAP_FIXED,
    INT_TRAP_ADAPT,
    INT_SIMP_FIXED,
    INT_SIMP_ADAPT,
    INT_ROMB_ADAPT,
    INT_TRAP_NL,
    INT_HYBRID,
    INT_CADRE,
    INT_GK49_ADAPT,
    INT_GK715_ADAPT,
    INT_GK1021_ADAPT,
    INT_UNDEFINED
};

enum class GAMMA_TYPE
{
    GAMMA_INF,
    GAMMA_FLAT,
    GAMMA_EVENMASS,
    GAMMA_STRUCT
};

struct fluxParams
{
    double theta;
    double phi;
    double cp;
    double ct;
    double st;
    double cto;
    double sto;

    double theta_obs;
    double t_obs;
    double nu_obs;
    double d_L;

    double E_iso;
    double n_0;
    double g_init;

    double p;
    double epsilon_E;
    double epsilon_B;
    double ksi_N;

    double theta_h;
    double E_iso_core;
    double theta_core;
    double theta_wing;
    double b;
    double E_tot;
    double g_core;
    double E_core_global;
    double theta_core_global;

    int envType;
    double As;
    double Rwind;

    double L0;
    double q;
    double ts;

    double current_theta_cone_hi;
    double current_theta_cone_low;
    double theta_obs_cur;
    int tRes;
    int latRes;
    int spread;
    int counterjet;

    INTEGRAL_TYPE int_type;
    double rtol_struct;
    double rtol_theta;
    double rtol_phi;
    int nmax_theta;
    int nmax_phi;

    double atol_theta;

    double Rt0;
    double Rt1;
    double ta;
    double tb;

    double C_BMsqrd;
    double C_STsqrd;

    double t_NR;

    std::vector<double> t_table;
    std::vector<double> R_table;
    std::vector<double> u_table;
    std::vector<double> th_table;
    std::vector<double> mu_table;
    int table_entries;

    std::vector<double> t_table_inner;
    std::vector<double> R_table_inner;
    std::vector<double> u_table_inner;
    std::vector<double> th_table_inner;
    std::vector<double> mu_table_inner;
    int table_entries_inner;

    int spec_type;
    GAMMA_TYPE gamma_type;

    std::function<double(double, void *)> f_e;

    std::vector<double> mask;
    int nmask;

    long nevals;

    int error;
    char *error_msg;
};

constexpr void err_chk_void(fluxParams &pars)
{
    if (pars.error)
    {
        return;
    }
}

constexpr int err_chk_int(fluxParams &pars)
{
    if (pars.error)
    {
        return 0;
    }
}

constexpr double err_chk_dbl(fluxParams &pars)
{
    if (pars.error)
    {
        return 0.0;
    }
}

namespace afterglowpy
{
    /**
     * @brief Distribution function for tophat blast wave
     *
     * @param theta  angle
     * @param params struct of simulation parameters
     * @return double
     */
    double f_E_tophat(double theta, void *params);

    /**
     * @brief Distribution function for tophat blast wave
     *
     * @param theta  angle
     * @param params struct of simulation parameters
     * @return double
     */
    double f_E_Gaussian(double theta, void *params);

    /**
     * @brief Distribution function for gaussian blast wave
     *
     * @param theta  angle
     * @param params struct of simulation parameters
     * @return double
     */
    double f_E_powerlaw(double theta, void *params);

    /**
     * @brief Distribution function for two-component blast wave
     *
     * @param theta  angle
     * @param params struct of simulation parameters
     * @return double
     */
    double f_E_twocomponent(double theta, void *params);

    /**
     * @brief Distribution function for exponential blast wave
     *
     * @param theta  angle
     * @param params struct of simulation parameters
     * @return double
     */
    double f_E_exponential(double theta, void *params);

    /**
     * @brief Energy distribution function for tophat blast wave
     *
     * @param theta  angle
     * @param params struct of simulation parameters
     * @return double
     */
    double f_Etot_tophat(void *params);

    /**
     * @brief Energy distribution function for Gaussian blast wave
     *
     * @param theta  angle
     * @param params struct of simulation parameters
     * @return double
     */
    double f_Etot_Gaussian(void *params);

    

    /**
     * @brief Energy distribution function for power law blast wave
     *
     * @param theta  angle
     * @param params struct of simulation parameters
     * @return double
     */
    double f_Etot_powerlaw(void *params);

    /**
     * @brief Make table for radial zones
     *
     * @param pars
     */
    void make_R_table(fluxParams &pars);

    /**
     * @brief make table for cos(theta) zones
     *
     * @param pars
     */
    void make_mu_table(fluxParams &pars);

    /**
     * @brief Check that the time in the emitter frame is valid
     *
     * @param t_e
     * @param mu
     * @param t_obs
     * @param mu_table
     * @param N
     * @return double
     */
    double check_t_e(double t_e, double mu, double t_obs, std::vector<double> &mu_table, int N);

    /**
     * @brief
     *
     * @param x
     * @param arr
     * @param N
     * @return int
     */
    int searchSorted(double x, std::vector<double> &arr, int N);

    /**
     * @brief Perform linear interpolation
     *
     * @param a
     * @param b
     * @param x
     * @param X
     * @param Y
     * @param N
     * @return double
     */
    double interpolateLin(int a, int b, double x, std::vector<double> &X, std::vector<double> &Y, int N);

    /**
     * @brief Perform logarithmic interpolation
     *
     * @param a
     * @param b
     * @param x
     * @param X
     * @param Y
     * @param N
     * @return double
     */
    double interpolateLog(int a, int b, double x, std::vector<double> &X, std::vector<double> &Y, int N);

    /**
     * @brief Compute the jet edge extremum
     *
     * @param phi
     * @param cto
     * @param sto
     * @param theta0
     * @param a_mu
     * @param a_thj
     * @param N
     * @return double
     */
    double find_jet_edge(double phi, double cto, double sto, double theta0,
                        std::vector<double> &a_mu, std::vector<double> &a_thj, int N);

    /**
     * @brief Cosine theta integrand in spherical coordinate jacobican
     *
     * @param a_theta
     * @param params
     * @return double
     */
    double costheta_integrand(double a_theta, void *params); // inner integral

    /**
     * @brief Phi integrand in spherical coordinate jacobian
     *
     * @param a_phi
     * @param params
     * @return double
     */
    double phi_integrand(double a_phi, void *params); // outer integral

    /**
     * @brief compute zone emissivity
     *
     * @param nu
     * @param R
     * @param mu
     * @param te
     * @param u
     * @param us
     * @param n0
     * @param p
     * @param epse
     * @param epsB
     * @param ksiN
     * @param specType
     * @return double
     */
    double emissivity(
        double nu, 
        double R, 
        double mu, 
        double te,
        double u, 
        double us, 
        double n0, 
        double p, 
        double epse,
        double epsB, 
        double ksiN, 
        int specType);

    /**
     * @brief Compute zone flux given t_obs
     *
     * @param pars
     * @param atol
     * @return double
     */
    double flux(fluxParams &pars, double atol);

    /**
     * @brief Compute observed flux for conical structure
     *
     * @param t_obs
     * @param nu_obs
     * @param E_iso
     * @param theta_h
     * @param theta_cone_low
     * @param theta_cone_hi
     * @param atol
     * @param pars
     * @return double
     */
    double flux_cone(double t_obs, double nu_obs, double E_iso, double theta_h,
                    double theta_cone_low, double theta_cone_hi,
                    double atol, fluxParams &pars);

    /**
     * @brief Compute intensity in a zone
     *
     * @param theta
     * @param phi
     * @param tobs
     * @param nuobs
     * @param theta_obs
     * @param theta_cone_hi
     * @param theta_cone_low
     * @param pars
     * @return double
     */
    double intensity(
        double theta, 
        double phi, 
        double tobs, 
        double nuobs,
        double theta_obs, 
        double theta_cone_hi, 
        double theta_cone_low,
        fluxParams &pars);

    /**
     * @brief Compute primitive variables in the shock
     *
     * @param theta
     * @param phi
     * @param tobs
     * @param t
     * @param R
     * @param u
     * @param thj
     * @param theta_obs
     * @param theta_cone_hi
     * @param theta_cone_low
     * @param pars
     */
    void shockVals(
        double theta, 
        double phi, 
        double tobs,
        std::vector<double> &t, 
        std::vector<double> &R, 
        std::vector<double> &u, 
        std::vector<double> &thj,
        double theta_obs, 
        double theta_cone_hi, 
        double theta_cone_low,
        fluxParams &pars);

    /**
     * @brief Compute intensity in a conical blast wave
     *
     * @param theta
     * @param phi
     * @param t
     * @param nu
     * @param I
     * @param N
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param pars
     */
    void intensity_cone(
        std::vector<double> &theta, 
        std::vector<double> &phi, 
        std::vector<double> &t, 
        std::vector<double> &nu,
        std::vector<double> &I, 
        int N, 
        double E_iso_core,
        double theta_h_core, 
        double theta_h_wing,
        fluxParams &pars);

    /**
     * @brief Compute  complete intensity profile
     *
     * @param theta
     * @param phi
     * @param t
     * @param nu
     * @param I
     * @param N
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param res_cones
     * @param f_E
     * @param pars
     */
    void intensity_struct(
        std::vector<double> &theta, 
        std::vector<double> &phi, 
        std::vector<double> &t, 
        std::vector<double> &nu,
        std::vector<double> &I, 
        int N,
        double E_iso_core,
        double theta_h_core, 
        double theta_h_wing,
        int res_cones, 
        std::function<double(double, void *)> func,
        fluxParams &pars);

    /**
     * @brief Compute intensity structure in a core
     *
     * @param theta
     * @param phi
     * @param t
     * @param nu
     * @param I
     * @param N
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param res_cones
     * @param f_E
     * @param pars
     */
    void intensity_structCore(
        std::vector<double> &theta, 
        std::vector<double> &phi, 
        std::vector<double> &t, 
        std::vector<double> &nu,
        std::vector<double> &I, 
        int N,
        double E_iso_core,
        double theta_h_core, 
        double theta_h_wing,
        int res_cones, 
        std::function<double(double, void *)> func,
        fluxParams &pars);

    /**
     * @brief Compute primitive shock variables in a core structure
     *
     * @param theta
     * @param phi
     * @param tobs
     * @param t
     * @param R
     * @param u
     * @param thj
     * @param N
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param pars
     */
    void shockVals_cone(
        std::vector<double> &theta, 
        std::vector<double> &phi, 
        std::vector<double> &tobs,
        std::vector<double> &t, 
        std::vector<double> &R, 
        std::vector<double> &u, 
        std::vector<double> &thj, 
        int N,
        double E_iso_core,
        double theta_h_core, 
        double theta_h_wing,
        fluxParams &pars);

    /**
     * @brief Compute shock primitive values in non-core structure
     *
     * @param theta
     * @param phi
     * @param tobs
     * @param t
     * @param R
     * @param u
     * @param thj
     * @param N
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param res_cones
     * @param f_E
     * @param pars
     */
    void shockVals_struct(
        std::vector<double> &theta, 
        std::vector<double> &phi, 
        std::vector<double> &tobs,
        std::vector<double> &t, 
        std::vector<double> &R, 
        std::vector<double> &u, 
        std::vector<double> &thj, 
        int N,
        double E_iso_core,
        double theta_h_core, 
        double theta_h_wing,
        int res_cones, 
        std::function<double(double, void *)> func,
        fluxParams &pars);

    /**
     * @brief Compute shock primitives in core strucrure
     *
     * @param theta
     * @param phi
     * @param tobs
     * @param t
     * @param R
     * @param u
     * @param thj
     * @param N
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param res_cones
     * @param f_E
     * @param pars
     */
    void shockVals_structCore(
        std::vector<double> &theta, 
        std::vector<double> &phi, 
        std::vector<double> &tobs,
        std::vector<double> &t, 
        std::vector<double> &R, 
        std::vector<double> &u, 
        std::vector<double> &thj, 
        int N,
        double E_iso_core,
        double theta_h_core, 
        double theta_h_wing,
        int res_cones, 
        std::function<double(double, void *)> func,
        fluxParams &pars);

    /**
     * @brief Compute light curve for a tophat blast wave
     *
     * @param t
     * @param nu
     * @param F
     * @param Nt
     * @param E_iso
     * @param theta_h
     * @param pars
     */
    void lc_tophat(
        std::vector<double> &t, 
        std::vector<double> &nu, 
        std::vector<double> &F, 
        int Nt,
        double E_iso,
        double theta_h, 
        fluxParams &pars);

    /**
     * @brief Compute light curve for conical blast wave
     *
     * @param t
     * @param nu
     * @param F
     * @param Nt
     * @param E_iso
     * @param theta_h
     * @param theta_wing
     * @param pars
     */
    void lc_cone(
        std::vector<double> &t, 
        std::vector<double> &nu, 
        std::vector<double> &F, 
        int Nt, 
        double E_iso,
        double theta_core, 
        double theta_wing, 
        fluxParams &pars);

    /**
     * @brief Compute light curve for power law blast wave in core
     *
     * @param t
     * @param nu
     * @param F
     * @param Nt
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param beta
     * @param theta_c_arr
     * @param E_iso_arr
     * @param res_cones
     * @param pars
     */
    void lc_powerlawCore(
        std::vector<double> &t, 
        std::vector<double> &nu, 
        std::vector<double> &F, 
        int Nt,
        double E_iso_core, 
        double theta_h_core,
        double theta_h_wing, 
        double beta,
        std::vector<double> &theta_c_arr, 
        std::vector<double> &E_iso_arr,
        int res_cones, fluxParams &pars);

    /**
     * @brief Compute light curve for power law blast wave
     *
     * @param t
     * @param nu
     * @param F
     * @param Nt
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param theta_c_arr
     * @param E_iso_arr
     * @param res_cones
     * @param pars
     */
    void lc_powerlaw(
        std::vector<double> &t, 
        std::vector<double> &nu, 
        std::vector<double> &F, 
        int Nt,
        double E_iso_core, 
        double theta_h_core,
        double theta_h_wing,
        std::vector<double> &theta_c_arr, 
        std::vector<double> &E_iso_arr,
        int res_cones, fluxParams &pars);

    /**
     * @brief Compyte light curve in Gaussian structure
     *
     * @param t
     * @param nu
     * @param F
     * @param Nt
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param theta_c_arr
     * @param E_iso_arr
     * @param res_cones
     * @param pars
     */
    void lc_Gaussian(
        std::vector<double> &t, 
        std::vector<double> &nu, 
        std::vector<double> &F, 
        int Nt,
        double E_iso_core,
        double theta_h_core, 
        double theta_h_wing,
        std::vector<double> &theta_c_arr, 
        std::vector<double> &E_iso_arr,
        int res_cones, fluxParams &pars);

    /**
     * @brief Compute light curve in Gaussian core structure
     *
     * @param t
     * @param nu
     * @param F
     * @param Nt
     * @param E_iso_core
     * @param theta_h_core
     * @param theta_h_wing
     * @param theta_c_arr
     * @param E_iso_arr
     * @param res_cones
     * @param pars
     */
    void lc_GaussianCore(std::vector<double> &t, std::vector<double> &nu, std::vector<double> &F, int Nt,
                        double E_iso_core,
                        double theta_h_core, double theta_h_wing,
                        std::vector<double> &theta_c_arr, std::vector<double> &E_iso_arr,
                        int res_cones, fluxParams pars);

    /**
     * @brief Calculate the total flux density
     *
     * @param jet_type
     * @param spec_type
     * @param t
     * @param nu
     * @param Fnu
     * @param N
     * @param fp
     */
    void calc_flux_density(
        int jet_type, 
        int spec_type, 
        std::vector<double> &t,
        std::vector<double> &nu, 
        std::vector<double> &Fnu, 
        int N,
        fluxParams &fp);

    /**
     * @brief Calculate the totqal intensity in a zone
     *
     * @param jet_type
     * @param spec_type
     * @param theta
     * @param phi
     * @param t
     * @param nu
     * @param Inu
     * @param N
     * @param fp
     */
    void calc_intensity(
        int jet_type, 
        int spec_type, 
        std::vector<double> &theta, 
        std::vector<double> &phi,
        std::vector<double> &t, 
        std::vector<double> &nu, 
        std::vector<double> &Inu, 
        int N,
        fluxParams &fp);

    /**
     * @brief Calculate the shock primitive values
     *
     * @param jet_type
     * @param theta
     * @param phi
     * @param tobs
     * @param t
     * @param R
     * @param u
     * @param thj
     * @param N
     * @param fp
     */
    void calc_shockVals(int jet_type, std::vector<double> &theta, std::vector<double> &phi, std::vector<double> &tobs,
                        std::vector<double> &t, std::vector<double> &R, std::vector<double> &u, std::vector<double> &thj, int N,
                        fluxParams &p);

    /**
     * @brief Setup the flux parameters for computation
     *
     * @param pars
     * @param d_L
     * @param theta_obs
     * @param E_iso_core
     * @param theta_core
     * @param theta_wing
     * @param b
     * @param L0
     * @param q
     * @param ts
     * @param n_0
     * @param p
     * @param epsilon_E
     * @param epsilon_B
     * @param ksi_N
     * @param g0
     * @param E_core_global
     * @param theta_core_global
     * @param ta
     * @param tb
     * @param tRes
     * @param latRes
     * @param int_type
     * @param rtol_struct
     * @param rtol_phi
     * @param rtol_theta
     * @param nmax_phi
     * @param nmax_theta
     * @param spec_type
     * @param mask
     * @param nmask
     * @param spread
     * @param counterjet
     * @param gamma_type
     */
    void setup_fluxParams(
        fluxParams &pars,
        double d_L,
        double theta_obs,
        double E_iso_core, 
        double theta_core, 
        double theta_wing,
        double b, 
        double L0, 
        double q, 
        double ts,
        double n_0,
        double p,
        double epsilon_E,
        double epsilon_B,
        double ksi_N,
        double g0,
        double E_core_global,
        double theta_core_global,
        double ta, 
        double tb,
        int tRes, 
        int latRes, 
        INTEGRAL_TYPE int_type,
        double rtol_struct,
        double rtol_phi, 
        double rtol_theta,
        int nmax_phi, 
        int nmax_theta,
        int spec_type,
        std::vector<double> &mask, 
        int nmask,
        int spread, int counterjet, 
        GAMMA_TYPE gamma_type);

    /**
     * @brief Set the jet params objects
     *
     * @param pars
     * @param E_iso
     * @param theta_h
     */
    void set_jet_params(fluxParams &pars, double E_iso, double theta_h);

    /**
     * @brief Set the obs params object
     *
     * @param pars
     * @param t_obs
     * @param nu_obs
     * @param theta_obs_cur
     * @param current_theta_cone_hi
     * @param current_theta_cone_low
     */
    void set_obs_params(fluxParams &pars, double t_obs, double nu_obs,
                        double theta_obs_cur, double current_theta_cone_hi,
                        double current_theta_cone_low);
    /**
     * @brief Check for any errors
     *
     * @param params
     * @return int
     */
    int check_error(void *params);

    /**
     * @brief Set the error object
     *
     * @param pars
     * @param msg
     */
    void set_error(fluxParams &pars, char msg[]);

    void free_fluxParams(fluxParams &pars);
    
} // namespace afterglowpy

#endif
