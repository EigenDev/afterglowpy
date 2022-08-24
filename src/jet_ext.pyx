# disutils: language = c++
from libcpp.vector cimport vector 
cimport numpy as np 
import numpy as np 

cdef extern from "../include/offaxis_struct.hpp":
    # cpdef enum class INTEGRAL_TYPE(int):
    #     NT_TRAP_FIXED,
    #     INT_TRAP_ADAPT,
    #     INT_SIMP_FIXED,
    #     INT_SIMP_ADAPT,
    #     INT_ROMB_ADAPT,
    #     INT_TRAP_NL,
    #     INT_HYBRID,
    #     INT_CADRE,
    #     INT_GK49_ADAPT,
    #     INT_GK715_ADAPT,
    #     INT_GK1021_ADAPT,
    #     INT_UNDEFINED,

    # cpdef enum class GAMMA_TYPE(int):
    #     GAMMA_INF,
    #     GAMMA_FLAT,
    #     GAMMA_EVENMASS,
    #     GAMMA_STRUCT,

    cdef cppclass fluxParams:
        double theta
        double phi
        double cp
        double ct
        double st
        double cto
        double sto

        double theta_obs
        double t_obs
        double nu_obs
        double d_L

        double E_iso
        double n_0
        double g_init

        double p
        double epsilon_E
        double epsilon_B
        double ksi_N

        double theta_h
        double E_iso_core
        double theta_core
        double theta_wing
        double b
        double E_tot
        double g_core
        double E_core_global
        double theta_core_global

        int envType
        double As
        double Rwind

        double L0
        double q
        double ts

        double current_theta_cone_hi
        double current_theta_cone_low
        double theta_obs_cur
        int tRes
        int latRes
        int spread
        int counterjet

        int int_type
        double rtol_struct
        double rtol_theta
        double rtol_phi
        int nmax_theta
        int nmax_phi

        double atol_theta

        double Rt0
        double Rt1
        double ta
        double tb

        double C_BMsqrd
        double C_STsqrd

        double t_NR

        vector[double] t_table
        vector[double] R_table
        vector[double] u_table
        vector[double] th_table
        vector[double] mu_table
        int table_entries

        vector[double] t_table_inner
        vector[double] R_table_inner
        vector[double] u_table_inner
        vector[double] th_table_inner
        vector[double] mu_table_inner
        int table_entries_inner

        int spec_type
        int gamma_type

        void *f_e

        vector[double] mask
        int nmask

        long nevals

        int error
        char *error_msg

    double calc_flux_density(
        int jet_type, 
        int spec_type,
        vector[double] t, 
        vector[double] nu, 
        vector[double] Fnu, 
        int N,
        fluxParams fp
    )


def calc_flux_density_wrapper(
    jet_type: int, 
    spec_type: int,
    t: np.ndarray,
    nu: np.ndarray,
    Fnu: np.ndarray,
    n: int,
    fp):
    pass