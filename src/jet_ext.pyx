# disutils: language = c++
from libcpp.vector cimport vector 
cimport numpy as np 
import numpy as np 

cdef extern from "../include/offaxis_struct.hpp":
    cdef cppclass INTEGRAL_TYPE:
        pass 

    cdef cppclass GAMMA_TYPE:
        pass

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

        INTEGRAL_TYPE int_type
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
        GAMMA_TYPE gamma_type

        void *f_e

        vector[double] mask
        int nmask

        long nevals

        int error
        char *error_msg
cdef extern from "../include/offaxis_struct.hpp" namespace "INTEGRAL_TYPE":
    cdef INTEGRAL_TYPE INT_TRAP_FIXED
    cdef INTEGRAL_TYPE INT_TRAP_ADAPT
    cdef INTEGRAL_TYPE INT_SIMP_FIXED
    cdef INTEGRAL_TYPE INT_SIMP_ADAPT
    cdef INTEGRAL_TYPE INT_ROMB_ADAPT
    cdef INTEGRAL_TYPE INT_TRAP_NL
    cdef INTEGRAL_TYPE INT_HYBRID
    cdef INTEGRAL_TYPE INT_CADRE
    cdef INTEGRAL_TYPE INT_GK49_ADAPT
    cdef INTEGRAL_TYPE INT_GK715_ADAPT
    cdef INTEGRAL_TYPE INT_GK1021_ADAPT
    cdef INTEGRAL_TYPE INT_UNDEFINED

cdef extern from "../include/offaxis_struct.hpp" namespace "GAMMA_TYPE":
    cdef GAMMA_TYPE GAMMA_INF
    cdef GAMMA_TYPE GAMMA_FLAT
    cdef GAMMA_TYPE GAMMA_EVENMASS
    cdef GAMMA_TYPE GAMMA_STRUCT

cdef extern from "../include/offaxis_struct.hpp" namespace "afterglowpy":
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
        vector[double] mask, 
        int nmask,
        int spread, 
        int counterjet, 
        GAMMA_TYPE gamma_type
    )

    void calc_flux_density(
        int jet_type, 
        int spec_type,
        vector[double] t, 
        vector[double] nu, 
        vector[double] Fnu, 
        int N,
        fluxParams fp
    )

cdef GAMMA_TYPE get_gamma_type(val: str):
    if val == 'gamma_inf':
        return GAMMA_INF 
    elif val == 'gamma_flat':
        return GAMMA_FLAT
    elif val == 'gamma_evenmass':
        return GAMMA_EVENMASS
    elif val == 'gamma_struct':
        return GAMMA_STRUCT

cdef INTEGRAL_TYPE get_intergal_type(val: str):
    if val == 'int_trap_fixed':
        return INT_TRAP_FIXED
    elif val == 'int_trap_adapt':
        return INT_TRAP_ADAPT
    elif val == 'int_simp_fixed':
        return INT_SIMP_FIXED
    elif val == 'int_simp_adapt':
        return INT_SIMP_ADAPT
    elif val == 'int_romb_adapt':
        return INT_ROMB_ADAPT
    elif val == 'int_trap_nl':
        return INT_TRAP_NL
    elif val == 'int_hybid':
        return INT_HYBRID
    elif val == 'int_cadre':
        return INT_CADRE
    elif val == 'int_gk49_adapt':
        return INT_GK49_ADAPT
    elif val == 'int_gk715_adapt':
        return INT_GK715_ADAPT
    elif val == 'int_gk1021_adapt':
        return INT_GK1021_ADAPT
    else:
        return INT_UNDEFINED

def calc_flux_density_wrapper(
    jet_type:  int, 
    spec_type: int,
    t:         vector[double],
    nu:        vector[double],
    fnu:       vector[double],
    n:         int,
    fp:        dict):
    cdef fluxParams flux_params
    cdef INTEGRAL_TYPE intergraL_type
    cdef GAMMA_TYPE gamma_type

    intergraL_type = get_intergal_type(fp['int_type'])
    gamma_type     = get_gamma_type(fp['gamma_type'])
    setup_fluxParams(
        flux_params,
        fp['d_L'],
        fp['theta_obs'],
        fp['E_iso_core'], 
        fp['thetaCore'], 
        fp['theta_wing'],
        fp['b'], 
        fp['L0'], 
        fp['q'], 
        fp['ts'],
        fp['n_0'],
        fp['p'],
        fp['epsilon_E'],
        fp['epsilon_B'],
        fp['ksi_N'],
        fp['g0'],
        fp['E_core_global'],
        fp['theta_core_global'],
        fp['ta'], 
        fp['tb'],
        fp['tRes'], 
        fp['latRes'], 
        intergraL_type,
        fp['rtol_struct'],
        fp['rtol_phi'], 
        fp['rtol_theta'],
        fp['nmax_phi'], 
        fp['nmax_theta'],
        fp['specType'],
        fp['mask'], 
        fp['nmask'],
        fp['spread'], 
        fp['counterjet'], 
        gamma_type
    )

    calc_flux_density(
        jet_type, 
        spec_type,
        t, 
        nu, 
        fnu, 
        n,
        flux_params
    )

    return fnu