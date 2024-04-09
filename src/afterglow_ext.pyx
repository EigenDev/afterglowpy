# disutils: language = c++
from libcpp.vector cimport vector 
cimport numpy as np 
import numpy as np 
import scipy.integrate as integrate
import math 

c      = 2.99792458e10
me     = 9.1093897e-28
mp     = 1.6726231e-24
h      = 6.6260755e-27
hbar   = 1.05457266e-27
ee     = 4.803e-10
sigmaT = 6.65e-25

Msun    = 1.98892e33
cgs2mJy = 1.0e26
mJy2cgs = 1.0e-26
deg2rad = np.pi/180.0
rad2deg = 180.0/np.pi
day2sec = 86400.0
sec2day = 1.0/day2sec
parsec  = 3.0857e18
Hz2eV   = 4.13566553853599e-15
eV2Hz   = 1.0/Hz2eV
 
cdef extern from "../include/offaxis_struct.hpp":
    cdef cppclass INTEGRAL_TYPE:
        pass 

    cdef cppclass GAMMA_TYPE:
        pass

    cdef cppclass fluxParams:
        pass 
        
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
        vector[double] &mask, 
        int nmask,
        int spread, 
        int counterjet, 
        GAMMA_TYPE gamma_type,
        int jetType
    )

    void calc_flux_density(
        int jet_type, 
        int spec_type,
        vector[double] &t, 
        vector[double] &nu, 
        vector[double] &Fnu, 
        int N,
        fluxParams &fp
    )

    void calc_intensity(
        int jet_type, 
        int spec_type, 
        vector[double] &theta, 
        vector[double] &phi,
        vector[double] &t, 
        vector[double] &nu, 
        vector[double] &Inu, 
        int N, 
        fluxParams &fp
    )

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
        int    specType
    )

cdef extern from "../include/shockEvolution.hpp" namespace "afterglowpy":
    void shockEvolveSpreadRK4(
        vector[double] &t, 
        vector[double] &R, 
        vector[double] &u, 
        vector[double] &th, 
        int N, 
        double R0, 
        double u0, 
        double th0, 
        void *args,
        int   spread
    )

    void shockEvolveRK4(
        vector[double] &t,
        vector[double] &R, 
        vector[double] &u, 
        int N, 
        double R0, 
        double u0, 
        void *args
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

def dP(
    double costheta, 
    np.ndarray[np.float64_t, ndim=1] amu, 
    np.ndarray[np.float64_t, ndim=1] ate, 
    np.ndarray[np.float64_t, ndim=1] au, 
    np.ndarray[np.float64_t, ndim=1] ar, 
    double nu, 
    double n0, 
    double p, 
    double epsE, 
    double epsB, 
    double ksiN, 
    specType: int) -> double:

    mu = costheta
    ib = np.searchsorted(amu, mu)
    N  = amu.shape[0]
    if ib <= 0:
        ib = 1
    elif ib >= N:
        ib = N-1
    ia = ib-1 

    te = ((mu-amu[ia])*ate[ib] + (amu[ib]-mu)*ate[ia]) / (amu[ib]-amu[ia])
    u  = au[ia]*math.pow(te/ate[ia], math.log(au[ib]/au[ia])
                        / math.log(ate[ib]/ate[ia]))
    r = ar[ia]*math.pow(te/ate[ia], math.log(ar[ib]/ar[ia])
                        / math.log(ate[ib]/ate[ia]))

    g = math.sqrt(u*u+1)

    us = 4*u*g / math.sqrt(8*u*u+9)

    em = emissivity(
        nu, r, mu, te, u, us, n0, p, epsE,
        epsB, ksiN, specType
    )

    return 2*math.pi * em
    
def spherical_flux_density_wrapper(t, nu, **kwargs):
    t           = np.asanyarray(t)
    specType    = kwargs['specType']
    uMax        = kwargs['uMax']   
    uMin        = kwargs['uMin']
    Er          = kwargs['Er']
    k           = kwargs['k']
    MFast_solar = kwargs['MFast_solar']
    n0          = kwargs['n0']
    p           = kwargs['p']
    epsilon_E   = kwargs['epsilon_E']
    epsilon_B   = kwargs['epsilon_B']
    ksiN        = kwargs['xi_N']
    dL          = kwargs['d_L']
    Einj        = kwargs['Einj']

    # Energy injection variables (off by default)
    L0  = kwargs['L0'] if 'L0' in kwargs else 0.0
    q   = kwargs['q']  if 'q'  in kwargs else 0.0
    ts  = kwargs['ts'] if 'ts' in kwargs else 0.0
   
    # Numerical integration variables
    rtol   = kwargs['rtol']   if 'rtol'   in kwargs else 1.0e-3
    tRes   = kwargs['tRes']   if 'tRes'   in kwargs else 1000
    latRes = kwargs['latRes'] if 'latRes' in kwargs else 0

    rho0 = mp * n0
    Mej  = MFast_solar * Msun
    u0   = uMax
    g0   = np.sqrt(1+u0*u0)
    bes0 = 4*u0*g0 / (4*u0*u0+3)
    Rd   = math.pow(9*g0*g0*Mej / (4*np.pi*(g0+1)*(4*u0*u0+3)*rho0), 1./3.)
    td   = Rd / (bes0 * c)

    t0 = min(1.0e-2*td, 5.0e-1 * g0*g0*t.min(), 5.0e-1 * t.min()/(1+bes0))
    t1 = 2. * g0*g0*t.max()

    NT = int(tRes * np.log10(t1/t0))
    # print("{0:.3e} {1:.3e} {2:.3e} {3:d}".format(t0, t1, t1/t0, NT))

    r0   = bes0*c*t0
    # Vej0 = 4.0/3.0*np.pi*r0*r0*r0

    cdef ate = np.geomspace(t0, t1, num=NT)
    cdef vector[double] arc
    cdef vector[double] auc 
    for i in range(ate.size):
        arc.push_back(0) 
        auc.push_back(0)

    cdef double shockArgs[9]
    shockArgs[:] = [u0, Mej, rho0, Einj, k, uMin, L0, q, ts]
    shockEvolveRK4(
        ate,
        arc, 
        auc, 
        NT, 
        r0, 
        u0, 
        shockArgs
    )
    ar = np.array(arc)
    au = np.array(auc)
    P = np.zeros_like(t)
    wopts = None

    for i in range(len(t)):
        amu = c * (ate - t[i]) / ar

        args = (amu, ate, au, ar, nu[i], n0, p, epsilon_E, epsilon_B, ksiN,
                specType)

        res = integrate.quad(dP, 0.0, 1.0, args, full_output=1, wopts=wopts,
                             epsrel=rtol)
        P[i] = res[0]

    Fnu = cgs2mJy * P / (4*np.pi*dL*dL)

    return Fnu

def jet_flux_density_wrapper(
    int jet_type, 
    int spec_type,
    vector[double] t,
    vector[double] nu,
    vector[double] fnu,
    int n,
    dict fp):
    cdef fluxParams flux_params
    cdef INTEGRAL_TYPE intergraL_type = get_intergal_type(fp['int_type'])
    cdef GAMMA_TYPE gamma_type = get_gamma_type(fp['gamma_type'])
    setup_fluxParams(
        flux_params,
        fp['d_L'],
        fp['thetaObs'],
        fp['E0'], 
        fp['thetaCore'], 
        fp['thetaWing'],
        fp['b'], 
        fp['L0'], 
        fp['q'], 
        fp['ts'],
        fp['n0'],
        fp['p'],
        fp['epsilon_E'],
        fp['epsilon_B'],
        fp['xi_N'],
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
        gamma_type,
        jet_type
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

def intensity_wrapper(
    int jet_type, 
    int spec_type,
    vector[double] t,
    vector[double] nu,
    vector[double] theta,
    vector[double] phi,
    vector[double] Inu,
    int n,
    dict fp): 

    cdef fluxParams flux_params
    cdef INTEGRAL_TYPE intergraL_type = get_intergal_type(fp['int_type'])
    cdef GAMMA_TYPE gamma_type =  get_gamma_type(fp['gamma_type'])

    setup_fluxParams(
        flux_params,
        fp['d_L'],
        fp['thetaObs'],
        fp['E0'], 
        fp['thetaCore'], 
        fp['thetaWing'],
        fp['b'], 
        fp['L0'], 
        fp['q'], 
        fp['ts'],
        fp['n0'],
        fp['p'],
        fp['epsilon_E'],
        fp['epsilon_B'],
        fp['xi_N'],
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
        gamma_type,
        jet_type
    )

    calc_intensity(
        jet_type, 
        spec_type, 
        theta, 
        phi,
        t, 
        nu, 
        Inu, 
        n, 
        flux_params
    )

    return Inu