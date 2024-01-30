import math 
import numpy as np 
from enum import IntEnum
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
deg2rad = math.pi/180.0
rad2deg = 180.0/math.pi
day2sec = 86400.0
sec2day = 1.0/day2sec
parsec  = 3.0857e18
Hz2eV   = 4.13566553853599e-15
eV2Hz   = 1.0/Hz2eV


class Jet(IntEnum):
    Cone           = -2,
    TopHat         = -1,
    Gaussian       = 0,
    PowerLawCore   = 1,
    GaussianCore   = 2,
    Spherical      = 3,
    PowerLaw       = 4,
    Exponential    = 5,
    Twocomponent   = 6,
    Exponential2   = 7,
    Ring           = 8,
    GaussianRing   = 9,
    
default_params = {
    "E_iso_core":  1.0e53,
    "theta_h_core":  0.1,
    "theta_h_wing":  0.4,
    "b": 0.0,
    "L0": 0.0,
    "q": 0.0,
    "ts": 0.0, 
    "p": 2.2,
    "d_L": 1.0e28,

    "latRes":  5,
    "tRes":  1000,
    "spread":  7,
    "counterjet":  0,
    "gamma_type":  "gamma_inf",
    "g0": -1.0,
    "E_core_global":  0.0,
    "theta_core_global": 0.0,

    "rtol_struct": 1.0e-2,
    "rtol_theta" : 1.0e-2,
    "rtol_phi": 1.0e-2,
    "int_type": "int_cadre",
    "nmax_phi": 1000,
    "nmax_theta": 1000,
    
    "jetType": Jet.Spherical,
    "specType": 0,
    "thetaObs": 0.05,
    "E0": 1e53,
    "thetaCore": 0.1,
    "n0": 1,
    "epsilon_E": 0.1,
    "epsilon_B": 0.1,
    "xi_N": 1.0,
    "thetaWing": 0.1,
    "mask": np.empty(1),
    "nmask": 0,   
    "z": 0.55,
}