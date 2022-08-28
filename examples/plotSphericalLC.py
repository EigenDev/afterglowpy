import afterglowpy as grb 
import time
import numpy as np 
import matplotlib.pyplot as plt 

params = {
    "jetType": grb.Jet.Spherical,
    "specType": 0,
    "thetaObs": 0.05,
    "E0": 1e53,
    "thetaCore": 0.1,
    "n0": 1,
    "epsilon_e": 0.1,
    "epsilon_b": 0.1,
    "xi_N": 1.0,
    "theta_wing": 0.1,
    "thetaWing": 0.1,
    "z": 0.55,
    "uMax": 1e3,
    "uMin": 10,
    "Er": 1e51,
    "MFast_solar": 1e-6,
    "k": 0,
    "Einj": 0,
}

t    = np.geomspace(1.0e3, 1.0e7, 300)
nu   = np.ones_like(t) * 1e9
tday = t * grb.sec2day
t1 = (time.time())
a    = grb.fluxDensity(t, nu, **params)
print(f"done in: {time.time() - t1:.2e} s")
zzz = input('')

plt.loglog(tday, a)
plt.show()