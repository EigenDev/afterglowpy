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
    "epsilon_E": 0.1,
    "epsilon_b": 0.1,
    "xi_N": 1.0,
    "thetaWing": 0.1,
    "z": 0.55,
    "uMax": 1e3,
    "uMin": 10,
    "Er": 1e51,
    "MFast_solar": 1e-6,
    "k": 0,
    "Einj": 0,
}
plt.rc('text', usetex=True)
nu0 = 1e9
t    = np.geomspace(1.0e3, 1.0e7, 300)
nu   = np.ones_like(t) * nu0
tday = t * grb.sec2day
print("Computing....")
a    = grb.fluxDensity(t, nu, **params)

print("Plotting...")
oom = int(np.floor(np.log10(nu0)))
fig, ax = plt.subplots(1, 1)
ax.loglog(tday, a)
ax.set_xlabel(r'$t [\rm{day}$')
ax.set_ylabel(r'$F_\nu[10^{%d} \rm Hz]$'%(oom))
plt.show()