
#! /usr/bin/env python
import numpy as np
import afterglowpy as grb
import matplotlib.pyplot as plt 

t0  = 1.0    # initial time in s
t1  = 1.0e8  # final time in s
g0  = 1.0e3  # initial Lorentz factor
E0  = 1.0e53 # Eiso in erg
n0  = 1.0    # ISM density in cm^{-3}
th0 = 0.1   # initial width of this jet sector in rad

b0 = np.sqrt(1.0 - 1.0 / (g0*g0))
u0 = b0 * g0
R0 = grb.c * b0 * t0

t = np.geomspace(t0, t1, 8000)  # array of lab times in s.
# the solver will take a single RK4 step between each of these times.
M0   = E0 / ((g0 - 1) * grb.c*grb.c)
rho0 = n0 * grb.mp

thC = th0        # angular size of jet. (controls when spreading starts)
spreadType = 7   # 0: no spreading, 7: spreading for tophat of width th0, 8: spreading for cone of with th0 in a jet with size thC

R, u, thj = grb.shock.shockEvolSpreadRK4(t, R0, u0, th0,
                                         M0, rho0, 
                                         0, 0, 0, 0, 0, 0,
                                         thC, spreadType)

plt.semilogx(t, u)
plt.show()