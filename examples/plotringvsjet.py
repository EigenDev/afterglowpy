import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb
from matplotlib.lines import Line2D

def order_of_mag(var: float) -> int:
    return int(np.floor(np.log10(var)))

def front_part(var: float) -> float:
    oom = order_of_mag(var)
    return var / 10**oom 
plt.rc('text', usetex=True)
eiso = 5e52
if front_part(eiso) == 1:
    eiso_label = r'$E_{{\rm iso}} = 10^{%d} \rm erg$'%(order_of_mag(eiso))
else:
    eiso_label = r'$E_{{\rm iso}} = %1.f \times 10^{%d} \rm erg$'%(front_part(eiso), order_of_mag(eiso))
# For convenience, place arguments into a dict.
Zjet = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
    'specType':    0,                  # Basic Synchrotron Emission Spectrum
    'spread':      False,
    'thetaObs':    0.00,   # Viewing angle in radians
    'E0':          eiso, # Isotropic-equivalent energy in erg
    'thetaCore':   0.1,    # Half-opening angle in radians
    'n0':          1.0,    # circumburst density in cm^{-3}
    'p':           2.5,    # electron energy distribution index
    'epsilon_e':   0.1,    # epsilon_e
    'epsilon_B':   0.1,   # epsilon_B
    'xi_N':        1.0,    # Fraction of electrons accelerated
    'd_L':         1.0e28, # Luminosity distance in cm
    'z':           0.0}   # redshift

Zring = {
    'jetType':     grb.jet.Cone,     # Top-Hat jet
    'specType':    0,                  # Basic Synchrotron Emission Spectrum
    'spread':      False,
    'thetaObs':    np.pi / 2.0,   # Viewing angle in radians
    'E0':          eiso, # Isotropic-equivalent energy in erg
    'thetaCore':   np.pi/2 - 0.1,    # Half-opening angle in radians
    'thetaWing':   np.pi / 2,
    'counterjet':  True,
    'n0':          1.0,    # circumburst density in cm^{-3}
    'p':           2.5,    # electron energy distribution index
    'epsilon_e':   0.1,    # epsilon_er'$E_{\rm iso} = 10^{53} \rm erg$'
    'epsilon_B':   0.1,   # epsilon_B
    'xi_N':        1.0,    # Fraction of electrons accelerated
    'd_L':         1.0e28, # Luminosity distance in cm
    'z':           0.0}   # redshift

Zring2 = {
    'jetType':     grb.jet.Ring,     # Top-Hat jet
    'specType':    0,                  # Basic Synchrotron Emission Spectrum
    'spread':      False,
    'thetaObs':    np.pi / 2.0,   # Viewing angle in radians
    'E0':          eiso, # Isotropic-equivalent energy in erg
    'thetaCore':   np.pi/2 - 0.1,    # Half-opening angle in radians
    'thetaWing':   np.pi / 2,
    'counterjet':  True,
    'n0':          1.0,    # circumburst density in cm^{-3}
    'p':           2.5,    # electron energy distribution index
    'epsilon_e':   0.1,    # epsilon_er'$E_{\rm iso} = 10^{53} \rm erg$'
    'epsilon_B':   0.1,   # epsilon_B
    'xi_N':        1.0,    # Fraction of electrons accelerated
    'd_L':         1.0e28, # Luminosity distance in cm
    'z':           0.0}   # redshift

ts = 1e-1 * grb.day2sec
te = 1e5 * grb.day2sec
nus = [1e9, 1e12, 1e15, 1e17]
nt = 300 
nf = len(nus)

# Space time points geometrically, from 10^3 s to 10^7 s
t = np.zeros(shape=(nt, nf))

# Calculate flux in a single X-ray band (all times have same frequency)
nu = np.empty(t.shape)
for idx, freq in enumerate(nus):
    nu[:, idx] = freq
    
t[:, :] =  np.geomspace(ts, te, num=nt)[:, None]
tday = t * grb.sec2day
# Calculate!

Fnu_jet   = grb.fluxDensity(t, nu, **Zjet)
Fnu_ring  = grb.fluxDensity(t, nu, **Zring)
# Fnu_ring2 = grb.fluxDensity(t, nu, **Zring2)
# Plot!

print("Plotting")
fig, ax     = plt.subplots(1, 1, figsize=(5,4))
vmin, vmax  = 0, 1
color_range = np.linspace(0, 1, len(nus))
color_scale = plt.get_cmap('viridis')
norm        = plt.Normalize(vmin, vmax)
colors      = color_scale(norm(color_range))
# plt.rc('font', family='serif')
reg_lines = [0] * len(nus)
for idx in range(len(nus)):
    order_of_mag = int(np.floor(np.log10(nus[idx])))
    reg_lines[idx], = ax.plot(tday[:, idx], Fnu_jet[:, idx], color=colors[idx], label = r'$10^{%s} \rm Hz$'%(order_of_mag))
    ax.plot(tday[:, idx], Fnu_ring[:, idx], linestyle='--', color=colors[idx])
    # ax.plot(tday[:, idx], Fnu_ring2[:, idx], linestyle=':', color=colors[idx])
ax.set(xscale='log', xlabel=r'$t [\rm s]$',
       yscale='log', ylabel=r'$F_\nu [\rm mJy]$')

ring_line  = Line2D([0],[0],linestyle='--', color='gray', label='conical ring')
ring_line2 = Line2D([0],[0],linestyle=':', color='gray', label='structured ring')
ax.set_title(eiso_label)
ax.legend(handles=[*reg_lines, ring_line])
# ax.legend()
fig.tight_layout()
print("Saving figure lc.png")
fig.savefig("lc.pdf")
plt.show()
plt.close(fig)
