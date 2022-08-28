#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb
import argparse
from matplotlib.lines import Line2D

def main():
    parser = argparse.ArgumentParser("arg parser for ring vs jet afterglow calculation")
    parser.add_argument('--eiso',  dest='eiso',  help='isotropic equivalent energy in erg', type=float, default=1e53)
    parser.add_argument('--sring', dest='sring', help='set if want structured ring', default=False,action=argparse.BooleanOptionalAction)
    parser.add_argument('--tdomain', dest='tdomain', help='time domain for afterglow in days', default=[1e-2,1e3], type=float, nargs='+')
    parser.add_argument('--fname', help='name of plot file to be saved', default='some_lc', type=str, dest='file_name')
    args = parser.parse_args()
    
    def order_of_mag(var: float) -> int:
        return int(np.floor(np.log10(var)))

    def front_part(var: float) -> float:
        oom = order_of_mag(var)
        return var / 10**oom 
    
    plt.rc('text', usetex=True)
    
    eiso = args.eiso
    
    if front_part(eiso) == 1:
        eiso_label = r'$E_{{\rm iso}} = 10^{%d} \rm erg$'%(order_of_mag(eiso))
    else:
        eiso_label = r'$E_{{\rm iso}} = %1.f \times 10^{%d} \rm erg$'%(front_part(eiso), order_of_mag(eiso))
    # For convenience, place arguments into a dict.
    Zjet = {'jetType':     grb.Jet.TopHat,     # Top-Hat jet
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
        'jetType':     grb.Jet.Cone,     # Top-Hat jet
        'specType':    0,                  # Basic Synchrotron Emission Spectrum
        'spread':      False,
        'thetaObs':    np.pi / 2.0,   # Viewing angle in radians
        'E0':          eiso, # Isotropic-equivalent energy in erg
        'thetaCore':   np.pi/2 - 0.17/2,    # Half-opening angle in radians
        'thetaWing':   np.pi / 2,
        'counterjet':  True,
        'n0':          1.0,    # circumburst density in cm^{-3}
        'p':           2.5,    # electron energy distribution index
        'epsilon_e':   0.1,    # epsilon_er'$E_{\rm iso} = 10^{53} \rm erg$'
        'epsilon_B':   0.1,   # epsilon_B
        'xi_N':        1.0,    # Fraction of electrons accelerated
        'd_L':         1.0e28, # Luminosity distance in cm
        'z':           0.0}   # redshift

    if args.sring:
        Zring2 = {
            'jetType':     grb.Jet.Ring,     # Top-Hat jet
            'specType':    0,                  # Basic Synchrotron Emission Spectrum
            'spread':      False,
            'thetaObs':    np.pi / 2.0,   # Viewing angle in radians
            'E0':          eiso, # Isotropic-equivalent energy in erg
            'thetaCore':   np.pi/2 - 0.17/2,    # Half-opening angle in radians
            'thetaWing':   np.pi / 2,
            'counterjet':  True,
            'n0':          1.0,    # circumburst density in cm^{-3}
            'p':           2.5,    # electron energy distribution index
            'epsilon_e':   0.1,    # epsilon_er'$E_{\rm iso} = 10^{53} \rm erg$'
            'epsilon_B':   0.1,   # epsilon_B
            'xi_N':        1.0,    # Fraction of electrons accelerated
            'd_L':         1.0e28, # Luminosity distance in cm
            'z':           0.0}   # redshift
    tdomain  = np.asanyarray(args.tdomain)
    ts, te = tdomain * grb.day2sec
    nus = [1e9, 8.5e9, 1e12, 1e15, 1e17]
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
    if args.sring:
        Fnu_ring2 = grb.fluxDensity(t, nu, **Zring2)
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
        if front_part(nus[idx]) == 1:
            nus_label = r'$10^{%d} \rm Hz$'%(order_of_mag(nus[idx]))
        else:
            nus_label = r'$%1.f \times 10^{%d} \rm Hz$'%(front_part(nus[idx]), order_of_mag(nus[idx]))
        reg_lines[idx], = ax.plot(tday[:, idx], Fnu_jet[:, idx], color=colors[idx], label = nus_label)
        ax.plot(tday[:, idx], Fnu_ring[:, idx], linestyle='--', color=colors[idx])
        if args.sring:
            ax.plot(tday[:, idx], Fnu_ring2[:, idx], linestyle=':', color=colors[idx])
            
    ax.set(xscale='log', xlabel=r'$t [\rm s]$',
        yscale='log', ylabel=r'$F_\nu [\rm mJy]$')

    ring_line1  = Line2D([0],[0],linestyle='--', color='gray', label='conical ring')
    ring_lines = [ring_line1]
    if args.sring:
        ring_line2 = Line2D([0],[0],linestyle=':', color='gray', label='structured ring')
        ring_lines += [ring_line2]
    ax.set_title(r"$\rm{Transient \ Afterglows}$, " + eiso_label)
    ax.legend(handles=[*reg_lines, *ring_lines])
    # ax.legend()
    fig.tight_layout()
    print(f"Saving figure {args.file_name + '.pdf'}")
    fig.savefig(args.file_name + ".pdf")
    plt.show()
    plt.close(fig)

if __name__ == '__main__':
    main()