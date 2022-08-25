import jet_ext
import numpy as np 
import matplotlib.pyplot as plt 
from enum import Enum 

global_n = 100

c = 2.99792458e10
me = 9.1093897e-28
mp = 1.6726231e-24
h = 6.6260755e-27
hbar = 1.05457266e-27
ee = 4.803e-10
sigmaT = 6.65e-25

Msun = 1.98892e33
cgs2mJy = 1.0e26
mJy2cgs = 1.0e-26
deg2rad = np.pi/180.0
rad2deg = 180.0/np.pi
day2sec = 86400.0
sec2day = 1.0/day2sec
parsec = 3.0857e18
Hz2eV = 4.13566553853599e-15
eV2Hz = 1.0/Hz2eV


class Jet(Enum):
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
    
default_params = {
    "jet_type":  -1,
    "spec_type":  0,
    "theta_obs":  0.0,
    "E_iso_core":  1.0e53,
    "theta_h_core":  0.1,
    "theta_h_wing":  0.4,
    "b": 0.0,
    "L0": 0.0,
    "q": 0.0,
    "ts": 0.0, 
    "n_0": 1.0,
    "p": 2.5,
    "epsilon_E": 0.1,
    "epsilon_B": 0.01,
    "ksi_N": 1.0, 
    "d_L": 1.0e27,

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
    
    "jetType": Jet.TopHat,
    "specType": 0,
    "thetaObs": 0,
    "E0": 1e53,
    "thetaCore": 0.1,
    "n0": 1,
    "epsilon_e": 0.1,
    "epsilon_b": 0.1,
    "xi_N": 1.0,
    "theta_wing": 0.1,
    "thetaWing": 0.1,
    "mask": np.zeros(global_n),
    "nmask": global_n,   
}

def fluxDensity(t, nu, *args, **kwargs):
    """
    Compute the flux density F_nu of a GRB afterglow.

    Utiliizes the single shell approximation described in Ryan et al 2020
    to compute the synchrotron emission from the forward shock of a blast
    wave at the specified times and frequencies.

    ``fluxDensity`` takes many model-dependent parameters which are specified
    via keyword arguments. It is recommended to collect these parameters in a
    dictionary ``Z`` and call fluxDensity() by::

        Fnu = fluxDensity(t, nu, **Z)

    Alternately model parameters may be specified as positional arguments,
    although this interface is **not recommended** and may be deprecated in
    the future. To call a jetted model with positional and keyword arguments::

        Fnu = fluxDensity(t, nu, jetType, specType, thetaObs, E0, thetaCore,
                          thetaWing, b, L0, q, ts, n0, p, epsilon_e, epsilon_B,
                          xi_N, d_L, **Z)

    To call a spherical refreshed shock model with positional arguments::

        Fnu = fluxDensity(t, nu, Jet.Spherical, specType, uMax, uMin, Er, k,
                          MFast_solar, L0, q, ts, n0, p, epsilon_e, epsilon_B,
                          xi_N, d_L, **Z)


    Parameters 
    ----------
    t : array_like or scalar
        Time since burst in observer frame, measured in seconds.
    nu : array_like or scalar
        Frequency of flux in observer frame, measured in Hz, same size as t.
    jetType : int
        Code for type of Jet. Model codes are available in ``afterglowpy.jet``
        and include: ``Jet.TopHat``, ``Jet.Cone``, ``Jet.Gaussian``,
        ``Jet.PowerLaw``, ``Jet.GaussianCore``, ``Jet.PowerLawCore``,
        ``Jet.Spherical`, and ``Jet.Ring''`.
    specType : int
        Code for type of spectrum.  Options are: 0 broken power law
        (Ryan+ 2020), 1 broken power law w/ inverse Compton cooling. Default: 0
    thetaObs : float
        Viewing angle in radians. Jet models only.
    E0: float
        Isotropic-equivalent energy along the jet axis in ergs. Jet models
        only.
    thetaCore: float
        Half opening angle of jet core in radians. Jet models only.
    thetaWing: float
        Outer truncation angle of the jet in radians. Ignored by
        ``Jet.TopHat``, jet models only.
    b: float
        Power law index of jet angular energy distribution. Only used by
        ``Jet.PowerLaw`` and ``Jet.PowerLawCore``.
    n0 : float
        Number density of protons in circumburst medium in cm^{-3}.
    p : float
        Power law index of relativistic electron energy distribution,
        generally p > 2.
    epsilon_e : float
        Fraction of thermal energy in relativistic electrons, epsilon_e <= 1.
    epsilon_B : float
        Fraction of thermal energy in magnetic field, epsilon_B <= 1.
    xi_N : float
        Fraction of electrons that get accelerated, xi_N <= 1.
    d_L : float
        Luminosity distance to burst, in cm.
    z : float, optional
        Redshift of burst, defaults to 0.

    L0: float, optional
        Luminosity of energy injection, in erg/s.  Default 0.0.
    q: float, optional
        Power law index of energy injection: L = L0 (t/t0)^{-q}, t0 = 1 ks.
        Default 0.0.
    ts: float, optional
        Time energy injection ends in burster frame, in seconds. Default 0.0.
    g0 : float, optional
        EXPERIMENTAL.  Initial Lorentz factor of outflow along jet axis,
        defaults to -1 (unset, jet has deceleration radius 0). Do not use with
        jet spreading enabled.

    uMax : float
        Maximum 4-velocity of outflow. Only for spherical models.
    uMin : float
        Minimum 4-velocity of outflow. Only for spherical models
    Er : float
        Normalization of outflow's energy distribution in ergs. Only for
        spherical models.
        E(u>U) = Er * U^{-k}
    k : float
        Power law index of outflow's energy distribution. Only for spherical
        models
    MFast_solar : float
        Mass of material at u_max in solar masses. Only for spherical models.

    Other Parameters
    ----------------

    spread : {True, False}, optional
        Whether to include jet spreading. Defaults to True.
    counterjet : {'True', 'False'}, optional
        Whether to include counterjet emission. Defaults to False.
    tRes : int, optional
        Time resolution, number of points per decade in t, for shock evolution.
        Defaults to 1000.
    latRes : int, optional
        Lateral resolution of structured jet, number of conical sections per
        thetaCore-sized interval. Defaults to 5.
    intType : int, optional
        Integration scheme to use when computing flux. Defaults to
        ``Jet.Cadre``. Changing this may result in longer run times or larger
        than expected numerical errors.
    rtolStruct : float, optional
        Overall relative tolerance of flux integration for structured jets.
        Defaults to 1.0e-2.
    rtolTheta : float, optional
        Relative tolerance of flux integration over theta within each
        conical section. Defaults to 1.0e-2.
    rtolPhi : float, optional
        Relative tolerance of flux integration over phi within
        each conical section. Defaults to 1.0e-2.
    NPhi : int, optional
        Maximum number of evaluations to perform in phi direction during
        numerical integration. Default 1000.
    NTheta : int, optional
        Maximum number of evaluations to perform in theta direction during
        numerical integration. Default 1000.

    Returns
    -------
    Fnu: array
        The flux density F_nu in the observer frame, same shape as t and nu.

    Raises
    ------

    ValueError
        If t, nu are the wrong shape or arguments take illegal values.
    """

    argsDict = parseArgs(args, kwargs)

    t, nu = checkTNu(t, nu)

    jetType = argsDict['jetType']
    checkJetArgs(**argsDict)
    if jetType == Jet.Spherical:
        checkCocoonArgs(**argsDict)
    else:
        checkJetArgs(**argsDict)

    # arguments are good, full steam ahead!
    z = argsDict.pop('z') if 'z' in argsDict else 0.0

    tz  = t / (1+z)
    nuz = nu * (1+z)

    # Default spreading method
    if 'spread' in argsDict:
        if argsDict['spread'] == True:
            if jetType == -2 and 'thetaCoreGlobal' in argsDict:
                argsDict['spread'] = 8
            else:
                argsDict['spread'] = 7

    # Intercept background luminosities
    # This was a bad idea to add to this function, but is kept for 
    # backwards compatibility. Please don't use these.
    LR = argsDict.pop('LR') if 'LR' in argsDict else 0.0
    LO = argsDict.pop('LO') if 'LO' in argsDict else 0.0
    LX = argsDict.pop('LX') if 'LX' in argsDict else 0.0
    tAdd = argsDict.pop('tAdd') if 'tAdd' in argsDict else 0.0
    
    # timeA = time.time()
    argsDict['ta'] = t.flat[0]
    argsDict['tb'] = t.flat[-1]
    Fnu = np.empty(tz.shape)
    Fnu.flat[:] = jet_ext.calc_flux_density_wrapper(
        jet_type  = -1, 
        spec_type = 1,
        t         = tz.flatten(),
        nu        = nuz.flatten(),
        fnu       = Fnu.flatten(),
        n         = global_n,
        fp        = argsDict
    )
    print(Fnu)
    zzz = input('')
    # if argsDict['jetType'] == Jet.Spherical:
    #     Fnu.flat[:] = cocoon.fluxDensity(tz.flat, nuz.flat, **argsDict)
    # else:
    #     Fnu.flat[:] = Jet.fluxDensity(tz.flat, nuz.flat, **argsDict)
    # timeB = time.time()
    # print("Eval took: {0:f} s".format(timeB - timeA))

    
    # Adding background luminosities.
    L_to_flux = cgs2mJy / ( 4 * np.pi * argsDict['d_L']**2)


    if LR > 0.0:
        rad = (nuz < 3.0e11) & (tz > tAdd)  # radio < 300 GHz
        Lnu = LR / 1.0e10  # 10 GHz bandwidth
        Fnu[rad] += Lnu*L_to_flux
    if LO > 0.0:
        # 300GHz<opt<100eV
        opt = (nuz >= 3.0e11) & (nuz < 100*cocoon.eV2Hz) & (tz > tAdd)
        Lnu = LO * 2.32478e-5 / cocoon.c # 2324.78 A bandwidth
        Fnu[opt] += Lnu*L_to_flux
    if LX > 0.0:
        xry = (nuz >= 100*cocoon.eV2Hz) & (tz > tAdd)  # xray > 100eV
        Lnu = LX / (9.7e3 * cocoon.eV2Hz)  # 9.7 keV bandwidth
        Fnu[xry] += Lnu*L_to_flux

    # K-correct the flux
    Fnu *= 1+z

    return Fnu

def checkTNu(t, nu):
    # Make sure t and nu are array_like or castable to an array.
    t = np.atleast_1d(t)
    nu = np.atleast_1d(nu)

    # Check shapes, if scalars make into right size array
    if t.shape != nu.shape:
        if t.shape == (1, ):
            T = t[0]
            t = np.empty(nu.shape)
            t[:] = T
        elif nu.shape == (1, ):
            NU = nu[0]
            nu = np.empty(t.shape)
            nu[:] = NU
        else:
            raise ValueError("t and nu must be same shape or scalars")

    return t, nu


def checkThetaPhiTNu(theta, phi, t, nu):
    # Make sure args are array_like or castable to an array.
    theta = np.atleast_1d(theta)
    phi = np.atleast_1d(phi)
    t = np.atleast_1d(t)
    nu = np.atleast_1d(nu)

    shape = (1, )
    if shape != theta.shape:
        shape = theta.shape
    elif shape != phi.shape:
        shape = phi.shape
    elif shape != t.shape:
        shape = t.shape
    else:
        shape = nu.shape

    if theta.shape != shape:
        if theta.shape == (1, ):
            TH = theta[0]
            theta = np.empty(shape)
            theta[:] = TH
        else:
            msg = "theta must be scalar or same shape as phi, t, nu"
            raise ValueError(msg)

    if phi.shape != shape:
        if phi.shape == (1, ):
            PH = phi[0]
            phi = np.empty(shape)
            phi[:] = PH
        else:
            msg = "phi must be scalar or same shape as theta, t, nu"
            raise ValueError(msg)

    if t.shape != shape:
        if t.shape == (1, ):
            T = t[0]
            t = np.empty(shape)
            t[:] = T
        else:
            msg = "t must be scalar or same shape as theta, phi, nu"
            raise ValueError(msg)

    if nu.shape != shape:
        if nu.shape == (1, ):
            NU = nu[0]
            nu = np.empty(shape)
            nu[:] = NU
        else:
            msg = "nu must be scalar or same shape as theta, phi, t"
            raise ValueError(msg)

    return theta, phi, t, nu


def checkJetArgs(**argsDict):

    # for _, x in argsDict.items():
    #     if not np.isfinite(x):
    #         raise ValueError("All parameters must be finite")


    jetType = argsDict['jetType']
    specType = argsDict['specType']

    theta_obs = argsDict['thetaObs']
    E0 = argsDict['E0']
    theta_c = argsDict['thetaCore']
    n0 = argsDict['n0']
    p = argsDict['p']
    epse = argsDict['epsilon_e']
    epsB = argsDict['epsilon_B']
    xiN = argsDict['xi_N']
    dL = argsDict['d_L']

    # More-or-less universal bounds
    if theta_obs < 0.0 or theta_obs > 0.5*np.pi:
        raise ValueError("theta_obs must be in [0.0, pi/2]")
    if E0 <= 0.0:
        raise ValueError("E0 must be positive")
    if jetType != Jet.Cone and (theta_c <= 0.0 or theta_c > 0.5*np.pi):
        raise ValueError("theta_c must be in (0.0, pi/2]")
    if jetType == Jet.Cone and (theta_c < 0.0 or theta_c > 0.5*np.pi):
        raise ValueError("theta_c must be in [0.0, pi/2]")
    if n0 <= 0.0:
        raise ValueError("n0 must be positive")
    if specType != 2 and p <= 2.0:
        raise ValueError("p must be in (2, inf)")
    if specType == 2 and p <= 1.0:
        raise ValueError("p must be in (1, inf)")
    if epse <= 0.0 or epse > 1.0:
        raise ValueError("epsilon_e must be in (0, 1]")
    if epsB <= 0.0 or epsB > 1.0:
        raise ValueError("epsilon_B must be in (0, 1]")
    if xiN <= 0.0 or xiN > 1.0:
        raise ValueError("xi_N must be in (0, 1]")
    if dL <= 0.0:
        raise ValueError("d_L must be positive")

    # Bounds on optional parameters

    if jetType != Jet.TopHat:
        if 'thetaWing' not in argsDict:
            raise KeyError('This jet type requires thetaWing')
        else:
            theta_w = argsDict['thetaWing']
            if (theta_w <= 0.0 or theta_w > 0.5*np.pi):
                raise ValueError("thetaWing must be in (0.0, pi/2]")

    if jetType == Jet.PowerLaw or jetType == Jet.PowerLawCore:
        if 'b' not in argsDict:
            raise KeyError('This jet type requires b')
        else:
            b = argsDict['b']
            if (b <= 0.0):
                raise ValueError("b must be positive")

    # Energy Injection
    if 'L0' in argsDict:
        L0 = argsDict['L0']
        if L0 < 0.0:
            raise ValueError("L0 must be non-negative")
    if 'ts' in argsDict:
        ts = argsDict['ts']
        if ts < 0.0:
            raise ValueError("ts must be non-negative")

    # Additional Luminosity
    if 'LR' in argsDict:
        LR = argsDict['LR']
        if LR < 0.0:
            raise ValueError("LR must be non-negative")
    if 'LO' in argsDict:
        LO = argsDict['LO']
        if LO < 0.0:
            raise ValueError("LO must be non-negative")
    if 'LX' in argsDict:
        LX = argsDict['LX']
        if LX < 0.0:
            raise ValueError("LX must be non-negative")

    if 'z' in argsDict:
        if argsDict['z'] < 0.0:
            raise ValueError("z must be non-negative")

    # Model Specific bounds
    if jetType == Jet.Cone and argsDict['thetaCore'] > argsDict['thetaWing']:
        raise ValueError("thetaWing must be larger than thetaCore"
                         "for cone model")

    return


def checkCocoonArgs(**argsDict):

    for _, x in argsDict.items():
        if not np.isfinite(x):
            raise ValueError("All parameters must be finite")

    specType = argsDict['specType']

    u_max = argsDict['uMax']
    u_min = argsDict['uMin']
    Er = argsDict['Er']
    MFast = argsDict['MFast_solar']
    n0 = argsDict['n0']
    p = argsDict['p']
    epse = argsDict['epsilon_e']
    epsB = argsDict['epsilon_B']
    xiN = argsDict['xi_N']
    dL = argsDict['d_L']

    if u_max <= 0.0:
        raise ValueError("u_max must be positive")
    if u_min <= 0.0:
        raise ValueError("u_min must be positive")
    if Er <= 0.0:
        raise ValueError("Er must be positive")
    if MFast <= 0.0:
        raise ValueError("MFast must be positive")
    if n0 <= 0.0:
        raise ValueError("n0 must be positive")
    if specType != 2 and p <= 2.0:
        raise ValueError("p must be in (2, inf)")
    if specType == 2 and p <= 1.0:
        raise ValueError("p must be in (1, inf)")
    if epse <= 0.0 or epse > 1.0:
        raise ValueError("epsilon_e must be in (0, 1]")
    if epsB <= 0.0 or epsB > 1.0:
        raise ValueError("epsilon_B must be in (0, 1]")
    if xiN <= 0.0 or xiN > 1.0:
        raise ValueError("xi_N must be in (0, 1]")
    if dL <= 0.0:
        raise ValueError("d_L must be positive")

    if 'z' in argsDict:
        if argsDict['z'] < 0.0:
            raise ValueError("z must be non-negative")

    # Energy Injection
    if 'L0' in argsDict:
        L0 = argsDict['L0']
        if L0 < 0.0:
            raise ValueError("L0 must be non-negative")
    if 'ts' in argsDict:
        ts = argsDict['ts']
        if ts < 0.0:
            raise ValueError("ts must be non-negative")

    # Additional Luminosity
    if 'LR' in argsDict:
        LR = argsDict['LR']
        if LR < 0.0:
            raise ValueError("LR must be non-negative")
    if 'LO' in argsDict:
        LO = argsDict['LO']
        if LO < 0.0:
            raise ValueError("LO must be non-negative")
    if 'LX' in argsDict:
        LX = argsDict['LX']
        if LX < 0.0:
            raise ValueError("LX must be non-negative")

    return


def parseArgs(args, kwargs):
    """
    Parse the arguments to fluxDensity() or intensity(). Supports both
    positional and keyword arguments for now.
    """

    argsDict = kwargs.copy()

    # If there were no extra positional arguments, things are easy.
    if len(args) == 0:
        return argsDict

    # Recommend users not to do this
    warnings.warn("Positional argument inferface to afterglowpy"
                  + " may not be supported in the future. Use keyword"
                  + " arguments instead, see documentation for details.",
                  FutureWarning)

    # Now for the fun part

    jetKeys = ['jetType', 'specType', 'thetaObs', 'E0', 'thetaCore',
               'thetaWing', 'b', 'L0', 'q', 'ts', 'n0', 'p', 'epsilon_e',
               'epsilon_B', 'xi_N', 'd_L', 'g0', 'LR', 'LO', 'LX', 'tAdd', 'z']
    sphKeys = ['jetType', 'specType', 'uMax', 'uMin', 'Er',
               'k', 'MFast_solar', 'L0', 'q', 'ts', 'n0', 'p', 'epsilon_e',
               'epsilon_B', 'xi_N', 'd_L', 'g0', 'LR', 'LO', 'LX', 'tAdd', 'z']

    jetType = args[0]

    # if jetType == Jet.Spherical:
    #     argsDict.update(zip(sphKeys, args))
    # else:
    #     argsDict.update(zip(jetKeys, args))

    return argsDict

t    = np.linspace(0.1, 1000, global_n) * day2sec
nu   = np.ones_like(t) * 1e9
tday = t * sec2day
a    = fluxDensity(t, nu, **default_params)

plt.loglog(tday, a)
plt.show()