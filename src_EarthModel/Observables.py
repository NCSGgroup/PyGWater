"""
Pluggable definition of loading observables.

Every linear elastic response to a surface load is described here ONCE, by
  (1) a per-degree factor X_n that maps the load's EWH spherical-harmonic coefficients [m] onto the observable, and
  (2) how the observable is evaluated from the resulting coefficients:
        'value' -> Y_nm                       (scalar observables)
        'north' -> -dY_nm/dtheta              (north component of a surface gradient; theta = colatitude)
        'east'  -> 1/sin(theta) dY_nm/dlambda (east component of a surface gradient)
Spherical harmonics (SHC.synthesis_observable / sh_synthesis) and Slepian functions (SlepianBasis) read this table,
so adding an observable here makes it available in both. Green's-function kernels (PointLoad, DiskLoad) need a
dedicated closed form/asymptotic treatment per observable and are implemented individually.

Notation: F_n = 3 rho_w / (rho_e (2n+1)),  N_n = F_n (1 + k_n) E_n is the geoid, V_n = g N_n the potential,
g = GM/a^2 of the model Earth (EarthConstant), so that g/a * F_n = 4 pi G rho_w / (2n+1) exactly.

  Observable (EnumClasses.Displacement)  factor X_n                                         unit   kind
  -------------------------------------  -------------------------------------------------  -----  -----------
  Vertical                               F_n h_n                                            m      value
  North / East (horizontal)              F_n l_n                                            m      north/east
  Geoheight (geoid)                      F_n (1 + k_n)                                      m      value
  Gravity (surface gravimeter)           g/a F_n [newton_n + (n+1) k_n - 2 h_n]             m/s^2  value
  GravityDisturbance (at height H)       g/a F_n (n+1) (1 + k_n) (a/(a+H))^(n+2)            m/s^2  value

Surface gravity (Farrell, 1972; positive = increase of gravity) = Newtonian attraction of the load + attraction of the
redistributed Earth masses ((n+1) k_n) + free-air effect of the vertical displacement (-2 h_n). The Newtonian term of
a thin surface layer depends on the side of the layer the sensor is on:
    newtonian='below' : newton_n = -n      (sensor beneath the load, e.g. atmosphere; Farrell's convention and the
                                            convention of the Green's functions of Wang et al., 2012)
    newtonian='above' : newton_n = n + 1   (sensor on top of the load, e.g. groundwater/soil moisture underneath)
    newtonian='none'  : newton_n = 0       (elastic/deformation part only)
'below' and 'above' differ only by the local term 4 pi G rho_w E at the sensor (about 83.9 microGal per metre of
water, i.e. twice the Bouguer-plate value 2 pi G rho_w = 41.9 microGal/m); away from the load (Green's functions at
theta > 0) they are identical.
1 microGal = 1e-8 m/s^2.

InSAR line of sight: line_of_sight() / look_vector() project (east, north, up) onto the ground-to-satellite unit
vector (positive = motion towards the satellite).
"""

import numpy as np

from src_auxiliary.constants import EarthConstant
from src_auxiliary.EnumClasses import Displacement, LLN_variable

'''alias: the enum lists observables, not only displacements'''
Observable = Displacement


def kind(variable: Displacement):
    """How the observable is evaluated from coefficients: 'value', 'north' or 'east'."""
    if variable in (Displacement.Vertical, Displacement.Geoheight, Displacement.Gravity,
                    Displacement.GravityDisturbance):
        return 'value'
    if variable == Displacement.North:
        return 'north'
    if variable == Displacement.East:
        return 'east'
    raise NotImplementedError('%s is not a supported observable (use North/East for horizontal displacement)'
                              % variable)


def love_factor(variable: Displacement, lln, n, newtonian='below', height=0.0):
    """
    Per-degree factor X_n from load EWH coefficients [m] to the observable (see module docstring).
    :param variable: EnumClasses.Displacement member
    :param lln: LoveNumber object (frame/Earth model chosen by the user)
    :param n: array of degrees
    :param newtonian: 'below' | 'above' | 'none', surface gravity only
    :param height: height above the reference sphere [m], GravityDisturbance only
    :return: array like n
    """
    n = np.asarray(n)
    LLN = lln.LLN
    if len(LLN[LLN_variable.h]) < n.max() + 1:
        raise ValueError('Love numbers are needed up to degree %d' % n.max())
    h, l, k = LLN[LLN_variable.h][n], LLN[LLN_variable.l][n], LLN[LLN_variable.k][n]
    a = EarthConstant.radiusm
    g = EarthConstant.GM / a ** 2  # gravity of the model Earth (consistent with rho_e, a): g/a * F_n = 4 pi G rho_w/(2n+1)
    F = 3 * EarthConstant.rhow / (EarthConstant.rhoear * (2 * n + 1))

    if variable == Displacement.Vertical:
        return F * h
    if variable in (Displacement.North, Displacement.East):
        return F * l
    if variable == Displacement.Geoheight:
        return F * (1 + k)
    if variable == Displacement.Gravity:
        newton = {'below': -n, 'above': n + 1, 'none': 0 * n}[newtonian]
        return g / a * F * (newton + (n + 1) * k - 2 * h)
    if variable == Displacement.GravityDisturbance:
        return g / a * F * (n + 1) * (1 + k) * (a / (a + height)) ** (n + 2)
    raise NotImplementedError(variable)


def sh_synthesis(load_cqlm, load_sqlm, lln, variable: Displacement, lat, lon, **kwargs):
    """
    Spherical-harmonic synthesis of an observable on a regular grid, from the LOAD (EWH [m]) coefficients.
    :param load_cqlm, load_sqlm: (L+1, L+1) or (ntime, L+1, L+1)
    :param lat, lon: 1d grid coordinates [degree]
    :param kwargs: newtonian=..., height=... (see love_factor)
    :return: (nlat, nlon) or (ntime, nlat, nlon)
    """
    from src_Basis.SphericalHarmonic.Harmonic import Harmonic
    from src_auxiliary.EnumClasses import PhysicalDimensions

    load_cqlm, load_sqlm = np.asarray(load_cqlm), np.asarray(load_sqlm)
    L = load_cqlm.shape[-1] - 1
    X = love_factor(variable, lln, np.arange(L + 1), **kwargs)[:, None]
    special = {'value': None, 'north': PhysicalDimensions.HorizontalDisplacementNorth,
               'east': PhysicalDimensions.HorizontalDisplacementEast}[kind(variable)]
    return Harmonic(np.asarray(lat), np.asarray(lon), L, option=1).synthesis(load_cqlm * X, load_sqlm * X,
                                                                           special_type=special)


def look_vector(incidence, heading, right_looking=True):
    """
    Ground-to-satellite unit vector (east, north, up) of a SAR acquisition.
    :param incidence: incidence angle [degree]
    :param heading: satellite heading, clockwise from north [degree]
    :param right_looking: True for right-looking SAR (Sentinel-1, ALOS-2, ...)
    :return: array (3,) or (3, N)
    """
    inc, hd = np.deg2rad(incidence), np.deg2rad(heading)
    beta = hd - np.pi / 2 if right_looking else hd + np.pi / 2  # azimuth of the ground -> satellite direction
    return np.array([np.sin(inc) * np.sin(beta), np.sin(inc) * np.cos(beta), np.cos(inc)])


def line_of_sight(east, north, up, look):
    """
    Line-of-sight displacement, positive towards the satellite (range decrease).
    :param east, north, up: displacement components (same shape)
    :param look: (3,) or (3, N) ground-to-satellite unit vector, see look_vector()
    """
    look = np.asarray(look)
    return look[0] * east + look[1] * north + look[2] * up
