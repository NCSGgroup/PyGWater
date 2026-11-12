import numpy as np


class EarthConstant:
    """
    Single source of truth for the physical constants used by ALL basis functions
    (spherical harmonics, finite/infinite point-mass LGF, D-LGF, CD-LGF).

    Consistency requirement (Yang et al., 2025, JGR Solid Earth, Eqs. 2, 7, 11, 21):
        M = 4/3 * pi * a^3 * rho_e
    so that  a/M * ... (Green's function)  and  3*rho_w/rho_e * ... (spherical harmonic)
    describe the same Earth. Previously the point-load Green's function used M = g*a^2/G
    (implied rho_e ~ 5500 kg/m3) while the SH and disk paths used rho_e = 5517 kg/m3,
    which produced a spurious 0.31% scale difference between the bases.

    Water density follows the paper (Sect. 3.1): fresh water, 1000 kg/m3. This is also
    the convention behind GRACE mascon/EWH products such as CSR-Mascon.
    If another value is wanted, change it HERE only.
    """
    ggg = 6.67384e-11  # Newton's constant (SI units)
    grav = 9.80665  # standard surface gravity (m/s/s), WMO
    rhow = 1000.0  # density of fresh water (kg/m^3), see Yang et al. (2025) Sect. 3.1
    rhoear = 5517.0  # mean Earth density (kg/m^3)
    radiusm = 6378136.3  # Earth radius (m), the GRACE Level-2 convention
    Mass = 4.0 / 3.0 * np.pi * radiusm ** 3 * rhoear  # Earth mass consistent with rhoear and radiusm (kg)
    GM = ggg * Mass  # (m^3/s^2), consistent with Mass


class EarthConstant_3:
    """
    Backward-compatible alias used by the spherical-harmonic conversion (ConvertSHC).
    All values are taken from EarthConstant to guarantee consistency.
    """
    density_earth = EarthConstant.rhoear  # unit[kg/m3]
    radius_earth = EarthConstant.radiusm  # unit[m]
    GM = EarthConstant.GM  # unit[m3/s2]
    g_wmo = EarthConstant.grav  # gravity constant g defined by WMO
    density_water = EarthConstant.rhow  # water density


def demo():
    print('M  =', EarthConstant.Mass)
    print('GM =', EarthConstant.GM)
    print('3/rho_e == 4*pi*a^3/M :', 3 / EarthConstant.rhoear,
          4 * np.pi * EarthConstant.radiusm ** 3 / EarthConstant.Mass)
    pass


if __name__ == '__main__':
    demo()
