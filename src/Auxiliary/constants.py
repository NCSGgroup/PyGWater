import numpy as np


class EarthConstant_2:
    ggg = 6.67384e-11  # % Newton's constant (SI units)
    radiusm = 6371 * 1e3  # % Radius of the Earth (m)
    grav = 9.8046961  # % Surface gravity (m/s/s)
    rhow = 1000  # % density of pure water(kg/m^3)
    rhoear = 3.0 * grav / 4.0 / ggg / np.pi / radiusm  # % Average Earth density (kg/m^3)
    Mass = grav * radiusm ** 2 / ggg


class EarthConstant:
    ggg = 6.67384e-11  # % Newton's constant (SI units)
    grav = 9.80665  # % Surface gravity (m/s/s)
    rhow = 1025.0  # % density of pure water(kg/m^3)
    rhoear = 5517  # unit[kg/m3]
    radiusm = 6378136.3  # unit[m]
    Mass = grav * radiusm ** 2 / ggg


class EarthConstant_3:
    density_earth = 5517  # unit[kg/m3]
    radius_earth = 6378136.3  # unit[m]
    # GM = 3.9860044150E+14  # unit[m3/s2]
    GM = 398940628222199.44  # unit[m3/s2]
    '''gravity constant g defined by WMO'''
    g_wmo = 9.80665
    # g_wmo = 9.7
    ''' water density'''
    density_water = 1025.0
    # density_water = 1000.0


def demo():
    print(EarthConstant.Mass*EarthConstant.ggg)
    print(EarthConstant.rhoear)
    EarthConstant.rhoear = 9
    print(EarthConstant.rhoear)
    pass


if __name__ == '__main__':
    demo()
