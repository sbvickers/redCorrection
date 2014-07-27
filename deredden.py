#! /usr/bin/env python3

# --------------------------------------------
#
#   function to deredden flux data using extinction
#   law of Cardelli, Clayton & Mathis 1989 (CCM 1989)
#   with updated coefficients from O'Donnell 1994
#
#   function takes wavelength, flux and A_v
#   and returns a dereddened flux vector
#
# --------------------------------------------

import numpy as np

def law_ir(a, b, wave):
    """
        returns two vectors of extinction law values
        for the ir regime 0.909um - ...
    """

    a[wave <= 1.1] = 0.574 * wave[wave <= 1.1]**(1.61)
    b[wave <= 1.1] = -0.527 * wave[wave <= 1.1]**(1.61)
    
    return a, b

def law_nir(a, b, wave):
    """ 
        returns two vectors of extinction law values
        for the fir regime 0.303um - 0.909um
    """

    y = wave[(wave > 1.1) & (wave <= 3.3)] - 1.82

# coefficients of CCM 1898
    #a[(wave > 1.1) & (wave <= 3.3)] = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
    #b[(wave > 1.1) & (wave <= 3.3)] = 0 + 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7

# coefficients of O'Donnell 1984
    a[(wave > 1.1) & (wave <= 3.3)] = 1 + 0.104*y - 0.609*y**2 + 0.701*y**3 + 1.137*y**4 - 1.718*y**5 - 0.827*y**6 + 1.647*y**7 - 0.505*y**8
    b[(wave > 1.1) & (wave <= 3.3)] = 0 + 1.952*y + 2.908*y**2 - 3.989*y**3 - 7.985*y**4 + 11.102*y**5 + 5.491*y**6 - 10.805*y**7 + 3.347*y**8
   
    return a, b

def law_uv(a, b, wave):
    """
       returns two vectors of extinction law values
       for the uv regime 0.125um - 0.303um
    """
    f_a = np.empty(len(wave[(wave >= 3.3) & (wave <= 8.0)]), dtype=float)
    f_a[:] = np.NaN
    
    f_b = np.empty(len(wave[(wave >= 3.3) & (wave <= 8.0)]), dtype=float)
    f_b[:] = np.NaN

    y = wave[(wave >=3.3) & (wave <= 8.0)]

    if len(y[y >= 5.9]) > 0:
        f_a[y >= 5.9] = -0.04473*(y[y >= 5.9] - 5.9)**2 - 0.009779*(y[y >= 5.9] - 5.9)**3
        f_b[y >= 5.9] = -0.2130*(y[y >= 5.9] - 5.9)**2 - 0.1207*(y[y >= 5.9] - 5.9)**3

    f_a[y < 5.9] = 0
    f_b[y < 5.9] = 0

    a[(wave >= 3.3) & (wave <= 8.0)] = 1.752 - 0.316 * y - 0.104/((y-4.67)**2 + 0.341) + f_a
    b[(wave>= 3.3) & (wave <= 8.0)] = -3.090 + 1.825 * y + 1.206/((y-4.62)**2 + 0.263) + f_b

    return a, b

def law_fuv(a, b, wave):
    """
        returns two vectors of extinction law values
        for the far-uv regime 0.1um - 0.125um
    """
    
    y = wave[(wave > 8.0) & (wave <= 10)] - 8

    a[(wave > 8.0) & (wave <= 10)] = -1.073 - 0.628*y + 0.137*y**2 - 0.070*y**3
    b[(wave > 8.0) & (wave <= 10)] = 13.670 + 4.257*y - 0.420*y**2 + 0.334*y**3

    return a, b

def dered(wave, flux, A_v, R_v = 3.1):
    """
        dereddens a spectra from 0.1 microns through
        to the far IR
    """
    wave = 1 / wave
   
    a = np.empty(len(wave), dtype=float)
    a[:] = np.NaN

    b = np.empty(len(wave), dtype=float)
    b[:] = np.NaN

# law_ir is for lambda > 0.909 microns

    a, b = law_uv(a, b, wave)
    a, b = law_fuv(a, b, wave)
    a, b = law_nir(a, b, wave)
    a, b = law_ir(a, b, wave)
    
    A_lambda = A_v * (a + b/R_v)

    return flux * 10**(0.4*A_lambda)

