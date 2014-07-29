# --------------------------------------------
#
#   function to deredden flux data using extinction
#   law of Cardelli, Clayton & Mathis 1989 (CCM 1989)
#   with updated near-IR coefficients from O'Donnell 1994
#
#   function takes wavelength, flux, E(B-V) and R_v (optional)
#   and returns a dereddened numpy array of fluxes.
#
# --------------------------------------------

import numpy as np

def law_ir(a, b, wave):
    """
        Returns a and b vectors of extinction law values
        for the IR regime (wavelength: 0.909um - ...; wavenumber: <= 1.1).

        This is just a power law extinction curve from 0.909um onwards.

        Parameters
        ----------
                a : ndarray
                Array of values used to correct for extinction. Only the section
                of the 'a' vector that lies in the IR and onwards are altered in
                this function.

                b : ndarray
                Array of values used to correct for extinction. Only the section
                of the 'b' vector that lies in the IR and onwards are altered in
                this function.

                wave : ndarray
                The wavelengths of the fluxes being corrected. Used to determine
                which 'a' and 'b' vector values to alter.

        Returns 
        ----------
                a : ndarray
                Modified 'a' vector values.

                b : ndarray
                Modified 'b' vector values.
    """

    a[wave <= 1.1] = 0.574 * wave[wave <= 1.1]**(1.61)
    b[wave <= 1.1] = -0.527 * wave[wave <= 1.1]**(1.61)
    
    return a, b

def law_nir(a, b, wave):
    """ 
        Returns a and b vectors of extinction law values for the far-IR 
        regime (wavelength: 0.303um - 0.909um; wavenumber: 1.1 - 3.3).

        Parameters
        ----------
                a : ndarray
                Array of values used to correct for extinction. Only the section
                of the 'a' vector that lies in the far-IR range are altered in
                this function.

                b : ndarray
                Array of values used to correct for extinction. Only the section
                of the 'b' vector that lies in the far-IR range are altered in
                this function.

                wave : ndarray
                The wavelengths of the fluxes being corrected. Used to determine
                which 'a' and 'b' vector values to alter.

        Returns 
        ----------
                a : ndarray
                Modified 'a' vector values.

                b : ndarray
                Modified 'b' vector values.
    """
    cond = '(wave > 1.1) & (wave <= 3.3)'

    y = wave[eval(cond)] - 1.82

    # ***********************************************
    # Original near-IR coefficients of CCM 1989
    # ***********************************************

    #a[(wave > 1.1) & (wave <= 3.3)] = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7

    #b[(wave > 1.1) & (wave <= 3.3)] = 0 + 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    # ***********************************************


    # ***********************************************
    # Updated near-IR coefficients of O'Donnell 1984
    # ***********************************************
    a[eval(cond)] = 1 + 0.104*y - 0.609*y**2 + 0.701*y**3 + 1.137*y**4 - 1.718*y**5 - 0.827*y**6 + 1.647*y**7 - 0.505*y**8

    b[eval(cond)] = 0 + 1.952*y + 2.908*y**2 - 3.989*y**3 - 7.985*y**4 + 11.102*y**5 + 5.491*y**6 - 10.805*y**7 + 3.347*y**8
    # ***********************************************

    return a, b

def law_uv(a, b, wave):
    """
       Returns a and b vectors of extinction law values for the UV regime
       (wavelength: 0.125um - 0.303um; wavenumber: 3.3 - 8.0).

        Parameters
        ----------
                a : ndarray
                Array of values used to correct for extinction. Only the section
                of the 'a' vector that lies in the far-UV range are altered in
                this function.

                b : ndarray
                Array of values used to correct for extinction. Only the section
                of the 'b' vector that lies in the far-UV range are altered in
                this function.

                wave : ndarray
                The wavelengths of the fluxes being corrected. Used to determine
                which 'a' and 'b' vector values to alter.

        Returns 
        ----------
                a : ndarray
                Modified 'a' vector values.

                b : ndarray
                Modified 'b' vector values.
            
    """
    cond = '(wave >=3.3) & (wave <= 8.0)'
    wave_cond = wave[eval(cond)]

    f_a = np.zeros(len(wave_cond))
    f_b = np.zeros(len(wave_cond))

    y = wave_cond

    if len(y[y >= 5.9]) > 0:
        f_a[y >= 5.9] = -0.04473*(y[y >= 5.9] - 5.9)**2 - 0.009779*(y[y >= 5.9] - 5.9)**3
        f_b[y >= 5.9] = -0.2130*(y[y >= 5.9] - 5.9)**2 - 0.1207*(y[y >= 5.9] - 5.9)**3

    a[eval(cond)] = 1.752 - 0.316 * y - 0.104/((y-4.67)**2 + 0.341) + f_a
    b[eval(cond)] = -3.090 + 1.825 * y + 1.206/((y-4.62)**2 + 0.263) + f_b

    return a, b

def law_fuv(a, b, wave):
    """
        Returns a and b vectors of extinction law values
        for the far-UV regime (wavelength: 0.1um - 0.125um; wavenumber: 8.0 -
        10.0).

        Parameters
        ----------
                a : ndarray
                Array of values used to correct for extinction. Only the section
                of the 'a' vector that lies in the far-UV range are altered in
                this function.

                b : ndarray
                Array of values used to correct for extinction. Only the section
                of the 'b' vector that lies in the far-UV range are altered in
                this function.

                wave : ndarray
                The wavelengths of the fluxes being corrected. Used to determine
                which 'a' and 'b' vector values to alter.

        Returns 
        ----------
                a : ndarray
                Modified 'a' vector values.

                b : ndarray
                Modified 'b' vector values.
    """
    
    y = wave[(wave > 8.0) & (wave <= 10)] - 8

    a[(wave > 8.0) & (wave <= 10)] = -1.073 - 0.628*y + 0.137*y**2 - 0.070*y**3
    b[(wave > 8.0) & (wave <= 10)] = 13.670 + 4.257*y - 0.420*y**2 + 0.334*y**3

    return a, b

def dered(wave, flux, e_bv, R_v = 3.089):
    """
        Corrects fluxes using the given reddening E(B-V), using the Cardelli,
        Clayton & Mathis (1989) formulation with updated IR coefficients.

        Parameters
        ----------
                wave : list, ndarray
                The wavelengths of the fluxes in microns, forces to ndarray 
                if type is list.  

                flux : list, ndarray
                The fluxes at each wavelength in spectral flux density units,
                forces to ndarray if type is list.

                e_bv : float, ufloat
                The B-V color excess for the object. If color excess is a ufloat
                then the resultant fluxes will have an associated uncertainty.

                R_v : float, ufloat
                The ratio of total to selective extinction in the UV bands B and
                V. Defaults to R_v = 3.089 if no value is given.

        Returns
        ----------
                flux : ndarray
                Returns a numpy array of extinction corrected fluxes, using the
                supplied E(B-V) and R_v values. If an uncertainty was defined
                for either the input fluxes or the E(B-V) value then the
                resultant ndarray will be an array of ufloat objects.
    """

    if (type(wave) != np.ndarray):
        wave = np.array(wave)

    if (type(flux) != np.ndarray):
        flux = np.array(flux)

    wave = 1 / wave

    a, b = np.zeros(len(wave)), np.zeros(len(wave))

    laws = [law_uv, law_fuv, law_nir, law_ir]

    for law in laws:
        a, b = law(a, b, wave)

    return flux * 10**(0.4 * ( (e_bv * R_v) * (a + b/R_v) ) )

