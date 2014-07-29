*********************************************** 
This is a function designed to take spectral
flux densities and a B-V color excess, E(B-V), 
and return extinction corrected fluxes using
the interstellar extinction law of Cardelli,
Clayton & Mathis (1989).
*********************************************** 

To use the function it can be imported via 

    import deredden

    or

    import deredden as ...

and then simply call the primary function via

    deredden.dered(wave, flux, e_bv)

where the arguments are lists/ndarray's of wavelengths 
and fluxes and E(B-V) and R_v (optional) floats or 
ufloats*.

If no R_v value is given then it will default to 

    R_v = 3.089

the average R_v for the diffuse interstellar medium.

*ufloats are objects from the uncertainties python
package.
