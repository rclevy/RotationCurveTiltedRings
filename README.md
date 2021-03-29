# RotationCurveTiltedRings

## FitTiltedRings.py
`FitTiltedRings.py` fits rotation curves from galaxy velocity fields using tilted rings and a first order harmonic decomposition. The Python implementation was developed by R. C. Levy (rlevy.astro@gmail.com) and is based on MATLAB scripts by A. D. Bolatto, J. D. Simon, and R. C. Levy. The MATLAB versions of this script were used to derive the rotation curves presented in:
 - [Simon et al. 2003](https://ui.adsabs.harvard.edu/abs/2003ApJ...596..957S/abstract)
 - [Simon et al. 2005](https://ui.adsabs.harvard.edu/abs/2005ApJ...621..757S/abstract)
 - [Levy et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...92L/abstract)


*Any use of these codes or data products must cite Cooke, Levy, et al. (2021, in prep.).*

### Input Parameter Descriptions
	gal_name : str
		Name of galaxy to fit
	header : astropy header
		astropy header of velocity field
	velfield : array
		2D numpy array containing the velocity map
	evelfield : array
		2D numpy array containing the errors on the velocity map
	RA : float
		J2000 right ascension of the kinematic center of the galaxy, in decimal degrees
	Dec : float
		J2000 declination of the kinematic center of the galaxy, in decimal degrees
	PA : float
		Position angle of the /approaching/ side of the major axis measured east of north, in degrees
	inc : float
		Inclination of the galaxy to the line of sight, in degrees
	Vsys : float
		Systemic (AKA recessional) velocity of the center of the galaxy, in km/s
	rmEndRings : int
		Number of end rings to remove from rotation curve fit, imposed after Rmax, can be set to 0 so no rings are removed
	Rmax : float
		Maximum radius of rings, if NaN this limit is not imposed, in arcsec
	save_dir : str
		Path to directory to save output plots
	plotON: flag
		If True or NOT GIVEN, will make and save intermediate plots
		If False, will NOT make intermediate plots
	rotmodel : flag
		If 'rotonly', fits only the rotation (cosine) component
		If 'full', other, or NOT GIVEN, fits full rotation model (rotation:cosine, radial:sine, systemic:constant)

### Returned Value Descriptions
	R : array
		numpy array containing the radii (center) of the fitted rings, in arcsec
	Vrot : array
		numpy array containing the fitted rotation velocity (cosine component), in km/s
	eVrot : array
		numpy array containing the uncertainty on Vrot, in km/s
	Vrad : array
		numpy array containting the fitted radial velocity (sine component), in km/s
	eVrad : array
		numpy array containing the uncertainty on Vrad, in km/s
	dVsys : array
		numpy array containing the fitted deviation from the systemic velocity (constant component), in km/s
	edVsys : array
		numpy array containing the uncertainty on Vsys, in km/s


## FitCORotationCurves.py
`FitCORotationCurves.py` is a convenience wrapper function to load the relevant data, plot the velocity fields, fit, save, and plot the rotation curves, and do a Monte Carlo to get uncertainties on the rotation curves. This wrapper was written by R. C. Levy (rlevy.astro@gmail.com).

*Any use of these codes or data products must cite Cooke, Levy, et al. (2021, in prep.).*

## CODataProducts_Cooke21.tar.gz
The ALMA data used for this study are part of program ADS/JAO.ALMA#2015.1.00820.S (PI L. Blitz) and are publicly available on the [ALMA Science Archive](https://almascience.nrao.edu/asax/). Data products published in Cooke, Levy et al. (2021, in prep.) are available here (`CODataProducts_Cooke2021.tar.gz`). See Section 2 of Cooke, Levy, et al. (2021, in prep.) for details on the data calibration and imaging. The products presented here are the results of Gaussian fits to the data cube to extract the peak intensity, velocity centroid, and line width (see [Levy et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...92L/abstract) for more details). These data products include:
 - Peak intensity: \*.gausspk.fits
 - Error in peak intensity: \*.egausspk.fits
 - Velocity centroid: \*.cgrad.vcen.fits
 - Error in velocity centroid: \*.evcen.fits
 - Linewidth (FWHM): \*.fwhm.fits
 - Error in linewidth: \*.efwhm.fits   

*Any use of these codes or data products must cite Cooke, Levy, et al. (2021, in prep.).*

## Notes
This repository's main branch is called [main, not master](https://www.cnet.com/news/microsofts-github-is-removing-coding-terms-like-master-and-slave/).
