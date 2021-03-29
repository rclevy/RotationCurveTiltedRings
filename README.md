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
		The RA is the same for all of the rings
	Dec : float
		J2000 declination of the kinematic center of the galaxy, in decimal degrees
		The Dec is the same for all of the rings
	PA : float
		Position angle of the /approaching/ side of the major axis measured east of north, in degrees
		The PA is the same for all of the rings
	inc : float
		Inclination of the galaxy to the line of sight, wher 0 deg is face-on and 90 deg is edge-on, in degrees
		The inc is the same for all of the rings
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
	eR : array
		numpy array containing the uncertainty on the radii given by the width of the rings, in arcsec
	Vrot : array
		numpy array containing the fitted rotation velocity (cosine component), in km/s
	eVrot : array
		numpy array containing the uncertainty on Vrot, in km/s
		This is the statistical fitting uncertaities, which are usually small and underestimate the true systematic uncertainties
	Vrad : array
		numpy array containting the fitted radial velocity (sine component), in km/s
	eVrad : array
		numpy array containing the uncertainty on Vrad, in km/s
		This is the statistical fitting uncertaities, which are usually small and underestimate the true systematic uncertainties
	dVsys : array
		numpy array containing the fitted deviation from the systemic velocity (constant component), in km/s
	edVsys : array
		numpy array containing the uncertainty on Vsys, in km/s
		This is the statistical fitting uncertaities, which are usually small and underestimate the true systematic uncertainties
	chisq :  array
		chi^2 statistic of fit per ring
	chisqr : array
		reduced chi^2 statistic of fit per ring
	rms :  array
		root-mean-square error of fit per ring

### Example
	# import modules
	import numpy as np
	from astropy.io import fits
	import pandas as pd
	from FitTiltedRings import fit_tilted_rings
	
	gal_name = 'NGC6106'
	
	# load galaxy parameters from file
	fname_gal_params = '../Data/Galaxy_Parameters.csv'
	gal_params = pd.read_csv(fname_gal_params) 
	idx = gal_params['Name'].values.tolist().index(gal_name)
	RA = gal_params['RA'][idx] #deg
	Dec = gal_params['Dec'][idx] #deg
	Inc = gal_params['Inc'][idx] #deg
	PA = gal_params['PA'][idx] #deg
	Vsys = float(gal_params['Vsys'][idx]) #km/s
	Rmax = np.nan #don't impose max ring radius
	rmEndRings = 0 #don't impose max ring radius
	rotmodel = 'full' #fit full rotation model (cosine+sine+constant)
	plotOn = True #save intermediate plots
	save_dir = '../Plots/ringfit/' 
	
	# open the velocity field data
	vel_fits = fits.open('../Data/gaussfit/'+gal_name+'.cgrad.vcen.fits',ignore_missing_end=True)
	evel_fits = fits.open('../Data/gaussfit/'+gal_name+'.evcen.fits',ignore_missing_end=True)
	hdr = vel_fits[0].header
	vel = vel_fits[0].data
	vel[vel==-999] = np.nan
	evel = evel_fits[0].data
	evel[evel==-999] = np.nan
	
	# run the ring fitting
	R,eR,Vrot,eVrot,Vrad,eVrad,dVsys,edVsys,chisq,chisqr,rms=fit_tilted_rings(gal_name,hdr,vel,evel,RA,Dec,PA,Inc,Vsys,rmEndRings,Rmax,save_dir,plotOn,rotmodel)
	
Example output plot of the fit in one ring: ![NGC6106_ringfit_21 292](https://user-images.githubusercontent.com/14076216/112869601-596b7300-908b-11eb-9458-7ecd29fb992c.png)
 - This is an example of a "good" fit. Things to look that indicate a good fit are:
 	- The data (blue points) track the full model (green) curve well. Data which are "flat topped" or "peaky" with respect to the curves usually indicates that the inclination is wrong. Data points that only sample one side of the galaxy (i.e., only positive or negative azimuthal angles) may produce unreliable fits.
 	- The full model (green) and rotation-only (red) curves overlap. A horizontal offset between these may mean that the PA is wrong or that there is a kinematic twist in the galaxy.
 	- The peak of the model curves is at an azimuthal angle of 0 deg. If instead there's a trough at 0 deg, the PA is wrong by 180 deg.

Summary plot of fits in all rings: ![NGC6106_ringfit](https://user-images.githubusercontent.com/14076216/112869631-61c3ae00-908b-11eb-9e74-c62dac2ca076.png)


Plot of the rotation curve components: ![NGC6106_COrotcurve](https://user-images.githubusercontent.com/14076216/112869682-66886200-908b-11eb-83ce-b925662cb9b6.png)
 - The uncertainties here reflect only the statistical fitting uncertainties. Other methods (bootstrapping, Monte Carlo) are useful to determine more representative uncertainties that account for systematic uncertainties. `FitCORotationCurves.py` (see below) gives an exmaple of doing a Monte Carlo over uncertainties in the center position, PA, and inc so that the errors reflect systematic uncertaintites in these parameters.


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
