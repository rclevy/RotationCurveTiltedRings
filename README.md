# RotationCurveTiltedRings
`FitTiltedRings.py` fits rotation curves from galaxy velocity fields using tilted rings and a first order harmonic decomposition. The Python implementation was developed by R. C. Levy (rlevy.astro@gmail.com) and is based on MATLAB scripts by A. D. Bolatto, J. D. Simon, and R. C. Levy. The MATLAB versions of this script were used to derive the rotation curves presented in:
 - [Simon et al. 2003](https://ui.adsabs.harvard.edu/abs/2003ApJ...596..957S/abstract)
 - [Simon et al. 2005](https://ui.adsabs.harvard.edu/abs/2005ApJ...621..757S/abstract)
 - [Levy et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...92L/abstract)   
 - *Any use of this code must cite Cooke, Levy, et al. (2021, in prep) and should cite the above three papers as well.*

`FitCORotationCurves.py` is a convenience wrapper function to load the relevant data, plot the velocity fields, fit, save, and plot the rotation curves, and do a Monte Carlo to get uncertainties on the rotation curves. This wrapper was written by R. C. Levy (rlevy.astro@gmail.com).

The ALMA data used for this study are part of program ADS/JAO.ALMA#2015.1.00820.S (PI L. Blitz) and are publicly available on the [ALMA Science Archive](https://almascience.nrao.edu/asax/). Data products published in Cooke, Levy et al. (2021, in prep.) are available here (`CODataProducts_Cooke2021.tar.gz`). See Section 2 of Cooke, Levy, et al. (2021, in prep.) for details on the data calibration and imaging. The products presented here are the results of Gaussian fits to the data cube to extract the peak intensity, velocity centroid, and line width (see [Levy et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...92L/abstract) for more details). These data products include:
 - Peak intensity: \*.gausspk.fits
 - Error in peak intensity: \*.egausspk.fits
 - Velocity centroid: \*.cgrad.vcen.fits
 - Error in velocity centroid: \*.evcen.fits
 - Linewidth (FWHM): \*.fwhm.fits
 - Error in linewidth: \*.efwhm.fits   
 - *Any use of these data products must cite Cooke, Levy, et al. (2021, in prep.).*

## Notes
This repository's main branch is called [main, not master](https://www.cnet.com/news/microsofts-github-is-removing-coding-terms-like-master-and-slave/).
