# RotationCurveTiltedRings
Script to fit rotation curves from galaxy velocity fields using tilted rings and a first order harmonic decomposition. Please see examples for use in the code header.

The Python implementation was developed by R. C. Levy (rlevy.astro@gmail.com) and is based on MATLAB scripts by A. D. Bolatto, J. D. Simon, and R. C. Levy.

The Python implementation will be published in a forthcoming paper (Cooke, Levy, Bolatto, et al. 2021, in prep.).

The MATLAB versions of this script were used to derive the rotation curves presented in:
- [Simon et al. 2003](https://ui.adsabs.harvard.edu/abs/2003ApJ...596..957S/abstract)
- [Simon et al. 2005](https://ui.adsabs.harvard.edu/abs/2005ApJ...621..757S/abstract)
- [Levy et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...92L/abstract)

*Any use of this code must cite Cooke, Levy, et al. (2021, in prep) and should cite the above three papers as well.*

### CO data products presented in Cooke, Levy, et al. (2021, in prep.) (CoDataProducts_Cooke2021.tar.gz)

The ALMA data used for this study are part of program ADS/JAO.ALMA#2015.1.00820.S (PI L. Blitz) and are publicly available on the [ALMA Science Archive](https://almascience.nrao.edu/asax/).

See Section 2 of Cooke, Levy, et al. (2021, in prep.) for details on the data calibration and imaging. The products presented here are the results of Gaussian fits to the data cube to extract the peak intensity, velocity centroid, and line width (see [Levy et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...92L/abstract) for more details).

*Any use of these data products must cite Cooke, Levy, et al. (2021, in prep.).*


#### Data Products
- Peak intensity: \*.gausspk.fits
- Error in peak intensity: \*.egausspk.fits
- Velocity centroid: \*.cgrad.vcen.fits
- Error in velocity centroid: \*.evcen.fits
- Linewidth (FWHM): \*.fwhm.fits
- Error in linewidth: \*.efwhm.fits

## Notes
This repository's main branch is called [main, not master](https://www.cnet.com/news/microsofts-github-is-removing-coding-terms-like-master-and-slave/).
