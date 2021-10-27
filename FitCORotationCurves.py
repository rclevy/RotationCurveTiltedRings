this_script = __file__

#wrapper to fit all CO rotation curves
import argparse
parser=argparse.ArgumentParser(
    description='''''',
    epilog='''''')
parser.add_argument('gal_name',type=str, help='Name of galaxy to fit')
parser.add_argument('--make_grid',type=str,help='Make the grid of models (True) or not (False). Default is True')
parser.add_argument('--rot_only',type=str,help='Fit only the rotation component (True) or rotation, radial, and systemic (False). Default is False.')
parser.add_argument('--plotOn',type=str,help='Save the individual ring fits (True) or not (False). Default is True.')

args=parser.parse_args()
gal_name = args.gal_name

def main(gal_name):
	r'''
	Main function that load data and calls functions for a single galaxy.

	Parameters
	----------
	gal_name : str
		name of galaxy to fit
	make_grid : str
		optional, if 'True' make a grid of models varying the geometric parameters, default is True, this is a slow process!
		
	Returns
	-------
	
	Notes
	-----
	Required modules: astropy, matplotlib, numpy, pandas, os (called functions have their own required modules)
	Based on WW2020_Data_Tutorial.ipynb
	Author: R. C. Levy (rlevy.astro@gmail.com)
	Last updated: 2021-10-27
	Change log:
		2020-05-28 : file created, RCL
		2021-03-09 : added support to fit rotation-only model, RCL
		2021-10-27 : better handling of uncertainties, make a grid of models varying the input parameters, RCL
	'''
	# Import some base modules we'll need throughout
	import numpy as np
	import pandas as pd
	from astropy.io import fits
	from astropy.wcs import WCS
	from astropy.visualization import (AsinhStretch,ImageNormalize,MinMaxInterval)
	import matplotlib.pyplot as plt
	import matplotlib.ticker as mticker
	from matplotlib.patches import Ellipse
	import os
	import scipy.ndimage as ndimage
	from FitTiltedRings import fit_tilted_rings
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['font.size'] = 14
	plt.rcParams['mathtext.rm'] = 'serif'
	plt.rcParams['mathtext.fontset'] = 'cm'

	import warnings
	from astropy.io.fits.verify import VerifyWarning
	warnings.simplefilter('ignore', category=VerifyWarning) #don't throw warnings about the FITS headers


	def arcsec2kpc(x):
		return x/206265*dist*1E3
	def kpc2arcsec(x):
		return x/(1E3*dist)*206265

	if args.make_grid:
		if args.make_grid=='True':
			make_grid=True
		else:
			make_grid=False
	else:
		make_grid=True

	if args.rot_only:
		if args.rot_only=='True':
			rot_only=True
			rotmodel='rotonly'
		else:
			rot_only=False
			rotmodel='full'
	else:
		rot_only=False
		rotmodel='full'

	if args.plotOn:
		if args.plotOn=='True':
			plotOn=True
		else:
			plotOn=False
	else:
		plotOn=True

	#load the geometric parameters
	fname_gal_params = '../Data/Galaxy_Parameters.csv'
	gal_params = pd.read_csv(fname_gal_params) 
	idx = gal_params['Name'].values.tolist().index(gal_name)
	RA = gal_params['RA'][idx] #deg
	Dec = gal_params['Dec'][idx] #deg
	dist = gal_params['Distance'][idx] #Mpc 
	Inc = gal_params['Inc'][idx] #deg
	PA = gal_params['PA'][idx] #deg
	Vsys = float(gal_params['Vsys'][idx]) #km/s, relativistic
	rmEndRings = gal_params['rmEndRings'][idx]

	#load the CO velocity data
	vel_fits = fits.open('../Data/gaussfit/'+gal_name+'.cgrad.vcen.fits',ignore_missing_end=True)
	evel_fits = fits.open('../Data/gaussfit/'+gal_name+'.evcen.fits',ignore_missing_end=True)
	hdr = vel_fits[0].header
	wcs = WCS(hdr)
	vel = vel_fits[0].data
	evel = evel_fits[0].data
	vel[vel == hdr['BLANK']]=np.nan
	evel[evel == hdr['BLANK']]=np.nan
	bmaj = hdr['BMAJ']/hdr['CDELT2'] #pixels
	bmin = hdr['BMIN']/hdr['CDELT2'] #pixels
	bpa = hdr['BPA']

	#convert velocity from radio to relativistic velocity frame
	c = 2.9979E5 #km/s
	vel = c*((1-(1-vel/c)**2)/(1+(1-vel/c)**2))
	evel = c*((1-(1-evel/c)**2)/(1+(1-evel/c)**2))


	RA_pix,Dec_pix = wcs.all_world2pix(RA,Dec,1,ra_dec_order=True)

	#mask the CO velocity field
	#open the gausspeak data for masking
	peak = fits.open('../Data/gaussfit/'+gal_name+'.gausspk.fits',ignore_missing_end=True)
	epeak = fits.open('../Data/gaussfit/'+gal_name+'.egausspk.fits',ignore_missing_end=True)
	data_peak = peak[0].data
	edata_peak = epeak[0].data
	#get snr
	snr = data_peak/edata_peak
	#remove NaNs
	data_peak[data_peak==-999]=np.nan
	edata_peak[edata_peak==-999]=np.nan
	#smooth the snr image
	snr[np.isnan(data_peak)==True]=0.
	snr = ndimage.gaussian_filter(snr,sigma=(bmaj/2.355/2,bmaj/2.355/2),order=0)
	vel[snr<3] = np.nan

	#plot the CO velocity field
	#crop to 1' square
	crop = 1.5/60/hdr['CDELT2']
	xcen = RA_pix
	ycen = Dec_pix
	xtext = (xcen-crop/2)+0.05*((xcen+crop/2)-(xcen-crop/2))
	ytext = (ycen-crop/2)+0.05*((ycen+crop/2)-(ycen-crop/2))
	xsb = (xcen-crop/2)+0.9*((xcen+crop/2)-(xcen-crop/2))
	ysb = (ycen-crop/2)+0.05*((ycen+crop/2)-(ycen-crop/2))
	#add a scale bar
	sb_kpc = 0.5 #kpc
	sb_arcsec = sb_kpc/(dist*1E3)*206265
	sb_pix = sb_arcsec/3600/hdr['CDELT2']
	sb_str = str(sb_kpc)+' kpc'


	#run the ringfitting
	save_dir='../Plots/ringfit/'
	print('Fitting with '+rotmodel+' model')
	R,eR,Vrot,eVrot,Vrad,eVrad,dVsys,edVsys,chisq,chisqr,rms=fit_tilted_rings(gal_name,hdr,vel,evel,RA,Dec,PA,Inc,Vsys,rmEndRings,np.nan,save_dir,plotOn,rotmodel)
	eVrot_fit = rms.copy() #error on Vrot from the rms of the fit in each ring
	# eVrot_stat = eVrot.copy() #statistical fitting error on Vrot, very small
	# eVrad_stat = eVrad.copy() #statistical fitting error on Vrad, very small
	# edVsys_stat = edVsys.copy() #statistical fitting error on dVsys, very small

	#plot the velocity field
	plt.figure(1)
	plt.clf()
	ax=plt.subplot(projection=wcs)
	im=ax.imshow(vel-Vsys,origin='lower',cmap='RdYlBu_r',vmin=-100,vmax=100)
	ax.plot(RA_pix,Dec_pix,'kx')
	cb=plt.colorbar(im)
	cb.set_label('V$_{\mathrm{CO}}$-V$_{\mathrm{sys}}$ (km s$^{-1}$)')
	ax.coords[0].set_major_formatter('hh:mm:ss.s')
	ax.coords[0].set_separator(('$^{\mathrm{h}}$','$^{\mathrm{m}}$','$^{\mathrm{s}}$'))
	ax.coords[0].display_minor_ticks(True)
	ax.coords[1].display_minor_ticks(True)
	ax.coords[0].set_minor_frequency(4)
	ax.coords[1].set_minor_frequency(4)
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.coords[1].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Decl. (J2000)')
	plt.text(0.03,0.97,gal_name,horizontalalignment='left',verticalalignment='top',transform=ax.transAxes,fontsize=plt.rcParams['font.size']+2)
	ax.contour(snr,levels=np.arange(3,4,1),colors='gray',linewidths=1.0) 
	ax.add_patch(Ellipse((xtext,ytext),bmaj,bmin,bpa+90, ec='k',fc='k'))
	ax.plot([xsb-sb_pix,xsb],[ysb,ysb],'k',lw=1.5)
	ax.text(np.mean([xsb-sb_pix,xsb]),ysb+2,sb_str,ha='center',va='bottom',fontsize=plt.rcParams['font.size']-2)
	ax.set_xlim(left=xcen-crop/2,right=xcen+crop/2)
	ax.set_ylim(bottom=ycen-crop/2,top=ycen+crop/2)
	plt.savefig('../Plots/Vfields/'+gal_name+'_vfield.pdf',bbox_inches='tight',metadata={'Creator': this_script})

	#add contours
	cstep = 20.
	levels = np.arange(-100,100+cstep,cstep)
	ax.contour(vel-Vsys,levels=levels,colors='gray',linewidths=0.5)
	plt.savefig('../Plots/Vfields/'+gal_name+'_vfield_contours.pdf',bbox_inches='tight',metadata={'Creator': this_script})

	#plot the rings over the velcity field
	plt.figure(1)
	R_pix = R[1:]/3600/hdr['CDELT2'] #ring center
	R_shift = np.concatenate([R_pix[1:],np.array([np.inf])])
	R_outer = np.nanmean([R_pix,R_shift],axis=0)[0:-1]
	[ax.add_patch(Ellipse((RA_pix,Dec_pix),2*ro,2*ro*np.cos(np.radians(Inc)),PA+90, ec='k',fc='None',zorder=5)) for ro in R_outer]
	plt.savefig('../Plots/Vfields/'+gal_name+'_vfield_rings.pdf',bbox_inches='tight',metadata={'Creator': this_script})


	#do a MonteCarlo over the geometric parameters
	print('Doing a Monte Carlo to get RC uncertainties...')
	ntrials = 50
	eevrot = np.zeros((ntrials,int(len(Vrot))))
	eevrad = np.zeros((ntrials,int(len(Vrot))))
	eedvsys = np.zeros((ntrials,int(len(Vrot))))
	eevrot_ro = np.zeros((ntrials,int(len(Vrot))))
	pa_off = 5. #deg
	inc_off = 4. #deg
	cen_off = 0.4/3600 #deg
	for j in range(ntrials):
		ra = RA + np.random.uniform(low=-cen_off,high=cen_off)
		dec = Dec + np.random.uniform(low=-cen_off,high=cen_off)
		pa = PA + np.random.uniform(low=-pa_off,high=pa_off)
		inc = Inc + np.random.uniform(low=-inc_off,high=inc_off)
		_,_,e1,_,e2,_,e0,_,_,_,_=fit_tilted_rings(gal_name,hdr,vel,evel,ra,dec,pa,inc,Vsys,rmEndRings,np.nan,save_dir,False,rotmodel)

		if (np.all(np.isnan(e1))==True) | (np.all(np.isnan(e2))==True) | (np.all(np.isnan(e0))==True):
			#throw away trials that result in NaNs
			eevrot[j,:]=np.nan
			eevrad[j,:]=np.nan
			eedvsys[j,:]=np.nan
		else:
			if len(e1) > len(Vrot):
				#if MC has more rings, trim extra rings
				e1=e1[0:len(Vrot)]
				e2=e2[0:len(Vrot)]
				e0=e0[0:len(Vrot)]
			elif len(e1) < len(Vrot):
				#if MC has fewer rings, pad with NaNs
				pad = np.zeros((int(len(Vrot)-len(e1)),))*np.nan
				e1 = np.concatenate((e1,pad))
				e2 = np.concatenate((e2,pad))
				e0 = np.concatenate((e0,pad))				
			
			eevrot[j,:]=e1
			eevrad[j,:]=e2
			eedvsys[j,:]=e0

	eVrot_MC = np.nanstd(eevrot,axis=0)[0:len(R)]
	eVrad_MC = np.nanstd(eevrad,axis=0)[0:len(R)]
	edVsys_MC = np.nanstd(eedvsys,axis=0)[0:len(R)]

	# #make a grid of RC models varying PA and Inc for the decompositon
	if make_grid==True:
		print('Making a grid of RCs in PA and Inc')
		grid_step = 1.
		PA_grid = np.arange(PA-pa_off,PA+pa_off+grid_step,grid_step)
		Inc_grid = np.arange(Inc-inc_off,Inc+inc_off+grid_step,grid_step)
		for j in range(len(PA_grid)):
			for k in range(len(Inc_grid)):
				r_g,_,vrot_g,_,_,_,_,_,_,_,evrot_g=fit_tilted_rings(gal_name,hdr,vel,evel,RA,Dec,PA_grid[j],Inc_grid[k],Vsys,rmEndRings,np.nan,save_dir,False,rotmodel)

				if len(vrot_g) > len(Vrot):
					#if MC has more rings, trim extra rings
					r_g = r_g[0:len(Vrot)]
					vrot_g=vrot_g[0:len(Vrot)]
					evrot_g=evrot_g[0:len(Vrot)]
				elif len(vrot_g) < len(Vrot):
					#if MC has fewer rings, pad with NaNs
					pad = np.zeros((int(len(Vrot)-len(vrot_g)),))*np.nan
					r_g = np.concatenate((r_g,pad))
					vrot_g = np.concatenate((vrot_g,pad))
					evrot_g = np.concatenate((evrot_g,pad))

				#save the grids to a file to be passed in the the decomposition later on
				if rotmodel == 'full':
					fname = '../Data/CORotationCurves/ParamGrid/'+gal_name+'_CORotCurve_PA'+str(int(PA_grid[j]))+'_Inc'+str(int(Inc_grid[k]))+'.txt'
				else:
					fname = '../Data/CORotationCurves/ParamGrid/'+gal_name+'_CORotCurve_'+rotmodel+'_PA'+str(int(PA_grid[j]))+'_Inc'+str(int(Inc_grid[k]))+'.txt'
			
				with open(fname,'w') as fp:
					fp.write('#galaxy = %s\n' %gal_name)
					fp.write('#PA = %.1f\n' %PA_grid[j])
					fp.write('#Inc = %.1f\n' %Inc_grid[k])
					fp.write('#Vsys = %.1f\n' %Vsys)
					fp.write('#RA = %f\n' %ra)
					fp.write('#Dec = %f\n' %dec)
					fp.write('#column [0] = radius (arcsec),[1] = radius (kpc), [2] = rotation velocity (km/s), [3] = uncertainty on rotation velocity from rms of fit in each ring (km/s)\n')
					for i in range(len(R)):
						fp.write('%.2f\t%.2f\t%.2f\t%.2f\n'
							%(r_g[i],arcsec2kpc(r_g[i]),vrot_g[i],evrot_g[i]))

			
	#remove the 0,0 point
	R=R[1:]
	eR=eR[1:]
	Vrot=Vrot[1:]
	eVrot_MC=eVrot[1:]
	eVrot_fit=eVrot_fit[1:]
	eevrot = eevrot[:,1:]
	Vrad=Vrad[1:]
	eVrad_MC=eVrad[1:]
	eevrad = eevrad[:,1:]
	dVsys=dVsys[1:]
	edVsys_MC=edVsys[1:]
	eedvsys = eedvsys[:,1:]
	chisq=chisq[1:]
	chisqr=chisqr[1:]
	rms=rms[1:]

	#remove all NaN points
	R=R[np.isnan(rms)==False]
	eR=eR[np.isnan(rms)==False]
	Vrot=Vrot[np.isnan(rms)==False]
	eVrot_MC=eVrot_MC[np.isnan(rms)==False]
	eVrot_fit=eVrot_fit[np.isnan(rms)==False]
	eevrot = eevrot[:,np.isnan(rms)==False]
	Vrad=Vrad[np.isnan(rms)==False]
	eVrad_MC=eVrad_MC[np.isnan(rms)==False]
	eevrad = eevrad[:,np.isnan(rms)==False]
	dVsys=dVsys[np.isnan(rms)==False]
	edVsys_MC=edVsys_MC[np.isnan(rms)==False]
	eedvsys = eedvsys[:,np.isnan(rms)==False]
	chisq=chisq[np.isnan(rms)==False]
	chisqr=chisqr[np.isnan(rms)==False]
	rms=rms[np.isnan(rms)==False]

	eVrot_MC_u = eVrot_MC.copy()
	eVrot_MC_l = eVrot_MC.copy()

	#make a second x-axis on the top in kpc
	R_kpc = arcsec2kpc(R) #kpc
	beam = hdr['BMAJ']*3600 #arcsec
	beam_kpc = arcsec2kpc(beam) #kpc

	#open the Ha rotation curve if it exists
	harc_file = '../Data/HalphaRotationCurves/'+gal_name+'.VrotHa.csv'
	if os.path.exists(harc_file)==True:
		# dist_file = pd.read_csv('../Data/Beta_starComparison/CommonGalaxies.csv')
		# idx = dist_file['Name'].values.tolist().index(gal_name)
		dist_R19 = gal_params['Distance_R19a'].values[idx] #Mpc
		vsys_R19_opt = gal_params['Vsys_R19a'].values[idx] #km/s, optical
		harc = pd.read_csv(harc_file)
		Rha = harc['Rkpc'].values/(1E3*dist_R19)*206265 #arcsec
		Vha = harc['Vrot'].values
		eVha = harc['eVrot'].values
		Rha_kpc = arcsec2kpc(Rha) #kpc

		#convert Ha rot curve to from optical to rel convention
		#see eqs 19 & 21 of Levy+2018
		vsys_r19 = c*(((vsys_R19_opt/c+1)**2-1)/((vsys_R19_opt/c+1)**2+1)) #km/s, rel
		Vha = Vha*(1-vsys_r19/c)



	#plot the rotation curve
	plt.figure(1)
	plt.clf()
	plt.axhline(0,color='gray',lw=0.75)
	plt.axvspan(xmin=0,xmax=beam_kpc,fc='lightgray',ec='gray',alpha=0.25,hatch='\\\\',label='CO Beam FWHM',zorder=2)
	plt.axvspan(0.3,0.8,ec='None',fc='k',alpha=0.1,zorder=1)
	for ii in range(eevrot.shape[0]):
		inx = np.where(np.isnan(eevrot[ii,:])==False)
		el1, = plt.plot(R_kpc[inx],np.squeeze(eevrot[ii,inx]),'b-',alpha=0.05,zorder=2,lw=1)
	if rotmodel=='full':
		for ii in range(eevrot.shape[0]):
			inx = np.where(np.isnan(eevrad[ii,:])==False)
			el2, = plt.plot(R_kpc[inx],np.squeeze(eevrad[ii,inx]),'r-',alpha=0.05,zorder=2,lw=1)
			inx = np.where(np.isnan(eedvsys[ii,:])==False)
			el3, = plt.plot(R_kpc[inx],np.squeeze(eedvsys[ii,inx]),'g-',alpha=0.05,zorder=2,lw=1)
	l1=plt.errorbar(R_kpc,Vrot,yerr=eVrot_fit,color='b',marker='o',linestyle='-',lw=1.5,markersize=4,capsize=3,label=r'CO V$_{\mathrm{rot}}$',zorder=4)
	if rotmodel=='full':
		l2,=plt.plot(R_kpc,Vrad,'r-',lw=1.5,label=r'CO V$_{\mathrm{rad}}$',zorder=4)
		l3,=plt.plot(R_kpc,dVsys,'g-',lw=1.5,label=r'CO $\Delta$V$_{\mathrm{sys}}$',zorder=4)
	xlim = plt.gca().get_xlim()
	if os.path.exists(harc_file)==True:
		l4=plt.errorbar(Rha_kpc,Vha,yerr=eVha,fmt='^-',color='rebeccapurple',markersize=4,capsize=3,lw=1.5,label=r'H$\alpha$ V$_{\mathrm{rot}}$')
	plt.minorticks_on()
	plt.xlabel('Radius (kpc)')
	plt.ylabel('Velocity (km s$^{-1}$)')
	plt.xlim(left=0,right=xlim[1])

	if os.path.exists(harc_file)==True:
		if rotmodel == 'full':
			plt.legend([(el1,l1),(el2,l2),(el3,l3),(l4)],[l1.get_label(),l2.get_label(),l3.get_label(),l4.get_label()],
				loc='lower right',fontsize=plt.rcParams['font.size']-2)
		else:
			plt.legend([(el1,l1),(l4)],[l1.get_label(),l4.get_label()],
				loc='lower right',fontsize=plt.rcParams['font.size']-2)
	else:
		if rotmodel=='full':
			plt.legend([(el1,l1),(el2,l2),(el3,l3)],[l1.get_label(),l2.get_label(),l3.get_label()],
				loc='lower right',fontsize=plt.rcParams['font.size']-2)
		else:
			plt.legend([(el1,l1)],[l1.get_label()],
				loc='lower right',fontsize=plt.rcParams['font.size']-2)

	plt.text(0.015,0.982,gal_name.replace('NGC','NGC '),transform=plt.gca().transAxes,fontsize=plt.rcParams['font.size']+2,va='top',ha='left')

	ax2 = plt.gca().secondary_xaxis('top',functions=(kpc2arcsec,arcsec2kpc))
	ax2.set_xlabel('Radius (arcsec)')
	ax2.minorticks_on()

	if rotmodel == 'full':
		plt.savefig('../Plots/CORotationCurves/'+gal_name+'_COrotcurve.pdf',bbox_inches='tight',metadata={'Creator': this_script})
	else:
		plt.savefig('../Plots/CORotationCurves/'+gal_name+'_COrotcurve_'+rotmodel+'.pdf',bbox_inches='tight',metadata={'Creator': this_script})

	plt.close('all')


	#write to file
	if rotmodel == 'full':
		fname = '../Data/CORotationCurves/'+gal_name+'_CORotCurve.txt'
	else:
		fname = '../Data/CORotationCurves/'+gal_name+'_CORotCurve_'+rotmodel+'.txt'
	
	with open(fname,'w') as fp:
		fp.write('#galaxy = %s\n' %gal_name)
		fp.write('#PA = %.1f\n' %PA)
		fp.write('#Inc = %.1f\n' %Inc)
		fp.write('#Vsys = %.1f\n' %Vsys)
		fp.write('#RA = %f\n' %RA)
		fp.write('#Dec = %f\n' %Dec)
		fp.write('#column [0] = radius (arcsec),[1] = uncertainty on radius, [2] = rotation velocity (km/s), [3] = upper uncertainty on rotation velocity from Monte Carlo, [4] = lower uncertainty on rotation velocity from Monte Carlo, [5] = radial velocity (km/s), [6] = uncertainty on radial velocity from Monte Carlo, [7] = devitation from systemic velocity (km/s), [8] = uncertainty on deviation from systemic velocity from Monte Carlo, [9] = reduced chi squared of fit, [10] = rms of fit\n')
		for i in range(len(R)):
			fp.write('%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'
				%(R[i],eR[i],Vrot[i],eVrot_MC_u[i],eVrot_MC_l[i],Vrad[i],eVrad_MC[i],dVsys[i],edVsys_MC[i],chisqr[i],rms[i]))


	#make a model velocity field
	import mkvelfield
	v_model = mkvelfield.mk_model_velfield(hdr,PA,Inc,RA,Dec,0.,0.,R,Vrot,Vrad,dVsys)


	#plot the model velocity field
	plt.figure(1)
	plt.clf()
	ax=plt.subplot(projection=wcs)
	im=ax.imshow(v_model,origin='lower',cmap='RdYlBu_r',vmin=-100,vmax=100)
	ax.plot(RA_pix,Dec_pix,'kx')
	cb=plt.colorbar(im)
	cb.set_label('V$_{\mathrm{model}}$ (km s$^{-1}$)')
	ax.coords[0].set_major_formatter('hh:mm:ss.s')
	ax.coords[0].set_separator(('$^{\mathrm{h}}$','$^{\mathrm{m}}$','$^{\mathrm{s}}$'))
	ax.coords[0].display_minor_ticks(True)
	ax.coords[1].display_minor_ticks(True)
	ax.coords[0].set_minor_frequency(4)
	ax.coords[1].set_minor_frequency(4)
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.coords[1].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Decl. (J2000)')
	plt.text(0.03,0.97,gal_name,horizontalalignment='left',verticalalignment='top',transform=ax.transAxes,fontsize=plt.rcParams['font.size']+2)
	ax.contour(snr,levels=np.arange(3,4,1),colors='gray',linewidths=1.0) 
	ax.add_patch(Ellipse((xtext,ytext),bmaj,bmin,bpa+90, ec='k',fc='k'))
	ax.plot([xsb-sb_pix,xsb],[ysb,ysb],'k',lw=1.5)
	ax.text(np.mean([xsb-sb_pix,xsb]),ysb+2,sb_str,ha='center',va='bottom',fontsize=plt.rcParams['font.size']-2)
	ax.set_xlim(left=xcen-crop/2,right=xcen+crop/2)
	ax.set_ylim(bottom=ycen-crop/2,top=ycen+crop/2)
	if rotmodel=='full':
		plt.savefig('../Plots/Vfields/'+gal_name+'_vfield_model.pdf',bbox_inches='tight',metadata={'Creator': this_script})
	else:
		plt.savefig('../Plots/Vfields/'+gal_name+'_vfield_model_'+rotmodel+'.pdf',bbox_inches='tight',metadata={'Creator': this_script})

	#find residual velocity field
	v_resid = (vel-Vsys)-v_model
	vmax = np.nanmedian(np.abs(v_resid))+3*np.nanstd(np.abs(v_resid))
	if vmax > np.nanmax(np.abs(v_resid)):
		vmax = np.nanmax(np.abs(v_resid))
	if vmax > 100.:
		vmax=100.

	#plot the residual velocity
	plt.figure(1)
	plt.clf()
	ax=plt.subplot(projection=wcs)
	im=ax.imshow(v_resid,origin='lower',cmap='RdYlBu_r',vmin=-vmax,vmax=vmax)
	ax.plot(RA_pix,Dec_pix,'kx')
	cb=plt.colorbar(im)
	cb.set_label('V$_{\mathrm{CO}}$-V$_{\mathrm{model}}$ (km s$^{-1}$)')
	ax.coords[0].set_major_formatter('hh:mm:ss.s')
	ax.coords[0].set_separator(('$^{\mathrm{h}}$','$^{\mathrm{m}}$','$^{\mathrm{s}}$'))
	ax.coords[0].display_minor_ticks(True)
	ax.coords[1].display_minor_ticks(True)
	ax.coords[0].set_minor_frequency(4)
	ax.coords[1].set_minor_frequency(4)
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.coords[1].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Decl. (J2000)')
	plt.text(0.03,0.97,gal_name,horizontalalignment='left',verticalalignment='top',transform=ax.transAxes,fontsize=plt.rcParams['font.size']+2)
	ax.contour(snr,levels=np.arange(3,4,1),colors='gray',linewidths=1.0) 
	ax.add_patch(Ellipse((xtext,ytext),bmaj,bmin,bpa+90, ec='k',fc='k'))
	ax.plot([xsb-sb_pix,xsb],[ysb,ysb],'k',lw=1.5)
	ax.text(np.mean([xsb-sb_pix,xsb]),ysb+2,sb_str,ha='center',va='bottom',fontsize=plt.rcParams['font.size']-2)
	ax.set_xlim(left=xcen-crop/2,right=xcen+crop/2)
	ax.set_ylim(bottom=ycen-crop/2,top=ycen+crop/2)
	plt.savefig('../Plots/Vfields/'+gal_name+'_vfield_residual.pdf',bbox_inches='tight',metadata={'Creator': this_script})

	#add velocity contours
	levels = np.arange(-100,100+cstep,cstep)
	level0 = np.arange(0,cstep/2,cstep)
	ax.contour(vel-Vsys,levels=levels,colors='gray',linewidths=0.5)
	ax.contour(vel-Vsys,levels=level0,colors='k',linewidths=0.5)
	if rotmodel=='full':
		plt.savefig('../Plots/Vfields/'+gal_name+'_vfield_residual_contours.pdf',bbox_inches='tight',metadata={'Creator': this_script})
	else:
		plt.savefig('../Plots/Vfields/'+gal_name+'_vfield_residual_contours_'+rotmodel+'.pdf',bbox_inches='tight',metadata={'Creator': this_script})


main(gal_name)



