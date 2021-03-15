#wrapper to fit all CO rotation curves
import argparse
parser=argparse.ArgumentParser(
    description='''''',
    epilog='''''')
parser.add_argument('gal_name',type=str, help='Name of galaxy to fit')

args=parser.parse_args()
gal_name = args.gal_name

def main(gal_name):
	r'''
	Main function that load data and calls functions for a single galaxy.

	Parameters
	----------
	gal_name : str
		name of galaxy to fit
		
	Returns
	-------
	
	Notes
	-----
	Required modules: astropy, matplotlib, numpy, pandas, os (called functions have their own required modules)
	Based on WW2020_Data_Tutorial.ipynb
	Author: R. C. Levy (rlevy.astro@gmail.com)
	Last updated: 2021-03-09
	Change log:
		2020-05-28 : file created, RCL
		2021-03-09 : added support to fit rotation-only model, RCL
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
	from FitTiltedRings import fit_tilted_rings
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['font.size'] = 14
	plt.rcParams['mathtext.rm'] = 'serif'
	plt.rcParams['mathtext.fontset'] = 'cm'

	#load the geometric parameters
	fname_gal_params = '../Data/Galaxy_Parameters.csv'
	gal_params = pd.read_csv(fname_gal_params) 
	idx = gal_params['Name'].values.tolist().index(gal_name)
	RA = gal_params['RA'][idx] #deg
	Dec = gal_params['Dec'][idx] #deg
	dist = gal_params['Distance'][idx] #Mpc 
	Inc = gal_params['Inc'][idx] #deg
	PA = gal_params['PA'][idx] #deg
	Vsys = float(gal_params['Vsys'][idx]) #km/s
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

	RA_pix,Dec_pix = wcs.all_world2pix(RA,Dec,1,ra_dec_order=True)

	#mask the CO velocity field
	#open the gausspeak data for masking
	peak = fits.open('../Data/gaussfit/'+gal_name+'.gausspk.fits',ignore_missing_end=True)
	epeak = fits.open('../Data/gaussfit/'+gal_name+'.egausspk.fits',ignore_missing_end=True)
	snr = peak[0].data/epeak[0].data
	vel[snr < 5] = np.nan


	#plot the CO velocity field
	#crop to 1' square
	crop = 1./60/hdr['CDELT2']
	xcen = hdr['CRPIX1']
	ycen = hdr['CRPIX2']
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
	rotmodel = 'full' #fit full model: rotation (cosine), radial (sine), systemic (constant)
	print('Fitting with '+rotmodel+' model')
	R,eR,Vrot,eVrot,Vrad,eVrad,dVsys,edVsys,chisq,chisqr,rms=fit_tilted_rings(gal_name,hdr,vel,evel,RA,Dec,PA,Inc,Vsys,rmEndRings,np.nan,save_dir,True,rotmodel)

	#now fit the rotation-only model
	R_ro,eR_ro,Vrot_ro,eVrot_ro,_,_,_,_,_,_,rms_ro=fit_tilted_rings(gal_name,hdr,vel,evel,RA,Dec,PA,Inc,Vsys,rmEndRings,np.nan,save_dir,False,'rotonly')

	#plot the velocity field
	plt.figure(1)
	plt.clf()
	ax=plt.subplot(projection=wcs)
	im=ax.imshow(vel-Vsys,origin='lower',cmap='RdYlBu_r',vmin=-100,vmax=100)
	ax.plot(RA_pix,Dec_pix,'kx')
	cb=plt.colorbar(im,shrink=0.7)
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
	ax.contour(snr,levels=np.arange(5,6,1),colors='gray',linewidths=1.0) 
	ax.add_patch(Ellipse((xtext,ytext),bmaj,bmin,bpa+90, ec='k',fc='k'))
	ax.plot([xsb-sb_pix,xsb],[ysb,ysb],'k',lw=1.5)
	ax.text(np.mean([xsb-sb_pix,xsb]),ysb+2,sb_str,ha='center',va='bottom',fontsize=plt.rcParams['font.size']-2)
	ax.set_xlim(left=xcen-crop/2,right=xcen+crop/2)
	ax.set_ylim(bottom=ycen-crop/2,top=ycen+crop/2)
	plt.savefig('../Plots/Vfields/'+gal_name+'_vfield.pdf',bbox_inches='tight')


	#plot the rings over the velcity field
	plt.figure(1)
	R_pix = R[1:]/3600/hdr['CDELT2'] #ring center
	R_shift = np.concatenate([R_pix[1:],np.array([np.inf])])
	R_outer = np.nanmean([R_pix,R_shift],axis=0)[0:-1]
	[ax.add_patch(Ellipse((RA_pix,Dec_pix),2*ro,2*ro*np.cos(np.radians(Inc)),PA+90, ec='k',fc='None',zorder=5)) for ro in R_outer]
	plt.savefig('../Plots/Vfields/'+gal_name+'_vfield_rings.pdf',bbox_inches='tight')


	#do a MonteCarlo over the geometric parameters
	ntrials = 50
	eevrot = np.zeros((ntrials,int(len(Vrot))))
	eevrad = np.zeros((ntrials,int(len(Vrot))))
	eedvsys = np.zeros((ntrials,int(len(Vrot))))
	eevrot_ro = np.zeros((ntrials,int(len(Vrot))))
	pa_off = 5. #deg
	inc_off = 5. #deg
	cen_off = 1.0/3600 #deg
	for j in range(ntrials):
		ra = RA + np.random.uniform(low=-cen_off,high=cen_off)
		dec = Dec + np.random.uniform(low=-cen_off,high=cen_off)
		pa = PA + np.random.uniform(low=-pa_off,high=pa_off)
		inc = Inc + np.random.uniform(low=-inc_off,high=inc_off)
		_,_,e1,_,e2,_,e0,_,_,_,_=fit_tilted_rings(gal_name,hdr,vel,evel,ra,dec,pa,inc,Vsys,rmEndRings,np.nan,save_dir,False,rotmodel)
		_,_,e1_ro,_,_,_,_,_,_,_,_=fit_tilted_rings(gal_name,hdr,vel,evel,ra,dec,pa,inc,Vsys,rmEndRings,np.nan,save_dir,False,'rotonly')

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


		if np.all(np.isnan(e1_ro))==True:
			#throw away trials that result in NaNs
			eevrot_ro[j,:]=np.nan
		else:
			if len(e1_ro) > len(Vrot_ro):
				#if MC has more rings, trim extra rings
				e1_ro=e1_ro[0:len(Vrot_ro)]
			elif len(e1_ro) < len(Vrot_ro):
				#if MC has fewer rings, pad with NaNs
				pad = np.zeros((int(len(Vrot_ro)-len(e1_ro)),))*np.nan
				e1_ro = np.concatenate((e1_ro,pad))				
			
			eevrot_ro[j,:]=e1_ro

			
	eVrot = np.nanstd(eevrot,axis=0)[0:len(R)]
	eVrad = np.nanstd(eevrad,axis=0)[0:len(R)]
	edVsys = np.nanstd(eedvsys,axis=0)[0:len(R)]
	eVrot_ro = np.nanstd(eevrot_ro,axis=0)[0:len(R_ro)]

	#remove the 0,0 point
	R=R[1:]
	eR=eR[1:]
	Vrot=Vrot[1:]
	eVrot=eVrot[1:]
	Vrad=Vrad[1:]
	eVrad=eVrad[1:]
	dVsys=dVsys[1:]
	edVsys=edVsys[1:]
	chisq=chisq[1:]
	chisqr=chisqr[1:]
	rms=rms[1:]
	R_ro=R_ro[1:]
	eR_ro=eR_ro[1:]
	Vrot_ro=Vrot_ro[1:]
	eVrot_ro=eVrot_ro[1:]
	rms_ro=rms_ro[1:]
	#remove all NaN points
	eR=eR[np.isnan(R)==False]
	Vrot=Vrot[np.isnan(R)==False]
	eVrot=eVrot[np.isnan(R)==False]
	Vrad=Vrad[np.isnan(R)==False]
	eVrad=eVrad[np.isnan(R)==False]
	dVsys=dVsys[np.isnan(R)==False]
	edVsys=edVsys[np.isnan(R)==False]
	chisq=chisq[np.isnan(R)==False]
	chisqr=chisqr[np.isnan(R)==False]
	rms=rms[np.isnan(R)==False]
	R=R[np.isnan(R)==False]
	eR_ro=eR_ro[np.isnan(R_ro)==False]
	Vrot_ro=Vrot_ro[np.isnan(R_ro)==False]
	eVrot_ro=eVrot_ro[np.isnan(R_ro)==False]
	rms_ro=rms_ro[np.isnan(R_ro)==False]
	R_ro=R_ro[np.isnan(R_ro)==False]

	eVrot = np.sqrt(eVrot**2+rms**2)
	eVrot_ro = np.sqrt(eVrot_ro**2+rms_ro**2)

	eVrot_u = np.nanmax([Vrot+eVrot,Vrot_ro+eVrot_ro],axis=0)-Vrot
	eVrot_l = Vrot-np.nanmin([Vrot-eVrot,Vrot_ro-eVrot_ro],axis=0)

	if gal_name=='NGC4150':
		#remove three rings where the gap is
		idx = [1,2,3]
		R = np.delete(R,idx)
		eR = np.delete(eR,idx)
		Vrot = np.delete(Vrot,idx)
		eVrot_u = np.delete(eVrot_u,idx)
		eVrot_l = np.delete(eVrot_l,idx)
		Vrad = np.delete(Vrad,idx)
		eVrad = np.delete(eVrad,idx)
		dVsys = np.delete(dVsys,idx)
		edVsys = np.delete(edVsys,idx)
		chisq = np.delete(chisq,idx)
		chisqr = np.delete(chisqr,idx)
		rms = np.delete(rms,idx)



	#make a second x-axis on the top in kpc
	def arcsec2kpc(x):
		return x/206265*dist*1E3
	def kpc2arcsec(x):
		return x/(1E3*dist)*206265

	R_kpc = arcsec2kpc(R) #kpc
	beam = hdr['BMAJ']*3600 #arcsec
	beam_kpc = arcsec2kpc(beam) #kpc

	#open the Ha rotation curve if it exists
	harc_file = '../Data/HalphaRotationCurves/'+gal_name+'.VrotHa.csv'
	if os.path.exists(harc_file)==True:
		dist_file = pd.read_csv('../Data/Beta_starComparison/CommonGalaxies.csv')
		idx = dist_file['Name'].values.tolist().index(gal_name)
		dist_R19 = dist_file['Distance_R19a'].values[idx] #Mpc
		harc = pd.read_csv(harc_file)
		Rha = harc['Rkpc'].values/(1E3*dist_R19)*206265 #arcsec
		Vha = harc['Vrot'].values
		eVha = harc['eVrot'].values
		Rha_kpc = arcsec2kpc(Rha) #kpc




	#plot the rotation curve
	plt.figure(1)
	plt.clf()
	plt.axhline(0,color='gray',lw=0.75)
	plt.axvspan(xmin=0,xmax=beam_kpc,fc='lightgray',ec='gray',alpha=0.25,hatch='\\\\',label='CO Beam FWHM',zorder=2)
	plt.axvspan(0.3,0.8,ec='None',fc='k',alpha=0.1,zorder=1)
	el1=plt.fill_between(R_kpc,Vrot+eVrot_u,Vrot-eVrot_l,fc='b',ec='None',alpha=0.3,zorder=3)
	if rotmodel=='full':
		el2=plt.fill_between(R_kpc,Vrad+eVrad,Vrad-eVrad,fc='r',ec='None',alpha=0.3,zorder=3)
		el3=plt.fill_between(R_kpc,dVsys+edVsys,dVsys-edVsys,fc='g',ec='None',alpha=0.3,zorder=3)
	l1,=plt.plot(R_kpc,Vrot,'bo-',markersize=4,label=r'CO V$_{\mathrm{rot}}$',zorder=4)
	if rotmodel=='full':
		l2,=plt.plot(R_kpc,Vrad,'r-',label=r'CO V$_{\mathrm{rad}}$',zorder=4)
		l3,=plt.plot(R_kpc,dVsys,'g-',label=r'CO $\Delta$V$_{\mathrm{sys}}$',zorder=4)
	xlim = plt.gca().get_xlim()
	if os.path.exists(harc_file)==True:
		el4=plt.fill_between(Rha_kpc,Vha+eVha,Vha-eVha,fc='rebeccapurple',ec='None',alpha=0.3,zorder=3)
		l4,=plt.plot(Rha_kpc,Vha,'^-',color='rebeccapurple',markersize=4,label=r'H$\alpha$ V$_{\mathrm{rot}}$')
	plt.minorticks_on()
	plt.xlabel('Radius (kpc)')
	plt.ylabel('Velocity (km s$^{-1}$)')
	plt.xlim(left=0,right=xlim[1])

	if os.path.exists(harc_file)==True:
		if rotmodel == 'full':
			plt.legend([(el1,l1),(el2,l2),(el3,l3),(el4,l4)],[l1.get_label(),l2.get_label(),l3.get_label(),l4.get_label()],
				loc='lower right',fontsize=plt.rcParams['font.size']-2)
		else:
			plt.legend([(el1,l1),(el4,l4)],[l1.get_label(),l4.get_label()],
				loc='lower right',fontsize=plt.rcParams['font.size']-2)
	else:
		if rotmodel=='full':
			plt.legend([(el1,l1),(el2,l2),(el3,l3)],[l1.get_label(),l2.get_label(),l3.get_label()],
				loc='lower right',fontsize=plt.rcParams['font.size']-2)
		else:
			plt.legend([(el1,l1)],[l1.get_label()],
				loc='lower right',fontsize=plt.rcParams['font.size']-2)

	plt.text(0.015,0.982,gal_name,transform=plt.gca().transAxes,fontsize=plt.rcParams['font.size']+2,va='top',ha='left')

	ax2 = plt.gca().secondary_xaxis('top',functions=(kpc2arcsec,arcsec2kpc))
	ax2.set_xlabel('Radius (arcsec)')
	ax2.minorticks_on()

	if rotmodel == 'full':
		plt.savefig('../Plots/CORotationCurves/'+gal_name+'_COrotcurve.pdf',bbox_inches='tight')
	else:
		plt.savefig('../Plots/CORotationCurves/'+gal_name+'_COrotcurve_'+rotmodel+'.pdf',bbox_inches='tight')

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
		fp.write('#column [0] = radius (arcsec),[1] = uncertainty on radius, [2] = rotation velocity (km/s), [3] = upper uncertainty on rotation velocity, [4] = lower uncertainty on rotation velocity, [5] = radial velocity (km/s), [6] = uncertainty on radial velocity, [7] = devitation from systemic velocity (km/s), [8] = uncertainty on deviation from systemic velocity, [9] = reduced chi squared of fit\n')
		for i in range(len(R)):
			fp.write('%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'
				%(R[i],eR[i],Vrot[i],eVrot_u[i],eVrot_l[i],Vrad[i],eVrad[i],dVsys[i],edVsys[i],chisqr[i]))


	#compare the rot-only and full fits
	#os.system('python CompareCORotationModels.py '+gal_name)

main(gal_name)



