this_script = __file__

def fit_tilted_rings(gal_name,header,velfield,evelfield,RA,Dec,PA,inc,Vsys,rmEndRings,Rmax,save_dir,plotOn):
	
	r'''
	Derive rotation curve using tilted rings and a first order harmonic decomposition.

	Parameters
	----------
	gal_name : str
		name of galaxy to fit
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
		Position angle of the approaching side of the major axis measured east of north, in degrees
	inc : float
		Inclination of the galaxy to the line of sight, in degrees
	Vsys : float
		Systemic (AKA recessional) velocity of the center of the galaxy, in km/s
	rmEndRings : int
		Number of end rings to remove from rotation curve fit
	Rmax : float
		Maximum radius of rings, if NaN this limit is not imposed, in arcsec
	save_dir : str
		Path to directory to save output plots
	plotON: flag
	    If true or NOT GIVEN, will make and save intermediate plots
		If false, will NOT make nor save intermediate plots

		
	Returns
	-------
	R : array
		numpy array containing the radii of the fitted rings
	Vrot : array
		numpy array containing the fitted rotation velocity (cosine component)
	eVrot : array
		numpy array containing the uncertainty on Vrot
	Vrad : array
		numpy array containting the fitted radial velocity (sine component)
	eVrad : array
		numpy array containing the uncertainty on Vrad
	dVsys : array
		numpy array containing the fitted deviation from the systemic velocity (constant component)
	edVsys : array
		numpy array containing the uncertainty on Vsys
	
	Notes
	-----
	Required packages: numpy, matplotlib, os
	Author: R. C. Levy (rlevy.astro@gmail.com)
	Based on bestfit.m and ringfit.m by A. D. Bolatto and bestgetrings.m by R. C. Levy
	Last updated: 2021-02-15
	Change log:
		2019-05-17 : file created, RCL
		2019-05-20 : finished writing code, RCL
		2020-05-28 : added save_dir arg, other minor cleaning up, RCL
		2021-01-29 : added parameter to remove some number of bad end rings from fit, RCL
		2021-02-25 : return rms, add script name to pdf metadata, added Rmax keyword, RCL

	Examples
	--------
	>>> from FitTiltedRings import fit_tilted_rings
	>>> R, eR, Vrot, eVrot, Vrad, eVrad, dVsys, edVsys, chisq, chisqr, rms = fit_tilted_rings(gal_name,header,velfield,evelfield,RA,Dec,PA,inc,Vsys,rmEndRings,Rmax,save_dir,plotON)
	
	'''
	if (len(locals()) < 11):
		plotOn = True
	#import modules
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import sys
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['mathtext.rm'] = 'serif'
	plt.rcParams['mathtext.fontset'] = 'cm'

	def mkcolist(header,velfield,evelfield):
		#make a flattened list of the co velocity field and coordinates
		#based on mkcolist.m by A. D. Bolatto
		#get arrays of RA and Dec coordinates at each pixel
		RA_start = header['CRVAL1']-header['CRPIX1']*header['CDELT1']
		RA_end = RA_start+header['NAXIS1']*header['CDELT1']
		Dec_start = header['CRVAL2']-header['CRPIX2']*header['CDELT2']
		Dec_end = Dec_start+header['NAXIS2']*header['CDELT2']
		x = (np.linspace(RA_start,RA_end,header['NAXIS1'],endpoint=False)-RA)*3600. #offset in arcsec
		y = (np.linspace(Dec_start,Dec_end,header['NAXIS2'],endpoint=False)-Dec)*3600. #offset in arcsec
		#get beam oversampling factor
		beam_osamp = header['BMAJ']*header['BMIN']/header['CDELT2']**2
		#make a grid of coordinates
		xx,yy = np.meshgrid(x,y)
		#find where velocity field is not NaN
		idx = np.where(np.isnan(velfield)==False)
		#get coordinates and velocities of not-NaN points
		xv = xx[idx]
		yv = yy[idx]
		vv = velfield[idx]
		evv = evelfield[idx]*np.sqrt(beam_osamp)
		#save flattened corrdinates and velocities in an array
		vel_list = np.array([xv,yv,vv,evv]).T 
		return vel_list, beam_osamp

	def getrings(vel_list,PA,inc,Vsys,bmaj,rmEndRings,Rmax):
		#find tilted rings
		#based on bestgetrings.m by R. C. Levy (based on getrings.m by A. D. Bolatto)
		#get galaxy axis ratio from inclination (simplistic)
		arat = np.cos(np.radians(inc))
		#get sin(inc) projection factor
		sini = np.sin(np.radians(inc))
		#convert pa from degrees to radian
		#par = np.radians(PA-90.)
		par = np.radians(270-PA)
		#get rotation matrix
		rot = np.array([[np.cos(par), np.sin(par)],[-np.sin(par), np.cos(par)]])
		#get coordinates
		x = vel_list[:,0] #arcsec
		y = vel_list[:,1] #arcsec
		xy = np.array([x,y])
		#rotate coordinates
		xyr = np.matmul(rot,xy)
		xr = xyr[0,:]
		yr = xyr[1,:]/arat
		vr = vel_list[:,2]
		evr = vel_list[:,3]
		#find radius of each pixel from the center in arcsec
		rr = np.sqrt(xr**2+yr**2) #arcsec
		#sort by increasing distance from the center
		idx_sort = np.argsort(rr)
		rr_sort = rr[idx_sort]
		vr_sort = vr[idx_sort]
		evr_sort = evr[idx_sort]
		xr_sort = xr[idx_sort]
		yr_sort = yr[idx_sort]
		#get rings every half beam
		if np.isnan(Rmax)==True:
			r_beamspace = np.arange(0.,np.max(rr),bmaj/2)
		else:
			r_beamspace = np.arange(0.,Rmax,bmaj/2)
		#unless there are fewer than min_ppr pixels in the ring
		min_ppr = 30 #minimum number of pixels in a ring
		nring = 0
		this_max_r = 0.0
		while (this_max_r < (np.max(rr_sort)-bmaj/2)) & (nring+1 < len(r_beamspace)):
			r_thisring = rr_sort[(rr_sort>=r_beamspace[nring]) & (rr_sort < r_beamspace[nring+1])]
			nring = nring+1
			if len(r_thisring)== 0:
				break
			this_max_r = np.max(r_thisring)
			
			idx = np.array([np.where(rr_sort==el)[0][0] for el in r_thisring])
			if len(idx) < min_ppr:
				dpix = min_ppr-len(idx) #number of pixels to add
				new_idx = idx[-1]+dpix
				if new_idx >= len(rr_sort):
					nring=nring-1
					this_max_r = np.inf
				else:
					dr = rr_sort[idx[-1]+dpix]-rr_sort[idx[-1]] #change in r to add to all future rings (to keep half beam spacing)
					#update radii of rings
					r_beamspace[nring+1:] = r_beamspace[nring+1:]+dr
					#re-find radii and indices of this ring
					r_thisring = rr_sort[(rr_sort>=r_beamspace[nring-1]) & (rr_sort < r_beamspace[nring])]
					this_max_r = np.max(r_thisring)
					idx = np.array([np.where(rr_sort==el)[0][0] for el in r_thisring])
				#else:
					#nring=nring-1
				#	this_max_r = np.inf
			#clean up
		r_rings = r_beamspace[0:nring+1] #inner radii of the rings, arcsec

		#remove bad end rings
		if rmEndRings > 0:
			r_rings = r_rings[0:-rmEndRings]

		return r_rings


	def ringfit(vel_list,r_rings,PA,inc,Vsys,beam_osamp,gal_name,plotOn,save_dir):
		#fit a first order harmonic decomposition using the previously derived rings
		#based on ringfit.m by A. D. Bolatto

		#get galaxy axis ratio (simplistic)
		arat = np.cos(np.radians(inc))
		#get sin(inc) projection factor
		sini = np.sin(np.radians(inc))
		#flag points > flagcut*sigma frim fit
		flatcut = 10.
		#get pa in radians
		#par = np.radians(PA-90.)
		par = np.radians(270-PA)
		#get rotation matrix
		rot = np.array([[np.cos(par), np.sin(par)],[-np.sin(par), np.cos(par)]])
		#get coordinates
		x = vel_list[:,0]
		y = vel_list[:,1]
		xy = np.array([x,y])
		#rotate coordinates
		xyr = np.matmul(rot,xy)
		xr = xyr[0,:]
		yr = xyr[1,:]/arat
		vr = vel_list[:,2]
		evr = vel_list[:,3]
		#find radius of each pixel from the center in arcsec
		rr = np.sqrt(xr**2+yr**2)
		#get number of rings to fit
		nrings = len(r_rings)
		#make arrays to hold values
		R = np.zeros(nrings)
		eR = np.zeros(nrings)
		Vrot = np.zeros(nrings)
		eVrot = np.zeros(nrings)
		Vrad = np.zeros(nrings)
		eVrad = np.zeros(nrings)
		dVsys = np.zeros(nrings)
		edVsys = np.zeros(nrings)
		chisq = np.zeros(nrings)
		chisqr = np.zeros(nrings)
		npt = np.zeros(nrings)
		rms = np.zeros(nrings)

		#set up summary plot
		ncols=5
		nrows=int(np.ceil((nrings-1)/ncols))
		fig,ax=plt.subplots(nrows,ncols,sharex=True,sharey=True,num=5,figsize=(10,10))
		ax = np.ravel(ax)
		nplt = ncols*nrows

		#loop over rings
		for i in range(1,nrings):
			r_thisring = rr[(rr>=r_rings[i-1]) & (rr < r_rings[i])]
			idx = np.array([np.where(rr==el)[0][0] for el in r_thisring])
			#print(len(idx))
			#print(idx.shape)
			if len(idx) > 2:
				npt[i]=len(idx)

			#grab coordinates and velocities for this ring
				xs = xr[idx]
				ys = yr[idx]
				vs = (vr[idx]-Vsys)/sini
				evs = evr[idx]/sini
				th = np.arctan2(ys,xs)
				ang = np.degrees(th)
				wt = 1./evs
			#fit rotation
				sol, err, rms[i] = fitrotn(th,vs,wt)
				Vrot[i] = sol[1]
				Vrad[i] = sol[2]
				dVsys[i] = sol[0]


			#get uncertainties
				#get "base" uncertainties from the fitting
				eVrot[i] = err[1]; eVrad[i] = err[2]; edVsys[i] = err[0]
				# #get uncertainies by simply multiplying "base" by sqrt(beam_osamp)
				# eVrot[i],eVrad[i],edVsys[i]=uncert_beamosamp(err,beam_osamp)
				# #get uncertainties by fitting each half of the galaxy seperately
				# eVrot[i],eVrad[i],edVsys[i]=uncert_fithalfgal(ang,th,vs,wt)
				# #get uncertainties through bootstrapping
				# frac = 0.6 #fraction of pixels to use in each ring
				# eVrot[i],eVrad[i],edVsys[i]=uncert_bootstrap(frac,th,vs,wt)
				# #get uncertainties by Monte Carlo over uncertainties on each point
				# eVrot[i],eVrad[i],edVsys[i]=uncert_montecarlo_all(th,vs,wt,evs)
				# #get uncertainties by Monte Carlo over average uncertainty in each ring
				# eVrot[i],eVrad[i],edVsys[i]=uncert_montecarlo_one(th,vs,wt,evs)

			#get mean radius of ring
				R[i] = np.mean(np.sqrt(xs**2+ys**2))
				eR[i] = np.std(np.sqrt(xs**2+ys**2))#/np.sqrt(npt[i])
			#get model rotation velocity
				mod = dVsys[i]+Vrot[i]*np.cos(th)+Vrad[i]*np.sin(th)
				chisq[i]=np.sum(((vs-mod)*wt)**2)
				chisqr[i] = chisq[i]/(npt[i]/np.sqrt(beam_osamp)-3.)
			#if chisqr is negative, remove this point
				if chisqr[i] <= 0:
					R[i]=np.nan
					eR[i]=np.nan
					Vrot[i]=np.nan
					eVrot[i]=np.nan
					Vrad[i]=np.nan
					eVrad[i]=np.nan
					dVsys[i]=np.nan
					edVsys[i]=np.nan
					chisq[i]=np.nan
					chisqr[i]=np.nan
					rms[i]=np.nan


				if (plotOn == True):
					print('Ring '+str(i)+' with r_max = '+str(np.round(r_rings[i],2))+
					      '" has '+str(int(npt[i]))+' pixels, RMS = '
					      +str(np.round(rms[i],2))+', Chi^2_r = '+str(np.round(chisqr[i],2)))
					#plot the fit for each ring
					plt.figure(4)
					plt.clf()
					ix = np.argsort(ang)
					allang = np.radians(np.linspace(-180.,180.,361))
					rota = dVsys[i]+Vrot[i]*np.cos(allang)#rotation only
					erota = np.sqrt(edVsys[i]**2+eVrot[i]**2*np.cos(allang)**2)
					moda = dVsys[i]+Vrot[i]*np.cos(allang)+Vrad[i]*np.sin(allang) #full model
					emoda = np.sqrt(edVsys[i]**2+eVrot[i]**2*np.cos(allang)**2+eVrad[i]**2*np.sin(allang)**2)
					plt.errorbar(ang[ix],vs[ix],yerr=evs[ix],fmt='bo',label='Data',alpha=0.3)
					plt.fill_between(np.degrees(allang),rota-erota,rota+erota,ec=None,fc='r',alpha=0.3)					
					plt.fill_between(np.degrees(allang),moda-emoda,moda+emoda,ec=None,fc='g',alpha=0.3)
					plt.plot(np.degrees(allang),rota,'r-',label='Rotation only model')
					plt.plot(np.degrees(allang),moda,'g-',label='Full model')
					plt.xlabel('Azimuthal Angle (degrees)')
					plt.ylabel('Velocity (km s$^{-1}$)')
					plt.legend(loc='upper right',fontsize=8)
					plt.text(0.025,0.975,'V$_\mathrm{rot}$ = %.1f $\pm$ %.1f km s$^{-1}$\nV$_\mathrm{rad}$ = %.1f $\pm$ %.1f km s$^{-1}$\n$\Delta$V$_\mathrm{sys}$ = %.1f $\pm$ %.1f km s$^{-1}$' %(Vrot[i],eVrot[i],Vrad[i],eVrad[i],dVsys[i],edVsys[i]),
						transform=plt.gca().transAxes,ha='left',va='top',fontsize=8)
					plt.minorticks_on()
					t_str = gal_name+': Radius = '+str(np.round(R[i],1))+'", RMS = '+str(np.round(rms[i],1))+', $\chi^2_r$ = '+str(np.round(chisqr[i],1))
					plt.title(t_str)
					fig_name = save_dir+gal_name+'_ringfit_'+str(np.round(R[i],3))+'.pdf'
					plt.savefig(fig_name,bbox_inches='tight',metadata={'Creator':this_script})
					plt.close()

					#plot a summary plot
					ax[i-1].plot(np.degrees(allang),rota,'g-',zorder=9)
					ax[i-1].plot(np.degrees(allang),moda,'r-',linewidth=2.,zorder=10)
					ax[i-1].errorbar(ang[ix],vs[ix],yerr=evs[ix],fmt='b',marker='None',linestyle='None',alpha=0.2)
					ax[i-1].plot(ang[ix],vs[ix],'b.')
					ax[i-1].set_xlim(-180., 180.)
					ax[i-1].set_ylim(-200., 200.)
					ax[i-1].text(0.05,0.95,str(i),fontsize=8.,verticalalignment='top',transform=ax[i-1].transAxes)
					ax[i-1].minorticks_on()

					if i==int((np.ceil(nrings/ncols)-1)*ncols+1):
						ax[i-1].set_xlabel('Azimuthal Angle',fontsize=8)
						ax[i-1].set_ylabel('Velocity',fontsize=8)
			else:
				R[i]=np.nan
				eR[i]=np.nan
				Vrot[i] = np.nan
				eVrot[i] = np.nan
				Vrad[i] = np.nan
				eVrad[i] = np.nan
				dVsys[i] = np.nan
				edVsys[i] = np.nan
				chisq[i] = np.nan
				chisqr[i] = np.nan
				rms[i] = np.nan

		if plotOn==True:
			di = nplt-i-1
			for axx in ax[nplt-di-1:]:
				axx.remove()
			#save summary plot
			plt.figure(5)
			fig_name = save_dir+gal_name+'_ringfit.pdf'
			plt.savefig(fig_name,bbox_inches='tight',metadata={'Creator':this_script})
			plt.close()

		return R,eR,Vrot,eVrot,Vrad,eVrad,dVsys,edVsys,chisq,chisqr,rms

	def fitrotn(th,vs,wt):
		#fit rotation with first order harmonic decomposition
		#based on fitrotn.m by A. D. Bolatto
		#needs exception for when N = 1 ************
		w2 = wt**2
		N = len(th)
		D = 1./N*np.sum(w2)
		X1 = np.mean(np.cos(th))
		X2 = np.mean(np.sin(th))
		Y=np.mean(vs)
		s1=1./(N-1)*np.sum(w2*(np.cos(th)-X1)**2)/D
		s12=1./(N-1)*np.sum(w2*(np.cos(th)-X1)*(np.sin(th)-X2))/D
		s2=1./(N-1)*np.sum(w2*(np.sin(th)-X2)**2)/D
		sy=1./(N-1)*np.sum(w2*(vs-Y)**2)/D
		s1y=1./(N-1)*np.sum(w2*(np.cos(th)-X1)*(vs-Y))/D
		s2y=1./(N-1)*np.sum(w2*(np.sin(th)-X2)*(vs-Y))/D
		r12=s12/np.sqrt(s1*s2);
		r11=1
		r22=1
		r21=r12
		r1y=s1y/np.sqrt(s1*sy)
		r2y=s2y/np.sqrt(s2*sy)
		r=np.array([[r11,r12],[r21,r22]])
		rm=np.linalg.inv(r)

		a1=np.sqrt(sy/s1)*(r1y*rm[0,0]+r2y*rm[0,1])
		a2=np.sqrt(sy/s2)*(r1y*rm[1,0]+r2y*rm[1,1])
		a0=np.sum(w2*(vs-a1*np.cos(th)-a2*np.sin(th)))/np.sum(w2)
		sol=np.array([a0,a1,a2])

		sa1=1/(N-1)/s1*rm[0,0]/D
		sa2=1/(N-1)/s2*rm[1,1]/D
		sa0=(1/N+1/(N-1)*(X1**2/s1*rm[0,0]+X2**2/s2*rm[0,1]+2*X1*X2/np.sqrt(s1*s2)*rm[0,1]))/D
		err=np.sqrt(np.array([sa0,sa1,sa2]))
		
		rms=np.std(vs-sol[1]*np.cos(th)-sol[2]*np.sin(th)-sol[0])

		return sol, err, rms

	def uncert_beamosamp(err,beam_osamp):
		#muliply uncertainties from the fit by sqrt(beam_osamp)
		eVrot = err[1]*np.sqrt(beam_osamp)
		eVrad = err[2]*np.sqrt(beam_osamp)
		edVsys = err[0]*np.sqrt(beam_osamp)
		return eVrot, eVrad, edVsys

	def uncert_fithalfgal(ang,th,vs,wt):
		#fit with data on half of the galaxy to get the uncertainties 
		idx_p = np.where(ang >= 0)[0]
		idx_n = np.where(ang <= 0)[0]
		if len(idx_p) > 3:
			sol_p, _, _ = fitrotn(th[idx_p],vs[idx_p],wt[idx_p])
		else:
			sol_p = np.array([np.nan,np.nan,np.nan])
		if len(idx_n > 3):
			sol_n, _, _ = fitrotn(th[idx_n],vs[idx_n],wt[idx_n])
		else:
			sol_n = np.array([np.nan,np.nan,np.nan])
		esol = np.nanmean([np.abs(sol-sol_p),np.abs(sol-sol_n)],axis=0)
		eVrot = esol[1]
		eVrad = esol[2]
		edVsys = esol[0]
		return eVrot,eVrad,edVsys

	def uncert_bootstrap(frac,th,vs,wt):
		#get uncertaities by bootstrapping
		ntrials = 100
		sol = np.zeros((3,ntrials))
		for j in range(ntrials):
			idx = np.random.uniform(low=0,high=len(th)-1,size=int(frac*len(th))).astype(int)
			sol[:,j], _, _ = fitrotn(th[idx],vs[idx],wt[idx])
		eVrot = np.std(sol[1,:])
		eVrad = np.std(sol[2,:])
		edVsys = np.std(sol[0,:])
		return eVrot,eVrad,edVsys

	def uncert_montecarlo_all(th,vs,wt,evs):
		#get uncertainties by allowing all points to vary randomly within their uncertainties
		ntrials = 100
		sol = np.zeros((3,ntrials))
		for j in range(ntrials):
			voff = np.random.uniform(low=-evs,high=evs)
			sol[:,j], _, _ = fitrotn(th,vs+voff,wt)
		eVrot = np.std(sol[1,:])
		eVrad = np.std(sol[2,:])
		edVsys = np.std(sol[0,:])
		return eVrot,eVrad,edVsys

	def uncert_montecarlo_one(th,vs,wt,evs):
		#get uncertainties by allowing all points to vary by a single value within their uncertainties
		ntrials = 100
		sol = np.zeros((3,ntrials))
		for j in range(ntrials):
			voff = np.random.uniform(low=-np.mean(evs),high=np.mean(evs))
			sol[:,j], _, _ = fitrotn(th,vs+voff,wt)
		eVrot = np.std(sol[1,:])
		eVrad = np.std(sol[2,:])
		edVsys = np.std(sol[0,:])
		return eVrot,eVrad,edVsys

	#get beam major axis size from header
	bmaj = header['BMAJ']*3600. #arcsec

	#make list of coordinates and velocities
	vel_list, beam_osamp = mkcolist(header,velfield,evelfield)

	#get radii of rings to use
	r_rings = getrings(vel_list,PA,inc,Vsys,bmaj,rmEndRings,Rmax)
	if len(r_rings) > 2:
		
		if plotOn==True:
			#remove existing ring fit figures
			ringfilename = save_dir+gal_name+'*.pdf'
			os.system('rm -f '+ringfilename)

	#fit the rotation curve
		R,eR,Vrot,eVrot,Vrad,eVrad,dVsys,edVsys,chisq,chisqr,rms=ringfit(vel_list,r_rings,PA,inc,Vsys,beam_osamp,gal_name,plotOn,save_dir)

	else:
		print('Too few rings for fitting.')
		R = np.nan
		eR = np.nan
		Vrot = np.nan
		eVrot = np.nan
		Vrad = np.nan
		eVrad = np.nan
		dVsys = np.nan
		edVsys = np.nan
		chisq = np.nan
		chisqr = np.nan
		rms=np.nan

	return R,eR,Vrot,eVrot,Vrad,eVrad,dVsys,edVsys,chisq,chisqr,rms






