def fit_tilted_rings(gal_name,header,velfield,evelfield,RA,Dec,PA,inc,Vsys):
	
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
	Author: R. C. Levy (rlevy@astro.umd.edu)
	Based on bestfit.m and ringfit.m by A. D. Bolatto and bestgetrings.m by R. C. Levy
	Last updated: 2019-05-20
	Change log:
		2019-05-17 : file created, RCL
		2019-05-20 : finished writing code, RCL
	
	Examples
	--------
	>>> from FitTiltedRings import fit_tilted_rings
	>>> R, eR, Vrot, eVrot, Vrad, eVrad, dVsys, edVsys, chisq, chisqr = fit_tilted_rings(gal_name,header,velfield,evelfield,RA,Dec,PA,inc,Vsys)
	
	'''

	#import modules
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	plt.rcParams['font.family'] = 'serif'

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
		idx = np.where(np.isnan(velfield)==0)
		#get coordinates and velocities of not-NaN points
		xv = xx[idx[0],idx[1]]
		yv = yy[idx[0],idx[1]]
		vv = velfield[idx[0],idx[1]]
		evv = evelfield[idx[0],idx[1]]*np.sqrt(beam_osamp)
		#save flattened corrdinates and velocities in an array
		vel_list = np.array([xv,yv,vv,evv]).T 
		return vel_list

	def getrings(vel_list,PA,inc,Vsys,bmaj):
		#find tilted rings
		#based on bestgetrings.m by R. C. Levy (based on getrings.m by A. D. Bolatto)
		#get galaxy axis ratio from inclination (simplistic)
		arat = np.cos(np.radians(inc))
		#get sin(inc) projection factor
		sini = np.sin(np.radians(inc))
		#convert pa from degrees to radian
		par = np.radians(PA+90.)
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
		#sort by increasing distance from the center
		idx_sort = np.argsort(rr)
		rr_sort = rr[idx_sort]
		vr_sort = vr[idx_sort]
		evr_sort = evr[idx_sort]
		xr_sort = xr[idx_sort]
		yr_sort = yr[idx_sort]
		#get rings every half beam
		r_beamspace = np.arange(0.,np.max(rr),bmaj/2)
		#unless there are fewer than min_ppr pixels in the ring
		min_ppr = 30 #minimum number of pixels in a ring
		nring = -1
		this_max_r = 0.0
		while this_max_r < np.max(rr_sort):
			nring = nring+1
			r_thisring = rr_sort[(rr_sort>=r_beamspace[nring]) & (rr_sort < r_beamspace[nring+1])]
			this_max_r = np.max(r_thisring)
			if len(r_thisring)==0:
				break
			idx = np.array([np.where(rr_sort==el)[0][0] for el in r_thisring])
			if len(idx) < min_ppr:
				dpix = min_ppr-len(idx) #number of pixels to add
				new_idx = idx[-1]+dpix
				if new_idx <= len(rr_sort):
					dr = rr_sort[idx[-1]+dpix]-rr_sort[idx[-1]] #change in r to add to all future rings (to keep half beam spacing)
					#update radii of rings
					r_beamspace[nring+1:] = r_beamspace[nring+1:]+dr
					#re-find radii and indices of this ring
					r_thisring = rr_sort[(rr_sort>=r_beamspace[nring]) & (rr_sort < r_beamspace[nring+1])]
					this_max_r = np.max(r_thisring)
					idx = np.array([np.where(rr_sort==el)[0][0] for el in r_thisring])
				else:
					nring=nring-1
					this_max_r = np.inf
			#clean up
		r_rings = r_beamspace[0:nring+1] #inner radii of the rings, arcsec
		return r_rings


	def ringfit(vel_list,r_rings,PA,inc,Vsys,gal_name):
		#fit a first order harmonic decomposition using the previously derived rings
		#based on ringfit.m by A. D. Bolatto

		#get galaxy axis ratio (simplistic)
		arat = np.cos(np.radians(inc))
		#get sin(inc) projection factor
		sini = np.sin(np.radians(inc))
		#flag points > flagcut*sigma frim fit
		flatcut = 10.
		#get pa in radians
		par = np.radians(PA+90.)
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

		#loop over rings
		for i in range(1,nrings):
			r_thisring = rr[(rr>=r_rings[i-1]) & (rr < r_rings[i])]
			idx = np.array([np.where(rr==el)[0][0] for el in r_thisring])
			npt[i]=len(idx)
			print('Ring '+str(i)+' with r_max = '+str(np.round(r_rings[i],2))+'" has '+str(npt[i])+' pixels')
			#grab coordinates and velocities for this ring
			xs = xr[idx]
			ys = yr[idx]
			vs = (vr[idx]-Vsys)/sini
			evs = evr[idx]/sini
			th = np.arctan2(ys,xs)
			ang = np.degrees(th)
			wt = 1./evs
			#fit rotation
			sol, err, rms = fitrotn(th,vs,wt)
			Vrot[i] = sol[1]
			Vrad[i] = sol[2]
			dVsys[i] = sol[0]
			eVrot[i] = err[1]
			eVrad[i] = err[2]
			edVsys[i] = err[0]
			#get mean radius of ring
			R[i] = np.mean(np.sqrt(xs**2+ys**2))
			eR[i] = np.std(np.sqrt(xs**2+ys**2))/np.sqrt(npt[i])
			#get model rotation velocity
			mod = dVsys[i]+Vrot[i]*np.cos(th)+Vrad[i]*np.sin(i)
			chisq[i]=np.sum(((vs-mod)*wt)**2)
			chisqr[i] = chisq[i]/(npt[i]-3.)

			#plot the fit for each ring
			plt.figure(4)
			plt.clf()
			ix = np.argsort(ang)
			allang = np.radians(np.linspace(-180.,180.,361))
			moda = dVsys[i]+Vrot[i]*np.cos(allang)+Vrad[i]*np.sin(allang) #full model
			rota = dVsys[i]+Vrot[i]*np.cos(allang)#rotation only
			plt.errorbar(ang[ix],vs[ix],yerr=evs[ix],fmt='bo',label='Data')
			plt.plot(np.degrees(allang),moda,'g-',label='Full model')
			plt.plot(np.degrees(allang),rota,'r-',label='Rotation only model')
			plt.xlabel('Azimuthal Angle (degrees)')
			plt.ylabel('Velocity (km s$^{-1}$)')
			plt.legend()
			plt.minorticks_on()
			t_str = gal_name+': Radius = '+str(np.round(R[i],2))+'", RMS = '+str(np.round(rms,2))+', $\chi^2_r$ = '+str(np.round(chisqr[i],2))
			plt.title(t_str)
			fig_name = 'ringfit/'+gal_name+'_ringfit_'+str(np.round(R[i],3))+'.pdf'
			plt.savefig(fig_name,bbox_inches='tight',overwrite=True)
			plt.close()

			#plot a summary plot
			plt.figure(5,figsize=(10,10))
			ax=plt.subplot(6,np.ceil(nrings/6),i)
			ax.plot(np.degrees(allang),rota,'g-')
			ax.plot(np.degrees(allang),moda,'r-',linewidth=2.)
			ax.errorbar(ang[ix],vs[ix],yerr=evs[ix],fmt='b.')
			ax.set_xlim(-180., 180.)
			ax.set_ylim(-150., 150.)
			ax.text(-180.,150.,str(i),fontsize=10.,verticalalignment='top')
			plt.minorticks_on()

		#save summary plot
		plt.figure(5)
		plt.tight_layout()
		fig_name = 'ringfit/'+gal_name+'_ringfit.pdf'
		plt.savefig(fig_name,overwrite=True)
		plt.close()

		return R,eR,Vrot,eVrot,Vrad,eVrad,dVsys,edVsys,chisq,chisqr

	def fitrotn(th,vs,wt):
		#fit rotation with first order harmonic decomposition
		#based on fitrotn.m by A. D. Bolatto
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


	#get beam major axis size from header
	bmaj = header['BMAJ']*3600. #arcsec

	#make list of coordinates and velocities
	vel_list = mkcolist(header,velfield,evelfield)

	#get radii of rings to use
	r_rings = getrings(vel_list,PA,inc,Vsys,bmaj)

	#remove existing ring fit figures
	ringfilename = 'ringfit/'+gal_name+'*.pdf'
	os.system('rm -f '+ringfilename)

	#fit the rotation curve
	R,eR,Vrot,eVrot,Vrad,eVrad,dVsys,edVsys,chisq,chisqr=ringfit(vel_list,r_rings,PA,inc,Vsys,gal_name)

	return R,eR,Vrot,eVrot,Vrad,eVrad,dVsys,edVsys,chisq,chisqr






