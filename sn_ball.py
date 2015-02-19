# -*- coding: utf-8 -*-
"""
This code is meant to give an estimated S/N for a single
balloon observation of a star, integrated over a specific
wavelength range.

Inputs: 
vmag: the v band magnitude of the star, G2V spectral type is
      assumed for now

lam0 - the blue edge of the desired bandpass, in micron
lam1 - the red edge of the desired bandpass, in micron
(a tophat banpass is assumed for now)

Optional Inputs: 
nd - this is an additional attenuation factor. 0.0
     would be no additinal, 0.90 would be an extra
     90%. 

Outputs: 
SN - the estimaged Signal and noise [S,N]. If there are saturated pixels,
     then it returns neg

Common Block: 
one common block is used to pass info to the function
used to integrate the PSF across a central pixel

Created by Avi Mandell, updated by Brett Morris
"""
import numpy as np
import pyfits
import os
from glob import glob
import cPickle
from urllib import urlopen
import time

def sgauss(x, sigma1):
    return (1./(2.*np.pi*sigma1**2))*np.exp(-(x**2)/(2*sigma1**2))

def deriv(y):
    return np.gradient(y)
    
def planck(wave, temp):
    '''
    Translated from: http://idlastro.gsfc.nasa.gov/ftp/pro/astro/planck.pro
    
 INPUT PARAMETERS: 
       WAVE   Scalar or vector giving the wavelength(s) in **Angstroms**
               at which the Planck function is to be evaluated.
       TEMP   Scalar giving the temperature of the planck function in degree K

 OUTPUT PARAMETERS:
       BBFLUX - Scalar or vector giving the blackbody flux (i.e. !pi*Intensity)
               in erg/cm^2/s/A in at the specified wavelength points.    
    
    '''
    bbflux = wave*0.
    # Gives the blackbody flux (i.e. PI*Intensity) ergs/cm2/s/a
    w = wave / 1.0e8                              #Angstroms to cm    
    # constants appropriate to cgs units.
    c1 =  3.7417749e-5                # =2*!DPI*h*c*c       
    c2 =  1.4387687                  # =h*c/k
    val =  c2/w/temp  
    bbflux = c1 / ( w**5 * ( np.exp(val)-1. ) )
    return bbflux*1.0e-8              # Convert to ergs/cm2/s/A

def downloadAndPickle():
    exodbPath = 'db'
    pklDatabaseName = os.path.join(exodbPath,'exoplanetDB.pkl')	 ## Name of exoplanet database C-pickle
    pklDatabasePaths = glob(pklDatabaseName)   ## list of files with the name pklDatabaseName in cwd
    csvDatabaseName = os.path.join(exodbPath,'exoplanets.csv')  ## Path to the text file saved from exoplanets.org
    csvDatabasePaths = glob(csvDatabaseName)

    '''First, check if there is an internet connection.'''

    '''If there's a previously archived database pickle in this current working 
        directory then use it, if not, grab the data from exoplanets.org in one big CSV file and make one.
        If the old archive is >14 days old, grab a fresh version of the database from exoplanets.org.
        '''
    if csvDatabasePaths == []:
        print 'No local copy of exoplanets.org database. Downloading one...'
        rawCSV = urlopen('http://www.exoplanets.org/csv-files/exoplanets.csv').read()
        saveCSV = open(csvDatabaseName,'w')
        saveCSV.write(rawCSV)
        saveCSV.close()
    else: 
        '''If the local copy of the exoplanets.org database is >14 days old, download a new one'''
        secondsSinceLastModification = time.time() - os.path.getmtime(csvDatabaseName) ## in seconds
        daysSinceLastModification = secondsSinceLastModification/(60*60*24*30)
        if daysSinceLastModification > 7:
            print 'Your local copy of the exoplanets.org database is >14 days old. Downloading a fresh one...'
            rawCSV = urlopen('http://www.exoplanets.org/csv-files/exoplanets.csv').read()
            saveCSV = open(csvDatabaseName,'w')
            saveCSV.write(rawCSV)
            saveCSV.close()
        else: print "Your local copy of the exoplanets.org database is <14 days old. That'll do."

    if len(pklDatabasePaths) == 0:
        print 'Parsing '+os.path.split(csvDatabaseName)[1]+', the CSV database from exoplanets.org...'
        rawTable = open(csvDatabaseName).read().splitlines()
        labels = rawTable[0].split(',')
        #labelUnits = rawTable[1].split(',')
        #rawTableArray = np.zeros([len(rawTable),len(labels)])
        exoplanetDB = {}
        planetNameColumn = np.arange(len(labels))[np.array(labels,dtype=str)=='NAME'][0]
        for row in range(1,len(rawTable)): 
            splitRow = rawTable[row].split(',')
            exoplanetDB[splitRow[planetNameColumn]] = {}	## Create dictionary for this row's planet
            for col in range(0,len(splitRow)):
                exoplanetDB[splitRow[planetNameColumn]][labels[col]] = splitRow[col]
        
        output = open(pklDatabaseName,'wb')
        cPickle.dump(exoplanetDB,output)
        output.close()
    else: 
        print 'Using previously parsed database from exoplanets.org...'
        inputFile = open(pklDatabaseName,'rb')
        exoplanetDB = cPickle.load(inputFile)
        inputFile.close()
    
    return exoplanetDB

db = downloadAndPickle()

def sn_ball(**kwargs):
    mag = kwargs.get('mag', 4.0)        # Object magnitude in a specific band
    band = kwargs.get('band', 'V').upper() # Band for magnitude entered (V,I,K)
    lam0 = kwargs.get('lam0', 1.29)     # Short-wavelength limit
    lam1 = kwargs.get('lam1', 1.51)     # Long-wavelength limit
    nd = kwargs.get('nd', 0)            # If no nuetral density, set it to 1.0
    verbose = kwargs.get('verbose', 1)  # If verbose not set, turn off report
    spt = kwargs.get('spt', 'G3').upper() # Spectral type
    ftime = kwargs.get('ftime', 7200.0) # Total observation time
    etime = kwargs.get('etime', 2.7)    # Integration exposure time
    pixel_lam = kwargs.get('pixlam', 0.0015) # pixel size in wavelength
    adapt = kwargs.get('adapt', 0)      # Use adaptive exposure time
    load = kwargs.get('load', 0)    

    opt = 1 if lam0 < 0.8 else 0
    mir = 1 if lam1 > 5.0 else 0
  
    # Define constants
    tel_diam   = 127.                   # in cm# TEEBall=127
    tel_area   = np.pi*(tel_diam/2.)**2 # telescope area in cm**2
    pixel_size = 0.18                   # pixel size in arcseconds TEEBall=0.18
    read_noise = 17.0                   # e-, H1RG=17.
    full_well  = 93000.                 # maximum number of electrons for 
                                        #    linearity H1RG=93K
    gain       = 2.15                   # e-/ADU, H1RG=2.15
    fixed_exptime=etime                 # [s] - H1RG min is 2.7s, next is 7.5
    dead_time  = 0.3                    # [s], reset time add to any exposure# 
                                        #    H1RG=0.3
    flight_time= ftime                  # [s]  total obseving time allocated
    ref_alt    = 36.6                   # [km] used for scintillation estimate
    airmass    = 1.5/143.	          #  Airmass for scintillation estimate 
                              #(denominator based on altitude see Bretts calcs)
    tel_eff    = 0.60*(1.0-nd)          # this is the total system efficiency, 
            # assumed not to depend on wavelength, effect of ND filter included
    Ttel = 220.0                        # Telescope temperature-TEEBall is 220K
    ntel = 0.1                          # Telescope emissivity after accounting 
               #for ohmic losses and spillover (0.9 for radio and up to 10 um)
    transname = '35K'		          # String for transmittance file - using 
                                        # 35K for balloon
    
    # wavelengths in microns for these 80% diameters
    lamD80      = np.arange(1, 6)
    # diameter of 80% of light in arcsec
    D80_of_lam  = 1.3*np.ones_like(lamD80)
    #  arcsec, size of PSF, assume it's a Gaussian
    #fwhm = np.interp(D80_of_lam, lamD80, 0.5*(lam0+lam1) )
    fwhm = np.interp(0.5*(lam0+lam1), lamD80, D80_of_lam)
    # here I am assuming that this "diameter of 80%" is about the FWHM
    
    phot_aper = 5.0*fwhm    # photometric aperture radius in arcseconds,
    #this should be close to optimal for background-limited observations,
    #see Naylor, 1998, MNRAS, 296, 339
    phot_area = phot_aper#**2
    phot_area_pix = (phot_aper/pixel_size)*((lam1-lam0)/pixel_lam)
    
    # If loading in saved arrays, skip file reading
    if not load:
        # first, let's read in the sky transmission spectrum. I have it
        # as a FITS binary table with wavelength in microns and [0,1] 
        # transmission
        skylam = np.ones(7000)*0.1 + 300.0
        if not opt:
            skytrans=pyfits.getdata(transname+'_0001.fits')
            skylam = skytrans[:, 0]*1000.   # Wavelength [nm]
            skyfrq = 1e7/skylam             # Frequency [cm-1]
            skyrad = skytrans[:, 1]   # Radiance [W/(cm2 sr cm-1)]
            skytrn = skytrans[:, 2]   # Transmittance [0-1]
        
        # This spectrum is really high resolution. 
        # We read in the stellar spectrum and convert it to photon units
        # The Kurucz models are in flux density ergs/s/cm**2/ster/nm,
        # wavelength is in nm
        DIR='kp00/'
        
        SP=['A0','A1','A2',
        'A3','A4','A5','A6','A7',
        'A8','A9','F0','F1','F2',
        'F3','F4','F5','F6','F7',
        'F8','F9','G0','G1','G2',
        'G3','G4','G5','G6','G7',
        'G8','G9','K0','K1','K2',
        'K3','K4','K5','K6','K7',
        'K8','K9','M0',
        'M1','M2','M3',
        'M4','M5']
        
        models=['kp00_9500.fits','kp00_9500.fits','kp00_9500.fits',
        'kp00_8250.fits','kp00_8250.fits','kp00_8250.fits','kp00_8250.fits',
        'kp00_8250.fits','kp00_7250.fits','kp00_7250.fits','kp00_7250.fits',
        'kp00_7250.fits','kp00_7250.fits','kp00_6500.fits','kp00_6500.fits',
        'kp00_6500.fits','kp00_6500.fits','kp00_6500.fits','kp00_6000.fits',
        'kp00_6000.fits','kp00_6000.fits','kp00_6000.fits','kp00_6000.fits',
        'kp00_5750.fits','kp00_5750.fits','kp00_5750.fits','kp00_5750.fits',
        'kp00_5750.fits','kp00_5250.fits','kp00_5250.fits','kp00_5250.fits',
        'kp00_5250.fits','kp00_5250.fits','kp00_4250.fits','kp00_4250.fits',
        'kp00_4250.fits','kp00_4250.fits','kp00_4250.fits','kp00_3750.fits',
        'kp00_3750.fits','kp00_3750.fits','kp00_3500.fits','kp00_3500.fits',
        'kp00_3500.fits','kp00_3500.fits','kp00_3500.fits']
        
        mtags=['G40','G40','G40',
        'G45','G45','G45','G45','G45',
        'G45','G45','G45','G45','G45',
        'G45','G45','G45','G45','G45',
        'G45','G45','G45','G45','G45',
        'G45','G45','G45','G45','G45',
        'G45','G45','G45','G45','G45',
        'G45','G45','G45','G45','G45',
        'G45','G45','G45',
        'G45','G45','G45',
        'G50','G50']
        
        # Read the selected spectral type
        model = models[spt == SP] 
        tag = mtags[spt == SP]
        star = pyfits.getdata(DIR+model)
        #tagnum = where(strcmp(tag,tag_names(star)))
        starflam = star[tag]*10.0            # Flux [ergs s-1 cm-2 nm-1]
        starlam  = star['WAVELENGTH']/10.0   # Wavelength [nm]
        h = 6.62606957e-27                 # Planck's constant [erg s photon-1]
        c = 2.99792458e17                    # Speed of light [nm/s]
        
        starfluxphotons=starflam*starlam/(h*c)  # Flux [photon s-1 cm-2 nm-1]
        
        # according to the gemini page, a star at V=0 is 9.71d3 
        # [photon s-1 cm-2 nm-1]  at a wavelength of 0.55 micron 
        # (http://www.gemini.edu/?q=node/10257). I=0 is 3.90d3, wavelength is 
        # 0.87 microns. K=0 is 4.5d2, wavelength is 2.2
        
        vband= (starlam > 545)*(starlam < 555)
        #iband= (starlam > 865)*(starlam < 875)
        jband= (starlam > 1105)*(starlam < 1349)
        kband= (starlam > 2195)*(starlam < 2205)
        
        if band == 'V':
            specnorm = np.mean(starfluxphotons[vband])/(9.71e3)
        #elif band == 'I':
        #    specnorm = np.mean(starfluxphotons[iband])/(3.90e3)
        elif band == 'J':
            specnorm = np.mean(starfluxphotons[jband])/(1.97e3)
        elif band in ['K', 'KS']:
            specnorm = np.mean(starfluxphotons[kband])/(4.50e2)
        
        # next, we will get the various sources of background photons tabulated
        # If observing in the optical, skip the thermal + OH contributions
        
        if not opt:
            # Telescope background
            telrad = planck(skylam*10.,Ttel)*ntel     # erg s-1 cm-2 sr-1 A-1
            telrad = telrad*np.abs(deriv(skylam*10.))  # erg s-1 cm-2 sr-1
            telrad = telrad*skylam/(h*c)              # ph  s-1 cm-2 sr-1
            telrad = telrad/((180.0/np.pi)**2)  # ph  s-1 cm-2 arcdeg-2
            telrad = telrad/(60*60.0)**2      # ph  s-1 cm-2 arcsec-2
            telrad = telrad*(100.0)**2         # ph  s-1  m-2 arcsec-2
            #telrad = 0.0
            
            # Sky Background - Generated by LBLRTM
            atmrad = skyrad.copy()          # W cm-2 sr-1 (cm-1)-1
            atmrad = 3.0e-9/skyfrq
            atmrad = 1.0e7*atmrad              # erg s-1 cm-2 sr-1 (cm-1)-1
            atmrad = atmrad*np.abs(deriv(skyfrq))    # erg s-1 cm-2 sr-1  
            atmrad = atmrad*skylam/(h*c)     # ph  s-1 cm-2 sr-1
            atmrad = atmrad/((180.0/np.pi)**2)  # ph  s-1 cm-2 arcdeg-2
            atmrad = atmrad/(60.0*60)**2      # ph  s-1 cm-2 arcsec-2
            atmrad = atmrad*(100.0)**2         # ph  s-1  m-2 arcsec-2
            
            thermalrad = telrad + atmrad     # ph  s-1  m-2 arcsec-2
            
        # Gemini Sky Background Data Files 
        # (http://www.gemini.edu/sciops/telescopes-and-sites/observing-
        #     condition-constraints/ir-background-spectra#Near-IR-short)
        # The files were manufactured starting from the sky transmission 
        # files generated by ATRAN (Lord, S. D., 1992, NASA Technical 
        # Memorandum 103957). These files were subtracted from unity to give 
        # an emissivity and then multiplied by a blackbody function of 
        # temperature 273 for Mauna Kea and 280 for Cerro Pachon. To these 
        # were added the OH emission spectrum (available from the European 
        # Southern Observatorys ISAAC web pages) a set of O2 lines near 1.3 
        # microns with estimated strengths based on observations at Mauna Kea, 
        # and the dark sky continuum (in part zodiacal light), approximated as 
        # a 5800K gray body times the atmospheric transmission and scaled to 
        # produce 18.2 mag/arcsec^2 in the H band, as observed on Mauna Kea by 
        # Maihara et al. (1993 PASP, 105, 940).
        
        # Im going to assume that OH lines are really only between 900 nm and 
        # 2300nm, so I can shift the spectrum to match the thermal component 
        # from our sky emission.
            
            mksky = pyfits.getdata('mk_skybg_zm_10_10_ph.fits')
            
            lamskyline = mksky['lam']
            skylines = mksky['lines']
            
            skynorm_mk = (lamskyline > 2200.)*(lamskyline < 2250.)
            skynorm_bal = (skylam  > 2200.)*(skylam < 2250.)
            skylines2 = skylines - (np.mean(skylines[skynorm_mk]) - 
                        np.mean(thermalrad[skynorm_bal])) 
                        # ph  s-1  m-2 arcsec-2 nm-1
            skylines2 = skylines2*np.abs(deriv(lamskyline[0])) 
                        # ph  s-1  m-2 arcsec-2
            
        elif opt:
            # Sky emission for optical wavelengths# from MODTRAN plots by Terry
            # Convert Rayleighs to ph s-1 -- 
            # from http://www.irf.se/~urban/avh/html/node13.html
            lamskyline = skylam
            skylines2 = 10.**(3.3+0.05*(50-ref_alt)-1.*skylam/450.) # kR nm-1
            skylines2 = skylines2*np.abs(deriv(skylam))   # kR
            skylines2 = skylines2*1.0e13/(4.*np.pi) # ph  s-1 m-2 sr-1
            skylines2 = skylines2/((180.0/np.pi)**2)  # ph  s-1 m-2 arcdeg-2
            skylines2 = skylines2/(60.0*60.0)**2      # ph  s-1 m-2 arcsec-2
            
    #loadskip:

    
    # define the part of the spectrum that we care about.    
    # The input bandpass definitions are in micron 
    
    inband= (skylam/1e3 >= lam0)*(skylam/1.0e3 <= lam1)
    inband = np.invert(inband)
    lambda_= skylam[inband]      # in nm
    if not opt:
        trans = skytrn[inband]      # Transmittance
    else:
        trans = lambda_*0.0+1.0        # If optical, than assume sky is clear
    
    # now we can scale the spectrum to the desired flux given the V magnitude
    
    starfluxphotons_mag=(starfluxphotons/specnorm)*10.0**(-mag*1.0/2.5)
    
    fjy = starfluxphotons_mag*(h*c/starlam) # ergs s-1 cm-2 nm-1
    fjy = fjy*np.abs(deriv(starlam)) # erg s-1 cm-2
    fjy = fjy/np.abs(deriv(c/starlam)) # erg s-1 cm-2 Hz-1
    fjy = fjy*1.0e23 # Jy
    
    # interpolate the stellar model onto the wavelength grid of
    # transmission spectrum, again keep microns and nm straight
    #starflux=interpol(starfluxphotons_mag, starlam, lambda_)
    #fjy = np.mean(interpol(fjy,starlam,lambda_))
    starflux=np.interp(lambda_, starlam, starfluxphotons_mag)
    fjy = np.mean(np.interp(lambda_, starlam, fjy))   
 
    # Flux as observed from observatory
    starflux_trn=starflux*trans # Flux [photon s-1 cm-2 nm-1] 
    
    # Now I integrate the flux density times the transmission and to get the
    # expected photons per second from the star. Starflux has units of [photon s-1 cm-2 nm-1]
    star_photons=np.sum(starflux_trn*np.abs(deriv(lambda_)))  # [photon s-1 cm-2]
    star_photons=star_photons*tel_area*tel_eff             # [photon s-1]
    
    # Sum up the thermal and sky backgrounds for the requested wavelengths
    if not opt: 
        thermalback = np.sum(thermalrad[inband]) 
    else:
        thermalback = 0.0
    
    if lam0 > 2.25:
        skyflux=0. 
    else: 
        inband_sky = (lamskyline/1000. >= lam0)*(lamskyline/1000. <= lam1)
        skyflux=np.sum(skylines2[inband_sky])
    
    totalsky=(thermalback+skyflux)*(tel_area/1.0e4)*tel_eff
    
    # this is now total background in photon/s/arcsec^2
    
    # now we need to figure out effective integration time
    
    # We assume that the PSF is a Gaussian with sigma = FWHM/2.36
    
    # sigma of gaussian, in pixels, assuming that the 80% diameter ~ FWHM
    sigma1=(fwhm/2.36)/pixel_size 
    
    central_pix_frac = sgauss(np.arange(101)-50, sigma1)[50] / \
                     (lam1-lam0)/pixel_lam
    
    # this is photons/s in central pixel, assuming PSF is aligned with center of pixel
    central_pix_flux=[central_pix_frac*star_photons, totalsky*pixel_size**2]
    
    # this is the exposure time - starts with initial value provided
    expt1=fixed_exptime
    
    # Loop to decrease exptime for saturation  
    adaptexp = True
    while adaptexp:
        # this is the number of exposures we get during flight
        nexposures=np.floor(flight_time/(expt1+dead_time))
        
        # this is the total effective integration time
        expt2=expt1*nexposures
        signal=expt2*star_photons
        
        # This is Nutzman & Charbonneau (2008) Eqn. 2 for uncertainty due to 
        # scintillation
        
        sigma_scint = signal*0.09*(airmass**(3/2.))/((tel_diam)**(2/3.)*\
                      np.sqrt(2.*expt2))*np.exp(-ref_alt/8.)
        
        noise1 = np.sqrt(signal+sigma_scint**2.+phot_area*expt2*totalsky + \
                 phot_area_pix*nexposures*(read_noise/np.sqrt(expt1/1.4))**2.)
        
        precision=[signal,noise1*np.sqrt(1.7)]
        
        # if we are non-linear within the exptime, then set the signal to negative
        # if using adaptive exptime, decrease exptime and recalculate
        cenpix = [expt1*np.array(central_pix_flux),full_well]
        decreasefactor = 0.75
        if np.sum(cenpix[0:1]) > full_well and expt1*decreasefactor > 1:
            #if not adapt:
            #    precision = -1.0*np.array(precision)
            #if adapt:
            expt1 *= decreasefactor 
            #    adaptexp = True
        else: 
            adaptexp = False
        
    # Print results
    if verbose:
        print '---------------------------------------------'
        print 'Spectral type: ', spt
        print 'Magnitude of the object: ', mag
        print 'Flux in Jy: ',fjy
        print 'Wavelength limits [um]: ', lam0, lam1
        print 'Neutral density filter: ', nd
        print '---------------------------------------------'
        print 'Object signal:  ', signal
        print 'Expected noise: ', noise1
        print 'Expected precision:   ', 1.0e6*noise1*np.sqrt(1.7)/signal
        print 'Fraction of Sat. Lim.:   ', np.sum(cenpix[0:1])/full_well
        print '---------------------------------------------'

def getplanetparams(planet, bandpass):
    '''
    Retrieve magnitude and spectral type for each target:
    
    Parameters
    ----------
    planet : str
        Name of planet
    '''
    return (('band', bandpass),
            ('mag', float(db[planet][bandpass])), 
            ('spt', temp2type(db[planet]['TEFF'])))

def temp2type(T_eff):
    '''
    Use lookup table, return spectral type for a given effective temperature
    '''
    lookup = np.genfromtxt('db/temp2type.txt', 
                           dtype=[('sptype','S2'), ('temp','i8')])
    closesttype = lookup['sptype'][np.argmin(np.abs(lookup['temp'] - 
                                   float(T_eff)))]
    return closesttype


kwargs = {
    'lam0': 1.29,
    'lam1': 1.51,
    'nd': 0,
    'verbose': 1,
    'ftime': 7200.0,
    'etime': 2.7,
    'pixlam': 0.0015,
    'adapt': 1,
    'load': 0
    }
    #'mag': 4.0,
    #'band': 'V',
    #'spt': 'G3',

def sn_ball_planet(planet, bandpass, **kwargs):
    '''
    For observations of `planet` in `bandpass`, run sn_ball()
    '''
    bandpasses = ['V', 'J', 'H', 'K', 'KS']
    if bandpass.upper() not in bandpasses:
        raise ValueError('Band "{0}" not in known bands: {1}'.format(bandpass,
                         ', '.join(bandpasses)))
    for key, value in getplanetparams(planet, bandpass):
        kwargs[key] = value
    return sn_ball(**kwargs)
    
print sn_ball_planet('WASP-6 b', 'K', **kwargs)
