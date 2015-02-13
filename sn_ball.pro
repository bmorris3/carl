;--------------------------------------------------------------------------------------------
;PURPOSE:  This code is meant to give an estimated S/N for a single
;          balloon observation of a star, integrated over a specific
;          wavelength range.
;
;INPUT: vmag: the v band magnitude of the star, G2V spectral type is
;       assumed for now
;
;       lam0 - the blue edge of the desired bandpass, in micron
;       lam1 - the red edge of the desired bandpass, in micron
;              (a tophat banpass is assumed for now)
;
;OPTIONAL INPUTS: nd - this is an additional attenuation factor. 0.0
;                      would be no additinal, 0.90 would be an extra
;                      90%. 
; 
;OUTPUTS: SN - the estimaged Signal and noise [S,N]. If there are saturated pixels,
;              then it returns neg
;
;COMMON BLOCK: one common block is used to pass info to the function
;used to integrate the PSF across a central pixel
;
;LAST UPDATED - mandell, Feb 2015
;--------------------------------------------------------------------------------------------
; Functions used
FUNCTION sgauss, x  & common gaussblock, sigma1 & return, (1./(2.*!pi*sigma1^2d))*exp(-(x^2.)/(2*sigma1^2d)) & end

;--------------------------------------------------------------------------------------------
PRO sn_ball, mag, lam0, lam1, precision, cenpix=cenpix, band=band, nd=nd, ftime=ftime, etime=etime, verbose=verbose, spt=spt, load=load, fjy=fjy, pixlam=pixlam, adapt=adapt

  ; Usage example
  ; sn_ball, 9.5, 1.1, 1.3, result, spt='G0'

  ; Define defaults
  if keyword_set(mag) eq 0 then mag=4.0      ; Object magnitude in a specific band
  if keyword_set(band) eq 0 then band = 'V'   ; Band for magnitude entered (V,I,or K)
  if keyword_set(lam0) eq 0 then lam0=1.29    ; Short-wavelength limit
  if keyword_set(lam1) eq 0 then lam1=1.51    ; Long-wavelength limit
  if keyword_set(nd) eq 0 then nd=0.          ; if the nuetral density is not set, then set it to 1.0
  if keyword_set(verbose) eq 0 then verbose=0 ; if verbose keyword is not set, turn off final report
  if keyword_set(spt) eq 0 then spt='G3'      ; Spectral type
  if keyword_set(ftime) eq 0 then ftime=7200. ; Total observation time
  if keyword_set(etime) eq 0 then etime=2.7   ; Integration exposure time
  if keyword_set(pixlam) eq 0 then pixlam=0.0015   ;pixel size in wavelength
  if keyword_set(adapt) eq 0 then adapt = 0    ; Use adaptive exposure time
  pixel_lam = double(pixlam)
  
  common gaussblock, sigma1
  common savblock, specnorm, starfluxphotons, skylam, skyfrq, skytrn, skyrad, starlam, thermalrad, lamskyline, skylines2, c, h

  if lam0 lt .8 then opt = 1 else opt = 0
  if lam1 gt 5. then mir = 1 else mir = 0

  ; Define constants
  tel_diam   = 127.                 ;in cm; TEEBall=127
  tel_area   = (tel_diam/2.)^(2.)*!pi ;telescope area in cm**2
  pixel_size = 0.18                  ;pixel size in arcseconds; TEEBall=0.18
  read_noise = 17.0                 ;e-, H1RG=17.
  full_well  = 93000.              ; maximum number of electrons for linearity; H1RG=93K
  gain       = 2.15                  ;e-/ADU, H1RG=2.15
  fixed_exptime=etime                ;in seconds - H1RG min is 2.7s, next is 7.5s
  dead_time  = 0.3d                ;seconds, reset time to add to any exposure; H1RG=0.3
  flight_time= ftime                ;seconds, total obseving time allocated
  ref_alt    = 36.6                  ;in km, used for scintillation estimate
  airmass    = 1.5d/143.	    ; Airmass for scintillation estimate (denominator based on altitude see Bretts calcs)
  tel_eff    = 0.60*(1.-nd)          ;this is the total system efficiency, assumed not to depend on wavelength, effect of ND filter included
  Ttel = 220d                       ; Telescope temperature - TEEBall is 220K
  ntel = 0.1d                       ; Telescope emissivity after accounting for ohmic losses and spillover (0.9 for radio and up to 10 um)
  transname = '35K'		    ; String for transmittance file - using 35K for balloon

  lamD80      =1.*[1,   2,   3,   4,   5]; wavelengths in microns for these 80% diameters
  D80_of_lam  = [1.3, 1.3, 1.3, 1.3, 1.3]; diameter of 80% of light in arcsec
  fwhm=interpol(D80_of_lam, lamD80, .5*(lam0+lam1) ); arcsec, size of PSF, assume it's a Gaussian
                               ;here I am assuming that this "diameter of 80%" is about the FWHM
  
  phot_aper=5.0*fwhm ;photometric aperture radius in arcseconds,
  ;this should be close to optimal for background-limited observations,
  ;see Naylor, 1998, MNRAS, 296, 339
  phot_area=phot_aper^2
  phot_area_pix=(phot_aper/pixel_size)*((lam1-lam0)/pixel_lam)

  ;If loading in saved arrays, skip file reading
  if keyword_set(load) then goto,loadskip
  
  ;first, let's read in the sky transmission spectrum. I have it
  ;as a FITS binary table with wavelength in microns and [0,1] transmission
  skylam = findgen(7000)*0.1 + 300.
  if not opt then begin
    skytrans=readfits(transname+'_0001.fits',/silent)
    skylam = reform(skytrans[0,*])*1000.   ; Wavelength [nm]
    skyfrq = 1d7/skylam              ; Frequency [cm-1]
    skyrad = reform(skytrans[1,*])   ; Radiance [W/(cm2 sr cm-1)]
    skytrn = reform(skytrans[2,*])   ; Transmittance [0-1]
  endif
    
  ;This spectrum is really high resolution. 
  ;We read in the stellar spectrum and convert it to photon units
  ;The Kurucz models are in flux density ergs/s/cm**2/ster/nm,
  ;wavelength is in nm
  DIR='kp00/'
  
  SP=['A0','A1','A2',$
      'A3','A4','A5','A6','A7',$
      'A8','A9','F0','F1','F2',$
      'F3','F4','F5','F6','F7',$
      'F8','F9','G0','G1','G2',$
      'G3','G4','G5','G6','G7',$
      'G8','G9','K0','K1','K2',$
      'K3','K4','K5','K6','K7',$
      'K8','K9','M0',$
      'M1','M2','M3',$
      'M4','M5']
    
  models=['kp00_9500.fits','kp00_9500.fits','kp00_9500.fits',$
          'kp00_8250.fits','kp00_8250.fits','kp00_8250.fits','kp00_8250.fits','kp00_8250.fits',$
          'kp00_7250.fits','kp00_7250.fits','kp00_7250.fits','kp00_7250.fits','kp00_7250.fits',$
          'kp00_6500.fits','kp00_6500.fits','kp00_6500.fits','kp00_6500.fits','kp00_6500.fits',$
          'kp00_6000.fits','kp00_6000.fits','kp00_6000.fits','kp00_6000.fits','kp00_6000.fits',$
          'kp00_5750.fits','kp00_5750.fits','kp00_5750.fits','kp00_5750.fits','kp00_5750.fits',$
          'kp00_5250.fits','kp00_5250.fits','kp00_5250.fits','kp00_5250.fits','kp00_5250.fits',$
          'kp00_4250.fits','kp00_4250.fits','kp00_4250.fits','kp00_4250.fits','kp00_4250.fits',$
          'kp00_3750.fits','kp00_3750.fits','kp00_3750.fits',$
          'kp00_3500.fits','kp00_3500.fits','kp00_3500.fits',$
          'kp00_3500.fits','kp00_3500.fits']

  mtags=['G40','G40','G40',$
        'G45','G45','G45','G45','G45',$
        'G45','G45','G45','G45','G45',$
        'G45','G45','G45','G45','G45',$
        'G45','G45','G45','G45','G45',$
        'G45','G45','G45','G45','G45',$
        'G45','G45','G45','G45','G45',$
        'G45','G45','G45','G45','G45',$
        'G45','G45','G45',$
        'G45','G45','G45',$
        'G50','G50']

  ; Read the selected spectral type
  model = models[where(strcmp(spt,SP))] & tag = (mtags[where(strcmp(spt,SP))])[0]
  star=mrdfits(dir+model[0],1,hdr,/silent) & tagnum = where(strcmp(tag,tag_names(star)))
  starflam = star.(tagnum)*10d            ; Flux [ergs s-1 cm-2 nm-1]
  starlam  = star.wavelength/10d          ; Wavelength [nm]
  h = 6.62606957d-27                      ; Planck's constant [erg s photon-1]
  c = 2.99792458d17                       ; Speed of light [nm/s]
  
  starfluxphotons=starflam*starlam/(h*c)  ; Flux [photon s-1 cm-2 nm-1]

  ;according to the gemini page, a star at V=0 is 9.71d3 [photon s-1 cm-2 nm-1] 
  ;at a wavelength of 0.55 micron (http://www.gemini.edu/?q=node/10257).
  ; I=0 is 3.90d3, wavelength is 0.87 microns. K=0 is 4.5d2, wavelength is 2.2
  
  vband=where(starlam gt 545 and starlam lt 555)
  iband=where(starlam gt 865 and starlam lt 875)
  kband=where(starlam gt 2195 and starlam lt 2205)

  if strcmp(band,'V') then specnorm=mean(starfluxphotons[vband])/(9.71d3)
  if strcmp(band,'I') then specnorm=mean(starfluxphotons[iband])/(3.90d3)
  if strcmp(band,'K') then specnorm=mean(starfluxphotons[kband])/(4.50d2)

  ;next, we will get the various sources of background photons tabulated
  ; If observing in the optical, we can skip the thermal and OH contributions

  if opt then goto,loadskip

  ; Telescope background
  telrad = planck(skylam*10.,Ttel)*ntel     ; erg s-1 cm-2 sr-1 A-1
  telrad = telrad*abs(deriv(skylam*10.))  ; erg s-1 cm-2 sr-1
  telrad = telrad*skylam/(h*c)              ; ph  s-1 cm-2 sr-1
  telrad = telrad/((180d/!PI)^2d)  ; ph  s-1 cm-2 arcdeg-2
  telrad = telrad/(60d*60d)^2      ; ph  s-1 cm-2 arcsec-2
  telrad = telrad*(100d)^2         ; ph  s-1  m-2 arcsec-2
  ;telrad = 0.d
  
  ; Sky Background - Generated by LBLRTM
  atmrad = reform(skyrad)          ; W cm-2 sr-1 (cm-1)-1
  atmrad = 3.d-9/skyfrq
  atmrad = 1d7*atmrad              ; erg s-1 cm-2 sr-1 (cm-1)-1
  atmrad = atmrad*abs(deriv(skyfrq))    ; erg s-1 cm-2 sr-1  
  atmrad = atmrad*skylam/(h*c)     ; ph  s-1 cm-2 sr-1
  atmrad = atmrad/((180d/!PI)^2d)  ; ph  s-1 cm-2 arcdeg-2
  atmrad = atmrad/(60d*60d)^2      ; ph  s-1 cm-2 arcsec-2
  atmrad = atmrad*(100d)^2         ; ph  s-1  m-2 arcsec-2
  
  thermalrad = telrad + atmrad     ; ph  s-1  m-2 arcsec-2

  ; Gemini Sky Background Data Files (http://www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/ir-background-spectra#Near-IR-short)
  ;
  ; The files were manufactured starting from the sky transmission 
  ; files generated by ATRAN (Lord, S. D., 1992, NASA Technical 
  ; Memorandum 103957). These files were subtracted from unity to give 
  ; an emissivity and then multiplied by a blackbody function of 
  ; temperature 273 for Mauna Kea and 280 for Cerro Pachon. To these 
  ; were added the OH emission spectrum (available from the European 
  ; Southern Observatorys ISAAC web pages) a set of O2 lines near 1.3 
  ; microns with estimated strengths based on observations at Mauna Kea, 
  ; and the dark sky continuum (in part zodiacal light), approximated as 
  ; a 5800K gray body times the atmospheric transmission and scaled to 
  ; produce 18.2 mag/arcsec^2 in the H band, as observed on Mauna Kea by 
  ; Maihara et al. (1993 PASP, 105, 940).
  
  ; Im going to assume that OH lines are really only between 900 nm and 
  ; 2300nm, so I can shift the spectrum to match the thermal component 
  ; from our sky emission.

  mksky=mrdfits('mk_skybg_zm_10_10_ph.fits',1,/silent)

  lamskyline=mksky.lam
  skylines=mksky.lines

  skynorm_mk = where(lamskyline gt 2200. and lamskyline lt 2250.)
  skynorm_bal = where(skylam  gt 2200. and skylam lt 2250.)
  skylines2 = skylines - (mean(skylines(skynorm_mk)) - mean(thermalrad(skynorm_bal))) ; ph  s-1  m-2 arcsec-2 nm-1
  skylines2 = skylines2*abs(deriv(lamskyline)) ; ph  s-1  m-2 arcsec-2
  
  if opt then begin
    ;Sky emission for optical wavelengths; from MODTRAN plots by Terry
    ; Convert Rayleighs to ph s-1 -- from http://www.irf.se/~urban/avh/html/node13.html
    lamskyline = skylam
    skylines2 = 10.^(3.3+0.05*(50-ref_alt)-1.*skylam/450.) ; kR nm-1
    skylines2 = skylines2*abs(deriv(skylam))   ; kR
    skylines2 = skylines2*1.d13/(4.*!PI) ; ph  s-1 m-2 sr-1
    skylines2 = skylines2/((180d/!PI)^2d)  ; ph  s-1 m-2 arcdeg-2
    skylines2 = skylines2/(60d*60d)^2      ; ph  s-1 m-2 arcsec-2
  endif
  
loadskip:

  ;define the part of the spectrum that we care about.    
  ; The input bandpass definitions are in micron 

  inband=where(skylam/1d3 ge lam0 and skylam/1d3 le lam1) & inband = reverse(inband)
  lambda=(skylam)[inband]      ; in nm
  if not opt then $
    trans =(skytrn)[inband] $     ; Transmittance
  else $
    trans = lambda*0.0+1.0        ; If optical, than assume sky is clear
    
  ;now we can scale the spectrum to the desired flux given the V magnitude
  	  
  starfluxphotons_mag=(starfluxphotons/specnorm)*10d^(-mag*1d/2.5)

  fjy = starfluxphotons_mag*(h*c/starlam) ; ergs s-1 cm-2 nm-1
  fjy = fjy*abs(deriv(starlam)) ; erg s-1 cm-2
  fjy = fjy/abs(deriv(c/starlam)) ; erg s-1 cm-2 Hz-1
  fjy = fjy*1.d23 ; Jy

  ;interpolate the stellar model onto the wavelength grid of
  ;transmission spectrum, again keep microns and nm straight
  starflux=interpol(starfluxphotons_mag, starlam, lambda)
  fjy = mean(interpol(fjy,starlam,lambda))
  
  ; Flux as observed from observatory
  starflux_trn=starflux*trans ; Flux [photon s-1 cm-2 nm-1] 
  
  ;Now I integrate the flux density times the transmission and to get the
  ;expected photons per second from the star. Starflux has units of [photon s-1 cm-2 nm-1]
  star_photons=total(starflux_trn*abs(deriv(lambda)))  ; [photon s-1 cm-2]
  star_photons=star_photons*tel_area*tel_eff             ; [photon s-1]

  ;Sum up the thermal and sky backgrounds for the requested wavelengths
  if not opt then thermalback = total(thermalrad[inband]) else thermalback = 0.0d

  if lam0 gt 2.25 then skyflux=0. else begin
  
    inband_sky=where(lamskyline/1000. ge lam0 and lamskyline/1000. le lam1)
    skyflux=total(skylines2[inband_sky])

  endelse

  totalsky=(thermalback+skyflux)*(tel_area/1d4)*tel_eff

  ;this is now total background in photon/s/arcsec^2

  ;now we need to figure out effective integration time

  ;We assume that the PSF is a Gaussian with sigma = FWHM/2.36
  
  sigma1=(fwhm/2.36d)/pixel_size ;sigma of gaussian, in pixels, assuming that the 80% diameter ~ FWHM

  central_pix_frac=(sgauss(indgen(101)-50))[50] / ((lam1-lam0)/pixel_lam)
  
   ;this is photons/s in central pixel, assuming PSF is aligned with center of pixel
  central_pix_flux=[central_pix_frac*star_photons, totalsky*pixel_size^2d]
  
  ;this is the exposure time - starts with initial value provided
  expt1=fixed_exptime

; Loop to decrease exptime for saturation  
adaptexp:

  ;this is the number of exposures we get during flight
  nexposures=floor(flight_time/(expt1+dead_time))

  ;this is the total effective integration time
  expt2=expt1*nexposures
  signal=expt2*star_photons

  ;This is Nutzman & Charbonneau (2008) Eqn. 2 for uncertainty due to scintillation

  sigma_scint=signal*0.09*(airmass^(3/2.))/((tel_diam)^(2/3.)*sqrt(2.*expt2))*exp(-ref_alt/8.)
  
  noise1=sqrt(signal+sigma_scint^2.+phot_area*expt2*totalsky+phot_area_pix*nexposures*(read_noise/sqrt(expt1/1.4))^2.)

  precision=[signal,noise1*sqrt(1.7)]

  ;if we are non-linear within the exptime, then set the signal to negative
  ;if using adaptive exptime, decrease exptime and recalculate
  cenpix = [expt1*central_pix_flux,full_well]
  if total(cenpix[0:1]) gt full_well then begin
  	if not adapt then precision = precision*(-1.d)
  	if adapt then begin & expt1 = expt1 - 1 & goto,adaptexp & endif
  endif  	
  	
  ; Print results
  if verbose then begin
    print,'---------------------------------------------'
    print,'Spectral type: ', spt
    print,'Magnitude of the object: ', mag
    print,'Flux in Jy: ',fjy
    print,'Wavelength limits [um]: ', lam0, lam1
    print,'Neutral density filter: ', nd
    print,'---------------------------------------------'
    print,'Object signal:  ', signal
    print,'Expected noise: ', noise1
    print,'Expected precision:   ', 1.d6*noise1*sqrt(1.7)/signal
    print,'Fraction of Sat. Lim.:   ', (total(cenpix[0:1]))/full_well
    print,'---------------------------------------------'
  endif
  
end
  
