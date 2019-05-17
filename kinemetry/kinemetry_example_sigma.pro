;################################################################
; This is an example routine which calls KINEMETRY 
; to analyse SAURON velocity dispersion (sigma) map of NGC2974 
; (Emsellem et al. 2004 MNRAS, 352, 271). It makes a plot
; of kinemetric coefficients: position angle and flattening
; of the ellipse, sigma, and 4th cosine coefficient. 
; 
; Note that for this example the range of PA and Q was limited
; to values expected from the photometry of NGC2974. An 
; unconstrained fit is not well constrained with respect to the
; PA and Q of the best fitting ellipse. In general, it is advised
; to run kinemetry on velocity diseprsion maps using either a 
; constant PA and Q (for example, the average values from the 
; photometry) or run kinemetry on set of rings with different PA 
; and Q, which were previously obtained by running kinemetry on 
; an image. In this latter case, it is also advised to run kinemetry
; on the iamge created from the data cube itself. 
; 
; This routine uses RDFLOAT.PRO and PLOTERROR.PRO avaiable from
; IDL Astronomer Unser's Library.
; 
; V1.0. 10.10.2013, Davor Krajnovic, Potsdam, 
;###############################################################

PRO kinemetry_example_sigma

;
; read in all data
;
file = '/Users/dkrajnov/WORK/idl/kinemetry/kinemetry_distribution_additional_stuff/NGC2974_SAURON_kinematics.dat'

rdfloat, file, num, xbin, ybin, velbin, er_velbin, sigbin, er_sigbin, SKIPLIN=1

;
; kinemetry on velocity map
;
t=systime(1)
KINEMETRY, xbin, ybin, sigbin, rad, pa, q, cf, ntrm=8, scale=0.8, $
  ERROR=er_sigbin, name='NGC2974',er_cf=er_cf, er_pa=er_pa, $
  er_q=er_q, /EVEN, /plot, /verbose, RANGEPA=[40,50]-90, RANGEQ=[0.5,0.99]
print, systime(1) -t, 'seconds'


;
; plot coeffs.
;
r = GET_SCREEN_SIZE()
window, 1, xsize=r[0]*0.3, ysize=r[1]*0.8
!p.charsize=3
!y.style=1
!p.multi=[0,1,4]
!Y.MARGIN=[0,0] ; show plots with shared X axis
!Y.OMARGIN=[5,3] ; allow for space for the axis labels
ploterror, rad, pa-180, er_pa, PSYM=-5, TITLE='NGC2974 sigma', xtickformat = '(A1)', YTITLE='!7C!X!N [degrees]', YRANGE=[30,70]
ploterror, rad, q, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', YTITLE='q'
ploterror, rad, cf[*,0], er_cf[*,0], PSYM=-5, xtickformat = '(A1)', YTITLE='a0 [km/s]', YRANGE=[50,245]
ploterror, rad, (cf[*,8]/cf[*,0])*100., er_cf[*,8], PSYM=-5, yrange=[-5,5], xtitle='RAD [arcsec]', ytitle='b4/a0 x 100'
oplot, [0,30], [0,0], color=0, linestyle=2
!P.MULTI=0
!Y.MARGIN=[4,2] ; back to default values
!Y.OMARGIN=[0,0]
!p.charsize=1

stop

END
