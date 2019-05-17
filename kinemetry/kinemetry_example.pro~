;################################################################
; This is an example routine which calls KINEMETRY 
; to analyse SAURON velocity map of NGC2974 as presented in
; (Emsellem et al. 2004 MNRAS, 352, 271). It makes a plot
; of kinemetric coefficients as in the first column of Fig.7 
; in Krajnovic et al. (2006).
;
; This routine uses RDFLOAT.PRO and PLOTERROR.PRO avaiable from
; IDL Astronomer Unser's Library.
; 
; Davor Krajnovic, Oxford, 07.12.2005.
;###############################################################

PRO kinemetry_example

;
; read in all data
;
file = 'NGC2974_SAURON_kinematics.dat'

rdfloat, file, num, xbin, ybin, velbin, er_velbin,  SKIPLIN=1

;
; kinemetry on velocity map
;
t=systime(1)
KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, ntrm=6, scale=0.8, $
	ERROR=er_velbin, name='NGC2974',er_cf=er_cf, er_pa=er_pa, $
	er_q=er_q, /plot, /verbose
print, systime(1) -t, 'seconds'

;
; kinemetry parameters as defined in Krajnovic et al. (2006)
; 
k0 = cf[*,0]
k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
k5 = SQRT(cf[*,5]^2 + cf[*,6]^2)
k51 = k5/k1
erk1 = (SQRT( (cf[*,1]*er_cf[*,1])^2 + (cf[*,2]*er_cf[*,2])^2 ))/k1
erk5 = (SQRT( (cf[*,5]*er_cf[*,5])^2 + (cf[*,6]*er_cf[*,6])^2 ))/k5
erk51 = ( SQRT( ((k5/k1) * erk1)^2 + erk5^2  ) )/k1 

;
; plot coeffs.
;
r = GET_SCREEN_SIZE()
window, 1, xsize=r[0]*0.3, ysize=r[1]*0.8
!y.style=1
!p.multi=[0,1,4]
!Y.MARGIN=[0,0] ; show plots with shared X axis
!Y.OMARGIN=[5,3] ; allow for space for the axis labels
ploterror, rad, pa, er_pa, PSYM=-5, TITLE='NGC2974', xtickformat = '(A1)', YTITLE='!7C!X!N [degrees]', YRANGE=[30,70]
ploterror, rad, q, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', YTITLE='q'
ploterror, rad, k1, erk1, PSYM=-5, xtickformat = '(A1)', YTITLE='k1 [km/s]',YRANGE=[0,245]
ploterror, rad, k51, erk51, PSYM=-5, XTITLE='R [arcsec]', YTITLE='k5/k1', YRANGE=[0,0.13]
!P.MULTI=0
!Y.MARGIN=[4,2] ; back to default values
!Y.OMARGIN=[0,0]

END
