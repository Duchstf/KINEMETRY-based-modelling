;################################################################
; This is an example routine which calls KINEMETRY 
; to analyse an SDSS image of NGC4473. It makes a plot
; of position angle, ellipticity, and the 4th cosine
; coefficient (the diskiness parameter as in e.g. Bender et al. 199)
;
; Kinemetry is run given the X0,Y0 (location of the centre) and 
; initial conditions for the position angle and ellipticity. The centre is 
; fixed in this example (comment out the last keyword "fixcen" to force kinemetry 
; to look for a centre of an each ring - it will be considerably slower)
; It is also possible to run kinemetry via a grid of PA and Q (then it is not
; necessary to set up the initial conditions for PA and Q), which could
; reult in a somewhat more precise values, but it will be considerably slower).
; 
;
; This routine uses MRDFITS.PRO, RDFLOAT.PRO, DIST_CIRCLE.PRO and PLOTERROR.PRO 
; avaiable from IDL Astronomer Unser's Library.
; 
; Davor Krajnovic, Potsdam, 10.10.2013.
; V2.0, DK, Potsdam 14.07.2014. Construct an error image and pass it to kinemetry.
;###############################################################

PRO kinemetry_example_phot

;
; read in SDSS image
;
file='/Users/dkrajnov/WORK/idl/kinemetry/kinemetry_distribution_additional_stuff/NGC2974r.fits'
img=mrdfits(file, 0, h)

s=size(img)  

;
; get pixel size
;
cd1_1=fxpar(h,'CD1_1')
cd1_2=fxpar(h,'CD1_2')
cd2_1=fxpar(h,'CD2_1')
cd2_2=fxpar(h,'CD2_2')
dx = sqrt(CD1_1^2 + CD2_1^2)
dy = sqrt(CD1_2^2 + CD2_2^2)
pixelSize = 3600d*(dx + dy)/2d
scale=pixelsize

;
; set the centre (read out from the image, but note, nuclear region 
; is actually saturated on this image)
; 
X0 = 900.
Y0 = 900.
angIn=45.
epsIn=0.4

;
; measure sky level 
;
sky, img, medSkytmp, sigSkytmp, /silent
img=img-medskytmp
 
;
; select bad pixels (mask the bright star 
;
dist_circle, r, s[1:2], 1076, 779
badpixels = WHERE(r LT 100)
        
; 
; determine size and make dummy arrays
;
n=s[1]*s[2]
yy=REPLICATE(1,s[1])#(indgen(s[2]))
xx=(indgen(s[1]))#REPLICATE(1,s[2])
x=REFORM(xx, n)
y=REFORM(yy, n)
flux=REFORM(img, n)

err_img= sqrt(img)
w=where(finite(err_img, /NAN))
err_img[w]=median(img)

;
; running kinemetry
;
print, 'start kinemetry'
t=systime(1)
KINEMETRY, X, Y, flux, rad, pa, q, cf, NTRM=10, /EVEN,  /verbose, /plot, $
           X0=x0, Y0=y0, XC=xc, YC=yc,  VELCIRC=model, VELKIN=modelF, $
           ER_CF=er_cf, ER_PA=er_pa, ER_q=er_q, ER_XC=er_xc, ER_YC=er_yc,$
           XELLIP=xellip, YELLIP=yellip, IMG=img, ERROR=err_img, /bmodel, badpix=badpixels, $
           /nogrid, PAQ=[angIn-90, 1-epsIn], sky=sigskytmp, /fixcen;
print, systime(1) -t, 'seconds'
print, 'end kinemetry'
;
; This is an alternative call to kinemetry when no intinal conditions for PA and Q
; are given and the minimization is first done via a grid (as in a typical ODD case). 
;
;t=systime(1)
;KINEMETRY, X, Y, flux, rad, pa, q, cf, NTRM=8, /EVEN, NPA=21, NQ=21, /verbose, $
;           X0=x0, Y0=y0, /plot, XC=xc, YC=yc,  VELCIRC=model, VELKIN=modelF, $
;           ER_CF=er_cf, ER_PA=er_pa, ER_q=er_q, ER_XC=er_xc, ER_YC=er_yc,$
;           XELLIP=xellip, YELLIP=yellip, IMG=img, /bmodel, badpix=badpixels, $
;           sky=sigskytmp, /fixcen
;print, systime(1) -t, 'seconds'


;
; plotting image and models
;
r = GET_SCREEN_SIZE()
window,1, xsize=r[0]*0.4, ysize=r[1]*0.95
!p.charsize=2
sauron_colormap
!p.multi=[0,1,3]
        
;
; set up levels for plotting
;
lev = ROTATE(ALOG10(max(Model)) - 3.5*FINDGEN(17)/17,2)
        
        
contour, alog10(img), levels=lev[2:16], title='NGC2974', xthick=1.5, ythick=1.5, /iso, /xstyle, /ystyle
;contour, alog10(model), levels=lev, color=200, /overplot, thick=1.8
;contour, alog10(modelf), levels=lev, color=50, /overplot, thick=1.8
XYOUTS, 1900, 1500, 'MODEL (intensity on ellipses)', color=200, charsize=2.
XYOUTS, 1900, 1300, 'FULL MODEL (all terms)', color=50, charsize=2.
XYOUTS, 1900, 1100, 'ELLIPSE', color=150, charsize=2.
FOR j=0, N_elements(rad)-1, 4 DO tvellipse, rad[j], rad[j]*q[j], 900,900, pa[j]+90, COLOR=150, THICK=1.5, /data


contour, alog10(img[x0-200:x0+200, y0-200:y0+200]), levels=lev, xthick=1.5, ythick=1.5, /iso, /xstyle, /ystyle
contour, alog10(model[x0-200:x0+200, y0-200:y0+200]),  levels=lev, color=200, /overplot, thick=1.5
contour, alog10(modelf[x0-200:x0+200, y0-200:y0+200]), levels=lev, color=50, /overplot, thick=1.8

contour, alog10(img[x0-80:x0+80, y0-80:y0+80]), levels=lev, xthick=1.5, ythick=1.5, /iso, /xstyle, /ystyle
contour, alog10(model[x0-80:x0+80, y0-80:y0+80]),  levels=lev, color=200, /overplot, thick=1.5
contour, alog10(modelf[x0-80:x0+80, y0-80:y0+80]), levels=lev, color=50, /overplot, thick=1.8

;
; ploting coefficients
;
r = GET_SCREEN_SIZE()
window, 2, xsize=r[0]*0.3, ysize=r[1]*0.8
!p.charsize=3
!y.style=1
!p.multi=[0,1,4]
!Y.MARGIN=[0,0] ; show plots with shared X axis
!Y.OMARGIN=[5,3] ; allow for space for the axis labels
ploterror, rad, pa-180, er_pa, PSYM=-5, TITLE='NGC2974', xtickformat = '(A1)', YTITLE='!7C!X!N [degrees]', YRANGE=[30,70], /xlog
ploterror, rad, q, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', YTITLE='q', /xlog
ploterror, rad, cf[*,0], er_cf[*,0], PSYM=-5, xtickformat = '(A1)', YTITLE='a0 [counts]', /ylog, /xlog

a_prime=cf[*,8]*100 & b=cf[*,0] & er_a_prime=er_cf[*,8]*100 & er_b=er_cf[*,0]
er = SQRT( (1./b^2) * er_a_prime^2 + (er_b^2)*(a_prime^2)/(b^4))
;ploterror, rad, (cf[*,8]/cf[*,0])*100., er_cf[*,8], PSYM=-5, yrange=[-5,5], xtitle='RAD [arcsec]', ytitle='b4/a0 x 100'
ploterror, rad, (a_prime/b), er, PSYM=-5, yrange=[-5,5], xtitle='RAD [arcsec]', ytitle='b4/a0 x 100', /xlog
!P.MULTI=0
!Y.MARGIN=[4,2] ; back to default values
!Y.OMARGIN=[0,0]
!p.charsize=1



stop
END