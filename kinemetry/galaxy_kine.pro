PRO galaxy_model_kine, galaxy
  ;Parameter: - name of the xdr file in string
  ;
  ;Action: Fit a model in KINEMETRY and write the returned values in a fits file (Specific details about returned values are in KINEMERTY's documentation)
  ; - The fits file parameters are described in detail in the python file "plotting_model" which analysses KINEMETRY's output in details. 
  
  restore, galaxy 
  
  ;process the velocity map + generate inputs for kenemetry. Save the error for later analysis
  velocity_map = windstr.e_vel['v%50c1','Halpha'] ;Change the second entry depending on whatever emission line the user want to use
  error = windstr.e_vel["v%50c1err","Halpha"]
  
  ; Take the center's coordinate
  windstr.center_nuclei -= 1 ;The coordinates in data cube is 1 larger than the actual index in the velocity map
  x_center = windstr.center_nuclei[0] ; might be reversed, not sure how array's indexes work
  y_center = windstr.center_nuclei[1]
  
  ; Take the size
  s=size(velocity_map)
  
  ;generate the inputs 
  n=s[1]*s[2]
  yy=REPLICATE(1,s[1])#(indgen(s[2]))
  xx=(indgen(s[1]))#REPLICATE(1,s[2])
  xbin=REFORM(xx, n)
  ybin=REFORM(yy, n)
  
  ;moment input is the velocity flattened to 1D array
  velbin = REFORM(velocity_map, n)
  
  ;Array to save good pixels
  xgood = []
  ygood = []
  velgood = []
  badpix = []
  
  ;Loop through the velbin to select the indexes of bad pixels
  for x = 0, n - 1 do begin
    if velbin[x] lt  500000000 then begin
      xgood = [xgood,xbin[x]]
      ygood = [ygood,ybin[x]]
      velgood = [velgood, velbin[x]]
    endif
    
    ;if the data is bad then save the index
  ;  if velbin[x] gt  500000000 then begin
  ;    badpix = [badpix,x]
  ;  endif 
  endfor
  
  ; Loop over each row to ignore ignore bad values
  for i  = 0, s[1]-1 do begin
    for j = 0, s[2]-1 do begin
      if velocity_map[i,j] gt  500000000 then begin ; ignore particularly big values
        velocity_map[i,j] = 0
      endif
      if error[i,j] eq 0 then begin
        error[i,j] = 10000000 ;So that to avoid dividing by zero in signal/noise analysis later 
      endif
    endfor
  endfor
  
  ;test output variables
  vel_model = []
  xellip = []
  yellip = []
  
  ;Adjust the center
  xgood -= x_center
  ygood -= y_center
  
  ;Call KINEMETRY
  KINEMETRY, xgood,ygood,velgood,rad,pa,q,cf,VELCIRC=vel_model,$
    /plot, /verbose, /bmodel
    
  ;Calculate the residual
  residual = velgood - vel_model
  
  ;clear bad data in vel_model and clear the values outside the range of model in residual
  si_vel = size(vel_model);
  for i = 0, si_vel[1]-1 do begin
      if vel_model[i] gt  500000 then begin ; ignore particularly big values
        vel_model[i] = 0
        residual[i] = 0 
       endif    
  endfor
  
  ;Add the center back to coordinate
  xgood += x_center
  ygood += y_center
  
  ;Create final arrays to display
  final_vel_mol= MAKE_ARRAY(s[1], s[2], /INTEGER, VALUE = 0)
  final_vel_residual= MAKE_ARRAY(s[1], s[2], /INTEGER, VALUE = 0)
  ;Populate the arrays 
  for i=0,si_vel[1]-1 do begin
    final_vel_mol[xgood[i],ygood[i]]= vel_model[i]
    final_vel_residual[xgood[i],ygood[i]]= residual[i] 
  endfor
  ;Calculate the signal to noise
  
  signal_to_noise = final_vel_residual/error
  
  ;Save the center's coordinate into fits as well (the fits file will be passed into python, which USE ROW, COL INDEXING INSTEAD OF USING COL,ROW LIKE IDL)
  center = [y_center,x_center] ; Flip the coordinate over for Python
  
  ;write rhe residual, model, original velocity to an fits file (this will later be processed in python)
  writefits, "velModel.fits", final_vel_mol, "Model"
  writefits, "velModel.fits", final_vel_residual, "Residual", /append
  writefits, "velModel.fits", velocity_map, "velOG",/append
  writefits, "velModel.fits", center, "Center Coordinate",/append
  writefits, "velModel.fits", rad, "Radii of best fitted ellipses",/append ; rad is still in spaniels
  writefits, "velModel.fits", q, "(1 - Eliipticity)",/append 
  writefits, "velModel.fits", pa, "Potisiotn Angle",/append
  writefits, "velModel.fits", signal_to_noise, "Signal/Noise",/append
  
  ;extract output parameteres
  k0 = cf[*,0]
  k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
  k5 = SQRT(cf[*,5]^2 + cf[*,6]^2)
  k51 = k5/k1
  
  ;change spaxel to kpc
  rad = (rad/0.8)*1.3
  
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
  plot, rad, k1, PSYM=-5, XTITLE='R [kpc]', YTITLE='Rotational Velocity [km/s]', YRANGE=[0,245]
  !P.MULTI=0
  !Y.MARGIN=[4,2] ; back to default values
  !Y.OMARGIN=[0,0]
  !p.charsize=1
  
END 