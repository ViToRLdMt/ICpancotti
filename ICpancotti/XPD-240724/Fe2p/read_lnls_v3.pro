function gaussn,X,Param
    ; This gaussian funtion was written to be used with the mpfitfun procedure to fit multiple gaussian peaks using the gauss1.pro function /idl_5.2/lib
    n=size(Param,/n_elements)/3
    p=dblarr(3)
    ;evaluate n times,the gaussian function with the appropriated parameters set in Param vector.
    ;Param=[xc0,sigma0,area0,xc1,sigma1,area1,...]

    functot=0.0D
    for i=0,n-1 do begin
        for j=0,2 do begin
           p(j)=Param(3*i+j)
        endfor
        functot=functot+gauss1(X,p)
    endfor
    return,functot
end

;start main program 
;This small program was developed by A. de Siervo. Last version December 2005.

pro read_lnls_v3,filein,thetai,thetaf,dtheta,phii,phif,dphi,thetatotal1,phitotal1,intenstotal1,channeltron
lixo=' '
phi0=0
nphi=abs((phif-phii)/dphi+1)
ntheta=abs((thetaf-thetai)/dtheta+1)
thetatotal=fltarr(ntheta*nphi) & phitotal=fltarr(ntheta*nphi) & intenstotal=fltarr(ntheta*nphi)
thetatotal1=fltarr(ntheta*nphi) & phitotal1=fltarr(ntheta*nphi) & intenstotal1=fltarr(ntheta*nphi)


;Reading the XPS from LNLS data format. Do sequence is that produced in the new manipulator XXXX_theta.phi


n=0 ; constant value. It should not be changed. 
region=2 ; spectral region to be used .example : 1-Pt 4f ; 2 In 3d ; ...
reg=1 ; should start always in 1
title=' '
var=fltarr(7)

nchanel=18 ; number of points in the spectrum. Should be changed for each region 

; variable definition
;for theta=thetai,thetaf,dtheta do begin
;j=0
;if theta le 60 then nchanel=33
inten3v=fltarr(ntheta,nphi,nchanel)
i0deltaaux=fltarr(nphi,1)
i0delta=fltarr(ntheta,nphi)
izero=fltarr(ntheta)
aux=fltarr(nphi,nchanel)
smoothcurve=fltarr(nchanel)
aux0=fltarr(1,nchanel)
aa=0.000 
auxiliar0=fltarr(nphi)
auxiliar=fltarr(nphi)
background=fltarr(nchanel)
background0=fltarr(nchanel)
fittingcurve=fltarr(nchanel)
ang=fltarr(ntheta*nphi) & ind=fltarr(ntheta*nphi)
Xvar=findgen(nchanel)
f=fltarr(nphi,2)
inversef=fltarr(nphi,2)
;j=j+1
;endif
;enfor

; erro: it is used to compute error analysis in the fitting procedure. I put it to 1, which means that every point will be considered with the same weigth 

erro=fltarr(nchanel)
for i=0,nchanel-1 do begin
 erro(i)=1.0 
endfor


; reading the files
for theta=thetai,thetaf,dtheta do begin
j=0
for phi=phii,phif,dphi do begin
;for file=1,nfiles do begin
;arq=filein+strcompress(Listprefix(n),/remove_all)+'.'+strcompress(string(file),/remove_all)

if theta eq 90 then begin 
     arq=filein+'.'+strcompress(string(phi),/remove_all)
endif else begin
	arq=filein+strcompress(string(theta),/remove_all)+'.'+strcompress(string(phi),/remove_all)
endelse

;print,arq
openr,1,arq
while reg le region do begin
        if reg eq 1 then begin 
            readf,1,title
        ;    print,title
        endif
        readf,1,var
        readf,1,lixo
       ; print,var
       ; print,lixo
  for k=0,var(5)-1 do begin
     readf,1,aa
     ;print,aa
      
   if reg eq region then inten3v(n,j,k)=aa
  endfor

reg=reg+1
endwhile
close,1
reg=1
j=j+1
endfor
n=n+1
endfor

;end reading files

; the Param vector is used to give the starting point for each parameter in the curve fitting model
; each internal bracket indicates a peak with the following parameters for a gaussian shape [centroide, sigma(width), peak area]

Param=[[9.5D,2.2D,9972.0D],[23.0D,2.9D,53922.7D],[37.05D,1.85D,6000.0D]]

openw,1,'fitresults.dat' ; all fitting results will be written in this file. Remember to change the file name to save previous results. Be careful with incongluences between the starting values and the limits fixed in the constrain vectors bellow.

;Shirley background
initback=1   ; you can change this number as you like. 
endback=nchanel-1-0 ; the last indice indicate the number of channels from the end to the fixed point. Zero means the last channel in the spectrum. Change it for a better background sub. of your particular data set.  
i=0
for theta=thetai,thetaf,dtheta do begin
;   if theta ge 40 then initback=2
;   if theta gt 20 then initback=2
;   if theta ge 25 then initback=5
;   if theta eq 35 then initback=3
;   if theta eq 40 then initback=4.0
;   if theta gt 40 then initback=3.0
   j=0
   for phi=phii,phif,dphi do begin 
      ;plot,inten3v(i,j,*)
      ;if inten3v(i,j,10) eq 0 then inten3v(i,j,*)=inten3v(i,j-1,*) 
      smoothcurve(*)=inten3v(i,j,*)
      smoothcurve=smooth(smoothcurve,3) 
      inten3v(i,j,*)=smoothcurve
      ;normalizing by i0delta
      ;inten3v(i,j,*)=inten3v(i,j,*)/i0delta(i,j)
   
      a=inten3v(i,j,initback)
      b=inten3v(i,j,endback)
      background(*)=0
      background0(*)=b
      for nint=1,6 do begin
        for k2=endback,initback,-1 do begin
             a=inten3v(i,j,initback)
	     b=inten3v(i,j,endback)
             sum1=0
             sum2=0
             for k=endback,k2,-1 do begin
               sum1=sum1+inten3v(i,j,k)-background0(k)
             endfor
             for k=endback,initback,-1 do begin
               sum2=sum2+inten3v(i,j,k)-background0(k)
             endfor
             background(k2)=(a-b)*sum1/sum2+b
        endfor
      for xx=0,initback do begin
         background(xx)=background(initback)
      endfor
      for xx=endback,nchanel-1 do begin
         background(xx)=background(endback)
      endfor
      background0=background 
      ;plot,inten3v(i,j,*)
      ;oplot,background,linestyle=2
     ; wait,0.01   ; only to see the background converging

    endfor
    ;here the background is subtracted directly from the spectrum, and the array    is substituted
    plot,inten3v(i,j,*),title='THETA='+string(theta)+'  Phi='+string(phi),psym=4,symsize=1
    inten3v(i,j,*)=inten3v(i,j,*)-background(*)
    oplot,inten3v(i,j,*),psym=5,symsize=1
    oplot,background,linestyle=2
  ;  read,lixo

fitting='no'  ; use only 'yes' or 'no'

if fitting eq'yes' then begin

	constrains=replicate({fixed:0,limited:[0,0],limits:[0.D,0.D]},size(Param,/n_elements))
        
	;area constrain
	constrains(2).limits(*)=[0.0D,55000.D] & constrains(2).limited(*)=[1,1]
	constrains(5).limits(*)=[500.D,90000.D] & constrains(5).limited(*)=[1,1]
	;constrains(8).limits(*)=[0.D,2000000.D] & constrains(8).limited(*)=[1,1]
	;constrains(11).limits(*)=[0.D,2000000.D] & constrains(11).limited(*)=[1,1]
        constrains(8).fixed=1 & Param(8)=1.5*Param(2)
        ;constrains(11).fixed=1 & Param(11)=0.5*Param(5)
        ;constrains(11).limits(*)=[0.0D,55000.D] & constrains(11).limited(*)=[1,1]
        ;constrains(14).limits(*)=[0.0D,55000.D] & constrains(14).limited(*)=[1,1]

	;peak position constrain;
	;constrains(0).limits(*)=[58.7D,60.3D] & constrains(0).limited(*)=[1,1]
        ;constrains(0).fixed=1 & Param(0)=58.8D
	constrains(0).limits(*)=[9.50D,10.0D] & constrains(0).limited(*)=[1,1]
        constrains(3).limits(*)=[23.0D,25.0D] & constrains(3).limited(*)=[1,1]
        constrains(6).limits(*)=[37.05D,38.0D] & constrains(6).limited(*)=[1,1]
        ;constrains(9).fixed=1 & Param(9)=Param(0)-6.0
	;constrains(12).fixed=1 & Param(12)=Param(6)-6.0 
	
        ;sigma (peak width);
	;constrains(1).limits(*)=[3.1D,6.1D] & constrains(1).limited(*)=[1,1]
	;constrains(4).limits(*)=[1.8D,2.1D] & constrains(4).limited(*)=[1,1]
        ;constrains(1).fixed=1 & Param(1)=2.0D
        ;constrains(4).fixed=1 & Param(4)=3.1D
	;constrains(7).fixed=1 & Param(7)=1.85D
        constrains(1).limits(*)=[1.85D,2.2D] & constrains(1).limited(*)=[1,1]
        constrains(4).limits(*)=[2.5D,2.90D] & constrains(4).limited(*)=[1,1]
        constrains(7).limits(*)=[1.85D,1.9D] & constrains(7).limited(*)=[1,1]
   	;constrains(10).limits(*)=[1.85D,2.2D] & constrains(10).limited(*)=[1,1]
        ;constrains(13).limits(*)=[1.85D,2.2D] & constrains(13).limited(*)=[1,1]

	results=mpfitfun('gaussn',Xvar,inten3v(i,j,*),erro,Param,PARINFO=constrains)
   
   for k=0, size(results,/n_elements)-1 do begin
 	inten3v(i,j,k)=results(k)
        
   endfor

   printf,1,theta,phi,results
   ;read,lixo
   ; plotting the fitting results
   oplot,gaussn(Xvar,results) ; total
   p=dblarr(3)
   for pp=0,(size(Param,/n_elements)/3)-1 do begin
       for ii=0,2 do begin
        p(ii)=Param(3*pp+ii)
       endfor
    oplot,gaussn(Xvar,p),linestyle=1; each component
   endfor
   Param=results
   ;wait,0.01
   read,lixo
endif ;end procedure fit

    j=j+1
   endfor
 Print,'Theta=',string(theta)
 i=i+1
endfor
close,1  ; close file fitresults.dat

;end shirley background and non-linear fit

;transforme from vector format to holoxwin format and calculate izero curves.
i=0
t=0
for theta=thetai,thetaf,dtheta do begin
   j=0
   izero(i)=0
      
   openw,2,'result_evap.txt' 
   for phi=phii,phif,dphi do begin
      intenstotal(t)=0
      k=channeltron
      if fitting ne 'yes' then begin
      	for k=0,nchanel-1 do begin
                ; if peak wasnt fitted, the intensity will be the integrated area.
        	intenstotal(t)=intenstotal(t)+inten3v(i,j,k);/-inten3v(i,j,nchanel-1)
      	endfor
      endif else begin 
        	intenstotal(t)=intenstotal(t)+inten3v(i,j,k);/-inten3v(i,j,nchanel-1)
      endelse
      printf,2,phi,intenstotal(t)
     
     ;plot,inten3v(i,j,*)
     ;wait,0.00
     inten3v(i,j,channeltron)=intenstotal(t)

     thetatotal(t)=theta*!pi/180. & phitotal(t)=(phi0+phi)*!pi/180.; & intenstotal(t)=(inten3v(i,j,channeltron)-0*inten3v(i,j,nchanel-1));/inten3v(i,j,nchanel-1)
      if (phitotal(t) gt 2.0*!pi) then begin
          phitotal(t)=phitotal(t)-2.0*!pi
      endif
      if (phitotal(t) lt 0.0) then begin
         phitotal(t)=phitotal(t)+2.0*!pi
      endif
      ; next 3 lines used for sort algorithm
      ind(t)=t
      ang(t)=10000*thetatotal(t)*180./!pi+phitotal(t)*180./!pi 

; Preparing to calculate izero using average of data
      izero(i)=izero(i)+intenstotal(t)
      j=j+1
      t=t+1
   endfor
   close,2
; smoothing the data
; or smoothing strongly for some angles
     smooth=4
     ;if theta eq 2 then wsmooth=6
     ;if theta eq 24 then wsmooth=6
     ;if theta eq 34 then wsmooth=6
     ;if theta eq 42 then wsmooth=6
     ;if theta eq 46 then wsmooth=6

     auxiliar=inten3v(i,*,channeltron)
     auxiliar=smooth(auxiliar,smooth,/NAN)

     ;smothdata=inten3v(i,*,channeltron)
     ;calculate fit_poly background
     x=findgen(nphi)
     results=poly_fit(x,auxiliar,3)
     print,results
     ;plot,i0delta(i,*)
     plot,inten3v(i,*,channeltron)
     oplot,auxiliar
     ;read,lixo

     ;substituting the original data by the smooth one.
     inten3v(i,*,channeltron)=auxiliar

     for kk=0,nphi-1 do begin
        auxiliar(kk)=results(0)+results(1)*x(kk)+results(2)*x(kk)^2+results(3)*x(kk)^3 ;+results(4)*x(kk)^4;results(5)*x(kk)^5
  ; calculationg chi using the fit polynominal curve
     inten3v(i,kk,channeltron)=(inten3v(i,kk,channeltron)-1.0*auxiliar(kk))/auxiliar(kk)
     endfor

     oplot,auxiliar
     wait,0.1 
     izero(i)=izero(i)/nphi
   i=i+1
endfor
  plot,izero,title='Polar scan'
  ;if you want the polar curve such a file, uncomment the following lines
 ; openw,1,'polar_scan.dat'
 ; printf,1,izero(*)
 ; close,1
  read,lixo
;calcule chi function using the average of data
;or use inten3v where the chi function (anizotropy) is already calculated.
i=0
t=0
for theta=thetai,thetaf,dtheta do begin
   j=0
   print,theta,izero(i)
   for phi=phii,phif,dphi do begin
        ;intenstotal(t)=(intenstotal(t)-izero(i))/izero(i)
        intenstotal(t)=inten3v(i,j,channeltron)
      j=j+1
      t=t+1
   endfor
   i=i+1
endfor
; end chi function

; simetrization
simetria='false'
i=0
t=0

; Fourier expansion --> This simetrization use Fourier exp. to take simetric part of the whole XPD data in each azimuth curve
t=0
t1=0
for theta=thetai,thetaf,dtheta do begin
        j=0
	for phi=phii,phif,dphi do begin 
		f(j,0)=intenstotal(t)
                f(j,1)=0
		t=t+1
		j=j+1	
        endfor
        plot, f(*,0),linestyle=1,title='THETA='+string(theta)
; calculate F(u)=FFT(f(x))
    for u=0,nphi-1 do begin
    	sumR=0
    	sumI=0
    	for x=0,nphi-1 do begin
       		SumR=sumR+f(x,0)*cos(2*!pi*u*x/nphi)/nphi
       		SumI=sumI-f(x,0)*sin(2*!pi*u*x/nphi)/nphi
    	endfor
    	inversef(u,0)=SumR
    	inversef(u,1)=SumI
     endfor
;taking only the c4 simetric componets

    for x=0,nphi-1 do begin
    	sumR=0
    	sumI=0
        for u=2,nphi-1,2 do begin 
      
	;for u=0,nphi-1,1 do begin  ; use this to see original data (all components of Fourien expansion.)
       		SumR=sumR+inversef(u,0)*cos(2*!pi*u*x/nphi)-inversef(u,1)*sin(2*!pi*u*x/nphi)
       		SumI=sumI+inversef(u,0)*sin(2*!pi*u*x/nphi)+inversef(u,1)*cos(2*!pi*u*x/nphi)
    	endfor
	f(x,0)=SumR
	f(x,1)=SumI
     endfor
        oplot, f(*,0),linestyle=2
	read, lixo
        wait,0.5
;retorning value to intenstotal vector
        j=0
	for phi=phii,phif,dphi do begin 
		intenstotal(t1)=f(j,0)
		t1=t1+1
		j=j+1	
        endfor
endfor

;end Fourier symmetrization

;Rotina auxiliar para ordenar
thetatotal1(*)=thetatotal(*)
phitotal1(*)=phitotal(*)
intenstotal1(*)=intenstotal(*)

max=0.0
min=9999999999.0
n=ntheta*nphi
for i=0,n-1 do begin
  for j=i,n-1 do begin
    if (ang(j) lt ang(i)) then begin
        min=ang(j)
        x=ang(i)
        ang(i)=ang(j)
        ang(j)=x
        x=ind(i)
        ind(i)=ind(j)
        ind(j)=x
     endif
   endfor
 endfor

for i=0,ntheta*nphi-1 do begin
      thetatotal(i)=thetatotal1(ind(i))
      phitotal(i)=phitotal1(ind(i))
      intenstotal(i)=intenstotal1(ind(i))
endfor
; end sort algorithm


mscdinp='yes'
If mscdinp eq 'yes' then begin
;preparing mscdinput
kphoton=19.59  ; change it for the correct k value each time
thetaerro=0.2  ; correcting experimental polar angle constant.
; Open a new file to write
OPENW, 1, 'exp_Fe2_GaO.out'
printf,1,"   323    17    0     datakind beginning-row linenumbers"
printf,1,'----------------------------------------------------------------'
printf,1,'MSCD Version 1.00 Yufeng Chen and Michel A Van Hove'
printf,1,'Lawrence Berkeley National Laboratory (LBNL), Berkeley, CA 94720'
printf,1,'Copyright (c) Van Hove Group 1997. All rights reserved'
printf,1,'--------------------------------------------------------------'
printf,1,' angle-resolved photoemission extended fine structure (ARPEFS)'
printf,1,' experimental data for Pd 3d3/2 from W(100)  excited with hv=1810eV'
printf,1,' '
printf,1,' provided by Pancotti et al. (LNLS in 12, June 2007) '
printf,1,'   intial angular momentum (l) = 2'
printf,1,'   photon polarization angle (polar,azimuth) = (  30.0,   0.0 ) (deg)'
printf,1,'   sample temperature = 300 K'
printf,1,' '
printf,1,'   photoemission angular scan curves'
printf,1,'     (curve point theta phi weightc weighte//k intensity chiexp)'
printf,1,ntheta, ntheta*nphi,1,ntheta,nphi,ntheta*nphi

i=0
t=0
for theta=thetai,thetaf,dtheta do begin
   j=0
   printf,1,i+1,nphi,kphoton,theta+thetaerro,1.0,0.0
   for phi=phii,phif,dphi do begin
       printf,1,phitotal(t)*180./!pi,intenstotal(t)*izero(i)+izero(i),izero(i),intenstotal(t)
      j=j+1
      t=t+1
   endfor
   i=i+1
endfor
; Close file unit 1:
CLOSE, 1
endif

;n-fold symmetry
n=1
thetatotal1=fltarr(ntheta*nphi*n) & phitotal1=fltarr(ntheta*nphi*n) & intenstotal1=fltarr(ntheta*nphi*n)
for i=0,nphi*ntheta-1 do begin
     for nn=0,n-1 do begin
      thetatotal1(i+nn*nphi*ntheta)=thetatotal(i)
      phitotal1(i+nn*nphi*ntheta)=phitotal(i)+(360.0*nn/n)*!pi/180.
      intenstotal1(i+nn*nphi*ntheta)=intenstotal(i)
     endfor
 endfor

return
end
