pro holoxwin, option,channeltron

;This program plots holograms resulting from unpolarized light experiments
;by adding the calculated intensities for the two components. This program
;is only useful if the hologram has 4 and 3 fold symetry and the hologram of each
;comnponent was calculated with theta and phi between 0 and 90 degrees for
;4-fold and between 0 and 120 degrees for p3. 
;The program automatically generates the whole hologram from symmetry 
;considerations.
;
;arqin1 = first hologram file 
;arqin2 = second hologram file
;option = define which kind of projection to be use in the plot
;       = 'stereo' -> x=2*tan(theta/2.)*cos(phi)
;       = 'theta'  -> x=(theta)*cos(phi)
;       = 'k'      -> x=k*sin(theta)*cos(phi)

;Device, pseudo_color=16
;Device, Get_Visual_Depth=thisDepth
;IF thisDepth GT 8 THEN Device, Decomposed=0



Device, true=24
Device, Get_Visual_Depth=thisDepth
IF thisDepth GT 8 THEN Device, Decomposed=0, retain=2



;IF thisDepth GT 8 THEN Device, 
;Device, Set_font = '/TIMES, /BOLD '  


;Device,true_color=16,
 surfdata = shift(dist(20),10,10)
  surfdata = bytscl(exp(-(surfdata/5)^2))

loadct, 8  ;maiores informacoes    (http://ham.space.umn.edu/johnd/ct/ct-names.html) 
    shade_surf, surfdata, shades=surfdata
    colour_surf = tvrd()



;read_lnls_v3,'DZ10_',18,72,3,3,120,3,thetatotal,phitotal,intensity,channeltron

;read_lnls_v3,'/home/fbernardi/Alex/data/DZ1208/Pd3d/DZ12_',30,72,3,0,120,3,thetatotal,phitotal,intensity,channeltron

read_LNLS_v3,'JL24_-',12,69,3,0,357,3,thetatotal,phitotal,intensity,channeltron


max_theta=max(thetatotal,min=min_theta)
print,'theta minimo=',min_theta
print,'theta maximo=',max_theta

pi2=!pi/2. & density=600

if (option eq 'stereo') then begin
 x=2*tan(thetatotal/2.)*cos(phitotal) & y=2*tan(thetatotal/2.)*sin(phitotal)
endif

if (option eq 'theta') then begin
 ;correction for 'const fild'
 ;phitotal=phitotal+thetatotal*sin(17*!pi/180)
 ;
 x=(thetatotal)*cos(phitotal) & y=(thetatotal)*sin(phitotal) 
endif

if (option eq 'k') then begin
 k=10.80
 x=k*sin(thetatotal)*cos(phitotal) & y=k*sin(thetatotal)*sin(phitotal)
endif

triangulate,x,y,triang
result=trigrid(x,y,intensity,triang,nx=density,ny=density)
final=smooth(result,15)


a=60./(max(final)-min(final)) & b=(max(final)-min(final))/59. - min(final)

ll=findgen(60)/a - b

xtt=findgen(density)/(density/(max(x)-min(x))) + min(x)
ytt=findgen(density)/(density/(max(y)-min(y))) + min(y)
window, 1, xsize=800,ysize=780,retain=2
contour, final, xtt,ytt, $
         levels=ll,$
         xrange=[-1.8,1.9], xst=1,$
         yrange=[-1.9,1.8], yst=1,$
         /fill
phitt=findgen(101)/100.*2.*!pi
xfac=(4.*!pi/180.0)+max_theta

;plotando a parte esterna do holograma
for i=0,1200 do begin
 k=float(i)
 plots,(k/100.+max_theta+0.5*!pi/180)*cos(phitt),(k/100.+max_theta+0.5*!pi/180)*sin(phitt),thick=4
endfor
print,min_theta
for i=0,100 do begin
 k=float(i)
 plots,(min_theta-2*!pi/180-k/100.*min_theta)*cos(phitt),(min_theta-2*!pi/180-k/100.*min_theta)*sin(phitt),thick=4
endfor
plots,max_theta*cos(phitt),max_theta*sin(phitt)

plots,xfac*cos(phitt),xfac*sin(phitt),color=1,thick=3; anel ao redor do padrao de difracao
;yfac=1.1
yfac=1.5; ; posicao vertical do eixo theta

plots,[-xfac,-xfac],[-yfac,0],color=1,thick=3; risco indo centro baixo lado esquerdo
plots,[xfac,xfac],[-yfac,0],color=1,thick=3; risco indo centro baixo lado direito
plots,[-xfac,xfac],[-yfac,-yfac],color=1,thick=3; risco da graduacao theta

print,max_theta
step=max_theta/3.
for i=0,3 do begin
 x=i*step
 print,x
 plots,[x,x],[-yfac,-yfac-0.1],color=1,thick=3; risco meio graduacao theta 
 plots,[-x,-x],[-yfac,-yfac-0.1],color=1,thick=3
 strx=strcompress(fix(x*180.0/!pi),/remove_all)
 if i eq 0 then begin
  xyouts,x-0.02,-yfac-0.22,strx,color=1,chars=3.0, CHARTHICK=2
 endif else begin
  xyouts,x-0.05,-yfac-0.22,strx,color=1,chars=3.0, CHARTHICK=2
  xyouts,-x-0.05,-yfac-0.22,strx,color=1,chars=3.0, CHARTHICK=2
 endelse
endfor
;xyouts,-0.1,-yfac-0.4,'!4h',color=1,chars=2.5, CHARTHICK=2; escreve o simbolo theta

plots,[xfac,xfac-0.1],[0,0],color=1,thick=3 ; risco do phi = 0
plots,[0,0],[xfac,xfac-0.1],color=1,thick=3 ; risco do phi = 0 
xyouts,xfac+0.06,0,'0',color=1,chars=3, CHARTHICK=2 ; escreve o 0 do phi
xyouts,-0.05,xfac+0.06,'90',color=1,chars=3,CHARTHICK=2 ; escreve o 90 do phi
 

;xyouts,xfac*cos(45.*!pi/180.)+0.1,xfac*sin(45.*!pi/180),'!4u',color=1,chars=3,CHARTHICK=3; escreve o phi

;printing anysotropy
Print,'max anisotropy=',max(intensity)
Print,'min anisotropy=',min(intensity)
anisotropy=' '
anisotropy=strcompress(string(round(100*(max(intensity)+abs(min(intensity))))),/remove_all)+'%'
xyouts,-yfac+0.3,-xfac-0.1,anisotropy,color=1,chars=3


;RA='!18R!Ia!N='+strmid(strcompress(string(rafactor),/remove_all),0,5)

;Device, Set_font = '/TIMES, /BOLD'
 
;xyouts, yfac-1.1,-xfac-0.11,RA, color=1,chars=3.5, CHARTHICK=2; escreve o Ra factor

xyouts,-0.06,-yfac-0.47,'!4h',color=1,chars=5,CHARTHICK=3 ;escreve o THETA
xyouts,xfac*cos(45.*!pi/180.)+0.1,xfac*sin(45.*!pi/180),'!4u',color=1,chars=5, CHARTHICK=3 ;escreve o PHI


end
