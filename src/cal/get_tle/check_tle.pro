;
fname='check_tle.dat'
get_lun,unit0
openr,unit0,fname
day=fltarr(1000)
apogee=fltarr(1000)
perigee=fltarr(1000)
semiax=fltarr(1000)
raan=fltarr(1000)
incl=fltarr(1000)
mean_anom=fltarr(1000)
arg_perig=fltarr(1000)
eccen=fltarr(1000)
d=0.0
a=0.0
p=0.0
sem=0.0
ra=0.0
inc=0.0
ma=0.0
ap=0.0
ec=0.0
j=0
while (not EOF(unit0)) do begin
   readf,unit0,d,a,p,sem,ra,inc,ma,ap,ec
   day[j]=d
   apogee[j]=a
   perigee[j]=p
   semiax[j]=sem
   raan[j]=ra
   incl[j]=inc
   mean_anom[j]=ma
   arg_perig[j]=ap
   eccen[j]=ec
   j=j+1
endwhile
day=day(0:j-1)
apogee=apogee(0:j-1)
perigee=perigee(0:j-1)
semiax=semiax(0:j-1)
raan=raan(0:j-1)
incl=incl(0:j-1)
mean_anom=mean_anom(0:j-1)
arg_perig=arg_perig(0:j-1)
eccen=eccen(0:j-1)
close,unit0
free_lun,unit0
!x.title='Day of Year 1999'
!y.title='Arbitrary units'
!p.title='FUSE Orbital Elements'
!x.style=1
!y.range=[0,1000]
!x.range=[170,day(0)+10]
!p.charsize=1.7
!p.linestyle=0
plot,day,apogee
oplot,day,perigee
oplot,day,semiax*100.0-714470.0+700.0
oplot,day,raan*2.0
oplot,day,incl*3500.0-87445.0+900.0
oplot,day,mean_anom*2.0
oplot,day,arg_perig*2.0
oplot,day,eccen/2E-6
xyouts,day(0),apogee(0),"Apogee",charsize=1.0
xyouts,day(0),perigee(0),"Perigee",charsize=1.0
xyouts,day(0),semiax(0)*100.0-714470.0+700.0,"Semimajor axis*100",charsize=1.0
xyouts,day(0),raan(0)*2.0,"RAAN*2",charsize=1.0
xyouts,day(0),incl(0)*3500.0-87445.0+900.0,"incl*3500",charsize=1.0
xyouts,day(0),mean_anom(0)*2.0,"mean anomaly*2",charsize=1.0
xyouts,day(0),arg_perig(0)*2.0,"argument of perigee*2",charsize=1.0
xyouts,day(0),eccen(0)/2E-6,"eccen/2E-6",charsize=1.0

end
