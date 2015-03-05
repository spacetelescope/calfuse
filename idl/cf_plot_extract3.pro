;+
;  cf_plot_extract3.pro is a procedure to plot the extraction windows over a 
;                   greyscale image for a given IDF to a standard image file. 
;
;  Author:  Edward M. Murphy
;  Written: 1999 November 24
;
;       Calling sequence: cf_plot_extract3, rootname
;
;       Inputs:    rootname will have 'idf.fit' added to it.
;
; 02/12/03 jch: Adapt from cf_plot_extract.pro to work for CalFUSE 3.0.
; 03/10/03 wvd: Silence calls to mrdfits.  Simplify image construction.
; 04/08/03 wvd: Shift apertures in X to match FPA positions.
; 06/02/03 wvd: Don't plot non-target apertures in HIST mode.
; 10/10/03 wvd: Exclude photons from bad times and those outside PHA limits.
; 12/09/03 wvd: Use extended apertures when appropriate.
; 02/26/04 wvd: Use strcmp to compare strings.
; 05/03/04 wvd: In HIST mode, don't plot non-target apertures.
;		Read expected spectral centroid from CHID_CAL file header.
; 05/14/04 wvd: Don't use strcmp to compare strings.
;		Modify boundaries between LiF and SiC channels.
; 06/14/04 wvd: Don't shift apertures in X to match FPA positions.
; 07/12/05 wvd: Clean up call to HIST_2D.
; 08/01/05 wvd: If IDL version > 6.0, generate GIF file.
; 08/30/05 wvd: If there are no good photons, plot all photons.
; 11/01/05 wvd: If DAYNIGHT = NIGHT, exclude daytime photons.
; 03/14/08 wvd: Require LOC_FLAG < 16B for good data.
;-

pro cf_plot_extract3,rootname

	!quiet=1
	filename=rootname+'idf.fit'
	a=mrdfits(filename,0,ah,/SILENT)

	; Read input file's primary header keywords
	wave_cal = strcompress(fxpar(ah,'WAVE_CAL'),/rem)
	chid_cal = strcompress(fxpar(ah,'CHID_CAL'),/rem)
	filedate = strcompress(fxpar(ah,'DATE'),/rem)
	obsdate = strcompress(fxpar(ah,'DATEOBS'),/rem)
	obstime = strcompress(fxpar(ah,'TIMEOBS'),/rem)
	cfvers = strcompress(fxpar(ah,'CF_VERS'),/rem)
	exptime = fxpar(ah,'EXPTIME')
	instmode = strcompress(fxpar(ah,'INSTMODE'),/rem)
	detector = strcompress(fxpar(ah,'DETECTOR'),/rem)
	aperture = strcompress(fxpar(ah,'APERTURE'),/rem)
	srctype = strmid(fxpar(ah,'SRC_TYPE'),0,1)
	daynight = strcompress(fxpar(ah,'DAYNIGHT'),/rem)

        ; Interpret aperture and srctype keywords.
        if (aperture EQ 'MDRS') then begin
		aplif = 2
		apsic = 6
        endif else if (aperture EQ 'HIRS') then begin
		aplif = 1
		apsic = 5
        endif else begin	; assume LWRS
		aplif = 3
		apsic = 7
        endelse
	if (srctype EQ 'P') then extend = 0 else extend = 8

	cal_path = getenv("CF_CALDIR")
	slitname=["HIRS","MDRS","LWRS","PINH","HIRS","MDRS","LWRS","PINH"]
	specname=["LiF","LiF","LiF","LiF","SiC","SiC","SiC","SiC"]    

	; Set the plot resolution.
	set_plot,'Z'
	device,set_resolution=[1124,612]

	; Define plot attributes
	!x.style=1
	!y.style=1
	!x.tickformat='(I5)'
	!p.title=filename
	!x.title="X pixel"
	!y.title="Y pixel"
	!p.charsize=1.2
	cs=0.8

        ; Read the X and Y columns from the first extension (photon list),
        ; exclude bad events, and generate an image from the data.
        ftab_ext, filename, 'X,Y,TIMEFLGS,LOC_FLGS', xarr, yarr, timeflag, loc_flag
	if (daynight EQ 'NIGHT') then begin
            good = where (timeflag eq 0 and loc_flag lt 16, n) 
	endif else begin
            good = where ((timeflag and not 1B) eq 0 and loc_flag lt 16, n) 
	endelse
        if (n gt 0) then begin
            xarr = xarr[good]
            yarr = yarr[good]
        endif
        a = HIST_2D(xarr,yarr,bin1=16,bin2=2,min1=0,max1=16383,min2=0,max2=1023)
        b = HIST_EQUAL(a)

	!p.charsize=1.3
	x0=85
	y0=65
	x1=1024+x0
	y1=512+y0
	dummy_x = indgen(2)*16383
	dummy_y = indgen(2)*1023

	; Plot the pixel axes and the photon image.
	plot,dummy_x,dummy_y,/nodata,position=[x0,y0,x1,y1],/DEVICE,$
		background=255,color=0,xticks=8
	tv,255-b,x0,y0
	plot,dummy_x,dummy_y,/nodata,/noerase,position=[x0,y0,x1,y1],/DEVICE,$
		background=255,color=0,xticks=8,xminor=4

	miny=10000
	maxy=0.0

	; Plot the extraction windows.
	for i=1,7 do begin
	    ; Skip the pinhole and any non-target apertures for HIST data.
	    if ((i ne 4) and ((instmode EQ 'TTAG') $
                or (i eq aplif) or (i eq apsic))) then begin

		; Read aperture limits from calibration file.
		j = 8; use extended apertures for non-target apertures
		if ((i eq aplif) or (i eq apsic)) then j = extend
		xw=mrdfits(cal_path+'/'+chid_cal,i+j,chidhdr,/SILENT)
		yhigh = xw.yhigh
		ylow  = xw.ylow 
		default_y_centroid = fxpar(chidhdr,'CENTROID')

		; Get the YCENT keywords and adjust the windows accordingly.
		key = 'YCENT' + strcompress(string(i),/rem)
		ycent = fxpar(ah, key)
		shift = nint(default_y_centroid - ycent)
		yhigh = yhigh - shift
		ylow  = ylow - shift

		; Plot the extraction windows.
		oplot,yhigh,color=0
		oplot,ylow,color=0
		miny=min([miny,ylow])
		maxy=max([maxy,yhigh])
		xyouts,15600,ylow[16350]+5,specname[i-1],$
			color=0,charsize=0.8
		xyouts,15600,ylow[16350]-20,slitname[i-1],$
			color=0,charsize=0.8
	    endif
	endfor

	; Write history information.
	xyouts,50,20,'Observation date='+obsdate+'T'+obstime,$
		charsize=0.90,color=0,/DEVICE
	xyouts,50,5,'IDF file date= '+filedate,charsize=0.90,color=0,/DEVICE
	xyouts,1124/2,5,'CALFUSE version '+cfvers,$
		charsize=0.90,color=0,/DEVICE,align=0.5
	xyouts,1124-50,20,'CHID calibration file='+chid_cal,$
		charsize=0.90,color=0,/DEVICE,align=1.0
	xyouts,1124-50,5,'Exposure time='+$
		strcompress(string(exptime,format="(F10.2)"),/rem),$
		charsize=0.90,color=0,/DEVICE,align=1.0

	;  Read in and plot the wavelength scale
	wl_lif=mrdfits(cal_path+'/'+wave_cal,4,ah,/SILENT)
	wl_sic=mrdfits(cal_path+'/'+wave_cal,8,ah,/SILENT)

	wlsc_lif=[wl_lif[0].wavelength,wl_lif[16383].wavelength]
	wlsc_sic=[wl_sic[0].wavelength,wl_sic[16383].wavelength]

	detector=strtrim(detector,2)

	; Determine locations of the wavelength axes.
	if (detector EQ '1A') then begin
		hip=maxy+25
		midp=425 
		lowp=miny-15
	endif else if (detector EQ '1B') then begin 
		hip=maxy+25
		midp=425
		lowp=miny-10
	endif else if (detector EQ '2A') then begin
		hip=maxy+100
		midp=500
		lowp=miny-100
	endif else if (detector EQ '2B') then begin
		hip=maxy+100
		midp=500
		lowp=miny-100
	endif else begin
		hip=maxy+25
		midp=512
		lowp=miny-15
	endelse

	; Draw wavelength axes.
	axis,0,hip,xaxis=1,xrange=wlsc_lif,xstyle=1,charsize=cs,$
		color=0,xtitle="Wavelength"
	axis,0,midp,xaxis=0,xrange=wlsc_lif,xstyle=1,charsize=0.001,$
		color=0,xtitle="Wavelength"

	axis,0,midp,xaxis=1,xrange=wlsc_sic,xs=1,charsize=0.001,$
		color=0,xtitle=""
	axis,0,lowp,xaxis=0,xrange=wlsc_sic,xs=1,charsize=cs,$
		color=0,xtitle="Wavelength"

	; Create a GIF/JPEG file: 
	ver = float(!version.release)
	if (ver ge 5.4 and ver le 6.0) then write_jpeg,rootname+'ext.jpg',TVRD() $
	else write_gif,rootname+'ext.gif',TVRD()

EXIT:
return
end
