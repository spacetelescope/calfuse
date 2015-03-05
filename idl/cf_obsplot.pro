;+
;  cf_obsplot.pro reads the final, combined *00000all* spectral file
;       for a given observation and produces quick-look plots for
;	each of the four FUSE detectors.  A fifth plot shows the 
;	best three channels that span the regions 900-1000,
;	1000-1100, and 1100-1200 A.
;
;	Author:  Ken Sembach
;	Last Update: February 23, 2001
;
;	V. Dixon, 07/02/2001 - Modified for inclusion in pipeline.
;       V. Dixon, 08/16/2001 - Read and write to current directory.
;       V. Dixon, 08/31/2001 - Add a sixth plot spanning FUSE band.
;				Mark region of worm in LiF 1B.
;	V. Dixon, 11/05/2001 - Compute mean flux with 1-sig errors.
;	V. Dixon, 11/20/2001 - Modify position of mean-flux plots.
;				Add sentence describing error bars.
;	V. Dixon, 01/17/2002 - Deal gracefully with missing data.
;	V. Dixon, 03/12/2002 - Create .gif files if possible.
;       B. Godard, 04/26/2004 - Print units and CalFUSE Ver
;                               Don't keep exposures with exptime=0
;                               or quality array=0
;                               Use quality array to compute min/max
;	V. Dixon, 05/04/2004 - Scale ymax by 1.1 in all plots.
;				Change output file names to match
;				MAST's expectations.
;	V. Dixon, 08/20/2004 - Change A to Angstroms at the bottom
;				of each page.
;				Exclude C III 977 when setting Y
;				scale.
;	V. Dixon, 04/14/2005 - Rewrite the program to use the new
;				*all* files produced by CalFUSE v3.1.
;	V. Dixon, 08/01/2005 - If IDL version > 6.0, generate GIF file.
;	V. Dixon, 03/21/2006 - If channel has no data, set mflux = -1.
;			       Use "Exposure" or "Exposures" as 
;				appropriate.
;			       Consider all four channels for the
;				1000-1100 A plot.
;			       Pass xrange and yscale to TOP_PLOT 
;				and BOTTOM_PLOT.
;	V. Dixon, 05/19/2006 - For BKGD targets, omit flux
;				comparison when selecting which
;				detector segments to use.
;			       Compute mean flux using same wavelength
;				intervals for each band.
;	V. Dixon, 05/24/2006 - For PC targets, don't include O VI in
;				calculation of mean flux.
;	V. Dixon, 12/12/2006 - If OBSTIME < 10 s for preferred segment
;				and > 100 s for other segment, use other.
;	V. Dixon, 12/19/2006 - Use the same sample regions in all plots.
;	V. Dixon, 04/06/2007 - Comment out calls to set_plot,'X' and 
;				CLEANPLOT, as they can cause problems
;				for OPUS.
;	V. Dixon, 04/11/2008 - If obstime < 1000, quote in seconds.
;	V. Dixon, 08/08/2008 - If airglow keyword is set, image files
;				are named *00900*
;
;       Calling sequence: OBSPLOT, obsname
;
;       Input:    obsname (e.g., P1030604)
;	Output:	  obsname00000lif1ttagf.jpg, obsname00000lif2ttagf.jpg, 
;	          obsname00000sic1ttagf.jpg, obsname00000sic2ttagf.jpg, 
;	          obsname00000specttagf.jpg
;		
;
;-

pro TOP_PLOT, data, header, xrange, yscale, aper, sample, mflux, title
   s=sample
   !x.range=xrange
   obstime = SXPAR(header,'OBSTIME')
   extname = SXPAR(header,'EXTNAME')
   num_exp = SXPAR(header,'NSPEC')
   if (num_exp eq 1) then exposures=' Exposure, ' else exposures=' Exposures, '
   IF (obstime gt 0) THEN BEGIN
      w = data.wave & f = data.flux 
      loc = WHERE(((w GT s[0]) AND (w LE s[1])) OR ((w GT s[2]) and (w LT s[3])) OR $
	((w GT s[4]) and (w LT s[5])))
      ftmp = SMOOTH(f[loc],6) & nftmp = N_ELEMENTS(ftmp)
      ymin = MIN(f[loc]) & ymax = MAX(ftmp)
      ydiv = FIX(1-ALOG10(ymax))
      ymax = yscale * ymax *10.^ydiv
      !y.range=[ymin,ymax]
      mflux = mean(f[loc])
      
      if  (obstime lt 1000) then time_str = STRING(round(obstime),'(I3)')+' sec)' else $
      time_str = STRING(obstime/1000,'(F4.1)')+' ksec)'

      PLOT,w,f*10.^ydiv,nsum=12,title=title,ytitle='Flux (10!u-'+STRTRIM(ydiv,2)+'!n)',xtitle=' '
      XYOUTS,(!x.range[1]-!x.range[0])*0.5+!x.range[0],(!y.range[1]-!y.range[0])*0.88+!y.range[0], $
       extname+' '+aper+' ('+STRTRIM(num_exp,2)+exposures+time_str,charsize=1.3,align=0.5
   ENDIF ELSE BEGIN
      !y.range=[0,1]
      mflux = -1.
      PLOT,[0,0],[0,0],title=title,ytitle='Flux',xtitle=' '
      XYOUTS,(!x.range[1]-!x.range[0])*0.5+!x.range[0],0.88, $
       extname+' '+aper+': No data found',charsize=1.5,align=0.5
   ENDELSE
   RETURN
END

pro BOTTOM_PLOT, data, header, xrange, yscale, aper, sample, mflux, vers
   s=sample
   !x.range=xrange
   obstime = SXPAR(header,'OBSTIME')
   extname = SXPAR(header,'EXTNAME')
   num_exp = SXPAR(header,'NSPEC')
   if (num_exp eq 1) then exposures=' Exposure, ' else exposures=' Exposures, '
   IF (obstime gt 0) THEN BEGIN
      w = data.wave & f = data.flux
      loc = WHERE(((w GT s[0]) AND (w LE s[1])) OR ((w GT s[2]) AND (w LT s[3])) OR $
	((w GT s[4]) AND (w LT s[5])))
      ftmp = SMOOTH(f(loc),6) & nftmp = N_ELEMENTS(ftmp)
      ymin = MIN(f[loc]) & ymax = MAX(ftmp)
      ydiv = FIX(1-ALOG10(ymax))
      ymax = yscale * ymax *10.^ydiv
      !y.range=[ymin,ymax]
      mflux = mean(f[loc])

      if  (obstime lt 1000) then time_str = STRING(round(obstime),'(I3)')+' sec)' else $
      time_str = STRING(obstime/1000,'(F4.1)')+' ksec)'

      PLOT,w,f*10.^ydiv,nsum=12,title=' ',ytitle='Flux (10!u-'+STRTRIM(ydiv,2)+'!n)', $
       xtitle='Wavelength ('+string("305B)+')',/noeras
      XYOUTS,(!x.range[1]-!x.range[0])*0.5+!x.range[0],(!y.range[1]-!y.range[0])*0.88+!y.range[0], $
        extname+' '+aper+' ('+STRTRIM(num_exp,2)+exposures+time_str, charsize=1.3,align=0.5
      IF (STRCMP(extname,'1BLIF',5, /FOLD_CASE) EQ 1) THEN $
          XYOUTS, 1160., (!y.range[1]-!y.range[0])*0.1+!y.range[0],'---- REGION OF WORM ----',align=0.5
   ENDIF ELSE BEGIN
      !y.range=[0,1]
      mflux = -1.
      PLOT,[0,0],[0,0],title=' ',ytitle='Flux', xtitle='Wavelength ('+string("305B)+')',/noeras
      XYOUTS,(!x.range[1]-!x.range[0])*0.5+!x.range[0],0.88, $
       extname+' '+aper+': No data found',charsize=1.5,align=0.5
   ENDELSE
   XYOUTS, 0.05, 0.01, 'Flux units are erg cm!E-2!N s!E-1!N '+string("305B)+'!E-1!N.', alignment=0., /NORMAL
   XYOUTS, 0.95, 0.01, 'CalFUSE v'+vers, alignment=1., /NORMAL
   RETURN
END

PRO CF_OBSPLOT, obsname, airglow=airglow

   ver = float(!version.release)

   if keyword_set(airglow) then fnumber = '00900' else fnumber = '00000'
   list = FINDFILE(obsname+fnumber+'all*fcal.fit', COUNT=count)
   lun = fxposit(list[0], 0, /SILENT)
   foo = MRDFITS(lun,0,h,/SILENT)

   root = SXPAR(h,'PRGRM_ID')+SXPAR(h,'TARG_ID')+SXPAR(h,'SCOBS_ID')
   root = STRCOMPRESS(root, /REMOVE_ALL)
   targ = STRCOMPRESS(SXPAR(h,'TARGNAME'), /REMOVE_ALL)
   objclass = SXPAR(h,'OBJCLASS')
   sp_type = STRCOMPRESS(SXPAR(h,'SP_TYPE'), /REMOVE_ALL)
   src_type = STRCOMPRESS(SXPAR(h,'SRC_TYPE'), /REMOVE_ALL)
   mode = STRLOWCASE(STRCOMPRESS(SXPAR(h,'INSTMODE'), /REMOVE_ALL))
   aper = STRCOMPRESS(SXPAR(h,'APERTURE'), /REMOVE_ALL)
   vers = STRCOMPRESS(SXPAR(h,'CF_VERS'), /REMOVE_ALL)
   ext = '4'
   IF aper EQ 'MDRS' THEN ext = '2'
   IF aper EQ 'HIRS' THEN ext = '3'
   if keyword_set(airglow) then title = ''+root+' [AIRGLOW]' else $
	title = ''+root+' ['+targ+']'
   if (src_type eq 'PE') then yscale = 1.0 else yscale = 1.3
   if (objclass eq 7 or sp_type eq 'BKGD') then bkgd_obs = 1 else bkgd_obs = 0

;
; Read each channel and its header.
;
   L1A = mrdfits(lun,0,HL1A)  
   L1B = mrdfits(lun,0,HL1B)  
   L2B = mrdfits(lun,0,HL2B)  
   L2A = mrdfits(lun,0,HL2A)  
   S1A = mrdfits(lun,0,HS1A)  
   S1B = mrdfits(lun,0,HS1B)  
   S2B = mrdfits(lun,0,HS2B)  
   S2A = mrdfits(lun,0,HS2A)  

;
; Use different sample regions for point and extended sources.
;
   S1B_SAMPLE = [912, 935, 955, 970, 980, 985]
   S2A_SAMPLE = [912, 935, 955, 970, 980, 985]

   if (src_type eq 'PC') then begin
	L1A_SAMPLE = [0, 0, 0, 0, 1045, 1070]
	S1A_SAMPLE = [0, 0, 0, 0, 1045, 1070]
	L2B_SAMPLE = [0, 0, 0, 0, 1045, 1070]
	S2B_SAMPLE = [0, 0, 0, 0, 1045, 1070]
   endif else begin
	L1A_SAMPLE = [0, 0, 1030, 1039, 1045, 1070]
	S1A_SAMPLE = [0, 0, 1030, 1039, 1045, 1070]
	L2B_SAMPLE = [0, 0, 1030, 1039, 1045, 1070]
	S2B_SAMPLE = [0, 0, 1030, 1039, 1045, 1070]
   endelse

   L1B_SAMPLE = [1095, 1130, 1140, 1165, 1170, 1190]
   L2A_SAMPLE = [1095, 1130, 1140, 1165, 1170, 1190]

;
; Plot individual channels.  Skip this for airglow files.
;
   ;set_plot,'X'
   ;CLEANPLOT
   set_plot,'Z'
   device,set_resolution=[650,580]

   !x.style = 1 & !y.style = 1
   !x.ticklen=0.045
   !x.charsize=1.5
   !y.charsize=1.5
   !y.minor = 2
   !p.title = ' '
   !x.title = ' '
   !y.title = ' '
   !p.color = 0
   !p.background = 255

;
; LiF 1A and 1B
;
   !p.position = [0.13,0.57,0.95,0.95]
   TOP_PLOT, L1A, HL1A, [985,1085], yscale, aper, L1A_SAMPLE, FL1A, title 

   !p.position = [0.13,0.1,0.95,0.48]
   BOTTOM_PLOT, L1B, HL1B, [1092,1190], yscale, aper, L1B_SAMPLE, FL1B, vers

   if (ver ge 5.4 and ver le 6.0) then write_jpeg,root+fnumber+'lif1'+mode+'f.jpg',TVRD() $
   else write_gif,root+fnumber+'lif1'+mode+'f.gif',TVRD()

;
; SiC 1A and 1B
;
   !p.position = [0.13,0.57,0.95,0.95]
   TOP_PLOT, S1A, HS1A, [1000,1094], yscale, aper, S1A_SAMPLE, FS1A, title

   !p.position = [0.13,0.1,0.95,0.48]
   BOTTOM_PLOT, S1B, HS1B, [910,995], yscale, aper, S1B_SAMPLE, FS1B, vers

   if (ver ge 5.4 and ver le 6.0) then write_jpeg,root+fnumber+'sic1'+mode+'f.jpg',TVRD() $
   else write_gif,root+fnumber+'sic1'+mode+'f.gif',TVRD()

;
; LiF 2A and 2B
;
   !p.position = [0.13,0.57,0.95,0.95]
   TOP_PLOT, L2A, HL2A, [1085,1185], yscale, aper, L2A_SAMPLE, FL2A, title

   !p.position = [0.13,0.1,0.95,0.48]
   BOTTOM_PLOT, L2B, HL2B, [975,1078], yscale, aper, L2B_SAMPLE, FL2B, vers

   if (ver ge 5.4 and ver le 6.0) then write_jpeg,root+fnumber+'lif2'+mode+'f.jpg',TVRD() $
   else write_gif,root+fnumber+'lif2'+mode+'f.gif',TVRD()

;
; SiC 2A and 2B
;
   !p.position = [0.13,0.57,0.95,0.95]
   TOP_PLOT, S2A, HS2A, [913,1009], yscale, aper, S2A_SAMPLE, FS2A, title

   !p.position = [0.13,0.1,0.95,0.48]
   BOTTOM_PLOT, S2B, HS2B, [1013,1106], yscale, aper, S2B_SAMPLE, FS2B, vers

   if (ver ge 5.4 and ver le 6.0) then write_jpeg,root+fnumber+'sic2'+mode+'f.jpg',TVRD() $
   else write_gif,root+fnumber+'sic2'+mode+'f.gif',TVRD()

;
; Now plot combined fluxes
;
   ;set_plot,'X'
   ;CLEANPLOT
   set_plot,'Z'
   device,set_resolution=[650,870]

   !x.style = 1 & !y.style = 1
   !x.ticklen=0.045
   !x.charsize=1.5
   !y.charsize=1.5
   !y.minor = 2
   !p.title = ' '
   !x.title = ' '
   !y.title = ' '
   !p.color = 0
   !p.background = 255

;
;        900 - 1000 A
;
   !p.position = [0.13,0.70,0.95,0.95]
   obs2a = SXPAR(HS2A, 'OBSTIME')
   obs1b = SXPAR(HS1B, 'OBSTIME')

   if ((bkgd_obs eq 1 and FS2A gt -1) or (FS2A gt 0.9*FS1B) or (obs2a gt 100 and obs1b lt 10)) then begin
      TOP_PLOT, S2A, HS2A, [913,1009], yscale, aper, S2A_SAMPLE, foo, title
   endif else begin
      TOP_PLOT, S1B, HS1B, [910,995], yscale, aper, S1B_SAMPLE, foo, title
   endelse

;
;       1000 - 1100 A
;
   !p.noerase = 1
   !p.position = [0.13,0.40,0.95,0.65]
   obl1a = SXPAR(HL1A, 'OBSTIME')
   obl2b = SXPAR(HL2B, 'OBSTIME')

   if ((bkgd_obs eq 1 and FL1A gt -1) or $
	((FL1A gt 0.9*FL2B and FL1A gt 0.7*FS1A and FL1A gt 0.7*FS2B) $
	  and not (obl2b gt 100 and obl1a lt 10))) then $
        TOP_PLOT, L1A, HL1A, [985,1085], yscale, aper, L1A_SAMPLE, foo, '' else $
   if ((bkgd_obs eq 1 and FL2B gt -1) or $ 
	((FL2B gt 0.7*FS1A and FL2B gt 0.7*FS2B) and not (obl2b lt 10 and obl1a gt 100))) then $
        TOP_PLOT, L2B, HL2B, [975,1078], yscale, aper, L2B_SAMPLE, foo, '' else $
   if ((bkgd_obs eq 1 and FS1A gt -1) or ((FS1A gt FS2B) and not (obl1a lt 10 and obl2b gt 100))) then $
        TOP_PLOT, S1A, HS1A, [1000,1094], yscale, aper, S1A_SAMPLE, foo, '' else $
        TOP_PLOT, S2B, HS2B, [1013,1106], yscale, aper, S2B_SAMPLE, foo, ''

;
;       1100 - 1200 A
;
   !p.position = [0.13,0.1,0.95,0.35]
   obl2a = SXPAR(HL2A, 'OBSTIME')
   obl1b = SXPAR(HL1B, 'OBSTIME')

   if ((bkgd_obs eq 1 and FL2A gt -1) or (FL2A gt 0.9*FL1B) or (obl2a gt 100 and obl1b lt 10)) then $
      BOTTOM_PLOT, L2A, HL2A, [1085,1182], yscale, aper, L2A_SAMPLE, foo, vers else $
      BOTTOM_PLOT, L1B, HL1B, [1092,1190], yscale, aper, L1B_SAMPLE, foo, vers

   if (ver ge 5.4 and ver le 6.0) then write_jpeg,root+fnumber+'spec'+mode+'f.jpg',TVRD() $
   else write_gif,root+fnumber+'spec'+mode+'f.gif',TVRD()

   RETURN
END

