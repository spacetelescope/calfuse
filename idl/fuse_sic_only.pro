pro fuse_sic_only, input
;+
; NAME:
;	FUSE_EXTRACT
;
;*PURPOSE:
;       To extract full-resolution spectra from the primary LiF and SiC
;	apertures of an intermediate data file.
;
;	THIS VERSION EXTRACTS ONLY A SIC SPECTRUM.
;
;*CATEGORY:
;	INPUT/OUTPUT
;
;*CALLING SEQUENCE:
;	FUSE_EXTRACT, input
;
;*INPUT:
;	input: intermediate data file
;
;*OUTPUT:
;	LiF and SiC spectral files in v2.x format.
;
;		
;*PROCEDURES USED:
;	FXADDPAR, FXHMAKE, FXBADDCOL, FXBCREATE, FXBWRITE, FXBFINISH, MRDFITS,
;	SXPAR
;
;*HISTORY:
;	Written by:	V. Dixon	March 2003
;	03/23/2003	V. Dixon	Correct X values for Doppler shift.
;	03/24/2003	V. Dixon	Fix bug in call to where().
;       06/27/2003      R. Robinson     Update to new version of IDF
;                                        - added fscale keyword to mrdfits call
;                                        - divide flux by exposure time
;	06/08/2004	V. Dixon	Use EXPTIME rather than RAWTIME.
;					Correct calculation of COUNTS
;	07/01/2004	V. Dixon	Don't read headers of extensions.
;					Employ faster scheme for calculating
;					Doppler shift.
;	07/02/2004	V. Dixon	Don't compute dl per pixel.  It suffers
;					from truncation errors.  Instead, use
;					value from file header.
;	07/15/2004	V. Dixon	Extract only SiC spectrum.
;-
;-----------------------------------------------------------------------------
;
; print calling sequence if no parameters supplied
;
	if n_params(0) lt 1 then begin
	    print,'CALLING SEQUENCE: fuse_sic_only, input'
	    return
	end
;
; Physical constants
;
	C = 2.99792458e5	; km/s
;
; Open input file and read keywords from header.
;
	dummy = mrdfits(input, 0, hdr, status=status)
	rootname = SXPAR(hdr, 'ROOTNAME')
	vhelio   = SXPAR(hdr, 'V_HELIO')
	detector = SXPAR(hdr, 'DETECTOR')
	aperture = SXPAR(hdr, 'APERTURE')
	wavecal  = SXPAR(hdr, 'WAVE_CAL')
        exptime  = SXPAR(hdr, 'EXPTIME')
	
	ap_array = ['HIRS', 'MDRS', 'LWRS', 'PINH']
	ap_code = ['3','2','4','1']
	for k = 0, 3 do if (strcmp(aperture, ap_array[k], 1)) then ap = k
	channel = ap+1

	GET_DATE, FITS_date
;
; Read photon information in the first extension of the IDF.
;
	a = mrdfits(input, 1, status=status, /fscale)
;
; Read timeline table from the third extension of the IDF.
;
	b = mrdfits(input, 3, status=status, /fscale)
	nseconds = n_elements(b.time)
	if (nint(b.time[nseconds-1] - b.time[0]) gt nseconds) then $
		print, 'WARNING: '+input+' timeline table not continuous.'

	; for k = channel, channel+4, 4 do begin
	k = channel + 4
            print, 'Analyzing extension ', k
	;
	; Read corresponding wavelength array from WAVE_CAL file.
	;
            wcalfile=strtrim('/data1/fuse/calfuse/v3.0/calfiles/'+wavecal,2)
            print, 'Reading wavelength data from ', wcalfile
	    w=mrdfits(wcalfile,k,w_hdr)
            wave = w.wavelength
	    wpc = SXPAR(w_hdr, 'DISPAPIX')
	;
	; Select photons from aperture of interest.
	;
	    g = where (a.channel eq k and (a.timeflgs and not 1B) eq 0 and $
			(a.loc_flgs and 16B) eq 0, n)
            print, 'Number of photons in spectrum = ', n 
	    x = a.x[g]
	;
	; Correct X values for Doppler shift.
	;
            print, 'Applying Doppler shift' 
	    t = nint(a.time[g] - b.time[0])
	    ovc = (b.ORBITAL_VEL + vhelio)/C
	    for i = 0L, n-1L do x[i] = x[i] + ovc[t[i]] * a.lambda[g[i]]/wpc
	;
	; Extract spectrum and populate output arrays.
	;
            print, 'Extracting spectrum' 
	    x = nint(x)
	    flux = fltarr(16384)
	    error = fltarr(16384)
	    cntserr = fltarr(16384)
	    for i = 0L, n-1L do flux[x[i]] = flux[x[i]] + a.ergcm2[g[i]]
	    flux = flux / exptime / abs(wpc)
	    counts = histogram(x, min=0, max=16383)
	    g = where (counts gt 0)
	    cntserr[g] = sqrt(counts[g])
	    error[g] = flux[g]/cntserr[g]
	    quality = indgen(16384)
	;
	; Create primary HDU of output file
	;
	    if (k lt 5) then ch = 'lif' else ch = 'sic'
	    filename = strcompress(rootname+strlowcase(detector)+ch+ap_code[ap]+'ttagfcal.fit',/remove_all)
	    FXADDPAR, hdr, 'DATE', FITS_date
	    FXADDPAR, hdr, 'FILENAME', filename
	    FXADDPAR, hdr, 'FILETYPE', 'CALIBRATED EXTRACTED SPECTRUM'
	    FXADDPAR, hdr, 'APER_ACT', strupcase(strcompress(ap_array[ap]+'_'+ch,/remove_all))
	    FXHMAKE, hdr, /EXTEND
	    FXWRITE, filename, hdr
	    
	    FXBHMAKE, HEADER, 1, 'FUSE 1D Spectrum'
	    FXBADDCOL, iwave,    HEADER, wave,    'WAVE',    TUNIT='ANGSTROMS'
	    FXBADDCOL, iflux,    HEADER, flux,    'FLUX',    TUNIT='ERG/CM2/S/A'
	    FXBADDCOL, ierror,   HEADER, error,   'ERROR',   TUNIT='ERG/CM2/S/A'
	    FXBADDCOL, iquality, HEADER, quality, 'QUALITY', TUNIT='UNITLESS'
	    FXBADDCOL, icounts,  HEADER, counts,  'COUNTS',  TUNIT='COUNTS'
	    FXBADDCOL, icntserr, HEADER, cntserr, 'CNTSERR', TUNIT='COUNTS'

	    FXBCREATE, out, filename, HEADER
	    FXBWRITE, out, wave,    iwave,    1
	    FXBWRITE, out, flux,    iflux,    1
	    FXBWRITE, out, error,   ierror,   1
	    FXBWRITE, out, quality, iquality, 1
	    FXBWRITE, out, counts,  icounts,  1
	    FXBWRITE, out, cntserr, icntserr, 1
	    FXBFINISH, out
	;endfor

	close, /all
	return
	end
