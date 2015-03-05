/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Description: Modifies X array to correct for instrumental astigmatism.
 *
 *		Since we have no astigmatism correction for extended 
 *		sources, this program does not attempt to apply one.
 *
 *              In keeping with the CalFUSE standard, the astigmatism 
 *		correction is added to the measured X position of each photon.
 *		Versions of the ASTG_CAL file prior to v009 have the opposite 
 *		sign and are incompatible with this version of the program.
 *                    
 * Returns:     0 upon successful completion.
 *
 * History:     12/02/02   1.1   jch	Adapt from cf_astig.c
 *		12/20/02   1.2   wvd	Install
 *		02/11/03   1.3   wvd	Apply astig correction only if 
 *					YCENT# > O for target aperture.
 *					Add ASTG_KEY keyword.  Comment field:
 *					Performed? 1=LiF, 2=SiC, 3=Both
 *		02/13/03   1.4   wvd	Change all indices to type long.
 *              02/24/03   1.5   peb    Changed include file to calfusettag.h
 *              03/04/03   1.6   peb    Removed unused variables and math.h.
 *                                      Corrected sprintf:185 argument error.
 *		03/11/03   1.7   wvd	Change channel to type char
 *		03/12/03   1.8   wvd	If wavelength undefined, set
 *					channel[k] = 0
 *		04/04/03   1.9   wvd	Test for overflow of astig array.
 *		04/07/03   1.10  wvd	Delete test for quality of YCENT.
 *					Delete ASTG_KEY keyword.
 *					Don't round yoff to nearest integer.
 *		04/09/03   1.11  wvd	Write name of ASTG_CAL file to trailer.
 *              05/20/03   1.12  rdr    Added call to cf_proc_check
 *              09/16/03   1.13  wvd    Change sign of ASTG_CAL file so that
 *					correction is additive.
 *              09/22/03   1.14  wvd    Initialize overflow counter.
 *					Print out number of overflows.
 *              10/17/03   1.15  wvd    Discard photons that fall outside of 
 *					the astigmatism-correction window.
 *              10/31/03   1.16  wvd    Change channel to unsigned char.
 *              11/11/03   1.17  wvd    Interpolate between tabulated
 *					wavelength values.  Use cf_read_col
 *					to read wavecal files.
 *              02/17/04   1.18  wvd    Replace cf_nint() with cf_nlong()
 *              03/02/04   1.19  wvd    Change name of assign_wavelengths
 *					to cf_x2lambda and make it
 *					generally accessible.
 *              07/16/04   1.20  wvd    Must change ASTG_COR to COMPLETE or
 *					SKIPPED by hand.
 *              02/18/05   1.21  wvd    Add some diagnostic output.
 *              11/30/05   1.22  wvd    Don't complain when photon events
 *					fall outside of astig images.
 *              05/15/06   1.23  wvd    Divide cf_astigmatism_and_dispersion
 *					into two separate routines.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_astigmatism";
static char CF_VER_NUM[] = "1.23";


/* 
 * Loop over the apertures, read appropriate astigmatism correction,
 * and add to the input X shifts of pixels within extraction window.
 */
static int
add_astig_shifts(fitsfile *infits, long nevents, float *x, float *y,
	unsigned char *channel)
{
    char  	astig_file[FLEN_VALUE], keyword[FLEN_VALUE];
    int   	status=0;
    int   	ap, active_ap[2], ast_cen, extno, i, overflow;
    long  	k, ii, jj, kk, npix, nx, ny;
    float 	centroid, yoff;
    float 	*astig=NULL;
    fitsfile 	*astigfits;

    /* Read target apertures from the file header */
    (void) cf_source_aper (infits, active_ap);

    /* Open the astigmatism file */
    FITS_read_key(infits, TSTRING, "ASTG_CAL", astig_file, NULL, &status);
    FITS_open_file(&astigfits, cf_cal_file(astig_file), READONLY, &status);

    /*   Since we presently have no astig correction for extended
     *   sources, we apply correction only to target aperture.
     */
    for (i = 0; i < 2; i++) {
	ap = active_ap[i];

	/* Read the centroid of the aperture from the input file. */
	sprintf(keyword, "YCENT%1d", ap); 
    	FITS_read_key(infits, TFLOAT, keyword, &centroid, NULL, &status);
        cf_verbose(3, "Target aperture: %d\t Y centroid: %f", ap, centroid);

	/* Read astigmatism correction for this aperture */
	extno = ap+1;
	cf_verbose(3, "Reading extension %d of %s", extno, astig_file);
	FITS_movabs_hdu(astigfits, extno, NULL, &status);
	FITS_read_key(astigfits, TINT, "SLIT_CEN", &ast_cen, NULL, &status);
	FITS_read_key(astigfits, TLONG, "NAXIS1", &nx, NULL, &status);
	FITS_read_key(astigfits, TLONG, "NAXIS2", &ny, NULL, &status);

	npix = nx * ny;
	astig = (float *) cf_malloc(sizeof(float) * npix);
	FITS_read_img(astigfits, TFLOAT, 1L, npix, 0, astig, NULL, &status);

	/* Determine offset between spectrum and astig correction file */
	yoff = centroid - ast_cen;
	cf_verbose(3, "Offset between spectrum and astig correction: %f", yoff);

	/*  Go through the event list and find the astigmatism correction
	 *  appropriate for the photon's (x, y) position.
	 */
	overflow = 0;
	for (k = 0; k < nevents; k++) {
	    if (channel[k] == ap) {
		if (((ii = cf_nlong(x[k])) < 0 || ii >= nx) ||
		    ((jj = cf_nlong(y[k] - yoff)) < 0 || jj >= ny) ||
		    ((kk = jj * nx + ii) >= npix)) {
			channel[k] = 0;
			overflow++;
			cf_verbose(3, "Cannot correct pixel %05ld (%f, %f)",
				k, x[k], y[k]);
		}
		else x[k] += astig[kk]; 
	    } 
	}
	/* If overflow flag is set, there is a problem either with
		the ASTG_CAL file or the input photon list.
	if (overflow)
    	    cf_if_warning("%d events fell outside of aperture %d ASTG_CAL window",
		overflow, ap);
	Comment this out, as it is always triggered by cf_bad_pixels. - wvd (11/30/05) */

        /* Space for astig array is allocated in each loop. */
	free(astig);
    	cf_verbose(3, "End of loop for aperture %d", ap);
    }

    FITS_close_file(astigfits, &status);
    return status;
}


/*
 *  If APERTURE = RFPT or target is not a point source, exit without
 *  applying an astigmatism correction.
 */
int cf_astigmatism(fitsfile *infits, long nevents, float *x,
	float *y, unsigned char *channel)
{
    char 	aperture[FLEN_VALUE];
    char	source_type[FLEN_VALUE];
    int 	errflg=0, fileok=TRUE, status=0;

    /* Enter a timestamp into the log. */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

    FITS_read_key(infits, TSTRING, "APERTURE", aperture, NULL, &status);
    if (!strncmp(aperture, "RFPT", 4)) {
        cf_verbose(1, "Aperture is RFPT. No astig correction applied.");
	fileok = FALSE;
    }

    /* Check source type */
    FITS_read_key(infits, TSTRING, "SRC_TYPE", source_type, NULL, &status);
    if (strncmp(source_type, "P", 1)) {
	cf_verbose(1, "Source type is %s. No astig correction applied.",
	    source_type);
	fileok = FALSE;
    }

    if (fileok) {
        add_astig_shifts(infits, nevents, x, y, channel);
	cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    }
    else
	cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");

    /* Update processing flags. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");

    return fileok;
}
