/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE 
 *****************************************************************************
 *
 * Synopsis:    cf_extraction_limits(fitsfile *infits, int aperture,
 *		int srctype, short **ylow, short **yhigh, short *xmin,
 *		short *xmax)
 *
 * Description: Read the extraction limits for a given aperture and source
 *		type (point or extended) from the CHID_CAL file.  Shift these
 *		limits to match the centroid of the observed spectrum, if
 *		that measurement has been performed.  For point sources
 *		observed in HIST mode, pad YLOW and YHIGH by SPECBINY pixels.
 *
 * Arguments:   fitsfile  *infits       Pointer to the location of the FITS
 *                                      header of the Intermediate Data File
 *              int       aperture	aperture designation (1 through 8)
 *		int	  srctype	0 = point source, 1 = extended
 *              short     **ylow	lower Y bound for the extraction window
 *              short     **yhigh	upper Y bound for the extraction window
 *              short     *xmin 	lower X bound for the extraction window
 *              short     *xmax  	upper X bound for the extraction window
 *
 * Calls:       None
 *
 * Return:      npts = length of arrays ylow and yhigh
 *
 * History:     08/22/03    1.1  BJG    Based on subroutine from v1.6 of
 *					cf_bad_pixels.c by RDR.
 *					Change coltype from char to int in
 *					cf_read_col.
 *		03/15/05    1.2  wvd	Read HIST_PAD from appropriate HDU
 *					of CHID_CAL file and pad apertures
 *					if data are in HIST mode.
 *		03/16/05    1.3  wvd	For extended sources, return large
 *					apertures.
 *		10/06/05    1.4  wvd	Always pad apertures by 8 pixels
 *					in HIST mode.
 *		05/03/06    1.5  wvd	Always pad apertures by SPECBINY
 *					pixels in HIST mode.
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calfuse.h"

long cf_extraction_limits(fitsfile *infits, int aperture, int srctype,
	short **ylow, short **yhigh, short *xmin, short *xmax) { 

    char CF_PRGM_ID[] = "cf_extraction_limits"; 
    char CF_VER_NUM[] = "1.5";

    fitsfile *chidfits ;
    int status=0, hdutype;
    int dcent, hdu ; 
    short histpad=0, specbiny;
    short xmint, xmaxt, *ylowt=NULL, *yhight=NULL;
    float ycent_data, ycent_tab ;
    long j, npts ;
    char ycentname[FLEN_VALUE]={'\0'}, chidfile[FLEN_VALUE] ;
    char instmode[FLEN_VALUE];

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Read the measured centroid from the input file header. */
    FITS_movabs_hdu(infits, 1, &hdutype, &status) ;
    sprintf(ycentname,"YCENT%1d",aperture) ;
    FITS_read_key(infits, TFLOAT, ycentname, &ycent_data, NULL, &status);
    cf_verbose(3, "Reading extraction limits for aperture %d ", aperture) ;
    cf_verbose(3, "ycent keyword = %s, value=%5.1f ", ycentname, ycent_data) ;

    /* Read instrument mode and binning factor from file header. */
    FITS_read_key(infits, TSTRING, "INSTMODE", instmode, NULL, &status);
    FITS_read_key(infits, TSHORT, "SPECBINY", &specbiny, NULL, &status);

    /* Which HDU to read depends on both channel and source type. */
    hdu = aperture + 1;
    if (srctype) hdu += 8;

    /* Open the CHID_CAL file. */
    FITS_read_key(infits, TSTRING, "CHID_CAL", chidfile, NULL, &status);
    cf_verbose(3, "spectral extraction file = %s, HDU = %d", chidfile, hdu);
    FITS_open_file(&chidfits, cf_cal_file(chidfile), READONLY, &status);
    FITS_movabs_hdu(chidfits, hdu, &hdutype, &status);

    /* Read X limits of extraction window. */
    FITS_read_key(chidfits, TSHORT, "XMIN", &xmint, NULL, &status);
    FITS_read_key(chidfits, TSHORT, "XMAX", &xmaxt, NULL, &status);

    /* Read centroid of tabulated extraction window. */
    FITS_read_key(chidfits,TFLOAT,"CENTROID",&ycent_tab,NULL,&status);

    /* For point-source data taken in HIST mode, pad window. */
    if ((!strncmp(instmode, "HIST", 4)) && !srctype) histpad = specbiny;
	/* FITS_read_key(chidfits,TSHORT,"HIST_PAD",&histpad,NULL,&status); */

    /* Read the upper and lower boundaries of the extraction slit. */
    npts=cf_read_col(chidfits, TSHORT, "YLOW",  (void **) &ylowt);
    npts=cf_read_col(chidfits, TSHORT, "YHIGH", (void **) &yhight);
    FITS_close_file(chidfits, &status) ;

    /* If the spectral centroid has been measured, shift the aperture
       limits so that the tabulated centroid equals the measured value. */
    if (ycent_data > 1 && ycent_data < NYMAX) {
	dcent = (short) cf_nint(ycent_data - ycent_tab) ;
        cf_verbose(3, "shifting extraction window: ycent_tab=%5.1f, delta=%d",
            ycent_tab, dcent) ;
        for (j = xmint; j <= xmaxt; j++) {
	    ylowt[j]  += dcent ;
            yhight[j] += dcent ; 
	}
    }

    /* If HIST_PAD > 0, pad the extraction window, top and bottom. */
    if (histpad > 0) {
        cf_verbose(3, "HIST point source.  "
	    "Padding extraction window by +/- %d pixels", histpad);
        for (j = xmint; j <= xmaxt; j++) {
	    ylowt[j]  -= histpad ;
	    yhight[j] += histpad ; 
	}
    }

    /* Return the extraction limits. */
    *xmin = xmint;
    *xmax = xmaxt;
    *ylow = ylowt ;
    *yhigh = yhight ;

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return npts ;
}
