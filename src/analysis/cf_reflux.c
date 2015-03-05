/*****************************************************************************
 *		Johns Hopkins University
 *		Center For Astrophysical Sciences
 *		FUSE
 *
 * Synopsis:	cf_reflux input_file output_file aeff_cal1 [aeff_cal2]
 *
 * Description: Program applies a new flux calibration to an extracted
 *		spectral file.  If two effective-area files are given,
 *		program interpolates between them based on the MDJ
 *		of the exposure.
 *		Program writes a history line to output file header.
 *
 * Arguments:	input_file	FUSE extracted spectral file (fcal.fit)
 *		output_file	Copy of input file with modified flux cal.
 *		aeff_cal1	Effective-area curve
 *		aeff_cal2	Optional effective-area curve, needed only
 *				if you wish to interpolate.
 *
 * Returns:     none
 *
 * History:     06/02/05 1.1    wvd	Based on programs cf_uninterp and
 *					cf_extract_spectra.
 *
 ****************************************************************************/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_reflux";
static char CF_VER_NUM[] = "1.1";

int main(int argc, char *argv[])
{
    char  aper_act[FLEN_VALUE], datestr[FLEN_VALUE];
    int   status=0, hdutype;
    int   aperture, timeref;
    long  i, frow=1, felement=1, nout;
    float exptime, wpc;
    float *wave, *flux, *error, *weights, *bkgd;
    unsigned char *channel;
    fitsfile *infits, *outfits;

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    if (argc < 4) {
	printf("Usage: cf_reflux input_file output_file aeff_cal1 [aeff_cal2]\n");
	return 0;
    }

    /* Copy input spectral file to output.  Close input file. */
    FITS_open_file(&infits, argv[1], READONLY, &status);
    FITS_create_file(&outfits, argv[2], &status);
    FITS_copy_hdu(infits, outfits, 0, &status);
    FITS_movrel_hdu(infits, 1, &hdutype, &status);
    FITS_copy_hdu(infits, outfits, 0, &status);
    FITS_close_file(infits, &status);

    /* Modify output file keywords. */
    FITS_movabs_hdu(outfits, 1, &hdutype, &status);
    FITS_update_key(outfits, TSTRING, "FILENAME", argv[2], NULL, &status);
    FITS_update_key(outfits, TSTRING, "AEFF1CAL", argv[3], NULL, &status);
    if (argc == 5)
	FITS_update_key(outfits, TSTRING, "AEFF2CAL", argv[4], NULL, &status);
    else
	FITS_update_key(outfits, TSTRING, "AEFF2CAL", argv[3], NULL, &status);

    /* Add HISTORY lines to the output file */
    fits_get_system_time(datestr, &timeref, &status);
    strcat(datestr, "   File recalibrated using cf_reflux");
    FITS_write_history(outfits, datestr, &status);

    /* Through which aperture was this spectrum obtained? */
    FITS_read_key(outfits, TSTRING, "APER_ACT", aper_act, NULL, &status);
         if (!strcmp(aper_act, "HIRS_LIF")) aperture = 1;
    else if (!strcmp(aper_act, "MDRS_LIF")) aperture = 2;
    else if (!strcmp(aper_act, "LWRS_LIF")) aperture = 3;
    else if (!strcmp(aper_act, "HIRS_SIC")) aperture = 5;
    else if (!strcmp(aper_act, "MDRS_SIC")) aperture = 6;
    else if (!strcmp(aper_act, "LWRS_SIC")) aperture = 7;
    else {
	printf("Unable to interpret header keyword APER_ACT = %s\n", aper_act);
        FITS_close_file(outfits, &status);
	return 0;
    }

    /* Read WAVE, FLUX, ERROR, WEIGHTS, and BKGD arrays */
    FITS_movabs_hdu(outfits, 2, &hdutype, &status);
    nout = cf_read_col(outfits, TFLOAT, "WAVE", (void **) &wave);
    nout = cf_read_col(outfits, TFLOAT, "FLUX", (void **) &flux);
    nout = cf_read_col(outfits, TFLOAT, "ERROR", (void **) &error);
    nout = cf_read_col(outfits, TFLOAT, "WEIGHTS", (void **) &weights);
    nout = cf_read_col(outfits, TFLOAT, "BKGD", (void **) &bkgd);
    channel = (unsigned char *) cf_malloc(sizeof(unsigned char) * nout);

    /* Compute inputs to flux-calibration routine */
    for (i = 0; i < nout; i++) {
	channel[i] = (unsigned char) aperture;
	if (flux[i] != 0.) error[i] /= flux[i];		/* Relative error */
	flux[i] = weights[i] - bkgd[i];			/* Units are counts */
    }
    
    /* Apply flux calibration */
    FITS_movabs_hdu(outfits, 1, &hdutype, &status);
    FITS_update_key(outfits, TSTRING, "FLUX_COR", "PERFORM", NULL, &status) ;
    cf_convert_to_ergs(outfits, nout, flux, flux, channel, wave);

    /* Convert flux and error arrays to units of erg/cm2/s/A. */
    FITS_movabs_hdu(outfits, 1, &hdutype, &status);
    FITS_read_key(outfits, TFLOAT, "EXPTIME", &exptime, NULL, &status);
    FITS_read_key(outfits, TFLOAT, "WPC", &wpc, NULL, &status);
    if (exptime < 1) {
	printf ("EXPTIME = %g.  Exiting.\n", exptime);
        FITS_close_file(outfits, &status);
	return 0;
    }
    for (i = 0; i < nout; i++) {
        flux[i]  /= exptime * wpc;
        error[i] *= flux[i];
    }
 
    /* Write new flux and error arrays to the output file. */
    FITS_movabs_hdu(outfits, 2, &hdutype, &status);
    FITS_write_col(outfits, TFLOAT, 2, frow, felement, nout, flux, &status);
    FITS_write_col(outfits, TFLOAT, 3, frow, felement, nout, error, &status);
 
    /* Close output file. */
    FITS_close_file(outfits, &status);

    /* Free memory. */
    free(wave);
    free(flux);
    free(error);
    free(weights);
    free(bkgd);
    free(channel);

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
