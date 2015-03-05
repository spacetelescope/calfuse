/*****************************************************************************
 *            Johns Hopkins University
 *            Center For Astrophysical Sciences
 *            FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_make_wave_array(fitsfile *header, int aperture, long *nout,
 *		float **wave_out)
 *
 * Description: Reads W0, WMAX, and WPC from the scrn*.fit file and returns
 *		output wavelength array and its size.  If not supplied by
 *		user, values are taken from the WAVE_CAL file.
 *
 *              The FITS wavelength-calibration file has 8 extensions in the
 *              following order:
 *              Ext# HDU# Channel Aperture
 *                1   2     LiF     HIRS
 *                2   3     LiF     MDRS
 *                3   4     LiF     LWRS
 *                4   5     LiF     PINH
 *                5   6     SiC     HIRS
 *                6   7     SiC     MDRS
 *                7   8     SiC     LWRS
 *                8   9     SiC     PINH
 *
 * Arguments:   fitsfile  *header       Pointer to IDF FITS file header
 *              int       aperture      Target aperture (values 1-8)
 *              long      *nout         Length of output arrays
 *              float     **wave_out    Pointer to output wavelength array
 *
 * Returns:     0 on success
 *
 * History:     02/24/03  1.1   wvd	Initial coding.
 *		03/02/03  1.2   wvd	Write W0 and WPC to IDF header.
 *		09/29/03  1.3   wvd	Change calfusettag.h to calfuse.h
 *		10/15/03  1.4   wvd	Read either LIF or SIC wavelength
 *					parameters from PARM_CAL file,
 *					depending upon aperture.
 *		03/16/04  1.5   wvd	Cast wave_out array as type float.
 *
 ****************************************************************************/

#include <stdio.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_make_wave_array";
static char CF_VER_NUM[] = "1.5";

int
cf_make_wave_array(fitsfile *header, int aperture, long *nout,
	float **wave_out)
{
    char  parm_file[FLEN_VALUE];	/* Name of parameter file */
    char  wave_file[FLEN_VALUE]; 	/* Name of WAVE_CAL file */
    char  w0_key[FLEN_KEYWORD];		/* Depends on aperture */
    char  wmax_key[FLEN_KEYWORD];	/* Depends on aperture */
    char  wpc_key[FLEN_KEYWORD];	/* Depends on aperture */

    int   status=0, hdunum;
    long  i;

    float default_w0, default_wmax, default_wpc; /* Default values */
    float w0, wmax, wpc;			 /* Requested values */

    fitsfile *parmfits, *wavefits;

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin processing");

    /* 
     *  Set keywords to read either LIF or SIC parameters, depending on
     *	the aperture number.
     */
    if (aperture < 5) {
	sprintf(w0_key,   "LIF_W0");
	sprintf(wmax_key, "LIF_WMAX");
	sprintf(wpc_key,  "LIF_WPC");
    } else {
	sprintf(w0_key,   "SIC_W0");
	sprintf(wmax_key, "SIC_WMAX");
	sprintf(wpc_key,  "SIC_WPC");
    }

    /* Read recommended values of W0, WMAX, and WPC from WAVE_CAL header */
    hdunum = aperture + 1;
    FITS_read_key(header, TSTRING, "WAVE_CAL", wave_file, NULL, &status);
    FITS_open_file(&wavefits, cf_cal_file(wave_file), READONLY, &status);
    FITS_movabs_hdu(wavefits, hdunum, NULL, &status);
    FITS_read_key(wavefits, TFLOAT, "W0", &default_w0, NULL, &status);
    FITS_read_key(wavefits, TFLOAT, "WMAX", &default_wmax, NULL, &status);
    FITS_read_key(wavefits, TFLOAT, "WPC", &default_wpc, NULL, &status);
    FITS_close_file(wavefits, &status);
    cf_verbose(3, "Recommended wavelength parameters:");
    cf_verbose(3, "\tW0 = %g, WMAX = %g, WPC = %g",
	default_w0, default_wmax, default_wpc);

    /* 
     * Read requested values of W0, WMAX, and WPC from PARM_CAL.
     * If no value requested, use default values.
     * Compare requested values with default; complain if weird.
     */
    FITS_read_key(header, TSTRING, "PARM_CAL", parm_file, NULL, &status);
    FITS_open_file(&parmfits, cf_parm_file(parm_file), READONLY, &status);
    if (fits_read_key(parmfits, TFLOAT, w0_key, &w0, NULL, &status)) {
	w0 = default_w0;
	status = 0;
    }
    else if (w0 < default_w0)
	cf_if_warning("Requested value of W0 is less than recommended value.");

    if (fits_read_key(parmfits, TFLOAT, wmax_key, &wmax, NULL, &status)) {
        wmax = default_wmax;
        status = 0;
    }   
    else if (wmax > default_wmax)
        cf_if_warning("Requested value of WMAX is greater than "
		      "recommended value.");

    if (fits_read_key(parmfits, TFLOAT, wpc_key, &wpc, NULL, &status)) {
        wpc = default_wpc;
        status = 0;
    }   
    else if (wpc < default_wpc)
        cf_if_warning("Requested value of WPC is less than "
		      "recommended value.");
    FITS_close_file(parmfits, &status);
    cf_verbose(3, "Using these wavelength parameters:");
    cf_verbose(3, "\tW0 = %g, WMAX = %g, WPC = %g", w0, wmax, wpc);

    /*  Compute length of output wavelength array, allocate memory,
     *  and fill array.
     */
    *nout = (long) ((wmax - w0)/wpc + 0.5) + 1L;
    *wave_out  = (float *) cf_calloc(*nout, sizeof(float));
    for (i = 0; i < *nout; i++)
	(*wave_out)[i] = w0 + wpc * (float) i;

    /*
     *  Write WO and WPC to output file header.
     */
    FITS_update_key(header, TFLOAT, "WPC", &wpc, NULL, &status);
    FITS_update_key(header, TFLOAT, "W0", &w0, NULL, &status);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
