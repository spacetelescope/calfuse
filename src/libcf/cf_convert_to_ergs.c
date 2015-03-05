/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_convert_to_ergs(fitsfile *header, long nevents,
 *              float *weight, float *ergcm2, unsigned char *channel,
 *		float *lambda);
 *
 * Description: Convert counts to ergs.
 *
 * Arguments:   *header                 Input FITS file pointer
 *              nevents                 number of points in the photon list
 *              *weight			weight of each event
 * 		*ergcm2			flux (returned)
 *              *channel                channel number in the photon list
 *              *lambda                 wavelength in the photon list
 *
 * Returns:    	0 (int) upon successful completion.
 *
 * History:    	01/02/03   1.0   jch    Initial coding
 *		02/10/03   1.1   wvd    Install
 *		02/10/03   1.2   wvd    Replace FLUX_CAL with AEFF_CAL
 *                                      Define HC in calfusettag.h
 *                                      Don't sort photons by wavelength
 *              03/04/03   1.3   peb    Fixed sprintf argument formating,
 *                                      deleted unused variable, localized
 *                                      scope of some variables, made sure
 *                                      that memory is allocated for aeff,
 *                                      and fixed memory leak with aeff,
 *                                      aeff1, and aeff2.
 *		03/11/03   1.4   wvd    Changed channel to type char
 *		03/12/03   1.5   wvd    Fixed bug in calculation of dl.
 *					If unable to calculate flux,
 *					set channel[k] = 0.
 *					Made aeff, aeff1, aeff2 doubles
 *					to match aeff*fit files.
 *		05/07/03   1.6   wvd	Change ERGCM2S to ERGCM2 throughout.
 *					Don't divide by EXPTIME.
 *					Use cf_verbose throughout.
 *		05/16/03   1.7   wvd	Update version number.
 *              05/20/03   1.8   rdr    Add call to cf_proc_check
 *              06/11/03   1.9   wvd    Pass datatype to cf_read_col.
 *					Change calfusettag.h to calfuse.h
 *              08/25/03   1.10  wvd    Change coltype from string to int
 *					in cf_read_col.
 *              09/17/03   1.11  wvd    Change aeff arrays to type float
 *					to match AEFF_CAL files.
 *              11/05/03   1.12  wvd    Change channel to unsigned char
 *              04/27/04   1.13  wvd    Print out relative weighting of
 *					the two effective-area files.
 *              11/18/05   1.14  wvd    Change aeff arrays to type double
 *					to match AEFF_CAL files.
 *              04/07/07   1.15  wvd	Clean up compiler warnings.
 *              04/07/07   1.16  wvd	Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_convert_to_ergs";
static char CF_VER_NUM[] = "1.16";


int
cf_convert_to_ergs(fitsfile *header, long nevents, float *weight, 
		   float *ergcm2, unsigned char *channel, float *lambda)
{
    char        aeff1file[FLEN_VALUE], aeff2file[FLEN_VALUE];
    int		ap, errflg=0, interp=0, status=0;
    float 	sig1=0, sig2=0;
    fitsfile    *aeff1fits, *aeff2fits;

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* Read the AEFF1CAL and AEFF2CAL keywords. */
    FITS_read_key(header, TSTRING, "AEFF1CAL", aeff1file, NULL, &status);
    FITS_read_key(header, TSTRING, "AEFF2CAL", aeff2file, NULL, &status);

    /* Open the calibration files. */
    FITS_open_file(&aeff1fits, cf_cal_file(aeff1file), READONLY, &status);

    /* Determine whether we need to interpolate getween AEFF_CAL files. */
    if (strcmp(aeff1file, aeff2file)) {
	float dt1, dt2, dt3, a1date, a2date, expstart, expend, expmjd;
        interp = 1;
        FITS_open_file(&aeff2fits, cf_cal_file(aeff2file), READONLY, &status);
        FITS_read_key(aeff1fits, TFLOAT, "EFFMJD", &a1date, NULL, &status);
        FITS_read_key(aeff2fits, TFLOAT, "EFFMJD", &a2date, NULL, &status);

	/* Linearly interpolate the aeff according to exposure time */
        FITS_read_key(header, TFLOAT, "EXPSTART", &expstart, NULL, &status);
        FITS_read_key(header, TFLOAT, "EXPEND", &expend, NULL, &status);
	expmjd = (expstart + expend) / 2.;
        dt1=a2date-a1date;
        dt2=a2date-expmjd;
        dt3=expmjd-a1date;
        if (dt1 < 0.01) {
            cf_if_error("Calibration file dates are identical: "
	    "cal1=%10.5f obs=%10.5f cal2=%10.5f", a1date, expmjd, a2date);
        }
        sig1=dt2/dt1;
        sig2=dt3/dt1;
	cf_verbose(3, "Interpolation: %0.2f x %s, %0.2f x %s", 
		sig1, aeff1file, sig2, aeff2file);
    }
    else
	cf_verbose(3, "No interpolation.  Using Aeff file %s", aeff1file);
	
    /* Step through the apertures, skipping aperture 4 (pinhole) */
    for (ap = 1; ap < 8; ap++) {
	long jmax=0,	/* size of wave and aeff arrays */
    	    j,		/* index through wave and aeff */
	    k;		/* index through photon array */

	double *aeff;
	float *wavelength, dl, w0;
	
        if (ap == 4)
	    continue;

	if (interp == 1) {
	    /* Interpolate the calibrated aeff. */
	    double *aeff1, *aeff2;

            FITS_movabs_hdu(aeff1fits, ap+1, NULL, &status);
	    jmax = cf_read_col(aeff1fits, TDOUBLE, "AREA", (void **) &aeff1);

            FITS_movabs_hdu(aeff2fits, ap+1, NULL, &status);
	    jmax = cf_read_col(aeff2fits, TDOUBLE, "AREA", (void **) &aeff2);

	    aeff = (double *) cf_calloc(jmax, sizeof(double));

	    for (j = 0; j < jmax; j++) {
                aeff[j] = sig1 * aeff1[j] + sig2 * aeff2[j];
	    }
	    free(aeff1);
	    free(aeff2);
	} else {
            FITS_movabs_hdu(aeff1fits, ap+1, NULL, &status);
	    jmax = cf_read_col(aeff1fits, TDOUBLE, "AREA", (void **) &aeff);
	}

	/* Read wavelength array in AEFF_CAL file */
	jmax = cf_read_col(aeff1fits, TFLOAT, "WAVE", (void **) &wavelength);

	/* Compute translation from wavelength to pixel within AEFF_CAL file. */
	w0 = wavelength[0];
	dl = (wavelength[jmax-1] - wavelength[0]) / (jmax - 1);

	/* Go through the photon list */
        for (k = 0; k < nevents; k++) {
            if (channel[k] == ap) {
		long j = (lambda[k] - w0) / dl + 0.5;
		if (j >= 0L && j < jmax)
		    ergcm2[k] = weight[k]*HC/lambda[k]/aeff[j];
		else
		    channel[k] = 0;
            }
        }
	free(wavelength);
	free(aeff);
    }

    /* Clean up. */
    FITS_close_file(aeff1fits, &status);
    if (interp == 1)
        FITS_close_file(aeff2fits, &status);

    /* Update processing flags. */
    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");
    return (status);
}
