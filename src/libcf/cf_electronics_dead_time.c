/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_electronics_dead_time(fitsfile *header, long nseconds,
 *                                      float *fec_rate, float *elec_dtc)
 *
 * Description: Scales the weighting factor to correct for the electronics
 *              dead time.
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              long      nseconds      The number of timeline values
 *              float     *fec_rate     An array of Front End Counter (FEC)
 *                                      rates
 *              float     *elec_dtc     Dead-time correction array (returned)
 * Calls:
 *
 * Return:      0  on success
 *
 * History:     10/27/02   1.1    peb   Begin work
 *              11/11/02   1.2    peb   Corrected function description and
 *                                      added cf_timestamp after ELEC_COR
 *                                      check.
 *		12/06/02   1.4    wvd   Calculate DT correction for each time
 *					step, then apply to photons.
 *					Set keyword DET_DEAD.
 *              05/20/03   1.5    rdr   Added proc_check call
 *              08/01/03   1.7    wvd   Just calculate correction; don't
 *                                      apply it.  Return elec_dtc array.
 *              08/04/03   1.8    wvd   Convert fec_rate to type short.
 *              11/26/03   1.9    wvd   Change fec_rate to type float.
 *              02/12/04   1.10   wvd   In verbose mode, write mean DTC.
 *              04/07/07   1.11   wvd   Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include <math.h>
#include "calfuse.h"

int
cf_electronics_dead_time(fitsfile *header, long nseconds, 
			 float *fec_rate, float *elec_dtc)
{
    char CF_PRGM_ID[] = "cf_electronics_dead_time";
    char CF_VER_NUM[] = "1.11";

    char  elecfile[FLEN_VALUE]={"\0"}, detector[FLEN_VALUE]={"\0"};
    char  keyword[FLEN_KEYWORD]={"\0"};
    int   errflg=0, status=0;
    long  k;
    float abort, clock, state;
    float mean_elec_dtc = 0.;
    fitsfile *elecfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    FITS_read_key(header, TSTRING, "DETECTOR", detector, NULL, &status);
    /*
     *  Get the electronics dead-time constants.
     */
    FITS_read_key(header, TSTRING, "ELEC_CAL", elecfile, NULL, &status);
    FITS_open_file(&elecfits, cf_cal_file(elecfile), READONLY, &status);
    sprintf(keyword, "ABORT_%s", detector);
    FITS_read_key(elecfits, TFLOAT, keyword, &abort, NULL, &status);
    sprintf(keyword, "CLOCK_%s", detector);
    FITS_read_key(elecfits, TFLOAT, keyword, &clock, NULL, &status);
    sprintf(keyword, "STATE_%s", detector);
    FITS_read_key(elecfits, TFLOAT, keyword, &state, NULL, &status);
    FITS_close_file(elecfits, &status);
    /*
     *  Calculate dead-time correction.
     */
    for (k=0; k<nseconds; k++) {
        float xa, x, x2, x4, x6, x10;

        xa = exp(-fec_rate[k]*abort);
        x  = exp(-fec_rate[k]*clock);
        x2 = x*x;
        x4 = x2*x2;
        x6 = x4*x2;
        x10 = x6*x4;

        elec_dtc[k] = (1. + fec_rate[k]*state + 4.*(x2-x) + x6 - x10)/xa;
        if (elec_dtc[k] < 1.0) elec_dtc[k] = 1.0;
	mean_elec_dtc += elec_dtc[k];
    }
    /*
     *  Calculate mean dead-time correction.  Write to file header.
     */
    mean_elec_dtc /= nseconds;
    cf_verbose(2, "Mean detector electronics dead-time correction: %f", mean_elec_dtc);
    FITS_update_key(header, TFLOAT, "DET_DEAD", &mean_elec_dtc, NULL, &status);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
