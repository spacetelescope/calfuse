/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_ids_dead_time(fitsfile *header, long nseconds,
 *                               float *aic_rate, float *ids_dtc)
 *
 * Description: Computes the weighting factor needed to correct
 *		for the IDS dead time.
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              long      nseconds      The number of timeline values
 *              float     *aic_rate     An array of Active Image Counter
 *                                      values
 *              float     *ids_dtc      Dead-time correction array (returned)
 *
 * Calls:
 *
 * Return:      0  on success
 *
 * History:     10/27/02   1.1    peb   Begin work
 *              11/11/02   1.2    peb   Correct function description and add
 *                                      cf_timestamp after IDS__COR check.
 *              12/09/02   1.4    wvd   Calculate DT correction for each time
 *                                      step, then apply to photons.
 *                                      Set keyword IDS_DEAD.
 *              05/20/03   1.5    rdr   Add proc_check call.
 *		08/01/03   1.8    wvd   Just calculate correction; don't
 *					apply it.  Return ids_dtc array.
 *		08/04/03   1.9    wvd   Convert aic_rate to type short.
 *		11/26/03   1.10   wvd   Convert aic_rate to type float.
 *		02/09/04   1.11   wvd   max_ids_rate depends on instrument
 *					mode. For TTAG data, include time
 *					stamps in aic_rate.
 *		04/07/07   1.12   wvd   Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include "calfuse.h"

int
cf_ids_dead_time(fitsfile *header, long nseconds, float *aic_rate,
	float *ids_dtc)
{
    char CF_PRGM_ID[] = "cf_ids_dead_time";
    char CF_VER_NUM[] = "1.12";

    char elecfile[FLEN_VALUE], instmode[FLEN_VALUE];
    int  errflg=0, status=0;
    long  k;
    float ids_rate, tstamps, ttperiod;
    float mean_ids_dtc = 0.;
    float max_ids_rate = 1.;
    fitsfile *elecfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;
    /*
     *  Read observing mode and interval between time stamps.
     */
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    FITS_read_key(header, TFLOAT, "TTPERIOD", &ttperiod, NULL, &status);
    /*
     *  Compute number of time stamps per second (0 for HIST data).
     */
    if (!strncmp(instmode, "TTAG", 4)) 		/* TTAG MODE */
	if (ttperiod > 1e-8) tstamps = 1./ttperiod;
	else tstamps = 1.;
    else					/* HIST MODE */
	tstamps = 0.;
    /*
     *  Maximum IDS rate depends on instrument mode.
     */
    FITS_read_key(header, TSTRING, "ELEC_CAL", elecfile, NULL, &status);
    FITS_open_file(&elecfits, cf_cal_file(elecfile), READONLY, &status);
    if (!strncmp(instmode, "TTAG", 4)) 
        FITS_read_key(elecfits, TFLOAT, "TTAG_BUS", &max_ids_rate, NULL, &status);
    else
        FITS_read_key(elecfits, TFLOAT, "HIST_BUS", &max_ids_rate, NULL, &status);
    FITS_close_file(elecfits, &status);
    cf_verbose(4, "max_ids_rate = %f", max_ids_rate);
    /*
     *  Calculate the IDS dead-time correction.
     */
    for (k=0; k<nseconds; k++) {
        ids_rate = aic_rate[k] + tstamps;

	if (ids_rate > max_ids_rate)
	    ids_dtc[k] = ids_rate/max_ids_rate;
	else
	    ids_dtc[k] = 1.0;
        mean_ids_dtc += ids_dtc[k];
    }
    /*
     *  Write mean IDS dead-time correction to file header.
     */
    mean_ids_dtc /= nseconds;
    cf_verbose(2, "Mean IDS dead-time correction: %f", mean_ids_dtc);
    FITS_update_key(header, TFLOAT, "IDS_DEAD", &mean_ids_dtc, NULL, &status);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
