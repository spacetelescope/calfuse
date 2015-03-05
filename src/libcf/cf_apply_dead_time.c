/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_apply_dead_time(fitsfile *header, long nevents, float *time,
 *                      float *weight, long nseconds, float *timeline, 
 *			float *ids_dtc, float *elec_dtc)
 *
 * Description: Applies dead-time correction to all photons.  For HIST data,
 *		uses mean of dead-time arrays.  Writes mean dead-time
 *		values to file header.  Issue warning if DET_DEAD > 1.5.
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              long      nevents       The number of events
 *              float     *time         An array of event times
 *              float     *weight       An array of event weights
 *              long      nseconds      The number of timeline values
 *              float     *timeline     An array of timeline times
 *              float     *aic_rate     Active Image Counter rate
 *              float     *ids_dtc      Output of cf_ids_dead_time
 *              float     *elec_dtc     Output of cf_elec_dead_time
 *
 * Calls:
 *
 * Return:      0  on success
 *
 * History:     08/01/03   1.1    wvd   Begin work
 *		08/04/03   1.2    wvd	Convert aic_rate to type short.
 *		11/25/03   1.3    wvd	For HIST data, require that
 *					mean_tot_dtc be >= 1.
 *		11/26/03   1.4    wvd	Change aic_rate to type float.
 *		02/09/04   1.5    wvd	Re-write with new logic.
 *              04/09/04   1.6    bjg   Include stdlib.h
 *		11/07/05   1.7    wvd	Add a test to catch bad DTC values.
 *		03/10/06   1.8    wvd	Clean up i/o.
 *		12/29/06   1.9    wvd	Check DET_DEAD rather than TOT_DEAD.
 *		04/07/07   1.10   wvd	Clean up compiler errors.
 *
 ****************************************************************************/
#include <stdlib.h>
#include <string.h>
#include "calfuse.h"

int
cf_apply_dead_time(fitsfile *header, long nevents, float *time, float *weight,
	long nseconds, float *timeline, float *ids_dtc, float *elec_dtc)
{
    char CF_PRGM_ID[] = "cf_apply_dead_time";
    char CF_VER_NUM[] = "1.10";

    char comment[FLEN_CARD], datestr[FLEN_CARD], instmode[FLEN_VALUE];
    char hkexists[FLEN_VALUE];
    double ctimeb, ctimee, eng_time;
    int  errflg=0, status=0, timeref;
    int  det_flag=FALSE, ids_flag=FALSE;
    long j, k;
    float *tot_dtc, det_dead, ids_dead, mean_tot_dtc=0.;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* Read header keywords.  Set flag if DET_DEAD > 1.5 */
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    FITS_read_key(header, TSTRING, "HKEXISTS", hkexists, NULL, &status);
    FITS_read_key(header, TFLOAT, "DET_DEAD", &det_dead, NULL, &status);
    if (det_dead > 1.5) det_flag = TRUE;

    /*
     *  If there is no housekeeping file, AND the engineering snapshot times
     *  are bad, AND the mean value of any DTC array is > 1.5, set the 
     *  combined DTC array to unity and issue a warning.
     */
    if(!strncasecmp(hkexists, "N", 1)) {
	FITS_read_key(header, TDOUBLE, "CTIME_B", &ctimeb, NULL, &status);
        FITS_read_key(header, TDOUBLE, "CTIME_E", &ctimee, NULL, &status);
        if ((eng_time = (ctimee-ctimeb) * 86400.) < 1.) {
            FITS_read_key(header, TFLOAT, "IDS_DEAD", &ids_dead, NULL, &status);
	    if (ids_dead > 1.5) ids_flag = TRUE;
	}
    }

    /*
     *  Compute total dead-time correction and its mean.
     */
    tot_dtc = (float *) cf_malloc(sizeof(float) * nseconds);
    if (det_flag || ids_flag) {
	for (k = 0; k < nseconds; k++) tot_dtc[k] = 1.0;
	mean_tot_dtc = 1.0;
    }
    else {
	for (k = 0; k < nseconds; k++) {
	    tot_dtc[k] = elec_dtc[k] * ids_dtc[k];
            mean_tot_dtc += tot_dtc[k];
	}
        mean_tot_dtc /= nseconds;
    }

    /* In TTAG mode, scale weights by tot_dtc as a function of time. */
    if (!strncmp(instmode, "TTAG", 4)) for (j=k=0; j<nevents; j++) {
        while(timeline[k+1]-FRAME_TOLERANCE < time[j] && k+1 < nseconds)
            k++;
        weight[j] *= tot_dtc[k];
    }

    /* In HIST mode, scale all weights by mean dead-time correction. */
    else if (mean_tot_dtc > 1.)
	for (j = 0; j < nevents; j++) 
            weight[j] *= mean_tot_dtc;
    else mean_tot_dtc = 1.;

    cf_verbose(2, "Mean total dead-time correction: %f", mean_tot_dtc);
    FITS_update_key(header, TFLOAT, "TOT_DEAD", &mean_tot_dtc, NULL, &status);

    if (det_flag || ids_flag) {
        FITS_write_comment(header, "  ", &status);
	if (ids_flag) {
	    sprintf(comment,
		"IDS_DEAD > 1.5.  Setting dead-time correction to unity.");
            cf_if_warning(comment);
            FITS_write_comment(header, comment, &status);
	}
	if (det_flag) {
	    sprintf(comment,
		"DET_DEAD > 1.5.  Setting dead-time correction to unity.");
            cf_if_warning(comment);
            FITS_write_comment(header, comment, &status);
	}
        fits_get_system_time(datestr, &timeref, &status);
        sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
        FITS_write_comment(header, comment, &status);
        FITS_write_comment(header, "  ", &status);
    }

    free (tot_dtc);
    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
