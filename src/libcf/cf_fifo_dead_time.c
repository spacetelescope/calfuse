/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_fifo_dead_time(fitsfile *header, long nevents,
 *		float *ptime, long nseconds, float *ttime, float *aic_rate,
 *		float *ids_dtc)
 *
 * Description: Searches high-count-rate observations for data drop-outs
 *		due to FIFO overflows.  If found, calculates time-dependent
 *		scale factor to correct for lost photon events.
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              long      nevents       Number of photon events in the file
 *              float     *ptime        Detection time for each photon
 *              long      nseconds      The number of timeline values
 *              float     *ttime        Tabulated times in the timeline
 *              float     *aic_rate     Active Image Counter array
 *              float     *ids_dtc      IDS dead-time correction array
 *					(modified)
 *
 * Calls:
 *
 * Return:      0  on success
 *
 * History of cf_screen_fifo_overflow:
 *
 *		02/10/04   1.1    wvd   Begin work
 *              02/17/05   1.2    wvd   Place parentheses around assignment
 *                                      used as truth value.
 *              11/09/05   1.3    wvd   In the current scheme, the rate_lif
 *					and rate_sic arrays contain values
 *					from the housekeeping file and do
 *					not reflect data loss due to FIFO
 *					overflows.  Now we calculate the
 *					count-rate array directly from the
 *					photon time array.
 *              02/23/06   1.4    wvd   Change length of data_rate array
 *					from nevents to nseconds.
 *
 * History:	12/29/06   1.1    wvd	Derived from cf_screen_fifo_overflow.
 *              01/15/07   1.2    wvd   When count rate is zero, correct one
 *					second before to two seconds after.
 *              02/08/07   1.3    wvd   If TTPERIOD = 0, set it to 1.
 *              04/07/07   1.4    wvd   Clean up compiler warnings.
 *
 ****************************************************************************/

#include <stdlib.h>
#include "calfuse.h"

int
compute_count_rate(long nevents, float *ptime, long nseconds,
    float *ttime, long **data_rate)
{
    float delta_max = 1.0 + FRAME_TOLERANCE;
    long j,k;

    *data_rate = (long *) cf_calloc(nseconds, sizeof(long));
    for (j=k=0; j<nevents; j++) {
        while(ttime[k+1]-FRAME_TOLERANCE < ptime[j] && k+1 < nseconds) k++;
        if ((ptime[j] - ttime[k]) < delta_max) (*data_rate)[k]++;
    }
    return 0;
}

int
cf_fifo_dead_time(fitsfile *header, long nevents, float *ptime,
    long nseconds, float *ttime, float *aic_rate, float *ids_dtc)
{
    char CF_PRGM_ID[] = "cf_fifo_dead_time";
    char CF_VER_NUM[] = "1.4";

    char elecfile[FLEN_VALUE];
    int  errflg=0, status=0;
    int  fifo_drain_rate;
    int  *flag_array, found_drop_out = FALSE;
    long *data_rate, k;
    float mean_ids_rate, ttperiod;
    float mean_ids_dtc;
    fitsfile *elecfits;

    /* Exit if data were obtained in HIST mode. */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* Read IDS drain rate from calibration file. */
    FITS_read_key(header, TSTRING, "ELEC_CAL", elecfile, NULL, &status);
    FITS_open_file(&elecfits, cf_cal_file(elecfile), READONLY, &status);
    FITS_read_key(elecfits, TINT, "FIFO_DRN", &fifo_drain_rate, NULL, &status);
    FITS_close_file(elecfits, &status);
    cf_verbose(4, "fifo_drain_rate = %d", fifo_drain_rate);

    /* Read interval between time stamps. */
    FITS_read_key(header, TFLOAT, "TTPERIOD", &ttperiod, NULL, &status);
    if (ttperiod < 1E-6) ttperiod = 1.0;

    /* Compute mean IDS count rate. */
    mean_ids_rate = 0.;
    for (k=0; k<nseconds; k++) 
        mean_ids_rate += aic_rate[k];
    mean_ids_rate /= nseconds;
    mean_ids_rate += 1./ttperiod;
    cf_verbose(4, "mean_ids_rate = %f", mean_ids_rate);

    /* If mean IDS count rate > 3200, look for data drop-outs. */
    if (mean_ids_rate > 3200 && nseconds > 4) {
	compute_count_rate(nevents, ptime, nseconds, ttime, &data_rate);
	flag_array = (int *) cf_calloc(nseconds, sizeof(int));

        /* If count rate goes to zero, set flag_array.  Where possible,
	   set flags from one second before to two seconds after. */
	for (k=0; k<1; k++) if (data_rate[k] == 0) 
	    flag_array[k] = flag_array[k+1] = flag_array[k+2] = TRUE;

	for ( ; k<nseconds-2; k++) if (data_rate[k] == 0) 
	    flag_array[k-1] = flag_array[k] = flag_array[k+1] = 
		flag_array[k+2] = TRUE;

	for ( ; k<nseconds-1; k++) if (data_rate[k] == 0) 
	    flag_array[k-1] = flag_array[k] = flag_array[k+1] = TRUE;

	for ( ; k<nseconds; k++) if (data_rate[k] == 0) 
	    flag_array[k-1] = flag_array[k] = TRUE;

	for (k=0; k<nseconds; k++) if (data_rate[k] == 0) {
	   found_drop_out = TRUE;
	   break;
	}

	/* If there are any drop-outs, modify IDS dead-time correction. */
	if (found_drop_out)
	     for ( ; k < nseconds; k++) if (flag_array[k])
		ids_dtc[k] = (aic_rate[k] - 1./ttperiod) / fifo_drain_rate;

	free(data_rate);
	free(flag_array);
    }

    /*
     *  Calculate mean dead-time correction.  Write to file header.
     */
    if (found_drop_out) {
	mean_ids_dtc = 0;
	for (k = 0; k < nseconds; k++) mean_ids_dtc += ids_dtc[k];
	mean_ids_dtc /= nseconds;
	cf_verbose(1, "IDS FIFO overflow detected.");
        cf_verbose(2, "Modified mean IDS dead-time correction: %f",
	    mean_ids_dtc);
        FITS_update_key(header, TFLOAT, "IDS_DEAD", &mean_ids_dtc, NULL,
	    &status);
    }

    /*
     *  Clean up and go home.
     */
    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
