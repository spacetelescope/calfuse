/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:	cf_screen_high_voltage(fitsfile *infits, long nseconds,
 *              unsigned char *timeline_status, short *timeline_hv);
 *
 * Description: Set the screening bit if detector voltage is low.
 *
 * Arguments:   fitsfile  *infits       Input FITS file pointer
 * 		long      nseconds      Number of points in the timeline table
 * 		unsigned char  *timeline_status   The status flag array
 * 		short     *timeline_hv  The high voltage array
 *
 * Calls:
 *
 * Returns:	0 on success
 *
 * History:	10/23/02   1.1   jch    Initial coding
 * 		12/20/02   1.3   jch    Make sure the flags are "unsigned" char
 *              01/24/03   1.4   rdr    Lower the voltage threshold required 
 *                                      to flag a time to 95% of the full 
 *                                      voltage
 * 		02/11/03   1.6   jch    Skip HV which is negative
 *              05/20/03   1.7   rdr    Add call to cf_proc_check
 *              05/22/03   1.8   wvd    Implement cf_verbose
 *              06/11/03   1.9   wvd    Change HV array to type short.
 *              06/13/03   1.10  wvd    Change calfusettag.h to calfuse.h
 *              09/10/03   1.11  wvd    Return errflg from cf_proc_check.
 *              07/21/04   1.12  wvd    Delete unused variable volt_cor.
 *              02/17/05   1.13  wvd    Place parentheses around assignment
 *					used as truth value.
 *              03/18/05   1.14  wvd    HV can be 0, so test for HV > -1
 *              03/24/05   1.15  wvd    Read HVGOODLM and HVBADLIM from 
 *					VOLT_CAL file.  If HV/FULLV drops
 *					below HVGOODLM, issue warning.
 *					If it falls below HVBADLM, set
 *					TEMPORAL_HV flag.
 *              06/22/06   1.16  wvd    Don't set EXP_STAT keyword.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include "calfuse.h"

int
cf_screen_high_voltage(fitsfile *infits, long nseconds,
	       unsigned char *timeline_status, short *timeline_hv)
{
    char CF_PRGM_ID[] = "cf_screen_high_voltage";
    char CF_VER_NUM[] = "1.16";

    char  file_name[FLEN_VALUE], mjd_str[FLEN_VALUE], full_str[FLEN_VALUE];
    char  datestr[FLEN_VALUE], comment[FLEN_COMMENT], instmode[FLEN_VALUE];
    int	  errflg=0, status=0, timeref, hvflag=FALSE, n, fullv;
    long  i;
    float expstart, hvfrac, hvgood, hvbad, mjdv;
    fitsfile *scrnfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

    /* Read exposure start time and instrument mode. */
    FITS_read_key(infits, TFLOAT, "EXPSTART", &expstart, NULL, &status);
    FITS_read_key(infits, TSTRING, "INSTMODE", instmode, NULL, &status);
    /*
     *  Open the high-voltage calibration file.
     */
    FITS_read_key(infits, TSTRING, "VOLT_CAL", file_name, NULL, &status);
    FITS_open_file(&scrnfits, cf_cal_file(file_name), READONLY, &status);
    /*
     *  Read the limits of the good and questionable voltage regimes.
     */
    FITS_read_key(scrnfits, TFLOAT, "HVGOODLM", &hvgood, NULL, &status);
    FITS_read_key(scrnfits, TFLOAT, "HVBADLIM", &hvbad,  NULL, &status);
    /*
     *  Read full-voltage level in use at time of observation.
     */
    n = 0;
    do {
	n++;
	sprintf(mjd_str, "MJD%d", n);
	FITS_read_key(scrnfits, TFLOAT, mjd_str, &mjdv, NULL, &status);
    } while (expstart > mjdv);
    n--;

    sprintf(full_str, "FULL%d", n);
    FITS_read_key(scrnfits, TINT, full_str, &fullv, NULL, &status);
    FITS_close_file(scrnfits, &status);

    /*
     *  Check for low values of the detector voltage.
     *  Set the high voltage bit if timeline_hv/fullv < hvbad.
     */
        cf_verbose(2, "Threshold voltage = %d", cf_nint(hvbad * fullv)) ; 

    for (i = 0; i < nseconds; i++) if (timeline_hv[i] > -1) {
	hvfrac = timeline_hv[i] / (float) fullv;
	if (hvfrac < hvgood) {
	    if (hvfrac < hvbad) timeline_status[i] |= TEMPORAL_HV; 
	    else hvflag = TRUE;
	}
    }

    /*
     *  If the detector voltage ever falls between the good and bad
     *  levels, issue a stern warning to the user.
     */
    if (hvflag) {
	/* Write warning to trailer file. */
	if (!strncmp(instmode, "HIST", 4))
	cf_if_warning("Detector voltage is < 90%% of optimal value.  "
	"Check photometry and wavelength scale.");
	else
	cf_if_warning("Detector voltage is < 90%% of optimal value.  "
	"Check pulse-height distribution.");

	/* Write warning to IDF file header. */
	FITS_write_comment(infits, " ", &status);
	FITS_write_comment(infits,
	"Detector voltage is less than 90%% of optimal value.", &status);
	if (!strncmp(instmode, "HIST", 4)) {
	FITS_write_comment(infits,
	"Carefully examine these data for photometric errors", &status);
	FITS_write_comment(infits, "and walk effects.", &status);
	}
	else {
	FITS_write_comment(infits,
	"Carefully examine the pulse height distribution of photons", &status);
	FITS_write_comment(infits,
	"in the target aperture before accepting these data.", &status);
	}
	fits_get_system_time(datestr, &timeref, &status);
	sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
	FITS_write_comment(infits, comment, &status);
	FITS_write_comment(infits, "  ", &status);
    }

    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
