/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_screen_limb_angle(fitsfile *infits, long nseconds,
 *                                   unsigned char *timeline_status,
 *                                   float *timeline_limb);
 *
 * Description: Set the screening bit due to limb angle constraint.
 *
 * Arguments:   fitsfile  *infits       Input FITS file pointer
 *              long      nseconds      Number of points in the timeline table
 *              unsigned char  *timeline_status   The status flag array
 * 		float     *timeline_limb The limb angle distance array
 *
 * Calls:
 *
 * Returns:	0 on success
 *
 * History:	10/21/02   1.1   jch    Initial coding
 * 		12/20/02   1.3   jch    Make sure the flags are unsigned char
 *              05/20/03   1.4   rdr    Added call to cf_proc_check
 *              09/10/03   1.5   wvd    Move test of LIMB_SCR to cf_fuv_init.
 *              07/21/04   1.7   wvd    Delete unused variable limb_cor.
 *              02/17/05   1.8   wvd    Place parentheses around assignment
 *                                      used as truth value.
 *              03/30/05   1.9   wvd    Write BRITLIMB and DARKLIMB keywords
 *					to IDF file header.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include "calfuse.h"

int
cf_screen_limb_angle(fitsfile *infits, long nseconds, 
		unsigned char *timeline_status, float *timeline_limb)
{
    char CF_PRGM_ID[] = "cf_screen_limb_angle";
    char CF_VER_NUM[] = "1.9";

    int	  errflg=0, status=0;
    long  i;
    char  file_name[FLEN_VALUE];
    double angle, britlimb, darklimb;
    fitsfile  *scrnfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

    /*
     *  Open the screening parameters file and read the parameters.
     */
    FITS_read_key(infits, TSTRING, "SCRN_CAL", file_name, NULL, &status);
    FITS_open_file(&scrnfits, cf_parm_file(file_name), READONLY, &status);

    FITS_read_key(scrnfits, TDOUBLE, "BRITLIMB", &britlimb, NULL, &status);
    FITS_read_key(scrnfits, TDOUBLE, "DARKLIMB", &darklimb, NULL, &status);
    /*
     *  Write BRITLIMB and DARKLIMB keywords to IDF file header.
     */
    FITS_update_key(infits, TDOUBLE, "BRITLIMB", &britlimb, NULL, &status);
    FITS_update_key(infits, TDOUBLE, "DARKLIMB", &darklimb, NULL, &status);
    /*
     *  If the limb angle is less than the limit, set the limb-angle bit.
     *  Use day or night limit as appropriate.
     */
    for (i = 0; i < nseconds; i++) {
        angle = darklimb;
        if ((timeline_status[i] & TEMPORAL_DAY) == TEMPORAL_DAY) {
	    angle = britlimb;
        }
        if (timeline_limb[i] < angle) {
	     timeline_status[i] |= TEMPORAL_LIMB; 
        }
    }
    FITS_close_file(scrnfits, &status);
    
    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
