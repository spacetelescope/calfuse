/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_modify_hist_pha (fitsfile *header, long nevents,
 *                          unsigned char *pha, unsigned char *channel)
 *
 * Description: For histogram data, set pulse heights of all photons to the
 *              most likely value, a function of aperture and obs date.
 *
 * Arguments:   fitsfile       *header      Input FITS file pointer
 *              long           nevents      Number of photon events
 *              unsigned char  *pha	    Pulse-height array (output)
 *              unsigned char  *channel     Channel assignment (input)
 *
 * Calls:
 *
 * Returns:     0 on success
 *
 * History:     03/01/05   1.1   wvd    Initial coding
 *		03/16/05   1.2   wvd	Read PHA array as type TBYTE, not TINT.
 *		04/28/05   1.3   wvd	Fix bug in loop through MJD array.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

int
cf_modify_hist_pha (fitsfile *header, long nevents, unsigned char *pha,
	unsigned char *channel)
{
    char CF_PRGM_ID[] = "cf_modify_hist_pha";
    char CF_VER_NUM[] = "1.3";

    char		phahfile[FLEN_FILENAME];
    unsigned char	pha_chan[8], *pha_mjd=NULL;
    int			errflg=0, status=0;
    long		i, j, jmax;
    double		expstart, *mjd=NULL;
    fitsfile		*phahfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* If data were not taken in HIST mode, exit now. */
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* Read header keywords. */
    FITS_read_key(header, TDOUBLE, "EXPSTART", &expstart, NULL, &status);

    /* Open the pulse-height calibration file. */
    FITS_read_key(header, TSTRING, "PHAH_CAL", phahfile, NULL, &status);
    cf_verbose(3, "Pulse-height calibration file = %s", phahfile);
    FITS_open_file(&phahfits, cf_cal_file(phahfile), READONLY, &status);

    /*
     * Loop through all channels, reading expected pulse heights from PHAH_CAL.
     * Use the values appropriate for the date of this exposure.
     */
    pha_chan[0] = 20;		/* Default for events not in a channel */
    for (i = 1; i < 8; i++) {
	if (i == 4) continue;
        FITS_movabs_hdu(phahfits, i+1, NULL, &status);
        jmax = cf_read_col(phahfits, TDOUBLE, "MJD", (void *) &mjd);
        jmax = cf_read_col(phahfits, TBYTE,   "PHA", (void *) &pha_mjd);
	j = 0;
	while (j < jmax-1 && mjd[j+1] < expstart) j++;
	pha_chan[i] = pha_mjd[j];
	cf_verbose(3, "channel = %d, j = %d, pha = %d", i, j, pha_chan[i]);
    }

    /* Close the pulse-height calibration file. */
    FITS_close_file(phahfits, &status);

    /* Set pulse height for each photon according to its channel number. */
    for (j = 0; j < nevents; j++)
	pha[j] = pha_chan[(int) channel[j]];

    free(mjd);
    free(pha_mjd);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");

    return status;
}
