/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_screen_pulse_height(fitsfile *infits, long nevents,
 *                              unsigned char *photon_ph, 
 *				unsigned char *photon_locflag);
 *
 * Description: Set the screening bit according to PHA limits.
 *
 * Arguments:   fitfile  *infits        Input FITS file pointer
 *              long     nevents        Number of points in the photon list
 *              unsigned char     *photon_ph     The photon pulse height array
 *              unsigned char     *photon_locflag The photon location array
 *
 * Returns:	0 on success
 *
 * History:	11/01/02   1.1   jch    Initial coding
 * 		12/20/02   1.3   jch    Make sure the flags are "unsigned" char
 *              02/14/03   1.4   rdr    Update phalow and phahigh keywords
 *                                       in the IDF
 *              05/20/03   1.5   rdr    Add call to cf_proc_check
 *              09/10/03   1.6   wvd    Write NBADPHA as type LONG.
 *              10/02/03   1.7   wvd    Move PHA bit from TEMPORAL to 
 *					LOCATION flags.
 *              05/04/04   1.9   bjg    Cosmetic change to prevent 
 *                                      compiler warning with gcc -Wall
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include "calfuse.h"

int
cf_screen_pulse_height(fitsfile *infits, long nevents, unsigned char *photon_ph,
		       unsigned char *photon_locflag)
{
    char CF_PRGM_ID[] = "cf_screen_pulse_height";
    char CF_VER_NUM[] = "1.9";

    char  scrnfile[FLEN_VALUE];
    unsigned char  phalow, phahigh;
    int	  errflg=0, status=0;
    long  i, nbadpha=0;
    fitsfile *scrnfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

     /*
     *  Open the screening parameters file and read the parameters.
     */
    FITS_read_key(infits, TSTRING, "SCRN_CAL", scrnfile, NULL, &status);
    FITS_open_file(&scrnfits, cf_parm_file(scrnfile), READONLY, &status);
    FITS_read_key(scrnfits, TBYTE, "PHALOW", &phalow, NULL, &status);
    FITS_read_key(scrnfits, TBYTE, "PHAHIGH", &phahigh, NULL, &status);
    FITS_close_file(scrnfits, &status);
    /*
     *  Set the pulse-height bit if pulse height is outside limits.
     */
    for (i = 0; i < nevents; i++) {
	if (photon_ph[i] < phalow || photon_ph[i] > phahigh) {
	    photon_locflag[i] |= LOCATION_PHA; 
	    nbadpha++;
	}
    }
    FITS_update_key(infits, TLONG, "NBADPHA", &nbadpha, NULL, &status);
    FITS_update_key(infits, TBYTE, "PHALOW",  &phalow,  NULL, &status);
    FITS_update_key(infits, TBYTE, "PHAHIGH", &phahigh, NULL, &status);
    
    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
