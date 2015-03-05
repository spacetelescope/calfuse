/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_set_user_gtis(fitsfile *infits, long nseconds,
 *                                  float *timeline_times,
 *                                  unsigned char *timeline_status);
 *
 * Description: Apply user-defined good time intervals to timeline table.
 *
 * Arguments:   fitsfile  *infits       Input FITS file pointer
 *              long      nseconds      Number of points in the timeline table
 * 		float     *timeline_times Array of times
 *              unsigned char  *timeline_status   The status flag array
 *
 * Calls:
 *
 * Returns:	0 on success
 *
 * History:	09.10.03   1.1   bjg    Initial coding
 * 		
 *              09.12.03   1.2   bjg    Bug fix. When nusergti was zero or less
 *                                      this function set all the photons to 
 *                                      bad. Now it doesn't change anything.
 *              04.20.04   1.3   bjg    Remove unused variables
 *                                      Change format to match arg type in
 *                                      sprintf.
 *
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include "calfuse.h"

int
cf_set_user_gtis(fitsfile *infits, long nseconds, 
	 float *timeline_times,	unsigned char *timeline_status)
{
    char CF_PRGM_ID[] = "cf_set_user_gtis";
    char CF_VER_NUM[] = "1.3";

    int	  errflg=0, status=0, nusergti;
    long  i,n;
    char  file_name[FLEN_VALUE];
    float limit1,limit2;
    char keywd1[FLEN_KEYWORD],keywd2[FLEN_KEYWORD];
    fitsfile  *scrnfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(infits, CF_PRGM_ID) )) return errflg;

    /*
     *  Open the screening parameters file and read the good time intervals.
     */
    FITS_read_key(infits, TSTRING, "SCRN_CAL", file_name, NULL, &status);
    FITS_open_file(&scrnfits, cf_parm_file(file_name), READONLY, &status);

    FITS_read_key(scrnfits, TINT, "NUSERGTI", &nusergti, NULL, &status);
    
    if (nusergti>0){
    
        i=0;
    
        for (n=1;n<=nusergti;n++){
           sprintf(keywd1, "GTIBEG%02ld",n);
           FITS_read_key(scrnfits, TFLOAT, keywd1, &limit1, NULL, &status);
           sprintf(keywd2, "GTIEND%02ld",n);
           FITS_read_key(scrnfits, TFLOAT, keywd2, &limit2, NULL, &status);
        
        
           while ((i<nseconds)&&(timeline_times[i]<limit1)) {
              timeline_status[i] |= TEMPORAL_USER;
              i++;
           }
        
           while ((i<nseconds)&&(timeline_times[i]<limit2)) {
              i++;
           }
    
        }


        while (i<nseconds){
             timeline_status[i] |= TEMPORAL_USER;
             i++;
        }
    
    }
    
    FITS_close_file(scrnfits, &status);
    
    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
