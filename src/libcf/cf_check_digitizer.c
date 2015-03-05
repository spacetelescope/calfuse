/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_check_digitizer(fitsfile *infits);
 *
 * Description: Check the digitizer values against a reference file.
 *
 * Arguments:   fitsfile  *infits       Input FITS file pointer
 *
 * Calls:       
 *
 * Returns:     0 on success
 *
 * History:     01/16/03   1.1   jch    Initial coding
 *		01/28/03   1.2   wvd    Check only values for this detector.
 *		02/29/03   1.3   rch    Correct error in making comparisons.
 *              01/29/03   1.4   rdr    Include math.h
 *              03/04/03   1.5   peb    Remove unused variables and header
 *              05/20/03   1.6   rdr    Add call to cf_proc_check
 *		07/29/03   1.8   wvd	If cf_proc_check fails, return errflg.
 *		01/28/05   1.9   wvd	On error, set EXP_STAT flag to -999 and
 *					issue a warning rather than an error.
 *		03/11/05   1.10  wvd	On error, write a warning to the file
 *					header using our standard format.
 *		03/30/07   1.11  wvd	On error, set EXP_STAT flag to -2,
 *					rather than -999.
 *
 ****************************************************************************/

#include <stdio.h>
#include <math.h>
#include "calfuse.h"

#define N 32

int
cf_check_digitizer(fitsfile *infits)
{
    char 	CF_PRGM_ID[] = "cf_check_digitizer";
    char 	CF_VER_NUM[] = "1.11";

    char  	comment[FLEN_COMMENT], datestr[FLEN_CARD], 
		file_name[FLEN_VALUE], detector[FLEN_VALUE];
    fitsfile	*digifits;
    int   	errflg=0, status=0, i, first, last, exp_stat = -2;
    int		timeref;
    float	val, ref_val;

    char    	key[N][9] = {"DET1ASCL", "DET1BSCL", "DET1AXOF", "DET1BXOF", "DET1AUCT", "DET1BUCT", "DET1ABWK", "DET1BBWK", "DET1AEWK", "DET1BEWK", "DET1ABSL", "DET1BBSL", "DET1ALCT", "DET1BLCT", "DET1ALTT", "DET1BLTT", "DET2ASCL", "DET2BSCL", "DET2AXOF", "DET2BXOF", "DET2AUCT", "DET2BUCT", "DET2ABWK", "DET2BBWK", "DET2AEWK", "DET2BEWK", "DET2ABSL", "DET2BBSL", "DET2ALCT", "DET2BLCT", "DET2ALTT", "DET2BLTT"};

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

    /*
     *  Determine which detector produced our data file.
     */
    FITS_read_key(infits, TSTRING, "DETECTOR", detector, 0, &status);
    if (detector[0] == '1')
	first = 0, last = N/2;
    else
	first = N/2, last = N;

    /*
     *  Read the digitizer reference file name and open it.
     */
    FITS_read_key(infits, TSTRING, "DIGI_CAL", file_name, 0, &status);
    FITS_open_file(&digifits, cf_cal_file(file_name), READONLY, &status);

    /*
     *  Read keywords.  Compare with reference values.
     */
    for (i = first; i < last; i++) {
        FITS_read_key(infits, TFLOAT, key[i], &val, NULL, &status);
        FITS_read_key(digifits, TFLOAT, key[i], &ref_val, NULL, &status);

        if (fabs(val - ref_val) > 2.0) {
            cf_if_warning("Digitizer keyword %s is out of bounds.", key[i]);
	    sprintf(comment, "Data suspect: keyword %s out of bounds", key[i]);
	    FITS_update_key(infits, TINT, "EXP_STAT", &exp_stat, comment,
		&status);    
            FITS_write_comment(infits, " ", &status);
	    sprintf(comment, "Data are suspect.  Keyword %s is out of bounds.",
		key[i]);
            FITS_write_comment(infits, comment, &status);
            fits_get_system_time(datestr, &timeref, &status);
            sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
            FITS_write_comment(infits, comment, &status);
            FITS_write_comment(infits, " ", &status);
	}
    }

    FITS_close_file(digifits, &status);
    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
