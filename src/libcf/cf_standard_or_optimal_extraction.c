/*****************************************************************************
 *            Johns Hopkins University
 *            Center For Astrophysical Sciences
 *            FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_standard_or_optimal_extraction (fitsfile *header,
 *			int *optimal)
 *
 * Description: Sets optimal to TRUE if optimal extraction is desired.
 *		FALSE if RUN_OPTI = 'NO' in parm*.fit file.
 *		FALSE if SRC_TYPE[0] = 'E' in file header.
 *
 * Arguments:   fitsfile  *header       Pointer to IDF FITS file header
 *              int       *optimal      TRUE if optimal extraction desired
 *
 * Calls:
 *
 * Returns:     0 on success
 *
 * History:     03/02/03  1.1   wvd	Initial coding.
 *		03/12/03  1.2   wvd	Write value of optimal to trailer file.
 *		03/19/03  1.3   wvd	Read SRC_TYPE, not SP_TYPE.
 *		03/19/03  1.4   wvd	Add call to cf_timestamp at start.	
 *					Change printf to cf_verbose
 *		09/30/03  1.5   wvd	Check for lower-case value of
 *					RUN_OPTI.
 *		03/15/05  1.6   wvd	No optimal extraction for HIST targets
 *					unless try_optimal = TRUE;
 *		03/18/05  1.7   wvd	v1.6 was a bad idea.  Returning to v1.5
 *
 ****************************************************************************/

#include <stdio.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_standard_or_optimal_extraction";
static char CF_VER_NUM[] = "1.7";

int
cf_standard_or_optimal_extraction(fitsfile *header, int *optimal)
{
    char  parm_file[FLEN_VALUE], instmode[FLEN_VALUE];
    char  run_opti[FLEN_VALUE], src_type[FLEN_VALUE];
    int   status=0;
    fitsfile *parmfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin processing");

    /* Read SRC_TYPE from IDF header. */
    FITS_read_key(header, TSTRING, "SRC_TYPE", src_type, NULL, &status);
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);

    /* Read RUN_OPTI from PARM_CAL. */
    FITS_read_key(header, TSTRING, "PARM_CAL", parm_file, NULL, &status);
    FITS_open_file(&parmfits, cf_parm_file(parm_file), READONLY, &status);
    FITS_read_key(parmfits, TSTRING, "RUN_OPTI", run_opti, NULL, &status);
    FITS_close_file(parmfits, &status);

    /* Determine value of optimal. */
    *optimal = TRUE;
    if (*run_opti != 'Y' && *run_opti != 'y') {
	cf_verbose (1, "RUN_OPTI = %s.  No optimal extraction.", run_opti);
	*optimal = FALSE;
    }
    if (*src_type != 'P') {
	cf_verbose (1, "SRC_TYPE = %s.  No optimal extraction.", src_type);
	*optimal = FALSE;
    }
    if (*optimal) cf_verbose (1, "Attempting optimal extraction.");

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
