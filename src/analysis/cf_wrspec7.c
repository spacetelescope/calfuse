/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Synopsis:    cf_wrspec7(fitsfile *outfits, long npts, float *wave, float *flux,
 *		           float *error, long *counts, float *weights,
 *                         float *bkgd,short *pothole)
 *
 * Description: Write a file containing a FUSE 2D spectrum.
 *
 * Arguments:   fitsfile *fname		Output FITS file structure
 *
 * Returns:     none
 *
 * History:     12/15/03        bjg     Begin work from cf_wrspec4.
 *              01/12/04        bjg     Call to fits_create_tbl with nrow=1
 *                                      instead of npts
 *              03/22/04        bjg     Change POTHOLE to QUALITY
 *              03/25/04        bjg     Moved time stamp writing to beginning 
 *                                      of routine
 *              04/07/05  v1.1  tc      Create a multi-rows table rather than
 *                                      a unique vector in a cell
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_wrspec7";
static char CF_VER_NUM[] = "1.1";

void cf_wrspec7(fitsfile *outfits, long npts, float *wave, float *flux,
		float *error, long *counts, float *weights,float *bkgd,short *pothole)
{
	char	*ttype[] = {"WAVE", "FLUX", "ERROR", "COUNTS", "WEIGHTS", "BKGD", "QUALITY" };
	char	*tform[] = {"1E", "1E", "1E", "1J", "1E", "1E", "1I"};
	char	*tunit[] = {"ANGSTROMS", "ERG/CM2/S/A", "ERG/CM2/S/A", "COUNTS", "COUNTS", "COUNTS", "UNITLESS"};
	char	extname[] = "FUSE 2D Spectrum";
	int	tfields = 7, status = 0;

	/* Write time stamp to log file. */
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "");
		
	fits_create_tbl(outfits, BINARY_TBL, 1, tfields, ttype,	tform, tunit, extname, &status);
     
	/* Write out the data. */
	FITS_write_col(outfits, TFLOAT, 1, 1L, 1L, npts, wave, &status);
	FITS_write_col(outfits, TFLOAT, 2, 1L, 1L, npts, flux, &status);
	FITS_write_col(outfits, TFLOAT, 3, 1L, 1L, npts, error, &status );
	FITS_write_col(outfits, TLONG, 4, 1L, 1L, npts, counts, &status);
	FITS_write_col(outfits, TFLOAT, 5, 1L, 1L, npts, weights, &status);
	FITS_write_col(outfits, TFLOAT, 6, 1L, 1L, npts, bkgd, &status);
	FITS_write_col(outfits, TSHORT, 7, 1L, 1L, npts, pothole, &status);
	
}
