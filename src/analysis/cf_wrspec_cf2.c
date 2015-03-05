/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Synopsis:    cf_wrspec_cf2(fitsfile *outfits, int npts, float *wave,
 *				float *spec, float *errs, short *qual,
 *                              float *counts, float *cntserr, int ncol,
 *				int areaflag);
 *
 * Description: Write a file containing a FUSE 1D spectrum, in a format
 *		compatible with Calfuse 2.4 and earlier.
 *
 * Arguments:   char    *fname		Output file name
 *		float	*wave		wavelengths
 *		float	*spec		spectrum
 *		float	*errs		associated 1-sigma error bars
 *		short	*qual		associated quality flags
 *		float	*counts		total counts in column
 *		float	*cntserr	1-sigma error bar on total counts
 *		int	ncol		specifies 4 or 6 column format
 *		int	areaflag	0 if spectrum column is FLUX
 *					1 if spectrum column is AREA
 *
 * Returns:     none
 *
 * History:     05/04/98        gak     Begin work.
 *              04/20/99        emm     Added FITS_ error checking routines,
 *                                      converted qual flags to BYTE.
 *              06/07/99 1.2    peb     Added reporting of version number.
 *              10/22/99 1.3    emm     Added total counts column.
 *              10/22/99 1.4    emm     Added total counts errors column.
 *              10/22/99 1.5    emm     totcnts and errors are just duplicates
 *                                      of flux and error for now
 *              10/29/99 1.6    emm     counts and cntserr are now read in
 *                                      from calling program.
 *              12/21/99 1.7    emm     Quality flags now short (1I) instead
 *                                      of byte char (1B).
 *              11/30/01 1.8    wvd     Don't write timestamp to log file.
 *		03/08/04 1.9    jwk     allow tfields=4 or 6, and option to
 *					change "FLUX" col to "AREA"; for use
 *					with cf_arith
 *              04/05/04        bjg     Include string.h
 *                                      Write timestamp
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_wrspec_cf2";
static char CF_VER_NUM[] = "1.9";

void cf_wrspec_cf2(fitsfile *outfits, int npts, float *wave, float *spec,
		float *errs, short *qual, float *counts, float *cntserr,
		int ncol, int areaflag)
{
	char	*ttype[] = {"WAVE", "FLUX", "ERROR", "QUALITY", "COUNTS","CNTSERR"};
	char	*tform[] = {"1E", "1E", "1E", "1I", "1E", "1E"};
	char	*tunit[] = {"ANGSTROMS", "ERG/CM2/S/A", "ERG/CM2/S/A", " ","COUNTS","COUNTS"};
	char	extname[] = "FUSE 1D Spectrum";
	int	tfields = 6,
		status = 0;
	long	felem = 1,
		frow = 1;	

/* Write time stamp to log file. */
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "");


/* Create the extension for the data in the output file.
   Columns will be
	WAVE	1E
	FLUX	1E
	ERROR	1E
	QUALITY	1B
	TOTAL COUNTS 1E      Same as FLUX column for now
	COUNTS ERRORS 1E     Same as ERROR column for now
*/

    if (ncol == 4)
	{
	tfields = 4;
	if (areaflag)
	   {
	   strcpy(ttype[1], "AREA");
	   strcpy(tunit[1], "CM^2");
	   strcpy(tunit[2], "CM^2");
	   }
	}

	fits_create_tbl (outfits, BINARY_TBL, (long)npts, tfields, ttype,
	     tform, tunit, extname, &status);

/* Write out the data. */
	FITS_write_col(outfits, TFLOAT, 1, frow, felem, (long)npts, 
		       wave, &status);

	FITS_write_col(outfits, TFLOAT, 2, frow, felem, (long)npts,
		       spec, &status);

	FITS_write_col(outfits, TFLOAT, 3, frow, felem, (long)npts,
		       errs, &status );

	FITS_write_col(outfits,  TSHORT, 4, frow, felem, (long)npts,
		       qual, &status);

	if (tfields == 6)
	  {
	  FITS_write_col(outfits, TFLOAT, 5, frow, felem, (long)npts,
		       counts, &status);

	  FITS_write_col(outfits, TFLOAT, 6, frow, felem, (long)npts,
		       cntserr, &status);
	  }

}
