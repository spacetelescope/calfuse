/*******************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *******************************************************************************
 *
 * Synopsis:	cf_cp_hdr(fitsfile *infits, int inhdu, fitsfile *outfits,
 *			int outhdu)
 *
 * Description: Copies the header, with no changes, from HDU inhdu of the
 *	 	input FITS file infits to HDU outhdu of the output FITS file
 *		outfits.
 *
 * Arguments:	fitsfile	*infits		Input FITS file pointer
 *		int		inhdu		HDU to be copied from
 *		fitsfile	*outfits	Output FITS file pointer
 *		int		outhdu		HDU to be copied to
 *
 * Returns:	none
 *
 * History:	04/02/98	gak	Begin work
 *		08/27/98	gak	Initialized status, hdutype;
 *					return status
 *
 ******************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

#define CF_PRGM_ID      "cf_cp_hdr"
#define CF_VER_NUM      "1.1"

int cf_cp_hdr(fitsfile *infits, int inhdu, fitsfile *outfits, int outhdu)
{
	char		buffer[FLEN_CARD];
	int 		status, hdutype; 
	int		nkeys, keynum;
	int		i;

	status = 0;
	hdutype = 0;

	if ( fits_movabs_hdu(infits, inhdu, &hdutype, &status) ) {
		cf_if_error(status);
		exit(status);
	}
	if ( fits_movabs_hdu(outfits, outhdu, &hdutype, &status) ) {
		cf_if_error(status);
		exit(status);
	}

	/* Get number of keywords in the input header, loop & copy. */
	if ( fits_get_hdrpos(infits, &nkeys, &keynum, &status) ) {
		cf_if_error(status);
		exit(status);
	}
	for ( i = 1; i <= nkeys; i++) {
		if ( fits_read_record(infits, i, buffer, &status) ) {
			cf_if_error(status);
			exit(status);
		}
		if ( fits_write_record(outfits, buffer, &status) ) {
			cf_if_error(status);
			exit(status);
		}
	}

	return status;
}
