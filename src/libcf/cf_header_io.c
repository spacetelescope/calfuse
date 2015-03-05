/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_copy_to_memory(*infile, **memp);
 *
 * Description: Returns an in-memory FITS file pointer
 *
 *
 * Arguments:   char     *infile        Input data file name
 *	        fitsfile **memp         In memory FITS file pointer
 *
 * Returns:     None
 *
 * History:     08/15/2002  1.0   jch   Initial coding
 *              12/18/03          bjg   Change calfusettag.h to calfuse.h
 *
 ****************************************************************************/


#include <string.h>
#include <stdio.h>
#include "calfitsio.h"
#include "calfuse.h"

void
cf_read_header(char *infile, fitsfile **memp)
{
    fitsfile *ip;       /* input file pointer */
    int	status=0;       /* status used in CFITSIO */

    /* open the input file */
    FITS_open_file(&ip, infile, READONLY, &status);

    /* create the temp file in memory */
    FITS_create_file(memp, "mem://", &status);

    /* copy the input primary HDU to the temp file */
    FITS_copy_hdu(ip, *memp, 0, &status);

    /* close the input file */
    FITS_close_file(ip, &status);
}


/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_copy_to_file(*memp, *outfile);
 *
 * Description: Copy the header in the in-memory FITS file to an output FITS
 *              file's primary HDU.
 *
 *
 * Arguments:   fitsfile *memp          In memory FITS file pointer
 *              char     *outfile       Output data file name
 *
 * Returns:     None
 *
 * History:     08/15/2002  jch   1.0   Initial coding
 *
 ****************************************************************************/

void
cf_write_header(fitsfile *memp, char *outfile)
{
    fitsfile *op;       /* output file pointer */
    int	status=0;       /* status used in CFITSIO */

    int	ncards_memp, ncards_op, morekeys=0;
    char card[FLEN_CARD];

    /* Open the output FITS file. */
    FITS_open_file(&op, outfile, READWRITE, &status);

    fits_get_hdrspace(memp, &ncards_memp, &morekeys, &status);
    fits_get_hdrspace(op,   &ncards_op,   &morekeys, &status);
	
    /* overwrite output file's primary HDU card by card */
    if (ncards_op >= ncards_memp) {
	int i;
	for (i=1; i<=ncards_memp; i++) {
	    FITS_read_record(memp, i, card, &status);
	    FITS_modify_record(op, i, card, &status);
	}
	/* delete from the bottom */
	for (i=ncards_op; i>=ncards_memp+1; i--) {
	    FITS_delete_record(op, i, &status);
	}
    }
    else {
	int i;
	for (i=1; i<=ncards_op; i++) {
	    FITS_read_record(memp, i, card, &status);
	    FITS_modify_record(op, i, card, &status);
	}
	for (i=ncards_op+1; i<=ncards_memp; i++) {
	    FITS_read_record(memp, i, card, &status);
	    FITS_write_record(op, card, &status);
	}
    }

    /* close output file. */
    FITS_close_file(op, &status);
}
