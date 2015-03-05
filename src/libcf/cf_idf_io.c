/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    npts = cf_read_col(fitsfile, coltype, colname, colval) 
 *
 * Description: Procedure to read the values from a column in a FUSE
 *              Intermediate Data File (IDF).
 *		Because floating-point data may be stored as shorts,
 *		must convert data to type expected by calling routine.
 *
 * Variables:   fitsfile  *fileptr      Pointer to the input file
 *              int       coltype       Data type of the column
 *					  (TBYTE, TDOUBLE, TSHORT, or TFLOAT)
 *              char      *colname      Name of the column to re read
 *                                        (e.g. TIME or XFARF)
 *              void      **colval      Address of the array containing
 *                                      the data (must be cast to void)
 *
 * Return:      long      npts          Number of points in the array
 *
 * Example:
 *              fitsfile *infits ;
 *              long npts ;
 *              float *xfarf ;
 *       
 *              FITS_open_file(&infits, argv[1], READWRITE, &status); 
 *              npts = cf_read_col(infits, TFLOAT, "XFARF", (void **) &xfarf); 
 *
 * History:     08/08/02    rdr		Begin work
 *		06/11/03    wvd  v1.4  	Pass coltype from calling routine.
 *		08/25/03    wvd  v1.5  	Change coltype from string to int.
 *		08/28/03    wvd  v1.6  	Use CASEINSEN in calls to
 *					fits_get_colnum
 *              12/18/03    bjg  v1.7   Change calfusettag.h to calfuse.h
 *              04/28/05    wvd  v1.8   Allow arrays to be read as doubles.
 *              08/24/07    wvd  v1.9   Read nrows as type TLONG.
 *
 *************************************************************************/

#include <string.h>
#include "calfuse.h"

long
cf_read_col(fitsfile *infits, int coltype, char *colname, void **colval)
{
    int  colnum, tcode, intnull=0, anynull, status=0;
    long nrows, width, repeat, nentries;

    /* Get the column number of the column to be extracted. */
    FITS_get_colnum(infits, CASEINSEN, colname, &colnum, &status);

    /* Determine the properties of the column. */
    fits_get_coltype(infits, colnum, &tcode, &repeat, &width, &status);

    /* Read the number of rows in the table. */
    FITS_read_key(infits, TLONG, "NAXIS2", &nrows, NULL, &status);

    /* Calculate the number of entries in the table. */
    nentries = nrows * repeat;

    /* Set tcode and width to values expected by calling routine. */
    if (coltype == TFLOAT) {
        tcode = TFLOAT;        
        width = (long) sizeof (float);
    }
    if (coltype == TDOUBLE) {
        tcode = TDOUBLE;        
        width = (long) sizeof (double);
    }

    /* Allocate space for the file. */
    *colval = (void *) cf_malloc(nentries * width);
    
    /* Read the data from the file. */
    FITS_read_col(infits, tcode, colnum, 1, 1, nentries, &intnull, *colval,
		  &anynull, &status); 

    return nentries;
}

/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_write_col(fitsfile, colname, colval, npts) 
 *
 * Description: Procedure to write an array into a specified column of a
 *              FUSE Intermediate Data File (IDF).
 *		Because floating-point data may be stored as shorts,
 *		data type of array may not match that of output file header.
 *
 * Variables:   fitsfile  *fileptr      Pointer to the input file
 *              int       coltype       Data type of the column
 *					  (TBYTE, TDOUBLE, TSHORT, or TFLOAT)
 *              char      *colname      Name of the column to re read
 *                                        (e.g. TIME or XFARF)
 *              void      **colval      Address of the array containing
 *                                      the data (must be cast to void)
 *              long      npts          Number of points in the array
 *
 * Reture:	status
 *
 * Example:
 *              fitsfile *infits ;
 *              long npts=100 ;
 *              float *xfarf ;
 *       
 *              FITS_open_file(&infits, argv[1], READWRITE, &status);
 *              npts = cf_read_col(infits, TFLOAT, "XFARF", (void **) &xfarf);
 *
 * History:     08/08/02    rdr         Begin work
 *		06/11/03    wvd  v1.4  	Pass coltype from calling routine.
 *					Check for overflow in X arrays.
 *		08/25/03    wvd  v1.5  	Change coltype from string to int.
 *		08/28/03    wvd  v1.6  	Use CASEINSEN in calls to
 *					fits_get_colnum
 *
 *************************************************************************/

int
cf_write_col(fitsfile *infits, int coltype, char *colname, void *colval,
	long npts)
{
    int  colnum, tcode, status=0;
    long i, width, nevents;

    /* Get the column number of the column to be written. */
    FITS_get_colnum(infits, CASEINSEN, colname, &colnum, &status);

    /* Determine the properties of the column. */
    fits_get_coltype(infits, colnum, &tcode, &nevents, &width, &status);

    /* 
     * If the calling routine and the IDF differ on the data type of this
     * array, use the data type from the calling routine.
     */
    if (coltype == TFLOAT) {
        tcode = TFLOAT;        
        width = (long) sizeof (float);
    }

    /* Check for overflow in X arrays. */
    if (!strncmp(colname, "X", 1)) {
	float *x;
	x = (float *) colval;
	for (i = 0; i < nevents; i++) {
            if (x[i] < 0) x[i] = 0.;
            if (x[i] > NXMAX - 1) x[i] = NXMAX - 1;
	}
    }

    /* Write the data to the relevant column. */
    FITS_write_col(infits, tcode, colnum, 1, 1, npts, colval, &status);

    return status;
}
