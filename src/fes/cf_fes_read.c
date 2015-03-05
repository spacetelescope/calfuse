/*******************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *******************************************************************************
 *
 * Synopsis:	cf_fes_read(fitsfile *infits, float **image, 
 *                                            int *naxis1, int *naxis2);
 *
 * Description: reads an FES image from the current hdu in the FITS file
 *              pointed to by *infits and stores it in the array pointed
 *              to by *image.  
 *              *infits may contain an FES image where each pixel is 
 *              stored as unsigned shorts or as floats. cf_fes_read
 *              determines this from BITPIX in the current hdu and proceeds 
 *              accordingly.
 *              cf_fes_read also takes care of any binning by determining 
 *              the image size from NAXIS1 and NAXIS2 in the current hdu.
 *              4/12/99  MLR These two values are now returned to the 
 *                           calling routine for reference.
 *              4/10/99  MLR removed int inhdu as an argument and removed
 *                       the call to fits_movabs_hdu.  This requires any
 *                       user of this subroutine to do these steps before
 *                       calling cf_fes_read. 
 *
 * Arguments:	fitsfile   *infits	Input FITS file pointer
 *		float      **image  	Pointer to array holding FES image
 *              int        *naxis1      Image size parameters
 *              int        *naxis2
 * 
 * Returns:	0 upon successfule completion
 *
 * History:	07/08/98	gak     calfes_design.070898 design documented
 *              07/16/98        mlr     started work
 *              04/12/99        mlr     modified the input arguments(see above)
 *              04/16/99        mlr     finished modifications to utilize
 *                                      libcf and FITSIO.h(error handling).
 ******************************************************************************/

#include "calfuse.h"

#define CF_PRGM_ID      "cf_fes_read"
#define CF_VER_NUM      "1.4"

int cf_fes_read(fitsfile *infits, float**image, int *naxis1, int *naxis2)
{
	char	buffer[FLEN_CARD];
	int 	status; 
	int	bitpix, nullval, anynull, npixels;
        int     i, j, hdutype;  
        short   *short_image;

        status = 0;
        hdutype=-1;

	FITS_read_key(infits, TINT, "NAXIS1", naxis1, buffer, &status);
 
        FITS_read_key(infits, TINT, "NAXIS2", naxis2, buffer, &status);

        npixels = (*naxis1) * (*naxis2);

        cf_if_memory_error(*image = (float *) malloc(sizeof(float) * npixels)); 

        FITS_read_key(infits, TINT, "BITPIX", &bitpix, buffer, &status);

        nullval = 0;
        if (bitpix == 16 || bitpix == -32)

             FITS_read_img(infits, TFLOAT, 1, npixels, 
                                     &nullval, *image, &anynull, &status);

        else 
             cf_if_error("INVALID data type for FES image. ");

        return(status);
}
