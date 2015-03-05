/*******************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *******************************************************************************
 *
 * Synopsis:	cf_fes_write(fitsfile *outfits, int hdu, float *image);
 *
 * Description: writes an FES image into HDU number hdu in the FITS file
 *              pointed to by *fptr.  Currently cf_fes_write extracts
 *              BITPIX, NAXIS1, and NAXIS2 from the specified HDU and 
 *              uses those values to write the image.  
 *
 * Arguments:	fitsfile   *outfits	output FITS file pointer
 *		int	   hdu 		HDU to be written to
 *		float      *image  	Pointer to array holding FES image
 *
 * Returns:	none
 *
 * History:	07/08/98    gak    calfes_design.070898 design documented
 *              07/20/98    mlr    started work
 *              04/19/99    mlr    finished modifications to utilize
 *                                 libcf and FITSIO.h(error handling).
 * ToDO:  I may at somepoint modify this subroutine to pass in the image
 *        type and size in as parameters... I already have the information
 *        I shouldn't read it in... or perhaps I should just verify that 
 *        they match???
 * ToDo:  I should check bitpix and write out the requested size.  At the
 *        moment I only write float.
 * 
 ******************************************************************************/

#include "calfuse.h"

#define CF_PRGM_ID      "cf_fes_write"

int cf_fes_write(fitsfile *outfits, int hdu, float *image)
{
	char	buffer[FLEN_CARD];
	int 	status, hdutype; 
	int	naxis, naxis1, naxis2;
        int     fpixel, npixels, bitpix;
        short   *short_image;

        hdutype = -1;       
        status = 0;

	FITS_movabs_hdu(outfits, hdu, &hdutype, &status);
        
        FITS_read_key(outfits, TINT, "NAXIS1", &naxis1, buffer, &status);
        FITS_read_key(outfits, TINT, "NAXIS2", &naxis2, buffer, &status);
        FITS_read_key(outfits, TINT, "BITPIX", &bitpix, buffer, &status);

        fpixel = 1;
        npixels = naxis1 * naxis2;

        FITS_write_img(outfits, TFLOAT, fpixel, npixels, image, &status);

        return(status);
}
