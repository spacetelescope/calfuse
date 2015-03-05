/*******************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *******************************************************************************
 *
 * Synopsis: cf_apply_fes_mask(fitsfile *fesfits, float **image, 
 *                                       float *mask, int axis1, int axis2)
 *
 * Description: Apply FES mask in *mask to the raw FES image in **image
 *     To apply the mask I simply loop through the 
 *     mask image pixel by pixel and 
 *    If the mask pixel value != 0.0
 *    then fes image pixel value = mask pixel value.
 *
 * Arguments:	fitsfile *fesfits       fits file for image, so that we can add
 *                                      history/comment lines to the HDU
 *              float   **image 	FES image
 *              float   *mask	 	FES mask
 *              int     axis1, axis2    size of each axis of the images.
 *                                      This presumes that both images are 
 *                                      the same size when they get here.
 *
 * Returns:	0 upon successful completion
 *
 * History:	07/08/98	gak     calfes_design.070898 design documented
 *              07/22/98        mlr     started work
 *              04/06/99        mlr     Broke subroutine out from cf_fes_mask.c
 *              08/10/99        mlr     added fesfits to the argument list
 *
 *
 ******************************************************************************/

#include "calfuse.h"
#include "cf_calfes.h"

#define CF_PRGM_ID      "cf_fes_apply_mask"
#define CF_VER_NUM      "1.4"

int cf_fes_apply_mask(fitsfile *fesfits, float **image,  
                               float *mask, int naxis1, int naxis2)
{
        int i, j;  /*loop counters */

        cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Applying FES pixel mask");

        for (j=0; j < naxis2; j++)
             for (i=0; i < naxis1; i++)
               if (mask[j*naxis1+i] != FES_GOOD_PIX )
                     (*image)[j*naxis1+i] = mask[j*naxis1+i];


       cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done applying FES  mask");

       return(0);
}
