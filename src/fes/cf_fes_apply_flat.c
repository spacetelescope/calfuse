/*******************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *******************************************************************************
 *
 * Synopsis: cf_apply_fes_flat(fitsfile *fesfits, float **image, 
 *                                       float *flat, int axis1, int axis2)
 *
 * Description: Apply FES flatfield in *flat to the raw FES image in **image
 *     To apply the flatfield, I simply loop through the 
 *     flatfield image pixel by pixel and 
 *    If the flatfield pixel value != 0.0
 *    then fes image[] pixel value =  image[]pixel value/mask[]pixel value.
 *    else (If the flatfield pixel value == 0.0)
 *    then fes image[] pixel value =  FES_BAD_PIX.
 *
 * Arguments:	fitsfile *fesfits       fits file for image..  so that we can add
 *                                      history/comment lines to the HDU
 *              float   **image 	FES image
 *              float   *mask	 	FES mask
 *              int     axis1, axis2    size of each axis of the images.
 *                                      This presumes that both images are 
 *                                      the same size before they get here.
 *
 * Returns:	0 upon successful completion
 *
 * History:	07/08/98	gak     calfes_design.070898 design documented
 *              07/22/98        mlr     started work
 *              04/06/99        mlr     Broke this actual subroutine out
 *                                      from cf_fes_mask.c
 *              04/20/99        mlr     instead of calculating the image row
 *                                      size each time -- pass it in
 *              08/10/99        mlr     added fesfits to argument list
 *
 ******************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"
#include "cf_calfes.h"

#define CF_PRGM_ID      "cf_fes_apply_flat"
#define CF_VER_NUM      "1.4"

int cf_fes_apply_flat(fitsfile *fesfits, float **image,  
                            float *flat, int naxis1, int naxis2)
{
        int i, j;  /*loop counters */
        int pixel;

        cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Applying FES flatfield");
        for (j=0; j < naxis2; j++)
             for (i=0; i < naxis1; i++)
             {
               pixel = j* naxis1+ i;            
               if ((*image)[pixel] != FES_BAD_PIX)
                   if (flat[pixel] != 0.0 )
                        (*image)[pixel] /= flat[pixel];
                   else 
                      (*image)[pixel] = FES_BAD_PIX;
             }
       cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done applying FES flat");

       return(0);
}
