/*******************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *******************************************************************************
 *
 * Synopsis: cf_apply_fes_bias(fitsfile *fesfits, float **image, 
 *                                       float *bias, int axis1, int axis2)
 *
 * Description: Apply FES bias in *bias to the raw FES image in **image
 *     To apply the bias, I simply loop through the 
 *     bias image pixel by pixel and subtract that 
 *     value from the given FES image.
 *
 * Arguments:	fitsfile *fesfits       fits file for image so that we can add
 *                                      history/comment lines to the HDU
 *              float   **image 	FES image
 *              float   *bias	 	FES bias image
 *              int     axis1, axis2    size of each axis of the images.
 *                                      This presumes that both images are 
 *                                      the same size when they get here.
 *
 * Returns:	0 upon successful completion
 *
 * History:	07/08/98	gak     calfes_design.070898 design documented
 *              07/22/98        mlr     started work
 *              04/06/99        mlr     Broke this actual subroutine out
 *                                      from cf_fes_bias.c
 *              08/10/99        mlr     added fesfits to the argument list
 *              04/25/02        wvd	initialize status to zero
 *
 ******************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"
#include "cf_calfes.h"

#define CF_PRGM_ID      "cf_fes_apply_bias"
#define CF_VER_NUM      "1.4"

int cf_fes_apply_bias(fitsfile *fesfits, float **image,  
                               float *bias, int naxis1, int naxis2)
{
        int i, j;  /*loop counters */
        int status=0;
        float os_bias, os_data;
        float factor;
        char buffer[80];

        cf_timestamp(CF_PRGM_ID,CF_VER_NUM, "Applying FES bias");
        
        os_bias=0.;
        os_data=0.;

	if (naxis1 == 520)  /* 1x1 binning */
   	{
           for (j=0; j<520; j++)
	      for (i=512; i<520; i++)
              {
	         os_bias += bias[j*naxis1+i];
	         os_data += (*image)[j*naxis1+1];
              }
        }
        else if (naxis1 == 260)  /* 2x2 binning */
        {
           for (j=0; j<260; j++)
               for (i=256; i<260; i++)
               {
	         os_bias += bias[j*naxis1+i];
	         os_data += (*image)[j*naxis1+1];
               }
        }
        else   /*bail out -- unsupported image size or binning */
            exit(-1);

      factor = os_data / os_bias;
      sprintf(buffer, "FES BIAS FACTOR = %lf",factor);
      cf_timestamp(CF_PRGM_ID, CF_VER_NUM, buffer);

/*      FITS_write_history(fesfits, buffer, &status); */

      for (j=0; j<naxis2; j++)
         for (i=0; i<naxis1; i++)
         {
/*     printf("pixel[%d, %d] = %lf\n",i,j, (*image)[j*naxis1+i]); */
                    if ( (*image)[j*naxis1+i] != FES_BAD_PIX )
                          (*image)[j*naxis1+i] -= bias[j*naxis1+i] * factor;
                    else
                        printf("bad pixel[%d, %d]\n",i,j);
          }

      cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished FES bias");

      return(status);

}
