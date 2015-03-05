/*******************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *******************************************************************************
 *
 * Synopsis:	int cf_fes_cal(int cal_type, fitsfile *fesfits, float **image, 
 *                                    int inaxis1, int inaxis2)
 *
 * Description: Get bias calibration filename for the raw FES image in
 *                  the current hdu pointed at by *infits.
 *              Call cf_fes_get_cal_image() to open the cal file and get
 *                                          the calibration image.
 *              Call cf_fes_apply_flat();
 *              return the modified image to the calling routine.
 *              4/10/99  MLR removed int inhdu as an argument and removed
 *                       the call to fits_movabs_hdu.  This requires any
 *                       user of this subroutine to do these steps before
 *                       calling cf_fes_flat.
 *              4/10/99  MLR removed (fitsfile *outfits, int outhdu) from 
 *                       the argument list.  We will leave writing the image
 *                       out to the calling routine.
 *   
 *
 * Arguments:	fitsfile   *infits	Input FITS file pointer
 *              float **image           already read in FES image
 *              int inaxis1             image size  of the input image
 *              int inaxis2
 *
 * Returns:	none
 *
 * History:	07/08/98	gak     calfes_design.070898 design documented
 *              07/22/98        mlr     started work
 *              04/10/99        mlt     pulled the actual work out of this
 *                                      routine and into cf_fes_apply_flat
 *              04/12/99        mlr     modified the input arguments(see above)
 *                                      use the now generic cf_fes_read to
 *                                      get the calibration image.
 *              04/14/99        mlr     pulled the reading out and put it 
 *                                      into cf_fes_get_cal_image.c
 *              04/16/99        mlr     finished modifications to utilize
 *                                      libcf and FITSIO.h(error handling).
 *
 * ToDo:  4/12/99 I have added the naxis1 and naxis2 parameters to the 
 *                argument list.  I need to use them to compare binning
 *                factors in the two images.
 *                I should combine cf_fes_mask, cf_fes_bias, cf_fes_flat
 *                into one generic cf_fes_do_cal("cf_fes_mask", *fesfits,...etc)
 *
 ******************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"
#include "cf_calfes.h"

#define CF_VER_NUM  "1.4"

int cf_fes_cal(int cal_type, fitsfile *fesfits, float **image, 
                              int inaxis1, int inaxis2)
{
	char	buffer[FLEN_CARD];
        char    calfilename[30];
        fitsfile *calfile;
	int 	status; 
        float   *calimage;
	char *prog_id, *msg, *keyword, *msg2;

	int (*func)(fitsfile *, float **, float *, int, int);

        status = 0;

        switch (cal_type) {
           case FES_FLAT:
		prog_id = "cf_fes_flat";
		msg = "Begin FES flatfield";
		msg2 = "Done  FES flatfield";
		keyword = "FFLT1FCL";
		func =  cf_fes_apply_flat;
		break;
           case FES_MASK:
		prog_id = "cf_fes_mask";
		msg = "Begin FES mask";
		msg2 = "Done  FES mask";
		keyword = "MASK_FCL";
		func =  cf_fes_apply_mask;
		break;
           case FES_BIAS:
		prog_id = "cf_fes_bias";
		msg = "Begin FES bias";
		msg2 = "Done  FES bias";
		keyword = "BIAS1FCL";
		func =  cf_fes_apply_bias;
		break;
	}


	cf_timestamp(prog_id, CF_VER_NUM, msg);

        FITS_read_key(fesfits, TSTRING, keyword, calfilename, buffer, &status);

        if (cf_fes_get_cal_image(calfilename, &calimage, inaxis1, inaxis2) == -1)
           return (-1);
       
        (*func)(fesfits, image, calimage, inaxis1, inaxis2);

        cf_fes_proc_update(fesfits, prog_id, "COMPLETE");

        free(calimage);

        cf_timestamp(prog_id, CF_VER_NUM, msg2);

        return(status);
}
