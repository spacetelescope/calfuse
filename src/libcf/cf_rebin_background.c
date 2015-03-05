/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences 
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:   cf_rebin_background(infits, aperture, nout, wave_out,
 *		  binx, biny, bnx, bny, bimage, barray) 
 *
 * Arguments:   *fitsfile  infits     :  Input IDF file
 *              int       aperture    :  Aperture ID for the analysis
 *              long      nout        :  Length of wave_out array
 *              float     wave_out    :  Wavelength vector
 *              int       binx, biny  :  Binning factors of the background array
 *              int       bnx, bny    :  X and Y dimension of the background array
 *              float     bimage      :  2-D background image (original) 
 *              float     barray      :  (OUTPUT) Binned background image.
 *
 * Description:
 *
 *	The input background image is binned in X.
 *	We need an image binned in wavelength.   
 *	To get it, we use the WAVE_CAL file (shifted to account
 *	for the FPA position) to derive a wavelength array for the
 *	background image.  For each row in the background image (bimage),
 *	we build a second image (barray) scaled to the size of the output 
 *	wavelength bins.
 *
 *  History     02/21/03   v1.1   rdr   Started work
 *              04/01/03   v1.2   rdr   Made aeff arrays double precision in 
 *                                      the cf_read_fluxcal routine to be 
 *                                      compatible with data in calfiles
 *		05/07/03   v1.4   wvd	Don't divide by EXPTIME.
 *					Read WPC from file header.
 *					Use cf_verbose.
 *              06/02/03   v1.5   rdr   Correct bug in filling output arrays
 *              06/11/03   v1.7   wvd	Pass datatype to cf_read_col.
 *					Change calfusettag.h to calfuse.h
 *              08/22/03   v1.8   bjg	Move get_extraction_limits to external
 *					subroutine.  Free orphan arrays.
 *					Change coltype from char to int in
 *					cf_read_col.
 *              08/28/03   v1.9   wvd	Debug
 *              09/29/03   v1.10  wvd	Shift background model to match
 *					position of FPA.
 *              03/16/04   v1.11  wvd	Correct error in calculation of dw_ave
 *					and dweff.
 *              04/05/04   v1.12  wvd	Modify i/o.
 *              04/26/04   v1.13  wvd	Do not flux-calibrate 2-D background
 *					model.  Do not extract 1-D background
 *					spectrum.
 *              05/13/04   v1.14  wvd	Correct value of CF_PRGM_ID.
 *              07/21/04   v1.16  wvd	Delete unused variables.
 *
 *****************************************************************************/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"


int cf_rebin_background(fitsfile *infits, int aperture, long nout,
	float *wave_out, int binx, int biny, int bnx, int bny,
	float *bimage, float **barray) {

    char CF_PRGM_ID[] = "cf_rebin_background";
    char CF_VER_NUM[] = "1.16";

    fitsfile  *wavefits ;
    int       status=0;
    int       dndx, dpix, nyout;
    long      npix, npix_out, nwave, ndx, ndx1, ndx2;
    long      i, j, k;
    float     dwb, dw_out, mf, ctot, dw_ave;
    float     dpix_fpa;
    float     *wave, *wavet, *outarr;
    char      wavefile[FLEN_VALUE], det[FLEN_VALUE] ;

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");

    /* Read wavelength spacing of the output array */
    FITS_read_key(infits, TFLOAT, "WPC", &dw_out, NULL, &status) ;
    cf_verbose(3, "Bin size of the output wavelength array = %5.3f A",
	dw_out) ;
    FITS_read_key(infits, TSTRING, "DETECTOR", det, NULL, &status) ;
    cf_verbose(3,"Detector = %s, aperture = %d", det, aperture) ;

    /* Read wavelength-calibration file */
    FITS_read_key(infits, TSTRING, "WAVE_CAL", wavefile, NULL, &status) ;
    cf_verbose(3, "Wavelength calibration file = %s",wavefile) ;
    FITS_open_file(&wavefits, cf_cal_file(wavefile), READONLY, &status) ;
    cf_verbose(3, "HDU to read = %d",aperture+1) ;
    FITS_movabs_hdu(wavefits, aperture +1, NULL, &status) ;
    nwave=cf_read_col(wavefits,TFLOAT,"WAVELENGTH",(void **) &wave);
    cf_verbose(3, "Points read = %d, min & max: %4.2f, %4.2f",
	nwave, wave[0], wave[nwave-1]);

    /* Read pixel shifts due to FPA positions */
    if (aperture < 5) 
	FITS_read_key(infits, TFLOAT, "FPADXLIF", &dpix_fpa, NULL, &status);
    else 
	FITS_read_key(infits, TFLOAT, "FPADXSIC", &dpix_fpa, NULL, &status);
    dpix = cf_nint(dpix_fpa);

    /* To correct for the position of the FPA, each photon has been shifted
	in X by dpix pixels.  Rather than shifting the background to match,
	we modify the wavelength array that maps the background onto the
	output spectrum. */
    cf_verbose(2, "Shifting aperture %d background by %d X pixels "
	"to match the FPA shift.", aperture, cf_nint(dpix_fpa));

    /*  Generate a new wavelength array binned to match the background array */
    cf_verbose(3, "X and Y binning factors for background array = %d, %d",
	binx, biny) ;

        wavet = (float *) cf_calloc(bnx, sizeof(float));
	for (i = 0; i < bnx; i++) {
          ndx=i * binx + (binx/2) + dpix;
	  if (ndx > nwave-1) ndx = nwave-1 ;
	  if (ndx < 0) ndx = 0;
	  wavet[i] = wave[ndx];
	 }

      /* Error check: average wavelength spacing */
      dw_ave = fabs(wavet[bnx-1]-wavet[0])/(bnx-1);
      cf_verbose(3,"Average wavelength spacing of background array = %f", dw_ave) ;

      if ((aperture < 4 && strncmp(det,"1",1) == 0) ||
	  (aperture > 4 && strncmp(det,"2",1) == 0)) {
           cf_verbose(3, "Maximum tabulated wavelength = %.2f", wavet[bnx-1]);
           cf_verbose(3, "Maximum of wave_out = %.2f", wave_out[nout-1] );
      }
      else {
           cf_verbose(3, "Maximum tabulated wavelength = %.2f", wave_out[nout-1] );
           cf_verbose(3, "Maximum of wave_out = %.2f", wave_out[nout-1] );
      }

  /* Define some sizes and allocate space for the output arrays */
  nyout = bny * biny ;		/* Y dimension of output 2-D bkgd array */
  npix_out = nout * nyout ;	/* Size of output 2-D bkgd array */
  npix = bnx * bny ;		/* Size of input 2-D bkgd array */
  outarr = (float *) cf_calloc(npix_out, sizeof(float) ); /* 2-D output array */

    /* Error check: total counts in the background image */
    ctot=0.;
    for (i=0; i<npix; i++) ctot+=bimage[i] ;
    cf_verbose(3, "Number of counts in the input bkgd image = %f",ctot) ;

  /* For each wavelength in the output array, find the nearest point
     in the input background array.  Scale its value by the relative
     sizes of the background and output wavelength bins.

     Remember: wavelengths for the Lif channels and the SiC channels run
               in opposite directions.  For side one, wavelengths increase 
               with index for Lif but decrease for SiC. The opposite is true
	       for side two. */
  ndx=0 ;
  dndx = 1 ; 
  if ((aperture < 4 && strncmp(det,"2",1) == 0)  ||
      (aperture > 4 && strncmp(det,"1",1) == 0 ) ) {
    ndx = bnx - 1 ;
    dndx = -1 ;
  }

  cf_verbose(3,"Filling output: starting ndx=%d, increment =%d ", ndx, dndx) ;

  dwb = dw_ave;
  for (i=0; i<nout; i++) {
    /* Locate the element of the bkgd array corresponding to this wavelength */
        while (wavet[ndx] < wave_out[i] && ndx < bnx && ndx >= 0) ndx += dndx ;

	/* compute the relative
           wavelength coverage of the background and output arrays */
        if (ndx < bnx-2) dwb = fabs((double) (wavet[ndx+1]-wavet[ndx])) ;
	   mf = dw_out/dwb ;
	/* correct the scale factor for binning in y */
           mf /= biny ;
	/* populate the column of the array */
	for (j=0; j< nyout ; j++) {
          k = j/biny ;
          ndx1 = j*nout + i ;
	  if (ndx1 > npix_out-1) ndx1 = npix_out-1 ;
	  ndx2 = k*bnx + ndx ;
	    if (ndx2 > npix-1) ndx2=npix-1 ;
	  outarr[ndx1] = mf * bimage[ndx2] ; 
	}
  }


    /* Error checking: total counts in the background image */
    ctot=0.;
    for (i=0; i<npix_out; i++) ctot+=outarr[i] ;
    cf_verbose(3, "Number of counts in the output bkgd image = %f", ctot) ;


 /* Return output arrays */
    *barray = outarr ;
     
 /* Close wavelength calibration file */
      FITS_close_file(wavefits, &status) ;

   free(wave);
   free(wavet);

   cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
   return status;
}
