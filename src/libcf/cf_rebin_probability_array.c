/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences  
 *              FUSE  
 *****************************************************************************
 *
 * Synopsis: cf_rebin_probability_array(infits, extended, aperture, nout,
 *		wave_out, pny, pycent, parray)
 *
 * Arguments:     *fitsfile  infits     :  Input IDF file
 *                int        extended   :  TRUE if target is extended
 *                int        aperture   :  Aperture ID for the analysis
 *                long       nout       :  number of output wavelength bins
 *                float      wave_out   :  output wavelength vector
 *                int        pny        :  y dimension of probability array
 *                float      pycent     :  y position of the peak of the 
 *                                           probability array
 *                float      parray     :  probability array
 *
 * History:       02/20/03   v1.1   rdr    Started work
 *                04/03/03   v1.2   rdr    Normalize so that integral along y
 *                                         at each sample point is 1
 *                06/02/03   v1.3   rdr    Correct error in filling output array
 *                09/30/03   v1.5   wvd    Use cf_read_col to read WAVE_CAL.
 *                03/22/04   v1.6   wvd    Change pycent from int to float.
 *                03/23/04   v1.7   wvd    Shift weights array to account for
 *					   FPA shifts.
 *                06/14/04   v1.8   wvd    Because the weights array is defined
 *					   wrt the FARF, we need not correct
 *					   for FPA shifts here.
 *                05/17/05   v1.9   wvd    Don't normalize regions where the
 *					   total probability is less than 0.1.
 *					   Set them to zero.
 *		  06/15/05   v1.10  wvd	   BUG FIX: Program always read the
 *					   probability arrays for point-source 
 *					   targets from the WGTS_CAL file.  Now
 *					   it distinguishes between point-source
 *					   and extended targets.
 *
 *****************************************************************************/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

int cf_rebin_probability_array(fitsfile *infits, int extended, int aperture,
	long nout, float *wave_out, int *pny, float *pycent, float **parray)
{

  char CF_PRGM_ID[] = "cf_rebin_probability_array";
  char CF_VER_NUM[] = "1.10";

  fitsfile *wgtsfits, *wavefits ;
  int status=0, fpixel=1, nullval=0, anynull=0 ;
  int hdu, nxwt, nywt, dndx, wtbin;
  long npix, npix_out, nwave, ndx, ndx1, ndx2, i, j ; 
  float *wave, *wavet, *wgts_arr, *ynorm, *outarr ;
  float total, wgt_cent ;
  char wgtsfile[FLEN_VALUE]={'\0'}, wavefile[FLEN_VALUE]={'\0'} ;
  char buffer[FLEN_CARD], det[FLEN_CARD] ;

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");

  /* Read calibration file names from the IDF file header. */
  FITS_movabs_hdu(infits, 1, NULL, &status) ;
  FITS_read_key(infits, TSTRING, "WGTS_CAL", wgtsfile, NULL, &status) ;
    cf_verbose(3,"weights cal file = %s ",wgtsfile) ;
  FITS_read_key(infits, TSTRING, "WAVE_CAL", wavefile, NULL, &status) ;
    cf_verbose(3, "wavelength cal file = %s ",wavefile) ;
  FITS_read_key(infits, TSTRING, "DETECTOR", det, NULL, &status) ;
    cf_verbose(3,"Detector = %s ",det) ;

  /* Open the weights file and read the array */
  FITS_open_file(&wgtsfits, cf_cal_file(wgtsfile), READONLY, &status) ;
  hdu = extended * 8 + aperture + 1;
  cf_verbose(3, "WGTS_CAL: extended = %d, aperture = %d, hdu = %d",
	extended, aperture, hdu);
  FITS_movabs_hdu(wgtsfits, hdu, NULL, &status) ;
  FITS_read_key(wgtsfits, TINT, "NAXIS1", &nxwt, buffer, &status) ;
  FITS_read_key(wgtsfits, TINT, "NAXIS2", &nywt, buffer, &status) ;
  FITS_read_key(wgtsfits, TFLOAT, "WGT_CENT", &wgt_cent, buffer, &status) ;
  npix=nxwt * nywt ;
  cf_verbose(3, "prob array: nx=%d, ny=%d, wgt_cent=%f ",nxwt, nywt, wgt_cent);
  wgts_arr = (float *) calloc(npix, sizeof(float)) ;
  FITS_read_img(wgtsfits, TFLOAT, fpixel, npix, &nullval, wgts_arr, 
		 &anynull, &status) ;
  FITS_close_file(wgtsfits, &status) ;

  /* Read the wavelength array */
  FITS_open_file(&wavefits, cf_cal_file(wavefile), READONLY, &status) ;
  FITS_movabs_hdu(wavefits, aperture +1, NULL, &status) ;
  nwave=cf_read_col(wavefits,TFLOAT,"WAVELENGTH", (void **) &wave);
  FITS_close_file(wavefits, &status) ;

  /* Determine the binning factor for the weights and generate a new
     wavelength array that corresponds to the binning factor used. */
	wtbin = NXMAX / nxwt ;
	cf_verbose(3, "Binning factor for weights array = %d ", wtbin) ;
	wavet = (float *) cf_calloc(nxwt, sizeof(float) ) ;
	for (i=0; i<nxwt; i++) {
	    ndx=i * wtbin + (wtbin/2);
	    if (ndx > nwave-1) ndx = nwave-1 ;
	    if (ndx < 0) ndx = 0 ;
	    wavet[i] = wave[ndx] ;
	}

  /* Allocate space for the output array */
  npix_out = nout * nywt ;
  if ((aperture < 4 && strncmp(det,"1",1) == 1) ||
      (aperture > 4 && strncmp(det,"2",1) == 1)) 
      cf_verbose(3, "maximum tabulated wavelength = %g, "
	"maximum of wave_out=%g", wavet[nxwt-1],wave_out[nout-1] );
  else
      cf_verbose(3, "maximum tabulated wavelength = %g, "
	"maximum of wave_out=%g", wavet[0],wave_out[nout-1] );
  outarr = (float *) cf_calloc(npix_out, sizeof(float)) ;

  /* Fill in the output array - use the weight profile which is closest
     to the individual wavelengths in the output wavelength array. 
     Remember : the Lif channels and the SiC channels have wavelengths
     that run in opposite directions, e.g. for side 1 wavelengths
     increase with index for Lif but decrease for SiC. The opposite is
     true for side 2. */

  ndx=0 ;
  dndx = 1 ; 
  if ((aperture < 4 && strncmp(det,"2",1) == 0) ||
      (aperture > 4 && strncmp(det,"1",1) == 0 ) ) {
    ndx = nxwt - 1 ;
    dndx = -1 ;
  }

  cf_verbose(3, "filling output: starting index = %d, increment = %d ",
	ndx, dndx) ;

  for (i=0; i<nout; i++) {
    /* locate the wavelength corresponding to this entry */
        while (wavet[ndx] < wave_out[i] && ndx <= nxwt-1 && ndx >= 0 )
		ndx += dndx ;
	/* populate the column of the array */
	for (j=0; j<nywt ; j++) {
          ndx1 = j*nout + i ;
	  if (ndx1 > npix_out-1) ndx1 = npix_out-1 ;
	  ndx2 = j*nxwt + ndx ;
	    if (ndx2 > npix-1) ndx2=npix-1 ;
	  outarr[ndx1] = wgts_arr[ndx2] ; 
	}
  }

  /* Adjust the probabilities so that the integral along y at every 
     wavelength position is 1 */

  /* First sum the points along y (ynorm) at each point */
  ynorm= (float *) cf_calloc(nout, sizeof(float)) ;
  for (i=0; i<nout; i++) {
    total=0. ;
    for (j=0; j<nywt; j++) total += outarr[j*nout+i] ;
    ynorm[i]=total ;
  }

  /* Normalize all points by dividing by the total */
  for (i=0; i<nout; i++) {
     if (ynorm[i] > 0.1)
	for (j=0; j<nywt; j++) outarr[j*nout+i] /= ynorm[i] ;
     else
	for (j=0; j<nywt; j++) outarr[j*nout+i] = 0. ;
  }
          
  /* Redirect the output pointers. */
  *pny = nywt ; 
  *pycent = wgt_cent ;
  *parray = outarr ; 

  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return status;
}
