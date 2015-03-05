/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_bad_pixels options input_file
 *
 * Description: Determines a map of the bad pixels on the detector, after
 *		correction for image motions during the observation. The
 *		procedure works by generating a series of pseudo photons
 *		corresponding to bad pixels within the area of the detector
 *		seen by the target aperture. This list of pseudo photons is
 *		repeated for each second of the observation where the data
 *		is expected to be acceptable (i.e., it excludes times of
 *		bursts, SAA passages, etc). The list is then run through
 *		various routines in the pipeline which correct for photon
 *		motions. A map is then generated from the final photon
 *		list.  
 *
 * Arguments:   input_file              FARF-corrected intermediate data
 *                                      file
 * Calibration files:  QUAL_CAL                 
 *
 * Returns:     file containing the bad pixel map 
 *
 *
 * History:     04/08/03   1.0   rdr   Begin work
 *		04/30/03   1.1   wvd   Install
 *		05/22/03   1.2   wvd   Modify call to cf_set_photon_flags
 *              05/28/03   1.3   rdr   Made adjustments for HIST data and 
 *                                      correct some bugs
 *		05/04/03   1.4   wvd   Reset JITR_COR and FLAG_COR to PERFORM.
 *				 	Pass weights to cf_set_photon_flags.
 *              09/06/03   1.5   rdr   Incorporate tscreen flag.
 *		06/11/03   1.6   wvd   Pass datatype to cf_read_col and
 *					cf_write_col.
 *		08/22/03   1.7   wvd   Change channel array to unsigned char.
 *					Delete GTI from call to
 *					cf_satellite_jitter.
 *		08/25/03   1.8   wvd   Change coltype from string to int in
 *					cf_read_col.
 *              10/17/03   1.9   wvd   If EXPTIME = 0, exit without
 *                                      generating a bad-pixel file.
 *              11/05/03   1.10  wvd   Change chan_pix to unsigned char
 *              11/05/03   1.11  wvd   Use fits_create_tbl alone to make
 *					binary table extension.
 *              12/01/03   1.12  bjg   Update FILENAME, FILETYPE and IDF_FILE
 *                                      keywords.
 *              12/01/03   1.13  bjg   Minor changes.
 *              12/21/03   1.14  wvd   Remove underscore from idf and bpm
 *					filenames.
 *              02/24/04   1.15  rdr   change maximum size of the output
 *                                     arrays in cf_combine_pothole_data
 *              03/03/04   1.16  rdr   make sure that lif_cnt and sic_cnt arrays 
 *                                     are non-zero
 *              04/06/04   1.19  bjg   Include ctype.h and strings.h
 *                                     Change formats to match arg types in 
 *                                     printf
 *                                     Change ( a<b & a>c ) to ((a<b)&&(a>c)) 
 *                                     in cf_combine_pothole_data
 *                                     Create potholes from HIST FILL_DATA
 *              06/07/04   1.20  wvd   Add exp_jitr to cf_satellite_jitter(),
 *				   	delete timeline from cf_apply_filters.
 *              07/21/04   1.21  wvd   Delete unused variable nmax in
 *					cf_generate_pseudo_photons
 *              10/29/04   1.22  bjg   Added check for selected potholes that
 *                                     fall out the extraction window
 *              12/02/04   1.23  wvd   If a processing step was skipped for
 *                                      the data file, skip it also for the
 *                                      bpm file.
 *              01/27/05   1.24  wvd   Ignore LIMB, SAA, and BRST when 
 *					determining good times for HIST data.
 *		02/01/05   1.25  wvd   In cf_generate_pseudo_photons, pad
 *					extraction windows by ten pixels
 *					when making pothole list.  
 *				       In cf_combine_pothole_data, toss any
 *					potholes that drift out of extraction
 *					window.  Call cf_get_extraction_limits
 *					only when the aperture changes.
 *		02/17/05   1.26  wvd   Increase estimated size of output 
 *					pixel arrays (nmax) by 20%.
 *					Add additional diagnostic output.
 *		03/15/05   1.27  wvd   Replace cf_get_extraction_limits with
 *					cf_extraction_limits
 *		03/16/05   1.28  wvd   Pass srctype to cf_extraction_limits.
 *              03/22/05   1.29  wvd   Change TIME_SUNRISE and TIME_SUNSET
 *                                      from floats to shorts.
 *              03/22/05   1.30  wvd   In cf_combine_pothole_data, use fabs()
 *					when looking for changes in the time
 *					array.
 *              04/12/05   1.31  wvd   Add -n flag, which forces creation
 *                                      of night-only BPM file.  Its argument
 *                                      is the name of the output BPM file.
 *					This file name is NOT written to the
 *					IDF file header.
 *              06/03/05   1.32  wvd   Include math.h
 *              11/25/05   1.33  wvd   Use cf_x2lambda and cf_get_potholes
 *					rather than writing new routines here.
 *					Don't add random numbers to output
 *					pixel coordinates.
 *              05/15/06   1.34  wvd    Divide cf_astigmatism_and_dispersion
 *                                      into two routines.  Note that wave-
 *                                      length calibration is performed even
 *                                      if astigmatism correction is not.
 *					Add ASTG_COR to list of updated
 *					keywords.
 *              06/13/06   1.35  wvd    In cf_combine_pothole_data,
 *					initialize nfill to 0 before each
 *					iteration and ignore bad pixels that
 *					lie outside of the extraction window.
 *              03/20/07   1.36  wvd    In cf_combine_pothole_data, if the
 *					entire dead spot falls outside of the
 *					aperture, stop and move on to the
 *					next one.
 *              03/21/07   1.37  wvd    Replace "break" with "continue" in
 *					cf_combine_pothole_data.
 *		04/03/07   1.38  wvd	Scale nmax by 1.5.  Call cf_error_init
 *					after each call to cf_extraction_limits.
 *					Make CF_PRGM_ID and CF_VER_NUM static
 *					variables.
 *		04/03/07   1.39  bot	In cf_generate_pseduo_photons, l.198
 *					and l.219 added a test on LOCATION_SHLD
 * 					for FILL DATA.
 *		07/18/08   1.40  wvd	If EXP_STAT = 2, target is bright
 *					earth/airglow. Mask out limb-angle flag.
 *
 ****************************************************************************/

#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_bad_pixels"; 
static char CF_VER_NUM[] = "1.40";


/****************************************************************************
 *
 *              CF_GENERATE_PSEUDO_PHOTONS
 *
 *   Procedure to generate an array of pseudo photons to represent the
 *   locations of the detector potholes
 *
 ****************************************************************************/

int cf_generate_pseudo_photons(fitsfile *header, long ngood, long *good_index, 
        float *timeline, unsigned char *statflag, float *lif_cnt,
	float *sic_cnt, long *nevents, float **time, float **weight,
	float **x, float **y, unsigned char **channel, float **rx,
	float **ry, unsigned char **timeflag, unsigned char **loc_flag) {

  char instmode[FLEN_CARD] ;
  long nevts, npot, nfill=0, npot_sel; 
  int *bpndx, npot2, status=0, hdutype ;
  int srctype, ap[2], bymax, bymin, pad=10 ;
  short *ylow, *yhigh; 
  long npts, i, j, k, xndx, num, num_tot, ndx, tndx ;
  float *xpot, *ypot, *rxpot, *rypot ;
  float *xpot2, *ypot2, *rxpot2, *rypot2 ;
  float *xval, *yval, *rxval, *ryval, *xout, *yout, *rxout, *ryout ;
  float *time_out, *weight_out, tval, ctot_lif, ctot_sic ;
  float wval_lif, wval_sic ;
  short binx, biny, xmin, xmax;
  unsigned char *chan, *chan_out;
  unsigned char *tflag_out, tflag_val ;
  char * hdu2_loc_flgs;
  float * hdu2_xfarf, * hdu2_yfarf;


  cf_verbose(3, "Generating pseudo-photons") ;

  /* read the instrument mode for the data */
    FITS_movabs_hdu(header, 1, &hdutype, &status);
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status) ;
    FITS_read_key(header, TSHORT, "SPECBINX", &binx, NULL, &status) ;
    FITS_read_key(header, TSHORT, "SPECBINY", &biny, NULL, &status) ;

  /* read the header and determine the target apertures for the exposure */
  srctype = cf_source_aper(header, ap) ;
  cf_verbose(3,"Target apertures = %d and %d ",ap[0],ap[1]) ;

  /* read in the pothole positions */
  npot2 = cf_get_potholes(header, &xpot2, &ypot2, &rxpot2, &rypot2) ; 
  cf_verbose(3, "number of potholes found = %d ",npot2) ;

  if (!(strncasecmp(instmode,"HIST",4))) {
    FITS_movabs_hdu(header, 2, &hdutype, &status);
    nevts=cf_read_col(header,TBYTE,"LOC_FLGS",(void **) &hdu2_loc_flgs);
    nevts=cf_read_col(header,TFLOAT,"XFARF",(void **) &hdu2_xfarf);
    nevts=cf_read_col(header,TFLOAT,"YFARF",(void **) &hdu2_yfarf);
    for (j=0;j<nevts;++j) {
      if (!((hdu2_loc_flgs[j] & LOCATION_FILL) == 0) &&
	   ((hdu2_loc_flgs[j] & LOCATION_SHLD) == 0)) ++nfill;
    }
    cf_verbose(3, "number of FILL DATA potholes found = %d ",nfill) ;
    npot=npot2+nfill;

    xpot=(float *) cf_calloc(npot, sizeof(float)) ;
    ypot=(float *) cf_calloc(npot, sizeof(float)) ;
    rxpot=(float *) cf_calloc(npot, sizeof(float)) ;
    rypot=(float *) cf_calloc(npot, sizeof(float)) ;

    for (j=0;j<npot2;++j){
      xpot[j]=xpot2[j];
      ypot[j]=ypot2[j];
      rxpot[j]=rxpot2[j];
      rypot[j]=rypot2[j];
    }

    k=npot2;

    for (j=0;j<nevts;++j){
      if (!((hdu2_loc_flgs[j] & LOCATION_FILL) == 0) &&
	   ((hdu2_loc_flgs[j] & LOCATION_SHLD) == 0)) {
	
	xpot[k]=hdu2_xfarf[j];
	ypot[k]=hdu2_yfarf[j];
	rxpot[k]=binx;
	rypot[k]=biny;
	k++;

      }
    }
    free(hdu2_loc_flgs) ;
    free(hdu2_xfarf) ;
    free(hdu2_yfarf) ;
  }
  else {
    npot=npot2;
    xpot=xpot2;
    ypot=ypot2;
    rxpot=rxpot2;
    rypot=rypot2;
  }  

  /* allocate space to hold the pothole list (for 1 second of data) */
  xval = (float *) cf_calloc(npot, sizeof(float)) ;
  yval = (float *) cf_calloc(npot, sizeof(float)) ;
  rxval = (float *) cf_calloc(npot, sizeof(float)) ;
  ryval = (float *) cf_calloc(npot, sizeof(float)) ;
  chan = (unsigned char *) cf_calloc(npot, sizeof(char)) ;
  num=-1 ;


    cf_verbose(3, "For initial pothole selection, "
	"we pad extraction windows by %d pixels in Y.", pad);

  /* do each aperture separately */

  for (j=0 ; j<2; j++) {

  /* get the extraction limits for the relevant apertures */
    npts = cf_extraction_limits(header, ap[j], srctype,
	&ylow, &yhigh, &xmin, &xmax) ;
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* select potholes located within the extraction window */
    bpndx = (int *) cf_calloc(npot, sizeof(int)) ;
    npot_sel = -1 ; 
    for (i=0; i<npot; i++) {
      xndx = (int) (xpot[i] + 0.5) ;
      bymax = (int) (ypot[i]+rypot[i]+0.5) ;
      bymin = (int) (ypot[i]-rypot[i]) ;
      if ( bymax >= ylow[xndx] - pad  && bymin <= yhigh[xndx] + pad ) {
        npot_sel += 1 ;
        bpndx[npot_sel] = i ;
	cf_verbose(3, "pothole %d selected: x=%d, y=%d, rx=%d, ry=%d ", i,
	cf_nint(xpot[i]),cf_nint(ypot[i]),cf_nint(rxpot[i]),cf_nint(rypot[i]));
      }
    }

    /* convert npot_sel from an array index to a number of potholes selected */
    npot_sel++ ;
    cf_verbose(3, "number of potholes selected for ap %d = %d ",
	ap[j], npot_sel) ;

    /* put the locations into the output arrays */
    for (i=0; i<npot_sel; i++) {
      num+=1 ;
      ndx = bpndx[i] ;
      xval[num] = xpot[ndx] ;
      yval[num] = ypot[ndx] ; 
      chan[num] = ap[j] ;
      rxval[num] = rxpot[ndx] ;
      ryval[num] = rypot[ndx] ;
    }
  }

  /* convert num from an index to the number of potholes found */
  num+=1 ;
  cf_verbose(3, "total number of potholes found = %d ",num) ;

  /* generate a timeline with only the good times included */
  num_tot = num * ngood ;

  /* copy the non-zero values into the output arrays */
  time_out = (float *) cf_calloc(num_tot, sizeof(float)) ;
  weight_out = (float *) cf_calloc(num_tot, sizeof(float)) ;
  xout = (float *) cf_calloc(num_tot, sizeof(float)) ;
  yout = (float *) cf_calloc(num_tot, sizeof(float)) ;
  rxout = (float *) cf_calloc(num_tot, sizeof(float)) ;
  ryout = (float *) cf_calloc(num_tot, sizeof(float)) ;
  chan_out = (unsigned char *) cf_calloc(num_tot, sizeof(char)) ;
  tflag_out = (unsigned char *) cf_calloc(num_tot, sizeof(char)) ;
  *loc_flag = (unsigned char *) cf_calloc(num_tot, sizeof(char)) ;

  tndx=-1 ;
  ctot_lif = 0. ;
  ctot_sic = 0. ;
  for (j=0; j<ngood; j++) {
    ndx=good_index[j] ;
      tval = timeline[ndx] ;
      ctot_lif += lif_cnt[ndx] ;
      ctot_sic += sic_cnt[ndx] ;
      wval_lif = lif_cnt[ndx] ;
      wval_sic = sic_cnt[ndx] ;
      tflag_val = statflag[ndx] ;
    for (i=0; i<num; i++) {
      tndx += 1;
      if (tndx > num_tot -1) {
          cf_verbose(2, "overflow tndx while filling pothole array ") ;
          tndx= num_tot-1 ; }
      time_out[tndx] = tval ;
        if (chan[i] < 4) weight_out[tndx] = wval_lif ;
        else weight_out[tndx] = wval_sic ;
      xout[tndx] = xval[i] ;
      yout[tndx] = yval[i] ;
      rxout[tndx] = rxval[i] ;
      ryout[tndx] = ryval[i] ;
      chan_out[tndx] = chan[i] ;
      tflag_out[tndx] = tflag_val ;
    }
  }

  /* normalize the weights for TTAG data and set to 1 for HIST data*/
  if (!strncmp(instmode,"T",1)) 
    for (i=0; i<num_tot; i++) {
      if (chan_out[i] < 4) weight_out[i] /= ctot_lif ;
      else weight_out[i] /= ctot_sic ;
  }
  else for (i=0; i<num_tot; i++) weight_out[i] = 1. ;

  /* assign pointers to output variables */
  *x = xout ;
  *y = yout ;
  *rx = rxout ;
  *ry = ryout ;
  *channel = chan_out ;
  *timeflag = tflag_out ;
  *nevents = num_tot ;
  *weight = weight_out ;
  *time = time_out ;

  free(xval) ;
  free(yval) ;
  free(rxval) ;
  free(ryval) ;
  free(chan) ;
  free(bpndx) ;

  return 0;

}


/*******************************************************************************
 *
 *                     CF_CORRECT_DOPPLER_MOTIONS
 *  
 *    procedure to correct the pothole centroid positions for doppler effects
 *
 ******************************************************************************/

 int cf_correct_doppler_motions(fitsfile *header, long nevents, float *time,
	float *x, unsigned char *channel, long nseconds, float *timeline,
	float *velocity)
{
    char  	wave_file[FLEN_FILENAME];
    float	*wavelength, v_helio, vel, lam, dlam, dx;
    fitsfile 	*wavefits;
    int   	status=0, anynull;
    int	  	fcol, targ_ap[2], ap, ii;
    long	i, k, ndx ;

    /* get the information needed in the analysis */

    /* first determine the target apertures */
    cf_source_aper(header, targ_ap) ;
    cf_verbose(3, "target apertures = %d and %d ",targ_ap[0], targ_ap[1]) ;

    /* read in the wavelength calibration */
    FITS_read_key(header, TSTRING, "WAVE_CAL", wave_file, NULL, &status);
    FITS_open_file(&wavefits, cf_cal_file(wave_file), READONLY, &status);
    wavelength = (float *) cf_calloc(NXMAX, sizeof(float));

    /* read the heliocentric velocity from the header */
     FITS_read_key(header, TFLOAT, "V_HELIO", &v_helio, NULL, &status);
    cf_verbose(3, "heliocentric velocity = %f ",v_helio) ;

    /* Go through each target aperture  */
    for (i = 0; i < 2; i++) {
      ap = targ_ap[i] ;

    	FITS_movabs_hdu(wavefits, ap+1, NULL, &status);
        FITS_get_colnum(wavefits, TRUE, "WAVELENGTH", &fcol, &status);
        fits_read_col(wavefits, TFLOAT, fcol, 1, 1, NXMAX, NULL, wavelength, 
			&anynull, &status);

	ndx = 0 ;
       for (k = 0; k < nevents; k++) 
	    if (channel[k] == ap) {
	      while (time[k] >= timeline[ndx] )	ndx++ ; 
	      if (ndx > nseconds-1 ) ndx=nseconds-1 ;
	        ii = (int) (x[k] + 0.5);
		if (ii >= 1 && ii < NXMAX) {
	            lam = wavelength[ii];
                    vel = velocity[ndx] + v_helio ;
	            dlam = vel * lam / C ;
                    dx = dlam / ( wavelength[ii] - wavelength[ii-1]) ;
		    x[k] += dx ;
		    if (vel > 100 || vel < -100 ) 
                         cf_if_warning ("at k=%ld, time=%f, lam=%f, vel=%f, "
				"dlam=%f, dx=%f \n",
				 k, time[k], lam, vel, dlam, dx) ;
		}
	    }

    }
    free(wavelength);
    FITS_close_file(wavefits, &status);

    return status;
}


/*******************************************************************************
 *
 *                     CF_COMBINE_POTHOLE_DATA
 *
 *    procedure to combine the centroid data as a function of time with the 
 *    pothole size information to generate a list of affected pixels on the
 *    detector
 *
 ******************************************************************************/

 long cf_combine_pothole_data(fitsfile *header, long nevents, float *time, 
        float *weight, float *x, float *rx, float *y, float *ry,
	unsigned char *channel, unsigned char *timeflag, float **xpix,
	float **ypix, unsigned char **chan_pix, float **wt_pix) {

   int srctype, expstat, aperture[2];
   int  npot, xndx, yndx, ap, nfill ;
   short sxmin, sxmax, *ylow, *yhigh ;
   int  xdim1, ydim1, xmin1, xmax1, ymin1, ymax1 ;
   int xminp, xmaxp, yminp, ymaxp ;
   long i, j, k, l, m, num, numt ;
   long nout, nmax, npts1, npts2, xdim, ydim, ndx, ndxp ;
   float *xpixt, *ypixt, *wt_pixt, xmin, xmax, ymin, ymax ;
   float *array1, *array2, wtot, wval, wmax ;
   float xv2, yv2, rx2, ry2, pos, xval, yval ;
   float rxp, ryp, xv, yv, time0 ;  
   unsigned char *chan_pixt;

   int status=0;
   char instmode[FLEN_VALUE];
   unsigned char TEMPORAL_MASK;

   cf_verbose(3, "Combining pothole data") ;

    /* 
     *  Read INSTMODE keyword.  If in HIST mode, mask out
     *  LIMB, SAA, and BRST flags.  TEMPORAL_DAY is the default.
     */
    TEMPORAL_MASK = TEMPORAL_DAY;
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    if (!strncmp(instmode, "HIST", 4)) {
	TEMPORAL_MASK |= TEMPORAL_LIMB;
	TEMPORAL_MASK |= TEMPORAL_SAA;
	TEMPORAL_MASK |= TEMPORAL_BRST;
    }
    /*
     * If EXP_STAT = 2, target is bright earth/airglow.  Mast out limb-angle flag.
     */
    FITS_read_key(header, TINT, "EXP_STAT", &expstat, NULL, &status);
    if (expstat == 2) TEMPORAL_MASK |= TEMPORAL_LIMB;

   /* Read source type from header. */
   srctype = cf_source_aper(header, aperture) ;

   /* determine the number of potholes by looking at the number
	of entries for time=time[0] */

   time0 = time[0] ;
   i = 1 ;
   while (fabs((time[i] - time0)) < 0.5) i++ ;
   cf_verbose(3, "Number of potholes found = %d ", npot = i) ;

   nmax = 0;
   /* estimate size of the array needed to store the data and open
      arrays to contain it */   
    for (i=0; i<npot; i++) nmax += (60.+2.*rx[i]) * (30.+2*ry[i]) ;
	nmax *= 1.5;
      cf_verbose(3, "Estimated storage space needed for output array = %d ",
	nmax) ;
      xpixt = (float *) cf_calloc(nmax, sizeof(float)) ;
      ypixt = (float *) cf_calloc(nmax, sizeof(float)) ;
      wt_pixt = (float *) cf_calloc(nmax, sizeof(float)) ;
      chan_pixt = (unsigned char *) cf_calloc(nmax, sizeof(char)) ;

      /* initialize nout = total number of pixels affected by potholes */
      nout = -1 ;

      /* process each pothole individually */ 

      ap = -1;
      for (i=0; i<npot; i++) {

	nfill = 0;

	/* If channel has changed, get the new extraction limits */
	if (ap != (int) channel[i]) {
	    ap = (int) channel[i] ;
            (void) cf_extraction_limits(header, ap, srctype,
		&ylow, &yhigh, &sxmin, &sxmax);
    	    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
	}
	
    	/* determine the range in X and Y for the pothole centroid positions */
	xmin = 17000. ;
        xmax = 0. ;
        ymin = 1024. ;
        ymax = 0. ;
        numt = 0 ;
        for (j=i; j<nevents ; j+= npot) { 
	  numt ++ ;
	     if ((x[j] > xmax) && (x[j] < 16383)) xmax = x[j] ;
             if ((x[j] < xmin) && (x[j] > 0)) xmin = x[j] ;
             if ((y[j] > ymax) && (y[j] < 1023)) ymax = y[j] ;
             if ((y[j] < ymin) && (y[j] > 0)) ymin = y[j] ;
	}
      cf_verbose(3,"For pothole %.1d, xmin=%.1d, "
	"xmax=%.1d, ymin=%.1d, ymax=%.1d",
	i, cf_nint(xmin), cf_nint(xmax), cf_nint(ymin), cf_nint(ymax)) ;

      /* set up a 2D array which contains all of the movement in the pothole */
      xdim = (long) (xmax - xmin + 1.5) ;
      ydim = (long) (ymax - ymin + 1.5) ;
      npts1 = xdim * ydim  ;
      cf_verbose(3, "Setting up array1: xdim=%d, ydim=%d, npts=%d ",
		xdim, ydim, npts1) ;
      array1 = (float *) cf_calloc(npts1, sizeof(float)) ;

      /* Determine the 2D distribution of positions for the pothole */
      num = 0 ;
      for (j=i ; j<nevents ; j+= npot) 
        if (!(timeflag[j] & ~TEMPORAL_MASK)) {
	  num++ ;
	xndx = cf_nint(x[j] - xmin) ;
	  if (xndx < 0 ) xndx = 0 ;
          if (xndx > xdim-1) xndx = xdim-1 ;
        yndx = cf_nint(y[j] - ymin) ;
    	  if (yndx < 0 ) yndx = 0 ;
          if (yndx > ydim-1) yndx = ydim-1 ;
        ndx = xndx + yndx*xdim ;
	  if (ndx > npts1-1) ndx=npts1-1 ;
        array1[ndx] += weight[j] ;
      }

      if (i == 0 ) cf_verbose(3, "%d good seconds out of a total of %d  ",
		num, numt) ;

      wtot=0. ;
      for(j=0; j<ydim; j++)
        for (k=0; k<xdim; k++)  wtot += array1[j*xdim + k] ;
      cf_verbose(3, "wtot for array 1= %d", cf_nint(wtot)) ;

      /* specify the size of the pothole */
         rxp = rx[i]  ;
         rx2 = rxp * rxp ;
         ryp = ry[i]  ;
         ry2 = ryp * ryp ;
	 cf_verbose(3,"pothole size: rx=%d, ry=%d ",cf_nint(rxp),cf_nint(ryp));

      /* set up a target array to contain the image of the full pothole, after
         reconstruction */
      xmax1 = xmax + rxp ;
      xmin1 = xmin - rxp ;
        xdim1 = (int) (xmax1 - xmin1 + 1.5) ;
        xndx = (long) x[i] ;
      ymax1 = ymax + ryp ;
      ymin1 = ymin - ryp ;
        cf_verbose(3,"area of the detector covered by pothole :") ;
	cf_verbose(3, "     xmin=%d, xmax=%d, ymin=%d, ymax=%d ",
		xmin1,xmax1,ymin1,ymax1) ;
        cf_verbose(3,"extraction window extends from %d to %d ",ylow[xndx], yhigh[xndx]) ;

	/* Include only pixels that lie within extraction window. */
	if (ymin1 < ylow[xndx])  ymin1 = ylow[xndx];
	if (ymax1 > yhigh[xndx]) ymax1 = yhigh[xndx];
        ydim1 = (int) (ymax1 - ymin1 + 1.5) ;

	/* If no pixels overlap the extraction window, skip to the next one. */
	if (ydim1 <= 0) {
	   cf_verbose(3, "This pothole does not overlap the extraction window.  We'll skip it.");
	   continue;
	}

      npts2 = xdim1 * ydim1 ;
      cf_verbose(3, "setting up array2: xdim = %d, ydim = %d, npts2 = %d ", 
         xdim1, ydim1, npts2) ;
      array2 = (float *) cf_calloc(npts2, sizeof(float)) ;

      /* fill in array2
         - scan all of the positions specified in array1  */
      for (j=0 ; j<ydim ; j++) 
	for (k=0 ; k<xdim ; k++){
	  /* specify the properties of the pothole at this location */
	  xval = k + xmin ;
          yval = j + ymin ;
	  ndx = j*xdim + k ;
	  if (ndx > npts1-1) ndx=npts1-1 ;
            wval = array1[ndx] ;

	if (wval > 0) {
	  /* set limits to the region of array2 affected by the pothole */
	  xminp = (int) (xval - rxp - xmin1 + 0.5)  ;
          xmaxp = (int) (xval + rxp - xmin1 + 0.5)  ;
          yminp = (int) (yval - ryp - ymin1 + 0.5) ;
          ymaxp = (int) (yval + ryp - ymin1 + 0.5)  ;
	    if (yminp < 0) yminp = 0 ;
            if (ymaxp > ydim1-1) ymaxp = ydim1-1 ;
            if (xminp < 0) xminp = 0 ;
            if (xmaxp > xdim1-1) xmaxp = xdim1-1 ;

	/* fill in the array for this pothole location */
	  for (l=yminp; l<=ymaxp; l++) 
            for (m=xminp; m <= xmaxp ; m++) {
              ndxp = l * xdim1 + m ;
	         if (ndxp > npts2-1) ndxp = npts2 -1 ;
	      xv = m + xmin1 - xval ;
              yv = l + ymin1 - yval ;
              xv2= xv * xv ;
              yv2 = yv * yv ;
              pos = (xv2/rx2) + (yv2/ry2)  ; 
	      if ( pos <= 1)  array2[ndxp] += wval ; 
	      nfill++;
	    }
	 } 
      }

	cf_verbose(3,"number of pixels affected by this pothole = %d ",nfill);

      /* determine the maximum weights across the pothole map */
      wmax=0. ;
      for(j=0; j<ydim1; j++)
        for (k=0; k<xdim1; k++) {
	  ndx = j * xdim1 + k ;
          if (array2[ndx] > wmax) wmax = array2[ndx] ;
	}
      cf_verbose(3, "wmax for array2 = %d ", cf_nint(wmax)) ;

      /* normalize the values of the affected pixels to the maximum and put 
         their locations into the output array */
      for(j=0; j<ydim1; j++)
        for (k=0; k<xdim1; k++) {
	  ndx = j * xdim1 + k ;
          if (array2[ndx] > 0) {
	    if (nout < nmax-1) {
              nout ++ ;
	      xpixt[nout] = k + xmin1;
              ypixt[nout] = j + ymin1;
              wt_pixt[nout] = array2[ndx]/wmax ;
              chan_pixt[nout] = channel[i] ;
            }
	    else {
	      cf_if_warning("Array overflow.  Truncating pseudo-photon list.");
	      break;
	    }
          }
	}


      cf_verbose(3, "%d points tabulated after processing pothole %d\n",nout,i);

      free(array1) ;
      free(array2) ;
      }

      /* assign output pointers */
      *xpix = xpixt ;
      *ypix = ypixt ;
      *wt_pix = wt_pixt ;
      *chan_pix = chan_pixt ;

      cf_verbose(3, "Total number of pixels tabulated = %d ", nout) ;

      return nout ;
}


/****************************************************************************/

int
main(int argc,  char *argv[])
{
    unsigned char  *timeflag=NULL, *statflag=NULL, *locflag=NULL;
    unsigned char *loc_flag=NULL, *channel=NULL ;
    int  grating_motion=1, mirror_motion=1, fpa_pos=1;
    int  jitter=1, doppler_motion=1, astig=1, tscreen=1, night_only=FALSE;
    int  status=0, hdutype, optc;
    long nevents, nseconds, i, nrows=1, npixels ;
    long ngood, *good_index, dtime, ntime ;
    float exptime, *time=NULL, *x=NULL, *y=NULL, *velocity=NULL;
    float *rx=NULL, *ry=NULL, *xpix, *ypix, *wt_pix ;
    float *lif_cnt=NULL, *sic_cnt=NULL, *lambda=NULL ;
    float *timeline=NULL, *weight=NULL;
    short *tsunrise=NULL, *tsunset=NULL ;
    fitsfile *header, *bpmfits, *memp;
    char rootname[FLEN_VALUE]={'\0'}, bpm_file[FLEN_FILENAME]={'\0'};
    char idf_file[FLEN_FILENAME]={'\0'}, det[FLEN_VALUE]={'\0'}; 
    char *keywords[] = {"GRAT_COR","FPA__COR","MIRR_COR","ASTG_COR","WAVE_COR"};
    char keyval[FLEN_VALUE];
    unsigned char *chan_pix ;

    char extname[]="POTHOLE_DATA", fmt_float[FLEN_VALUE], fmt_byte[FLEN_VALUE] ;
    char instmode[FLEN_VALUE] ;
    int tfields=5 ;
    char *ttype[]={ "X", "Y", "CHANNEL", "WEIGHT", "LAMBDA"} ;
    char *tform[5] ;  /* we will define the tform once the number of 
                         array elements is known */
    char *tunit[]={ "pixels", "pixels",  "unitless", "unitless", "Angstroms"} ;
    
    char opts[]  = "hgmfjdasn:v:";
    char usage[] = 
	"Usage:\n"
	"  cf_bad_pixels [-hvgmfjdas] [-n bpm_filename] [-v level] idffile\n";
    char option[] =
	"Options:\n"
	"  -h:  this help message\n"
        "  -v:  verbosity level (=1; 0 is silent)\n" 
        "  -g:  no grating motion correction \n"
        "  -m:  no mirror motion correction \n" 
        "  -n:  Create BPM file for nighttime fraction of exposure.\n" 
	"       Argument is the name of the output BPM file.\n"
        "  -f:  no fpa motion correction \n" 
        "  -j:  no jitter correction \n" 
        "  -d:  no doppler correction \n" 
        "  -a:  no astigmatism correction \n"
        "  -s:  no screening on time flags (if 0)" ;

    verbose_level = 1;

    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
	case 'v':
	   verbose_level = atoi(optarg);
	    break;
        case 'g':
	    grating_motion = 0 ;
            break;
        case 'm':
	  mirror_motion = 0 ;
          break ;
        case 'n':
            sprintf(bpm_file,"%s", optarg) ;
            night_only = TRUE;
            break;
        case 'f':
          fpa_pos = 0 ;
          break ;
        case 'j':
          jitter = 0 ;
          break ;
        case 'a':
          astig = 0 ;
          break;
        case 'd':
          doppler_motion = 0 ;
          break ;
        case 's':
          tscreen = 0 ;
          break ;

	}
    }

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    if (argc <= optind) {
        printf("%s", usage);
	return -1;
    }

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Open the input IDF file. */
    FITS_open_file(&header, argv[optind], READWRITE, &status);
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status) ;

    /* If EXPTIME < 1, exit without generating a bad-pixel file. */
    FITS_read_key(header, TFLOAT, "EXPTIME", &exptime, NULL, &status);
    if (night_only)
	FITS_read_key(header, TFLOAT, "EXPNIGHT", &exptime, NULL, &status);
    if (exptime < 1.) {
        cf_verbose(1, "EXPTIME = %g.  "
	    "Will not attempt to generate bad-pixel map.", exptime);
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
        return 0;
    }

    /* Read timeline data from the input file */
    FITS_movabs_hdu(header, 4, &hdutype, &status);
    nseconds = cf_read_col(header, TFLOAT, "TIME",         (void **) &timeline);
    nseconds = cf_read_col(header, TBYTE,  "STATUS_FLAGS", (void **) &statflag);
    nseconds = cf_read_col(header, TSHORT, "TIME_SUNRISE", (void **) &tsunrise);
    nseconds = cf_read_col(header, TSHORT, "TIME_SUNSET",  (void **) &tsunset) ;
    nseconds = cf_read_col(header, TFLOAT, "ORBITAL_VEL",  (void **) &velocity);
    nseconds = cf_read_col(header, TFLOAT, "LIF_CNT_RATE", (void **) &lif_cnt) ;
    nseconds = cf_read_col(header, TFLOAT, "SIC_CNT_RATE", (void **) &sic_cnt) ;
    FITS_movabs_hdu(header, 1, NULL, &status);

   /*  
    *  If night-only spectrum is requested from the command line, copy 
    *  IDF header into memory, close the IDF, modify copy, and pass it
    *  to subsequent routines.
    */
    if (night_only) {
        cf_verbose(3, "Night-only spectrum requested from command line.");
	FITS_create_file(&memp, "mem://", &status);
	FITS_copy_hdu(header, memp, 0, &status);
	if (!strncmp(instmode, "H",1)) {
	    FITS_movabs_hdu(header, 2, NULL, &status);
	    FITS_copy_hdu(header, memp, 0, &status);
	}
	FITS_close_file(header, &status);
	header = memp;

	FITS_movabs_hdu(header, 1, NULL, &status);
        FITS_update_key(header, TFLOAT, "EXPTIME",  &exptime, NULL, &status);
        FITS_update_key(header, TSTRING, "DAYNIGHT", "NIGHT", NULL, &status);
    }
    else { 		/* Generate name of output file */
	FITS_read_key(header, TSTRING, "ROOTNAME", rootname, NULL, &status) ;
	FITS_read_key(header, TSTRING, "DETECTOR", det, NULL, &status) ;
	det[1]=tolower(det[1]) ;
	if (!strncmp(instmode, "T",1)) 
	    sprintf(bpm_file,"%11s%2sttagfbpm.fit",rootname,det) ;
	else
	    sprintf(bpm_file,"%11s%2shistfbpm.fit",rootname,det) ;
    }

    /* Write name of BPM file to header of IDF. */
    cf_verbose(3,"Output bad pixel map to file %s", bpm_file) ;
    FITS_update_key(header, TSTRING, "BPM_CAL", bpm_file, NULL, &status);

    /* Create the output file and copy header from the input */
    FITS_create_file(&bpmfits, bpm_file, &status) ;
    FITS_copy_hdu(header, bpmfits, 0, &status) ;

    FITS_update_key(bpmfits, TSTRING, "FILENAME", bpm_file, 
		    NULL, &status) ;
    FITS_update_key(bpmfits, TSTRING, "FILETYPE", "BAD PIXEL MAP", 
		    NULL, &status) ; 

    FITS_read_key(header, TSTRING, "FILENAME", idf_file, NULL, &status) ;
    FITS_update_key(bpmfits, TSTRING, "IDF_FILE", idf_file, 
		    NULL, &status) ;
     
    /* Reset analysis keywords in the top level of the output file */
    FITS_movabs_hdu(bpmfits, 1, &hdutype, &status);
    for (i = 0; i < 5; i++) {
        FITS_read_key(bpmfits, TSTRING, keywords[i], keyval, NULL, &status);
        if (!strncmp(keyval, "COMPLETE", 8))
            FITS_update_key(bpmfits, TSTRING, keywords[i], "PERFORM", NULL,
		&status) ;
        else
            FITS_update_key(bpmfits, TSTRING, keywords[i], "OMIT", NULL,
		&status) ;
    }

    /* Make sure that the count rates are non-zero */
    for (i=0; i<nseconds; i++) {
      if (lif_cnt[i] < 0.01) lif_cnt[i] = 0.01 ;
      if (sic_cnt[i] < 0.01) sic_cnt[i] = 0.01 ;
    }

    /* generate a null array for locflag */
    locflag = (unsigned char *) cf_calloc(nseconds, sizeof(unsigned char)) ;

    /* select only the good times */ 
    FITS_movabs_hdu(header, 1, &hdutype, &status);
    cf_apply_filters(header, tscreen, nseconds, statflag, locflag, nseconds,
		     statflag, &dtime, &ntime, &ngood, &good_index) ;

    free(locflag) ;

    cf_verbose(3, "nseconds=%d, dtime=%d, ntime=%d, ngood=%d ",
	   nseconds, dtime, ntime, ngood) ;

    /* Generate the pseudo photons */
    cf_generate_pseudo_photons(header, ngood, good_index, timeline, statflag,
			       lif_cnt, sic_cnt, &nevents, &time, &weight, &x,
			       &y, &channel, &rx, &ry, &timeflag, &loc_flag) ;

    if (night_only)
	FITS_delete_file(header, &status);
    else
	FITS_close_file(header, &status);

    cf_verbose(3, "number of pseudo photon events = %d ",nevents) ;

    /* Call routines to remove motions */
    if (grating_motion) {
      cf_verbose(4,"Correcting for grating motions") ;
       cf_grating_motion(bpmfits, nevents, time, x, y, channel,  
			 nseconds, timeline, tsunrise) ; }
    if (fpa_pos) {
      cf_verbose(4,"Correcting for fpa motions ") ;
      cf_fpa_position(bpmfits, nevents, x, channel) ;}
    if (mirror_motion) {
      cf_verbose(4,"Correcting for mirror motions") ;
      cf_mirror_motion(bpmfits, nevents, time, x, y, channel,
		       nseconds, timeline, tsunset) ; }
    if (jitter) {
      cf_verbose(4,"Correcting for satellite jitter ") ;
      cf_satellite_jitter(bpmfits, nevents, time, x, y, channel, nseconds, 
	 timeline, statflag) ;
      }

     /* correct for the Doppler motions */
     if (doppler_motion)
        cf_correct_doppler_motions(bpmfits, nevents, time, x, channel, 
				   nseconds, timeline, velocity) ;

    /* combine the pothole centroid information and generate a pothole map by
	including the pothole sizes and shapes */
    npixels = cf_combine_pothole_data(bpmfits, nevents, time, weight,
	x, rx, y, ry, channel, timeflag, &xpix, &ypix, &chan_pix, &wt_pix) ;

    /* generate a blank wavelength array */
    lambda = (float *) cf_calloc( npixels, sizeof(float) ) ;

    /* correct for astigmatism and assign wavelengths  */
    if (astig) cf_astigmatism(bpmfits, npixels, xpix, ypix, chan_pix) ; 
    cf_dispersion(bpmfits, npixels, xpix, chan_pix, lambda);
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
 
    /* create a table and output the results to the file */
    /* first set tform values */
    cf_verbose(3, "Writing pseudo-pixel information to output file.");
    sprintf(fmt_float, "%1ldE", npixels) ;
    sprintf(fmt_byte, "%1ldB", npixels) ;
    tform[0]=fmt_float ;
    tform[1]=fmt_float ;
    tform[2]=fmt_byte ;
    tform[3]=fmt_float ;
    tform[4]=fmt_float ;

    /* put the data into the table */
    FITS_create_tbl(bpmfits, BINARY_TBL, nrows, tfields, ttype, tform,
		     tunit, extname, &status) ;
    cf_verbose(3, "    Photon list created.");
    FITS_write_col(bpmfits, TFLOAT, 1, 1L, 1L, npixels, xpix, &status) ;
    cf_verbose(3, "    X array has been written.");
    FITS_write_col(bpmfits, TFLOAT, 2, 1L, 1L, npixels, ypix, &status) ;
    cf_verbose(3, "    Y array has been written.");
    FITS_write_col(bpmfits, TBYTE, 3, 1L, 1L, npixels, chan_pix, &status) ;
    cf_verbose(3, "    Channel array has been written.");
    FITS_write_col(bpmfits, TFLOAT, 4, 1L, 1L, npixels, wt_pix, &status) ;
    cf_verbose(3, "    Weights array has been written.");
    FITS_write_col(bpmfits, TFLOAT, 5, 1L, 1L, npixels, lambda, &status) ; 
    cf_verbose(3, "    Lambda array has been written.");

    FITS_close_file(bpmfits, &status);

    free(statflag);
    free(timeline);
    free(timeflag);
    free(rx) ;
    free(ry) ;
    free(y);
    free(x);
    free(time);
    free(xpix) ;
    free(ypix) ;
    free(chan_pix) ;
    free(wt_pix) ;
    free(lambda) ;

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
