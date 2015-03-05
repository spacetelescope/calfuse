/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE 
 *****************************************************************************
 *
 * Synopsis:    cf_identify_channel(fitsfile *header, long nevents,
 *              float *xfarf, float *yfarf, unsigned char *chan, 
 *		unsigned char *locflags, int pad, int final_call)
 *
 * Description: Assigns an aperture to each photon in the IDF.
 *
 *		Change of logic in v1.11: This routine may be called
 *		multiple times.  All photons are examined and assigned to
 *		a channel unless final_call = TRUE.  On the final call, only
 *		photons with channel > 0 are considered.  Because only photons
 *		assigned to a channel are motion corrected, the spectrum can
 *		move on top of previously-rejected photons.  We don't want
 *		to include them in the final spectrum, so leave their channel
 *		values unchanged.
 *
 *		Change of logic in v1.16: We've decided to use extended-
 *		source apertures for the non-target channels.  These
 *		apertures can overlap, so we assign disputed photons to
 *		their most likely source.  In decreasing order, these are
 *		the target aperture, the larger non-target aperture, and
 *		the smaller non-target aperture.
 *		
 * Arguments:   fitsfile  *header       Pointer to the location of the FITS
 *                                      header of the Intermediate Data File
 *              long      nevents       Number of photons in the file
 *              float     *xfarf        X positions for each photon event
 *              float     *yfarf        Y positions for each photon event
 *     unsigned char      *chan         Pointer to pointer to the array 
 *                                      containing the channel information
 *     unsigned char      *locflags     Location flags for each photon 
 *              int       pad           Switch saying whether to expand the
 *                                      width of the tabulated extraction
 *                                      window (yes if > 0). The number refers
 *                                      to the amount of padding (e.g. pad=5
 *                                      will add 5 pixels to each side of the
 *                                      extraction window.
 *		int	  final_call	This subroutine may be called
 *					many times.  The last time, set
 *					final_call = TRUE.
 *
 * Return:      0 on success
 *
 * History:     08/09/02    1.0  RDR    Begin work 
 *              09/04/02    1.1  rdr    Change program so that it just fills
 *                                      the channel array rather than creating
 *                                      and filling it.  Decreased the size of
 *                                      the mask.  Use the value of pad to
 *                                      specify degree of window padding.
 *                                      Check for extended source.
 *              12/11/02    1.2  RDR    Adjust aperture limits by measured
 *                                      y centroid (if available)
 *              02/11/03    1.3  rdr    changed test for ycent values
 *              02/12/03    1.6  wvd    Don't crash if YCENT# is undefined
 *              03/25/03    1.7  rdr    Ignore pinhole aperture 
 *		03/31/03    1.10 wvd    Ignore regions where YMAX = YMIN + 1.
 *					They are not valid channel regions.
 *		04/07/03    1.11 wvd    Shift channel boundaries to match FPA
 *					positions.  Implement use of cf_verbose.
 *					Update CHID_COR if final_call = TRUE.
 *              05/20/03    1.12 rdr    Added call to cf_proc_check
 *              07/23/03    1.14 rdr    Correct procedure which reads CHID
 *                                      calibration file
 *		08/21/03    1.15 wvd    Change channel to unsigned char.
 *		08/25/03    1.16 wvd    Use extended windows for non-target
 *					apertures.  Replace aperture mask
 *					with a grid of aperture limits.
 *					Use cf_nint() to round floats.
 *              08/25/03    1.17 wvd    Change coltype from string to int
 *					in cf_read_col.
 *              09/08/03    1.18 wvd    Default value is channel[i] = 0.
 *					Identify low and high values of x 
 *					by setting y limits to -NYMAX.
 *              09/18/03    1.19 wvd    Consider only photons falling on the
 *					active area of the detector.
 *              11/12/03    1.20 wvd    Don't bother with non-target channels
 *					in HIST mode.
 *              06/14/04    1.21 wvd    Don't bother shifting extraction
 *					windows to correct for FPA position.
 *					If known, shift is not needed.
 *              06/17/04    1.22 wvd    Because the extraction windows don't
 *					quite match the WGTS arrays, we pad
 *					windows by XPAD pixels in X.
 *              07/21/04    1.23 wvd    Change i from long to int.
 *              03/15/05    1.24 wvd    Call cf_extraction_limits to get limits
 *					of extraction windows, rather than
 *					reading CHID_CAL file directly.
 *					Do not pad windows in X.
 *              03/16/05    1.25 wvd    Pass srctype to cf_extraction_limits.
 *
 ****************************************************************************/

#include <string.h>
#include "calfuse.h"

int
cf_identify_channel(fitsfile* header, long npts, float* x, float* y,
	unsigned char* channel, unsigned char* locflags, int pad,
	int final_call)
{
    char CF_PRGM_ID[] = "cf_identify_channel";
    char CF_VER_NUM[] = "1.25";

    char  instmode[FLEN_VALUE], zero=0;
    int   errflg=0, extended, n_chan=6, status=0;
    int   active_ap[2], chan, srctype[8];
    int   i;	/* Always steps through targ_ap */
    int   *targ_ap,
	  lwrs[]={3,7,2,6,1,5},	/* LWRS, MDRS, HIRS (LiF & SiC) */
	  mdrs[]={2,6,3,7,1,5},	/* MDRS, LWRS, HIRS (LiF & SiC) */
	  hirs[]={1,5,3,7,2,6};	/* HIRS, LWRS, MDRS (LiF & SiC) */
    long  j,	/* Always steps through extraction window arrays */
	  k;	/* Always steps through photon list */
    short xint, yint, xmin, xmax, *ymin[6], *ymax[6];

    /* Initialize error checking */ 
     cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a time stamp into the log */
     cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Check whether routine is appropriate for this data file. */
     if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* If data were taken in HIST mode, ignore non-target channels. */
     FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
     if (!strncmp(instmode, "HIST", 4)) n_chan = 2;

    /* Determine source type and active apertures */
     extended = cf_source_aper(header, active_ap);

    /* Set the order in which we will assign photons to channels. */
     switch (active_ap[0]) {
	case 1: targ_ap = hirs;
		break;
	case 2: targ_ap = mdrs;
		break;
	default: targ_ap = lwrs;
		break;
     }

    /* For target apertures, srctype reflects target.
	Use extended apertures for the rest. */
     for (i = 0; i < 2;  i++) srctype[i] = extended;
     for ( ; i < n_chan; i++) srctype[i] = 1;

    /* While we're at it, here's one more bit of i/o. */
     if (pad > 0)
        cf_verbose(3, "Padding extraction windows by %d Y pixels", pad);
 
    /* 
     * Loop through all channels, read extraction limits,
     * and apply any requested padding.
     */
     for (i = 0; i < n_chan; i++) {
	chan = targ_ap[i];

       /* Read the extraction limits for this aperture. */
	(void) cf_extraction_limits(header, chan, srctype[i],
	    (ymin+i), (ymax+i), &xmin, &xmax);
 
	/* Add any padding requested by the user.  */
	if (pad > 0) {
	    for (j = xmin; j <= xmax; j++) {
		ymin[i][j] -= pad;
		ymax[i][j] += pad;
	    }
	}
    }

    /* Loop through photon list.  Assign a channel ID to each. */
    cf_verbose(3, "Entering loop to assign channel ID's.");
    for (k = 0; k < npts; k++) {
	if (final_call && channel[k] == zero) continue;
	channel[k] = zero;
	if (!(locflags[k] & LOCATION_SHLD)) {
	    xint = (short) cf_nint(x[k]);
	    yint = (short) cf_nint(y[k]);
	    if (xint < xmin || xint > xmax) continue;
	    for (i = 0; i < n_chan; i++) {
	        if (yint >= ymin[i][xint] && yint <= ymax[i][xint]) {
		    channel[k] = (char) targ_ap[i];
		    break;
		}
	    }
	}
    } 

  /* If final_call = TRUE, update IDF header keyword for this subroutine. */
     if (final_call) 
	cf_proc_update(header, CF_PRGM_ID, "COMPLETE");

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");

    return status;
}
