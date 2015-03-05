/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Description: Correct XFARF array for instrumental astigmatism.
 *
 * Returns:    0 upon successful completion.
 *
 * History:       12/12/02   1.1   jch  Begin work
 *                02/24/03   1.2   peb  Using calfusettag.h include file
 *                03/04/03   1.3   peb  Removed math.h - unnecessary
 *                03/11/03   1.4   wvd  Changed channel to type *char
 *                04/04/03   1.5   wvd  Test for overflow of astig array.
 *                05/20/03   1.6   rdr  Add call to cf_proc_check
 *                09/16/03   1.7   wvd  Remove calls to astig_target_aperture
 *					and astig_read_file.  Change sign of
 *					astigmatism correction to match
 *					the new astg**009.fit files.
 *					Make yoff a float.
 *                11/05/03   1.8   wvd  Change channel to unsigned char.
 *                04/07/07   1.9   wvd  Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_astig_farf";
static char CF_VER_NUM[] = "1.9";


/* 
 * Add astigmatism correction to XFARF for each photon within target
 * aperture.  Note that spectrum centroid moves with time.
 */
static int
farf_astig_shifts(fitsfile *header, long nevents, float *xfarf, 
	float *yfarf, unsigned char *channel, float *photon_time, 
	long nseconds, float *timeline_time, float *ycentl, float *ycents)
{
    char        astig_file[FLEN_VALUE];
    int   	status=0;
    int   	active_ap[2], ap, ast_cen, extno, overflow=0;
    long	i, j, k;
    long  	npix, nx, ny, ii, jj, kk;
    float 	yoff;
    float       *astig=NULL, *centroid=NULL;
    fitsfile    *astigfits;

    /* Read target apertures from the file header */
    (void) cf_source_aper (header, active_ap);

    /* Open the astigmatism file */
    FITS_read_key(header, TSTRING, "ASTG_CAL", astig_file, NULL, &status);
    FITS_open_file(&astigfits, cf_cal_file(astig_file), READONLY, &status);

    /*   Since we presently have no astig correction for extended
     *   sources, we apply correction only to target aperture.
     */
    for (i = 0; i < 2; i++) {

        ap = active_ap[i];
	if (i == 0) 
	    centroid = ycentl;
	else
	    centroid = ycents;

        /* Read astigmatism correction for this aperture */
        extno = ap+1;
        cf_verbose(3, "Reading extension %d of %s", extno, astig_file);
        FITS_movabs_hdu(astigfits, extno, NULL, &status);
        FITS_read_key(astigfits, TINT, "SLIT_CEN", &ast_cen, NULL, &status);
        FITS_read_key(astigfits, TLONG, "NAXIS1", &nx, NULL, &status);
        FITS_read_key(astigfits, TLONG, "NAXIS2", &ny, NULL, &status);

        npix = nx * ny;
        astig = (float *) cf_malloc(sizeof(float) * npix);
        FITS_read_img(astigfits, TFLOAT, 1L, npix, 0, astig, NULL, &status);

	/* Go through the photon list, find the astigmatism (delta(X)) 
         * corresponding to the (xfarf, yfarf) position and channel,
	 * and add it to xfarf[k].
	 */
	k = 0;
        for (j = 0; j < nevents; j++) {
            if (channel[j] == ap) {
                while ((photon_time[j] > timeline_time[k+1] - FRAME_TOLERANCE)
		    && (k < nseconds-1)) {
                    k++;
                }

        	/* Determine offset between spectrum and astig correction file. */
        	yoff = centroid[k] - ast_cen;

		/* Find the astigmatism correction appropriate 
		   to the photon's (xfarf, yfarf) position. */
                if ((ii = (xfarf[j] + 0.5)) >= nx) {overflow = 1; continue;};
                if ((jj = (yfarf[j] - yoff + 0.5)) >= ny) {overflow = 1; continue;}
                if ((kk = jj * nx + ii) >= npix) {overflow = 1; continue;}
                xfarf[j] += astig[kk];
            }
        }

        /* If overflow flag is set, there is a problem with the ASTG_CAL file. */
	if (overflow)
    	    cf_if_warning("Overflow of ASTG_CAL file in aperture %d", ap);

        /* Space for astig array is allocated in each loop. */
        free(astig);
    }

    FITS_close_file(astigfits, &status);
    return status;
}

int cf_astig_farf(fitsfile *header, long nevents, float *xfarf, float *yfarf, 
		  unsigned char *channel, float *photon_time, long nseconds, 
		  float *timeline_time, float *ycentl, float *ycents)
{
    int 	errflg=0, status=0;

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* Check that astigmatism correction is appropriate for input image */
    if (astig_check_input_image(header) == 1) {

        /* Read astigmatism corrections for each aperture, add to XFARF */
        farf_astig_shifts(header, nevents, xfarf, yfarf, channel,
		photon_time, nseconds, timeline_time, ycentl, ycents);
    }

    /* Update processing flags. */
    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");
    return (status);
}
