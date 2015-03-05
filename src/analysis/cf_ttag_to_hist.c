/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:   cf_ttag_to_hist input_file output_file
 *
 * Description: Creates a FITS 2D image from a time-tagged data file (FITS 
 *              binary table of photon address data).  The 2D SIA table is
 *              written to the first extension (HDU=1) and the other four
 *              extensions contain the two apertures and stim-lamp pulses.
 *
 * 
 * Arguments:   input_file      Input raw time-tagged FITS file name
 *              output_file     Output raw histogram FITS file name
 *
 * Returns:     none
 *
 * History:     03/28/00 1.01   peb     Begin work
 *              12/21/00 1.03   peb     Fixed type mismatched in cfitsio
 *                                      calls
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

#define  STIMX 2048

static char CF_PRGM_ID[] = "cf_ttag_to_hist";
static char CF_VER_NUM[] = "1.03";

int main(int argc, char *argv[])
{
    fitsfile *infits, *outfits;
    char   buffer[FLEN_CARD], *siabuf;
    int    status = 0, intnull, anynull, hdutype, fcol, begkey=5, endkey;
    int    j, k1, kx1, ky1, k2, kx2, ky2, k3, kx3, ky3, k4, kx4, ky4;
    int    naxis, naxis1, naxis2, naxis3, naxis4;
    int    *posx, *posy, npixel, fpixel=1, frow, felem, nrow;
    int    begx1, begy1, begx2, begy2, begx3, begy3, begx4, begy4;
    int    rcnt1, rcnt2, rcnt3, rcnt4, min1, min2, min3, min4;
    int    max1, max2, max3, max4, binx, biny;
    int    keynum, nextend=4, extver=0, bzero=32768, bscale=1;
    short  *outbuf1, *outbuf2, *outbuf3, *outbuf4;
    long   naxes[2], naxes1[2], naxes2[2], naxes3[2], naxes4[2];
    float  fzero=0.0;
    double ra_targ, dec_targ, pa_aper, equinox;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    /*
     *  Open input ttag data file and create output file.
     */
    FITS_open_file(&infits, argv[1], READONLY, &status);
    FITS_create_file(&outfits, argv[2], &status);
    /*
     *  Create primary image header
     */
    naxis = 2; naxes[0] = 8; naxes[1] = 64;
    FITS_create_img(outfits, BYTE_IMG, naxis, naxes, &status);
    /*
     *  Copy header from primary HDU of the input file to the output file,
     *  but leave the primary HDU data array empty.
     */
    FITS_read_key(infits, TINT, "SPECBINX", &binx, buffer, &status);
    FITS_read_key(infits, TINT, "SPECBINY", &biny, buffer, &status);
    FITS_get_hdrpos(infits, &endkey, &keynum, &status);
    for (j = begkey; j <= endkey; j++) {
	FITS_read_record(infits, j, buffer, &status);
	FITS_write_record(outfits, buffer, &status);
    }
    FITS_update_key(outfits, TINT, "NEXTEND", &nextend, NULL, &status);
    FITS_update_key(outfits, TSTRING, "FILETYPE", "RAW HISTOGRAM",
		    NULL, &status);
    FITS_update_key(outfits, TSTRING, "INSTMODE", "HIST", NULL, &status);
    /*
     *  Populate primary data or SIA table.
     */
    npixel = naxes[0]*naxes[1];
    siabuf = cf_malloc(npixel);
    memset(siabuf, 1, npixel);
    FITS_write_img(outfits, TBYTE, fpixel, npixel, siabuf, &status);
    free(siabuf);
    /*
     *  Move to the 1st HDU of the input file and read the binary table.
     */
    FITS_movabs_hdu(infits, 2, &hdutype, &status);
    if(hdutype != BINARY_TBL) {
	cf_if_error("Error: This HDU is not a binary table!");
    }
    FITS_read_key(infits, TINT,    "NAXIS2",  &nrow, buffer, &status);
    FITS_read_key(infits, TDOUBLE, "RA_TARG", &ra_targ, buffer, &status);
    FITS_read_key(infits, TDOUBLE, "DEC_TARG", &dec_targ, buffer, &status);
    FITS_read_key(infits, TDOUBLE, "PA_APER", &pa_aper, buffer, &status);
    FITS_read_key(infits, TDOUBLE, "EQUINOX", &equinox, buffer, &status);
    /*
     *  Allocate memory for X, Y photon event arrays
     */
    posx = cf_malloc(sizeof(int) * nrow);
    posy = cf_malloc(sizeof(int) * nrow);
    /*
     *  Read X and Y from table
     */
    frow = felem = 1; intnull = anynull = 0;
    
    FITS_get_colnum(infits, TRUE, "X", &fcol, &status);
    FITS_read_col(infits, TINT, fcol, frow, felem, nrow, &intnull,
		  posx, &anynull, &status);
    FITS_get_colnum(infits, TRUE, "Y", &fcol, &status);
    FITS_read_col(infits, TINT, fcol, frow, felem, nrow, &intnull,
		  posy, &anynull, &status);
    /*
     *  Create buffers (outbufn) to hold binned images, and
     *  initialize them.
     */
    naxis1 = 2; naxes1[0] = NXMAX; naxes1[1] = 20;
    begx1 = 0; begy1 = 48; rcnt1 = 56; min1 = max1 = 0;
    outbuf1 = cf_calloc(naxes1[0]*naxes1[1], sizeof(short));

    naxis2 = 2; naxes2[0] = NXMAX; naxes2[1] = 20;
    begx2 = 0; begy2 = 76; rcnt2 = 57; min2 = max2 = 0;
    outbuf2 = cf_calloc(naxes2[0]*naxes2[1], sizeof(short));

    naxis3 = 2; naxes3[0] = STIMX; naxes3[1] = 2;
    begx3 = 0; begy3 = 76; rcnt3 = 1; min3 = max3 = 0;
    outbuf3 = cf_calloc(naxes3[0]*naxes3[1], sizeof(short));

    naxis4 = 2; naxes4[0] = STIMX; naxes4[1] = 2;
    begx4 = 14336; begy4 = 78; rcnt4 = 1; min4 = max4 = 0;
    outbuf4 = cf_calloc(naxes4[0]*naxes4[1], sizeof(short));
    
    for(j = 0; j < nrow; j++) {
	kx1 = posx[j]/binx-begx1;
	ky1 = posy[j]/biny-begy1;
	if (kx1 >= 0 && kx1 < naxes1[0] && ky1 >= 0 && ky1 < naxes1[1]) {
	    k1 = naxes1[0]*ky1 + kx1;
	    outbuf1[k1] += 1;
	    if (outbuf1[k1] < min1)
		min1 = outbuf1[k1];
	    if (outbuf1[k1] > max1)
		max1 = outbuf1[k1];
	}
	kx2 = posx[j]/binx-begx2;
	ky2 = posy[j]/biny-begy2;
	if (kx2 >= 0 && kx2 < naxes2[0] && ky2 >= 0 && ky2 < naxes2[1]) {
	    k2 = naxes2[0]*ky2 + kx2;
	    outbuf2[k2] += 1;
	    if (outbuf2[k2] < min2)
		min2 = outbuf2[k2];
	    if (outbuf2[k2] > max2)
		max2 = outbuf2[k2];
	}
	kx3 = posx[j]/binx-begx3;
	ky3 = posy[j]/biny-begy3;
	if (kx3 >= 0 && kx3 < naxes3[0] && ky3 >= 0 && ky3 < naxes3[1]) {
	    k3 = naxes3[0]*ky3 + kx3;
	    outbuf3[k3] += 1;
	    if (outbuf3[k3] < min3)
		min3 = outbuf3[k3];
	    if (outbuf3[k3] > max3)
		max3 = outbuf3[k3];
	}
	kx4 = posx[j]/binx-begx4;
	ky4 = posy[j]/biny-begy4;
	if (kx4 >= 0 && kx4 < naxes4[0] && ky4 >= 0 && ky4 < naxes4[1]) {
	    k4 = naxes4[0]*ky4 + kx4;
	    outbuf4[k4] += 1;
	    if (outbuf4[k4] < min4)
		min4 = outbuf4[k4];
	    if (outbuf4[k4] > max4)
		max4 = outbuf4[k4];
	}
    }
    /*
     *  Extension 1
     */
    FITS_create_img(outfits, SHORT_IMG, naxis1, naxes1, &status);
    extver = 1;
    
    FITS_write_key(outfits, TSTRING, "EXTNAME", "HISTOGRAM",
		   "name of this extension", &status);
    FITS_write_key(outfits, TINT, "EXTVER" , &extver,
		   "extension version number", &status);
    FITS_write_key(outfits, TINT, "BZERO"  , &bzero,
		   "image brightness offset", &status);
    FITS_write_key(outfits, TINT, "BSCALE" , &bscale,
		   "image brightness scale", &status);
    
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits,
		       "    World Coordinate System and Related Parameters",
		       &status);
    FITS_write_comment(outfits, "        ", &status);
    
    FITS_write_key(outfits, TFLOAT,  "CRPIX1",  &fzero,
		   "x-coordinate of reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRPIX2",  &fzero, 
		   "y-coordinate of reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRVAL1",  &fzero, 
		   "first axis value at reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRVAL2",  &fzero,
		   "second axis value at reference pixel", &status);
    FITS_write_key(outfits, TSTRING, "CTYPE1", "LAMBDA", 
		   "the coordinate type for the first axis", &status);
    FITS_write_key(outfits, TSTRING, "CYTPE2", "ANGLE", 
		   "the coordinate type for the second axis", &status);
    FITS_write_key(outfits, TFLOAT,  "CD1_1",   &fzero,
		   "partial of first axis coordinate w.r.t. x", &status);
    FITS_write_key(outfits, TFLOAT,  "CD1_2",   &fzero,
		   "partial of first axis coordinate w.r.t. y", &status);
    FITS_write_key(outfits, TFLOAT,  "CD2_1",   &fzero,
		   "partial of second axis coordinate w.r.t. x", &status);
    FITS_write_key(outfits, TFLOAT,  "CD2_2",   &fzero,
		   "partial of second axis coordinate w.r.t. y", &status);
    FITS_write_key(outfits, TDOUBLE,  "RA_TARG", &ra_targ,
		   "RA of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "DEC_TARG", &dec_targ,
		   "Declination of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "PA_APER", &pa_aper,
		   "Position Angle of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "EQUINOX", &equinox,
		   "equinox of celestial coord. system", &status);
    
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits, "    HISTOGRAM REGION PARAMETERS", &status);
    FITS_write_comment(outfits, "        ", &status);
    
    FITS_write_key(outfits, TINT, "XORIGIN", &begx1,
		   "offset of this region in the x dimension", &status);
    FITS_write_key(outfits, TINT, "YORIGIN", &begy1,
		   "offset of this region in the y dimension", &status);
    FITS_write_key(outfits, TINT, "RCOUNT",  &rcnt1,
		   "count of rectangles in this region", &status);
    FITS_write_key(outfits, TINT, "MINVAL",  &min1,
		   "minimum value within rectangles", &status);
    FITS_write_key(outfits, TINT, "MAXVAL",  &max1,
		   "maximum value within rectangles", &status);

    FITS_write_img(outfits, TSHORT, fpixel, naxes1[0]*naxes1[1],
		   outbuf1, &status);
    /*
     *  Extension 2
     */
    FITS_create_img(outfits, SHORT_IMG, naxis2, naxes2, &status);
    extver = 2;
    
    FITS_write_key(outfits, TSTRING, "EXTNAME", "HISTOGRAM",
		   "name of this extension", &status);
    FITS_write_key(outfits, TINT, "EXTVER" , &extver,
		   "extension version number", &status);
    FITS_write_key(outfits, TINT, "BZERO"  , &bzero,
		   "image brightness offset", &status);
    FITS_write_key(outfits, TINT, "BSCALE" , &bscale,
		   "image brightness scale", &status);
    
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits,
		       "    World Coordinate System and Related Parameters",
		       &status);
    FITS_write_comment(outfits, "        ", &status);
    
    FITS_write_key(outfits, TFLOAT,  "CRPIX1",  &fzero,
		   "x-coordinate of reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRPIX2",  &fzero, 
		   "y-coordinate of reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRVAL1",  &fzero, 
		   "first axis value at reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRVAL2",  &fzero,
		   "second axis value at reference pixel", &status);
    FITS_write_key(outfits, TSTRING, "CTYPE1", "LAMBDA", 
		   "the coordinate type for the first axis", &status);
    FITS_write_key(outfits, TSTRING, "CYTPE2", "ANGLE", 
		   "the coordinate type for the second axis", &status);
    FITS_write_key(outfits, TFLOAT,  "CD1_1",   &fzero,
		   "partial of first axis coordinate w.r.t. x", &status);
    FITS_write_key(outfits, TFLOAT,  "CD1_2",   &fzero,
		   "partial of first axis coordinate w.r.t. y", &status);
    FITS_write_key(outfits, TFLOAT,  "CD2_1",   &fzero,
		   "partial of second axis coordinate w.r.t. x", &status);
    FITS_write_key(outfits, TFLOAT,  "CD2_2",   &fzero,
		   "partial of second axis coordinate w.r.t. y", &status);
    FITS_write_key(outfits, TDOUBLE,  "RA_TARG", &ra_targ,
		   "RA of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "DEC_TARG", &dec_targ,
		   "Declination of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "PA_APER", &pa_aper,
		   "Position Angle of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "EQUINOX", &equinox,
		   "equinox of celestial coord. system", &status);
    
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits, "    HISTOGRAM REGION PARAMETERS", &status);
    FITS_write_comment(outfits, "        ", &status);
    
    FITS_write_key(outfits, TINT, "XORIGIN", &begx2,
		   "offset of this region in the x dimension", &status);
    FITS_write_key(outfits, TINT, "YORIGIN", &begy2,
		   "offset of this region in the y dimension", &status);
    FITS_write_key(outfits, TINT, "RCOUNT",  &rcnt2,
		   "count of rectangles in this region", &status);
    FITS_write_key(outfits, TINT, "MINVAL",  &min2,
		   "minimum value within rectangles", &status);
    FITS_write_key(outfits, TINT, "MAXVAL",  &max2,
		   "maximum value within rectangles", &status);

    FITS_write_img(outfits, TSHORT, fpixel, naxes2[0]*naxes2[1],
		   outbuf2, &status);
    /*
     *  Extension 3
     */
    FITS_create_img(outfits, SHORT_IMG, naxis3, naxes3, &status);
    extver = 3;
    
    FITS_write_key(outfits, TSTRING, "EXTNAME", "HISTOGRAM",
		   "name of this extension", &status);
    FITS_write_key(outfits, TINT, "EXTVER" , &extver,
		   "extension version number", &status);
    FITS_write_key(outfits, TINT, "BZERO"  , &bzero,
		   "image brightness offset", &status);
    FITS_write_key(outfits, TINT, "BSCALE" , &bscale,
		   "image brightness scale", &status);
    
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits,
		       "    World Coordinate System and Related Parameters",
		       &status);
    FITS_write_comment(outfits, "        ", &status);
    
    FITS_write_key(outfits, TFLOAT,  "CRPIX1",  &fzero,
		   "x-coordinate of reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRPIX2",  &fzero, 
		   "y-coordinate of reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRVAL1",  &fzero, 
		   "first axis value at reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRVAL2",  &fzero,
		   "second axis value at reference pixel", &status);
    FITS_write_key(outfits, TSTRING, "CTYPE1", "LAMBDA", 
		   "the coordinate type for the first axis", &status);
    FITS_write_key(outfits, TSTRING, "CYTPE2", "ANGLE", 
		   "the coordinate type for the second axis", &status);
    FITS_write_key(outfits, TFLOAT,  "CD1_1",   &fzero,
		   "partial of first axis coordinate w.r.t. x", &status);
    FITS_write_key(outfits, TFLOAT,  "CD1_2",   &fzero,
		   "partial of first axis coordinate w.r.t. y", &status);
    FITS_write_key(outfits, TFLOAT,  "CD2_1",   &fzero,
		   "partial of second axis coordinate w.r.t. x", &status);
    FITS_write_key(outfits, TFLOAT,  "CD2_2",   &fzero,
		   "partial of second axis coordinate w.r.t. y", &status);
    FITS_write_key(outfits, TDOUBLE,  "RA_TARG", &ra_targ,
		   "RA of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "DEC_TARG", &dec_targ,
		   "Declination of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "PA_APER", &pa_aper,
		   "Position Angle of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "EQUINOX", &equinox,
		   "equinox of celestial coord. system", &status);
    
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits, "    HISTOGRAM REGION PARAMETERS", &status);
    FITS_write_comment(outfits, "        ", &status);
    
    FITS_write_key(outfits, TINT, "XORIGIN", &begx3,
		   "offset of this region in the x dimension", &status);
    FITS_write_key(outfits, TINT, "YORIGIN", &begy3,
		   "offset of this region in the y dimension", &status);
    FITS_write_key(outfits, TINT, "RCOUNT",  &rcnt3,
		   "count of rectangles in this region", &status);
    FITS_write_key(outfits, TINT, "MINVAL",  &min3,
		   "minimum value within rectangles", &status);
    FITS_write_key(outfits, TINT, "MAXVAL",  &max3,
		   "maximum value within rectangles", &status);

    FITS_write_img(outfits, TSHORT, fpixel, naxes3[0]*naxes3[1],
		   outbuf3, &status);
    /*
     *  Extension 4
     */
    FITS_create_img(outfits, SHORT_IMG, naxis4, naxes4, &status);
    extver = 4;
    
    FITS_write_key(outfits, TSTRING, "EXTNAME", "HISTOGRAM",
		   "name of this extension", &status);
    FITS_write_key(outfits, TINT, "EXTVER" , &extver,
		   "extension version number", &status);
    FITS_write_key(outfits, TINT, "BZERO"  , &bzero,
		   "image brightness offset", &status);
    FITS_write_key(outfits, TINT, "BSCALE" , &bscale,
		   "image brightness scale", &status);
    
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits,
		       "    World Coordinate System and Related Parameters",
		       &status);
    FITS_write_comment(outfits, "        ", &status);
    
    FITS_write_key(outfits, TFLOAT,  "CRPIX1",  &fzero,
		   "x-coordinate of reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRPIX2",  &fzero, 
		   "y-coordinate of reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRVAL1",  &fzero, 
		   "first axis value at reference pixel", &status);
    FITS_write_key(outfits, TFLOAT,  "CRVAL2",  &fzero,
		   "second axis value at reference pixel", &status);
    FITS_write_key(outfits, TSTRING, "CTYPE1", "LAMBDA", 
		   "the coordinate type for the first axis", &status);
    FITS_write_key(outfits, TSTRING, "CYTPE2", "ANGLE", 
		   "the coordinate type for the second axis", &status);
    FITS_write_key(outfits, TFLOAT,  "CD1_1",   &fzero,
		   "partial of first axis coordinate w.r.t. x", &status);
    FITS_write_key(outfits, TFLOAT,  "CD1_2",   &fzero,
		   "partial of first axis coordinate w.r.t. y", &status);
    FITS_write_key(outfits, TFLOAT,  "CD2_1",   &fzero,
		   "partial of second axis coordinate w.r.t. x", &status);
    FITS_write_key(outfits, TFLOAT,  "CD2_2",   &fzero,
		   "partial of second axis coordinate w.r.t. y", &status);
    FITS_write_key(outfits, TDOUBLE,  "RA_TARG", &ra_targ,
		   "RA of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "DEC_TARG", &dec_targ,
		   "Declination of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "PA_APER", &pa_aper,
		   "Position Angle of reference aperture center", &status);
    FITS_write_key(outfits, TDOUBLE,  "EQUINOX", &equinox,
		   "equinox of celestial coord. system", &status);
    
    FITS_write_comment(outfits, "        ", &status);
    FITS_write_comment(outfits, "    HISTOGRAM REGION PARAMETERS", &status);
    FITS_write_comment(outfits, "        ", &status);
    
    FITS_write_key(outfits, TINT, "XORIGIN", &begx4,
		   "offset of this region in the x dimension", &status);
    FITS_write_key(outfits, TINT, "YORIGIN", &begy4,
		   "offset of this region in the y dimension", &status);
    FITS_write_key(outfits, TINT, "RCOUNT",  &rcnt4,
		   "count of rectangles in this region", &status);
    FITS_write_key(outfits, TINT, "MINVAL",  &min4,
		   "minimum value within rectangles", &status);
    FITS_write_key(outfits, TINT, "MAXVAL",  &max4,
		   "maximum value within rectangles", &status);

    FITS_write_img(outfits, TSHORT, fpixel, naxes4[0]*naxes4[1],
		   outbuf4, &status);
    /*
     *  Free memory and close the input and output files
     */
    free(outbuf1);
    free(outbuf2);
    free(outbuf3);
    free(outbuf4);
    free(posx);
    free(posy);
    
    FITS_close_file(outfits, &status);
    FITS_close_file(infits, &status);

    return 0;
}
