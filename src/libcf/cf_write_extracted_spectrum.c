/**************************************************************************
 *              Johns Hopkins University
 *              Center for Astrophysical Sciences
 *              FUSE
 *************************************************************************
 *
 * synopsis:    cf_write_extracted_spectrum(infile, aperture, nout, wave_out,
 *                    flux_out, var_out, counts_out, weights_out, bkgd_out,
 *                    bpix_out, outrootname)
 *
 * Description: Write the results of the data reduction to a file.
 *              Information from the header is used to generate the name
 *              of the output file.  Aperture coding is as follows:
 *
 *              1-8 -> LiF[HIRS,MDRS,LWRS,PINH], SiC[HIRS,MDRS,LWRS,PINH]
 *
 * Arguments:   fitsfile  infile   :   pointer to Intermediate Data File
 *              int       aperture :   Target aperture (values are 1 to 8)
 *		int  valid_spectrum:   True if spectrum contains photons.
 *              long      nout     :   number of tabulated points
 *              float     wave_out :   Output wavelength array
 *              float     flux_out :   Output flux array
 *              float     var_out  :   Output variance array
 *              long    counts_out :   Output counts array (in counts)
 *              float  weights_out :   Output weights array (in counts)
 *              float     bkgd_out :   Output background array (in counts)
 *              short     bpix_out :   Output bad pixel spectrum
 *              char    outrootname:   Optional output file name
 *
 * Returns:	0 if successful
 *                                
 * HISTORY:     02/10/03   1.1   RDR    started work
 *              02/28/03   1.2   peb    Added header files to file and tidied
 *                                      code using lint.
 *		03/12/03   1.3   wvd    Output ERROR = sqrt(var_out)
 *		03/31/03   1.4   wvd	Update keywords FILENAME, FILETYPE,
 *                                      and APER_ACT.
 *		05/22/03   1.5   wvd	Direct cf_error_init to stderr.
 *              05/29/03   1.6   rdr    Added bad pixel column
 *              06/03/03   1.7   rdr    Change naming convention to allow
 *					hist data
 *              07/17/03   1.9   wvd    Test var_out before taking sqrt.
 *					Change calfusettag.h to calfuse.h
 *              10/01/03   1.10  wvd    Modify verbose levels.
 *              11/10/03   1.11  wvd    Correct numbering scheme for HIRS
 *					and MDRS apertures.
 *              02/13/04   1.12  wvd    Correct error in COMMENT lines.
 *              03/03/04   1.13  rdr    Populate archive search keywords
 *              03/05/04   1.14  wvd    Change name of POTHOLE column 
 *					to QUALITY.
 *              03/24/04   1.15  wvd    If valid_spectrum=FALSE, write
 *					warning message to file header.
 *              04/06/04   1.16  bjg    Include string.h and ctype.h
 *              04/09/04   1.17  wvd    If valid_spectrum=FALSE, set
 *					EXPTIME = 0 in file header.
 *              04/27/04   1.18  wvd    Move conversion from variance
 *					to sigma from this routine to
 *					cf_optimal_extraction.
 *              06/11/04   1.19  peb    Added the -r option to override the
 *                                      default rootname
 *		02/01/05   1.20  wvd	Revert to standard FITS binary
 *					table format for output file.
 *		04/11/05   1.21  wvd	In HIST mode, don't write SIA
 *					table to output spectral file.
 *		06/03/05   1.22  wvd	Return status from cf_copy_hist_header.
 *		06/03/05   1.23  wvd	Change cf_copy_hist_header to type void.
 *		08/11/05   1.24  wvd	Delete FITS_write_date from 
 *					cf_copy_hist_header.
 *
 **************************************************************************/

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "calfuse.h"

static char *ap_code[] = {"", "lif3", "lif2", "lif4", "lif1",
                              "sic3", "sic2", "sic4", "sic1"};

static void
cf_copy_hist_header(fitsfile *infits, fitsfile *outfits, int *status)
{
	char	card[FLEN_CARD];
	int	i, nkeys, naxis=0, nextend=1;
	long	naxes[2];

	naxes[0] = naxes[1] = 0;

	FITS_create_img(outfits, SHORT_IMG, naxis, naxes, status);
	FITS_delete_key(outfits, "COMMENT", status);
	FITS_delete_key(outfits, "COMMENT", status);
	FITS_write_key(outfits, TINT, "NEXTEND", &nextend,
		"Number of standard extensions", status);
	fits_get_hdrspace(infits, &nkeys, NULL, status);
	cf_if_fits_error(*status);
	for (i = 9; i <= nkeys; i++) {
	    FITS_read_record(infits, i, card, status);
	    FITS_write_record(outfits, card, status);
	}
}

int
cf_write_extracted_spectrum(fitsfile *infits, int aperture,
	int valid_spectrum, long nout, float *wave_out,
	float *flux_out, float *var_out, long *counts_out,
	float *weights_out, float *bkgd_out, short *bpix_out,
	char *outrootname)
{
    char CF_PRGM_ID[] = "cf_write_extracted_spectrum";
    char CF_VER_NUM[] = "1.24";

    fitsfile *outfits;
    int nextend=1, status=0, tfields=7, timeref;
    long i, frow=1, felem=1;
    float bw=0., center=0., wmin=1200., wmax=0. ;
    char comment[FLEN_CARD], datestr[FLEN_CARD];
    char rootname[FLEN_VALUE]={'\0'}, instmode[FLEN_VALUE] ;
    char det[FLEN_VALUE]={'\0'}, outfile[FLEN_VALUE]={'\0'};
    char extname[]="SPECTRUM";
    char *aperact[] = {"HIRS_LIF", "MDRS_LIF", "LWRS_LIF", "PINH_LIF",
                      "HIRS_SIC", "MDRS_SIC", "LWRS_SIC", "PINH_SIC"};
    char *ttype[]={"WAVE","FLUX","ERROR","COUNTS","WEIGHTS","BKGD","QUALITY"};
    char *tform[]={"1E", "1E", "1E", "1J", "1E", "1E", "1I"};
    char *tunit[]={"ANGSTROMS","ERG CM^-2 S^-1 ANG^-1","ERG CM^-2 S^-1 ANG^-1",
		"COUNTS", "COUNTS", "COUNTS","UNITLESS"};

    /* Initialize error checking */ 
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a time stamp into the log */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /*
     *  Read the observation ID from the input file and create an
     *  output file name.
    */
    FITS_movabs_hdu(infits, 1, NULL, &status);
    FITS_read_key(infits, TSTRING, "ROOTNAME", rootname, NULL, &status);
    FITS_read_key(infits, TSTRING, "DETECTOR", det, NULL, &status);
    FITS_read_key(infits, TSTRING, "INSTMODE", instmode, NULL, &status);

    /* Override the default rootname if the -r option is given. */
    if (outrootname != NULL)
        strcpy(rootname, outrootname);

    /* Convert detector name to lower case */
    det[1] = tolower(det[1]);
    
    if (!strncmp(instmode,"T",1))
        sprintf(outfile, "%11s%2s%4sttagfcal.fit",
                rootname, det, ap_code[aperture]);
    else if (!strncmp(instmode,"H",1))
        sprintf(outfile, "%11s%2s%4shistfcal.fit",
                rootname, det, ap_code[aperture]);
    else {
      cf_if_error("Cannot create name for the output files.") ;
      return 1;
    }
    cf_verbose(3, "Output file name = %s", outfile);

    /* Open the output file and copy the header info from the input file */
    FITS_create_file(&outfits, outfile, &status);
    if (!strncmp(instmode,"T",1))
        FITS_copy_header(infits, outfits, &status);
    else
	cf_copy_hist_header(infits, outfits, &status);
    FITS_update_key(outfits, TINT, "NEXTEND", &nextend, NULL, &status);
    FITS_update_key(outfits, TSTRING, "FILENAME", outfile, NULL, &status);
    FITS_update_key(outfits, TSTRING, "FILETYPE",
                 "CALIBRATED EXTRACTED SPECTRUM", NULL, &status);
    FITS_update_key(outfits, TSTRING, "APER_ACT", aperact[aperture-1],
                 NULL, &status);

    /* If valid_spectrum = FALSE, write warning to file header. */
    if (valid_spectrum == FALSE) {
	float exptime = 0.;
	FITS_update_key(outfits, TFLOAT, "EXPTIME", &exptime, NULL, &status);
        FITS_write_comment(outfits, "  ", &status);
        FITS_write_comment(outfits,
		"No good data.  Output arrays set to zero.", &status);
        fits_get_system_time(datestr, &timeref, &status);
        sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
        FITS_write_comment(outfits, comment, &status);
        FITS_write_comment(outfits, "  ", &status);
    }

    /* Populate archive search keywords */
    for (i=0 ; i<nout ; i++) {
      if (wave_out[i] < wmin) wmin=wave_out[i] ;
      if (wave_out[i] > wmax) wmax=wave_out[i] ;
    }
       bw = wmax-wmin ;
       center = (wmax + wmin) / 2. ;
    FITS_update_key(outfits, TFLOAT, "BANDWID", &bw, NULL, &status);
    FITS_update_key(outfits, TFLOAT, "CENTRWV", &center, NULL, &status);
    FITS_update_key(outfits, TFLOAT, "WAVEMIN", &wmin, NULL, &status);
    FITS_update_key(outfits, TFLOAT, "WAVEMAX", &wmax, NULL, &status);

    /* Put the data into the table */
    FITS_create_tbl(outfits, BINARY_TBL, nout, tfields, ttype, tform,
		    tunit, extname, &status);
    FITS_write_col(outfits, TFLOAT, 1, frow, felem, nout, wave_out, &status);
    FITS_write_col(outfits, TFLOAT, 2, frow, felem, nout, flux_out, &status);
    FITS_write_col(outfits, TFLOAT, 3, frow, felem, nout, var_out, &status);
    FITS_write_col(outfits, TLONG,  4, frow, felem, nout, counts_out, &status);
    FITS_write_col(outfits, TFLOAT, 5, frow, felem, nout, weights_out, &status);
    FITS_write_col(outfits, TFLOAT, 6, frow, felem, nout, bkgd_out, &status);
    FITS_write_col(outfits, TSHORT, 7, frow, felem, nout, bpix_out, &status);
    
    /* Write some comments into the header */
    FITS_write_comment(outfits, " ", &status);
    FITS_write_comment(outfits, "RAW COUNTS = COUNTS ",
		       &status);
    FITS_write_comment(outfits, "TARGET COUNTS = WEIGHTS - BKGD ",
		       &status);
    FITS_write_comment(outfits, "TARGET FLUX = TARGET COUNTS * HC / LAMBDA / "
		       "AEFF / EXPTIME / WPC", &status);
    FITS_write_comment(outfits, " ", &status);

    FITS_close_file(outfits, &status);

    /* Enter a time stamp into the log */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");

    return (status);
}
