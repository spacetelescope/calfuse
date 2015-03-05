/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Usage:    cf_pack combined_exposure_list out_file_name.fit
 *
 * Description: Write the 8 segments (1alif, 1blif etc...) into a single fits
 *              file with 8 extensions (in addition to the primary array).
 *
 * History:     04/08/05        tc      v1.0	First release
 *		05/02/05	wvd	v1.1	Copy SPEC and WOFF keywords to
 *						output file.
 *		06/03/05	wvd	v1.2	Delete unused variables.
 *		06/29/05	wvd	v1.3	Set keyword COMB_COR = COMPLETE.
 *		09/30/05	wvd	v1.4	Update archive search keywords.
 *		03/23/06	wvd	v1.5	If a segment is missing, write
 *						an empty image extension.
 *		05/19/06	wvd	v1.6	Write an index to the file
 *						extensions to the primary HDU.
 *		05/24/06	wvd	v1.7	Remove SPEC* keywords from 
 *						primary HDU.
 *		05/29/06	wvd	v1.8	Use while construction to 
 *						delete WOFF* keywords.
 *		06/22/06	wvd	v1.9	Truncate APER_ACT from (for
 *						example) MDRS_LIF to MDRS.
 *		08/14/06	wvd	v1.10	Tinker with code to copy SPEC
 *						keywords to new file header.
 *		04/07/07	wvd	v1.11	Clean up compiler warnings.
 *		08/27/07	bot	v1.12	Changed TINT to TLONG l.163,214,
 *						and 215 ; added L to every number
 * 						occurence related to i and j ;
 * 						Changed i to k l.123 to 131 to
 * 						keep an int in this case
 *
 ******************************************************************************/
 
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

/* Calfuse variables */
static char CF_PRGM_ID[] = "cf_pack";
static char CF_VER_NUM[] = "1.12";

int main(int argc, char *argv[])
{
	/* Variables */
	float obstime, obstime_max, woffset;
	float w0, wpc;
	int wave_cnum, flux_cnum, error_cnum;
	int n_spec, n_best_spec;
	long i, j, num_exp, n_rows;
	int k;
	
	char *ttype[] = { "WAVE", "FLUX", "ERROR" };
	char *tform[] = { "1E", "1E", "1E" };
	char *tunit[] = { "ANGSTROMS", "ERG/CM2/S/A", "ERG/CM2/S/A" };
	char *extname[] = { "1ALIF", "1BLIF", "2BLIF", "2ALIF", "1ASIC", "1BSIC", "2BSIC", "2ASIC" };
	char spec_list[][80] = { "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL" };

	/* New values for archive search keywords */
	float bandwid = 285.;
	float centrwv = 1046.;
	float wavemin = 904.;
	float wavemax = 1189.;
	
	FILE *pt_file;
	char comment[FLEN_COMMENT], spec_lis[80], out_name[80], spec_name[80];
	char aper_act[80], detector[80], segment[80];
	char err_str[80], aper_key[8], combmeth_key[9];
	char keyname[FLEN_KEYWORD], keyvalue[FLEN_VALUE];
	fitsfile *pt_fits, *pt_fits_out;
	int status = 0, hdutype = 0;
	int nextend_key = 8;
	float obstime_key = 0;

	    
	/***********************************
	** Enter a timestamp into the log **
	***********************************/
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
	cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);


	/***********************************************************	
	** Check for proper command-line usage and read arguments **
	***********************************************************/
	if (argc != 3)
		cf_if_error("Usage: cf_pack combined_exposure_list out_file_name.fit\n");
	
	strcpy(spec_lis, argv[1]);
	strcpy(out_name, "!");
	strcat(out_name, argv[2]); /* Force overwrite */
	
	/********************	
	** Open input list **
	********************/
	if ((pt_file = fopen(spec_lis, "r")) == NULL)
	{
		sprintf(err_str, "Unable to open file %s\n", spec_lis);
		cf_if_error(err_str);
	}
	
	/* ***********************************************
	** Determine which spectrum has maximum OBSTIME **
	** Write input file list to spec_list array.	**
	*************************************************/
	n_best_spec = -1;
	n_spec = 0;
	obstime_max = -1;
    
	while (fscanf(pt_file, "%80s", spec_name) != EOF)
	{
		FITS_open_file(&pt_fits, spec_name, READONLY, &status);
		FITS_read_key(pt_fits, TFLOAT, "OBSTIME" , &obstime, NULL, &status);
		FITS_read_key(pt_fits, TSTRING, "DETECTOR" , &detector, NULL, &status);
		FITS_read_key(pt_fits, TSTRING, "APER_ACT" , &aper_act, NULL, &status);
		FITS_close_file(pt_fits, &status);

		sprintf(segment, "%2s%3s", detector, aper_act+5);
		for (k=0; k<8; k++)
			if (!strcmp(segment, extname[k])) {
				strcpy(spec_list[k], spec_name); 
				break;
			}
        	if (obstime_max < obstime) {
			n_best_spec = k;
			obstime_max = obstime;
		}
		n_spec++;
	}
	fclose(pt_file);
	if (n_spec == 0) cf_if_error("No valid segments provided");
	if (n_spec != 8) cf_if_warning("Combining fewer than 8 segments");

	/* ***************************************************
	** Copy header of reference spectrum to output file **
	** Update a few header keywords.	 	    **
	*****************************************************/

	FITS_open_file(&pt_fits, spec_list[n_best_spec], READONLY, &status);
	FITS_create_file(&pt_fits_out, out_name, &status);
	FITS_copy_header(pt_fits, pt_fits_out, &status);
	FITS_close_file(pt_fits, &status);

	FITS_update_key(pt_fits_out, TSTRING, "FILENAME", (out_name + 1), NULL, &status);
	FITS_update_key(pt_fits_out, TSTRING, "DETECTOR", "ALL", NULL, &status);
	FITS_update_key(pt_fits_out, TSTRING, "COMB_COR", "COMPLETE", NULL, &status);
	
	sprintf(aper_key, "%4s", aper_act);
	FITS_update_key(pt_fits_out, TSTRING, "APER_ACT", aper_key, NULL, &status);
	FITS_update_key(pt_fits_out, TINT, "NEXTEND", &nextend_key, NULL, &status);

	FITS_update_key(pt_fits_out, TFLOAT, "BANDWID", &bandwid, NULL, &status);
	FITS_update_key(pt_fits_out, TFLOAT, "CENTRWV", &centrwv, NULL, &status);
	FITS_update_key(pt_fits_out, TFLOAT, "WAVEMIN", &wavemin, NULL, &status);
	FITS_update_key(pt_fits_out, TFLOAT, "WAVEMAX", &wavemax, NULL, &status);

	FITS_read_key(pt_fits_out, TLONG, "NSPEC", &num_exp, NULL, &status);
	FITS_delete_key(pt_fits_out, "NSPEC", &status);
	for (i = 0L; i < num_exp; i++) {
	    sprintf(keyname, "SPEC%03ld", i+1L);
	    FITS_delete_key(pt_fits_out, keyname, &status);
            sprintf(keyname, "WOFF*%03ld", i+1L);
	    FITS_delete_key(pt_fits_out, keyname, &status);
	    fits_delete_key(pt_fits_out, keyname, &status);
	    if (status != 0) status = 0;
	}


	/* ***************************************************
	** Step through detector segments.  
	** Copy spectra into extensions of output file.
	*****************************************************/
	for (i = 0L; i < 8L; i++) {

	    if (!strncmp(spec_list[i], "NULL", 4)) {
		n_rows = 0L;
		FITS_create_tbl(pt_fits_out, BINARY_TBL, n_rows, 3, ttype, tform,
		tunit, extname[i], &status);
		FITS_write_comment(pt_fits_out, "Extension is empty.  No data available.", &status);
	    }
	    else {
		FITS_open_file(&pt_fits, spec_list[i], READONLY, &status);

		FITS_movabs_hdu(pt_fits, 2, &hdutype, &status);
        	FITS_get_colnum(pt_fits, TRUE, "WAVE", &wave_cnum, &status);
		FITS_get_colnum(pt_fits, TRUE, "FLUX", &flux_cnum, &status);
		FITS_get_colnum(pt_fits, TRUE, "ERROR", &error_cnum, &status);
		FITS_get_num_rows(pt_fits, &n_rows, &status);
		
		FITS_create_tbl(pt_fits_out, BINARY_TBL, n_rows, 3, ttype, tform,
			tunit, extname[i], &status);
		
		FITS_copy_col(pt_fits, pt_fits_out, wave_cnum, 1, FALSE, &status);
		FITS_copy_col(pt_fits, pt_fits_out, flux_cnum, 2, FALSE, &status);
		FITS_copy_col(pt_fits, pt_fits_out, error_cnum, 3, FALSE, &status);
		
		/* Copy various header keywords */
		FITS_read_key(pt_fits, TSTRING, "COMBMETH", combmeth_key, NULL, &status);
		FITS_update_key(pt_fits_out, TSTRING, "COMBMETH", combmeth_key, NULL, &status);

		FITS_movabs_hdu(pt_fits, 1, &hdutype, &status);
		FITS_read_key(pt_fits, TFLOAT, "W0", &w0, NULL, &status);
		FITS_update_key(pt_fits_out, TFLOAT, "W0", &w0, NULL, &status);
		FITS_read_key(pt_fits, TFLOAT, "WPC", &wpc, NULL, &status);
		FITS_update_key(pt_fits_out, TFLOAT, "WPC", &wpc, NULL, &status);
		FITS_read_key(pt_fits, TFLOAT, "OBSTIME", &obstime_key, NULL, &status);
		FITS_update_key(pt_fits_out, TFLOAT, "OBSTIME", &obstime_key, NULL, &status);
		FITS_read_key(pt_fits, TLONG, "NSPEC", &num_exp, NULL, &status);
		FITS_update_key(pt_fits_out, TLONG, "NSPEC", &num_exp, NULL, &status);

		for (j = 0L; j < num_exp; j++) {
		    sprintf(keyname, "SPEC%03ld", j+1L);
		    FITS_read_key(pt_fits, TSTRING, keyname, keyvalue, NULL, &status);
		    FITS_update_key(pt_fits_out, TSTRING, keyname, keyvalue, NULL, &status);
		    sprintf(keyname, "WOFF%03ld", j+1L);
		    fits_read_key(pt_fits, TFLOAT, keyname, &woffset, comment, &status);
		    if (status) {
			status = 0;
			woffset = 0.;
			sprintf(comment, "[A]");
		    }
		    FITS_update_key(pt_fits_out, TFLOAT, keyname, &woffset, comment, &status);
		}
	        FITS_close_file(pt_fits, &status);
	    }
	}
	if (n_spec != 8) {
		char datestr[FLEN_CARD];
		int timeref;
		FITS_movabs_hdu(pt_fits_out, 1, &hdutype, &status);
		FITS_write_comment(pt_fits_out, "  ", &status);
		sprintf(comment, "Segments");
		for (i=0L; i<8L; i++) if (!strncmp(spec_list[i], "NULL", 4))
			sprintf(comment, "%s %s", comment, extname[i]);
		sprintf(comment, "%s are missing.", comment);
		FITS_write_comment(pt_fits_out, comment, &status);
		sprintf(comment, "Extensions");
		for (i=0L; i<8L; i++) if (!strncmp(spec_list[i], "NULL", 4))
			sprintf(comment, "%s %ld", comment, i+1L);
		sprintf(comment, "%s of this file contain no data.", comment);
		FITS_write_comment(pt_fits_out, comment, &status);
		fits_get_system_time(datestr, &timeref, &status);
		sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
		FITS_write_comment(pt_fits_out, comment, &status);
		FITS_write_comment(pt_fits_out, "  ", &status);
 	}

	/* Add an index to the file extensions. */
	FITS_movabs_hdu(pt_fits_out, 1, &hdutype, &status);
        FITS_write_comment(pt_fits_out, "  ", &status);
        sprintf(comment, "Index to file extensions:");
        FITS_write_comment(pt_fits_out, comment, &status);
        sprintf(comment, "Extension 0   This Header");
        FITS_write_comment(pt_fits_out, comment, &status);
        for (i=0L; i<8L; i++) {
                sprintf(comment, "Extension %ld   %s", i+1L, extname[i]);
        	FITS_write_comment(pt_fits_out, comment, &status);
	}
        FITS_write_comment(pt_fits_out, "  ", &status);

	FITS_close_file(pt_fits_out, &status);
	
	return(0);
}
