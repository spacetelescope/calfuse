/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_proc_update(fitsfile *fptr, char *prgm_id, char *key_value)
 *
 * Description: cf_proc_update will update the FITS header keyword 
 *              which corresponds to the given prgm_id in the file 
 *              fptr to the value given in key_value.
 *
 * Arguments:   fitsfile    *fptr	Pointer to input file
 *              char        *prgm_id 	Procedure name
 *              char        *key_value 	Updated procedure status
 *
 * History:     05/21/98        emm     Begin work.
 *              05/21/98        emm     finished
 *              06/17/98        emm     Modified so that the keyword_tab
 *                                      structure is read from ed_calfuse.h
 *                                      file, which simplifies revisions in
 *                                      the order of processing.
 *              09/08/99        peb     Added lines to update NEXTEND keyword
 *              01/31/00        emm     Added update to DATE keyword.
 *		04/01/03	wvd	Changed cf_errmsg to cf_if_warning,
 *					fits_modify_key_str to FITS_update_key,
 *					fits_read_key_str to FITS_read_key
 *              12/18/03        bjg     Change calfusettag.h to calfuse.h
 *              04/07/07  1.5   wvd     Delete CF_PRGM_ID, as it is not used.
 *
 ****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "calfuse.h"

int cf_proc_update(fitsfile *fptr, char *prgm_id, char *key_value) 
{

    char  comment[FLEN_CARD], instmode[FLEN_CARD];
    int i, status=0, hdutype=0, nextend;
    /*
     *  The calfuse.h file contains the definitions of keyword_tab
     *  NUM_PROC_STEPS, and CALIBRATION_STEP_KEYS.
     */
    struct keyword_tab keytab[NUM_PROC_STEPS]=CALIBRATION_STEP_KEYS;

    FITS_movabs_hdu(fptr, 1, &hdutype, &status);
    FITS_get_num_hdus(fptr, &nextend, &status);
    nextend -= 1;
    FITS_update_key(fptr, TINT, "NEXTEND", &nextend, NULL, &status);
    
    FITS_read_key(fptr, TSTRING, "INSTMODE", instmode, comment, &status);

    FITS_write_date(fptr, &status);

    switch (instmode[0]) {
    case 'H':
	/*  Find the keyword associated with prgm_id by looping 
	 *  through keytab[i].hist_proc */
	i=0;
	while ((strncmp(keytab[i].hist_proc,prgm_id,
			strlen(keytab[i].hist_proc)) != 0) &&
	       (i < NUM_PROC_STEPS)) i++;
	if (i < NUM_PROC_STEPS) {
	    /* We found a match to prgm_id, so change the associated
	     * keyword in the header. */
	    FITS_update_key(fptr, TSTRING,keytab[i].name,key_value, NULL, &status);
	} else {
	    /* The given prgm_id did not match any of the known
	     * keytab[i].hist_proc, so return 1. */
	    cf_if_warning("No PROCESSING STEP keyword is defined for %s", prgm_id);
	    return 1;
	}
	break;
    case 'T':
	/*  Timetag mode */
	/*  Find the keyword associated with prgm_id by looping 
	 *  through keytab[i].ttag_proc */
	i=0;
	while ((strncmp(keytab[i].ttag_proc,prgm_id,
			strlen(keytab[i].ttag_proc)) != 0) &&
	       (i < NUM_PROC_STEPS)) i++;
	if (i < NUM_PROC_STEPS) {
	    /* We found a match to prgm_id, so change the associated
	     * keyword in the header. */
	    FITS_update_key(fptr, TSTRING,keytab[i].name,key_value, NULL, &status);
	} else {
	    /* The given prgm_id did not match any of the known
	     * keytab[i].ttag_proc, so return 1. */
	    cf_if_warning("No PROCESSING STEP keyword is defined for %s", prgm_id);
	    return 1;
	}
	break;
    }
    return status;
}
