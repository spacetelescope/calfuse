/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_proc_check(fitsfile *fptr, char *prog_id)
 *
 * Description: cf_proc_check will determine whether a given calibration
 *              step is to be performed on the data.
 *
 * Arguments:   fitsfile    *fptr	Pointer to input file
 *              char        *prog_id 	Procedure name
 *
 * History:     05/21/98        emm     Begin work.
 *              05/21/98        emm     finished
 *              06/17/98        emm     Modified so that the keyword_tab
 *                                      structure is read from ed_calfuse.h
 *                                      file, which simplifies revisions in
 *                                      the order of processing.
 *              03/04/99        emm     I added use of errflg instead of
 *                                      return so that it checks all parms.
 *		04/01/03	wvd	Change fits_read_key_str to 
 *					FITS_read_key, delete fits_print_err
 *              05/20/03        rdr     remove part which checks previous steps
 *              05/22/03 v1.6   wvd     Modify i/o; implement cf_verbose
 *              07/23/03 v1.7   wvd     Modify i/o
 *              09/10/03 v1.8   wvd     Bug fix: was printing HIST values
 *					whenever called with verbose option.
 *					Now returns TTAG when appropriate.
 *              03/17/04 v1.9   wvd     Call cf_verbose rather than
 *					cf_if_warning when program not
 *					appropriate for HIST data.
 *              06/02/04 v1.10  wvd     Call cf_verbose rather than
 *					cf_if_warning when program not
 *					appropriate for TTAG data.
 *              11/26/04 v1.11  wvd     If header keyword is set to OMIT, skip
 *					this step and set keyword to SKIPPED.
 *              03/30/05 v1.12  wvd     If header keyword is set to SKIPPED,
 *					skip this step.
 *
 ****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "calfuse.h"


int cf_proc_check(fitsfile *fptr, char *prog_id) 
{
    char  comment[FLEN_CARD], instmode[FLEN_CARD];
    char  key_value[FLEN_CARD];
    int   j=0, errflg=0, status=0;

    /*  The calfuse.h file contains the definitions of keyword_tab,
     *  NUM_PROC_STEPS, and CALIBRATION_STEP_KEYS.
     */
    struct keyword_tab keytab[NUM_PROC_STEPS]=CALIBRATION_STEP_KEYS;

    FITS_read_key(fptr, TSTRING, "INSTMODE", instmode, comment, &status);
   
    switch (instmode[0]) {
    case 'H':
	/*  Histogram mode */
	/*  First, determine if this procedure is even supposed to be
	 *  run on this data. */

	while ((strncmp(keytab[j].hist_proc,prog_id,
			strlen(keytab[j].hist_proc)) != 0) &&
	       (j < NUM_PROC_STEPS)) j++;
	if (j < NUM_PROC_STEPS) 
	    cf_verbose(3, "cf_proc_check: prog_id=%s, key=%s, j=%d, value=%s",
               prog_id, keytab[j].hist_proc, j, keytab[j].hist_value) ;
	if ((j >= NUM_PROC_STEPS) || 
	    (strncmp(keytab[j].hist_value, "PERFORM", 7))) {
	    cf_verbose(1, "Not appropriate for HIST data.  Exiting.");
	    errflg=1;
	}

	break;

    case 'T':
	/*  Time-tag mode */
	/*  First, determine if this procedure is even supposed to be
	 *  run on this data. */
	while ((strncmp(keytab[j].ttag_proc,prog_id,
			strlen(keytab[j].ttag_proc)) != 0) &&
	       (j < NUM_PROC_STEPS)) j++;
	if (j < NUM_PROC_STEPS) 
	    cf_verbose(3, "cf_proc_check: prog_id=%s, key=%s, j=%d, value=%s",
               prog_id, keytab[j].ttag_proc, j, keytab[j].ttag_value) ;
	if ((j >= NUM_PROC_STEPS) || 
	    (strncmp(keytab[j].ttag_value, "PERFORM", 7))) {
	    cf_verbose(1, "Not appropriate for TTAG data.  Exiting.");
	    errflg=1;
	}

	break;

    default:
	errflg = 1;
    }

    if (!errflg) {
	/* Check to see whether user wants to skip this step. */
	FITS_read_key(fptr, TSTRING, keytab[j].name, key_value, comment,
		&status);
	if (!strncmp(key_value, "OMIT", 4) || !strncmp(key_value, "SKIPPED", 7)) {
	    cf_verbose(1, "Header keyword %s = %s.  Exiting.",
		keytab[j].name, key_value);
	    FITS_update_key(fptr, TSTRING, keytab[j].name, "SKIPPED", NULL,
		&status);
	    errflg=1;
	}

	/* Now check to see if the step has already been completed. */
	if (!strncmp(key_value, "COMPLETE", 7)) {
	    cf_if_warning("Exiting.\n\t%s has already been run on this file.", 
		prog_id);
	    errflg=1;
	}
    }
    return errflg;
}
