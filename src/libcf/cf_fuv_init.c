/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_fuv_init(fitsfile *fptr)
 *
 * Description: cf_fuv_init performs the following functions:
 *              1.  Populates the data processing step keywords based upon
 *                  the type of data (TTAG or HIST).  The processing steps
 *                  and the associated keywords are defined in the 
 *                  structure CALIBRATION_STEP_KEYS in calfuse.h.
 *              2.  Populates the orbital parameters in the header by
 *                  calling read_tle.  The orbital elements come from
 *                  the file FUSE.TLE, which is in standard NORAD two-line
 *                  element format.
 *              3.  Populates the calibration file keywords based on the 
 *                  time and date of observation, the detector segment, 
 *                  and the interpolation status of the calibration file.
 *                  The relation between processing steps and associated
 *                  calibration file keywords is defined in the structure
 *                  CALIBRATION_FILE_KEYS in calfuse.h.  The names and
 *                  effective dates of the calibration files are found in
 *                  the ASCII-format file master_calib_file.dat in the 
 *                  CF_PARMDIR directory.
 *
 * Input:   fitsfile    *fptr	Pointer to input file
 *
 * History:     05/11/98        emm     Begin work.
 *              05/15/98        emm     finished
 *              06/04/98        emm     Modified keywords to reflect updated
 *                                      calfuse design
 *              11/12/98        emm     Added ability to read orbital elements
 *                                      from file FUSE.TLE.
 *              03/18/99        emm     Added call to cf_cal_file.  master_
 *                                      calib_file.dat must now be in the
 *                                      CF_CALDIR directory.  (on 08/25/99
 *                                      this became CF_PARMDIR).
 *              03/18/99        emm     Edited header only.
 *              04/14/99        peb     Increased strings sizes to prevent
 *                                        character overflow.
 *                                      Implemented FITS wrapper functions.
 *              05/19/99        peb     Implemented cf_error_msg functions.
 *              06/07/99        peb     Added reporting of version number.
 *              07/29/99        emm     Added CF_VERS keyword update.
 *              08/25/99        emm     master_calib_file.dat now in CF_PARMDIR
 *              12/06/99 v1.8   emm     Removed cf_if_error call when master
 *                                      calib file keyword not found.
 *		10/13/00 v1.9	jwk	change CF_VER_NUM from "1.5"
 *		08/23/01 v1.10  wvd     Read PARM file; if requested, set
 *					BRST_COR, ASTG_COR, PHAX_COR, 
 *					MKBK_COR, or BKGD_COR to "OMIT"
 *		10/23/01 v1.11  hsu     Check OPUS version
 *		10/26/01 v1.12  wvd     Install subroutine.
 *		02/11/02 v1.13  wvd     Check detector bias levels.
 *		02/13/02 v1.14  wvd     Modified bias logic; write warning
 *					to data header as COMMENT lines.
 *		02/13/02 v1.15  wvd     Don't crash if bias keywords missing.
 *		02/19/02 v1.16  wvd     Check whether voltage too high
 *		02/20/02 v1.17  wvd     If keyword RAWTIME does not exist,
 *					create with value of EXPTIME
 *              10/30/02  1.18  peb     Converted header to calfusettag.h
 *		03/26/03  1.4   wvd	Removed "Warning" from warnings.
 *              04/16/03 v1.5   wvd     Test for corrupted bias values.
 *                                      Update comment field for HV_FLAG.
 *		04/28/03 v1.6   wvd	Change logic of assigning cal files
 *                                      never to extrapolate forward in time.
 *		07/23/03 v1.7   wvd	Write HSKP_CAL to file header.
 *		07/23/03 v1.8   wvd	Change calfusettag.h to calfuse.h
 *		09/09/03 v1.9   wvd	Read SAA_SCR and LIMB_SCR from
 *					SCRN_CAL and set RUN_SAA and RUN_LIMB.
 *		04/27/04 v1.10  wvd	Keywords MKBK_COR and BKGD_COR are
 *					gone, so we don't populate them.
 *		06/04/04 v1.11  wvd	Allow for the possibility that the
 *					time associated with a calibration
 *					file exactly equals EXPSTART.
 *		06/14/04 v1.12  wvd	Check value of BRIT_OBJ keyword and
 *					modify header accordingly.
 *					Write CF_VERS to trailer file.
 *		07/16/04 v1.13  wvd	Must set ASTG_COR to PERFORM, as
 *					pipeline doesn't do it automatically.
 *		11/26/04 v1.14  wvd	Check whether SAA_SCR and LIMB_SCR
 *					are ON or OFF, not YES or NO.
 *		12/02/04 v1.15  wvd	Read value of RUN_JITR from PARM file.
 *		03/24/05 v1.16  wvd	Use HVGOODLM to determine whether
 *					detector voltage is low.
 *		06/12/06 v1.17  wvd	Call cf_verbose rather than
 *					cf_if_warning for non-optimal voltage.
 *              07/18/08 v1.18  wvd     For bright-earth and airglow exposures, 
 *					change EXP_STAT to 2, 
 *					change SRC_TYPE to EE, and
 *					write warning to file header.
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calfuse.h"

#define  MAXCHARS 120
#define  MASTER_CAL_FILE "master_calib_file.dat"

int cf_fuv_init(fitsfile *fptr)
{
    char CF_PRGM_ID[] = "cf_fuv_init";
    char CF_VER_NUM[] = "1.18";

    int   i, status=0, expstat, interp_in, timeref;
    int   n, hv_flag, lowv, highv, saav, fullv;
    char  comment[FLEN_CARD], instmode[FLEN_CARD], detector[FLEN_CARD];
    char  opus_str[FLEN_CARD], keyword_in[FLEN_CARD], segment_in[FLEN_CARD];
    char  datestr[FLEN_CARD];
    char  filename_in[]="                    ", rootname[FLEN_CARD];
    char  hkexists[FLEN_CARD];
    char  keyword_out1[FLEN_CARD], keyword_out2[FLEN_CARD];
    char  linin[MAXCHARS];
    char  run_limb[FLEN_VALUE], run_saa[FLEN_VALUE], src_type[FLEN_VALUE];
    char  run_brst[FLEN_VALUE], run_walk[FLEN_VALUE], run_astg[FLEN_VALUE];
    /* char  run_mkbk[FLEN_VALUE], run_bkgd[FLEN_VALUE], brit_obj[FLEN_VALUE]; */
    char  run_jitr[FLEN_VALUE], brit_obj[FLEN_VALUE];
    char  hv_str[FLEN_CARD], lv_str[FLEN_CARD], mjd_str[FLEN_CARD];
    char  saa_str[FLEN_CARD], full_str[FLEN_CARD];
    float hvgood, mjdv, expstart, aftermjd_in;
    fitsfile *parmfits, *scrnfits;

    struct keyword_tab  keytab[NUM_PROC_STEPS]=CALIBRATION_STEP_KEYS;
    struct cal_file_tab calkey[NUMCALKEYS]=CALIBRATION_FILE_KEYS;
    FILE  *fpin;

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin initializing headers");

    /* Check the OPUS version number */
    fits_read_key(fptr, TSTRING, "OPUSVERS", opus_str, NULL, &status);
    if (status) {
        status = 0;
	cf_add_header_keywords(fptr);
	cf_verbose(2, "Data file generated by old version of OPUS.");
    }
    else if (atof(opus_str) < OPUS_VERSION) {
	cf_add_header_keywords(fptr);
	cf_verbose(2, "Data file generated by old version (%s) of OPUS.", opus_str);
    }
    FITS_update_key(fptr, TSTRING, "CF_VERS", CALFUSE_VERSION, NULL, &status);
    cf_verbose (1, "This is CalFUSE v%s", CALFUSE_VERSION);

    /*
     *  Populate data processing keywords based on observing mode found
     *  in the INSTMODE keyword.  NOTE:  We check only the first 
     *  letter of the INSTMODE keyword.
     */
    FITS_read_key(fptr, TSTRING, "INSTMODE", instmode, NULL, &status);

    switch (instmode[0]) {
    case 'H':
	/*  Histogram mode */
	for (i=0; i<NUM_PROC_STEPS; i++)
	    FITS_update_key(fptr, TSTRING, keytab[i].name,
			    keytab[i].hist_value, 0, &status);
	break;
    case 'T':
	/*  Time tag mode */
	for (i=0; i<NUM_PROC_STEPS; i++)
	    FITS_update_key(fptr, TSTRING, keytab[i].name,
			    keytab[i].ttag_value, 0, &status);
	break;
    default:
	/*  Invalid mode */
	cf_if_error("Invalid INSTMODE value");
	FITS_update_key(fptr, TSTRING, "INIT_COR", "ERROR", 0, &status);
	break;
    } /* End switch (instmode[0]) */

    /*  Populate the orbital parameters from a file of two-line elements. */
    read_tle(fptr);

    /*  Populate calibration file keywords based on the date of the
     *  observation and the segment/channel number. */

    FITS_read_key(fptr, TFLOAT, "EXPSTART", &expstart, NULL, &status);
    FITS_read_key(fptr, TSTRING, "DETECTOR", detector, NULL, &status);

    /*  Open the Master calibration database file. */

    fpin = fopen(cf_parm_file(MASTER_CAL_FILE), "r");
    if (fpin == NULL)
	cf_if_error("Master calibration database file not found");

    /*  Cycle through the master cal file until we run out of lines. */

    /*  Here's how this section works.  We need to keep track of 3 calibration
     *  files for each keyword: the cal file which immediately preceedes the 
     *  date of observation, the cal file which immediately preceedes the cal
     *  file which immediately preceedes the observation, and the cal file 
     *  immediately after the observation.  These are stored in strut calkey
     *  in:
     *     calkey[i].filenames[0] => second cal file preceeding observation
     *     calkey[i].filenames[1] => cal file immediately preceeding observa-
     *                                 tion
     *     calkey[i].filenames[2] => cal file immediately after observation
     *
     *  Here is the plan.  Cycle through the master cal file, reading each
     *  line.  Determine the keyword for each line, and use that to determine
     *  the appropriate index [i] for calkey.  Then, if the detector segment 
     *  is correct, determine where the new file fits into the structure above.
     *  If it is more recent than one or both of the preceeding cal files, 
     *  but is older than the observation (expstart) replace either 
     *  filenames[0] or filenames[1].  If filenames[1] is replaced, shift
     *  the old value down into filenames[0].  If the file is more recent 
     *  than the observation, and closer to the observation than the file 
     *  in filenames[2], replace filenames[2].  The values 
     *  calkey[i].aftermjd[] and calkey[i].interp[] store the date and 
     *  interpolation status of the best files.  This organization allows
     *  the master calibration database file to be written in any order.
     */

    while (fgets(linin, MAXCHARS, fpin) != NULL) {
        /* Check for comment lines */

	if ((linin[0] != '#') && (linin[0] != '\n')) {
	    sscanf(linin, "%4c%*2c%2c%*2c%13c%*2c%9f%*2c%1d",keyword_in,
		   segment_in,filename_in,&aftermjd_in,&interp_in);
	    /*
	     *  Determine which keyword we just read in.
	     */
	    i=0;
	    while ((strncmp(calkey[i].name,keyword_in,4) != 0) &&
		   (i < NUMCALKEYS)) i++;
	    if (i >= NUMCALKEYS) {
		cf_if_warning("Keyword %s in %s is unknown", keyword_in, MASTER_CAL_FILE);
	    } else {
		if ((strncmp(segment_in,"  ",2) == 0) ||
		    (strncmp(segment_in,detector,2) == 0)) {
		    if ((aftermjd_in <= expstart) &&
			(aftermjd_in > calkey[i].aftermjd[0]) &&
			(aftermjd_in > calkey[i].aftermjd[1])) {
			strncpy(calkey[i].filenames[0],
				calkey[i].filenames[1],18);
			calkey[i].aftermjd[0]=calkey[i].aftermjd[1];
			calkey[i].interp[0]=calkey[i].interp[1];
			strncpy(calkey[i].filenames[1],
				filename_in,18);
			calkey[i].aftermjd[1]=aftermjd_in;
			calkey[i].interp[1]=interp_in;
		    } else if ((aftermjd_in < expstart) &&
			       (aftermjd_in > calkey[i].aftermjd[0])) {
			strncpy(calkey[i].filenames[0],
				filename_in,18);
			calkey[i].aftermjd[0]=aftermjd_in;
			calkey[i].interp[0]=interp_in;
		    }
		    if ((aftermjd_in > expstart) &&
			(aftermjd_in < calkey[i].aftermjd[2])) {
			strncpy(calkey[i].filenames[2],
				filename_in,18);
			calkey[i].aftermjd[2]=aftermjd_in;
			calkey[i].interp[2]=interp_in;
		    }
		} /* Endif segment match */
	    } /* end else unrecognized keyword */		   
	}  /* end while linin not comment */
    } /* end while(linin != NULL) */

    fclose(fpin);

/*************************************************************************
        I have changed Rule 3 so that we NEVER extrapolate calibration
        files forward in time.  Instead, we just use the most recent file.
        This reduces from 8 to 4 the number of possible cases, since we
        no longer care about the value of interp[0].
        This section and the previous one ought to be re-written,
        but I will simply change 0 to 1 for Case 110 below.
        WVD - 04/28/03
*************************************************************************/

    /*  OK, now that we have identified the closest two files before
     *  the observation, and the file after the observation, we can
     *  decide what to do with them.  First, the easy part:  if 
     *  calkey[i].numfiles=1, no interpolation is possible, so we
     *  write out calkey[i].filenames[1], which is the file immediately
     *  preceeding the observation.  If calkey[i].numfiles=2, interpolation
     *  is possible, and we will have to write two keywords into the header.
     *  The decision as to which two keywords to write follows these rules:
     *  Rule 1: If an observation occurs between two files which can be 
     *          interpolated, then an interpolation must be done. 
     *  Rule 2: If the calibration file immediately preceeding an 
     *          observation is to be used in a stepwise manner, then said
     *          file shall be used in stepwise manner. 
     *  Rule 3: If the two calibration files immediately preceeding an 
     *          observation can be interpolated, but the following file 
     *          cannot, then the two previous files will be extrapolated 
     *          forward in time. 
     *  Rule 4: If the closest preceeding calibration file has an interpolation
     *          flag, but the file preceeding it has a stepwise flag, the the 
     *          closest preceeding file will be used in a sterpwise manner.
     *  Rule 5: Never interpolate backward in time.
     *
     *  In fact, there are only 8 possible cases:
     *
     *  Remember, the files were sorted such that the observation 
     *  always occurs between [1] and [2].  Also remember that 
     *  interp[i]=0 => stepwise):
     *
     *    Case               Action
     *
     *  interp[0] = 0
     *  interp[1] = 0    filenames[1] will be used in a stepwise manner
     *  interp[2] = 0
     *
     *  interp[0] = 0
     *  interp[1] = 0    filenames[1] will be used in a stepwise manner
     *  interp[2] = 1
     *
     *  interp[0] = 0
     *  interp[1] = 1    filenames[1] will be used in a stepwise manner
     *  interp[2] = 0
     *
     *  interp[0] = 1
     *  interp[1] = 0    filenames[1] will be used in a stepwise manner
     *  interp[2] = 0
     *
     *  interp[0] = 1
     *  interp[1] = 0    filenames[1] will be used in a stepwise manner
     *  interp[2] = 1
     *
     *  interp[0] = 0
     *  interp[1] = 1    filenames[1] and filenames[2] will be interpolated
     *  interp[2] = 1
     *
     *  interp[0] = 1
     *  interp[1] = 1    filenames[0] and filenames[1] will be extrapolated
     *  interp[2] = 0
     *
     *  interp[0] = 1
     *  interp[1] = 1    filenames[1] and filenames[2] will be interpolated
     *  interp[2] = 1
     *
     */

    /*  Write the calibration files into the header keywords. */
    /*  Cycle through the keywords. */
    for (i=0; i<NUMCALKEYS; i++) {
	if (calkey[i].numfiles == 1) {
	    /*  This is the easy case, no interpolation is possible. */

	    sprintf(keyword_out1,"%4.4s%1.1s%3.3s",calkey[i].name,"_",
		    calkey[i].extension);

	    FITS_update_key(fptr, TSTRING, keyword_out1,
			    calkey[i].filenames[1], 0, &status);
	}
	else {
	    /*  OK, interpolation is possible, now decide which two files
	     *  to write out. */

	    /* Create the two output keywords */
	    sprintf(keyword_out1,"%4.4s%1.1s%3.3s",calkey[i].name,"1",
                    calkey[i].extension);
	    sprintf(keyword_out2,"%4.4s%1.1s%3.3s",calkey[i].name,"2",
                    calkey[i].extension);

	    if (calkey[i].interp[1] == 0) {
		/*  Easy, preceeding file is a 0, so write filenames[1] to
		 *  both keywords.  Takes care of cases 000, 001, 100, 101.
		 */
		FITS_update_key(fptr, TSTRING, keyword_out1,
				calkey[i].filenames[1], 0, &status);
		FITS_update_key(fptr, TSTRING, keyword_out2,
				calkey[i].filenames[1], 0, &status);
	    }
	    else if ((calkey[i].interp[1] == 1) &&
		     (calkey[i].interp[2] == 1)) {
		/* Case 011 and 111 */
		FITS_update_key(fptr, TSTRING, keyword_out1,
				calkey[i].filenames[1], 0, &status);
		FITS_update_key(fptr, TSTRING, keyword_out2,
				calkey[i].filenames[2], 0, &status);
	    } else if ((calkey[i].interp[0] == 1) &&
		       (calkey[i].interp[1] == 1) &&
		       (calkey[i].interp[2] == 0)){
		/* Case 110 */
	    /* CHANGED 0 to 1 in the following line - wvd 04/28/03 */
		FITS_update_key(fptr, TSTRING, keyword_out1,
				calkey[i].filenames[1], 0, &status);
		FITS_update_key(fptr, TSTRING, keyword_out2,
				calkey[i].filenames[1], 0, &status);
	    }
	    else if ((calkey[i].interp[0] == 0) && (calkey[i].interp[1] == 1)
		     && (calkey[i].interp[2] == 0)){
		/* Case 010 */
		FITS_update_key(fptr, TSTRING, keyword_out1,
				calkey[i].filenames[1], 0, &status);
		FITS_update_key(fptr, TSTRING, keyword_out2,
				calkey[i].filenames[1], 0, &status);
	    }
	    else
                cf_if_error("The interpolation/stepwise logic has"
                            "failed.  I am confused.");
	} /* end else calkey.numfiles == 1 */
    } /* end for i=0; i<NUMCALKEYS */

    /*
     *  If housekeeping file exists,
     *	add HSKP_CAL and JITR_CAL filenames to the header.
     */
    FITS_read_key(fptr, TSTRING, "HKEXISTS", hkexists, NULL, &status);
    if(!strncasecmp(hkexists, "Y", 1)) {
	FITS_read_key(fptr, TSTRING, "ROOTNAME", rootname, NULL, &status);
	sprintf(keyword_out1, "%s%s", rootname, "hskpf.fit");
	FITS_update_key(fptr, TSTRING, "HSKP_CAL", keyword_out1, 0, &status);
	sprintf(keyword_out2, "%s%s", rootname, "jitrf.fit");
	FITS_update_key(fptr, TSTRING, "JITR_CAL", keyword_out2, 0, &status);
    }
    /*
     *  Check min and max voltage levels and set HV_FLAG accordingly.
     */
    sprintf(hv_str, "DET%cHV%cH", detector[0], detector[1]);
    sprintf(lv_str, "DET%cHV%cL", detector[0], detector[1]);

    fits_read_key(fptr, TINT, hv_str, &highv, NULL, &status);
    if (status) {
        status = 0;
        hv_flag = -1;
        sprintf(comment, "Detector bias keywords not found");
        FITS_update_key(fptr, TINT, "HV_FLAG", &hv_flag, comment, &status);
        cf_verbose(2, comment);
    }
    else {
        FITS_read_key(fptr, TINT, lv_str, &lowv, NULL, &status);
        FITS_read_key(fptr, TSTRING, "VOLT_CAL", filename_in, NULL, &status);
        FITS_open_file(&parmfits, cf_cal_file(filename_in), READONLY, &status);
        FITS_read_key(parmfits, TFLOAT, "HVGOODLM", &hvgood, NULL, &status);

        n = 0;
        do {
            n++;
            sprintf(mjd_str, "MJD%d", n);
            FITS_read_key(parmfits, TFLOAT, mjd_str, &mjdv, NULL, &status);
        } while (expstart > mjdv);
        n--;

        sprintf(saa_str, "SAA%d", n);
        sprintf(full_str, "FULL%d", n);
        FITS_read_key(parmfits, TINT, saa_str, &saav, NULL, &status);
        FITS_read_key(parmfits, TINT, full_str, &fullv, NULL, &status);
        FITS_close_file(parmfits, &status);

             if (lowv > fullv + 3) {
                hv_flag = 5;
                sprintf(comment, "Detector voltage high throughout exposure");
            }
        else if ((float) lowv / fullv > hvgood) {
                hv_flag = 4;
                sprintf(comment, "Detector voltage full throughout exposure");
            }
        else if (lowv > saav  + 3) {
                hv_flag = 3;
                sprintf(comment, "Detector voltage low throughout exposure");
            }
        else if (lowv > saav  - 3) {
                hv_flag = 2;   
                sprintf(comment, "Detector at SAA voltage throughout exposure");
            }
        else {
                hv_flag = 1;
                sprintf(comment, "Detector below SAA voltage throughout exposure");
            }

             if (highv - lowv > 3) {
                hv_flag = 0;
                sprintf (comment, "Detector voltage changed during exposure");
            }
        else if (lowv > highv) {
                hv_flag = -1;
                sprintf(comment, "Detector bias keywords are corrupted");
            }

        FITS_update_key(fptr, TINT, "HV_FLAG", &hv_flag, comment, &status);

        if (hv_flag != 4) {
            if (hv_flag == -1) cf_if_warning(comment);
            else cf_verbose(1, comment);
            FITS_write_comment(fptr, "  ", &status);
            FITS_write_comment(fptr, comment, &status);
	    fits_get_system_time(datestr, &timeref, &status);
	    sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
            FITS_write_comment(fptr, comment, &status);
            FITS_write_comment(fptr, "  ", &status);
        }
    }

    /*  Read BRIT_OBS keyword and modify header accordingly. */

    FITS_read_key(fptr, TSTRING, "SRC_TYPE", src_type, NULL, &status);
    FITS_read_key(fptr, TSTRING, "BRIT_OBJ", brit_obj, NULL, &status);
    if (!strcmp(brit_obj, "DEFOCUS") && !strncmp(src_type, "P", 1)) {
	*src_type = 'E';
	FITS_update_key(fptr, TSTRING, "SRC_TYPE", src_type,
		"DEFOCUS OBS ==> Treat as extended source", &status);
    }

    /* Flag bright-earth and airglow (M106, S100, and 900+) exposures. */

    FITS_read_key(fptr, TSTRING, "ROOTNAME", rootname, NULL, &status);
    if (!strncmp(rootname, "M106", 4) || !strncmp(rootname, "S100", 4)) {

        cf_verbose(1, "Bright-earth exposure.");

	/* Set EXP_STAT to TEMPORAL_LIMB. */
	expstat = (int) TEMPORAL_LIMB;
	FITS_update_key(fptr, TINT, "EXP_STAT", &expstat, 
		"Bright-earth exposure", &status);
    }

   if (!strncmp(rootname+8, "9", 1)) {

        cf_verbose(1, "Airglow exposure.  Will treat as extended source.");

	/* Set EXP_STAT to TEMPORAL_LIMB. */
	expstat = (int) TEMPORAL_LIMB;
	FITS_update_key(fptr, TINT, "EXP_STAT", &expstat, 
		"Airglow exposure", &status);

	/* Change SRC_TYPE to EE. */
	(void) strcpy(src_type, "EE");
	FITS_update_key(fptr, TSTRING, "SRC_TYPE", src_type,
		"Airglow exposure ==> Treat as extended source", &status);

        /* Write warning to trailer file and IDF header. */
	FITS_write_comment(fptr, " ", &status);
	sprintf(comment, "Airglow exposure. Not an astrophysical target.");
	FITS_write_comment(fptr, comment, &status);
	fits_get_system_time(datestr, &timeref, &status);
	sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
	FITS_write_comment(fptr, comment, &status);
	FITS_write_comment(fptr, " ", &status);
    }

    /*  Read the SCRN file and see if user wants to skip any steps. */

    FITS_read_key(fptr, TSTRING, "SCRN_CAL", filename_in, NULL, &status);
    FITS_open_file(&scrnfits, cf_parm_file(filename_in), READONLY, &status);
    FITS_read_key(scrnfits, TSTRING, "SAA_SCR", run_saa, NULL, &status);
    FITS_read_key(scrnfits, TSTRING, "LIMB_SCR", run_limb, NULL, &status);
    FITS_close_file(scrnfits, &status);

    /*  Now read the PARM file and see if user wants to skip any steps. */

    FITS_read_key(fptr, TSTRING, "PARM_CAL", filename_in, NULL, &status);
    FITS_open_file(&parmfits, cf_parm_file(filename_in), READONLY, &status);
    FITS_read_key(parmfits, TSTRING, "RUN_BRST", run_brst, NULL, &status);
    FITS_read_key(parmfits, TSTRING, "RUN_PHAX", run_walk, NULL, &status);
    /* FITS_read_key(parmfits, TSTRING, "RUN_MKBK", run_mkbk, NULL, &status); */
    /* FITS_read_key(parmfits, TSTRING, "RUN_BKGD", run_bkgd, NULL, &status); */
    FITS_read_key(parmfits, TSTRING, "RUN_ASTG", run_astg, NULL, &status);
    FITS_read_key(parmfits, TSTRING, "RUN_JITR", run_jitr, NULL, &status);
    FITS_close_file(parmfits, &status);

    if (!strncmp(run_limb, "OFF", 3) || !strncmp(run_limb, "off", 3))
	FITS_update_key(fptr, TSTRING, "LIMB_COR", "OMIT", NULL, &status);

    if (!strncmp(run_saa, "OFF", 3) || !strncmp(run_saa, "off", 3))
	FITS_update_key(fptr, TSTRING, "SAA__COR", "OMIT", NULL, &status);

    if (!strncmp(run_brst, "N", 1) || !strncmp(run_brst, "n", 1))
	FITS_update_key(fptr, TSTRING, "BRST_COR", "OMIT", NULL, &status);

    if (!strncmp(run_walk, "N", 1) || !strncmp(run_walk, "n", 1))
	FITS_update_key(fptr, TSTRING, "PHAX_COR", "OMIT", NULL, &status);

    /*
    if (!strncmp(run_mkbk, "N", 1) || !strncmp(run_mkbk, "n", 1)) 
	FITS_update_key(fptr, TSTRING, "MKBK_COR", "OMIT", NULL, &status);

    if (!strncmp(run_bkgd, "N", 1) || !strncmp(run_bkgd, "n", 1)) 
	FITS_update_key(fptr, TSTRING, "BKGD_COR", "OMIT", NULL, &status);
    */

    if (!strncmp(run_astg, "N", 1) || !strncmp(run_astg, "n", 1))
	FITS_update_key(fptr, TSTRING, "ASTG_COR", "OMIT", NULL, &status);
    else
	FITS_update_key(fptr, TSTRING, "ASTG_COR", "PERFORM", NULL, &status);

    if (!strncmp(run_jitr, "N", 1) || !strncmp(run_jitr, "n", 1))
	FITS_update_key(fptr, TSTRING, "JITR_COR", "OMIT", NULL, &status);
    else
	FITS_update_key(fptr, TSTRING, "JITR_COR", "PERFORM", NULL, &status);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done initializing headers");
    return status;
}
