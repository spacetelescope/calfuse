/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_gainmap input_files
 *
 * Description: Adds the counts from a list of input ttag idf files into
 *              the gain map.  The accumulated gain map cube is
 *              determined by reading the GAIN_CAL keyword in the
 *              input file's header.  To prevent multiple copies
 *              of reprocessed input files from being added to
 *              the total exposure map, the gain map contains an
 *              index list of the files which have been added.
 *              The new file is skipped if it is in the list.
 * 
 * 
 * Arguments:   input_files     Input ttag IDF FITS file name
 *
 * Calibration files required:  None
 *
 * Returns:     int  0 = successful completion
 *
 * Note:  This routine requires significantly more memory than your
 *        average routine.  Its design has not been optimized to 
 *        reduce memory usage.
 *
 * History:     03/31/03    1.1  peb  Finished work.
 *                                    Based on cf_ttag_countmap.
 *		06/11/03    1.2  wvd  Pass datatype to cf_read_col
 *		08/25/03    1.3  wvd  Change coltype from string to int in
 *					cf_read_col.
 *              01/09/04    1.4  bjg  Change CHIDCOR1 to CHID_COR
 *              02/17/04    1.5  bjg  Added scaling option 
 *                                    histogram files given in the argument 
 *                                    list are ignored
 *              02/26/04    1.6  bjg  Bug Fix. Routine name changed from 
 *                                    cf_ttag_gainmap to cf_gainmap.
 *                                    Added option to select coordinates
 *                                    to use: RAW, FARF (default), final
 *                                    Print ignored HIST files
 *              04/08/04    1.7  bjg  Remove unused function get_xy_columns
 *                                    Use cf_verbose instead of printf
 *              04/27/04    1.8  bjg  Added lower and upper limits option
 *                                    for PHA 
 *                                    Change author name of gainfile
 *              08/26/04    1.9  bjg  In force mode, issue a warning 
 *                                    if file already included. 
 *              03/02/05    1.10 wvd  Use FLOAT_IMG to insure that output
 *					images are written as floats.
 *              03/02/05    1.11 wvd  Don't write AUTHOR keyword to header.
 *              04/11/05    1.12 wvd  Force detector string to be upper case.
 *					Add -t option, which prevents creation
 *					of detector image for PHA=31.  It's
 *					required if you don't bin the data.
 *              06/03/05    1.13 wvd  Include ctype.h
 *              09/09/05    1.14 wvd  Update list of allowed options.
 *              08/24/07    1.15 bot  Changed nx and ny to long and called
 *				      as TLONG ; changed binx and biny to int
 *				      and called as TINT.
 *
 ****************************************************************************/

#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "calfuse.h"

#define ROOT_LENGTH 13
#define MAXCHARS 120
#define MASTER_CAL_FILE "master_calib_file.dat"
#define IMAGE_EXT 1
#define INDEX_EXT 34

static char CF_PRGM_ID[] = "cf_gainmap";
static char CF_VER_NUM[] = "1.15";


static int
new_gainfile(char *newfile, char *detector, float effmjd, int binx, int biny, int truncate)
{
    char  date[FLEN_CARD]={'\0'};
    char  *ttype[] = {"ROOTNAME", "OBSDATE", "COUNTS", "EXPTIME", "ADDDATE", "SCALE", "PHAMIN", "PHAMAX"};
    char  *tform[] = {"13A",      "19A",     "1D",     "1E",      "19A",     "1E",     "1I",    "1I"};
    char  *tunit[] = {" ",        " ",       "counts", "seconds", " ",       " ",      " ",     " "};
    int   naxis=2, status=0, version=1, tref=0, tfields=8, phaext;
    long  axis[2], naxis2=0;
    float exptime=0.;
    double ncounts=0.;
    fitsfile *gainfits;

    cf_verbose(3, "Entering subroutine new_gainfile");
    cf_verbose(3, "Creating file %s", newfile);

    FITS_create_file(&gainfits, newfile, &status);

    FITS_create_img(gainfits, 8, 0, NULL, &status);
    FITS_write_key(gainfits, TSTRING, "CALFTYPE", "GAIN",
		   "Calibration file type", &status);
    FITS_write_key(gainfits, TINT,    "CALFVERS", &version,
		   "Calibration file version", &status);
    FITS_write_key(gainfits, TSTRING, "DETECTOR", detector,
		   "detector (1A, 1B, 2A, 2B)", &status);
    FITS_write_key(gainfits, TFLOAT,  "EFFMJD",   &effmjd,
		   "Date on which file should be applied (MJD)", &status);
    fits_get_system_time(date, &tref, &status);
    FITS_write_key(gainfits, TSTRING, "DATE",     date,
		   "file creation date (YYYY-MM-DDThh:mm:ss UTC)", &status);
    FITS_write_key(gainfits, TDOUBLE, "NEVENTS",  &ncounts,
		   "Accumulated weighted events", &status);
    FITS_write_key(gainfits, TFLOAT,  "EXPTIME",  &exptime,
		   "Accumulated exposure time", &status);

    axis[0] = NXMAX/binx;
    axis[1] = NYMAX/biny;
    for (phaext=0; phaext<32-truncate; phaext++) {
	
	cf_verbose(3, "Creating extension %d", phaext);
	FITS_create_img(gainfits, FLOAT_IMG, naxis, axis, &status);
	FITS_write_key(gainfits, TINT,    "PHA",         &phaext,
		       "PHA value of image", &status);
	FITS_write_key(gainfits, TINT,    "SPECBINX",    &binx,
		       "Binsize in detector X coordinate", &status);
	FITS_write_key(gainfits, TINT,    "SPECBINY",    &biny,
		       "Binsize in detector Y coordinate", &status);

    }

    cf_verbose(3, "Creating binary table extension");
    FITS_create_tbl(gainfits, BINARY_TBL, naxis2, tfields, ttype, tform, tunit,
		    "PROCESSED FILES", &status);

    cf_verbose(3, "Closing file %s", newfile);
    FITS_close_file(gainfits, &status);
    cf_verbose(3, "Exiting new_gainfile");
    return status;
}

static int
_check_index(fitsfile *gainfits, char *rootname, float *scale, int *phamin,
	int *phamax, int truncate)
{
    int   status=0, anynull=0, hdutype;
    float *scal; 
    int *phami, *phama;
    long  j, nrows, frow=1, felem=1;
    char  **filename, strnull[] = "          ";

    FITS_movabs_hdu(gainfits, INDEX_EXT-truncate, &hdutype, &status);
    FITS_read_key(gainfits, TLONG, "NAXIS2", &nrows, NULL, &status);

    filename = (char **) cf_malloc(nrows * sizeof(char *));
    filename[0] = (char *) cf_malloc(nrows*(ROOT_LENGTH+1)*sizeof(char));
    scal = (float *) cf_malloc(nrows * sizeof(float));
    phami = (int *) cf_malloc(nrows * sizeof(int));
    phama = (int *) cf_malloc(nrows * sizeof(int));
    for (j=1; j<nrows; j++)
	filename[j] = filename[j-1]+(ROOT_LENGTH+1);

    FITS_read_col(gainfits, TSTRING, 1, frow, felem, nrows, strnull,
		  filename, &anynull, &status);
    FITS_read_col(gainfits, TFLOAT,  6, frow, felem, nrows, strnull,
		  scal, &anynull, &status);
    FITS_read_col(gainfits, TINT,  7, frow, felem, nrows, strnull,
		  phami, &anynull, &status);
    FITS_read_col(gainfits, TINT,  8, frow, felem, nrows, strnull,
		  phama, &anynull, &status);

    for (j=0; j<nrows; j++) {
	if (strncmp(filename[j], rootname, strlen(rootname)) == 0) {
	    *scale = scal[j];
	    *phamin = phami[j];
	    *phamax = phama[j];
	    free(filename[0]);
	    free(filename);
	    free(scal);
	    free(phami);
	    free(phama);
            return j+1;
	}
    }
    free(filename[0]);
    free(filename);
    free(scal);
    free(phami);
    free(phama);
    return 0;
}

static int
_update_index(fitsfile *gainfits, char *rootname, char *obsdate,
     double ncounts, float exptime, float scale, int phamin, int phamax,
     int truncate)
{
    char  *strptr[1], adddate[FLEN_CARD]={'\0'};
    int   status=0, hdutype, tref=0;
    long  nrows;

    FITS_movabs_hdu(gainfits, INDEX_EXT-truncate, &hdutype, &status);
    FITS_read_key(gainfits, TLONG, "NAXIS2", &nrows, NULL, &status);

    strptr[0] = rootname;
    FITS_write_col(gainfits, TSTRING, 1, nrows+1, 1, 1, strptr,   &status);
    strptr[0] = obsdate;
    FITS_write_col(gainfits, TSTRING, 2, nrows+1, 1, 1, strptr,   &status);
    FITS_write_col(gainfits, TDOUBLE, 3, nrows+1, 1, 1, &ncounts, &status);
    FITS_write_col(gainfits, TFLOAT,  4, nrows+1, 1, 1, &exptime, &status);
    fits_get_system_time(adddate, &tref, &status);
    strptr[0] = adddate;
    fits_write_col(gainfits, TSTRING, 5, nrows+1, 1, 1, strptr,   &status);
    FITS_write_col(gainfits, TFLOAT,  6, nrows+1, 1, 1, &scale,   &status);
    fits_write_col(gainfits, TINT,    7, nrows+1, 1, 1, &phamin,  &status);
    FITS_write_col(gainfits, TINT,    8, nrows+1, 1, 1, &phamax,  &status);
    return 0;
}


static double
add_events(fitsfile *infits, char *xcolumn, char *ycolumn, float scale,
	   int phamin, int phamax, fitsfile *gainfits)
{
    char   *pha;
    int    binx, biny, phaext, status=0, hdutype, anynull=0;
    long   nevents, j;
    long   nx, ny;
    float  *xpos, *ypos, *weight, *image;
    double ncounts=0.;

    cf_verbose(3,"Entering add_events");

    FITS_movabs_hdu(infits, 2, &hdutype, &status);
    nevents = cf_read_col(infits, TFLOAT, xcolumn,  (void **) &xpos);
    nevents = cf_read_col(infits, TFLOAT, ycolumn,  (void **) &ypos);
    nevents = cf_read_col(infits, TBYTE,  "PHA",    (void **) &pha);
    nevents = cf_read_col(infits, TFLOAT, "WEIGHT", (void **) &weight);

    for (phaext=phamin; phaext<=phamax; phaext++) {
	FITS_movabs_hdu(gainfits, phaext+2, &hdutype, &status);
	if (hdutype == BINARY_TBL) {
	    cf_if_warning("Can't find image extension for PHA = %d", phaext);
	    break;
	}
	FITS_read_key(gainfits, TLONG, "NAXIS1",   &nx,   NULL, &status);
	FITS_read_key(gainfits, TLONG, "NAXIS2",   &ny,   NULL, &status);
	FITS_read_key(gainfits, TINT, "SPECBINX", &binx, NULL, &status);
	FITS_read_key(gainfits, TINT, "SPECBINY", &biny, NULL, &status);
    
	image = (float *) cf_malloc(sizeof(float)*nx*ny);
	FITS_read_img(gainfits, TFLOAT, 1, nx*ny, NULL, image,
		      &anynull, &status);

	for (j=0; j<nevents; j++) {
	    if (pha[j] == phaext) {
		int x, y;
		x = xpos[j];
		y = ypos[j];
		if ((x >= 0) && (x < nx*binx) && (y >= 0) && (y < ny*biny)) {
		    image[(y/biny)*nx + (x/binx)] += weight[j]*scale;
		    ncounts += weight[j]*scale;
		}
	    }
	}
	FITS_write_img(gainfits, TFLOAT, 1, nx*ny, image, &status);
	free(image);
    }

    free(weight);
    free(pha);
    free(ypos);
    free(xpos);
    return ncounts;
}

static double
subtract_events(fitsfile *infits, char *xcolumn, char *ycolumn, 
		float scale, int phamin, int phamax, fitsfile *gainfits)
{
    char   *pha;
    long   nx, ny;
    int    phaext, status=0, hdutype, anynull=0, binx, biny;
    long   nevents, j;
    float  *xpos, *ypos, *weight, *image;
    double ncounts=0.;

    FITS_movabs_hdu(infits, 2, &hdutype, &status);
    nevents = cf_read_col(infits, TFLOAT, xcolumn,  (void **) &xpos);
    nevents = cf_read_col(infits, TFLOAT, ycolumn,  (void **) &ypos);
    nevents = cf_read_col(infits, TBYTE, "PHA",    (void **) &pha);
    nevents = cf_read_col(infits, TFLOAT, "WEIGHT", (void **) &weight);

    for (phaext=phamin; phaext<=phamax; phaext++) {
	FITS_movabs_hdu(gainfits, phaext+2, &hdutype, &status);
	FITS_read_key(gainfits, TLONG, "NAXIS1", &nx, NULL, &status);
	FITS_read_key(gainfits, TLONG, "NAXIS2", &ny, NULL, &status);
	FITS_read_key(gainfits, TINT, "SPECBINX", &binx, NULL, &status);
	FITS_read_key(gainfits, TINT, "SPECBINY", &biny, NULL, &status);

	image = (float *) cf_malloc(sizeof(float) * nx*ny);
	FITS_read_img(gainfits, TFLOAT, 1, nx*ny, NULL, image,
		      &anynull, &status);

	for (j=0; j<nevents; j++) {
	    if (pha[j] == phaext) {
		int x, y;
		x = xpos[j];
		y = ypos[j];
		if ((x >=0) && (x < nx*binx) && (y >= 0) && (y < ny*biny)) {
		    image[(y/biny)*nx + (x/binx)] -= weight[j]*scale;
		    ncounts += weight[j]*scale;
		}
	    }
	}
	FITS_write_img(gainfits, TFLOAT, 1, nx*ny, image, &status);
	free(image);
    }

    free(weight);
    free(pha);
    free(ypos);
    free(xpos);
    return ncounts;
}


int main(int argc, char *argv[])
{
    char  *newfile = NULL, *detector = NULL, buffer[FLEN_CARD]={'\0'}, instmode[FLEN_CARD]={'\0'};
    int   subtract=0, force=0, newbinx=8, newbiny=4, status=0, optc;
    int   truncate=0;
    float effmjd = 0.,scale=1., old_scale;
    int phamin=0, phamax=31; 
    int old_phamin, old_phamax;

    char   xcolumn[FLEN_CARD]={'\0'}, ycolumn[FLEN_CARD]={'\0'};

    char opts[] = "hv:fsrtcn:d:j:w:x:y:l:u:";
    char usage[] = 
	"Usage:\n"
	"  cf_gainmap [-hfsrtc] [-v level] [-n filename] [-d segment]\n"
	"    [-w scale] [-j effmjd] [-x xbin] [-y ybin] \n"
        "    [-l phamin] [-u phamax] idffiles\n";
    char option[] =
	"Options:\n"
	"  -h:  this help message\n"
        "  -v:  verbosity level (=1; 0 is silent)\n"
	"  -f:  force the files\n"
	"  -s:  subtract files\n"
	"  -n:  gain map file\n"
	"  -d:  gain map segment (no default for new file)\n"
        "  -j:  effective MJD (=0.0)\n"
        "  -w:  weighting factor for image (=1.0)\n"
	"  -x:  X-bin size (=8)\n"
	"  -y:  Y-bin size (=4)\n"
        "  -r:  use raw XY coordinates instead of farf\n"
        "  -t:  truncate map at PHA = 30\n"
        "  -c:  use final XY coordinates instead of farf\n"
        "  -l:  lower limit for pha (=0) \n"
        "  -u:  upper limit for pha (=31) \n";

    verbose_level = 1;

    strcpy(xcolumn, "XFARF");
    strcpy(ycolumn, "YFARF");

   

    /* Check number of options and arguments */
    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {
	    
	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
	case 'v':
	    verbose_level = atoi(optarg);
	    break;
	case 'f':
	    force = 1;
	    break;
	case 's':
	    subtract = 1;
	    break;
	case 'n':
	    newfile = optarg;
	    break;
	case 'd':
	    detector = optarg;
	    break;
	case 'j':
	    effmjd = atof(optarg);
	    break;
	case 'w':
	    scale = atof(optarg);
	    break;
	case 'x':
	    newbinx = atoi(optarg);
	    break;
	case 'y':
	    newbiny = atoi(optarg);
	    break;
	case 'r':
	    strcpy(xcolumn, "XRAW");
	    strcpy(ycolumn, "YRAW");
	    break;
	case 't':
	    truncate = 1;
	    break;
	case 'c':
	    strcpy(xcolumn, "X");
	    strcpy(ycolumn, "Y");
	    break;
	case 'l':
	    phamin = atoi(optarg);
	    break;
	case 'u':
	    phamax = atoi(optarg);
	    break;
	}
    }
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    if (argc <= optind)
        cf_if_error("%s\nIncorrect number of program arguments", usage);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin processing");
    
    while (optind < argc) {

	char   gainfile[FLEN_CARD], rootname[FLEN_CARD]={'\0'};
	char   segment[FLEN_CARD] = {'\0'};
	int    hdutype;
	FILE   *fptr;
	fitsfile *infits, *gainfits;

	FITS_open_file(&infits, argv[optind], READONLY, &status);
	FITS_movabs_hdu(infits, 1, &hdutype, &status);

	FITS_read_key(infits, TSTRING, "INSTMODE", instmode, NULL,
			  &status);

	if (strncmp(instmode,"TTAG",4)){
	  FITS_close_file(infits, &status);
	  cf_verbose(1,"Skipped HIST file: %s\n",argv[optind]);
	  optind ++;
	  continue;
	}

	if (newfile) {
	    cf_verbose(3, "newfile is TRUE");
	    strcpy(gainfile, newfile);

	    if ((fptr = fopen(gainfile, "r")))
		fclose(fptr);
	    else {
		if (!detector)
		    cf_if_error("No detector segment specified: "
				"1A, 1B, 2A, 2B");
		detector[1] = toupper(detector[1]);
		status = new_gainfile(gainfile, detector, effmjd,
				      newbinx, newbiny, truncate);
	    }
	    FITS_open_file(&gainfits, gainfile, READWRITE, &status);
	    FITS_movabs_hdu(gainfits, 1, &hdutype, &status);
	    detector = buffer;
	    FITS_read_key(gainfits, TSTRING, "DETECTOR", detector, NULL,
			  &status);
	} else {
	    cf_verbose(3, "newfile is FALSE");
	    FITS_read_key(infits, TSTRING, "GAIN_MAP", gainfile, NULL,
			  &status);
	    detector = buffer;
	    FITS_read_key(infits, TSTRING, "DETECTOR", detector, NULL,
			  &status);

	    if ((fptr = fopen(cf_hist_file(gainfile), "r"))) {
		fclose(fptr);
	    } else {
		char  keyword[FLEN_CARD]={'\0'}, segment[FLEN_CARD]={'\n'};
		char  linin[MAXCHARS]={'\0'}, filename[FLEN_CARD]={'\0'};
		int   interp;
		float aftermjd;
		FILE  *master;
		/*
		 * Get the effective MJD from the master calibration file
		 */
		master = fopen(cf_parm_file(MASTER_CAL_FILE), "r");
		if (master == NULL)
		    cf_if_error("Master calibration database file not found");
		while (fgets(linin, MAXCHARS, master) != NULL) {
		    /*
		     * Check for comment lines
		     */
		    if ((linin[0] != '#') && (linin[0] != '\n')) {
			sscanf(linin, "%4c%*2c%2c%*2c%13c%*2c%9f%*2c%1d",
			       keyword, segment, filename, &aftermjd, &interp);
			if (strcmp(filename, gainfile) == 0) {
			    effmjd = aftermjd;
			    break;
			}
		    }
		}
		fclose(master);

		status = new_gainfile(cf_hist_file(gainfile), detector,
				      effmjd, newbinx, newbiny, truncate);
	    }
	    FITS_open_file(&gainfits, cf_hist_file(gainfile), READWRITE,
			   &status);
	}

	strncpy(rootname, argv[optind], ROOT_LENGTH);
	FITS_read_key(infits, TSTRING,  "DETECTOR", segment, NULL, &status);
	if (strcmp(detector, segment) != 0) {
	    cf_if_warning("%s segment is not the same as %s",
			  argv[optind], gainfile);
	}
	else if (subtract) {
	    int  index;
	    if ((index = _check_index(gainfits, rootname, &scale, &phamin,
		&phamax, truncate))) {
		
		int    hdutype;
		float  exptime=0., acctime=0.;
		double ncounts=0, accevnt=0;

		FITS_read_key(infits, TFLOAT,  "EXPTIME", &exptime, NULL,
			      &status);

		
		ncounts = subtract_events(infits, xcolumn, ycolumn, scale, phamin,
			phamax-truncate, gainfits);

		FITS_movabs_hdu(gainfits, IMAGE_EXT, &hdutype, &status);
		FITS_read_key(gainfits, TDOUBLE, "NEVENTS",  &accevnt,
			      NULL, &status);
		FITS_read_key(gainfits, TFLOAT, "EXPTIME", &acctime,
			      NULL, &status);
		accevnt -= ncounts;
		acctime -= exptime;
		FITS_update_key(gainfits, TDOUBLE,  "NEVENTS", &accevnt,
				NULL, &status);
		FITS_update_key(gainfits, TFLOAT, "EXPTIME", &acctime,
				NULL, &status);

		FITS_movabs_hdu(gainfits, INDEX_EXT-truncate, &hdutype, &status);
		FITS_delete_rows(gainfits, index, 1, &status);
		cf_verbose(1, "%s is subtracted from %s",
			   argv[optind], gainfile);
	    } else {
		cf_if_warning("%s is not in %s.", argv[optind], gainfile);
	    }
	} else {
	    cf_verbose(3, "Adding new data file.");
	    if (!force && _check_index(gainfits, rootname, &old_scale, 
		&old_phamin, &old_phamax, truncate)) {
		cf_if_warning("%s is already in %s with scale = %f, "
		"phamin = %d, phamax=%d.",
		argv[optind], gainfile, old_scale, old_phamin, old_phamax);
	    } else {
		char   obsdate[FLEN_CARD]={'\0'}, dateobs[FLEN_CARD]={'\0'};
		char   timeobs[FLEN_CARD]={'\0'};
		
		float  exptime=0., acctime=0.;
		double ncounts=0, accevnt=0;

		FITS_read_key(infits, TSTRING, "DATEOBS", dateobs, NULL,
			      &status);
		FITS_read_key(infits, TSTRING, "TIMEOBS", timeobs, NULL,
			      &status);
		obsdate[0] = '\0';
		sprintf(obsdate, "%s%s%s", dateobs, "T", timeobs);
		FITS_read_key(infits, TFLOAT,  "EXPTIME", &exptime, NULL,
			      &status);

		if ((force) && _check_index(gainfits, rootname, &old_scale,
		   &old_phamin, &old_phamax, truncate)) 
		   cf_if_warning("FORCE MODE - %s is already in %s with scale = %f, "
			"phamin = %d, phamax=%d.",
			argv[optind], gainfile, old_scale, old_phamin, old_phamax);
		
		ncounts = add_events(infits, xcolumn, ycolumn, scale, phamin,
			phamax-truncate, gainfits);

		FITS_movabs_hdu(gainfits, IMAGE_EXT, &hdutype, &status);
		FITS_read_key(gainfits, TDOUBLE, "NEVENTS",  &accevnt,
			      NULL, &status);
		FITS_read_key(gainfits, TFLOAT, "EXPTIME", &acctime,
			      NULL, &status);
		accevnt += ncounts;
		acctime += exptime;
		FITS_update_key(gainfits, TDOUBLE,  "NEVENTS", &accevnt,
				NULL, &status);
		FITS_update_key(gainfits, TFLOAT, "EXPTIME", &acctime,
				NULL, &status);

		status = _update_index(gainfits, rootname, obsdate, ncounts,
				       exptime, scale, phamin, phamax, truncate);
		cf_verbose(1, "%s is added to %s", argv[optind], gainfile);
	    }
	}
	FITS_close_file(gainfits, &status);
	FITS_close_file(infits, &status);
	optind ++;
    }
    
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
    return 0;
}
