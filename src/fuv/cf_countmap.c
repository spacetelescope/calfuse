/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_countmap input_files
 *
 * Description: Adds the counts from a list of input idf files into
 *              the count map. The accumulated count map cube is
 *              determined by reading the COUNT_CAL keyword in the
 *              input file's header. To prevent multiple copies
 *              of reprocessed input files from being added to
 *              the total exposure map, if the input file has been
 *              previously added, it is skipped over and a message is
 *              printed out.
 *
 * 
 * Arguments:   input_files     Input initialized ttag FITS file names
 *
 * Calibration files required:  None
 *
 * Returns:     int  0 = successful completion
 *
 * Note:  This routine requires significantly more memory than your
 *        average routine.  Its design has not been optimized to 
 *        reduce memory usage.
 *
 * History:     03/25/03    1.1  peb  Begin and finish work. Based on
 *                                    cf_ttag_to_image
 *              03/31/03    1.2  peb  Removed filtering of photons and fixed
 *                                    exposure time bug when subtracting data
 *		06/11/03    1.3  wvd  Pass datatype to cf_read_col
 *		08/25/03    1.4  wvd  Change coltype from string to int in
 *					cf_read_col.
 *		08/25/03    1.5  wvd  Change CHIDCOR1 to CHID_COR
 *              02/17/04    1.6  bjg  Added scaling option
 *                                    Routine name changed from 
 *                                    cf_ttag_countmap to cf_countmap.
 *              02/26/04    1.7  bjg  Added option to select coordinates
 *                                    to use: RAW, FARF (default), final
 *              04/08/04    1.8  bjg  Remove unused function get_xy_columns
 *              04/27/04    1.9  bjg  Change author name of countfile
 *              08/25/04    1.10 bjg  Added lower and upper limits option
 *                                    for PHA. In force mode, issue a warning 
 *                                    if file already included. 
 *              03/02/05    1.11 wvd  Don't write AUTHOR keyword to header.
 *              08/24/07    1.12 bot  Changed nx and ny to long and called
 *				      as TLONG ; changed binx and biny to int
 *				      and called as TINT.
 *
 ****************************************************************************/

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>  
#include "calfuse.h"

#define ROOT_LENGTH 13
#define MAXCHARS 120
#define MASTER_CAL_FILE "master_calib_file.dat"
#define IMAGE_EXT 1
#define INDEX_EXT 2

static char CF_PRGM_ID[] = "cf_countmap";
static char CF_VER_NUM[] = "1.12";


static int
new_imagfile(char *newfile, char *detector, float effmjd, int binx, int biny)
{
    char  date[FLEN_CARD]={'\0'};
    char  *ttype[] = {"ROOTNAME", "OBSDATE", "COUNTS", "EXPTIME", "ADDDATE", "SCALE", "PHAMIN", "PHAMAX"};
    char  *tform[] = {"13A"     , "19A"    , "1J"    , "1E"     , "19A"    , "1E"   , "1I"    , "1I"    };
    char  *tunit[] = {" "       , " "      , "counts", "seconds", " "      , " "    , " "     , " "     };
    int   status=0, version=1, tref = 0, tfields = 8;
    long  axis[2];
    float exptime=0.;
    double nevents=0.;
    fitsfile *imagfits;

    FITS_create_file(&imagfits, newfile, &status);
    
    axis[0] = NXMAX/binx;
    axis[1] = NYMAX/biny;
    
    FITS_create_img(imagfits, FLOAT_IMG, 2, axis, &status);
    FITS_write_key(imagfits, TINT,    "SPECBINX",    &binx,
		   "Binsize in detector X coordinate", &status);
    FITS_write_key(imagfits, TINT,    "SPECBINY",    &biny,
		   "Binsize in detector Y coordinate", &status);
    FITS_write_key(imagfits, TSTRING, "CALFTYPE", "COUNT",
		   "Calibration file type", &status);
    FITS_write_key(imagfits, TINT,    "CALFVERS", &version,
		   "Calibration file version", &status);
    FITS_write_key(imagfits, TSTRING, "DETECTOR", detector,
		   "detector (1A, 1B, 2A, 2B)", &status);
    FITS_write_key(imagfits, TFLOAT,  "EFFMJD",   &effmjd,
		   "Date on which file should be applied (MJD)", &status);
    fits_get_system_time(date, &tref, &status);
    FITS_write_key(imagfits, TSTRING, "DATE",     date,
		   "file creation date (YYYY-MM-DDThh:mm:ss UTC)", &status);
    FITS_write_key(imagfits, TDOUBLE,   "NEVENTS",  &nevents,
		   "Accumulated events", &status);
    FITS_write_key(imagfits, TFLOAT,  "EXPTIME",  &exptime,
		   "Accumulated exposure time", &status);

    FITS_create_tbl(imagfits, BINARY_TBL, 0, tfields, ttype, tform, tunit,
		    "PROCESSED FILES", &status);

    FITS_close_file(imagfits, &status);
    return status;
}

static int
_check_index(fitsfile *imagfits, char *rootname, float *scale, int *phamin, int *phamax)
{
    int   status=0, anynull=0, hdutype;
    float *scal;
    int *phami, *phama;
    long  j, nrows, frow=1, felem=1;
    char  **filename, strnull[] = "          ";

    FITS_movabs_hdu(imagfits, INDEX_EXT, &hdutype, &status);
    FITS_read_key(imagfits, TLONG, "NAXIS2", &nrows, NULL, &status);

    filename = (char **) cf_malloc(nrows * sizeof(char *));
    filename[0] = (char *) cf_malloc(nrows*(ROOT_LENGTH+1)*sizeof(char));
    scal = (float *) cf_malloc(nrows * sizeof(float));
    phami = (int *) cf_malloc(nrows * sizeof(int));
    phama = (int *) cf_malloc(nrows * sizeof(int));
    for (j=1; j<nrows; j++)
	filename[j] = filename[j-1]+(ROOT_LENGTH+1);

    FITS_read_col(imagfits, TSTRING, 1, frow, felem, nrows, strnull,
		  filename, &anynull, &status);
    FITS_read_col(imagfits, TFLOAT,  6, frow, felem, nrows, strnull,
		  scal, &anynull, &status);
    FITS_read_col(imagfits, TINT,  7, frow, felem, nrows, strnull,
		  phami, &anynull, &status);
    FITS_read_col(imagfits, TINT,  8, frow, felem, nrows, strnull,
		  phama, &anynull, &status);

    for (j=0; j<nrows; j++) {
        if (strncmp(filename[j], rootname, strlen(rootname)) == 0) {
	    *scale = scal[j];
	    *phamin = phami[j];
	    *phamax = phama[j];
            free(filename[0]);
            free(filename);
	    free(scal);
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
_update_index(fitsfile *imagfits, char *rootname, char *obsdate,
	      double nevents, float exptime, float scale, int phamin, int phamax)
{
    char  *strptr[1], adddate[FLEN_CARD]={'\0'};
    int   status=0, hdutype, tref=0;
    long  nrows;

    FITS_movabs_hdu(imagfits, INDEX_EXT, &hdutype, &status);
    FITS_read_key(imagfits, TLONG, "NAXIS2", &nrows, NULL, &status);

    strptr[0] = rootname;
    FITS_write_col(imagfits, TSTRING, 1, nrows+1, 1, 1, strptr, &status);
    strptr[0] = obsdate;
    FITS_write_col(imagfits, TSTRING, 2, nrows+1, 1, 1, strptr, &status);
    FITS_write_col(imagfits, TDOUBLE, 3, nrows+1, 1, 1, &nevents, &status);
    FITS_write_col(imagfits, TFLOAT,  4, nrows+1, 1, 1, &exptime, &status);
    fits_get_system_time(adddate, &tref, &status);
    strptr[0] = adddate;
    fits_write_col(imagfits, TSTRING, 5, nrows+1, 1, 1, strptr, &status);
    FITS_write_col(imagfits, TFLOAT,  6, nrows+1, 1, 1, &scale, &status); 
    FITS_write_col(imagfits, TINT,    7, nrows+1, 1, 1, &phamin,  &status);
    FITS_write_col(imagfits, TINT,    8, nrows+1, 1, 1, &phamax,  &status);
    return 0;
}



static double
add_events(fitsfile *infits, char *xcolumn, char *ycolumn, float scale,
	   int phamin, int phamax, fitsfile *imagfits)
{
    int    binx, biny, status=0, hdutype, anynull=0;
    long   nx, ny, nevents, j;
    float  *xpos, *ypos, *weight, *image;
    char   *pha;
    double ncounts=0.;
    
    FITS_movabs_hdu(infits, 2, &hdutype, &status);
    nevents = cf_read_col(infits, TFLOAT, xcolumn,  (void **) &xpos);
    nevents = cf_read_col(infits, TFLOAT, ycolumn,  (void **) &ypos);
    nevents = cf_read_col(infits, TFLOAT, "WEIGHT", (void **) &weight);
    nevents = cf_read_col(infits, TBYTE, "PHA", (void **) &pha);

    FITS_movabs_hdu(imagfits, IMAGE_EXT, &hdutype, &status);
    FITS_read_key(imagfits, TLONG, "NAXIS1",   &nx,   NULL, &status);
    FITS_read_key(imagfits, TLONG, "NAXIS2",   &ny,   NULL, &status);
    FITS_read_key(imagfits, TINT, "SPECBINX", &binx, NULL, &status);
    FITS_read_key(imagfits, TINT, "SPECBINY", &biny, NULL, &status);
    
    image = (float *) cf_malloc(sizeof(float)*nx*ny);
    FITS_read_img(imagfits, TFLOAT, 1, nx*ny, NULL, image, &anynull, &status);

    for (j=0; j<nevents; j++) {
      if ((pha[j]>=phamin) && (pha[j]<=phamax)){
	int x, y;
	x = xpos[j];
	y = ypos[j];
	if ((x >= 0) && (x < nx*binx)  &&  (y >= 0) && (y < ny*biny)) {
	  image[(y/biny)*nx + (x/binx)] += weight[j]*scale;
	  ncounts += weight[j]*scale;
	}
      }
    }
    FITS_write_img(imagfits, TFLOAT, 1, nx*ny, image, &status);

    free(image);
    free(weight);
    free(ypos);
    free(xpos);
    return ncounts;
}

static double
subtract_events(fitsfile *infits, char *xcolumn, char *ycolumn,
		float scale, int phamin, int phamax, fitsfile *imagfits)
{
    int    binx, biny, status=0, hdutype, anynull=0;
    long   nx, ny, nevents, j;
    float  *xpos, *ypos, *weight, *image;
    double ncounts=0.;
    char    *pha;
    
    FITS_movabs_hdu(infits, 2, &hdutype, &status);
    nevents = cf_read_col(infits, TFLOAT, xcolumn,  (void **) &xpos);
    nevents = cf_read_col(infits, TFLOAT, ycolumn,  (void **) &ypos);
    nevents = cf_read_col(infits, TFLOAT, "WEIGHT", (void **) &weight);
    nevents = cf_read_col(infits, TBYTE, "PHA",    (void **) &pha);

    FITS_movabs_hdu(imagfits, IMAGE_EXT, &hdutype, &status);
    FITS_read_key(imagfits, TLONG, "NAXIS1",   &nx,   NULL, &status);
    FITS_read_key(imagfits, TLONG, "NAXIS2",   &ny,   NULL, &status);
    FITS_read_key(imagfits, TINT, "SPECBINX", &binx, NULL, &status);
    FITS_read_key(imagfits, TINT, "SPECBINY", &biny, NULL, &status);
    
    image = (float *) cf_malloc(sizeof(float)*nx*ny);
    FITS_read_img(imagfits, TFLOAT, 1, nx*ny, NULL, image, &anynull, &status);

    for (j=0; j<nevents; j++) {
      if ((pha[j]>=phamin) && (pha[j]<=phamax)){
	int x, y;
	x = xpos[j];
	y = ypos[j];
	if ((x >= 0) && (x < nx*binx)  &&  (y >= 0) && (y < ny*biny)) {
	  image[(y/biny)*nx + (x/binx)] -= weight[j]*scale;
	  ncounts += weight[j]*scale;
	}
      }
    }
    FITS_write_img(imagfits, TFLOAT, 1, nx*ny, image, &status);

    free(image);
    free(weight);
    free(ypos);
    free(xpos);
    return ncounts;
}


int main(int argc, char *argv[])
{
    char  *newfile = NULL, *detector = NULL, buffer[FLEN_CARD]={'\0'};
    int   subtract=0, force=0, newbinx=1, newbiny=1, status=0, optc;
    float effmjd = 0.,scale=1.,old_scale;
    int phamin=0, phamax=31;
    int old_phamin, old_phamax;
    
    char   xcolumn[FLEN_CARD]={'\0'}, ycolumn[FLEN_CARD]={'\0'};

    char opts[] = "hv:fsrcn:d:j:w:x:y:l:u:";
    char usage[] =
        "Usage:\n"
	"  cf_countmap [-hfsrc] [-v level] [-n filename] [-d segment]\n"
        "    [-w scale] [-j effmjd] [-x xbin] [-y ybin] \n"
        "    [-l phamin] [-u phamax] idffiles\n";
    char option[] =
	"Options:\n"
	"  -h:  this help message\n"
        "  -v:  verbosity level (=1; 0 is silent)\n"
	"  -f:  force processing\n"
	"  -s:  subtract files\n"
	"  -n:  image map file\n"
	"  -d:  image map segment (no default for new file)\n"
	"  -j:  effective MJD (=0.0)\n"
	"  -w:  weighting factor for image (=1.0)\n"
	"  -x:  X-bin size (=1)\n"
	"  -y:  Y-bin size (=1)\n"
        "  -r:  use raw XY coordinates instead of farf\n"
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
    if (argc <= optind) {
        cf_if_error("%s\nIncorrect number of program arguments", usage);
    }

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin processing");
    
    while (optind < argc) {

	char   imagfile[FLEN_CARD], rootname[FLEN_CARD]={'\0'};
	char   segment[FLEN_CARD] = {'\0'};
	int    hdutype;
	FILE   *fptr;
	fitsfile *infits, *imagfits;

	FITS_open_file(&infits, argv[optind], READONLY, &status);
	FITS_movabs_hdu(infits, 1, &hdutype, &status);

	if (newfile) {
	    strcpy(imagfile, newfile);

	    if ((fptr = fopen(imagfile, "r")))
		fclose(fptr);
	    else {
		if (!detector)
		    cf_if_error("No detector segment specified: "
				"1A, 1B, 2A, 2B");
		status = new_imagfile(imagfile, detector, effmjd,
				       newbinx, newbiny);
	    }
	    FITS_open_file(&imagfits, imagfile, READWRITE, &status);
	    FITS_movabs_hdu(imagfits, IMAGE_EXT, &hdutype, &status);
	    detector = buffer;
	    FITS_read_key(imagfits, TSTRING, "DETECTOR", detector, NULL,
			  &status);
	}
	else {
	    FITS_read_key(infits, TSTRING, "TCNT_MAP", imagfile, NULL,
			  &status);
	    detector = buffer;
	    FITS_read_key(infits, TSTRING, "DETECTOR", detector, NULL,
			  &status);

	    if ((fptr = fopen(cf_hist_file(imagfile), "r"))) {
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
			if (strcmp(filename, imagfile) == 0) {
			    effmjd = aftermjd;
			    break;
			}
		    }
		}
		fclose(master);

		status = new_imagfile(cf_hist_file(imagfile), detector,
				      effmjd, newbinx, newbiny);
	    }
	    FITS_open_file(&imagfits, cf_hist_file(imagfile), READWRITE,
			   &status);
	}

	strncpy(rootname, argv[optind], ROOT_LENGTH);
	FITS_read_key(infits, TSTRING,  "DETECTOR", segment, NULL,
			  &status);
	if (strcmp(detector, segment) != 0) {
	    cf_if_warning("%s segment is not the same as %s.",
			  argv[optind], imagfile);
	}
	else if (subtract) {
	    int  index;
	    if ((index = _check_index(imagfits, rootname, &scale, &phamin, &phamax))) {
		
		int    hdutype;
		float  exptime=0., acctime=0.;
		double nevents, accevnt=0.;

		FITS_read_key(infits, TFLOAT,  "EXPTIME", &exptime, NULL,
			      &status);

		
		nevents = subtract_events(infits, xcolumn, ycolumn, scale, phamin, phamax, imagfits);

		FITS_movabs_hdu(imagfits, IMAGE_EXT, &hdutype, &status);
		FITS_read_key(imagfits, TDOUBLE, "NEVENTS",  &accevnt,
			      NULL, &status);
		FITS_read_key(imagfits, TFLOAT, "EXPTIME", &acctime,
			      NULL, &status);
		accevnt -= nevents;
		acctime -= exptime;
		FITS_update_key(imagfits, TDOUBLE,  "NEVENTS", &accevnt,
				NULL, &status);
		FITS_update_key(imagfits, TFLOAT, "EXPTIME", &acctime,
				NULL, &status);

		FITS_movabs_hdu(imagfits, INDEX_EXT, &hdutype, &status);
		FITS_delete_rows(imagfits, index, 1, &status);
		cf_verbose(1, "%s is subtracted from %s",
			   argv[optind], imagfile);
	    } else {
	      cf_if_warning("%s is not in %s.", argv[optind], imagfile);
	    }
	} else {
	    if (!force && _check_index(imagfits, rootname, &old_scale, &old_phamin, &old_phamax)) {
	      cf_if_warning("%s is already in %s with scale = %f, phamin = %d, phamax=%d.",
			    argv[optind], imagfile, old_scale, old_phamin, old_phamax);
	    } else {
		char   obsdate[FLEN_CARD]={'\0'}, dateobs[FLEN_CARD]={'\0'};
		char   timeobs[FLEN_CARD]={'\0'};
		
		float  exptime=0., acctime=0.;
		double nevents, accevnt=0;

		FITS_movabs_hdu(imagfits, IMAGE_EXT, &hdutype, &status);
		FITS_read_key(infits, TSTRING, "DATEOBS", dateobs, NULL,
			      &status);
		FITS_read_key(infits, TSTRING, "TIMEOBS", timeobs, NULL,
			      &status);
		obsdate[0] = '\0';
		sprintf(obsdate, "%s%s%s", dateobs, "T", timeobs);
		FITS_read_key(infits, TFLOAT,  "EXPTIME", &exptime, NULL,
			      &status);

		if ((force) && _check_index(imagfits, rootname, &old_scale, &old_phamin, &old_phamax)) 
		  cf_if_warning("FORCE MODE - %s is already in %s with scale = %f, phamin = %d, phamax=%d.",
				argv[optind], imagfile, old_scale, old_phamin, old_phamax);

		nevents = add_events(infits, xcolumn, ycolumn, scale, phamin, phamax, imagfits);

		FITS_movabs_hdu(imagfits, IMAGE_EXT, &hdutype, &status);
		FITS_read_key(imagfits, TDOUBLE, "NEVENTS",  &accevnt,
			      NULL, &status);
		FITS_read_key(imagfits, TFLOAT, "EXPTIME", &acctime,
			      NULL, &status);
		accevnt += nevents;
		acctime += exptime;
		FITS_update_key(imagfits, TDOUBLE, "NEVENTS", &accevnt,
				NULL, &status);
		FITS_update_key(imagfits, TFLOAT, "EXPTIME", &acctime,
				NULL, &status);

		status = _update_index(imagfits, rootname, obsdate, nevents,
				       exptime, scale, phamin, phamax);
		cf_verbose(1, "%s is added to %s", argv[optind], imagfile);
	    }
	}
	FITS_close_file(imagfits, &status);
	FITS_close_file(infits, &status);
	optind ++;
    }
    
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
    return 0;
}
