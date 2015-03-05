/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Synopsis:	cf_arith infile1 op infile2 outfile
 *
 * Description: This is a generic utility program that performs arithmetic 
 *		operations on two FUSE FITS files.  The files must
 *		contain extracted FUSE spectra.
 *		Input files are assumed to have 6 columns in a binary
 *		table in the first extension (wavelength, flux, error, quality,
 *		counts, counts_error).  Files with 4 columns (no counts or
 *		counts_error) are also handled for compatibility with
 *		pipeline output written prior to 12/1/99.
 *
 *		(The next paragraph is included for completeness, but is
 *		not expected to be applicable to most users).
 *		If infile1 and infile2 have the same number of extensions,
 *		the specified operation is performed for each extension
 *		(outfile_ext_i = infile1_ext_i op infile2_ext_i).
 *		If infile2 has a single extension, the operation is performed
 *		for each extension in infile1 with the extension in infile2
 *		(outfile_ext_i = infile1_ext_i op infile2_ext_1).
 *
 *		infile1 and infile2 must have the same number of elements 
 *		in each binary table (ie. binned by the same factor).
 *		infile2 may be a numeric string instead of a file.
 *		The header in outfile will be the same as for infile1,
 *		plus some HISTORY lines.
 *		No check is performed that the wavelengths in infile2 match
 *		those in infile1.
 *		It is assumed that error arrays contain errors, not variances.
 *
 *		If infile2 is the string "EXPTIME" (without the quotes)
 *		operand 2 is set to the exposure time in the header of infile1.
 *
 *		valid choices for the operation "op" are:
 *		+, -, *, /, avg_t, avg_s, repl_d, repl_e, boxcar, delta, 
 *		min, max, shift, bin, deriv, addnoise
 *
 *		avg_t = average weighted by exposure time
 *		avg_s = average weighted by sigma (errors)
 *		  avg_t, avg_s valid only if infile2 is a fits file.
 *		repl_d  = replace data in operand1 with data from operand2
 *		repl_e  = replace errors in operand1 with errors from operand2
 *		boxcar  = smooth data in operand1 with moving boxcar
 *			 of width operand2.  If operand2 is a FITS file,
 *			 the "FLUX" column is used as the width.  This
 *			 permits application of variable smoothing.
 *		delta = adds delta functions to operand 1 every 100 pixels;
 *			height is equal to operand 2 (must be a scalar).
 *			Used only for creating test spectra.
 *		shift = shift data by the specified number of pixels.  Operand2
 *			must be an integer.  A positive shift moves the data 
 *			to higher pixel numbers.
 *		bin   = bin data by the specified number of pixels.  Operand2
 *			must be an integer.  The flux column is averaged, the
 *			counts column is summed.
 *		deriv = take derivative of input: out[i] = in[i+1] - in[i]
 *			operand 2 is ignored for this!
 *		addnoise = a uniformly distributed random number in the range
 *			+/-operand2 is added pixel-by-pixel to the flux 
 *			Used only for creating test spectra.
 *
 *		If op = "+" and the second operand is a file (not a scalar), 
 *			or if op = "avg_t" or if op = "avg_s", the exposure
 *			time in the output is the sum of the exposure times
 *			in the input files.
 *
 *		For op1 min op2 out:  op1 is replaced by op2 if op1 < op2.
 *		  (op2 sets a min value for the flux)
 *		For op1 max op2 out:  op1 is replaced by op2 if op1 > op2.
 *		  (op2 sets a max value for the flux)
 *
 *		The same operations are performed on the counts/countserr
 *		columns as on the flux/fluxerr columns.  In most cases this
 *		leads to reasonable results.  However, in a few cases, such
 *		as addition or subtraction of a scalar, the result is likely
 *		to be meaningless for either the flux or counts columns.
 *		Exceptions:  MIN, MAX, REPL_D, REPL_E, DELTA, only affect
 *		the flux and fluxerr columns.  The counts and countserror
 *		columns are unaffected.
 *
 *		Note that avg_t and avg_s also average the counts; in some
 *		cases the user may have really wanted the counts to be summed!
 *
 *		NOTE:  Must enclose * in quotes (ie. "*") or the shell
 *			will expand it into filenames in the pwd!
 *
 *
 * History:     08/15/99        jwk     started work
 *		10/25/99 v1.1	jwk	avg_t,_s sum exp times in output
 *		11/30/99	jwk	handle either 4 or 6-column files
 *		01/06/00 v1.2	jwk	add bin option
 *		05/01/00 v1.3	jwk	handle byte or short quality flags
 *		09/20/00 v1.4   jwk	add derivative operation
 *		10/11/00 v1.5   jwk     update FILENAME keyword in output file
 *		12/11/00 v1.6   jwk	improve numeric accuracy of error prop.
 *		01/23/01 v1.7	jwk	smooth quality rather than sum when
 *					doing BOXCAR, to avoid overflow of
 *					short datasize.
 *		10/10/01 v1.8	jwk	add handling of eff area files
 *              05/30/03 v1.9   jwk     change timeval from long to time_t; put
 *                                      ifdef SOLARIS/endif around ieee_retro
 *                                      spec field width for final timestamp
 *                                      (make consistent with calfuse version)
 *		02/04/04 v2.0   jwk     implement addnoise option
 *		02/25/04 v2.1   jwk     check for zero errors in avg_s code
 *		03/06/04 v3.0   jwk     add support for calfuse 3.0 files
 *              04/05/04 v3.1   bjg     Declarations for cf_wrspec7 and
 *                                      cf_wrspec_cf2
 *                                      Include ctype.h
 *                                      Correct formats to match argument
 *                                      types in printf
 *                                      Remove unused variable long npix
 *                                      Casting to char * arg of setstate
 *                                      Remove local definitions of
 *                                      initstate() and setstate()
 *                                      Correct use of cf_wrspec7
 *              08/26/05 v3.2   wvd     Change name from cf_arith3 to cf_arith.
 *					If no arguments, don't return error.
 *					Just print message and exit.\
 *              08/27/07 v3.3   bot 	Added L to long l.399 and 537
 *
 ****************************************************************************/

#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calfuse.h"
#include <ctype.h>

#define	PLUS	1
#define	MINUS	2
#define	MULT	3
#define	DIV	4
#define	AVG_T	5
#define	AVG_S	6
#define REPL_D	7
#define REPL_E	8
#define	BOXCAR	9
#define DELTA	10
#define MIN	11
#define	MAX	12
#define SHIFT	13
#define BIN	14
#define DERIV	15
#define ADDNOISE 16

/* arbitrary value for error when no data is present */
#define DEF_ERROR 1.e-3

#define MAX_RAND 0x7fffffff

static char CF_PRGM_ID[] = "cf_arith";
static char CF_VER_NUM[] = "3.3";


void cf_wrspec7(fitsfile *outfits, long npts, float *wave, float *flux, 
		       float *error, long *counts, float *weights, 
		       float *bkgd,short *pothole);
void cf_wrspec_cf2(fitsfile *outfits, int npts, float *wave,
			  float *spec, float *errs, short *qual,
			  float *counts, float *cntserr, int ncol, 
			  int areaflag);




int main(int argc, char *argv[])
{
    char  *timestr;
    time_t  timeval;
    char  infile1[80], infile2[80], outfile[80];
    char  comment[FLEN_CARD];
    char  expt_comment[FLEN_CARD];
    char  errstr[80];
    char  datatype[20];
    int   areaflag;
    int	  nbin;
    int   next_1, next_2;
    int   status=0, anynull=0, intnull=0;
    int	  hdutype_1, hdutype_2, hdutype_o;
    int	  hdunum_1, hdunum_2, hdunum_o;
    long  npix_1, npix_2, npix_o;
    long  naxis2_1, naxis2_2;
    int   fluxcol, wavecol, errorcol, qualcol, cntscol, cntserrcol;
    int   bkgdcol, wgtscol;
    int	  tfields_1, tfields_2;
    int	  opcode;
    int   cf_version;
    long  i, j, k, k1, k2;
    int   arg2_scalar;
    int	  p_shift;
    float arg2_val;
    float swave, sflux, sferr, scnts, scerr, spix;
    float sbkgd, swgts;
    long   lqual, slcts;
    short  squal;
    short  *qual1, *qual2, *qual_o;
    float *wave1, *flux1, *flxerr1, *cnts1, *cntserr1;
    float *flux2, *flxerr2, *cnts2, *cntserr2;
    float *flux_o, *flxerr_o, *cnts_o, *cntserr_o;
    float *wave_o;
    float *wgts1, *wgts2, *wgts_o;
    float *bkgd1, *bkgd2, *bkgd_o;
    long  *lcts1, *lcts2, *lcts_o;
    char  tformstr[FLEN_CARD];
    char  key_str[FLEN_CARD];
    char  qual_lbl[FLEN_CARD];
    float f1, f2, fn;
    float tmp1, tmp2, tmp3, tmp4;
    float exptime_1, exptime_2, exptime_o;
    fitsfile *infits_1, *infits_2, *outfits;
    int	  chkop();
    int   check_num();
    void  padstr();
    /* following items used by random number generator */
    clock_t    clock();        /* system clock function */
    int     ns;             /* used by random # generator */
    unsigned long   seed;   /* random # generator seed */
 /*   char    *initstate(), *setstate(); */
    double  rand1();

/* Initialize random number generator state array. */
    static unsigned long state1[32]  =  {       3,       0x9a319039,
          0x32d9c024,  0x9b663182,  0x5da1f342,       0x7449e56b,
          0xbeb1dbb0,  0xab5c5918,  0x946554fd,       0x8c2e680f,
          0xeb3d799f,  0xb11ee0b7,  0x2d436b86,       0xda672e2a,
          0x1588ca88,  0xe369735d,  0x904f35f7,       0xd7158fd6,
          0x6fa6f051,  0x616e6b96,  0xac94efdc,       0xde3b81e0,
          0xdf0a6fb5,  0xf103bc02,  0x48f340fb,       0x36413f93,
          0xc622c298,  0xf5a42ab8,  0x8a88d77b,       0xf5ad9d0e,
          0x8999220b, 0x27fb47b9      };


    /* perform an initial call to clock() - clock() returns the elapsed time
     * since the first call, so in order to use it for a random # generator
     * seed later we have to get it started here.  Subsequent calls follow
     * opening of the first file, so that should provide some run-to-run
     * variations.
     */
    seed = clock();

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);


    /* check command line */
    if (argc != 5) {
	printf("Usage: cf_arith infile1 op infile2 outfile\n");
	return 1;
    }

    strcpy (infile1, argv[1]);
    strcpy (infile2, argv[3]);
    strcpy (outfile, argv[4]);

    if ((opcode = chkop (argv[2])) == 0)
	{
	sprintf (comment,"Undefined operation: %s\n", argv[2]);
	cf_if_error (comment);
	}

   /* Open input file(s) */

   FITS_open_file(&infits_1, infile1, READONLY, &status);

   /* get number of extensions in operand1 */
   FITS_get_num_hdus (infits_1, &next_1, &status);
   next_1--;

   /* get exposure time in case we need it later. */
   ffgky(infits_1, TFLOAT, "EXPTIME", &exptime_1, expt_comment, &status);
   if (status)
	{
	fprintf (stderr, "Keyword EXPTIME not found; set exptime_1=1.0\n");
	exptime_1 = 1.0;
	status = 0;
	}
   exptime_o = exptime_1;

   /* check nature of second operand:  is it a number? */
   /* can have leading sign, and may contain E,e,D,d and additional sign */
   arg2_scalar = check_num (argv[3]);
   if (arg2_scalar)
	arg2_val = atof (argv[3]);

   /* Also allow the possibility of having argv[3] = "EXPTIME": */
   if (!strcmp (argv[3], "EXPTIME"))
	{
	arg2_scalar = TRUE;
	arg2_val = exptime_1;
	}

   /* check for valid operand if operation is SHIFT: */
   p_shift = 0;
   if (opcode == SHIFT) 
	{
	if (!arg2_scalar)
	   cf_if_error ("SHIFT requires scalar argument!\n");

	/* convert argument to an integer */
	p_shift = (int) strtol (argv[3], (char **)NULL, 10);
	}

   /* check for valid operand if operation is BIN: */
   nbin = 1;
   if (opcode == BIN) 
	{
	if (!arg2_scalar)
	   cf_if_error ("BIN requires scalar argument!\n");

	/* convert argument to an integer */
	nbin = (int) strtol (argv[3], (char **)NULL, 10);
	if (nbin < 1)
	   {
	   sprintf (errstr, "bin factor = %d; must be >0!\n", nbin);
	   cf_if_error (errstr);
	   }
	}

   /* Check for valid operand if operation is ADDNOISE: 
    * I suppose we could relax this and take a pixel-by-pixel amplitude
    * from a spectrum, but let's keep it simple for now.
    * Also initialize random number generator.
    */
   if (opcode == ADDNOISE) 
	{
	if (!arg2_scalar)
	   cf_if_error ("ADDNOISE requires scalar argument!\n");

	/* initialize random number generator */
	/* the local clock function doesn't change rapidly enough to be
	 * useful, even though it is nominally at microsecond resolution */
        /* seed  =  2 * (clock() & 0x0ffff) + 1; */
	time (&timeval);
	seed = 2 * timeval + 1;
	
        ns  =  128;
        initstate (seed, (char *) state1, ns);
        setstate ((char *)state1);
	}

   /* bypass operand2 if op= derivative; we ignore op2 in this case*/
   if (opcode == DERIV)
	arg2_scalar = TRUE;

   if (!arg2_scalar)
	{
	FITS_open_file (&infits_2, infile2, READONLY, &status);
	FITS_get_num_hdus (infits_2, &next_2, &status);
	next_2--;

	if ((next_2 != next_1) && (next_2 != 1))
	   {
	   fprintf (stderr, "file %s has %d extensions.\n", infile1, next_1);
	   fprintf (stderr, "file %s has %d extensions.\n", infile2, next_2);
	   cf_if_error ("Illegal combination, aborting!\n");
	   }
	ffgky(infits_2, TFLOAT, "EXPTIME", &exptime_2, expt_comment, &status);
	if (status)
	  {
	  fprintf (stderr, "Keyword EXPTIME not found; set exptime_2=1.0\n");
	  exptime_2 = 1.0;
	  status = 0;
	  }
	}

   /* open output file */

   remove(outfile);	/*In case the old file name exists*/
   FITS_create_file (&outfits, outfile, &status);

   /* Copy header from first operand.  The output file will have the same
    * number of extensions as operand1.
    */
   /* cf_cp_hdr(infits_1, 1, outfits, 1);  */
   FITS_copy_header (infits_1, outfits, &status);

    /* Loop over extensions */
    for ( i = 1; i <= next_1; i++ ) 
	{
	hdunum_1=i+1;

	FITS_movabs_hdu(infits_1, hdunum_1, &hdutype_1, &status);
	/* for Calfuse 2 and before, array length is specified by NAXIS2.
	 * For Calfuse 3, need to get it from TFORM1.
	 */
	FITS_read_key(infits_1, TLONG, "NAXIS2", &naxis2_1, NULL, &status);
	if (naxis2_1 > 1L)
	   npix_1 = naxis2_1;
	else
	   {
	   FITS_read_key(infits_1, TSTRING, "TFORM1", tformstr, NULL, &status);
	   sscanf (tformstr, "%ldE", &npix_1);
	   }


	/* need to figure out what kind of data we have.
	 * tfields=4 implies very old calfuse files, or effective area files;
	 * tfields=6 implies calfuse 1.8 or 2.x
	 * tfields=7 implies calfuse 3.x
	 * Each has different types of data columns.
	 */
	FITS_read_key(infits_1, TINT, "TFIELDS", &tfields_1, NULL, &status);

	if ((tfields_1 == 4) || (tfields_1 == 6))
	   cf_version = 2;
	else if (tfields_1 == 7)
	   cf_version = 3;
	else
	   {
	   sprintf (errstr, "Expecting 4,6 or 7 cols, found %d\n", tfields_1);
	   cf_if_error (errstr);
	   }

        /* Allocate data arrays */
	wave1 = (float *) cf_malloc(sizeof(float) * npix_1);
	flux1 = (float *) cf_malloc(sizeof(float) * npix_1);
	flxerr1 = (float *) cf_malloc(sizeof(float) * npix_1);
	qual1 = (short *) cf_malloc(sizeof(short) * npix_1);

	if (cf_version == 2)
	   {
	   cnts1 = (float *) cf_malloc(sizeof(float) * npix_1);
	   cntserr1 = (float *) cf_malloc(sizeof(float) * npix_1);
	   }
	else
	   {
	   lcts1 = (long *) cf_malloc(sizeof(long) * npix_1);
	   bkgd1 = (float *) cf_malloc(sizeof(float) * npix_1);
	   wgts1 = (float *) cf_malloc(sizeof(float) * npix_1);
	   }

	/* Get the column numbers from file 1 */
	FITS_read_key(infits_1, TSTRING, "TTYPE2", datatype, NULL, &status);
	if (!strcmp(datatype, "FLUX"))
	   {
	   areaflag = FALSE;
	   FITS_get_colnum(infits_1, TRUE, "FLUX", &fluxcol, &status);
	   }
	else 
	   {
	   areaflag = TRUE;
	   FITS_get_colnum(infits_1, TRUE, "AREA", &fluxcol, &status);
	   }

	FITS_get_colnum(infits_1, TRUE, "WAVE", &wavecol, &status);
	FITS_get_colnum(infits_1, TRUE, "ERROR", &errorcol, &status);

	/* don't know at present if calfuse 3.x will use "QUALITY" or
	 * "POTHOLE", so check both.
	 */
	if (cf_version == 2)
	   FITS_get_colnum(infits_1, TRUE, "QUALITY", &qualcol, &status);
	else
	   {
	   qualcol = 0;
	   for (j=1; j<=tfields_1; j++)
	      {
	      sprintf (key_str, "TTYPE%1ld", j);
	      FITS_read_key(infits_1, TSTRING, key_str,datatype, NULL, &status);
	      if (!strcmp (datatype, "QUALITY") || !strcmp(datatype,"POTHOLE"))
		{
		strcpy (qual_lbl, datatype);
		qualcol = j;
		}
	      }
	   if (qualcol == 0)
	      cf_if_error ("QUALITY column not found!");
	   }

	if (tfields_1 == 6)
	   {
	   FITS_get_colnum(infits_1, TRUE, "COUNTS", &cntscol, &status);
	   FITS_get_colnum(infits_1, TRUE, "CNTSERR", &cntserrcol, &status);
	   }

	if (tfields_1 == 7)
	   {
	   FITS_get_colnum(infits_1, TRUE, "COUNTS", &cntscol, &status);
	   FITS_get_colnum(infits_1, TRUE, "WEIGHTS", &wgtscol, &status);
	   FITS_get_colnum(infits_1, TRUE, "BKGD", &bkgdcol, &status);
	   }

	FITS_read_col(infits_1, TFLOAT, fluxcol, 1L, 1L, npix_1,
		&intnull, flux1, &anynull, &status);
	FITS_read_col(infits_1, TFLOAT, wavecol, 1L, 1L, npix_1,
		&intnull, wave1, &anynull, &status);
	FITS_read_col(infits_1, TFLOAT, errorcol, 1L, 1L, npix_1,
		&intnull, flxerr1, &anynull, &status);
	FITS_read_col(infits_1, TSHORT, qualcol, 1L, 1L, npix_1,
		&intnull, qual1, &anynull, &status);

	if (tfields_1 == 6)
	   {
	   FITS_read_col(infits_1, TFLOAT, cntscol, 1L, 1L, npix_1,
		&intnull, cnts1, &anynull, &status);
	   FITS_read_col(infits_1, TFLOAT, cntserrcol, 1L, 1L, npix_1,
		&intnull, cntserr1, &anynull, &status);
	   }
	else if (tfields_1 == 7)
	   {
	   FITS_read_col(infits_1, TLONG, cntscol, 1L, 1L, npix_1,
		&intnull, lcts1, &anynull, &status);
	   FITS_read_col(infits_1, TFLOAT, wgtscol, 1L, 1L, npix_1,
		&intnull, wgts1, &anynull, &status);
	   FITS_read_col(infits_1, TFLOAT, bkgdcol, 1L, 1L, npix_1,
		&intnull, bkgd1, &anynull, &status);
	   }
	else
	   {
	   for (j=0; j<npix_1; j++)
		{
		cnts1[j] = 0.;
		cntserr1[j] = 0.001;
		}
	   }

	if (!arg2_scalar)
	   {
	   /* for file2, stop at first extension if that is all there is */
	   if (i <= next_2)
		{
		hdunum_2 = i + 1;
		FITS_movabs_hdu(infits_2, hdunum_2, &hdutype_2, &status);
		FITS_read_key(infits_2, TLONG,"NAXIS2",&naxis2_2,NULL,&status);
		if (naxis2_2 > 1L)
		   npix_2 = naxis2_2;
		else
		   {
		   FITS_read_key(infits_2, TSTRING, "TFORM1", tformstr, NULL, 
			&status);
		   sscanf (tformstr, "%ldE", &npix_2);
		   }

		FITS_read_key(infits_2, TINT, "TFIELDS", &tfields_2, NULL, 
				&status);

		if (tfields_2 != tfields_1)
		   {
		   sprintf (errstr, "File1 has %d cols; File 2 has %d cols!\n", 
			tfields_1, tfields_2);
		   cf_if_error (errstr);
		   }

		if (npix_2 != npix_1)
		  {
		  fprintf (stderr, "npix1 = %ld, npix2 = %ld\n", npix_1, npix_2);
		  cf_if_error ("Aborting");
		  }

        	/* Allocate data arrays */
		flux2 = (float *) cf_malloc(sizeof(float) * npix_2);
		flxerr2 = (float *) cf_malloc(sizeof(float) * npix_2);
		qual2 = (short *) cf_malloc(sizeof(short) * npix_2);
		if (cf_version == 2)
		  {
		  cnts2 = (float *) cf_malloc(sizeof(float) * npix_2);
		  cntserr2 = (float *) cf_malloc(sizeof(float) * npix_2);
		  }
		else
		  {
		  lcts2 = (long *) cf_malloc(sizeof(long) * npix_2);
		  wgts2 = (float *) cf_malloc(sizeof(float) * npix_2);
		  bkgd2 = (float *) cf_malloc(sizeof(float) * npix_2);
		  }

		/* Read the data from file 2 */

		if (areaflag)
		   FITS_get_colnum(infits_2, TRUE, "AREA", &fluxcol, &status);
		else
		   FITS_get_colnum(infits_2, TRUE, "FLUX", &fluxcol, &status);

		FITS_get_colnum(infits_2, TRUE, "ERROR", &errorcol, &status);

		/* don't know at present if calfuse 3.x will use "QUALITY" or
		 * "POTHOLE", so check both.
		 * (there may also be a mix of files present!)
		 */
		if (cf_version == 2)
		   FITS_get_colnum(infits_2, TRUE, "QUALITY", &qualcol,&status);
		else
		   {
		   qualcol = 0;
		   for (j=1; j<=tfields_2; j++)
		      {
		      sprintf (key_str, "TTYPE%1ld", j);
		      FITS_read_key(infits_2, TSTRING, key_str,datatype, NULL, 
				&status);
		      if (!strcmp (datatype, "QUALITY") || 
			  !strcmp(datatype,"POTHOLE"))
			{
			strcpy (qual_lbl, datatype);
			qualcol = j;
			}
		      }
		   if (qualcol == 0)
		      cf_if_error ("QUALITY column not found!");
		   }

		if (tfields_2 == 6)
		   {
		   FITS_get_colnum(infits_2, TRUE, "COUNTS", &cntscol, &status);
		   FITS_get_colnum(infits_2, TRUE, "CNTSERR", &cntserrcol, 
				&status);
		   }

		if (tfields_2 == 7)
		   {
		   FITS_get_colnum(infits_2,TRUE, "COUNTS", &cntscol, &status);
		   FITS_get_colnum(infits_2,TRUE, "WEIGHTS", &wgtscol, &status);
		   FITS_get_colnum(infits_2,TRUE, "BKGD", &bkgdcol, &status);
		   }

		FITS_read_col(infits_2, TFLOAT, fluxcol, 1L, 1L, npix_2,
			&intnull, flux2, &anynull, &status);
		FITS_read_col(infits_2, TFLOAT, errorcol, 1L, 1L, npix_2,
			&intnull, flxerr2, &anynull, &status);
	 	FITS_read_col(infits_2, TSHORT, qualcol, 1L, 1L, npix_2,
			&intnull, qual2, &anynull, &status);

		if (tfields_2 == 6)
		   {
		   FITS_read_col(infits_2, TFLOAT, cntscol, 1L, 1L, npix_2,
			&intnull, cnts2, &anynull, &status);
		   FITS_read_col(infits_2, TFLOAT, cntserrcol, 1L, 1L, npix_2,
			&intnull, cntserr2, &anynull, &status);
		   }
		else if (tfields_1 == 7)
		   {
		   FITS_read_col(infits_2, TLONG, cntscol, 1L, 1L, npix_2,
			&intnull, lcts2, &anynull, &status);
		   FITS_read_col(infits_2, TFLOAT, wgtscol, 1L, 1L, npix_2,
			&intnull, wgts2, &anynull, &status);
		   FITS_read_col(infits_2, TFLOAT, bkgdcol, 1L, 1L, npix_2,
			&intnull, bkgd2, &anynull, &status);
		   }
		else
		   {
		   for (j=0; j<npix_2; j++)
			{
			cnts2[j] = 0.;
			cntserr2[j] = 0.001;
			}
		   }

		} /* end of if (i<=next_2) */
	   } /* end of if(!arg2_scalar) */

	if ((opcode == BIN) && (nbin > 1))
	   {
	   if (nbin > npix_1)
	      {
	      sprintf (errstr, "bin factor = %d; must be < npix!\n", nbin);
	      cf_if_error (errstr);
	      }
	   npix_o = npix_1 / nbin;
	   }
	else
	   npix_o = npix_1;

	/* make space for output arrays */
	flux_o = (float *) cf_malloc(sizeof(float) * npix_o);
	flxerr_o = (float *) cf_malloc(sizeof(float) * npix_o);
	qual_o = (short *) cf_malloc(sizeof(short) * npix_o);

	if (cf_version == 2)
	  {
	  cnts_o = (float *) cf_malloc(sizeof(float) * npix_o);
	  cntserr_o = (float *) cf_malloc(sizeof(float) * npix_o);
	  }
	else
	  {
	  lcts_o = (long *) cf_malloc(sizeof(long) * npix_o);
	  wgts_o = (float *) cf_malloc(sizeof(float) * npix_o);
	  bkgd_o = (float *) cf_malloc(sizeof(float) * npix_o);
	  }

	if (opcode == BIN)
	   wave_o =  (float *) cf_malloc(sizeof(float) * npix_o);
	else
	   wave_o = wave1;

	/********* perform arithmetic operations here *********/
	switch (opcode)
	   {
	   case PLUS:
	     if (arg2_scalar)
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] + arg2_val;
		   flxerr_o[j] = flxerr1[j];
		   cnts_o[j] = cnts1[j]; 
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] + arg2_val;
		   flxerr_o[j] = flxerr1[j];
		   lcts_o[j] = lcts1[j]; 
		   qual_o[j] = qual1[j];
		   bkgd_o[j] = bkgd1[j];
		   wgts_o[j] = wgts1[j];
		   }
		  }
		} /* end if(arg2_scalar) */
	     else
		{
		exptime_o = exptime_1 + exptime_2;
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] + flux2[j];
		   tmp1 = flxerr1[j]*flxerr1[j];
		   tmp2 = flxerr2[j]*flxerr2[j];
		   flxerr_o[j] = sqrt (tmp1 + tmp2);
		   cnts_o[j] = cnts1[j] + cnts2[j];
		   tmp1 = cntserr1[j]*cntserr1[j];
		   tmp2 = cntserr2[j]*cntserr2[j];
		   cntserr_o[j] = sqrt (tmp1 + tmp2);
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] + flux2[j];
		   tmp1 = flxerr1[j]*flxerr1[j];
		   tmp2 = flxerr2[j]*flxerr2[j];
		   flxerr_o[j] = sqrt (tmp1 + tmp2);
		   lcts_o[j] = lcts1[j] + lcts2[j];
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = (qual1[j] + qual2[j])/2;
		   bkgd_o[j] = bkgd1[j] + bkgd2[j];
		   wgts_o[j] = wgts1[j] + wgts2[j];
		   }
		  }
		}
	     break;

	   case MINUS:
	      if (arg2_scalar)
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] - arg2_val;
		   flxerr_o[j] = flxerr1[j];
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] - arg2_val;
		   flxerr_o[j] = flxerr1[j];
		   lcts_o[j] = lcts1[j];
		   qual_o[j] = qual1[j];
		   bkgd_o[j] = bkgd1[j];
		   wgts_o[j] = wgts1[j];
		   }
		  }
		}
	      else
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] - flux2[j];
		   tmp1 = flxerr1[j]*flxerr1[j];
		   tmp2 = flxerr2[j]*flxerr2[j];
		   flxerr_o[j] = sqrt (tmp1 + tmp2);
		   cnts_o[j] = cnts1[j] - cnts2[j];
		   tmp1 = cntserr1[j]*cntserr1[j];
		   tmp2 = cntserr2[j]*cntserr2[j];
		   cntserr_o[j] = sqrt (tmp1 + tmp2);
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] - flux2[j];
		   tmp1 = flxerr1[j]*flxerr1[j];
		   tmp2 = flxerr2[j]*flxerr2[j];
		   flxerr_o[j] = sqrt (tmp1 + tmp2);
		   lcts_o[j] = lcts1[j] - lcts2[j];
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = (qual1[j] + qual2[j])/2;
		   bkgd_o[j] = bkgd1[j] - bkgd2[j];
		   wgts_o[j] = wgts1[j] - wgts2[j];
		   }
		  }
		}
		break;

	   case MULT:
	      if (arg2_scalar)
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] * arg2_val;
		   flxerr_o[j] = flxerr1[j] * arg2_val;
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] * arg2_val;
		   flxerr_o[j] = flxerr1[j] * arg2_val;
		   lcts_o[j] = lcts1[j];
		   qual_o[j] = qual1[j];
		   bkgd_o[j] = bkgd1[j];
		   wgts_o[j] = wgts1[j];
		   }
		  } /* end cf_version == 3 */
		} /* end if arg2_scalar */
	      else
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] * flux2[j];
		   if (flux1[j] == 0.)
			{
			flxerr_o[j] = fabs (flux2[j] * flxerr1[j]);
			}
		   else if (flux2[j] == 0.)
			{
			flxerr_o[j] = fabs (flux1[j] * flxerr2[j]);
			}
		   else
			{
			tmp1 = flxerr1[j] / flux1[j];
			tmp2 = flxerr2[j] / flux2[j];
			tmp3 = tmp1 * tmp1 + tmp2 * tmp2;
			flxerr_o[j] = fabs (flux_o[j]) * sqrt (tmp3);
			}
		   cnts_o[j] = cnts1[j] * cnts2[j];
		   tmp1 = cnts2[j] * cnts2[j] * cntserr1[j]*cntserr1[j];
		   tmp2 = cnts1[j] * cnts1[j] * cntserr2[j]*cntserr2[j];
		   cntserr_o[j] = sqrt (tmp1 + tmp2);
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }  /* end loop over j */
		  }   /* end cf_version == 2 */
		else
		  {   /* start cf_version == 3 */
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] * flux2[j];
		   if (flux1[j] == 0.)
			{
			flxerr_o[j] = fabs (flux2[j] * flxerr1[j]);
			}
		   else if (flux2[j] == 0.)
			{
			flxerr_o[j] = fabs (flux1[j] * flxerr2[j]);
			}
		   else
			{
			tmp1 = flxerr1[j] / flux1[j];
			tmp2 = flxerr2[j] / flux2[j];
			tmp3 = tmp1 * tmp1 + tmp2 * tmp2;
			flxerr_o[j] = fabs (flux_o[j]) * sqrt (tmp3);
			}
		   lcts_o[j] = lcts1[j] * lcts2[j];
		   wgts_o[j] = wgts1[j] * wgts2[j];
		   bkgd_o[j] = bkgd1[j] * bkgd2[j];
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = (qual1[j] + qual2[j])/2;
		   }  /* end loop over j */
		  }  /* end cf_version == 3 */
		} /* end !arg2_scalar */
		break;

	   case DIV:
	      if (arg2_scalar)
		{
		if (arg2_val == 0.)
		   cf_if_error ("Commanded division by 0!");
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] / arg2_val;
		   flxerr_o[j] = flxerr1[j] / arg2_val;
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] / arg2_val;
		   flxerr_o[j] = flxerr1[j] / arg2_val;
		   lcts_o[j] = lcts1[j];
		   qual_o[j] = qual1[j];
		   bkgd_o[j] = bkgd1[j];
		   wgts_o[j] = wgts1[j];
		   }
		  }
		}
	      else
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   if (flux2[j] !=0.)
			{
			flux_o[j] = flux1[j] / flux2[j];
			if (flux1[j] != 0.)
			   {
			   tmp1 = flxerr1[j] / flux1[j];
			   tmp2 = flxerr2[j] / flux2[j];
			   tmp3 = tmp1 * tmp1;
			   tmp4 = tmp2 * tmp2;
			   flxerr_o[j] = fabs (flux_o[j]) * sqrt (tmp3 + tmp4);
			   }
			else
			   flxerr_o[j] = flxerr1[j] / flux2[j];
			if ((qual1[j] == 0) || (qual2[j] == 0))
			   qual_o[j] = 0;
			else
			   qual_o[j] = qual1[j];
			}
		   else
			{
			flux_o[j] = 0.;
			flxerr_o[j] = flxerr1[j];
			qual_o[j] = 0;
			}
		   if (cnts2[j] != 0.)
			{
			cnts_o[j] = cnts1[j] / cnts2[j];
			if (cnts1[j] != 0.)
			   tmp1 = cntserr1[j] / cnts1[j];
			else
			   tmp1 = 0.;
			tmp2 = cntserr2[j] / cnts2[j];
			tmp3 = tmp1 * tmp1;
			tmp4 = tmp2 * tmp2;
			cntserr_o[j] = fabs (cnts_o[j]) * sqrt (tmp3 + tmp4);
			}
		   else
			{
			cnts_o[j] = 0.;
			cntserr_o[j] = cntserr1[j];
			}
		   } /* end loop over j */
		  } /* end cf_version == 2 */
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   if (flux2[j] !=0.)
			{
			flux_o[j] = flux1[j] / flux2[j];
			if (flux1[j] != 0.)
			   {
			   tmp1 = flxerr1[j] / flux1[j];
			   tmp2 = flxerr2[j] / flux2[j];
			   tmp3 = tmp1 * tmp1;
			   tmp4 = tmp2 * tmp2;
			   flxerr_o[j] = fabs (flux_o[j]) * sqrt (tmp3 + tmp4);
			   }
			else
			   flxerr_o[j] = flxerr1[j] / flux2[j];
			if ((qual1[j] == 0) || (qual2[j] == 0))
			   qual_o[j] = 0;
			else
			   qual_o[j] = (qual1[j] + qual2[j])/2;
			}
		   else
			{
			flux_o[j] = 0.;
			flxerr_o[j] = flxerr1[j];
			qual_o[j] = 0;
			}

		   if (lcts2[j] != 0.)
			lcts_o[j] = lcts1[j] / lcts2[j];
		   else
			lcts_o[j] = 0.;

		   if (wgts2[j] != 0.)
			wgts_o[j] = wgts1[j] / wgts2[j];
		   else
			wgts_o[j] = 0.;

		   if (bkgd2[j] != 0.)
			bkgd_o[j] = bkgd1[j] / bkgd2[j];
		   else
			bkgd_o[j] = 0.;

		   } /* end loop over j */
		  } /* end cf_version == 3 */
		} /* end !arg2_scalar */
		break;

	   case AVG_T:
	      if (arg2_scalar)
		{
		cf_if_error ("Illegal operation with scalar operand!");
		}
	      else
		{
		if (cf_version == 2)
		  {
		f1 = exptime_1 / (exptime_1 + exptime_2);
		f2 = exptime_2 / (exptime_1 + exptime_2);
		exptime_o = exptime_1 + exptime_2;
		for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = f1 * flux1[j] + f2 * flux2[j];
		   tmp1 = f1 * f1 * flxerr1[j]*flxerr1[j];
		   tmp2 = f2 * f2 * flxerr2[j]*flxerr2[j];
		   flxerr_o[j] = sqrt (tmp1 + tmp2);
		   cnts_o[j] = f1 * cnts1[j] + f2 * cnts2[j];
		   tmp1 = f1 * f1 * cntserr1[j]*cntserr1[j];
		   tmp2 = f2 * f2 * cntserr2[j]*cntserr2[j];
		   cntserr_o[j] = sqrt (tmp1 + tmp2);
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  f1 = exptime_1 / (exptime_1 + exptime_2);
		  f2 = exptime_2 / (exptime_1 + exptime_2);
		  exptime_o = exptime_1 + exptime_2;
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = f1 * flux1[j] + f2 * flux2[j];
		   tmp1 = f1 * f1 * flxerr1[j]*flxerr1[j];
		   tmp2 = f2 * f2 * flxerr2[j]*flxerr2[j];
		   flxerr_o[j] = sqrt (tmp1 + tmp2);
		   lcts_o[j] = lcts1[j] + lcts2[j];
		   wgts_o[j] = wgts1[j] + wgts2[j];
		   bkgd_o[j] = bkgd1[j] + bkgd2[j];
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = (qual1[j] + qual2[j])/2;
		   }
		  }
		}
		break;

	   case AVG_S:
	      if (arg2_scalar)
		{
		cf_if_error ("Illegal operation with scalar operand!");
		}
	      else
		{
		exptime_o = exptime_1 + exptime_2;
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   if ((flxerr1[j] > 0.) && (flxerr2[j] > 0.))
		      {
		      f1 = 1. / (flxerr1[j] * flxerr1[j]);
		      f2 = 1. / (flxerr2[j] * flxerr2[j]);
		      }
		   else
		      f1 = f2 = 1.;
		   fn = f1 + f2;
		   flux_o[j] = (f1 * flux1[j] + f2 * flux2[j]) / fn;
		   flxerr_o[j] = sqrt (1./fn);
		   if ((cntserr1[j] > 0.) && (cntserr2[j] > 0.))
		      {
		      f1 = 1. / (cntserr1[j] * cntserr1[j]);
		      f2 = 1. / (cntserr2[j] * cntserr2[j]);
		      }
		   else
		      f1 = f2 = 1.;
		   fn = f1 + f2;
		   cnts_o[j] = (f1 * cnts1[j] + f2 * cnts2[j]) / fn;
		   cntserr_o[j] = sqrt (1./fn);
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		  }  /* end if cf_version == 2 */
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   if ((flxerr1[j] > 0.) && (flxerr2[j] > 0.))
		      {
		      f1 = 1. / (flxerr1[j] * flxerr1[j]);
		      f2 = 1. / (flxerr2[j] * flxerr2[j]);
		      }
		   else
		      f1 = f2 = 1.;
		   fn = f1 + f2;
		   flux_o[j] = (f1 * flux1[j] + f2 * flux2[j]) / fn;
		   flxerr_o[j] = sqrt (1./fn);
		   lcts_o[j] = (f1 * lcts1[j] + f2 * lcts2[j]) / fn;
		   wgts_o[j] = (f1 * wgts1[j] + f2 * wgts2[j]) / fn;
		   bkgd_o[j] = (f1 * bkgd1[j] + f2 * bkgd2[j]) / fn;
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = (f1*qual1[j] + f2*qual2[j]) / fn;
		   }
		  }
		}
		break;

	   case REPL_D:
	      if (arg2_scalar)
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = arg2_val;
		   flxerr_o[j] = flxerr1[j];
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = arg2_val;
		   flxerr_o[j] = flxerr1[j];
		   lcts_o[j] = lcts1[j];
		   qual_o[j] = qual1[j];
		   wgts_o[j] = wgts1[j];
		   bkgd_o[j] = bkgd1[j];
		   }
		  }
		}
	      else
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux2[j];
		   flxerr_o[j] = flxerr1[j];
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux2[j];
		   flxerr_o[j] = flxerr1[j];
		   lcts_o[j] = lcts1[j];
		   wgts_o[j] = wgts1[j];
		   bkgd_o[j] = bkgd1[j];
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		  }
		}
		break;

	   case REPL_E:
	      if (arg2_scalar)
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j];
		   flxerr_o[j] = arg2_val;
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j];
		   flxerr_o[j] = arg2_val;
		   lcts_o[j] = lcts1[j];
		   qual_o[j] = qual1[j];
		   wgts_o[j] = wgts1[j];
		   bkgd_o[j] = bkgd1[j];
		   }
		  }
		}
	      else
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j];
		   flxerr_o[j] = flxerr2[j];
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j];
		   flxerr_o[j] = flxerr2[j];
		   lcts_o[j] = lcts1[j];
		   wgts_o[j] = wgts1[j];
		   bkgd_o[j] = bkgd1[j];
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		  }
		}
		break;

	   case DELTA:
	      /* start by copying input file to output */
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j];
		   flxerr_o[j] = flxerr1[j];
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j];
		   flxerr_o[j] = flxerr1[j];
		   lcts_o[j] = lcts1[j];
		   qual_o[j] = qual1[j];
		   wgts_o[j] = wgts1[j];
		   bkgd_o[j] = bkgd1[j];
		   }
		  }
	      if (arg2_scalar)
		{
		for (j=100; j<npix_1; j+=100)
		   {
		   flux_o[j] += arg2_val;
		   if (flux1[j] != 0) 
			{
			flxerr_o[j] *= fabs (flux_o[j] / flux1[j]);
			}
		   else
			{
			flxerr_o[j] += 0.1 * fabs (flux_o[j]);
			}
		   }
		}
	      else
		{
		for (j=100; j<npix_1; j+=100)
		   {
		   flux_o[j] += flux2[j];
		   if (flux1[j] != 0) 
			{
			flxerr_o[j] *= fabs (flux_o[j] / flux1[j]);
			}
		   else
			{
			flxerr_o[j] += 0.1 * fabs (flux_o[j]);
			}
		   }
		}
		break;

	   case BOXCAR:
	      if (arg2_scalar)
		{
		/* if width is <= 1, just copy input to output */
		if (arg2_val <= 1.)
		  {
		  if (cf_version == 2)
		    {
		    for (j=0; j<npix_1; j++)
		     {
		     flux_o[j] = flux1[j];;
		     flxerr_o[j] = flxerr1[j];
		     cnts_o[j] = cnts1[j];
		     cntserr_o[j] = cntserr1[j];
		     qual_o[j] = qual1[j];
		     }
		    } /* end if_cf_version == 2 */
		  else
		    {  /* cf_version == 3 */
		    for (j=0; j<npix_1; j++)
		     {
		     flux_o[j] = flux1[j];;
		     flxerr_o[j] = flxerr1[j];
		     lcts_o[j] = lcts1[j];
		     qual_o[j] = qual1[j];
		     wgts_o[j] = wgts1[j];
		     bkgd_o[j] = bkgd1[j];
		     }
		    } /* end if cf_version == 3 */
		  } /* end if arg2_val <=1 */
		else
		  {
		  if (cf_version == 2)
		   {
		   for (j=0; j<npix_1; j++)
		     {
		     k1 = j - arg2_val / 2. + 0.5;
		     k2 = j + arg2_val / 2. + 0.5;
		     if (k1 < 0)
			k1 = 0;
		     if (k2 > npix_1-1)
			k2 = npix_1-1;
		     sflux = 0.;
		     sferr = 0.;
		     scnts = 0.;
		     scerr = 0.;
		     lqual = 0;
		     for (k=k1; k<=k2; k++)
			{
			sflux += flux1[k];
			sferr += flxerr1[k] * flxerr1[k];
			scnts += cnts1[k];
			scerr += cntserr1[k] * cntserr1[k];
			lqual += (long) qual1[k];
			}
		     spix = (float)(k2 - k1 + 1);
		     flux_o[j] = sflux / spix;
		     flxerr_o[j] = sqrt (sferr) / spix;
		     cnts_o[j] = scnts / spix;
		     cntserr_o[j] = sqrt (scerr) / spix;
		     qual_o[j] = (short) (lqual / spix);
		     } /* end of for j... */
		   } /* end if cf_version == 2 */
		  else
		   { /* start cf_version == 3 for scalar smoothing > 1 */
		   for (j=0; j<npix_1; j++)
		     {
		     k1 = j - arg2_val / 2. + 0.5;
		     k2 = j + arg2_val / 2. + 0.5;
		     if (k1 < 0)
			k1 = 0;
		     if (k2 > npix_1-1)
			k2 = npix_1-1;
		     sflux = 0.;
		     sferr = 0.;
		     slcts = 0;
		     sbkgd = 0.;
		     swgts = 0.;
		     lqual = 0;
		     for (k=k1; k<=k2; k++)
			{
			sflux += flux1[k];
			sferr += flxerr1[k] * flxerr1[k];
			slcts += lcts1[k];
			swgts += wgts1[k];
			sbkgd += bkgd1[k];
			lqual += (long) qual1[k];
			}
		     spix = (float)(k2 - k1 + 1);
		     flux_o[j] = sflux / spix;
		     flxerr_o[j] = sqrt (sferr) / spix;
		     lcts_o[j] = slcts / spix;
		     bkgd_o[j] = sbkgd / spix;
		     wgts_o[j] = swgts / spix;
		     qual_o[j] = (short) (lqual / spix);
		     } /* end of for j... */
		   } /* end cf_version == 3 */
		  } /* end of arg2_scalar > 1 */
		} /* end of if (arg2_scalar) */
	      else
		{
		/* allow the smoothing width to vary pixel-by-pixel */
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   arg2_val = flux2[j];
		   if (arg2_val < 1.)
		     {
		     flux_o[j] = flux1[j];
		     flxerr_o[j] = flxerr1[j];
		     cnts_o[j] = cnts1[j];
		     cntserr_o[j] = cntserr1[j];
		     qual_o[j] = qual1[j];
		     }
		   else
		     {
		     k1 = j - arg2_val / 2. + 0.5;
		     k2 = j + arg2_val / 2. + 0.5;
		     if (k1 < 0)
			k1 = 0;
		     if (k2 > npix_1-1)
			k2 = npix_1-1;
		     sflux = 0.;
		     sferr = 0.;
		     scnts = 0.;
		     scerr = 0.;
		     lqual = 0;
		     for (k=k1; k<=k2; k++)
			{
			sflux += flux1[k];
			sferr += flxerr1[k] * flxerr1[k];
			scnts += cnts1[k];
			scerr += cntserr1[k] * cntserr1[k];
			lqual += (long) qual1[k];
			}
		     spix = (float)(k2 - k1 + 1);
		     flux_o[j] = sflux / spix;
		     flxerr_o[j] = sqrt (sferr) / spix;
		     cnts_o[j] = scnts / spix;
		     cntserr_o[j] = sqrt (scerr) / spix;
		     qual_o[j] = (short) (lqual / spix);
		     } /* end of arg2_val > 1 */
		   } /* end of for j... */
		  } /* end if cf_version == 2 */
		else
		  { /* cf_version == 3 */
		  for (j=0; j<npix_1; j++)
		   {
		   arg2_val = flux2[j];
		   if (arg2_val < 1.)
		     {
		     flux_o[j] = flux1[j];
		     flxerr_o[j] = flxerr1[j];
		     lcts_o[j] = lcts1[j];
		     wgts_o[j] = wgts1[j];
		     bkgd_o[j] = bkgd1[j];
		     qual_o[j] = qual1[j];
		     }
		   else
		     {
		     k1 = j - arg2_val / 2. + 0.5;
		     k2 = j + arg2_val / 2. + 0.5;
		     if (k1 < 0)
			k1 = 0;
		     if (k2 > npix_1-1)
			k2 = npix_1-1;
		     sflux = 0.;
		     sferr = 0.;
		     slcts = 0;
		     sbkgd = 0.;
		     swgts = 0.;
		     lqual = 0;
		     for (k=k1; k<=k2; k++)
			{
			sflux += flux1[k];
			sferr += flxerr1[k] * flxerr1[k];
			slcts += lcts1[k];
			swgts += wgts1[k];
			sbkgd += bkgd1[k];
			lqual += (long) qual1[k];
			}
		     spix = (float)(k2 - k1 + 1);
		     flux_o[j] = sflux / spix;
		     flxerr_o[j] = sqrt (sferr) / spix;
		     lcts_o[j] = slcts / spix;
		     wgts_o[j] = swgts / spix;
		     bkgd_o[j] = sbkgd / spix;
		     qual_o[j] = (short) (lqual / spix);
		     } /* end of arg2_val > 1 */
		   } /* end of for j... */
		  } /* end if cf_version == 3 */
		} /* end of if (!arg2_scalar) */
		break; /* end of case BOXCAR */

	   case MIN:
	      if (arg2_scalar)
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   if (flux1[j] < arg2_val)
			flux_o[j] = arg2_val;
		   else
			flux_o[j] = flux1[j];
		   flxerr_o[j] = flxerr1[j];
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   if (flux1[j] < arg2_val)
			flux_o[j] = arg2_val;
		   else
			flux_o[j] = flux1[j];
		   flxerr_o[j] = flxerr1[j];
		   lcts_o[j] = lcts1[j];
		   wgts_o[j] = wgts1[j];
		   bkgd_o[j] = bkgd1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		}
	      else
		{
		for (j=0; j<npix_1; j++)
		   {
		   if (flux1[j] < flux2[j])
			{
			flux_o[j] = flux2[j];
			flxerr_o[j] = flxerr2[j];
			}
		   else
			{
			flux_o[j] = flux1[j];
			flxerr_o[j] = flxerr1[j];
			}
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		    {
		    cnts_o[j] = cnts1[j];
		    cntserr_o[j] = cntserr1[j];
		    }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		    {
		    lcts_o[j] = lcts1[j];
		    wgts_o[j] = wgts1[j];
		    bkgd_o[j] = bkgd1[j];
		    }
		  }
		}
		break;

	   case MAX:
	      if (arg2_scalar)
		{
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   if (flux1[j] > arg2_val)
			flux_o[j] = arg2_val;
		   else
			flux_o[j] = flux1[j];
		   flxerr_o[j] = flxerr1[j];
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   if (flux1[j] > arg2_val)
			flux_o[j] = arg2_val;
		   else
			flux_o[j] = flux1[j];
		   flxerr_o[j] = flxerr1[j];
		   lcts_o[j] = lcts1[j];
		   wgts_o[j] = wgts1[j];
		   bkgd_o[j] = bkgd1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		}
	      else
		{
		for (j=0; j<npix_1; j++)
		   {
		   if (flux1[j] > flux2[j])
			{
			flux_o[j] = flux2[j];
			flxerr_o[j] = flxerr2[j];
			}
		   else
			{
			flux_o[j] = flux1[j];
			flxerr_o[j] = flxerr1[j];
			}
		   if ((qual1[j] == 0) || (qual2[j] == 0))
			qual_o[j] = 0;
		   else
			qual_o[j] = qual1[j];
		   }
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		    {
		    cnts_o[j] = cnts1[j];
		    cntserr_o[j] = cntserr1[j];
		    }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		    {
		    lcts_o[j] = lcts1[j];
		    wgts_o[j] = wgts1[j];
		    bkgd_o[j] = bkgd1[j];
		    }
		  }
		}
		break;

	   case SHIFT:
		if (cf_version == 2)
		  {
	          if (p_shift > 0)
		    {
		    for (j=0; j<p_shift; j++)
		       {
		       flux_o[j] = 0.;
		       flxerr_o[j] = DEF_ERROR;
		       cnts_o[j] = 0.;
		       cntserr_o[j] = DEF_ERROR;
		       qual_o[j] = 0;
		       }
		    for (j=p_shift; j<npix_1; j++)
		       {
		       k = j - p_shift;
		       flux_o[j] = flux1[k];
		       flxerr_o[j] = flxerr1[k];
		       cnts_o[j] = cnts1[k];
		       cntserr_o[j] = cntserr1[k];
		       qual_o[j] = qual1[k];
		       }
		    }
	          if (p_shift <= 0)
		    {
		    for (j=npix_1+p_shift; j<npix_1; j++)
		       {
		       flux_o[j] = 0.;
		       flxerr_o[j] = DEF_ERROR;
		       cnts_o[j] = 0.;
		       cntserr_o[j] = DEF_ERROR;
		       qual_o[j] = 0;
		       }
		    for (j=0; j<npix_1+p_shift; j++)
		       {
		       k = j - p_shift;
		       flux_o[j] = flux1[k];
		       flxerr_o[j] = flxerr1[k];
		       cnts_o[j] = cnts1[k];
		       cntserr_o[j] = cntserr1[k];
		       qual_o[j] = qual1[k];
		       }
		    }
		  } /* end if cf_version == 2 */
		else
		  {
	          if (p_shift > 0)
		    {
		    for (j=0; j<p_shift; j++)
		       {
		       flux_o[j] = 0.;
		       flxerr_o[j] = DEF_ERROR;
		       lcts_o[j] = 0;
		       qual_o[j] = 0;
		       wgts_o[j] = 0.;
		       bkgd_o[j] = 0.;
		       }
		    for (j=p_shift; j<npix_1; j++)
		       {
		       k = j - p_shift;
		       flux_o[j] = flux1[k];
		       flxerr_o[j] = flxerr1[k];
		       lcts_o[j] = lcts1[k];
		       wgts_o[j] = wgts1[k];
		       bkgd_o[j] = bkgd1[k];
		       qual_o[j] = qual1[k];
		       }
		    }
	          if (p_shift <= 0)
		    {
		    for (j=npix_1+p_shift; j<npix_1; j++)
		       {
		       flux_o[j] = 0.;
		       flxerr_o[j] = DEF_ERROR;
		       lcts_o[j] = 0;
		       wgts_o[j] = 0.;
		       bkgd_o[j] = 0.;
		       qual_o[j] = 0;
		       }
		    for (j=0; j<npix_1+p_shift; j++)
		       {
		       k = j - p_shift;
		       flux_o[j] = flux1[k];
		       flxerr_o[j] = flxerr1[k];
		       lcts_o[j] = lcts1[k];
		       wgts_o[j] = wgts1[k];
		       bkgd_o[j] = bkgd1[k];
		       qual_o[j] = qual1[k];
		       }
		    }
		  }
		break;

	   case BIN:
		/* if nbin = 1, just copy input to output */
		if (nbin == 1)
		   {
		   if (cf_version == 2)
		     {
		      for (j=0; j<npix_o; j++)
		        {
		        wave_o[j] = wave1[j];
		        flux_o[j] = flux1[j];;
		        flxerr_o[j] = flxerr1[j];
		        cnts_o[j] = cnts1[j];
		        cntserr_o[j] = cntserr1[j];
		        qual_o[j] = qual1[j];
		        }
		     }
		   else
		     {
		      for (j=0; j<npix_o; j++)
		        {
		        wave_o[j] = wave1[j];
		        flux_o[j] = flux1[j];;
		        flxerr_o[j] = flxerr1[j];
		        lcts_o[j] = lcts1[j];
		        wgts_o[j] = wgts1[j];
		        bkgd_o[j] = bkgd1[j];
		        qual_o[j] = qual1[j];
		        }
		     } /* end cf_version == 3 */
		   } /* end nbin == 1 */
		else
		   {
		   if (cf_version == 2)
		     {
		      for (j=0; j<npix_o; j++)
		        {
		        k1 = j * nbin;
		        k2 = k1 + nbin;
		        if (k2 > npix_1)
			   k2 = npix_1;
		        swave = 0.;
		        sflux = 0.;
		        sferr = 0.;
		        scnts = 0.;
		        scerr = 0.;
		        squal = 0;
		        spix = 0;
		        for (k=k1; k<k2; k++)
			   {
			   swave += wave1[k];
			   if (qual1[k] > 0)
			      {
			      sflux += flux1[k];
			      sferr += flxerr1[k] * flxerr1[k];
			      scnts += cnts1[k];
			      scerr += cntserr1[k] * cntserr1[k];
			      squal += qual1[k];
			      spix++;
			      }
			   }
		        if (spix == 0)
			   spix = 1;
		        wave_o[j] = swave / (float)(k2 - k1);
		        flux_o[j] = sflux / spix;
		        flxerr_o[j] = sqrt (sferr) / spix;
		        cnts_o[j] = scnts;
		        cntserr_o[j] = sqrt (scerr);
		        qual_o[j] = squal;
		        } /* end of for j... */
		     } /* end cf_version == 2 */
		   else
		     {
		      for (j=0; j<npix_o; j++)
		        {
		        k1 = j * nbin;
		        k2 = k1 + nbin;
		        if (k2 > npix_1)
			   k2 = npix_1;
		        swave = 0.;
		        sflux = 0.;
		        sferr = 0.;
		        slcts = 0;
		        swgts = 0.;
		        sbkgd = 0.;
		        squal = 0;
		        spix = 0;
		        for (k=k1; k<k2; k++)
			   {
			   swave += wave1[k];
			   if (qual1[k] > 0)
			      {
			      sflux += flux1[k];
			      sferr += flxerr1[k] * flxerr1[k];
			      slcts += lcts1[k];
			      swgts += wgts1[k];
			      sbkgd += bkgd1[k];
			      squal += qual1[k];
			      spix++;
			      }
			   }
		        if (spix == 0)
			   spix = 1;
		        wave_o[j] = swave / (float)(k2 - k1);
		        flux_o[j] = sflux / spix;
		        flxerr_o[j] = sqrt (sferr) / spix;
		        lcts_o[j] = slcts;
		        qual_o[j] = squal / spix;
		        wgts_o[j] = swgts / spix;
		        bkgd_o[j] = sbkgd;
		        } /* end of for j... */
		     } /* end cf_version == 3 */
		   } /* end of nbin > 1 */
		break;

	   case DERIV:
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1-1; j++)
		   {
		   flux_o[j] = flux1[j+1] - flux1[j];
		   tmp1 = flxerr1[j+1]*flxerr1[j+1];
		   tmp2 = flxerr1[j]*flxerr1[j];
		   flxerr_o[j] = sqrt (tmp1 + tmp2);
		   tmp1 = cntserr1[j+1]*cntserr1[j+1];
		   tmp2 = cntserr1[j]*cntserr1[j];
		   cnts_o[j] = cnts1[j+1] - cnts1[j];
		   cntserr_o[j] = sqrt (tmp1 + tmp2);
		   qual_o[j] = qual1[j];
		   }
		  flux_o[npix_1-1] = flux_o[npix_1-2];
		  flxerr_o[npix_1-1] = flxerr_o[npix_1-2];
		  cnts_o[npix_1-1] = cnts_o[npix_1-2];
		  cntserr_o[npix_1-1] = cntserr_o[npix_1-2];
		  } /* end cf_version == 2 */
		else
		  {
		  for (j=0; j<npix_1-1; j++)
		   {
		   flux_o[j] = flux1[j+1] - flux1[j];
		   tmp1 = flxerr1[j+1]*flxerr1[j+1];
		   tmp2 = flxerr1[j]*flxerr1[j];
		   flxerr_o[j] = sqrt (tmp1 + tmp2);
		   lcts_o[j] = lcts1[j+1] - lcts1[j];
		   qual_o[j] = qual1[j];
		   wgts_o[j] = wgts1[j];
		   bkgd_o[j] = bkgd1[j];
		   }
		  flux_o[npix_1-1] = flux_o[npix_1-2];
		  flxerr_o[npix_1-1] = flxerr_o[npix_1-2];
		  lcts_o[npix_1-1] = lcts_o[npix_1-2];
		  } /* end cf_version == 3 */
		break;

	   case ADDNOISE:
		if (cf_version == 2)
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] + arg2_val * rand1();
		   flxerr_o[j] = flxerr1[j];
		   cnts_o[j] = cnts1[j];
		   cntserr_o[j] = cntserr1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		else
		  {
		  for (j=0; j<npix_1; j++)
		   {
		   flux_o[j] = flux1[j] + arg2_val * rand1();
		   flxerr_o[j] = flxerr1[j];
		   lcts_o[j] = lcts1[j];
		   wgts_o[j] = wgts1[j];
		   bkgd_o[j] = bkgd1[j];
		   qual_o[j] = qual1[j];
		   }
		  }
		break;

	   default:
		cf_if_error ("Illegal opcode; should never get here!");
		break;

	   } /* end switch (opcode) */
 
	/* write the extension containing the results as a binary table */
	hdunum_o=i+1;

	/* cf_write_spec creates the new extension and sets the pointer
	 * to it.
	 */
	if (tfields_1 < 7)
	   cf_wrspec_cf2 (outfits, npix_o, wave_o, flux_o, flxerr_o, qual_o,
		cnts_o, cntserr_o, tfields_1, areaflag);
	else
	   cf_wrspec7 (outfits, npix_o, wave_o, flux_o, flxerr_o,
		lcts_o, wgts_o, bkgd_o, qual_o);

 
	if (wave_o != wave1)
	  free(wave_o);
	free(wave1);
	free(flux1);
	free(flxerr1);
	free(qual1);
	free(flux_o);
	free(flxerr_o);
	free(qual_o);
	if (cf_version == 2)
	   {
	   free(cnts1);
	   free(cntserr1);
	   free(cnts_o);
	   free(cntserr_o);
	   }
	else
	   {
	   free(lcts1);
	   free(bkgd1);
	   free(wgts1);
	   free(lcts_o);
	   free(bkgd_o);
	   free(wgts_o);
	   }

	/* is there only 1 extension in operand 2, so that we need to keep
	 * the data we've read in, or do we need to make space to read the
	 * next extension?
	 */
	if (!arg2_scalar && (i < next_2))
	   {
	   free(flux2);
	   free(flxerr2);
	   free(qual2);
	   if (cf_version == 2)
		{
		free(cnts2);
		free(cntserr2);
		}
	   else
		{
		free(lcts2);
		free(wgts2);
		free(bkgd2);
		}
	   }

	} /* end of loop over next_1 */

   /* add HISTORY lines to output file */
   hdunum_o = 1;
   FITS_movabs_hdu (outfits, hdunum_o, &hdutype_o, &status);

   sprintf (comment, "cf_arith:  operation = %s", argv[2]);
   FITS_write_history (outfits, comment, &status);
   sprintf (comment, "cf_arith:  operand1 = %s", infile1);
   FITS_write_history (outfits, comment, &status);
   sprintf (comment, "cf_arith:  operand2 = %s", infile2);
   FITS_write_history (outfits, comment, &status);

   /* update exposure time in output header if necessary */
   if (exptime_o != exptime_1)
	{
	FITS_update_key (outfits, TFLOAT, "EXPTIME", &exptime_o, expt_comment,
			&status);
	sprintf (comment, "cf_arith: EXPTIME = EXPTIME1 + EXPTIME2");
	FITS_write_history (outfits, comment, &status);
	}

   /* update FILENAME keyword in output file */
   FITS_update_key (outfits, TSTRING, "FILENAME", outfile, NULL, &status);

   time (&timeval);
   timestr = ctime (&timeval);
   sprintf (comment, "cf_arith:  completed at %.24s", timestr);
   /* over-write the trailing \n in timestr */
/* following not needed when field width is specified as 24 above 
   comment[strlen(comment)-1] = '\0';
   padstr (comment, FLEN_CARD);
*/
   FITS_write_history (outfits, comment, &status);

   /* Close the files. */

    FITS_close_file(infits_1, &status);
    if (!arg2_scalar)
	FITS_close_file(infits_2, &status);
    FITS_close_file(outfits, &status);

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
#ifdef SOLARIS
    ieee_retrospective(stdout);
#endif
    return 0;
}

int	chkop (opstr)
char	*opstr;
{
	if (!strcmp (opstr, "+"))
	   return (PLUS);
	else if (!strcmp (opstr, "-"))
	   return (MINUS);
	else if (!strcmp (opstr, "*"))
	   return (MULT);
	else if (!strcmp (opstr, "/"))
	   return (DIV);
	else if (!strcmp (opstr, "avg_t"))
	   return (AVG_T);
	else if (!strcmp (opstr, "avg_s"))
	   return (AVG_S);
	else if (!strcmp (opstr, "repl_d"))
	   return (REPL_D);
	else if (!strcmp (opstr, "repl_e"))
	   return (REPL_E);
	else if (!strcmp (opstr, "boxcar"))
	   return (BOXCAR);
	else if (!strcmp (opstr, "delta"))
	   return (DELTA);
	else if (!strcmp (opstr, "min"))
	   return (MIN);
	else if (!strcmp (opstr, "max"))
	   return (MAX);
	else if (!strcmp (opstr, "shift"))
	   return (SHIFT);
	else if (!strcmp (opstr, "bin"))
	   return (BIN);
	else if (!strcmp (opstr, "deriv"))
	   return (DERIV);
	else if (!strcmp (opstr, "addnoise"))
	   return (ADDNOISE);
	else
	   return (0);
}

/* check_num (in_str) returns TRUE if in_str is a string representation of
 * a valid number, FALSE otherwise.
 * This is a simple-minded check:  just see if first char is +,-,digit.
 */
int	check_num (in_str)
char	*in_str;
{
	if ((in_str[0] == '+') || (in_str[0] == '-'))
	   {
	   if (isdigit (in_str[1]))
		return (TRUE);
	   else
		return (FALSE);
	   }
	if (isdigit (in_str[0]))
	   return (TRUE);

	return (FALSE);
}

void	padstr (pstr, len)
char	*pstr;
int	len;
{
	int	i, slen;

	slen = strlen (pstr);
	for (i=slen; i<len; i++)
	   *(pstr + i) = ' ';

	return;
}
/*      return a random number in range -1 -> 1 */
double  rand1 ()
{
        long    random();
        double  r;

        r =  2. * (double) random() / (double)MAX_RAND - 1.;
        return (r);
}

