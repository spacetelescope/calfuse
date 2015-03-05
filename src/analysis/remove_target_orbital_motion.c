
/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * 
 *
 * Usage:
 *     remove_target_orbital_motion [-hb] [-v level]
 *                                   input_idf_file output_idf_file
 *                                   hjd0 period vrmax
 *                                   [ra_h ra_m ra_s dec_d dec_m dec_s]
 *
 * Arguments:
 *       mjd0                              :  Heliocentric or Geocentric
 *                                            Julian Day 
 *                                            of closest approch on orbit.
 *       period                            :  Period in seconds.
 *       vrmax                             :  Maximum radial velocity in km/s
 *       ra_h ra_m ra_s dec_d dec_m dec_s  :  (Optional) RA and DEC of target.
 *                                            If present mjd0 is interpreted 
 *                                            as Heliocentric.
 *
 * Options:
 *     -h:  this help message
 *     -v:  verbosity level (=1; 0 is silent)
 *     -b:  update timeline doppler information (required for BPM)
 *
 *
 * History:	12/01/03	bjg	v1.0	Begin work.
 *              12/02/03        bjg             Changed calling 
 *                                              parameters format.        	
 *              12/04/03        bjg             Number of elements in HDU2 
 *                                              is now determined from HDU2 
 *                                              header instead of HDU1 
 *                                              header NEVENTS keyword
 *              12/10/03        bjg		Updated header.
 *              04/05/04        bjg             Remove unused variables  
 *              06/10/04        bjg             Some code cleaning    
 *              10/07/04        bjg     v1.1    Added heliocentric
 *                                              time correction (optional)
 *                                              Added timeline doppler
 *                                              correction information
 *                                              (optional) for cf_bad_pixels 
 *              10/26/04        bjg     v1.2    Cosmetic change
 *
 ****************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "calfuse.h"




#define LIGHT_SPEED 299792.458 /* speed of light in kilometers per second */
#define SECONDS_PER_DAY (24.0*3600.0)

static char CF_PRGM_ID[]= "remove_target_orbital_motion";
static char CF_VER_NUM[]= "1.2";

struct key {
  char keyword[FLEN_KEYWORD];
  float value;
};


double gethmjd(double, int, int, float, int, int, float);

int main(int argc,char *argv[]){
  
 
  char     fmt_byte[FLEN_CARD],fmt_float[FLEN_CARD],fmt_short[FLEN_CARD];  
  
  fitsfile * infits;
  fitsfile * outfits; 

  char       input_filename[FLEN_CARD],output_filename[FLEN_CARD];  
  char *     corrected_filename;
  char *     string_pointer;

  int        intnull=0,anynull;
  int        status=0;
  int        hdutype=0;
  int        ncol;
  char       tempstring[FLEN_CARD]; 
  char       tempchar;
  char       keyword[FLEN_CARD],card[FLEN_CARD];
 

  struct key hdu4_tscal[16];
  struct key hdu4_tzero[16];
  
  char hdu4_extname[]="TIMELINE";
  int hdu4_tfields=16;  /* output table will have 16 columns */
  char *hdu4_ttype[]={"TIME", "STATUS_FLAGS", "TIME_SUNRISE", "TIME_SUNSET",
		      "LIMB_ANGLE", "LONGITUDE", "LATITUDE", "ORBITAL_VEL",
		      "HIGH_VOLTAGE", "LIF_CNT_RATE", "SIC_CNT_RATE",
		      "FEC_CNT_RATE", "AIC_CNT_RATE", "BKGD_CNT_RATE",
		      "YCENT_LIF","YCENT_SIC"};

  char *hdu4_tform[16]; /* we will define tform later, when the number
			   of photons is known */
  
  char *hdu4_tunit[]={"seconds", "unitless", "seconds", "seconds", "degrees",
		      "degrees", "degrees",  "km/s", "unitless", "counts/sec",
		      "counts/sec", "counts/sec", "counts/sec", "counts/sec",
		      "pixels","pixels" };


  long       nevents, nrecords;
  long       i;

  float *    floatlist;
  short *    shortlist;
  char  *    charlist;

  double *    hdu4_time;
  double *    hdu4_orbvel;

  double *    hdu2_time;
  double *    hdu2_lambda;

  double      delta_t;
  double      period, vrmax, mjd0;
  double      expstart, mjd;
 
  int        rah,ram,decd,decm;
  float      ras,decs;

  int        optc;

  char opts[] = "hbv:";

  char usage[] = 
    "Usage:\n"
    "  remove_target_orbital_motion [-hb] [-v level] \n"
    "    input_idf_file output_idf_file\n" 
    "    mjd0 period vrmax\n"
    "    [ra_h ra_m ra_s dec_d dec_m dec_s]\n";
  char argument[] =
    "Arguments:\n"
    "    mjd0      :  (Geocentric or Heliocentric) Modified Julian Day\n"
    "                 of closest approch on orbit.\n"
    "    period    :  Period in seconds.\n"
    "    vrmax     :  Maximum radial velocity in km/s\n"
    "    ra_h ra_m ra_s dec_d dec_m dec_s  :  Optional RA and DEC of target.\n" 
    "                                         If present mjd0 is interpreted\n"
    "                                         as Heliocentric. \n";

  char option[] =
    "Options:\n"
    "  -h:  this help message\n"
    "  -v:  verbosity level (=1; 0 is silent)\n"
    "  -b:  update timeline doppler information (required for BPM)\n";
  
  
  
  int        bigtime=1; /* How to store TSUNRISE and TSUNSET" */
  int        do_timeline=0;
  int        helio_corr=0;
  
  verbose_level = 1;
 
  
  /* Check number of options and arguments */
  while ((optc = getopt(argc, argv, opts)) != -1) {
    switch(optc) {
    case 'h':
      printf("%s\n%s\n%s", usage, argument, option);
      return EXIT_SUCCESS;
    case 'v':
      verbose_level = atoi(optarg);
      break; 
    case 'b':
      do_timeline=1;
      break;
    }
  }

  /* Initialize error checking. */
  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
  if ((argc != optind+11)&&(argc != optind+5))
    cf_if_error("%s\nIncorrect number of program arguments", usage);
  
  if (argc==optind+11) helio_corr=1;
  
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  
  strcpy(input_filename,argv[optind]);
  strcpy(output_filename,argv[optind+1]);

  sscanf(argv[optind+2],"%lf",&mjd0);
  sscanf(argv[optind+3],"%lf",&period);
  sscanf(argv[optind+4],"%lf",&vrmax);  
  
  cf_verbose(1,"Target Orbital Motion Doppler Correction");
  cf_verbose(1,"Processing file %s",input_filename);
  
  if (helio_corr) {
    sscanf(argv[optind+5],"%d",&rah);
    sscanf(argv[optind+6],"%d",&ram);  
    sscanf(argv[optind+7],"%f",&ras);
    sscanf(argv[optind+8],"%d",&decd);
    sscanf(argv[optind+9],"%d",&decm);
    sscanf(argv[optind+10],"%f",&decs);
    cf_verbose(1,"Right Ascension of target: %d h %d m %f s",rah,ram,ras);
    cf_verbose(1,"Declination of target: %d deg %d m %f s",decd,decm,decs);
    cf_verbose(1,"Heliocentric Modified Julian Day of closest approch on orbit: %lf",mjd0);
  } else {
    cf_verbose(1,"Geocentric Modified Julian Day of closest approch on orbit: %lf",mjd0);
  }
  
  cf_verbose(1,"Period in seconds: %lf",period);
  cf_verbose(1,"Maximum radial velocity in kilometers per second: %lf",vrmax);
  
  

  FITS_open_file(&infits,input_filename,READONLY,&status);
  FITS_read_key(infits,TDOUBLE,"EXPSTART",&(expstart),NULL,&status);
   
  FITS_create_file(&outfits,output_filename,&status);
  FITS_copy_hdu(infits,outfits,0,&status);
  
  string_pointer=strrchr(output_filename,'/');
  if (string_pointer==NULL) corrected_filename=output_filename;
  else corrected_filename=&(string_pointer[1]);
    
  FITS_update_key(outfits,TSTRING,"FILENAME",corrected_filename,NULL,&status);

  strcpy(tempstring,"");

  FITS_update_key(outfits,TSTRING,"BPM_CAL",tempstring,NULL,&status);

  
  
  FITS_movabs_hdu(infits,2,&hdutype,&status);
  FITS_read_key(infits,TSTRING,"TFORM1",tempstring,NULL,&status);
  sscanf(tempstring,"%ld%c",&nevents,&tempchar);

  FITS_create_hdu(outfits,&status);
  FITS_copy_hdu(infits,outfits,0,&status);
  
  
  hdu2_time=(double *)malloc(nevents*sizeof(double));
  hdu2_lambda=(double *)malloc(nevents*sizeof(double));

    
   
  FITS_get_colnum(infits, TRUE, "TIME", &ncol, &status);
  FITS_read_col(infits, TDOUBLE, ncol, 1, 1, nevents, &intnull,
		hdu2_time, &anynull, &status);
 
  FITS_get_colnum(infits, TRUE, "LAMBDA", &ncol, &status);
  FITS_read_col(infits, TDOUBLE, ncol, 1, 1, nevents, &intnull,
		hdu2_lambda, &anynull, &status);

 
		
  
  cf_verbose(1,"Geocentric Modified Julian Day of start of observation/exposure is: %f",expstart);
  if (helio_corr) {
    mjd=gethmjd(expstart,rah,ram,ras,decd,decm,decs);
    cf_verbose(1,"Heliocentric Modified Julian Day of start of observation/exposure is: %f",mjd);
  }

  for (i=0;i<nevents;i++){
    
    mjd=expstart+(hdu2_time[i])/SECONDS_PER_DAY;

    if (helio_corr) mjd=gethmjd(mjd,rah,ram,ras,decd,decm,decs);
      
    delta_t=(mjd-mjd0)*SECONDS_PER_DAY;

    hdu2_lambda[i]=hdu2_lambda[i]*(1-(vrmax/LIGHT_SPEED)*sin(2*M_PI*delta_t/period));        

  }

  FITS_get_colnum(outfits, TRUE, "LAMBDA", &ncol, &status);
  FITS_write_col(outfits, TDOUBLE, ncol, 1, 1, nevents, hdu2_lambda, &status);


  free(hdu2_time);
  free(hdu2_lambda);

   
  FITS_movabs_hdu(infits,3,&hdutype,&status);
  FITS_create_hdu(outfits,&status);
  FITS_copy_hdu(infits,outfits,0,&status);
  

  FITS_movabs_hdu(infits,4,&hdutype,&status);
  
  if (do_timeline) {	

    FITS_read_key(infits,TSTRING,"TFORM1",tempstring,NULL,&status);
    sscanf(tempstring,"%ld%c",&nrecords,&tempchar);		

    /* Generate the tform array */
    sprintf(fmt_byte,  "%ldB", nrecords);
    sprintf(fmt_float, "%ldE", nrecords);
    sprintf(fmt_short, "%ldI", nrecords);
  
  
  
    hdu4_tform[0] = fmt_float;
    hdu4_tform[1] = fmt_byte;

    for (i=2; i<hdu4_tfields; i++)	hdu4_tform[i] = fmt_short;
  
    /* Set TSCALE and TZERO values */
    /* Not all of these will be used; see below. */
    for (i=2; i<hdu4_tfields; i++) {
      sprintf(hdu4_tscal[i].keyword, "TSCAL%ld", i+1);
      hdu4_tscal[i].value = 0.1;
      sprintf(hdu4_tzero[i].keyword, "TZERO%ld", i+1);
      hdu4_tzero[i].value = 0.;
    }
  
    /* Need TZEROn = 3000. for tsunrise and tsunset. */
    hdu4_tzero[2].value = hdu4_tzero[3].value = 3000.;

    /* Let's keep three significant figures for vorb. */
    hdu4_tscal[7].value = 0.01;

    /* FEC_CNT_RATE and AIC_CNT_RATE can get big, so we
       (essentially) store them as unsigned shorts. */
    hdu4_tscal[11].value = 1;
    hdu4_tscal[12].value = 1;
    hdu4_tzero[11].value = 32768;
    hdu4_tzero[12].value = 32768;
    
  
    if (bigtime) { 
      hdu4_tform[2]=hdu4_tform[3]=fmt_float; 
    }

    hdu4_tform[7]=fmt_float;
	
    /* Append a new empty binary table to the output file */
    FITS_create_tbl(outfits, BINARY_TBL, 1, hdu4_tfields, hdu4_ttype, hdu4_tform,
		    hdu4_tunit, hdu4_extname, &status);

  
	
    /* Write TSCALE and TZERO entries to header */
    for (i=2; i<hdu4_tfields; i++) {

      /* If bigtime is TRUE, store tsunrise and tsunset as floats.*/
      if ((i == 2 || i == 3) && bigtime) continue;
    
      /* Omit TSCALE and TZERO for these four arrays. */
      if (i ==  8) continue;	/* HIGH_VOLTAGE */
      if (i ==  9) continue;	/* LIF_CNT_RATE */
      if (i == 10) continue;	/* SIC_CNT_RATE */
      if (i == 13) continue;	/* BKGD_CNT_RATE */


      /* Store ORBITAL_VEL as float */
      if (i == 7) continue;


      sprintf(keyword, "TUNIT%ld", i+1);
      if (fits_read_keyword(outfits, keyword, card, NULL, &status))
	cf_if_fits_error(status);
      FITS_insert_key_flt(outfits, hdu4_tscal[i].keyword, hdu4_tscal[i].value,
			  -1, NULL, &status);
      FITS_insert_key_flt(outfits, hdu4_tzero[i].keyword, hdu4_tzero[i].value,
			  -5, NULL, &status);
    
    } 
  
    floatlist=(float *)malloc(nrecords*sizeof(float));
    shortlist=(short *)malloc(nrecords*sizeof(short));
    charlist=(char *)malloc(nrecords*sizeof(char));
    hdu4_time=(double *)malloc(nevents*sizeof(double));
    hdu4_orbvel=(double *)malloc(nevents*sizeof(double));
  
    FITS_get_colnum(infits, TRUE, "TIME", &ncol, &status);
    FITS_read_col(infits, TDOUBLE, ncol, 1, 1, nrecords, &intnull,
		  hdu4_time, &anynull, &status);	
    FITS_write_col(outfits, TDOUBLE, 1, 1, 1, nrecords, hdu4_time, &status);
  
    FITS_get_colnum(infits, TRUE, "STATUS_FLAGS", &ncol, &status);
    FITS_read_col(infits, TBYTE, ncol, 1, 1, nrecords, &intnull,
		  charlist, &anynull, &status);
    FITS_write_col(outfits, TBYTE, 2, 1, 1, nrecords, charlist, &status);
  
    FITS_get_colnum(infits, TRUE, "TIME_SUNRISE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 3, 1, 1, nrecords, floatlist, &status);
  
    FITS_get_colnum(infits, TRUE, "TIME_SUNSET", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 4, 1, 1, nrecords, floatlist, &status);
  
    FITS_get_colnum(infits, TRUE, "LIMB_ANGLE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 5, 1, 1, nrecords, floatlist, &status);
  
    FITS_get_colnum(infits, TRUE, "LONGITUDE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 6, 1, 1, nrecords, floatlist, &status);
  
    FITS_get_colnum(infits, TRUE, "LATITUDE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 7, 1, 1, nrecords, floatlist, &status);
    
    FITS_get_colnum(infits, TRUE, "ORBITAL_VEL", &ncol, &status);
    FITS_read_col(infits, TDOUBLE, ncol, 1, 1, nrecords, &intnull,
		  hdu4_orbvel, &anynull, &status);

    for (i=0;i<nrecords;i++){
    
      mjd=expstart+(hdu4_time[i])/SECONDS_PER_DAY;

      if (helio_corr) mjd=gethmjd(mjd,rah,ram,ras,decd,decm,decs);
      
      delta_t=(mjd-mjd0)*SECONDS_PER_DAY;

      hdu4_orbvel[i]=hdu4_orbvel[i]-vrmax*sin(2*M_PI*delta_t/period);        

    }


    FITS_write_col(outfits, TDOUBLE, 8, 1, 1, nrecords, hdu4_orbvel, &status);
    
    FITS_get_colnum(infits, TRUE, "HIGH_VOLTAGE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 9, 1, 1, nrecords, floatlist, &status);
    
    FITS_get_colnum(infits, TRUE, "LIF_CNT_RATE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 10, 1, 1, nrecords, floatlist, &status);
    
    FITS_get_colnum(infits, TRUE, "SIC_CNT_RATE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 11, 1, 1, nrecords, floatlist, &status);
		
    FITS_get_colnum(infits, TRUE, "FEC_CNT_RATE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 12, 1, 1, nrecords, floatlist, &status);
    
    FITS_get_colnum(infits, TRUE, "AIC_CNT_RATE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 13, 1, 1, nrecords, floatlist, &status);
    
    FITS_get_colnum(infits, TRUE, "BKGD_CNT_RATE", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);	      
    FITS_write_col(outfits, TFLOAT, 14, 1, 1, nrecords, floatlist, &status);
    
    FITS_get_colnum(infits, TRUE, "YCENT_LIF", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 15, 1, 1, nrecords, floatlist, &status);
    
    FITS_get_colnum(infits, TRUE, "YCENT_SIC", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecords, &intnull,
		  floatlist, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 16, 1, 1, nrecords, floatlist, &status);
    
    
    free(hdu4_time);
    free(hdu4_orbvel);
    
    free(floatlist);
    free(charlist);
    free(shortlist);

  } else {
    FITS_create_hdu(outfits,&status);
    FITS_copy_hdu(infits,outfits,0,&status);
  }
       
  FITS_close_file(infits,&status);
  FITS_close_file(outfits,&status);	
    

  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return EXIT_SUCCESS;
  
}

