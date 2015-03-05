
/*************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *************************************************************************
 *
 *
 * Usage:
 *     ttag_lightcurve [-hf] [-v level] input_file output_file windows_file 
 *                                     bins LiForSiC
 * 
 *
 * Arguments:
 *   input_file     : ttag IDF file                                         
 *   output_file    : ASCII file with 2 columns : 
 *                        - time (MJD)         
 *                        - countrate                       
 *
 *   windows_file   : input ASCII file with 2 columns : 
 *                        - start of spectral windows 
 *                        - stop of spectral windows  
 *
 *   bins           : bin size in seconds                                         
 *   LiForSiC       : 1=LiF, 2=SiC
 *
 *
 *
 *
 * Options:
 *     -h:  this help message
 *     -f:  output calibrated flux (erg/cm2/s) instead of countrate
 *     -v:  verbosity level (=1; 0 is silent)
 *          
 *
 *
 *
 * History:	10/06/04	bjg	v1.0              
 *		06/03/05	wvd	v1.1 Delete unused variables.
 *   
 ***************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "calfuse.h"




#define N_MAX_WIN 100   /* Maximum number of spectral windows */

typedef struct {
  float start;
  float stop;
} interval;             /* type to store a spectral window  */ 


int is_inside(float val,interval * interval_ptr, int nintervals);
int is_good_time(char timeflag,char * daynight);
int is_in_lif_channel(char channel, char * aperture);
int is_in_sic_channel(char channel, char * aperture);





static char CF_PRGM_ID[]= "ttag_lightcurve";
static char CF_VER_NUM[]= "1.1";




int main(int argc,char *argv[]){
  
 
  
  char       input_filename[FLEN_CARD];
  char       output_filename[FLEN_CARD]; 
  char       windefs_filename[FLEN_CARD];
  
  FILE     * output_file;
  FILE     * windefs_file;

  fitsfile * infits;

  interval   swindows[N_MAX_WIN];   /* the spectral windows are stored 
                                           in an array of intervals */
  int        nwindows = 0;          /* number of spectral windows TBD */ 

  float      bins;    

  int        intnull=0,anynull;
  int        status=0;
  int        hdutype=0;

  long       nevents;
  long       i;
  int        ncol;

  char       tempchar;
  char       tempstring[FLEN_CARD];

  char       daynight[FLEN_CARD];
  char       instmode[FLEN_CARD];
  char       aperture[FLEN_CARD];

  double *   hdu2_time;
  double *   hdu2_lambda;
  double *   hdu2_weight;
  char   *   hdu2_channel;
  char   *   hdu2_timeflgs;
  char   *   hdu2_loc_flgs;

  float  *   gti_start;
  float  *   gti_stop;
  int        ngtis;
  
  interval * gtis;

  double t, tmax, ctot;

  double expstart;

  int badtime;

  int LiForSiC;
  

  int optc;
  
  int do_flux=0;

  char opts[] = "hfv:";
  char usage[] = 
    "Usage:\n"
    "  ttag_lightcurve [-hf] [-v level] input_file output_file windows_file\n"
    "                                  bins LiForSiC\n\n"
    "Arguments:\n"
    "  input_file     : ttag IDF file\n\n"
    "  output_file    : ASCII file with 2 columns:\n" 
    "                        - time (MJD)\n"
    "                        - countrate corrected for deadtime\n\n"
    "  windows_file   : input ASCII file with 2 columns :\n"
    "                        - start of spectral windows\n"
    "                        - stop of spectral windows\n\n"
    "  bins           : bin size in seconds\n"
    "  LiForSiC       : 1=LiF, 2=SiC\n";
  char option[] =
    "Options:\n"
    "  -h:  this help message                 \n"
    "  -f:  output calibrated flux (erg/cm2/s) instead of countrate  \n"
    "  -v:  verbosity level (=1; 0 is silent) \n";

  verbose_level = 1;
  
  
  /* Check number of options and arguments */
  while ((optc = getopt(argc, argv, opts)) != -1) {
    switch(optc) {    
    case 'h':
      printf("%s\n%s", usage, option);
      return EXIT_SUCCESS;
    case 'v':
      verbose_level = atoi(optarg);
      break;
    case 'f':
      do_flux=1;
      break;
    }
  }

  /* Initialize error checking. */
  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
  if (argc != optind+5)
    cf_if_error("%s\nIncorrect number of program arguments", usage);
  
  
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  
  strcpy(input_filename,argv[optind]); 
  strcpy(output_filename,argv[optind+1]);  
  strcpy(windefs_filename,argv[optind+2]);   
  sscanf(argv[optind+3],"%f",&bins);  
  sscanf(argv[optind+4],"%d",&LiForSiC);

  if ((LiForSiC>2)||(LiForSiC<1)) {
    cf_if_error("%s\nLiForSiC value must be 1 (LiF) or 2 (SiC)", usage);
  }
  
  
  if (bins<=0) {
    cf_if_error("Time interval is less than or equal to zero.");
  }

    
  /* Open the spectral windows definition file */
  if ((windefs_file = fopen (windefs_filename, "r")) == NULL){
    cf_if_error ("Unable to open file %s\n", windefs_filename);
  }


  /* Read the spectral windows definition file and store the spectral windows */
  while (fscanf(windefs_file,"%f %f",&((swindows[nwindows]).start),
                                     &((swindows[nwindows]).stop))!=EOF)
    {
      nwindows++;
    }

  


  for (i=0;i<nwindows;i++){
    printf("Window %ld : %f - %f\n",i+1,swindows[i].start,swindows[i].stop);
  }

  

  /* Create the output file */
  if ((output_file = fopen(output_filename,"w")) == NULL){
    fclose(windefs_file);
    cf_if_error ("Unable to create file %s\n", output_filename);    
  }
  


  FITS_open_file(&infits,input_filename,READONLY,&status);

  FITS_read_key(infits,TSTRING,"INSTMODE",instmode,NULL,&status);
  
  if (strncasecmp(instmode,"TTAG",4)){
    FITS_close_file(infits,&status);
    fclose(windefs_file);
    fclose(output_file);
    cf_if_error("Input file is not TTAG.");
  }
  
  FITS_read_key(infits,TDOUBLE,"EXPSTART",&expstart,NULL,&status);
  
  FITS_read_key(infits,TSTRING,"DAYNIGHT",daynight,NULL,&status);

  FITS_read_key(infits,TSTRING,"APERTURE",aperture,NULL,&status);


  FITS_movabs_hdu(infits,2,&hdutype,&status);
  FITS_read_key(infits,TSTRING,"TFORM1",tempstring,NULL,&status);
  sscanf(tempstring,"%ld%c",&nevents,&tempchar);

 
  
  
  hdu2_time=(double *)malloc(nevents*sizeof(double));
  hdu2_lambda=(double *)malloc(nevents*sizeof(double));
  hdu2_weight=(double *)malloc(nevents*sizeof(double));
  hdu2_channel=(char *)malloc(nevents*sizeof(char));
  hdu2_timeflgs=(char *)malloc(nevents*sizeof(char));
  hdu2_loc_flgs=(char *)malloc(nevents*sizeof(char));
 
  FITS_get_colnum(infits, TRUE, "TIME", &ncol, &status);
  FITS_read_col(infits, TDOUBLE, ncol, 1, 1, nevents, &intnull,
		hdu2_time, &anynull, &status);
 
  FITS_get_colnum(infits, TRUE, "LAMBDA", &ncol, &status);
  FITS_read_col(infits, TDOUBLE, ncol, 1, 1, nevents, &intnull,
		hdu2_lambda, &anynull, &status);

  if (do_flux) {
    FITS_get_colnum(infits, TRUE, "ERGCM2", &ncol, &status);
    FITS_read_col(infits, TDOUBLE, ncol, 1, 1, nevents, &intnull,
		  hdu2_weight, &anynull, &status);
  } else {
    FITS_get_colnum(infits, TRUE, "WEIGHT", &ncol, &status);
    FITS_read_col(infits, TDOUBLE, ncol, 1, 1, nevents, &intnull,
		  hdu2_weight, &anynull, &status);
  }

  FITS_get_colnum(infits, TRUE, "CHANNEL", &ncol, &status);
  FITS_read_col(infits, TBYTE, ncol, 1, 1, nevents, &intnull,
		hdu2_channel, &anynull, &status);

  FITS_get_colnum(infits, TRUE, "TIMEFLGS", &ncol, &status);
  FITS_read_col(infits, TBYTE, ncol, 1, 1, nevents, &intnull,
		hdu2_timeflgs, &anynull, &status);
  
  FITS_get_colnum(infits, TRUE, "LOC_FLGS", &ncol, &status);
  FITS_read_col(infits, TBYTE, ncol, 1, 1, nevents, &intnull,
		hdu2_loc_flgs, &anynull, &status);


  FITS_movabs_hdu(infits,3,&hdutype,&status);

  ngtis=cf_read_col(infits,TFLOAT,"START",(void **) &gti_start);
  ngtis=cf_read_col(infits,TFLOAT,"STOP",(void **) &gti_stop);

  gtis=(interval *)malloc(ngtis*sizeof(interval));

  for (i=0;i<ngtis;i++){
    gtis[i].start=gti_start[i];
    gtis[i].stop=gti_stop[i];
  }



  t=0.0;
  tmax=bins;
  i=0;
  ctot=0.0;
  badtime=0;  
  while (i<nevents){

    while (hdu2_time[i]>=tmax) {
      if (!(is_inside(t,gtis,ngtis))) badtime=1;
      if (!(is_inside(tmax,gtis,ngtis))) badtime=1;
      if (!badtime) fprintf(output_file,"%15lf %15lf\n",expstart+t/(3600.0*24.0),ctot/bins);
      t=tmax;
      tmax=t+bins;
      ctot=0.0;
      badtime=0;
    }

    if (!(is_good_time(hdu2_timeflgs[i],daynight))) badtime=1;

    if (is_inside(hdu2_lambda[i],swindows,nwindows)){
      if (hdu2_loc_flgs[i]==0){	
	if (LiForSiC==1){
	  if (is_in_lif_channel(hdu2_channel[i],aperture)){
	    ctot+=hdu2_weight[i];
	  }
	}
	else {
	  if (is_in_sic_channel(hdu2_channel[i],aperture)){
	    ctot+=hdu2_weight[i];
	  }
	}
      }       
    }
      
    i++;
  }
    
  free(hdu2_time);
  free(hdu2_lambda);
  free(hdu2_weight);
  free(hdu2_channel);
  free(hdu2_timeflgs);
  free(hdu2_loc_flgs);
 

       
  FITS_close_file(infits,&status);
  fclose(windefs_file);
  fclose(output_file);

  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return EXIT_SUCCESS;
  
}

int is_inside(float val,interval * interval_ptr, int nintervals){
  int i;
  i=0;
  while (i<nintervals){
    if ((val<=interval_ptr[i].stop)&&(val>=interval_ptr[i].start)){
      return 1;
    }
    i++;
  }
  return 0;
}

int is_good_time(char timeflag,char * daynight){

    if (!(strncasecmp(daynight,"DAY"  ,3))) return (timeflag  == 1);
    if (!(strncasecmp(daynight,"NIGHT",5))) return (timeflag  == 0);
    
    return ((timeflag&254) == 0);

}

int is_in_lif_channel(char channel, char * aperture){
  if (!strncasecmp(aperture,"LWRS",4)) return (channel == 3);
  if (!strncasecmp(aperture,"MDRS",4)) return (channel == 2);		 
  if (!strncasecmp(aperture,"HIRS",4)) return (channel == 1);
  return 0;
}

int is_in_sic_channel(char channel, char * aperture){
  if (!strncasecmp(aperture,"LWRS",4)) return (channel == 7);
  if (!strncasecmp(aperture,"MDRS",4)) return (channel == 6);		 
  if (!strncasecmp(aperture,"HIRS",4)) return (channel == 5);
  return 0;
}
