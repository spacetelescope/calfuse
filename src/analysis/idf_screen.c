
/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Synopsis:	idf_screen infits outfits timeflag value
 *
 * 
 *
 * History:	10/05/04	bjg   1.0     	Begin work.
 *	        10/26/04        bjg   1.1     	Updates important header 
 *                                              keywords.
 *
 ****************************************************************************/


#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "calfuse.h"


typedef char filename[FLEN_CARD];

struct key {
  char keyword[FLEN_KEYWORD];
  float value;
};

static char CF_PRGM_ID[]= "idf_screen";
static char CF_VER_NUM[]= "1.1";

int main(int argc,char *argv[]){

  char       tempstring[FLEN_CARD]; 
  char       daynight[FLEN_CARD];
  char *     corrected_filename;
  char *     string_pointer;

  char time_test;
  
  int status=0; 
  int hdutype=0;

  long i, npts;

  char * tflags, * lflags;

  fitsfile *infits     , *outfits;

  filename input_fname , output_fname;

  char selstr1[FLEN_CARD];
  char selstr2[FLEN_CARD];

  char sel1;
  int sel2;
 
  long nbadevnt=0;
  long nbadpha=0;
 
  double exptime=0.0;
  double exp_bad=0.0;
  double exp_brst=0.0;
  double exp_hv=0.0;
  double exp_jitr=0.0;
  double exp_lim=0.0;
  double exp_saa=0.0;
  double expnight=0.0;

  int optc;

  char opts[] = "hv:";
  char usage[] = 
    "Usage:\n"
    "  idf_screen [-h]  [-v level] input_idf_file output_idf_file \n"
    "                              timeflag value\n";
  char argument[] = 
    "Arguments :\n"
    "  timeflag :  USER/JITR/OPUS/BRST/HV/SAA/LIMB/DAY \n"
    "  value    :  GOOD/BAD/EITHER           \n";
  char option[] =
    "Options:\n"
    "  -h:  this help message\n"
    "  -v:  verbosity level (=1; 0 is silent)\n";

  verbose_level = 1;

  /* Check number of options and arguments */
  while ((optc = getopt(argc, argv, opts)) != -1) {
    switch(optc) { 
    case 'h':
      printf("%s\n%s\n%s", usage, argument, option);
      return 0;
    case 'v':
      verbose_level = atoi(optarg);
      break;
    }
  }

 
  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

  if (argc != optind+4) {
    printf("%s", usage);
    cf_if_error("Incorrect number of arguments");
  }
    
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  

  strcpy(input_fname,argv[optind]);
  strcpy(output_fname,argv[optind+1]);
  strcpy(selstr1,argv[optind+2]);
  strcpy(selstr2,argv[optind+3]);  

  

  if (!strcasecmp(selstr2,"GOOD")) {
    sel2=0;
  } else if (!strcasecmp(selstr2,"BAD")) {
    sel2=1;
  } else if (!strcasecmp(selstr2,"EITHER")) {
    sel2=2;
  } else { 
    printf("%s", usage);
    cf_if_error("VALUE not recognized.");
  }


  if (!strcasecmp(selstr1,"USER")) {
    sel1=TEMPORAL_USER;
  } else if (!strcasecmp(selstr1,"JITR")) {
    sel1=TEMPORAL_JITR;
  } else if (!strcasecmp(selstr1,"OPUS")) {
    sel1=TEMPORAL_OPUS;
  } else if (!strcasecmp(selstr1,"BRST")) {
    sel1=TEMPORAL_BRST;
  } else if (!strcasecmp(selstr1,"HV")) {
    sel1=TEMPORAL_HV;
  } else if (!strcasecmp(selstr1,"SAA")) {
    sel1=TEMPORAL_SAA;
  } else if (!strcasecmp(selstr1,"LIMB")) {
    sel1=TEMPORAL_LIMB;
  } else if (!strcasecmp(selstr1,"DAY")) {
    sel1=TEMPORAL_DAY;
  } else { 
    printf("%s", usage);
    cf_if_error("TIMEFLAG not recognized.");
  }

  cf_verbose(1,"Input  Filename : %s",input_fname);
  cf_verbose(1,"Output Filename : %s",output_fname);
  cf_verbose(1,"FLAG : %s - %d",selstr1,sel1);
  cf_verbose(1,"VALUE : %s - %d",selstr2,sel2);

  if ((sel1>1)&&(sel2==0)){
    cf_verbose(1,"Nothing to do.");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
    return EXIT_SUCCESS;
  }

  cf_verbose(1,"OPENING INPUT FILE");
  FITS_open_file(&infits,input_fname,READONLY,&status);

  cf_verbose(1,"CREATING OUTPUT FILE");
  FITS_create_file(&outfits,output_fname,&status);
  FITS_copy_hdu(infits,outfits,0,&status);

  string_pointer=strrchr(output_fname,'/');
  if (string_pointer==NULL) corrected_filename=output_fname;
  else corrected_filename=&(string_pointer[1]);

  FITS_update_key(outfits,TSTRING,"FILENAME",corrected_filename,NULL,&status);

  strcpy(tempstring,"");

  FITS_update_key(outfits,TSTRING,"BPM_CAL",tempstring,NULL,&status);

  if (sel1 == 1) {
    switch(sel2) {
    case 0:
      strcpy(daynight,"NIGHT");
      break;
    case 1:
      strcpy(daynight,"DAY");      
      break;
    case 2:
      strcpy(daynight,"BOTH");     
      break;
    }
    FITS_update_key(outfits,TSTRING,"DAYNIGHT",daynight,NULL,&status);
  } else {
    FITS_read_key(infits,TSTRING,"DAYNIGHT",daynight,NULL,&status);
  }


  FITS_movabs_hdu(infits,2,&hdutype,&status);
  FITS_create_hdu(outfits,&status);
  FITS_copy_hdu(infits,outfits,0,&status);

  
  npts = cf_read_col(outfits, TBYTE, "TIMEFLGS", (void **) &tflags);  
  npts = cf_read_col(outfits, TBYTE, "LOC_FLGS", (void **) &lflags); 
  for (i=0;i<npts;i++){
    if (sel1>1){
      if (sel2==1){
	tflags[i]=(tflags[i])^sel1;
      } else if (sel2==2){
	tflags[i]=(tflags[i])&(~sel1);
      }
    }
	
    switch(daynight[0]) {
    case 'D':
      time_test=tflags[i]^TEMPORAL_DAY;
      break;
    case 'N':
      time_test=tflags[i];
      break;
    default:
      time_test=tflags[i]&(~TEMPORAL_DAY);
    }
      
    if (lflags[i]&LOCATION_PHA){
      nbadpha++;
    } 
    
    if ((lflags[i]&~LOCATION_AIR) || (time_test)) {
      nbadevnt++;
    }   
  }

  if ((sel1>1)&&(sel2>0)){
    cf_write_col(outfits, TBYTE, "TIMEFLGS", (void *) tflags, npts);  
  }

  free(tflags);
  free(lflags);


  FITS_movabs_hdu(infits,3,&hdutype,&status);
  FITS_create_hdu(outfits,&status);
  FITS_copy_hdu(infits,outfits,0,&status);
  

  FITS_movabs_hdu(infits,4,&hdutype,&status);
  FITS_create_hdu(outfits,&status);
  FITS_copy_hdu(infits,outfits,0,&status);

 
  npts = cf_read_col(outfits, TBYTE, "STATUS_FLAGS", (void **) &tflags);  
    
  for (i=0;i<npts;i++){
    if (sel1>1){
      if (sel2==1){
	tflags[i]=(tflags[i])^sel1;
      } else if (sel2==2){
	tflags[i]=(tflags[i])&(~sel1);
      }
    }

    switch(daynight[0]) {
    case 'D':
      time_test=tflags[i]^TEMPORAL_DAY;
      break;
    case 'N':
      time_test=tflags[i];
      break;
    default:
      time_test=tflags[i]&(~TEMPORAL_DAY);
    }

    
    if (time_test==0) {
      
      exptime+=1.0;
      
      if (!(tflags[i] & TEMPORAL_DAY))
	expnight+=1.0;
      
    } else {
      
      exp_bad+=1.0;
      
      if (tflags[i] & TEMPORAL_BRST)
	exp_brst+=1.0;
      
      if (tflags[i] & TEMPORAL_HV)
	exp_hv=+1.0;	
      
      if (tflags[i] & TEMPORAL_JITR)
	exp_jitr=+1.0;
      
      if (tflags[i] & TEMPORAL_LIMB)
	exp_lim=+1.0;
      
      if (tflags[i] & TEMPORAL_SAA)
	exp_saa=+1.0;
    }

  }

 if ((sel1>1)&&(sel2>0)){
  cf_write_col(outfits, TBYTE, "STATUS_FLAGS", (void *) tflags, npts);  
 }
 
 free(tflags);
  
 /* Update Main Header Keyword */
 FITS_movabs_hdu(outfits,1,&hdutype,&status);  
 
 FITS_update_key(outfits,TDOUBLE ,"EXPTIME"  ,&exptime   ,NULL,&status);	
 FITS_update_key(outfits,TLONG   ,"NBADEVNT" ,&nbadevnt  ,NULL,&status);
 FITS_update_key(outfits,TLONG   ,"NBADPHA"  ,&nbadpha   ,NULL,&status);
 FITS_update_key(outfits,TDOUBLE ,"EXP_BAD"  ,&exp_bad   ,NULL,&status);
 FITS_update_key(outfits,TDOUBLE ,"EXP_SAA"  ,&exp_saa   ,NULL,&status);
 FITS_update_key(outfits,TDOUBLE ,"EXP_LIM"  ,&exp_lim   ,NULL,&status);
 FITS_update_key(outfits,TDOUBLE ,"EXP_BRST" ,&exp_brst  ,NULL,&status);		
 FITS_update_key(outfits,TDOUBLE ,"EXP_JITR" ,&exp_jitr  ,NULL,&status);
 FITS_update_key(outfits,TDOUBLE ,"EXPNIGHT" ,&expnight  ,NULL,&status);

 cf_verbose(1,"CLOSING INPUT FILE");
 FITS_close_file(infits,&status);
 
 cf_verbose(1,"CLOSING OUTPUT FILE");
 FITS_close_file(outfits,&status);
  
 
 cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
 return EXIT_SUCCESS;
  
}


