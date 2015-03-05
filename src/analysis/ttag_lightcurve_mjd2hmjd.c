
/*************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *************************************************************************
 *
 *
 * Usage:
 *     ttag_lightcurve_mjd2hmjd [-h] [-v level] output_file input_file
 *                                      hh mm ss dd mm ss     
 *
 *
 * Arguments:
 *   input_file     : ASCII file with 2 columns : 
 *                        - time (MJD)        
 *                        - signal
 *                                        
 *   output_file    : ASCII file with 2 columns : 
 *                        - time (HMJD)      
 *                        - signal                       
 *
 *   hh mm ss dd mm ss : RA and DEC of target 
 *
 *
 *
 *
 * Options:
 *     -h:  this help message
 *     -v:  verbosity level (=1; 0 is silent)
 *          
 *
 *
 *
 * History:     10/06/04  bjg   v1.0
 *		06/03/05  wvd	v1.1	Delete unused variables.
 *   
 ***************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "calfuse.h"

double gethmjd(double, int, int, float, int, int, float);

static char CF_PRGM_ID[]= "ttag_lightcurve_mjd2hmjd";
static char CF_VER_NUM[]= "1.1";


int main(int argc,char *argv[]){
  
 
  
  char       input_filename[FLEN_CARD]; 
  char       output_filename[FLEN_CARD]; 

  
  FILE     * output_file;
  FILE     * input_file;


  double t,t2,v;
  int    rah, ram, decd, decm;
  float    ras, decs;


  int optc;
   
  char opts[] = "hv:";
  char usage[] = 
    "Usage:\n"
    "  ttag_lightcurve_mjd2hmjd [-h] [-v level] output_file input_file\n" 
    "                                   hh mm ss dd mm ss\n\n"
    "Arguments:\n"
    "  input_file     : ASCII file with 2 columns :\n"  
    "                     - time (MJD) \n"        
    "                     - signal\n\n"       
    "  output_file    : ASCII file with 2 columns : \n" 
    "                     - time (HMJD) \n"        
    "                     - signal \n\n"
    "  hh mm ss dd mm ss : RA and DEC of target\n";
  char option[] =
    "Options:\n"
    "  -h:  this help message                 \n"
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
    }
  }

  /* Initialize error checking. */
  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
  if (argc != optind+8)
    cf_if_error("%s\nIncorrect number of program arguments", usage);
  
  
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  
  strcpy(output_filename,argv[optind]); 
  strcpy(input_filename,argv[optind+1]);  
  sscanf(argv[optind+2],"%d",&rah);
  sscanf(argv[optind+3],"%d",&ram);
  sscanf(argv[optind+4],"%f",&ras);
  sscanf(argv[optind+5],"%d",&decd);
  sscanf(argv[optind+6],"%d",&decm);
  sscanf(argv[optind+7],"%f",&decs);


  if ((input_file = fopen(input_filename,"r")) == NULL){
    cf_if_error ("Unable to open file %s\n", input_filename);    
  }

  /* Create the output file */
  if ((output_file = fopen(output_filename,"w")) == NULL){
    fclose(input_file);
    cf_if_error ("Unable to create file %s\n", output_filename);    
  }
  
  while(fscanf(input_file,"%lf %lf",&t,&v)!=EOF) {
    t2=gethmjd(t, rah, ram, ras, decd, decm, decs);
   
    fprintf(output_file,"%15lf  %15lf\n",t2,v);
  }

  fclose(input_file);
  fclose(output_file);

  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return EXIT_SUCCESS;
  
}


