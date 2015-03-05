
/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * 
 *
 * Usage:
 *     mjd2hjd [-h] [-v level] mjd ra_h ra_m ra_s dec_d dec_m dec_s
 *
 *
 * Arguments:
 *       mjd                               :  Modified Julian Day
 *       ra_h ra_m ra_s dec_d dec_m dec_s  :  RA and DEC of target : hh mm ss.s dd mm ss.s
 * 
 *
 * Options:
 *     -h:  this help message
 *     -v:  verbosity level (=1; 0 is silent)
 *
 *
 *
 *
 *
 * History:	10/06/04	bjg	v1.0	      
 *       	10/26/04	bjg	v1.1    Cosmetic change
 *
 ****************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "calfuse.h"



#define SECONDS_PER_DAY (24.0*3600.0)

static char CF_PRGM_ID[]= "mjd2hjd";
static char CF_VER_NUM[]= "1.1";


double gethmjd(double, int, int, float, int, int, float);

int main(int argc,char *argv[]){
  
  int optc;



  double     lighttime;
  double     jd, mjd, hjd, hmjd;

  int        rah,ram,decd,decm;
  float      ras,decs;
 

  char opts[] = "hv:";
  char usage[] = 
    "Usage:\n"
    "  mjd2hjd [-h] [-v level] mjd ra_h ra_m ra_s dec_d dec_m dec_s\n";
  char argument[] =
    "Arguments:\n"  
    "    mjd                               :  Modified Julian Day.\n"
    "    ra_h ra_m ra_s dec_d dec_m dec_s  :  RA and DEC of target.\n";
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
      return EXIT_SUCCESS;
    case 'v':
      verbose_level = atoi(optarg);
      break;
    }
  }

  /* Initialize error checking. */
  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
  if (argc != optind+7)
    cf_if_error("%s\nIncorrect number of program arguments", usage);
  
  
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  
  sscanf(argv[optind],"%lf",&mjd);
  sscanf(argv[optind+1],"%d",&rah);
  sscanf(argv[optind+2],"%d",&ram);  
  sscanf(argv[optind+3],"%f",&ras);
  sscanf(argv[optind+4],"%d",&decd);
  sscanf(argv[optind+5],"%d",&decm);
  sscanf(argv[optind+6],"%f",&decs);

 
  hmjd=gethmjd(mjd,rah,ram,ras,decd,decs,decm);

  jd   = mjd  + 2400000.5;
  hjd  = hmjd + 2400000.5;

  lighttime=(hmjd-mjd)*SECONDS_PER_DAY;

  cf_verbose(1,"Right Ascension of target: %d h %d m %f s",rah,ram,ras);
  cf_verbose(1,"Declination of target: %d deg %d m %f s",decd,decm,decs);
  cf_verbose(1,"jd= %lf",jd);
  cf_verbose(1,"mjd= %lf",mjd);
  cf_verbose(1,"hjd= %lf",hjd);
  cf_verbose(1,"hmjd= %lf",hmjd);
  cf_verbose(1,"heliocentric light time correction = %lf s",lighttime);
  cf_verbose(1,"---positive when Earth is closer to target than Sun");
  return EXIT_SUCCESS;

}

