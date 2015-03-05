
/*************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *************************************************************************
 *
 *
 * Usage:
 *     extract_jitter [-h] [-v level] input_file
 *                                    
 * 
 *
 * Arguments:
 *   input_file     : jitter file                                                    
 *
 *
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
 * History:  	11/04/2004  v1.0 bjg    
 *		06/03/2005  v1.1 wvd	Fix bug that set expstart=expend.
 *					Delete unused variables.
 *   
 ***************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "calfuse.h"


static char CF_PRGM_ID[]= "extract_jitter";
static char CF_VER_NUM[]= "1.0";


int main(int argc,char *argv[]){
  
 
  
  char       input_filename[FLEN_CARD];
  char       output_filename_x[FLEN_CARD]; 
  char       output_filename_y[FLEN_CARD];
  
  FILE     * output_file_x;
  FILE     * output_file_y;

  fitsfile * infits;

 

  int        intnull=0,anynull;
  int        status=0;
  int        hdutype=0;

  long       i, N;
  int        ncol;

  double *   hdu2_time;
  double *   hdu2_dx;
  double *   hdu2_dy;
  int    *   hdu2_trkflg;
  
 

  double mjd,dx,dy;

  int trkflg;

  double expstart, expend;

  int optc;
  
  char opts[] = "hfv:";
  char usage[] = 
    "Usage:\n"
    "  extract_jitter [-h] [-v level] input_file\n\n"
    "Arguments:\n"
    "  input_file     : jitter file (FITS)\n\n";
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
  if (argc != optind+1)
    cf_if_error("%s\nIncorrect number of program arguments", usage);
  
  
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  
  strcpy(input_filename,argv[optind]); 
  


  strcpy(output_filename_x,input_filename);
  strcpy(output_filename_y,input_filename);
  strcat(output_filename_x,".jtx");
  strcat(output_filename_y,".jty");
  

  /* Create the output file */
  output_file_x = fopen(output_filename_x,"w");
  output_file_y = fopen(output_filename_y,"w");
  


  FITS_open_file(&infits,input_filename,READONLY,&status);

  FITS_read_key(infits,TDOUBLE,"EXPSTART",&expstart,NULL,&status);
  FITS_read_key(infits,TDOUBLE,"EXPEND",&expend,NULL,&status);

  FITS_movabs_hdu(infits,2,&hdutype,&status);

  FITS_read_key(infits,TLONG,"NAXIS2",&N,NULL,&status);
 
  hdu2_time=(double *)malloc(N*sizeof(double));
  hdu2_dx=(double *)malloc(N*sizeof(double));
  hdu2_dy=(double *)malloc(N*sizeof(double));
  hdu2_trkflg=(int *)malloc(N*sizeof(int));
 
  FITS_get_colnum(infits, TRUE, "TIME", &ncol, &status);
  FITS_read_col(infits, TDOUBLE, ncol, 1, 1, N, &intnull,
		hdu2_time, &anynull, &status);
 
  FITS_get_colnum(infits, TRUE, "DX", &ncol, &status);
  FITS_read_col(infits, TDOUBLE, ncol, 1, 1, N, &intnull,
		hdu2_dx, &anynull, &status);

  FITS_get_colnum(infits, TRUE, "DY", &ncol, &status);
  FITS_read_col(infits, TDOUBLE, ncol, 1, 1, N, &intnull,
		hdu2_dy, &anynull, &status);
  
  FITS_get_colnum(infits, TRUE, "TRKFLG", &ncol, &status);
  FITS_read_col(infits, TINT, ncol, 1, 1, N, &intnull,
		hdu2_trkflg, &anynull, &status);
  

  for (i=0;i<N-100;i++){
    mjd=expstart+(hdu2_time[i])/(3600.0*24.0);
    dx=hdu2_dx[i];
    dy=hdu2_dy[i];
    trkflg=hdu2_trkflg[i];

    if ((mjd>=expstart) && (mjd<=expend) && (trkflg==5)){
      fprintf(output_file_x,"%15lf %15lf\n",mjd,dx);
      fprintf(output_file_y,"%15lf %15lf\n",mjd,dy);
    }
  }


    
  free(hdu2_time);
  free(hdu2_dx);
  free(hdu2_dy);
  free(hdu2_trkflg);
 

       
  FITS_close_file(infits,&status);
 
  fclose(output_file_x);
  fclose(output_file_y);


  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return EXIT_SUCCESS;
  
}
