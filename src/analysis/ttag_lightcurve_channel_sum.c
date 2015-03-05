
/*************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *************************************************************************
 *
 *
 * Usage:
 *     ttag_lightcurve_channel_sum [-h] [-v level] output_file 
 *                     input_file1 input_file2
 * 
 *
 * Arguments:
 *   input_files     : ASCII files with 2 columns : 
 *                        - time         
 *                        - signal
 *                                        
 *   output_file    : ASCII file with 2 columns : 
 *                        - time         
 *                        - signal                       
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
 * History:           10/06/04	bjg	v1.0
 *                    11/04/04	bjg	v1.1  time equality criteria changed
 *   
 ***************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "calfuse.h"


#define MAX_INPUT_FILE 20

static char CF_PRGM_ID[]= "ttag_lightcurve_channel_sum";
static char CF_VER_NUM[]= "1.1";


int is_equal(double a, double b, double delta){

  return ( ((a-b) <= delta) && ((a-b) > -delta) );

}


int main(int argc,char *argv[]){
  
 
  
  char       input_filename1[FLEN_CARD]; 
  char       input_filename2[FLEN_CARD];
  char       output_filename[FLEN_CARD]; 

  
  FILE     * output_file;
  FILE     * input_file1;
  FILE     * input_file2;

  double t1,v1,t2,v2, delta_t;

  int end_of_file;

  int optc;
   
  char opts[] = "hv:";
  char usage[] = 
    "Usage:\n"
    "  ttag_lightcurve_channel_sum [-h] [-v level] output_file\n" 
    "                  input_file1 input_file2\n\n"
    "Arguments:\n"
    "  input_files    : ASCII files with 2 columns :\n"  
    "                     - time\n"        
    "                     - signal\n\n"                                        
    "  output_file    : ASCII file with 2 columns : \n" 
    "                     - time \n"        
    "                     - signal\n";
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
  if (argc != optind+3)
    cf_if_error("%s\nIncorrect number of program arguments", usage);
  
  
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  
  strcpy(output_filename,argv[optind]); 
  strcpy(input_filename1,argv[optind+1]);  
  strcpy(input_filename2,argv[optind+2]);   



  if ((input_file1 = fopen(input_filename1,"r")) == NULL){
    cf_if_error ("Unable to open file %s\n", input_filename1);    
  }

  if ((input_file2 = fopen(input_filename2,"r")) == NULL){
    fclose(input_file1);
    cf_if_error ("Unable to open file %s\n", input_filename2);    
  }


  /* Create the output file */
  if ((output_file = fopen(output_filename,"w")) == NULL){
    fclose(input_file1);
    fclose(input_file2);
    cf_if_error ("Unable to create file %s\n", output_filename);    
  }
  
  fscanf(input_file1,"%lf %lf",&t1,&v1);
  fscanf(input_file1,"%lf %lf",&t2,&v2);
  delta_t=t2-t1;

  fseek(input_file1,0,SEEK_SET);
 
  
  end_of_file=(fscanf(input_file2,"%lf %lf",&t2,&v2)==EOF)||(fscanf(input_file1,"%lf %lf",&t1,&v1)==EOF);
  

  while (!end_of_file){
    if (is_equal(t2,t1,delta_t/2)){
      fprintf(output_file,"%15lf  %15lf\n",t1,v1+v2);
      end_of_file=(fscanf(input_file2,"%lf %lf",&t2,&v2)==EOF)||(fscanf(input_file1,"%lf %lf",&t1,&v1)==EOF);
    } 
    else if (t2<t1) {
      end_of_file=(fscanf(input_file2,"%lf %lf",&t2,&v2)==EOF);
    }
    else {
      end_of_file=(fscanf(input_file1,"%lf %lf",&t1,&v1)==EOF);
    }



  }
  

  fclose(input_file1);
  fclose(input_file2);
  fclose(output_file);

  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return EXIT_SUCCESS;
  
}


