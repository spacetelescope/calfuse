
/****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 ****************************************************************************
 *
 *
 * Usage:
 *    ttag_lightcurve_periodogram [-hs] [-v level] input_file output_file 
 *                                                minf maxf stepf 
 * 
 *
 * Arguments:
 *    input_file  : an ASCII file with 2 columns : 
 *                     - time (days) 
 *                     - countrate
 *
 *    output_file : an ASCII file with 2 columns : 
 *                     - frequency (Hz)
 *                     - normalized estimated power spectral density
 *   
 *    minf        : start frequency (Hz)
 *    maxf        :   end frequency (Hz)
 *    stepf       :  frequency step (Hz)
 *
 *
 *
 * Options:
 *     -h:  this help message
 *     -v:  verbosity level (=1; 0 is silent)
 *     -s:  time is in seconds instead of days
 *
 * 
 **************************************************************************
 *
 *
 *
 * Description:
 *
 *     Computes the Lomb periodogram of signal h
 *     from samples (hj) j=1..N at times (tj) j=1..N 
 * 
 *
 *      w=2*pi*f
 *
 *
 *
 *              1      /  (sum_j[(hj-H)*cos(w(tj-tau))])^2 
 * Pn(w) = ----------- | --------------------------------- 
 *          2*sigma2   \     sum_j[(cos(w(tj-tau)))^2]           
 *
 *
 *                              (sum_j[(hj-H)*sin(w(tj-tau))])^2 \
 *                        +    --------------------------------- | 
 *                                 sum_j[(cos(w(tj-tau)))^2]     /
 *
 *
 *
 *
 * H=sum_i[hi]/N       sigma2=sum_i[(hi-H)^2]/(N-1)
 * 
 *
 *                 sum_j[sin(2*w*tj)]  
 * tan (2*w*tau) = ------------------ 
 *                 sum_j[cos(2*w*tj)]
 *
 *
 * PAINFULLY SLOW
 * NEED TO BE REWRITTEN TO DO A FAST TRANSFORM (RECURSIVE)
 *
 *
 *
 *****************************************************************************
 *
 * History:	10/07/04	bjg	v1.0              
 *              11/04/04        bjg     v1.1 
 * 
 *                         ----v1.1 CHANGES WERE CANCELLED ---
 *              12/06/04        bjg     v1.2   
 *              12/13/04        bjg     v1.3 Changed N and i to long   
 *                                           Initialize N to 0  
 *		06/03/05	wvd	v1.4 Delete unused variables.
 * 
 *
 ***************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>


#include "calfuse.h"


static char CF_PRGM_ID[]= "ttag_lightcurve_periodogram";
static char CF_VER_NUM[]= "1.1";



int main(int argc,char *argv[]){
  
 
  
  char       input_filename[FLEN_CARD],output_filename[FLEN_CARD];
  
  FILE     * output_file;
  FILE     * input_file;

  float      minf,maxf;
  float      stepf;


  double     H;
  double     sigma2;
  double     tau;


  double     val1,val2;
  double   * t;
  double   * h;
  double     total1,total2,total3,total4;
  double     aux;
  double     Pn;
  double     f;
  double     j0=0;

  long        i, N;
  long ndx;

  int        time_unit=0;

  int optc;

  char opts[] = "hsv:";

  char usage[] = 
    "Usage:\n"
    "  ttag_lightcurve_periodogram [-hs] [-v level] input_file output_file\n"
    "                                              minf, maxf, stepf\n"
    "Arguments:\n"
    "  input_file  : an ASCII file with 2 columns :\n"
    "                   - time (JD or MJD, Helio- or Geo-centric)\n"
    "                   - countrate\n\n"
    "  output_file : an ASCII file with 2 columns :\n" 
    "                   - frequency (Hz)\n"
    "                   - normalized estimated power spectral density\n\n"
    "  minf        : start frequency (Hz)\n"
    "  maxf        :   end frequency (Hz)\n"
    "  stepf       :  frequency step (Hz)\n";

  char option[] =
    "Options:\n"
    "  -h:  this help message\n"
    "  -v:  verbosity level (=1; 0 is silent)\n"
    "  -s:  time is in seconds instead of days\n";

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
    case 's':
      time_unit=1;
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
  sscanf(argv[optind+2],"%f",&minf);   
  sscanf(argv[optind+3],"%f",&maxf);
  sscanf(argv[optind+4],"%f",&stepf);          
  
    
  
  if ((input_file = fopen (input_filename, "r")) == NULL){
    cf_if_error ("Unable to open file %s\n", input_filename);
  }


  N=0;
  while (fscanf(input_file,"%lf %lf",&val1,&val2)!=EOF){
    N++;
  }

  t=(double *)malloc(N*sizeof(double));
  h=(double *)malloc(N*sizeof(double));


  fseek(input_file,0,SEEK_SET);

  i=0;
  total1=0;
  while (fscanf(input_file,"%lf %lf",&val1,&val2)!=EOF){
    if (i==0) { 
      t[i]=0.0; j0=val1;
    }
    else {
      if (time_unit) t[i]=val1-j0;
      else t[i]=(val1-j0)*3600.0*24;
    }

    h[i]=val2;

    total1+=h[i];

    i++;

  }
  
  H=total1/N;

  total1=0;
  for (i=0;i<N;i++){
    aux=(h[i]-H);
    total1+=aux*aux;
  }

  sigma2=total1/(N-1);



  if ((output_file = fopen(output_filename,"w")) == NULL){
    fclose(input_file);
    cf_if_error ("Unable to create file %s\n", output_filename);    
  }
  
  ndx=0;
  for (f=minf; f<=maxf; f+=stepf){
    cf_verbose(4,"n=%ld;minf=%f;maxf=%f;freq=%lf",ndx,minf,maxf,f);
    ndx++;
    total1=0;
    total2=0;
    for (i=0;i<N;i++){
      total1+=sin(4*M_PI*f*t[i]);
      total2+=cos(4*M_PI*f*t[i]);
    }


    tau=atan(total1/total2)/(4*M_PI*f);



    total1=0;
    total2=0;
    total3=0;
    total4=0;

    for (i=0;i<N;i++){
      aux=cos(2*M_PI*f*(t[i]-tau));
      total1+=(h[i]-H)*aux;
      total2+=aux*aux;
      aux=sin(2*M_PI*f*(t[i]-tau));
      total3+=(h[i]-H)*aux;
      total4+=aux*aux;
    }

    Pn=(total1*total1/total2+total3*total3/total4)/(2*sigma2);
   
      
    fprintf(output_file,"%10lf  %10lf\n",f,Pn);

  }

  fclose(input_file);
  fclose(output_file);
  
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return EXIT_SUCCESS;
  
}


