/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Synopsis:	getshift spec1 spec2 w1 w2
 *
 * Description:  Computes the wavelength shift between the spectra from two
 *               different exposures spec1 and spec2 using a cross-correlation
 *               in the wavelength window defined by w1 and w2  
 *
 *
 * History:   01/08/04  bjg  1.0  begin work 
 *            04/05/04  bjg  1.1  Return EXIT_SUCCESS
 *                                Correct formats to match
 *                                argument types in printf
 *            08/26/05  wvd  1.2  If given no arguments, don't return error.
 *				  Just print message and exit.
 *
 ****************************************************************************/

#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"


#define MAXSHIFT 25
#define SMOOTH_PAR 4


static char CF_PRGM_ID[] = "get_shift";
static char CF_VER_NUM[] = "1.2";


int wave_search(float *table, int i, int j, float elmt){

  int k;
  
  if (i>=j) return i;
  if (elmt==table[i]) return i;
  if (elmt==table[j]) return j;
  k=(i+j)/2;
  if (elmt==table[k]) return k;
  if ((elmt-table[i])*(elmt-table[k])<0) return wave_search(table,i,k,elmt);
  return wave_search(table,k+1,j,elmt);
}

int main(int argc, char *argv[])
{
   
  int      status=0;
  int      hdutype=0;

  long n,n1,n2,laux,i,j,k;
  
  float *wave, *flux1, *flux2, *flux1s, *flux2s;
  
  float w1,w2;
  
  fitsfile *infits1, *infits2;
  
  float correl[2*MAXSHIFT+1];
  
  
  /* Enter a timestamp into the log. */
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

  /* Initialize error checking */
  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
  
  
  /* check command line */
  if (argc != 5) {
    printf("Usage: getshift spec1 spec2 w1 w2 \n");
    return 1;
  }
  
  
  FITS_open_file(&infits1, argv[1], READONLY, &status);
  FITS_open_file(&infits2, argv[2], READONLY, &status);
  
  sscanf(argv[3],"%f",&w1);
  sscanf(argv[4],"%f",&w2);
  
	
  FITS_movabs_hdu(infits1, 2, &hdutype, &status);
  FITS_movabs_hdu(infits2, 2, &hdutype, &status);
  

  n=cf_read_col(infits1,TFLOAT,"WAVE",(void **) &wave);
  n=cf_read_col(infits1,TFLOAT,"FLUX",(void **) &flux1);
  n=cf_read_col(infits2,TFLOAT,"FLUX",(void **) &flux2);
  
  
  flux1s=(float*)malloc(n*sizeof(float));
  flux2s=(float*)malloc(n*sizeof(float));
  
  
  for (i=0;i<n;i++){
        flux1s[i]=0;
        flux2s[i]=0;
        k=0;
  	for (j=i-SMOOTH_PAR;j<=i+SMOOTH_PAR;j++){
  		if ((j>=0)&&(j<n)) {
  			flux1s[i]+=flux1[j];
                        flux2s[i]+=flux2[j];
                        k=k+1;
  		}
  	}
  	flux1s[i]/=k;
  	flux2s[i]/=k;
  	
  	
  }
  
 

  if ((w1-wave[0])*(w1-wave[n-1])>0) 
    cf_if_error("Wavelength parameters out of range\n");  

  if ((w2-wave[0])*(w2-wave[n-1])>0) 
    cf_if_error("Wavelength parameters out of range\n");
  
  n1=wave_search(wave,0,n-1,w1);
  n2=wave_search(wave,0,n-1,w2);
  
  
  if (n1>n2) {
    laux=n2;
    n2=n1;
    n1=laux;
  }
  
  if ((n1-MAXSHIFT<0)||(n2+MAXSHIFT>=n)) 
    cf_if_error("Wavelength parameters too close to bounds\n");
  
  k=0;
  
  for (i=0;i<=2*MAXSHIFT;i++){
    correl[i]=0;
    for (j=n1;j<=n2;j++){
      correl[i]=correl[i]+flux1s[j]*flux2s[j+i-MAXSHIFT];
    }
    if (correl[i]>correl[k]) k=i;    
  }
  
  
  
  FITS_close_file(infits1, &status);
  FITS_close_file(infits2, &status);

  /* Enter a timestamp into the log. */
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");

  printf("Shift is %ld steps.\n",k-MAXSHIFT);
  
  return EXIT_SUCCESS;
 

}
