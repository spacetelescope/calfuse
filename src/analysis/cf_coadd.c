/******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 ******************************************************************************
 *
 * Synopsis:  cf_coadd root_name
 *
 *	      This program calls the routines that produce the exposure-level
 *	      all files and the observation-level all, ano, and nvo files.
 *
 *            03/29/04  v0.1  bjg Begin work.
 *            04/05/04  v0.2  bjg Remove unused variables
 *                                Include unistd.h
 *            06/07/04  v0.3  bjg Updated some keywords
 *            08/12/04  v0.4  bjg Keywords QIK_COR and COMB_COR                      
 *                                now set to COMPLETE.
 *            08/26/04  v0.5  bjg Handle case aperture=RFPT 
 *            03/17/05  v0.6  tc  Rename EXPTIME to OBSTIME in all extensions
 *            08/22/05  v1.0  wvd Completely re-write for CalFUSE v3.1.7.
 *				  Calls cf_make_all_exp for each exposure and
 *				  cf_make_all_obs.csh only once.
 *            08/27/07  v1.1  bot Changed i from int to long to match n
 *				  Added L to corresponding values
 *
 *****************************************************************************/



#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_coadd";
static char CF_VER_NUM[] = "1.1";


/*********************************************************************/

int main(int argc, char *argv[])
{

  char asnf_file[120], filename[120], instmode[100], aper_string[100], aper_num[100];
  fitsfile   *asnf_fits_ptr, *infits;
  long n;
  int status = 0;
  long i;
  int anynull;
  FILE * ffile;
  

  char command0[256];


  char    memname[120], *memname_ptr;
  char    memtype[120], *memtype_ptr;
  char ** memname_ptr2;
  char ** memtype_ptr2;

  char NULL_STR[] = {'\0'};

  char file_present;
  char expo_num_string[10];
  int expo_num;

  char know=0;

  /* Enter a timestamp into the log. */
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);


  memname_ptr=memname;
  memtype_ptr=memtype;

  memname_ptr2=&memname_ptr;
  memtype_ptr2=&memtype_ptr;

  
  /* Input Parameters */
  if (argc != 2)
    cf_if_error("Usage: cf_coadd rootname");
  
  /*Open association file*/
  
  strcpy(asnf_file, argv[1]);
  strcat(asnf_file, "asnf.fit");
  FITS_open_file(&asnf_fits_ptr, asnf_file, READONLY, &status);
  FITS_movabs_hdu(asnf_fits_ptr, 2, NULL, &status);
  
  
  FITS_read_key(asnf_fits_ptr, TLONG, "NAXIS2", &n, NULL, &status);

  
  for (i=0L;i<n;i++){
    
    cf_if_fits_error(fits_read_col_str(asnf_fits_ptr, 1, i+1L, 1,1, NULL_STR,
				       memname_ptr2, &anynull, &status));
 
    cf_if_fits_error(fits_read_col_str(asnf_fits_ptr, 2, i+1L, 1,1, NULL_STR,
				       memtype_ptr2, &anynull, &status));

    cf_if_fits_error(fits_read_col_log(asnf_fits_ptr, 3, i+1L, 1, 1, 0,
				       &file_present, &anynull, &status));
    
    if (file_present == FALSE) continue;
    if (strncasecmp(memtype,"EXPOSURE",8)) continue;

    strncpy(expo_num_string,&memname[8],3);
    expo_num=atoi(expo_num_string);

    if (expo_num==701) continue;
        
    strcpy(filename,memname);
    strcat(filename,"1ahistfraw.fit");
    strcpy(instmode,"hist");
    ffile=fopen(filename,"r");
    if (ffile == NULL) {
      strcpy(instmode,"ttag");
      strcpy(filename,memname);
      strcat(filename,"1attagfraw.fit");
      ffile=fopen(filename,"r");
      if (ffile == NULL) {
	printf("Warning: Expo %s not present\n",expo_num_string);
	continue;
      }
      fclose(ffile);
    }
    else {
      fclose(ffile);
    }
    
    
    if (!know) {

      FITS_open_file(&infits, filename, READONLY, &status);
      FITS_read_key(infits, TSTRING, "APERTURE", aper_string, NULL, &status);

      if (!strncasecmp(aper_string,"LWRS",4)) strcpy(aper_num,"4");
      if (!strncasecmp(aper_string,"MDRS",4)) strcpy(aper_num,"2");
      if (!strncasecmp(aper_string,"HIRS",4)) strcpy(aper_num,"3");
      if (!strncasecmp(aper_string,"PINH",4)) strcpy(aper_num,"1");
      if (!strncasecmp(aper_string,"RFPT",4)) strcpy(aper_num,"4");

      FITS_close_file(infits, &status);
      know=1;
    }


    strcpy(command0,"cf_make_all_exp ");
    strcat(command0,memname);
    strcat(command0,"00all");
    strcat(command0,aper_num);
    strcat(command0,instmode);
    strcat(command0,"fcal.fit ");
    strcat(command0,memname);
    strcat(command0,"1alif");
    strcat(command0,aper_num);
    strcat(command0,instmode);
    strcat(command0,"fcal.fit ");
    strcat(command0,memname);
    strcat(command0,"1blif");
    strcat(command0,aper_num);
    strcat(command0,instmode);
    strcat(command0,"fcal.fit ");
    strcat(command0,memname);
    strcat(command0,"2alif");
    strcat(command0,aper_num);
    strcat(command0,instmode);
    strcat(command0,"fcal.fit ");
    strcat(command0,memname);
    strcat(command0,"2blif");
    strcat(command0,aper_num);
    strcat(command0,instmode);
    strcat(command0,"fcal.fit ");

    system(command0);

  }

  FITS_close_file(asnf_fits_ptr,&status);

  sprintf(command0,"cf_make_all_obs.csh %s", asnf_file);
  system(command0);

  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");

  return(EXIT_SUCCESS);

}
