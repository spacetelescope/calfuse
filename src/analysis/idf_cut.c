
/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Usage:
 *   idf_cut [-hm]  [-v level] idf_file RefTime Period Nout
 * Arguments:
 *   RefTime :  Reference Time in seconds since EXPSTART
 *   Period  :  Period in seconds
 *   Nout    :  Number of output files.
 * Options:
 *   -m      :  interprets RefTime as MJD
 *   -h      :  this help message
 *   -v      :  verbosity level (=1; 0 is silent)
 *
 * Description: 
 *
 *     Cuts the IDF file into several smaller IDF files, sorting the 
 *   events according to their phase. The GTIs table and timeline get
 *   also cut. 
 *     From the input IDF file, idf_cut generates Nout IDF files named
 *   {input_idf_filename}.p{X}.fit with X=0..Nout-1. An event in the 
 *   input file is sorted in output file X if it happens in the interval
 *   [RefTime+(k+X/Nout)*Period,RefTime+(k+(X+1)/Nout)*Period[     
 *   with k an integer.
 *
 * History:	10/08/04	bjg  v1.0
 *              10/26/04        bjg  v1.1 Now updates the GTIs table
 *                                        Fixed EXPNIGHT computation
 *                                        Added -m switch.
 *                       		  + some cosmetic changes
 *              12/01/04        bjg  v1.2 Free allocated memory in the 
 *                                        subroutines
 *              12/02/04        bjg  v1.3 Fixed bug in make_gtitab
 *              06/03/05        wvd  v1.4 Change nevents and nseconds to longs.
 *              01/18/06        wvd  v1.5 Read arrays with cf_read_col.
 *
 ****************************************************************************/


#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "calfuse.h"

#define _MIN(a,b) (((a)<(b))?(a):(b))

typedef char filename[FLEN_CARD];

struct key {
  char keyword[FLEN_KEYWORD];
  float value;
};

static char CF_PRGM_ID[]= "idf_cut";
static char CF_VER_NUM[]= "1.1";


static int _isinphase(double t, double RefTime, double Period, double vinf, double vsup){
  
  float v;

  v=fmod(t-RefTime,Period);
  if (v<0) v+=Period; 

  return ((v>=vinf)&&(v<vsup));
}

static double _getphasenextstart(double t, double RefTime, double Period, double vinf, double vsup){

  float v;
  
  v=fmod(t-RefTime,Period);
  if (v<0) v+=Period; 
  
  if (v>vinf) return t-v+vinf+Period; 
  else return t-v+vinf;

}
static double _getphasenextstop(double t, double RefTime, double Period, double vinf, double vsup){

  float v;
  
  v=fmod(t-RefTime,Period);
  if (v<0) v+=Period; 
 
  if (v>vsup) return t-v+vsup+Period; 
  else return t-v+vsup;
}

struct _gti{
  double start;
  double stop;
  struct _gti *next;
};

typedef struct{
  struct _gti * gtilist;
  long N;
} _tempgtis;

static void _pushgti(double a, double b, _tempgtis * tempgtis){
  struct _gti *aux;
  
  aux=tempgtis->gtilist;
  
  tempgtis->gtilist=malloc(sizeof(struct _gti));
  tempgtis->N++;
  tempgtis->gtilist->start=a;
  tempgtis->gtilist->stop=b;
  tempgtis->gtilist->next=aux; 
} 

static void _popgti(double * a, double * b, _tempgtis * tempgtis){
  struct _gti *aux;
  
  *a=tempgtis->gtilist->start;
  *b=tempgtis->gtilist->stop;

  aux=tempgtis->gtilist;
  tempgtis->gtilist=tempgtis->gtilist->next;
  free(aux);
  tempgtis->N--;
} 

static _tempgtis * _makegtis(){
  _tempgtis * tempgtis;

  tempgtis=malloc(sizeof(_tempgtis));
  tempgtis->gtilist=NULL;
  tempgtis->N=0;
  return tempgtis;
}

static void _deletegtis(_tempgtis * tempgtis){
  double a,b;
  while (tempgtis->gtilist!=NULL) _popgti(&a,&b,tempgtis);
  free(tempgtis);
  tempgtis=NULL;
}

static long _ngtis(_tempgtis * tempgtis){
  return tempgtis->N;
}


/****************************************************************************/

long make_evlist(fitsfile *infits,fitsfile *outfits,
		 double RefTime,double Period, long Nout, 
		 long phase,char * daynight,
		 long * nbadevnt,long * nbadpha){
 
  
  int status=0;

  char pha;
  char loc_flgs;
  char timeflgs;
  char time_test;

  char     keyword[FLEN_CARD],card[FLEN_CARD];

  char     fmt_byte[FLEN_CARD],fmt_float[FLEN_CARD],fmt_short[FLEN_CARD]; 
  
  int  tfields=14;
  
  char extname[]="TTAG DATA";   /* Name of this extension */
  
  char *ttype[]={"TIME", "XRAW", "YRAW", "PHA", "WEIGHT", "XFARF", 
		 "YFARF",  "X", "Y", "CHANNEL", "TIMEFLGS",
		 "LOC_FLGS", "LAMBDA", "ERGCM2" };
  
  char *tform[14];  /* We'll assign values when we know 
		       the number of elements in the data set. */
  
  
  char *tunit[]={"SECONDS", "PIXELS", "PIXELS", "UNITLESS", "UNITLESS",
		 "PIXELS", "PIXELS", "PIXELS", "PIXELS", "UNITLESS",
		 "UNITLESS", "UNITLESS", "ANGSTROMS", "ERG CM^-2"};
  
  struct key tscal[] = {{"TSCAL6", 0.25},  {"TSCAL7", 0.1}, 
			{"TSCAL8", 0.25},  {"TSCAL9", 0.1}};
  
  struct key tzero[] = {{"TZERO6", 8192.}, {"TZERO7", 0.}, 
			{"TZERO8", 8192.}, {"TZERO9", 0.}};

  float *i_time;
  short *i_xraw;
  short *i_yraw;
  char *i_pha;
  float *i_weight;
  float *i_xfarf;
  float *i_yfarf;
  float *i_x;
  float *i_y;
  char *i_channel;
  char *i_timeflgs;
  char *i_loc_flgs;
  float *i_lambda;
  float *i_ergcm2;


  float *o_time;
  short *o_xraw;
  short *o_yraw;
  char *o_pha;
  float *o_weight;
  float *o_xfarf;
  float *o_yfarf;
  float *o_x;
  float *o_y;
  char *o_channel;
  char *o_timeflgs;
  char *o_loc_flgs;
  float *o_lambda;
  float *o_ergcm2;

  long npts;

  long i;

  float vinf,vsup;

  long nevents=0;

  *nbadevnt=0;
  *nbadpha=0;


  npts = cf_read_col(infits, TFLOAT, "TIME", (void **) &i_time);
  npts = cf_read_col(infits, TSHORT, "XRAW", (void **) &i_xraw);
  npts = cf_read_col(infits, TSHORT, "YRAW", (void **) &i_yraw);
  npts = cf_read_col(infits, TBYTE, "PHA", (void **) &i_pha);
  npts = cf_read_col(infits, TFLOAT, "WEIGHT", (void **) &i_weight);
  npts = cf_read_col(infits, TFLOAT, "XFARF", (void **) &i_xfarf);
  npts = cf_read_col(infits, TFLOAT, "YFARF", (void **) &i_yfarf);
  npts = cf_read_col(infits, TFLOAT, "X", (void **) &i_x);
  npts = cf_read_col(infits, TFLOAT, "Y", (void **) &i_y); 
  npts = cf_read_col(infits, TBYTE, "CHANNEL", (void **) &i_channel);
  npts = cf_read_col(infits, TBYTE, "TIMEFLGS", (void **) &i_timeflgs);
  npts = cf_read_col(infits, TBYTE, "LOC_FLGS", (void **) &i_loc_flgs);
  npts = cf_read_col(infits, TFLOAT, "LAMBDA", (void **) &i_lambda);
  npts = cf_read_col(infits, TFLOAT, "ERGCM2", (void **) &i_ergcm2); 


  o_time     = (float *) malloc(npts*sizeof(float));
  o_xraw     = (short *) malloc(npts*sizeof(short));
  o_yraw     = (short *) malloc(npts*sizeof(short));
  o_pha      = (char *)  malloc(npts*sizeof(char));
  o_weight   = (float *) malloc(npts*sizeof(float));
  o_xfarf    = (float *) malloc(npts*sizeof(float));
  o_yfarf    = (float *) malloc(npts*sizeof(float));
  o_x        = (float *) malloc(npts*sizeof(float));
  o_y        = (float *) malloc(npts*sizeof(float));
  o_channel  = (char *)  malloc(npts*sizeof(char));
  o_timeflgs = (char *)  malloc(npts*sizeof(char));
  o_loc_flgs = (char *)  malloc(npts*sizeof(char));
  o_lambda   = (float *) malloc(npts*sizeof(float));
  o_ergcm2   = (float *) malloc(npts*sizeof(float)); 

  vinf=(phase*Period)/((float)Nout);
  vsup=((phase+1)*Period)/((float)Nout);

  
  for (i=0;i<npts;i++){
 
    if (_isinphase(i_time[i], RefTime, Period, vinf, vsup)) {
      
      o_time[nevents]     = i_time[i];
      o_xraw[nevents]     = i_xraw[i];
      o_yraw[nevents]     = i_yraw[i]; 
      o_pha[nevents]      = i_pha[i];
      o_weight[nevents]   = i_weight[i];
      o_xfarf[nevents]    = i_xfarf[i];
      o_yfarf[nevents]    = i_yfarf[i];
      o_x[nevents]        = i_x[i];
      o_y[nevents]        = i_y[i];
      o_channel[nevents]  = i_channel[i]; 
      o_timeflgs[nevents] = i_timeflgs[i]; 
      o_loc_flgs[nevents] = i_loc_flgs[i];
      o_lambda[nevents]   = i_lambda[i];
      o_ergcm2[nevents]   = i_ergcm2[i];  

      pha=i_pha[i];
      loc_flgs=i_loc_flgs[i];
      timeflgs=i_timeflgs[i];

      switch(daynight[0]) {
      case 'D':
	time_test=timeflgs^TEMPORAL_DAY;
	break;
      case 'N':
	time_test=timeflgs;
	break;
      default:
	time_test=timeflgs&(~TEMPORAL_DAY);
      }

      if (loc_flgs&LOCATION_PHA){
	*nbadpha++;
      } 

      if ((loc_flgs&~LOCATION_AIR) || (time_test)) {
	*nbadevnt++;
      }

	
      nevents++;
     
    }
    

  }

  


  /* Generate the tform array */
  sprintf(fmt_byte,  "%ldB", nevents);
  sprintf(fmt_float, "%ldE", nevents);
  sprintf(fmt_short, "%ldI", nevents);
  
  tform[0] = fmt_float;
  tform[1] = fmt_short;
  tform[2] = fmt_short;
  tform[3] = fmt_byte;
  tform[4] = fmt_float;
  for (i=5; i<9; i++)
    tform[i] = fmt_short;
  for ( ; i<12; i++)
    tform[i] = fmt_byte;
  tform[12] = fmt_float;
  tform[13] = fmt_float;
  

  /* Append a new empty binary table to the output file */
  FITS_create_tbl(outfits, BINARY_TBL, 1, tfields, ttype, tform,
		  tunit, extname, &status);
  
  /* Write TSCALE and TZERO entries to header */
  for (i=0; i<4; i++) {
    sprintf(keyword, "TUNIT%ld", i+6);
    if (fits_read_keyword(outfits, keyword, card, NULL, &status))
      cf_if_fits_error(status);
    FITS_insert_key_flt(outfits, tscal[i].keyword, tscal[i].value,
			-2, NULL, &status);
    FITS_insert_key_flt(outfits, tzero[i].keyword, tzero[i].value,
			-4, NULL, &status);
  } 

  if (nevents > 0) {
   
    FITS_write_col(outfits, TFLOAT,  1, 1, 1, nevents, o_time, &status);
    FITS_write_col(outfits, TSHORT,  2, 1, 1, nevents, o_xraw, &status);
    FITS_write_col(outfits, TSHORT,  3, 1, 1, nevents, o_yraw, &status);
    FITS_write_col(outfits, TBYTE ,  4, 1, 1, nevents, o_pha, &status);
    FITS_write_col(outfits, TFLOAT,  5, 1, 1, nevents, o_weight, &status);
    FITS_write_col(outfits, TFLOAT,  6, 1, 1, nevents, o_xfarf, &status);
    FITS_write_col(outfits, TFLOAT,  7, 1, 1, nevents, o_yfarf, &status);
    FITS_write_col(outfits, TFLOAT,  8, 1, 1, nevents, o_x, &status);
    FITS_write_col(outfits, TFLOAT,  9, 1, 1, nevents, o_y, &status);
    FITS_write_col(outfits, TBYTE , 10, 1, 1, nevents, o_channel, &status);
    FITS_write_col(outfits, TBYTE , 11, 1, 1, nevents, o_timeflgs, &status);
    FITS_write_col(outfits, TBYTE , 12, 1, 1, nevents, o_loc_flgs, &status);
    FITS_write_col(outfits, TFLOAT, 13, 1, 1, nevents, o_lambda, &status);
    FITS_write_col(outfits, TFLOAT, 14, 1, 1, nevents, o_ergcm2, &status);
  }

  free(o_time);
  free(o_xraw );   
  free(o_yraw);    
  free(o_pha);     
  free(o_weight);   
  free(o_xfarf);    
  free(o_yfarf);    
  free(o_x);        
  free(o_y);        
  free(o_channel); 
  free(o_timeflgs);
  free(o_loc_flgs); 
  free(o_lambda);   
  free(o_ergcm2);

  free(i_time);
  free(i_xraw );   
  free(i_yraw);    
  free(i_pha);     
  free(i_weight);   
  free(i_xfarf);    
  free(i_yfarf);    
  free(i_x);        
  free(i_y);        
  free(i_channel); 
  free(i_timeflgs);
  free(i_loc_flgs); 
  free(i_lambda);   
  free(i_ergcm2); 
  
  return nevents;
}

/****************************************************************************/
long make_gtitab(fitsfile * infits, fitsfile * outfits, 
		 double RefTime, double Period, long Nout, 
		 long phase){
  
  

  int status=0;
  
  int  tfields=2;
  char extname[]="GTI";   /* Name of this extension */
  char *ttype[]={"START", "STOP" };
  char *tform[2]={"1D", "1D" };  
  char *tunit[]={"seconds", "seconds"};

  double *i_start;
  double *i_stop;
  double *o_start;
  double *o_stop;
  
  _tempgtis * tempgtis;

  long npts=0;
  long ngtis=0;
  long i;

  float vinf,vsup;
  
  double a,b;
  double _start, _stop;

  npts = cf_read_col(infits, TDOUBLE, "START", (void **) &i_start);
  npts = cf_read_col(infits, TDOUBLE, "STOP",  (void **) &i_stop);

  vinf=(phase*Period)/((float)Nout);
  vsup=((phase+1)*Period)/((float)Nout);

  tempgtis=_makegtis();
  for (i=0;i<npts;i++){
    cf_verbose(4,"Input GTI %ld = %lf/%lf",i,i_start[i],i_stop[i]);
    a=i_start[i];
    b=i_stop[i];
    if (_isinphase(a,RefTime, Period, vinf, vsup)) _start=a;
    else _start=_MIN(b,_getphasenextstart(a,RefTime, Period, vinf, vsup));
   
    while (_start<b){
      _stop=_MIN(b,_getphasenextstop(_start,RefTime, Period, vinf, vsup));
      _pushgti(_start,_stop,tempgtis);
      cf_verbose(4,"Output GTI = %lf,%lf",_start,_stop);
      _start=_getphasenextstart(_stop,RefTime, Period, vinf, vsup);
    }
  }

  ngtis=_ngtis(tempgtis);

  o_start=(double *)malloc(ngtis*sizeof(double));
  o_stop=(double *)malloc(ngtis*sizeof(double)); 

  for (i=ngtis-1;i>=0;i--){
    _popgti(&o_start[i],&o_stop[i],tempgtis);
  } 
  
  _deletegtis(tempgtis);
  
  /* Append a new empty binary table to the output file */
  FITS_create_tbl(outfits, BINARY_TBL, ngtis, tfields, ttype, tform,
		  tunit, extname, &status);
 

  if (ngtis>0){
    FITS_write_col(outfits, TDOUBLE, 1, 1, 1, ngtis, o_start, &status);
    FITS_write_col(outfits, TDOUBLE, 2, 1, 1, ngtis, o_stop, &status);
  }
  
  return ngtis;

}
/****************************************************************************/
long make_timeline(fitsfile *infits, fitsfile *outfits, 
		   double RefTime, double Period, long Nout, long phase, 
		   char *daynight, 
		   double *exptime, double *exp_bad, double *exp_brst, double *exp_hv, 
		   double *exp_jitr, double *exp_lim, double *exp_saa, double *expnight){

  char timeflgs;
  char time_test;
  int status=0;

  int bigtime=1;
  int big_vel=1;

  char     keyword[FLEN_CARD],card[FLEN_CARD];

  char     fmt_byte[FLEN_CARD],fmt_float[FLEN_CARD],fmt_short[FLEN_CARD]; 

  struct key tscal[16];
  struct key tzero[16];
  
  char extname[]="TIMELINE";

  int tfields=16;  /* output table will have 16 columns */

  char *ttype[]={"TIME", "STATUS_FLAGS", "TIME_SUNRISE", "TIME_SUNSET",
		 "LIMB_ANGLE", "LONGITUDE", "LATITUDE", "ORBITAL_VEL",
		 "HIGH_VOLTAGE", "LIF_CNT_RATE", "SIC_CNT_RATE",
		 "FEC_CNT_RATE", "AIC_CNT_RATE", "BKGD_CNT_RATE",
		 "YCENT_LIF","YCENT_SIC"};

  char *tform[16]; /* we will define tform later, when the number
		      of photons is known */
  
  char *tunit[]={"seconds", "unitless", "seconds", "seconds", "degrees",
		 "degrees", "degrees",  "km/s", "unitless", "counts/sec",
		 "counts/sec", "counts/sec", "counts/sec", "counts/sec",
		 "pixels","pixels" };  

  float *i_time;
  char  *i_status_flags;
  float *i_time_sunrise;
  float *i_time_sunset;
  float *i_limb_angle;
  float *i_longitude;
  float *i_latitude;
  float *i_orbital_vel;
  short *i_high_voltage;
  short *i_lif_cnt_rate;
  short *i_sic_cnt_rate;
  short *i_fec_cnt_rate;
  short *i_aic_cnt_rate;
  short *i_bkgd_cnt_rate;
  float *i_ycent_lif;
  float *i_ycent_sic;

  float *o_time;
  char  *o_status_flags;
  float *o_time_sunrise;
  float *o_time_sunset;
  float *o_limb_angle;
  float *o_longitude;
  float *o_latitude;
  float *o_orbital_vel;
  short *o_high_voltage;
  short *o_lif_cnt_rate;
  short *o_sic_cnt_rate;
  short *o_fec_cnt_rate;
  short *o_aic_cnt_rate;
  short *o_bkgd_cnt_rate;
  float *o_ycent_lif;
  float *o_ycent_sic;

  long npts;

  long i;

  float vinf,vsup;

  long nseconds=0;

  *exptime=0;
  *expnight=0;
  *exp_bad=0;
  *exp_brst=0;
  *exp_hv=0;
  *exp_jitr=0;
  *exp_lim=0;
  *exp_saa=0;
  *expnight=0;



  npts = cf_read_col(infits, TFLOAT, "TIME"          , (void **) &i_time);
  npts = cf_read_col(infits, TBYTE , "STATUS_FLAGS"  , (void **) &i_status_flags);
  npts = cf_read_col(infits, TFLOAT, "TIME_SUNRISE"  , (void **) &i_time_sunrise);
  npts = cf_read_col(infits, TFLOAT, "TIME_SUNSET"   , (void **) &i_time_sunset);
  npts = cf_read_col(infits, TFLOAT, "LIMB_ANGLE"    , (void **) &i_limb_angle);
  npts = cf_read_col(infits, TFLOAT, "LONGITUDE"     , (void **) &i_longitude);
  npts = cf_read_col(infits, TFLOAT, "LATITUDE"      , (void **) &i_latitude);
  npts = cf_read_col(infits, TFLOAT, "ORBITAL_VEL"   , (void **) &i_orbital_vel);
  npts = cf_read_col(infits, TSHORT, "HIGH_VOLTAGE"  , (void **) &i_high_voltage); 
  npts = cf_read_col(infits, TSHORT, "LIF_CNT_RATE"  , (void **) &i_lif_cnt_rate);
  npts = cf_read_col(infits, TSHORT, "SIC_CNT_RATE"  , (void **) &i_sic_cnt_rate);
  npts = cf_read_col(infits, TSHORT, "FEC_CNT_RATE"  , (void **) &i_fec_cnt_rate);
  npts = cf_read_col(infits, TSHORT, "AIC_CNT_RATE"  , (void **) &i_aic_cnt_rate);
  npts = cf_read_col(infits, TSHORT, "BKGD_CNT_RATE" , (void **) &i_bkgd_cnt_rate); 
  npts = cf_read_col(infits, TFLOAT, "YCENT_LIF"     , (void **) &i_ycent_lif);
  npts = cf_read_col(infits, TFLOAT, "YCENT_SIC"     , (void **) &i_ycent_sic);


  o_time            = (float *) malloc(npts*sizeof(float));
  o_status_flags    = (char *)  malloc(npts*sizeof(char));
  o_time_sunrise    = (float *) malloc(npts*sizeof(float));
  o_time_sunset     = (float *) malloc(npts*sizeof(float));
  o_limb_angle      = (float *) malloc(npts*sizeof(float));
  o_longitude       = (float *) malloc(npts*sizeof(float));
  o_latitude        = (float *) malloc(npts*sizeof(float));
  o_orbital_vel     = (float *) malloc(npts*sizeof(float));
  o_high_voltage    = (short *) malloc(npts*sizeof(short));
  o_lif_cnt_rate    = (short *) malloc(npts*sizeof(short));
  o_sic_cnt_rate    = (short *) malloc(npts*sizeof(short));
  o_fec_cnt_rate    = (short *) malloc(npts*sizeof(short));
  o_aic_cnt_rate    = (short *) malloc(npts*sizeof(short));
  o_bkgd_cnt_rate   = (short *) malloc(npts*sizeof(short));
  o_ycent_lif       = (float *) malloc(npts*sizeof(float));
  o_ycent_sic       = (float *) malloc(npts*sizeof(float));
 

  vinf=(phase*Period)/((float)Nout);
  vsup=((phase+1)*Period)/((float)Nout);


  for (i=0;i<npts;i++){
    
    if (_isinphase(i_time[i], RefTime, Period, vinf, vsup)) {

      o_time[nseconds]            = i_time[i];
      o_status_flags[nseconds]    = i_status_flags[i];
      o_time_sunrise[nseconds]    = i_time_sunrise[i];
      o_time_sunset[nseconds]     = i_time_sunset[i];
      o_limb_angle[nseconds]      = i_limb_angle[i];
      o_longitude[nseconds]       = i_longitude[i];
      o_latitude[nseconds]        = i_latitude[i];
      o_orbital_vel[nseconds]     = i_orbital_vel[i];
      o_high_voltage[nseconds]    = i_high_voltage[i];
      o_lif_cnt_rate[nseconds]    = i_lif_cnt_rate[i];
      o_sic_cnt_rate[nseconds]    = i_sic_cnt_rate[i];
      o_fec_cnt_rate[nseconds]    = i_fec_cnt_rate[i];
      o_aic_cnt_rate[nseconds]    = i_aic_cnt_rate[i];
      o_bkgd_cnt_rate[nseconds]   = i_bkgd_cnt_rate[i];
      o_ycent_lif[nseconds]       = i_ycent_lif[i];
      o_ycent_sic[nseconds]       = i_ycent_sic[i];

      timeflgs=i_status_flags[i];

      switch(daynight[0]) {
      case 'D':
	time_test=timeflgs^TEMPORAL_DAY;
	break;
      case 'N':
	time_test=timeflgs;
	break;
      default:
	time_test=timeflgs&(~TEMPORAL_DAY);
      }

    
      if (time_test==0) {

	*exptime+=1.0;

	if (!(timeflgs & TEMPORAL_DAY))
	  *expnight+=1.0;

      } else {

	*exp_bad+=1.0;
	
	if (timeflgs & TEMPORAL_BRST)
	  *exp_brst+=1.0;

	if (timeflgs & TEMPORAL_HV)
	  *exp_hv=+1.0;	

	if (timeflgs & TEMPORAL_JITR)
	  *exp_jitr=+1.0;

	if (timeflgs & TEMPORAL_LIMB)
	  *exp_lim=+1.0;

	if (timeflgs & TEMPORAL_SAA)
	  *exp_saa=+1.0;

      }
     
      nseconds++;
    }
  }
    
  /* Generate the tform array */
  sprintf(fmt_byte,  "%ldB", nseconds);
  sprintf(fmt_float, "%ldE", nseconds);
  sprintf(fmt_short, "%ldI", nseconds);

    
  tform[0] = fmt_float;
  tform[1] = fmt_byte;

  for (i=2; i<tfields; i++)	tform[i] = fmt_short;
  
  /* Set TSCALE and TZERO values */
  /* Not all of these will be used; see below. */
  for (i=2; i<tfields; i++) {
    sprintf(tscal[i].keyword, "TSCAL%ld", i+1);
    tscal[i].value = 0.1;
    sprintf(tzero[i].keyword, "TZERO%ld", i+1);
    tzero[i].value = 0.;
  }
  
  /* Need TZEROn = 3000. for tsunrise and tsunset. */
  tzero[2].value = tzero[3].value = 3000.;

  /* Let's keep three significant figures for vorb. */
  tscal[7].value = 0.01;

  /* FEC_CNT_RATE and AIC_CNT_RATE can get big, so we
     (essentially) store them as unsigned shorts. */
  tscal[11].value = 1;
  tscal[12].value = 1;
  tzero[11].value = 32768;
  tzero[12].value = 32768;
    
  
  if (bigtime) { 
    tform[2]=tform[3]=fmt_float; 
  }

  if (big_vel) tform[7]=fmt_float; 
    
  /* Append a new empty binary table to the output file */
  FITS_create_tbl(outfits, BINARY_TBL, 1, tfields, ttype, tform,
		  tunit, extname, &status);

  
	
  /* Write TSCALE and TZERO entries to header */
  for (i=2; i<tfields; i++) {

    /* If bigtime is TRUE, store tsunrise and tsunset as floats.*/
    if ((i == 2 || i == 3) && bigtime) continue;

    if ((i==7)&&(big_vel)) continue;
    
    /* Omit TSCALE and TZERO for these four arrays. */
    if (i ==  8) continue;	/* HIGH_VOLTAGE */
    if (i ==  9) continue;	/* LIF_CNT_RATE */
    if (i == 10) continue;	/* SIC_CNT_RATE */
    if (i == 13) continue;	/* BKGD_CNT_RATE */

      
    sprintf(keyword, "TUNIT%ld", i+1);
    if (fits_read_keyword(outfits, keyword, card, NULL, &status))
      cf_if_fits_error(status);
    FITS_insert_key_flt(outfits, tscal[i].keyword, tscal[i].value,
			-1, NULL, &status);
    FITS_insert_key_flt(outfits, tzero[i].keyword, tzero[i].value,
			-5, NULL, &status);
  } 

  if (nseconds > 0) { 

    FITS_write_col(outfits, TFLOAT,  1, 1, 1, nseconds, o_time          , &status);
    FITS_write_col(outfits, TBYTE ,  2, 1, 1, nseconds, o_status_flags  , &status);
    FITS_write_col(outfits, TFLOAT,  3, 1, 1, nseconds, o_time_sunrise  , &status);
    FITS_write_col(outfits, TFLOAT,  4, 1, 1, nseconds, o_time_sunset   , &status);
    FITS_write_col(outfits, TFLOAT,  5, 1, 1, nseconds, o_limb_angle    , &status);
    FITS_write_col(outfits, TFLOAT,  6, 1, 1, nseconds, o_longitude     , &status);
    FITS_write_col(outfits, TFLOAT,  7, 1, 1, nseconds, o_latitude      , &status);
    FITS_write_col(outfits, TFLOAT,  8, 1, 1, nseconds, o_orbital_vel   , &status);
    FITS_write_col(outfits, TSHORT,  9, 1, 1, nseconds, o_high_voltage  , &status);
    FITS_write_col(outfits, TSHORT, 10, 1, 1, nseconds, o_lif_cnt_rate  , &status);
    FITS_write_col(outfits, TSHORT, 11, 1, 1, nseconds, o_sic_cnt_rate  , &status);
    FITS_write_col(outfits, TSHORT, 12, 1, 1, nseconds, o_fec_cnt_rate  , &status);
    FITS_write_col(outfits, TSHORT, 13, 1, 1, nseconds, o_aic_cnt_rate  , &status);
    FITS_write_col(outfits, TSHORT, 14, 1, 1, nseconds, o_bkgd_cnt_rate , &status);
    FITS_write_col(outfits, TFLOAT, 15, 1, 1, nseconds, o_ycent_lif     , &status);
    FITS_write_col(outfits, TFLOAT, 16, 1, 1, nseconds, o_ycent_sic     , &status);

  } 

  free(o_time);
  free(o_status_flags);
  free(o_time_sunrise);  
  free(o_time_sunset);    
  free(o_limb_angle);    
  free(o_longitude);   
  free(o_latitude); 
  free(o_orbital_vel);
  free(o_high_voltage);
  free(o_lif_cnt_rate);
  free(o_sic_cnt_rate);
  free(o_fec_cnt_rate);
  free(o_aic_cnt_rate);
  free(o_bkgd_cnt_rate);
  free(o_ycent_lif);
  free(o_ycent_sic);

  free(i_time);
  free(i_status_flags);
  free(i_time_sunrise);  
  free(i_time_sunset);    
  free(i_limb_angle);    
  free(i_longitude);   
  free(i_latitude); 
  free(i_orbital_vel);
  free(i_high_voltage);
  free(i_lif_cnt_rate);
  free(i_sic_cnt_rate);
  free(i_fec_cnt_rate);
  free(i_aic_cnt_rate);
  free(i_bkgd_cnt_rate);
  free(i_ycent_lif);
  free(i_ycent_sic);

  return nseconds;
}
/****************************************************************************/


int main(int argc,char *argv[]){
  time_t   vtime;
  char     stime[FLEN_CARD];

  fitsfile * infits;
  fitsfile * outfits; 
 
  filename input_filename,output_filename, ifname;

  char string0[FLEN_CARD], daynight[FLEN_CARD];
  

  	
  long nevents;
  long nseconds;
  long ngtis;

  long nbadevnt;
  long nbadpha;
 
  double exptime;
  double exp_bad;
  double exp_brst;
  double exp_hv;
  double exp_jitr;
  double exp_lim;
  double exp_saa;
  double expnight;
 
  double RefTime;
  double Period;
  long Nout;

  double expstart;

  long i;

  int optc;
  
  char opts[] = "hmv:";
  char usage[] = 
    "Usage:\n"
    "  idf_cut [-hm]  [-v level] idf_file RefTime Period Nout\n";
  char argument[] =
    "Arguments:\n"
    "  RefTime :  Reference Time in seconds since EXPSTART\n"
    "  Period  :  Period in seconds\n"
    "  Nout    :  Number of output files.\n";
  char option[] =
    "Options:\n"
    "  -m:  interprets RefTime as MJD\n"
    "  -h:  this help message\n"
    "  -v:  verbosity level (=1; 0 is silent)\n";

  int hdutype=0;
  int status=0;
  int rtime_mjd=0;
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
    case 'm':
      rtime_mjd=1;
      break;
    }
  }

  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

  if (argc != optind+4) {
    printf("%s", usage);
    cf_if_error("Incorrect number of arguments");
  }

  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution."); 

  /* get and display time */
  vtime = time(NULL) ;
  strcpy(stime,ctime(&vtime));

  strcpy(input_filename,argv[optind]);

  sscanf(argv[optind+1],"%lf",&RefTime);
  sscanf(argv[optind+2],"%lf",&Period);
  sscanf(argv[optind+3],"%ld",&Nout);
  
  FITS_open_file(&infits,input_filename,READONLY,&status);
  FITS_read_key(infits,TSTRING,"DAYNIGHT",daynight,NULL,&status);
  FITS_read_key(infits,TSTRING,"FILENAME",ifname,NULL,&status);

  if (rtime_mjd){
    FITS_read_key(infits,TDOUBLE,"EXPSTART",&expstart,NULL,&status);
    RefTime=(RefTime-expstart)*3600.0*24.0;
  }

  for (i=0;i<Nout;i++){
    FITS_movabs_hdu(infits,1,&hdutype,&status);
    
    
    strcpy(output_filename,ifname);
    sprintf(string0,".p%ld.fit",i);
    strcat(output_filename,string0);

    FITS_create_file(&outfits,output_filename,&status);
    FITS_copy_hdu(infits,outfits,0,&status);

    FITS_update_key(outfits,TSTRING,"FILENAME",output_filename,NULL,&status);

    cf_verbose(1,"Phase: %ld/%ld",i+1,Nout);

    cf_verbose(2,"Creating Events List");
    /* Create Events List        */
    FITS_movabs_hdu(infits,2,&hdutype,&status);
    nevents=make_evlist(infits, outfits, RefTime, Period, Nout, i, daynight,
			&nbadevnt, &nbadpha);
    cf_verbose(2,"Number of Events: %ld",nevents);

    cf_verbose(2,"Creating GTI table");
    /* Create GTIs Table         */
    FITS_movabs_hdu(infits,3,&hdutype,&status);
    ngtis=make_gtitab(infits, outfits, RefTime, Period, Nout, i);
    cf_verbose(2,"Number of GTIs: %ld",ngtis);

    cf_verbose(2,"Creating Timeline");
    /* Create Timeline */
    FITS_movabs_hdu(infits,4,&hdutype,&status);
    nseconds=make_timeline(infits, outfits, RefTime, Period, Nout, i, daynight,
			   &exptime, &exp_bad, &exp_brst, &exp_hv, &exp_jitr,
			   &exp_lim, &exp_saa, &expnight);
    cf_verbose(2,"Number of records in timeline: %ld",nseconds);


    /* Update Main Header Keyword */
    FITS_movabs_hdu(outfits,1,&hdutype,&status);

   
    FITS_update_key(outfits,TDOUBLE ,"EXPTIME"  ,&exptime   ,NULL,&status);
    FITS_update_key(outfits,TLONG   ,"RAWTIME"  ,&nseconds  ,NULL,&status);	
    FITS_update_key(outfits,TLONG   ,"NEVENTS"  ,&nevents   ,NULL,&status);	
    FITS_update_key(outfits,TLONG   ,"NBADEVNT" ,&nbadevnt  ,NULL,&status);
    FITS_update_key(outfits,TLONG   ,"NBADPHA"  ,&nbadpha   ,NULL,&status);
    FITS_update_key(outfits,TDOUBLE ,"EXP_BAD"  ,&exp_bad   ,NULL,&status);
    FITS_update_key(outfits,TDOUBLE ,"EXP_SAA"  ,&exp_saa   ,NULL,&status);
    FITS_update_key(outfits,TDOUBLE ,"EXP_LIM"  ,&exp_lim   ,NULL,&status);
    FITS_update_key(outfits,TDOUBLE ,"EXP_BRST" ,&exp_brst  ,NULL,&status);		
    FITS_update_key(outfits,TDOUBLE ,"EXP_JITR" ,&exp_jitr  ,NULL,&status);
    FITS_update_key(outfits,TDOUBLE ,"EXPNIGHT" ,&expnight  ,NULL,&status);

    FITS_close_file(outfits,&status);
  }
  FITS_close_file(infits,&status);

  return (status);
}
