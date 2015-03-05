
/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Synopsis:	idf_combine outfile file1 file2 file3 file4 file5 ...
 *
 * Description: Reads multiple contiguous idf files and concatenates them 
 *              into a single file.
 *
 * History:	08/11/03    bjg		Begin work from ttag_combine.c
 *                                      sort the input files in time order
 *              08/12/03    bjg	 v1.0	First usable version
 *              08/13/03    bjg		Some minor changes
 *              08/14/03    bjg		Fixed some important keywords  
 *					in the main HDU
 *		08/20/03    bjg	 v1.1	Changed the way the binary tables
 *					are read, created and written
 *					(same as cf_ttag_init)
 *              09/09/03    bjg  v1.2   Fixed some keywords
 *              10/29/03    bjg         Changed the scaling format for
 *                                      output timeline to conform the 
 *                                      new IDF formatting.
 *              11/17/03    bjg  v1.3   Now warns users if the Focal
 *                                      Plane Assembly has moved between 
 *                                      exposures.
 *              11/22/03    bjg         Now prints the name of the
 *                                      included files in the main  
 *                                      header.
 *              12/10/03    bjg         Updated NSPEC keyword and 
 *					SPECxxx, WOFFLxxx, WOFFSxxx keywords
 *					Removed path in filenames written
 *					into header.
 *              02/25/04    bjg         Bug fix 
 *              03/25/04    bjg         Added TZERO and TSCALE values for 
 *                                      AIC and FEC countrates
 *              04/05/04    bjg         Remove unused variables
 *                                      Change formats to match arg types
 *                                      in printf
 *                                      Remove static keyword in struct key 
 *                                      definition
 *                                      Add braces in TSCAL and TZERO
 *                                      definitions
 *              06/08/04    bjg         Changes to handle EXP_STAT keyword
 *              08/13/04    bjg  v1.4   Recalculates Y centroid on demand
 *              10/06/04    bjg  v1.5   Added option -b to store 
 *                                      ORBITAL_VEL in a float 
 *              10/26/04    bjg  v1.6   Cosmetic change 
 *              01/28/2005  bjg  v1.7   Fixed crash within CFITSIO with
 *                                      with IDF files having empty tables
 *              02/16/2005  wvd  v1.8   Double TSCALE and TZERO 
 *                                      for the AIC_CNT_RATE array.
 *              03/22/2005  wvd  v1.9	Read and write TIME_SUNRISE
 *                                      and TIME_SUNSET as shorts.
 *              04/15/2005  tc   v1.10  Update PLANTIME keyword
 *              06/08/2005  tc   v1.11  Comment out warning if FPA has
 *					shifted.  Pipeline corrects for this.
 *              09/19/2005  wvd  v1.12  Don't even read FPA positions.
 *              05/22/2006  wvd  v1.13  Change most header keyword
 *					types from float to long.
 *              06/22/2006  wvd  v1.14  Add -z flag to force creation of
 *					an output file, even if empty.
 *              11/10/2006  wvd  v1.15  Add comment fields to WOFFS and
 *					WOFFL keywords.
 *              08/27/2007  bot  v1.16  Changed long for int l.138
 *		07/25/2008  wvd  v1.17  Set output EXP_STAT keyword to lowest
 *					non-negative value among input files.
 *              08/22/2008  wvd  v1.18  If -a flag is set, include files with
 *					positive values of EXP_STAT, but not
 *					those with negative values.
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

static char CF_PRGM_ID[]= "idf_combine";
static char CF_VER_NUM[]= "1.18";

int main(int argc,char *argv[]){

  char     date[FLEN_CARD]={'\0'};      
  char     rootname[FLEN_CARD];
  char     comment[FLEN_COMMENT];

  int      felem_hdu2,felem_hdu4,frow_hdu3;

  char     stime[FLEN_CARD],keyword[FLEN_CARD],card[FLEN_CARD];

  char     fmt_byte[FLEN_CARD],fmt_float[FLEN_CARD],fmt_short[FLEN_CARD]; 

  char *corrected_filename;
  char *string_pointer;
  char string0[100], string1[100];

  time_t   vtime;

  fitsfile *infits,*outfits;
	
  filename *filelist;
  double   *expstartlist;
  double   *expendlist;
  long     *neventslist;
  long     *nrecordslist;
  long     *ngtislist;
  
  double   delta_t;
  
  float    zero=0.0;
  
  filename tempstring;
  double   tempdouble;
  long     templong;
  float    tempfloat;
  char     tempchar;
  
  double   minexpstart;
  int      minindex;
  
  int      nfiles;

  int      intnull=0,anynull;
  int      ncol;
  
  int      status=0;
  int      hdutype=0;
  int      tref = 0;

  int n2;

  int exp_stat, exp_stat_out=99;

  long     nevents=0, nrecords=0, ngtis=0;
  long     n_real_events=0;
  
  long     i,j;
  
  float   *floatlist;
  short   *shortlist;
  char    *charlist;
  double  *doublelist;
  
  float    exptime=0, rawtime=0;
  long     neventscreened=0, neventscreenedpha=0, plantime=0, timehv=0;
  long     timescreened=0, timesaa=0, timelowlimbangle=0, timeburst=0,
		timejitter=0,timenight=0;

  /* float    fpasx0, fpasx, fpalx0, fpalx, dfpasx, dfpalx ; */
  /* int      fpa_split = FALSE; */

  int  hdu2_tfields=14;
  
  char hdu2_extname[]="TTAG DATA";   /* Name of this extension */
  
  char *hdu2_ttype[]={"TIME", "XRAW", "YRAW", "PHA", "WEIGHT", "XFARF", 
		      "YFARF",  "X", "Y", "CHANNEL", "TIMEFLGS",
		      "LOC_FLGS", "LAMBDA", "ERGCM2" };

  char *hdu2_tform[14];  /* We'll assign values when we know 
			    the number of elements in the data set. */


  char *hdu2_tunit[]={"SECONDS", "PIXELS", "PIXELS", "UNITLESS", "UNITLESS",
		      "PIXELS", "PIXELS", "PIXELS", "PIXELS", "UNITLESS",
		      "UNITLESS", "UNITLESS", "ANGSTROMS", "ERG CM^-2"};

  struct key hdu2_tscal[] = {{"TSCAL6", 0.25},  {"TSCAL7", 0.1}, 
			     {"TSCAL8", 0.25},  {"TSCAL9", 0.1}};
	
  struct key hdu2_tzero[] = {{"TZERO6", 8192.}, {"TZERO7", 0.}, 
			     {"TZERO8", 8192.}, {"TZERO9", 0.}};
  
  


  int  hdu3_tfields=2;
	
  char hdu3_extname[]="GTI";   /* Name of this extension */
  
  char *hdu3_ttype[]={"START", "STOP" };
  
  char *hdu3_tform[2]={"1D", "1D" };  /* We'll assign values when we know 
					 the number of elements in the data set. */


  char *hdu3_tunit[]={"seconds", "seconds"};

  struct key hdu4_tscal[16];
  struct key hdu4_tzero[16];
  
  char hdu4_extname[]="TIMELINE";
  int hdu4_tfields=16;  /* output table will have 16 columns */
  char *hdu4_ttype[]={"TIME", "STATUS_FLAGS", "TIME_SUNRISE", "TIME_SUNSET",
		      "LIMB_ANGLE", "LONGITUDE", "LATITUDE", "ORBITAL_VEL",
		      "HIGH_VOLTAGE", "LIF_CNT_RATE", "SIC_CNT_RATE",
		      "FEC_CNT_RATE", "AIC_CNT_RATE", "BKGD_CNT_RATE",
		      "YCENT_LIF","YCENT_SIC"};

  char *hdu4_tform[16]; /* we will define tform later, when the number
			   of photons is known */
  
  char *hdu4_tunit[]={"seconds", "unitless", "seconds", "seconds", "degrees",
		      "degrees", "degrees",  "km/s", "unitless", "counts/sec",
		      "counts/sec", "counts/sec", "counts/sec", "counts/sec",
		      "pixels","pixels" };


  int ignore_exp_stat=0;
  int do_ycent=0;
  int big_vel=0;
  int do_empty=FALSE, empty_output=FALSE;

  float    *x=NULL, *y=NULL, *weight=NULL;
  unsigned char *channel=NULL, *timeflags=NULL, *locflags=NULL;

  int optc;

  char opts[] = "hacbv:z";
  char usage[] = 
    "Usage:\n"
    "  idf_combine [-ahbcz]  [-v level] output_idf_file input_idf_files\n";
  char option[] =
    "Options:\n"
    "  -h:  this help message\n"
    "  -v:  verbosity level (=1; 0 is silent)\n"
    "  -a:  ignore EXP_STAT keyword \n"
    "  -c:  recalculates Y centroids (use target events) \n"
    "  -b:  store ORBITAL_VEL in a float \n"
    "  -z:  if no good data, generate empty output file\n";

  verbose_level = 1;

  /* Check number of options and arguments */
  while ((optc = getopt(argc, argv, opts)) != -1) {
    switch(optc) {
      
    case 'h':
      printf("%s\n%s", usage, option);
      return 0;
    case 'a':
      ignore_exp_stat=1;
      break;
    case 'v':
      verbose_level = atoi(optarg);
      break;
    case 'c':
      do_ycent=1;
      break; 
    case 'b':
      big_vel=1;
      break;
    case 'z':
      do_empty=TRUE;
      break;
    }
  }

  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

  if (argc < optind+2) {
    printf("%s", usage);
    return 1;
  }
   
    
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  
  /* get and display time */
  vtime = time(NULL) ;
  strcpy(stime,ctime(&vtime));
  
      
  nfiles=argc-optind-1;
  
  filelist     = (filename *)cf_calloc(nfiles, sizeof(filename));
  expstartlist = (double *)cf_calloc(nfiles, sizeof(double));
  expendlist   = (double *)cf_calloc(nfiles, sizeof(double));
  neventslist  = (long *)cf_calloc(nfiles, sizeof(long));
  nrecordslist = (long *)cf_calloc(nfiles, sizeof(long));
  ngtislist    = (long *)cf_calloc(nfiles, sizeof(long));
  
  if (!ignore_exp_stat)
    cf_verbose(1,"Will include only files with EXP_STAT=0") ;
  else
    cf_verbose(1,"Will include only files with EXP_STAT >= 0") ;

  cf_verbose(1,"GETTING INFORMATION ON INPUT FILES") ;

  n2=0;

  for (i=0; i<nfiles;i++) {
    
    FITS_open_file(&infits,argv[optind+i+1],READONLY,&status);

    FITS_read_key(infits,TINT,"EXP_STAT",&exp_stat,NULL,&status);

    strcpy(filelist[n2],argv[optind+i+1]);
    FITS_read_key(infits,TDOUBLE,"EXPSTART",&(expstartlist[n2]),NULL,&status);
    FITS_read_key(infits,TDOUBLE,"EXPEND",&(expendlist[n2]),NULL,&status);
    
    /* Set EXP_STAT to lowest non-negative value in input files. */
    if (exp_stat_out > exp_stat && exp_stat >= 0) exp_stat_out = exp_stat;

    /* Skip files with bad values of EXP_STAT. */
    if (!ignore_exp_stat && exp_stat != 0) {
	cf_if_warning("File %s rejected. EXP_STAT not equal to zero. Use -a flag to override.",
		argv[optind+i+1]) ;
	FITS_close_file(infits,&status);
	continue;
    }
    /* Never include files with negative values of EXP_STAT. */
    else if (exp_stat < 0) {
	cf_if_warning("File %s rejected. EXP_STAT less than zero.", argv[optind+i+1]) ;
	FITS_close_file(infits,&status);
	continue;
    }

    FITS_read_key(infits,TFLOAT,"EXPTIME",&(tempfloat),NULL,&status);
    exptime+=tempfloat;
    FITS_read_key(infits,TLONG,"NEVENTS",&templong,NULL,&status);
    n_real_events+=templong;
    FITS_read_key(infits,TFLOAT,"RAWTIME",&(tempfloat),NULL,&status);
    rawtime+=tempfloat;
    FITS_read_key(infits,TLONG,"PLANTIME",&(templong),NULL,&status);
    plantime+=tempfloat;
		
    FITS_read_key(infits,TLONG,"NBADEVNT",&templong,NULL,&status);
    neventscreened+=templong;
    FITS_read_key(infits,TLONG,"NBADPHA",&templong,NULL,&status);
    neventscreenedpha+=templong;
    
    FITS_read_key(infits,TLONG,"EXP_BAD",&templong,NULL,&status);
    timescreened+=templong;
    FITS_read_key(infits,TLONG,"EXP_BRST",&templong,NULL,&status);
    timeburst+=templong;
    FITS_read_key(infits,TLONG,"EXP_HV",&templong,NULL,&status);
    timehv+=templong;
    FITS_read_key(infits,TLONG,"EXP_JITR",&templong,NULL,&status);
    timejitter+=templong;
    FITS_read_key(infits,TLONG,"EXP_LIM",&templong,NULL,&status);
    timelowlimbangle+=templong;
    FITS_read_key(infits,TLONG,"EXP_SAA",&templong,NULL,&status);
    timesaa+=templong;
    FITS_read_key(infits,TLONG,"EXPNIGHT",&templong,NULL,&status);
    timenight+=templong;
    
    FITS_movabs_hdu(infits,2,&hdutype,&status);
    FITS_read_key(infits,TSTRING,"TFORM1",&tempstring,NULL,&status);
    sscanf(tempstring,"%ld%c",&neventslist[n2],&tempchar);
    nevents+=neventslist[n2];
    
    FITS_movabs_hdu(infits,3,&hdutype,&status);
    FITS_read_key(infits,TLONG,"NAXIS2",&(ngtislist[n2]),NULL,&status);
    ngtis+=ngtislist[n2];
    
    FITS_movabs_hdu(infits,4,&hdutype,&status);
    FITS_read_key(infits,TSTRING,"TFORM1",&tempstring,NULL,&status);
    sscanf(tempstring,"%ld%c",&nrecordslist[n2],&tempchar);
    nrecords+=nrecordslist[n2];
    
    FITS_close_file(infits,&status);
		
    n2++;
  }

  nfiles=n2;

  if (nfiles==0){
    if (do_empty) {
	cf_verbose(1,"No files to combine. Generating empty output file.");
	do_ycent = FALSE;
	empty_output = TRUE;
	expendlist[0]=expstartlist[0];
	nfiles = 1;
    }
    else {
	cf_verbose(1,"No files to combine. Exiting.");
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
	exit(0);
    }
  }
  
  cf_verbose(1,"SORTING INPUT FILES IN TIME ORDER") ;
  
  for (i=0; i<nfiles-1;i++) {
    minexpstart=expstartlist[i];
    minindex=i;
    for (j=i+1; j<nfiles;j++) {
      if (expstartlist[j]<minexpstart){
	minexpstart=expstartlist[j];
	minindex=j;
      }
    }
    
    strcpy(tempstring,filelist[minindex]);
    strcpy(filelist[minindex],filelist[i]);
    strcpy(filelist[i],tempstring);
    
    tempdouble=expstartlist[minindex];
    expstartlist[minindex]=expstartlist[i];
    expstartlist[i]=tempdouble;
    
    tempdouble=expendlist[minindex];
    expendlist[minindex]=expendlist[i];
    expendlist[i]=tempdouble;
    
    templong=neventslist[minindex];
    neventslist[minindex]=neventslist[i];
    neventslist[i]=templong;
		
    templong=nrecordslist[minindex];
    nrecordslist[minindex]=nrecordslist[i];
    nrecordslist[i]=templong;
    
    templong=ngtislist[minindex];
    ngtislist[minindex]=ngtislist[i];
    ngtislist[i]=templong;
  }

  cf_verbose(2,
	"N: Filename                  ExpStart ExpEnd  NEvents NGTIs NSeconds");

  for (i=0; i<nfiles;i++) {
    cf_verbose(2,"%ld: %s %7.1f  %7.1f %4ld  %2ld    %5ld",i,filelist[i],
	expstartlist[i],expendlist[i],neventslist[i],ngtislist[i],
	nrecordslist[i]);
    
  }
  
  cf_verbose(1,"CREATING OUTPUT FILE");
  
  FITS_open_file(&infits,filelist[0],READONLY,&status);
  FITS_create_file(&outfits,argv[optind],&status);
	
	
  cf_verbose(1,"WRITING MAIN HEADER");
  
  FITS_copy_hdu(infits,outfits,0,&status);
	
  FITS_read_key(infits,TSTRING,"ROOTNAME",rootname,NULL,&status);
  	
  rootname[8]='9';
  rootname[9]='9';
  rootname[10]='9';
  	
  FITS_update_key(outfits,TSTRING,"ROOTNAME",rootname,NULL,&status);

  string_pointer=strrchr(argv[optind],'/');
  if (string_pointer==NULL) corrected_filename=argv[optind];
  else corrected_filename=&(string_pointer[1]); 
  
  FITS_update_key(outfits,TINT,"EXP_STAT",&exp_stat_out,NULL,&status);

  FITS_update_key(outfits,TSTRING,"FILENAME",corrected_filename,NULL,&status);
	
  FITS_update_key(outfits,TSTRING,"EXP_ID","999",NULL,&status);
	
  fits_get_system_time(date, &tref, &status);
	
  FITS_update_key(outfits,TSTRING,"DATE",date,NULL,&status);
  	
  FITS_update_key(outfits,TDOUBLE,"EXPEND",&expendlist[nfiles-1],NULL,&status);	
  FITS_update_key(outfits,TFLOAT,"EXPTIME",&exptime,NULL,&status);
  FITS_update_key(outfits,TLONG,"NEVENTS",&n_real_events,NULL,&status);	
  FITS_update_key(outfits,TFLOAT,"RAWTIME",&rawtime,NULL,&status);
  FITS_update_key(outfits,TLONG,"PLANTIME",&plantime,NULL,&status);
  FITS_update_key(outfits,TLONG,"NBADEVNT",&neventscreened,NULL,&status);
  FITS_update_key(outfits,TLONG,"NBADPHA",&neventscreenedpha,NULL,&status);
  FITS_update_key(outfits,TLONG,"EXP_BAD",&(timescreened),NULL,&status);
  FITS_update_key(outfits,TLONG,"EXP_BRST",&(timeburst),NULL,&status);		
  FITS_update_key(outfits,TLONG,"EXP_HV",&(timehv),NULL,&status);
  FITS_update_key(outfits,TLONG,"EXP_JITR",&(timejitter),NULL,&status);
  FITS_update_key(outfits,TLONG,"EXP_LIM",&(timelowlimbangle),NULL,&status);
  FITS_update_key(outfits,TLONG,"EXP_SAA",&(timesaa),NULL,&status);
  FITS_update_key(outfits,TLONG,"EXPNIGHT",&(timenight),NULL,&status);
  
  if (do_ycent)   {
    strcpy(string0,"PERFORM");
    FITS_update_key(outfits, TSTRING, "YCNT_COR", string0, NULL, &status);

    for (i = 1; i < 8; i++) {
      if (i == 4) continue;
      sprintf(string1, "YQUAL%1ld", i);
      strcpy(string0,"HIGH");
      FITS_update_key(outfits, TSTRING, string1, string0, NULL, &status);
    }
  } 

  /* reset BPM_CAL */
  fits_read_key(outfits,TSTRING,"BPM_CAL",string0,NULL,&status);
  if (status!=0) {
    status = 0;
  }
  else {
    strcpy(string0,"");
    FITS_update_key(outfits,TSTRING,"BPM_CAL",string0,NULL,&status);
  }


  if (empty_output) {
     fits_write_history(outfits,"FILE CONTAINS NO DATA.",&status);
	n2=0;
     FITS_update_key(outfits,TINT,"NSPEC",&n2,NULL,&status);
  }
  else
     FITS_update_key(outfits,TINT,"NSPEC",&nfiles,NULL,&status);
  
  fits_write_history(outfits," COMBINED WITH IDF_COMBINE ",&status);
    
  for (i=0;i<nfiles;i++){

    string_pointer=strrchr(filelist[i],'/');
    if (string_pointer==NULL) corrected_filename=filelist[i];
    else corrected_filename=&(string_pointer[1]);  


    sprintf(keyword, "SPEC%.3ld", i+1);
    FITS_update_key(outfits,TSTRING,keyword,corrected_filename,NULL,&status);
    sprintf(keyword, "WOFFL%.3ld", i+1);
    sprintf(comment, "[Angstroms] Shift applied to LiF spectrum");
    FITS_update_key(outfits,TFLOAT,keyword,&zero,comment,&status);
    sprintf(keyword, "WOFFS%.3ld", i+1);
    sprintf(comment, "[Angstroms] Shift applied to SiC spectrum");
    FITS_update_key(outfits,TFLOAT,keyword,&zero,comment,&status);
  }
  
  FITS_close_file(infits,&status);
  
  cf_verbose(1,"PREPARING EVENTS LIST HDU");
  
  /* Generate the tform array */
  sprintf(fmt_byte,  "%ldB", nevents);
  sprintf(fmt_float, "%ldE", nevents);
  sprintf(fmt_short, "%ldI", nevents);
  
  hdu2_tform[0] = fmt_float;
  hdu2_tform[1] = fmt_short;
  hdu2_tform[2] = fmt_short;
  hdu2_tform[3] = fmt_byte;
  hdu2_tform[4] = fmt_float;
  for (i=5; i<9; i++)
    hdu2_tform[i] = fmt_short;
  for ( ; i<12; i++)
    hdu2_tform[i] = fmt_byte;
  hdu2_tform[12] = fmt_float;
  hdu2_tform[13] = fmt_float;
  

  /* Append a new empty binary table to the output file */
  FITS_create_tbl(outfits, BINARY_TBL, 1, hdu2_tfields, hdu2_ttype, hdu2_tform,
		  hdu2_tunit, hdu2_extname, &status);
  
  /* Write TSCALE and TZERO entries to header */
  for (i=0; i<4; i++) {
    sprintf(keyword, "TUNIT%ld", i+6);
    if (fits_read_keyword(outfits, keyword, card, NULL, &status))
      cf_if_fits_error(status);
    FITS_insert_key_flt(outfits, hdu2_tscal[i].keyword, hdu2_tscal[i].value,
			-2, NULL, &status);
    FITS_insert_key_flt(outfits, hdu2_tzero[i].keyword, hdu2_tzero[i].value,
			-4, NULL, &status);
  } 
  
  
  
  cf_verbose(1,"PREPARING GTIs LIST HDU");
  
  /* Append a new empty binary table to the output file */
  FITS_create_tbl(outfits, BINARY_TBL, ngtis, hdu3_tfields, hdu3_ttype, hdu3_tform,
		  hdu3_tunit, hdu3_extname, &status);
			
			
  cf_verbose(1,"PREPARING TIMELINE HDU");	
		
  /* Generate the tform array */
  sprintf(fmt_byte,  "%ldB", nrecords);
  sprintf(fmt_float, "%ldE", nrecords);
  sprintf(fmt_short, "%ldI", nrecords);
  

  
  hdu4_tform[0] = fmt_float;
  hdu4_tform[1] = fmt_byte;

  for (i=2; i<hdu4_tfields; i++) hdu4_tform[i] = fmt_short;
  
  /* Set TSCALE and TZERO values */
  /* Not all of these will be used; see below. */
  for (i=2; i<hdu4_tfields; i++) {
    sprintf(hdu4_tscal[i].keyword, "TSCAL%ld", i+1);
    hdu4_tscal[i].value = 0.1;
    sprintf(hdu4_tzero[i].keyword, "TZERO%ld", i+1);
    hdu4_tzero[i].value = 0.;
  }
  
  /* Let's keep three significant figures for vorb. */
  hdu4_tscal[7].value = 0.01;

    /* FEC_CNT_RATE and AIC_CNT_RATE can get big, so we
        use TZERO AND TSCALE to compress them. */
    hdu4_tscal[11].value = 1;                
    hdu4_tscal[12].value = 2;              
    hdu4_tzero[11].value = 32768;          
    hdu4_tzero[12].value = 65536;
  
  if (big_vel) hdu4_tform[7]=fmt_float; 
	
  /* Append a new empty binary table to the output file */
  FITS_create_tbl(outfits, BINARY_TBL, 1, hdu4_tfields, hdu4_ttype, hdu4_tform,
		  hdu4_tunit, hdu4_extname, &status);

	
  /* Write TSCALE and TZERO entries to header */
  for (i=2; i<hdu4_tfields; i++) {

    if ((i==7)&&(big_vel)) continue;
    
    /* Omit TSCALE and TZERO for these six arrays. */
    if (i ==  2) continue;	/* TIME_SUNRISE */
    if (i ==  3) continue;	/* TIME_SUNSET  */
    if (i ==  8) continue;	/* HIGH_VOLTAGE */
    if (i ==  9) continue;	/* LIF_CNT_RATE */
    if (i == 10) continue;	/* SIC_CNT_RATE */
    if (i == 13) continue;	/* BKGD_CNT_RATE */

    
    sprintf(keyword, "TUNIT%ld", i+1);
    if (fits_read_keyword(outfits, keyword, card, NULL, &status))
      cf_if_fits_error(status);
    FITS_insert_key_flt(outfits, hdu4_tscal[i].keyword, hdu4_tscal[i].value,
			-1, NULL, &status);
    FITS_insert_key_flt(outfits, hdu4_tzero[i].keyword, hdu4_tzero[i].value,
			-5, NULL, &status);
    
  } 
  
  
  
  
  
  
  felem_hdu2=1; frow_hdu3=1; felem_hdu4=1;
  
  
  cf_verbose(1,"COMBINING HDU");
  
  for (i=0;i<nfiles;i++){
    
    cf_verbose(2,"PROCESSING FILE %ld EVENTS LIST",i);
    
    delta_t=(expstartlist[i]-expstartlist[0])*3600*24;
    
    FITS_open_file(&infits,filelist[i],READONLY,&status);
    
    FITS_movabs_hdu(infits, 2, &hdutype, &status);
    FITS_movabs_hdu(outfits, 2, &hdutype, &status);

    if (neventslist[i]>0){
      floatlist=(float *)malloc(neventslist[i]*sizeof(float));
      shortlist=(short *)malloc(neventslist[i]*sizeof(short));
      charlist=(char *)malloc(neventslist[i]*sizeof(char));

      FITS_get_colnum(infits, TRUE, "TIME", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		    floatlist, &anynull, &status);
      for (j=0;j<neventslist[i];j++) floatlist[j]+=delta_t;	
      FITS_write_col(outfits, TFLOAT, 1, 1, felem_hdu2, neventslist[i], floatlist, &status);
    
    
    
      FITS_get_colnum(infits, TRUE, "XRAW", &ncol, &status);
      FITS_read_col(infits, TSHORT, ncol, 1, 1, neventslist[i], &intnull,
		    shortlist, &anynull, &status);
      FITS_write_col(outfits, TSHORT, 2, 1, felem_hdu2, neventslist[i], shortlist, &status);
    
    
    
      FITS_get_colnum(infits, TRUE, "YRAW", &ncol, &status);
      FITS_read_col(infits, TSHORT, ncol, 1, 1, neventslist[i], &intnull,
		    shortlist, &anynull, &status);
      FITS_write_col(outfits, TSHORT, 3, 1, felem_hdu2, neventslist[i], shortlist, &status);
		
    
    
      FITS_get_colnum(infits, TRUE, "PHA", &ncol, &status);
      FITS_read_col(infits, TBYTE, ncol, 1, 1, neventslist[i], &intnull,
		    charlist, &anynull, &status);
      FITS_write_col(outfits, TBYTE, 4, 1, felem_hdu2, neventslist[i], charlist, &status);
		

    
      FITS_get_colnum(infits, TRUE, "WEIGHT", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 5, 1, felem_hdu2, neventslist[i], floatlist, &status);
    

    
      FITS_get_colnum(infits, TRUE, "XFARF", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 6, 1, felem_hdu2, neventslist[i], floatlist, &status);
		

				
      FITS_get_colnum(infits, TRUE, "YFARF", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 7, 1, felem_hdu2, neventslist[i], floatlist, &status);
    
    
    
      FITS_get_colnum(infits, TRUE, "X", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 8, 1, felem_hdu2, neventslist[i], floatlist, &status);
    
    
    
      FITS_get_colnum(infits, TRUE, "Y", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 9, 1, felem_hdu2, neventslist[i], floatlist, &status);
    
    
    
      FITS_get_colnum(infits, TRUE, "CHANNEL", &ncol, &status);
      FITS_read_col(infits, TBYTE, ncol, 1, 1, neventslist[i], &intnull,
		    charlist, &anynull, &status);
      FITS_write_col(outfits, TBYTE, 10, 1, felem_hdu2, neventslist[i], charlist, &status);
    
    
    
      FITS_get_colnum(infits, TRUE, "TIMEFLGS", &ncol, &status);
      FITS_read_col(infits, TBYTE, ncol, 1, 1, neventslist[i], &intnull,
		    charlist, &anynull, &status);
      FITS_write_col(outfits, TBYTE, 11, 1, felem_hdu2, neventslist[i], charlist, &status);
    
    
    
      FITS_get_colnum(infits, TRUE, "LOC_FLGS", &ncol, &status);
      FITS_read_col(infits, TBYTE, ncol, 1, 1, neventslist[i], &intnull,
		    charlist, &anynull, &status);
      FITS_write_col(outfits, TBYTE, 12, 1, felem_hdu2, neventslist[i], charlist, &status);
    
    
    
      FITS_get_colnum(infits, TRUE, "LAMBDA", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 13, 1, felem_hdu2, neventslist[i], floatlist, &status);
    

    
      FITS_get_colnum(infits, TRUE, "ERGCM2", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 14, 1, felem_hdu2, neventslist[i], floatlist, &status);
    
    
      free(floatlist);
      free(charlist);
      free(shortlist);
    }
    
    cf_verbose(2,"PROCESSING FILE %ld GTIs LIST",i);
    
    FITS_movabs_hdu(infits, 3, &hdutype, &status);
    FITS_movabs_hdu(outfits, 3, &hdutype, &status);

    if (ngtislist[i]>0){
      doublelist=(double *)malloc(ngtislist[i]*sizeof(double));
    
      FITS_get_colnum(infits, TRUE, "START", &ncol, &status);
      FITS_read_col(infits, TDOUBLE, ncol, 1, 1, ngtislist[i], &intnull,
		    doublelist, &anynull, &status);
      for (j=0;j<ngtislist[i];j++) doublelist[j]+=delta_t;
      FITS_write_col(outfits, TDOUBLE, 1, frow_hdu3, 1, ngtislist[i], doublelist, &status);
    
      FITS_get_colnum(infits, TRUE, "STOP", &ncol, &status);
      FITS_read_col(infits, TDOUBLE, ncol, 1, 1, ngtislist[i], &intnull,
		    doublelist, &anynull, &status);
      for (j=0;j<ngtislist[i];j++) doublelist[j]+=delta_t;
      FITS_write_col(outfits, TDOUBLE, 2, frow_hdu3, 1, ngtislist[i], doublelist, &status);
    
    
      free(doublelist);
    }

    cf_verbose(2,"PROCESSING FILE %ld TIMELINE",i);
    
    FITS_movabs_hdu(infits, 4, &hdutype, &status);
    FITS_movabs_hdu(outfits, 4, &hdutype, &status);
    
    if (nrecordslist[i]>0){
      floatlist=(float *)malloc(nrecordslist[i]*sizeof(float));
      shortlist=(short *)malloc(nrecordslist[i]*sizeof(short));
      charlist=(char *)malloc(nrecordslist[i]*sizeof(char));
    
    
      FITS_get_colnum(infits, TRUE, "TIME", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      for (j=0;j<nrecordslist[i];j++) floatlist[j]+=delta_t;	
      FITS_write_col(outfits, TFLOAT, 1, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "STATUS_FLAGS", &ncol, &status);
      FITS_read_col(infits, TBYTE, ncol, 1, 1, nrecordslist[i], &intnull,
		    charlist, &anynull, &status);
      FITS_write_col(outfits, TBYTE, 2, 1, felem_hdu4, nrecordslist[i], charlist, &status);
    
      FITS_get_colnum(infits, TRUE, "TIME_SUNRISE", &ncol, &status);
      FITS_read_col(infits, TSHORT, ncol, 1, 1, nrecordslist[i], &intnull,
		    shortlist, &anynull, &status);
      FITS_write_col(outfits, TSHORT, 3, 1, felem_hdu4, nrecordslist[i], shortlist, &status);
    
      FITS_get_colnum(infits, TRUE, "TIME_SUNSET", &ncol, &status);
      FITS_read_col(infits, TSHORT, ncol, 1, 1, nrecordslist[i], &intnull,
		    shortlist, &anynull, &status);
      FITS_write_col(outfits, TSHORT, 4, 1, felem_hdu4, nrecordslist[i], shortlist, &status);
    
      FITS_get_colnum(infits, TRUE, "LIMB_ANGLE", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 5, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "LONGITUDE", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 6, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "LATITUDE", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 7, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "ORBITAL_VEL", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 8, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "HIGH_VOLTAGE", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 9, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "LIF_CNT_RATE", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 10, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "SIC_CNT_RATE", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 11, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
		
      FITS_get_colnum(infits, TRUE, "FEC_CNT_RATE", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 12, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "AIC_CNT_RATE", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 13, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "BKGD_CNT_RATE", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);	      
      FITS_write_col(outfits, TFLOAT, 14, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "YCENT_LIF", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 15, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
      FITS_get_colnum(infits, TRUE, "YCENT_SIC", &ncol, &status);
      FITS_read_col(infits, TFLOAT, ncol, 1, 1, nrecordslist[i], &intnull,
		    floatlist, &anynull, &status);
      FITS_write_col(outfits, TFLOAT, 16, 1, felem_hdu4, nrecordslist[i], floatlist, &status);
    
    
    
    
      free(floatlist);
      free(charlist);
      free(shortlist);
    }
    
    
    felem_hdu2+=neventslist[i];
    frow_hdu3+=ngtislist[i];
    felem_hdu4+=nrecordslist[i];
    FITS_close_file(infits,&status);
    
    
  }

/*
  if (fpa_split)
    cf_if_warning("FPA motion detected. Resolution may be compromised.");
*/

  if (do_ycent){
    FITS_movabs_hdu(outfits, 2, &hdutype, &status);
    nevents = cf_read_col(outfits, TFLOAT, "X",  (void **) &x);
    nevents = cf_read_col(outfits, TFLOAT, "Y",  (void **) &y);
    nevents = cf_read_col(outfits, TFLOAT, "WEIGHT", (void **) &weight);
    nevents = cf_read_col(outfits, TBYTE,  "TIMEFLGS", (void **) &timeflags);
    nevents = cf_read_col(outfits, TBYTE,  "LOC_FLGS", (void **) &locflags);
    nevents = cf_read_col(outfits, TBYTE,  "CHANNEL", (void **) &channel);
    
    FITS_movabs_hdu(outfits, 1, &hdutype, &status);
    
    /* Compute centroids. */
    cf_calculate_y_centroid(outfits, nevents, weight, x, y, channel, 
			    timeflags, locflags);
  } 
  
  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
  cf_verbose(1,"CLOSING OUTPUT FILE");
  
  FITS_close_file(outfits,&status);
  
 
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return EXIT_SUCCESS;
  
}


