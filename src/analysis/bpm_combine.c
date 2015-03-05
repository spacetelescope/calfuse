
/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSEbpm_combine
 *****************************************************************************
 *
 * Synopsis:	bpm_combine output_bpm_file combined_idf_file
 *		
 *
 * Description: Creates a bpm file associated with a combined idf file
 *              It gets the names of the idf files from the combined idf 
 *		file header. It gets the name of the bpm files from the idf 
 *		files header. These bpm files are then combined using the   
 *		offset information provided in the combined idf file header.  
 *		The BPM_CAL keyword is updated in the combined idf file 
 *		header.
 *
 * 		WARNING: all the single exposure IDF files (who took part in 
 *		the creation of the idf_file) and their associated BPM files
 *		must be in the working directory.
 *		
 *
 * History:	12/03/03	bjg	v1.0	Begin work.
 *              12/05/03        bjg             First version that compiles
 *              
 *		12/10/03 	bjg 	        Accept now only idf as parameter
 *                                              Updated NSPEC keyword and 
 *						SPECxxx, WOFFLxxx, WOFFSxxx
 *						keywords
 *						Removed paths in filenames 
 *						written into header.
 *              04/05/04        bjg             Remove unused variables
 *                                              Change formats to match arg
 *                                              types in printf
 *              05/25/04        bjg             Skip when BPM_CAL unpopulated 
 *                                              or not present.
 *              04/13/05	wvd	v2.0	combined_idf_file may be
 *						replaced by a file containing
 *						a list of BPM files.  The first
 *						line of this file must contain
 *						the number of entries that 
 *						follow.
 *
 ****************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "calfuse.h"

typedef char filename[FLEN_CARD];



static char CF_PRGM_ID[]= "bpm_combine";
static char CF_VER_NUM[]= "2.0";

int main(int argc,char *argv[]){

  char  date[FLEN_CARD]={'\0'};      
  char     rootname[FLEN_CARD];
  char *corrected_filename;
  char *string_pointer;

  int      felem_hdu2;

  char     stime[FLEN_CARD],keyword[FLEN_CARD];

 

  time_t   vtime;

  fitsfile *infits,*outfits,*idffits;
  char     *has_bpm_list;
  filename *filelist;
  double   *expstartlist;
  double   *expendlist;
  long     *neventslist;
  double   *sicshiftlist;
  double   *lifshiftlist;
  double   *exptimelist;
  
  double delta_t;
  
  filename tempstring,tempstring2;
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

  long     nevents=0;
  long     n_real_events=0;	
  long     i,j,istart;

  
  float   * xfield,*yfield,*weightfield,*lambdafield;
  char   *channelfield;
  
  char     fmt_byte[FLEN_CARD],fmt_float[FLEN_CARD],fmt_short[FLEN_CARD];  
  
  double   totalexptime=0, rawtime=0;
  long     neventscreened=0, neventscreenedpha=0;
  float    timescreened=0, timesaa=0, timelowlimbangle=0, timeburst=0, timejitter=0,timenight=0;

  

  int  hdu2_tfields=5;
  
  char hdu2_extname[]="POTHOLE_DATA";   /* Name of this extension */
  
  char *hdu2_ttype[]={"X", "Y", "CHANNEL", "WEIGHT", "LAMBDA"};

  char *hdu2_tform[5];  /* We'll assign values when we know 
			    the number of elements in the data set. */


  char *hdu2_tunit[]={"PIXELS", "PIXELS", "UNITLESS", "UNITLESS", "ANGSTROMS"};

  
  FILE *fp=NULL;
  char line[FLEN_FILENAME];
  int  maxline=FLEN_FILENAME;

  
  if (argc != 3) {
    printf("Incorrect number of arguments.\n");
    printf("Calling sequence:bpm_combine bpm_file combined_idf_file\n");
    printf("Final argument may be the name of a file containing a list of BPM files.\n");
    printf("First line must be number of BPM files in list.\n");
    exit(1);
  }
  
  /* Initialize error checking. */
  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  
  /* get and display time */
  vtime = time(NULL) ;
  strcpy(stime,ctime(&vtime));
  
  fits_open_file(&idffits,argv[2],READONLY,&status);
  if (status) {
	status = 0;
	if ((fp = fopen(argv[2], "r")) == NULL) {
	    printf("Can't open file %s\n", argv[2]);
	    return 1;
	}
	if ((fgets(line, maxline, fp)) == NULL) {
	    printf("Error reading file %s\n", argv[2]);
	    return 1;
	}
        sscanf(line,"%d",&nfiles);
  }
  else FITS_read_key(idffits,TINT,"NSPEC",&nfiles,NULL,&status);
  
  filelist     = (filename *)malloc(nfiles*sizeof(filename));
  neventslist  = (long *)malloc(nfiles*sizeof(long));     
  sicshiftlist  = (double *)calloc((size_t)nfiles,sizeof(double)); 
  lifshiftlist  = (double *)calloc((size_t)nfiles,sizeof(double)); 
  expstartlist = (double *)malloc(nfiles*sizeof(double));
  expendlist   = (double *)malloc(nfiles*sizeof(double));
  exptimelist   = (double *)malloc(nfiles*sizeof(double));
  has_bpm_list     = (char *)malloc(nfiles*sizeof(char));

  for (i=0; i<nfiles;i++) {

    has_bpm_list[i]=1;
    
    if (fp != NULL) {
	if((fgets(line, maxline, fp)) == NULL) {
	   printf("Error reading %s\n", argv[2]);
           return 1;
	}
        sscanf(line, "%s", filelist[i]);
    }
    else {
	/* Read name of individual IDF file from combined file header. */
      sprintf(tempstring, "SPEC%.3ld", i+1);
      FITS_read_key(idffits,TSTRING,tempstring,tempstring2,NULL,&status);
	/* Open individual IDF file and read name of associated BPM file. */
      FITS_open_file(&infits,tempstring2,READONLY,&status);
      fits_read_key(infits,TSTRING,"BPM_CAL",filelist[i],NULL,&status);
      if (status) {
        status=0;
        FITS_close_file(infits,&status);
        printf("BPM_CAL keyword not found in IDF %s; skipping.\n",tempstring2);
        has_bpm_list[i]=0;
        filelist[i][0]='\0';
        continue;
      }
	/* Close individual IDF file. */
      FITS_close_file(infits,&status);

	/* These keywords don't apply to a user-provided list of BPM files. */
      sprintf(tempstring, "WOFFL%.3ld", i+1);
      FITS_read_key(idffits,TDOUBLE,tempstring,&lifshiftlist[i],NULL,&status);
    
      sprintf(tempstring, "WOFFS%.3ld", i+1);
      FITS_read_key(idffits,TDOUBLE,tempstring,&sicshiftlist[i],NULL,&status);
    }

    /* Open individual BPM file. */
    fits_open_file(&infits,filelist[i],READONLY,&status);
    if (status) {
      status=0;
      printf("Could not open BPM %s for IDF %s; BPM_CAL keyword may be incorrectly populated; skipping.\n",filelist[i],tempstring2);
      has_bpm_list[i]=0;
      filelist[i][0]='\0';
      continue;
    }

    FITS_read_key(infits,TDOUBLE,"EXPSTART",&(expstartlist[i]),NULL,&status);
    FITS_read_key(infits,TDOUBLE,"EXPEND",&(expendlist[i]),NULL,&status);
   
    FITS_read_key(infits,TDOUBLE,"EXPTIME",&(exptimelist[i]),NULL,&status); 
    totalexptime+=exptimelist[i];

    FITS_read_key(infits,TDOUBLE,"RAWTIME",&(tempdouble),NULL,&status);
    rawtime+=tempdouble;
    
    FITS_read_key(infits,TLONG,"NBADEVNT",&templong,NULL,&status);
    neventscreened+=templong;
    FITS_read_key(infits,TLONG,"NBADPHA",&templong,NULL,&status);
    neventscreenedpha+=templong;
    
    FITS_read_key(infits,TFLOAT,"EXP_BAD",&tempfloat,NULL,&status);
    timescreened+=tempfloat;
    FITS_read_key(infits,TFLOAT,"EXP_SAA",&tempfloat,NULL,&status);
    timesaa+=tempfloat;
    FITS_read_key(infits,TFLOAT,"EXP_LIM",&tempfloat,NULL,&status);
    timelowlimbangle+=tempfloat;
    FITS_read_key(infits,TFLOAT,"EXP_BRST",&tempfloat,NULL,&status);
    timeburst+=tempfloat;
    FITS_read_key(infits,TFLOAT,"EXP_JITR",&tempfloat,NULL,&status);
    timejitter+=tempfloat;
    FITS_read_key(infits,TFLOAT,"EXPNIGHT",&tempfloat,NULL,&status);
    timenight+=tempfloat;
    
    FITS_read_key(infits,TFLOAT,"NEVENTS",&templong,NULL,&status);
    n_real_events+=templong;

    
    FITS_movabs_hdu(infits,2,&hdutype,&status);
    
    FITS_read_key(infits,TSTRING,"TFORM1",&tempstring,NULL,&status);
    sscanf(tempstring,"%ld%c",&neventslist[i],&tempchar);
    nevents+=neventslist[i];
    
    FITS_close_file(infits,&status);
    
  }
  
  if (fp != NULL) fclose(fp);
  else FITS_close_file(idffits,&status); 
  
  
  
  

  
  printf("---SORTING INPUT FILES IN TIME ORDER---\n") ;
  
  for (i=0; i<nfiles-1;i++) {
    if (!(has_bpm_list[i])) continue;
    minexpstart=expstartlist[i];
    minindex=i;
    for (j=i+1; j<nfiles;j++) {
      if (!(has_bpm_list[j])) continue;
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

    tempdouble=exptimelist[minindex];
    exptimelist[minindex]=exptimelist[i];
    exptimelist[i]=tempdouble;
    
    templong=neventslist[minindex];
    neventslist[minindex]=neventslist[i];
    neventslist[i]=templong;
		    
    tempdouble=lifshiftlist[minindex];
    lifshiftlist[minindex]=lifshiftlist[i];
    lifshiftlist[i]=tempdouble;
    
    tempdouble=sicshiftlist[minindex];
    sicshiftlist[minindex]=sicshiftlist[i];
    sicshiftlist[i]=tempdouble;
    
  }

  istart=-1;

  for (i=0; i<nfiles;i++) {
    if (!(has_bpm_list[i])) continue; 
      printf("%s %7.1f %7.1f %4ld %f %f\n",filelist[i],expstartlist[i],expendlist[i],neventslist[i],lifshiftlist[i],sicshiftlist[i]);
      if (istart==-1) istart=i;
  }
  printf("\n");

  if (istart==-1) {
    printf("No BPM files found. Exiting\n");
    exit(0);
  }
 
 
  
  
  printf("--------CREATING OUTPUT FILE-----\n");
  
  FITS_open_file(&infits,filelist[istart],READONLY,&status);
  FITS_create_file(&outfits,argv[1],&status);
	
  
	
	
  printf("--------WRITING MAIN HEADER-----\n");
  
  FITS_copy_hdu(infits,outfits,0,&status);
	
  FITS_read_key(infits,TSTRING,"ROOTNAME",rootname,NULL,&status);
	
  rootname[8]='9';
  rootname[9]='9';
  rootname[10]='9';
	
  FITS_update_key(outfits,TSTRING,"ROOTNAME",rootname,NULL,&status);

  string_pointer=strrchr(argv[1],'/');
  if (string_pointer==NULL) corrected_filename=argv[1];
    else corrected_filename=&(string_pointer[1]); 
	
  FITS_update_key(outfits,TSTRING,"FILENAME",corrected_filename,NULL,&status);
	
  FITS_update_key(outfits,TSTRING,"EXP_ID","999",NULL,&status);
	
  string_pointer=strrchr(argv[2],'/');
  if (string_pointer==NULL) corrected_filename=argv[2];
    else corrected_filename=&(string_pointer[1]); 
	

  if (fp == NULL) 
  FITS_update_key(outfits,TSTRING,"IDF_FILE",corrected_filename,NULL,&status);




  fits_get_system_time(date, &tref, &status);
	
  FITS_update_key(outfits,TSTRING,"DATE",date,NULL,&status);
	
  FITS_update_key(outfits,TDOUBLE,"EXPEND",&expendlist[nfiles-1],NULL,&status);	
  FITS_update_key(outfits,TDOUBLE,"EXPTIME",&totalexptime,NULL,&status);
  FITS_update_key(outfits,TDOUBLE,"RAWTIME",&rawtime,NULL,&status);	
  FITS_update_key(outfits,TLONG,"NEVENTS",&n_real_events,NULL,&status);	
  FITS_update_key(outfits,TLONG,"NBADEVNT",&neventscreened,NULL,&status);
  FITS_update_key(outfits,TLONG,"NBADPHA",&neventscreenedpha,NULL,&status);
  FITS_update_key(outfits,TFLOAT,"EXP_BAD",&(timescreened),NULL,&status);
  FITS_update_key(outfits,TFLOAT,"EXP_SAA",&(timesaa),NULL,&status);
  FITS_update_key(outfits,TFLOAT,"EXP_LIM",&(timelowlimbangle),NULL,&status);
  FITS_update_key(outfits,TFLOAT,"EXP_BRST",&(timeburst),NULL,&status);		
  FITS_update_key(outfits,TFLOAT,"EXP_JITR",&(timejitter),NULL,&status);
  FITS_update_key(outfits,TFLOAT,"EXPNIGHT",&(timenight),NULL,&status);

    fits_write_history(outfits," COMBINED WITH BPM_COMBINE ",&status);
    
  FITS_update_key(outfits,TINT,"NSPEC",&nfiles,NULL,&status);
  
  for (i=0;i<nfiles;i++){

    if (!(has_bpm_list[i])) continue; 
    
    string_pointer=strrchr(filelist[i],'/');
    if (string_pointer==NULL) corrected_filename=filelist[i];
    else corrected_filename=&(string_pointer[1]); 
    
    sprintf(keyword, "SPEC%.3ld", i+1);
    FITS_update_key(outfits,TSTRING,keyword,corrected_filename,NULL,&status);
    sprintf(keyword, "WOFFL%.3ld", i+1);
    FITS_update_key(outfits,TFLOAT,keyword,&lifshiftlist[i],NULL,&status);
    sprintf(keyword, "WOFFS%.3ld", i+1);
    FITS_update_key(outfits,TFLOAT,keyword,&sicshiftlist[i],NULL,&status);
  
      

  }



  
  FITS_close_file(infits,&status);
  
  printf("--------PREPARING EVENTS LIST HDU-----\n");
  
  /* Generate the tform array */
  sprintf(fmt_byte,  "%ldB", nevents);
  sprintf(fmt_float, "%ldE", nevents);
  sprintf(fmt_short, "%ldI", nevents);
  
  hdu2_tform[0] = fmt_float;
  hdu2_tform[1] = fmt_float;
  hdu2_tform[2] = fmt_byte;
  hdu2_tform[3] = fmt_float;
  hdu2_tform[4] = fmt_float;
  
  

  /* Append a new empty binary table to the output file */
  FITS_create_tbl(outfits, BINARY_TBL, 1, hdu2_tfields, hdu2_ttype, hdu2_tform,
		  hdu2_tunit, hdu2_extname, &status);
  
 
  
  
  
 
  
  printf("--------COMBINING HDU-----\n");
  felem_hdu2=1;
  for (i=0;i<nfiles;i++){
    if (!(has_bpm_list[i])) continue;
    
    
    delta_t=(expstartlist[i]-expstartlist[0])*3600*24;
    
    FITS_open_file(&infits,filelist[i],READONLY,&status);

   
    FITS_movabs_hdu(infits, 2, &hdutype, &status);
    FITS_movabs_hdu(outfits, 2, &hdutype, &status);
    printf("Processing file %s\n",filelist[i]);
    
    xfield=(float *)malloc(neventslist[i]*sizeof(float));
    yfield=(float *)malloc(neventslist[i]*sizeof(float));
    channelfield=(char *)malloc(neventslist[i]*sizeof(char));
    weightfield=(float *)malloc(neventslist[i]*sizeof(float));
    lambdafield=(float *)malloc(neventslist[i]*sizeof(float));
    
    
    FITS_get_colnum(infits, TRUE, "X", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		  xfield, &anynull, &status);
    
    
    FITS_get_colnum(infits, TRUE, "Y", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		  yfield, &anynull, &status);
   

    FITS_get_colnum(infits, TRUE, "CHANNEL", &ncol, &status);
    FITS_read_col(infits, TBYTE, ncol, 1, 1, neventslist[i], &intnull,
		  channelfield, &anynull, &status);
    

    FITS_get_colnum(infits, TRUE, "WEIGHT", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		  weightfield, &anynull, &status);
   
		
    FITS_get_colnum(infits, TRUE, "LAMBDA", &ncol, &status);
    FITS_read_col(infits, TFLOAT, ncol, 1, 1, neventslist[i], &intnull,
		  lambdafield, &anynull, &status);
    
    

    
    for (j=0;j<neventslist[i];j++){
	if (channelfield[j]<5) lambdafield[j]=lambdafield[j]+lifshiftlist[i];  
        if (channelfield[j]>4) lambdafield[j]=lambdafield[j]+sicshiftlist[i]; 
	weightfield[j]=weightfield[j]*exptimelist[i]/totalexptime;
    }


    
    FITS_write_col(outfits, TFLOAT, 1, 1, felem_hdu2, neventslist[i], xfield, &status);
    FITS_write_col(outfits, TFLOAT, 2, 1, felem_hdu2, neventslist[i], yfield, &status);
    FITS_write_col(outfits, TBYTE, 3, 1, felem_hdu2, neventslist[i], channelfield, &status);
    FITS_write_col(outfits, TFLOAT, 4, 1, felem_hdu2, neventslist[i], weightfield, &status);
    FITS_write_col(outfits, TFLOAT, 5, 1, felem_hdu2, neventslist[i], lambdafield, &status);
    
    
    free(xfield);
    free(yfield);
    free(channelfield);
    free(weightfield);
    free(lambdafield);
    
    felem_hdu2+=neventslist[i];
    
    FITS_close_file(infits,&status);
    
    
  }
  
  printf("--------CLOSING OUTPUT FILE-----\n");

  
  FITS_close_file(outfits,&status);
  
  
  if (fp == NULL) {
    FITS_open_file(&idffits,argv[2],READWRITE,&status);

    string_pointer=strrchr(argv[1],'/');
    if (string_pointer==NULL) corrected_filename=argv[1];
    else corrected_filename=&(string_pointer[1]);  

    FITS_update_key(idffits,TSTRING,"BPM_CAL",corrected_filename,NULL,&status);
    FITS_close_file(idffits,&status);
  }
  

  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return EXIT_SUCCESS;
  
}


