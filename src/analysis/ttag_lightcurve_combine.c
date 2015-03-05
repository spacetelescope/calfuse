
/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Synopsis:	ttag_lightcurve_combine outfile file1 file2 file3 file4 ...
 *
 *
 * History:	09/07/04   bjg   v1.0
 *
 *
 ****************************************************************************/


#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "calfuse.h"


typedef char filename[FLEN_CARD];

static char CF_PRGM_ID[]= "ttag_lightcurve_combine";
static char CF_VER_NUM[]= "1.0";

int main(int argc,char *argv[]){

  char     stime[FLEN_CARD];
  time_t   vtime;

	
  filename *filelist;
  double   *expstartlist;
 
  filename tempstring;
  double   tempdouble;
  
  double   minexpstart;
  int      minindex;
  
  long     nfiles, n2;
  long     i,j;
  
  double   t,v;
  FILE   * fin;
  FILE   * fout;

  int optc;

  char opts[] = "hv:";
  char usage[] = 
    "Usage:\n"
    "  ttag_lightcurve_combine [-h]  [-v level] output_file input_files\n";
  char option[] =
    "Options:\n"
    "  -h:  this help message\n"
    "  -v:  verbosity level (=1; 0 is silent)\n";

  verbose_level = 1;

  /* Check number of options and arguments */
  while ((optc = getopt(argc, argv, opts)) != -1) {
    switch(optc) {
    case 'h':
      printf("%s\n%s", usage, option);
      return 0;
    case 'v':
      verbose_level = atoi(optarg);
      break;
    }
  }

  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

  if (argc < optind+2) {
    printf("%s", usage);
    cf_if_error("Incorrect number of arguments");
  }
   
    
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");
  
  /* get and display time */
  vtime = time(NULL) ;
  strcpy(stime,ctime(&vtime));
  
      
  nfiles=argc-optind-1;
  
  filelist     = (filename *)malloc(nfiles*sizeof(filename));
  expstartlist = (double *)malloc(nfiles*sizeof(double));


  cf_verbose(1,"GETTING INFORMATION ON INPUT FILES") ;

  n2=0;

  for (i=0; i<nfiles;i++) {

    fin=fopen(argv[optind+i+1],"r");

    if (fin==NULL){
      cf_if_warning("Could not open file %s. Skipped.",argv[optind+i+1]);
      continue;
    }
      
    fscanf(fin,"%lf %lf",&t,&v);

    expstartlist[n2]=t;
    strcpy(filelist[n2],argv[optind+i+1]);

    fclose(fin);
		
    n2++;
  }

  nfiles=n2;

  if (nfiles==0){
    cf_verbose(1,"No files to combine. Exiting.");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
    exit(0);
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
    
  }

  cf_verbose(2,"N: Filename                  ExpStart ");

  for (i=0; i<nfiles;i++) {
    cf_verbose(2,"%ld: %s %7.1f",i,filelist[i],expstartlist[i]);
  }
  
  
  cf_verbose(1,"CREATING OUTPUT FILE");
  if ((fout = fopen(argv[optind],"w")) == NULL){
    cf_if_error ("Unable to create file %s\n", argv[optind]);    
  }

  
  
  for (i=0; i<nfiles;i++) {
    fin=fopen(filelist[i],"r");
    while (fscanf(fin, "%lf %lf", &t, &v)!=EOF){
      fprintf(fout, "%15lf %15lf\n", t, v);
    }
    fclose(fin);
  }



  fclose(fout);
  
  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
  return EXIT_SUCCESS;
  
}


