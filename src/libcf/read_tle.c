/******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 ******************************************************************************
 *
 * Synopsis:    read_tle(fitsfile *fptr)
 *
 * Description: read_tle will read the standard tle format file FUSE.TLE 
 *              and place the orbital elements into the header of the
 *              fitsfile.  The FUSE.TLE file will contain multiple orbital 
 *              element sets, so the set closest in time to the observation
 *              will be used.
 *
 * Arguments:   fitsfile    *fptr    Pointer to input file
 *
 * History:     11/11/98   emurphy   Begin work.
 *              03/18/98   emurphy   Added To do list.
 *              04/08/99       peb   Tidied code and added necessary
 *                                   include files
 *              06/07/99       peb   Added reporting of version control.
 *              06/22/99       peb   Added FITS_ wrappers.
 *              08/05/99       emm   Added statements to print date
 *                                   of obs and TLE if dtime > 5.
 *                                   Also modified program to work with
 *                                   FES keywords TEXPSTRT and TEXPEND
 *              05/31/00       peb   Implemented cfortran.h calls for slalib
 *                                   functions.
 *		02/13/03 v1.4  wvd   Change 0 to NULL in FITS_update_key
 *              03/01/03 v1.5  wvd   Correct use of pointer in FITS_read_key
 *              04/01/03 v1.6  wvd   Replace cf_errmsg with cf_if_warning,
 *					printf with cf_verbose
 *              12/18/03 v1.7  bjg   Change calfusettag.h to calfuse.h
 *              04/07/07 v1.8  wvd   Initialize min_mjd to zero to silence
 *					compiler warnings.
 *
 *****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef CFORTRAN
#include "cfortran.h"
PROTOCCALLSFSUB5(SLA_CLDJ, sla_cldj, INT, INT, INT, PDOUBLE, PINT)
#define slaCldj(IY, IM, ID, DJM, J) \
     CCALLSFSUB5(SLA_CLDJ, sla_cldj, INT, INT, INT, PDOUBLE, PINT, \
		 IY, IM, ID, DJM, J)
#else
#include "slalib.h"
#include "slamac.h"
#endif

#include "calfuse.h"

#define  MAXCHARS 120
static char CF_PRGM_ID[] = "read_tle";
static char CF_VER_NUM[] = "1.8";

void read_tle(fitsfile *fptr)
{
    char   line1[MAXCHARS], line2[MAXCHARS], line3[MAXCHARS];
    char   sat_num[10], security_class[10], international_num[10];
    char   inchar[15][15], sp[15][2], sgn[4][2];
    char   comment[FLEN_CARD], eccen_str[20];
    char   n6_mant_str[20], drag_mant_str[20], instrument[FLEN_CARD];
    int    status=0, hdutype=0, card_num, epoch_year, n6_exponent;
    int    drag_exponent, ephemeris, elset_num, checksum1, checksum2; 
    int    rev_num, e_month, e_day;
    double n6_mantissa, drag_mantissa, epoch_day, n2, inclin, raan;
    double aop, mean_anom, mean_mot, n6, drag, semimaj, revspersec;
    double expstart, expend, expmiddle, eccentr, dtime=1.0e9;
    double i_day, f_day, mjd, mjd_f, mjd_d, min_mjd=0;
    FILE   *ftle;

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin reading TLE file");

    /*  Open the file with the standard two-line elements. */
    ftle = NULL;
    ftle = fopen(cf_cal_file("FUSE.TLE"),"r");
    if (ftle == NULL)
	cf_if_error("Cannot find file FUSE.TLE, which contains the"
		  " orbital elements.");

    /*  Rewind the fits file and get the start time of the exposure */
    FITS_movabs_hdu(fptr, 1, &hdutype, &status);
    FITS_read_key(fptr, TSTRING, "INSTRUME", instrument, comment, &status);

    /* In the primary header, the FES files use keywords with a different
       name.  This is because the FES may have many exposures in a file.
       The TEXPSTRT is the start of the first exposure and TEXPEND is the
       end of the last exposure. */
    if (!strncmp(instrument,"FES",3)) {
         /* Read in FES keywords */
         FITS_read_key(fptr, TDOUBLE, "TEXPSTRT", &expstart, comment, &status);
         FITS_read_key(fptr, TDOUBLE, "TEXPEND", &expend, comment, &status);
    } else {
         /* Read in FUV keywords */
         FITS_read_key(fptr, TDOUBLE, "EXPSTART", &expstart, comment, &status);
         FITS_read_key(fptr, TDOUBLE, "EXPEND", &expend, comment, &status);
    }

    /* expmiddle is in mjd */
    expmiddle = (expstart+expend)/2.0;

    /*  Read in FUSE.TLE file by getting one line at a time  */
    while (fgets(line1, MAXCHARS,ftle) != NULL) {
	fgets(line2, MAXCHARS, ftle);
	fgets(line3, MAXCHARS, ftle);

#ifdef DEBUG
	printf("\n\n%-80.80s\n",line1);
	printf("%-80.80s\n",line2);
	printf("%-80.80s\n",line3);
#endif
	/*
	 *  This may look like a stupid way to do this, but it has to be
	 *  done this way.  You see, the fields can contain leading blanks
	 *  which should be interpreted as zeros.  However, c skips over
	 *  leading blanks when reading %f and %d.
	 */
	sscanf(line2,"%1c%1c%5c%1c%1c%8c%1c%2c%12c%1c%1c%9c%1c%1c%5c%2c%1c"
	       "%1c%5c%2c%1c%1c%1c%4c%1c",
	       inchar[1],sp[1],inchar[2],inchar[3],sp[2],inchar[4],sp[3],
	       inchar[5],inchar[6],sp[4],sgn[1],inchar[7],sp[5],sgn[2],
	       inchar[8],inchar[9],sp[6],sgn[3],inchar[10],inchar[11],
	       sp[7],inchar[12], sp[8],inchar[13],inchar[14]);

	inchar[1][1]='\0';
	card_num=atoi(inchar[1]);

	inchar[2][5]='\0';
	strncpy(sat_num,inchar[2],5);

	inchar[3][1]='\0';
	strncpy(security_class,inchar[3],1);

	inchar[4][8]='\0';
	strncpy(international_num,inchar[4],8);

	inchar[5][2]='\0';
	epoch_year=atoi(inchar[5]);
	if (epoch_year < 50)
	    epoch_year+=2000;
	else
	    epoch_year+=1900;

	inchar[6][12]='\0';
	epoch_day=atof(inchar[6]);

	inchar[7][9]='\0';
	n2=atof(inchar[7]);
	if (sgn[1][0] == '-')
	    n2*=-1.0;

	inchar[8][5]='\0';
	strcpy(n6_mant_str,"0.");
	strncat(n6_mant_str,inchar[8],7);
#ifdef DEBUG
	printf("n6=%-10.10s\n",n6_mant_str);
#endif
	n6_mantissa=atof(n6_mant_str);
	if (sgn[2][0] == '-')
	    n6_mantissa*=-1.0;

	inchar[9][2]='\0';
	n6_exponent=atoi(inchar[9]);

	inchar[10][5]='\0';
	strcpy(drag_mant_str,"0.");
	strncat(drag_mant_str,inchar[10],7);
#ifdef DEBUG
	printf("Drag=%-10.10s\n",drag_mant_str);
#endif
	drag_mantissa=atof(drag_mant_str);
	if (sgn[3][0] == '-')
	    drag_mantissa*=-1.0;

	inchar[11][2]='\0';
	drag_exponent=atoi(inchar[11]);

	inchar[12][1]='\0';
	ephemeris=atoi(inchar[12]);

	inchar[13][4]='\0';
	elset_num=atoi(inchar[13]);

	inchar[14][1]='\0';
	checksum1=atoi(inchar[14]);

	n6=n6_mantissa*pow(10.0,(double) n6_exponent);
	drag=drag_mantissa*pow(10.0,(double) drag_exponent);

#ifdef DEBUG
	printf("n6=>%10.5f %5d %15.6E\n", n6_mantissa, n6_exponent, n6);
	printf("drag=>%10.5f %5d %15.6E\n", drag_mantissa,
	       drag_exponent, drag);
	for(i=1; i<=14; i++)
	    printf("%2d %-15.15s\n",i,inchar[i]);
	for(i=1; i<=3; i++)
	    printf("%2d %-5.5s\n",i,sgn[i]);
#endif
	sscanf(line3,"%1c%1c%5c%1c%8c%1c%8c%1c%7c%1c%8c%1c%8c%1c%11c%5c%1c",
	       inchar[1],sp[1],inchar[2],sp[2],inchar[3],sp[3],inchar[4],
	       sp[4],inchar[5],sp[5],inchar[6],sp[6],inchar[7],sp[7],
	       inchar[8],inchar[9],inchar[10]);
	inchar[1][1]='\0';
	card_num=atoi(inchar[1]);

	inchar[2][5]='\0';
	strncpy(sat_num,inchar[2],5);

	inchar[3][8]='\0';
	inclin=atof(inchar[3]);

	inchar[4][8]='\0';
	raan=atof(inchar[4]);

	inchar[5][7]='\0';
	/* 
	 *  strcpy(eccen_str,"0.");
	 *  strncat(eccen_str,inchar[5],8);
	 */
	sprintf(eccen_str,"0.%-8s",inchar[5]);
#ifdef DEBUG
	printf("Ecstr=%-20.20s\n",eccen_str);
#endif
	eccentr=atof(eccen_str);

	inchar[6][8]='\0';
	aop=atof(inchar[6]);

	inchar[7][8]='\0';
	mean_anom=atof(inchar[7]);

	inchar[8][11]='\0';
	mean_mot=atof(inchar[8]);

	inchar[9][5]='\0';
	rev_num=atof(inchar[9]);

	inchar[10][1]='\0';
	checksum2=atoi(inchar[10]);

#ifdef DEBUG
	for(i=1; i<=10; i++)
	    printf("%2d %-15.15s\n",i,inchar[i]);
	printf("%4d\n%12.8f\n%10.8f\n%10.5E\n%10.5E%5d%5d%5d\n",
	       epoch_year, epoch_day, n2, n6, drag, ephemeris,
	       elset_num, checksum1);
	printf("\n%8.4f\n%8.4f\n%10.8f\n%8.4f\n%8.4f\n%11.8f\n%5d\n%1d\n",
	       inclin, raan, (float) eccentr, aop, mean_anom,
	       mean_mot, rev_num, checksum2);  
#endif
	/* Convert the year and day of year into a year month day */

	f_day = modf(epoch_day, &i_day);

	month_day(epoch_year, (int) i_day, &e_month, &e_day);

	/* Convert year month day into a Julian date. */
#ifdef CFORTRAN
	slaCldj(epoch_year, e_month, e_day, mjd, status);
#else
	slaCldj(epoch_year, e_month, e_day, &mjd, &status);
#endif

	mjd += f_day;

	if (fabs(expmiddle-mjd) < dtime) {
	    dtime = fabs(expmiddle-mjd);
            min_mjd=mjd;
	    revspersec = mean_mot/(24.0*60.0*60.0);
	    semimaj = pow(MU/(4.0*PI*PI*revspersec*revspersec),1.0/3.0);
            mjd_f=modf(mjd, &mjd_d);
            FITS_update_key(fptr, TDOUBLE, "EPCHTIMD", &mjd_d, NULL, &status);
            FITS_update_key(fptr, TDOUBLE, "EPCHTIMF", &mjd_f, NULL, &status);
	    FITS_update_key(fptr, TDOUBLE, "INCLINAT", &inclin, NULL, &status);
	    FITS_update_key(fptr, TDOUBLE, "ECCENTRY", &eccentr, NULL, &status);
	    FITS_update_key(fptr, TDOUBLE, "MEANANOM", &mean_anom, NULL,&status);
	    FITS_update_key(fptr, TDOUBLE, "ARGPERIG", &aop, NULL, &status);
	    FITS_update_key(fptr, TDOUBLE, "RASCASCN", &raan, NULL, &status);
	    FITS_update_key(fptr, TDOUBLE, "SEMIMAJR", &semimaj, NULL, &status);
	    FITS_update_key(fptr, TDOUBLE, "MEANMOTN", &mean_mot, NULL, &status);
	    FITS_update_key(fptr, TDOUBLE, "FDM2COEF", &n2, NULL, &status);
	    FITS_update_key(fptr, TDOUBLE, "SDM6COEF", &n6, NULL, &status);
	    FITS_update_key(fptr, TDOUBLE, "DRAGCOEF", &drag, NULL, &status);
	    FITS_update_key(fptr, TSTRING, "PROPMODL", "SGP4", NULL, &status);
	}
    }  /* endwhile */

    if (dtime > 5.0) {
        cf_verbose(2, "MJD OBS=%12.5f MJD TLE=%12.5f DT=%12.5f\n", 
                     expmiddle, min_mjd, dtime);
	cf_if_warning("Orbital elements are more than 5 days old.");
        }

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done reading TLE file");
}

