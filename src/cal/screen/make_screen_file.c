/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Synopsis:    make_screen_file
 *
 * Description: Make the screening parameters file.
 *
 * History:     07/14/98        emm     Begin work.
 *              07/08/99        peb     Added PHAMIN and PHAMAX keywords.
 *              07/29/99        emm     Added SAA_SCR and LIMB_SCR keywords
 ******************************************************************************/

#include <stdio.h>
#include "calfuse.h"
#define VERSION 3
#define EFFMJD 50000.0

int main()
{
char *channel[]={"LiF","SiC"};
char *segment[]={"1A","1B","2A","2B"};
char *segments[]={"1a","1b","2a","2b"};
char *aperture[]={"HIRS","MDRS","LWRS","PINH"};
char dummy[30];
char extname[]="SCREENING PARAMETERS";
char filename[80];


    int   i, status=0, vers;
    long  numbl, dumar[2]={0, 0};
    float numb, effmjd;
    fitsfile *parmfits;

for (i=0; i<4; i++) {  /* 1A,1B,2A,2B */

    sprintf(filename,"scrn%2.2s%03d.fit", segments[i],VERSION);
    printf("%20.20s\n",filename); 

    FITS_create_file(&parmfits,filename,&status);
    FITS_create_img(parmfits, SHORT_IMG, 0, dumar, &status);

    FITS_write_key(parmfits, TSTRING, "CALFTYPE", "SCRN", 
		   "Calibration file type", &status);

    vers=VERSION;
    FITS_write_key(parmfits,TINT,"CALFVERS",&vers,
		   "Calibration file version", &status);

    FITS_write_key(parmfits,TSTRING,"DETECTOR",segment[i],
                   "detector (1A, 1B, 2A, 2B", &status);

    effmjd=EFFMJD;
    FITS_write_key(parmfits,TFLOAT,"EFFMJD",&effmjd,
		   "Date on which file should be applied (MJD)", &status);

    FITS_write_date(parmfits, &status);

    FITS_write_key(parmfits,TSTRING,"AUTHOR","EDWARD MURPHY",
		   "Author of file", &status);

    numb=10.0;
    FITS_write_key(parmfits,TFLOAT,"TIMESTEP",&numb,
		   "[sec] Screening time step", &status);

    FITS_write_key(parmfits,TSTRING,"SAA_SCR","ON",
		   "SAA screening ON/OFF", &status);

    FITS_write_key(parmfits,TSTRING,"LIMB_SCR","ON",
		   "Limb angle screening ON/OFF", &status);

    FITS_write_key(parmfits,TSTRING,"DAYNIGHT","BOTH",
		   "Use only DAY, NIGHT or BOTH", &status);

    numbl = 0;
    FITS_write_key(parmfits, TLONG, "PHALOW",&numbl,
		   "Minimum acceptable PHA value", &status);

    numbl = 31;
    FITS_write_key(parmfits, TLONG, "PHAHIGH",&numbl,
		   "Maximum acceptable PHA value", &status);

    numb = 1.0;
    FITS_write_key(parmfits, TFLOAT, "PHA_BKGD",&numb,
		   "Background scaling factor for PHA screening", &status);

    numb=15.0;
    FITS_write_key(parmfits,TFLOAT,"BRITLIMB",&numb,
		   "[deg] Bright limb avoidance angle", &status);

    numb=10.0;
    FITS_write_key(parmfits,TFLOAT,"DARKLIMB",&numb,
		   "[deg] Dark limb avoidance angle", &status);

    numbl=0;
    FITS_write_key(parmfits,TLONG,"NUSERGTI",&numbl,
		   "Number of user defined good time intervals",&status);

    numb=0.0;
    FITS_write_key(parmfits,TFLOAT,"GTIBEG01",&numb,
		   "[sec] Beginning good time interval", &status);

    numb=0.0;
    FITS_write_key(parmfits,TFLOAT,"GTIEND01",&numb,
		   "[sec] Ending good time interval", &status);

    numb=0.0;
    FITS_write_key(parmfits,TFLOAT,"GTIBEG02",&numb,
		   "[sec] Beginning good time interval", &status);

    numb=0.0;
    FITS_write_key(parmfits,TFLOAT,"GTIEND02",&numb,
		   "[sec] Ending good time interval", &status);

    FITS_close_file(parmfits, &status);

}

    return 0;
}
