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
char extname[]="PARAMETER FILE";
char filename[80];


    int   i, status=0, vers;
    long  numbl, dumar[2]={0, 0};
    float numb, effmjd;
    fitsfile *parmfits;

for (i=0; i<4; i++) {  /* 1A,1B,2A,2B */

    sprintf(filename,"parm%2.2s%03d.fit", segments[i],VERSION);
    printf("%20.20s\n",filename); 

    FITS_create_file(&parmfits,filename,&status);
    FITS_create_img(parmfits, SHORT_IMG, 0, dumar, &status);

    FITS_write_key(parmfits, TSTRING, "CALFTYPE", "PARM", 
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

    FITS_write_comment(parmfits, "  ", &status);
    FITS_write_comment(parmfits, "The SPEX_SIC and SPEX_LIF keywords allow "
        "the user to center the extraction windows along a given row rather "
        "than having the pipeline try to centroid the spectrum.  A number "
        "less than 0 or greater than 1023 will default to pipeline "
        "centroiding.  A value between 0 and 1023 will define the center "
        "of the extraction window.  Since the pipeline expands histogram "
        "images out to 1023 pixels in Y, the number must be given for a "
        "full 1023 image height.  That is, if the center is determined "
        "from binned images, the center must be multiplied by the binning "
        "factor before entering it in this file.", &status);
    FITS_write_comment(parmfits, "  ", &status);

    numbl = -1;
    FITS_write_key(parmfits, TLONG, "SPEX_SIC",&numbl,
		   "SiC extraction window Y center (0-1023)", &status);

    numbl = -1;
    FITS_write_key(parmfits, TLONG, "SPEX_LIF",&numbl,
		   "LiF extraction window Y center (0-1023)", &status);

    FITS_write_comment(parmfits, "  ", &status);
    FITS_write_comment(parmfits, "The EMAX_SIC and EMAX_LIF keywords allow "
        "the user to limit the amount that the pipeline can shift the "
        "extraction windows based on the centroid of the spectrum.  In "
        "cases where the source is very faint, the centroid routine "
        "may find bright detector artifacts or other apertures instead "
        "of the desired aperture.  If the calculated centroid differs from "
        "the predicted centroid by more than EMAX_SIC or EMAX_LIF, the "
        "pipeline uses the default extraction window location.", &status);
    FITS_write_comment(parmfits, "  ", &status);

    numbl = 20;
    FITS_write_key(parmfits, TLONG, "EMAX_SIC",&numbl,
		   "SiC extraction window maximum Y movement", &status);

    numbl = 20;
    FITS_write_key(parmfits, TLONG, "EMAX_LIF",&numbl,
		   "LiF extraction window maximum Y movement", &status);

    FITS_close_file(parmfits, &status);

}

    return 0;
}
