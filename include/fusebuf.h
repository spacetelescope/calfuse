/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Synopsis:	#include	"fusebuf.h"
 *
 * Description: Structure definitions and routine prototyping for multi-line
 *		file buffering as handled by fusebuf.c.
 *
 * History:	08/17/98	gak	Begin work
 *		08/18/98	gak	Tested and working
 *					(as called by cf_make_ff)
 *
 ******************************************************************************/

typedef	struct {
	fitsfile *fits;		/* Pointer to open FITS structure.            */
	int	hdu;		/* HDU containing image to buffer.            */
	int	nx, ny;		/* X and Y dimensions of the image.           */
	float	**buf,		/* Array of pointers to buffered lines.       */
		**y;		/* Pointers to buffered lines in y order.     */
	int	nl,		/* Number of lines buffered.                  */
		yfirst, ylast,	/* Row numbers of buffered lines              */
		znext;		/* Ordinal of next line ptr to use.           */
} imgbuf;

int	cf_openextn(fitsfile *, int, imgbuf *, int, int *);

int	cf_closeextn(imgbuf *, int *);

float	getpixf(imgbuf *, int, int);

