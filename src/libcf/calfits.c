/*****************************************************************************
 *            Johns Hopkins University
 *            Center For Astrophysical Sciences
 *            FUSE
 *****************************************************************************
 *
 * Description: Wrappers for cfitsio routines which includes error checking.
 *
 * History:     04/14/99        peb   Finished work on subset of routines
 *                                    from cfitsio v2.031.  
 *              04/15/99        emm   Replaced #include FITSIO.h with 
 *                                    #include calfuse.h
 *              04/15/99    barrett   Added FITS_get_rowsize.
 *              04/16/99    barrett   Added FITS_get_hdu_num.
 *              04/20/99    barrett   Added FITS_insert_key_[str, log, lng,
 *                                      fixflt, flt, fixdbl, dbl, fixcmp,
 *                                      cmp, fixdblcmp, dblcmp]
 *              04/23/99    emurphy   Deleted error checking from
 *                                    cf_get_hdu_num
 *		04/04/01    kruk      Force casesen to FALSE in FITS_get_colnum
 *              12/18/03    bjg       Change calfusettag.h to calfuse.h
 *
 ****************************************************************************/

#include "fitsio.h"
#include "calfuse.h"

int FITS_open_file(fitsfile **fptr, const char *filename, int iomode,
                   int *status)
{
    ffopen(fptr, filename, iomode, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_reopen_file(fitsfile *openfptr, fitsfile **newfptr, int *status)
{
    ffreopen(openfptr, newfptr, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_create_file(fitsfile **fptr, const char *filename, int *status)
{
    ffinit(fptr, filename, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_flush_file(fitsfile *fptr, int *status)
{
    ffflus(fptr, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_close_file(fitsfile *fptr, int *status)
{
    ffclos(fptr, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_delete_file(fitsfile *fptr, int *status)
{
    ffdelt(fptr, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_file_name(fitsfile *fptr, char *filename, int *status)
{
    ffflnm(fptr, filename, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_file_mode(fitsfile *fptr, int *filemode, int *status)
{
    ffflmd(fptr, filemode, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_record(fitsfile *fptr, const char *card, int *status)
{
    ffprec(fptr, card, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_key(fitsfile *fptr, int datatype, char *keyname, void *value,
                   char *comm, int *status)
{
    ffpky(fptr, datatype, keyname, value,
          comm, status);
    if (*status) fprintf(stderr, "Error writing keyword %16.16s\n",keyname);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_comment(fitsfile *fptr, const char *comm, int *status)
{
    ffpcom(fptr, comm, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_history(fitsfile *fptr, const char *history, int *status)
{
    ffphis(fptr, history, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_date(fitsfile *fptr, int *status)
{
    ffpdat(fptr, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_copy_key(fitsfile *infptr,fitsfile *outfptr,int incol,int outcol,
                  char *rootname, int *status)
{
    ffcpky(infptr,outfptr,incol,outcol,
           rootname, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_imghdr(fitsfile *fptr, int bitpix, int naxis, long naxes[],
                      int *status)
{
    ffphps( fptr, bitpix, naxis, naxes, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_grphdr(fitsfile *fptr, int simple, int bitpix, int naxis,
                      long naxes[], long pcount, long gcount, int extend,
                      int *status)
{
    ffphpr( fptr, simple, bitpix, naxis, naxes,
            pcount, gcount, extend, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_btblhdr(fitsfile *fptr, long naxis2, int tfields, char **ttype,
                       char **tform, char **tunit, char *extname, long pcount,
                       int *status)
{
    ffphbn(fptr, naxis2, tfields, ttype,
          tform, tunit, extname, pcount, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_hdrpos(fitsfile *fptr, int *nexist, int *position, int *status)
{
    ffghps(fptr, nexist, position, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_read_record(fitsfile *fptr, int nrec, char *card, int *status)
{
    ffgrec(fptr, nrec, card, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_read_card(fitsfile *fptr, char *keyname, char *card, int *status)
{
    ffgcrd(fptr, keyname, card, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_read_key(fitsfile *fptr, int datatype, char *keyname, void *value,
                  char *comm, int *status)
{
    ffgky( fptr, datatype, keyname, value,
           comm, status);
    if (*status) fprintf(stderr, "Error reading keyword %16.16s\n",keyname);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_read_imghdr(fitsfile *fptr, int maxdim, int *simple, int *bitpix,
                     int *naxis, long naxes[], long *pcount, long *gcount,
                     int *extend, int *status)
{
    ffghpr(fptr, maxdim, simple, bitpix, naxis,
          naxes, pcount, gcount, extend, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_read_btblhdr(fitsfile *fptr, int maxfield, long *naxis2, int *tfields,
                      char **ttype, char **tform, char **tunit, char *extname,
                      long *pcount, int *status)
{
    ffghbn(fptr, maxfield, naxis2, tfields,
           ttype, tform, tunit, extname,
           pcount, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_update_card(fitsfile *fptr, char *keyname, char *card, int *status)
{
    ffucrd(fptr, keyname, card, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_update_key(fitsfile *fptr, int datatype, char *keyname, void *value,
                    char *comm, int *status)
{
    ffuky(fptr, datatype, keyname, value,
          comm, status);
    if (*status) fprintf(stderr, "Error updating keyword %16.16s\n",keyname);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_modify_record(fitsfile *fptr, int nkey, char *card, int *status)
{
    ffmrec(fptr, nkey, card, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_modify_card(fitsfile *fptr, char *keyname, char *card, int *status)
{
    ffmcrd(fptr, keyname, card, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_modify_comment(fitsfile *fptr, char *keyname, char *comm, int *status)
{
    ffmcom(fptr, keyname, comm, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_record(fitsfile *fptr, int nkey, char *card, int *status)
{
    ffirec(fptr, nkey, card, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_str(fitsfile *fptr, char *keyname, char *value,
			char *comment, int *status)
{
    ffikys(fptr, keyname, value, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_log(fitsfile *fptr, char *keyname, int value,
			char *comment, int *status)
{
    ffikyl(fptr, keyname, value, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_lng(fitsfile *fptr, char *keyname, long value,
			char *comment, int *status)
{
    ffikyj(fptr, keyname, value, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_flt(fitsfile *fptr, char *keyname, float value,
			int decimals, char *comment, int *status)
{
    ffikye(fptr, keyname, value, decimals, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_fixflt(fitsfile *fptr, char *keyname, float value,
			   int decimals, char *comment, int *status)
{
    ffikyf(fptr, keyname, value, decimals, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_dbl(fitsfile *fptr, char *keyname, double value,
			int decimals, char *comment, int *status)
{
    ffikyd(fptr, keyname, value, decimals, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_fixdbl(fitsfile *fptr, char *keyname, double value,
			   int decimals, char *comment, int *status)
{
    ffikyg(fptr, keyname, value, decimals, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_cmp(fitsfile *fptr, char *keyname, float *value,
			int decimals, char *comment, int *status)
{
    ffikyc(fptr, keyname, value, decimals, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_dblcmp(fitsfile *fptr, char *keyname, double *value,
			   int decimals, char *comment, int *status)
{
    ffikym(fptr, keyname, value, decimals, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_fixcmp(fitsfile *fptr, char *keyname, float *value,
			   int decimals, char *comment, int *status)
{
    ffikfc(fptr, keyname, value, decimals, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_fixdblcmp(fitsfile *fptr, char *keyname, double *value,
			      int decimals, char *comment, int *status)
{
    ffikfm(fptr, keyname, value, decimals, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_key_null(fitsfile *fptr, char *keyname, char *comment,
			 int *status)
{
    ffikyu(fptr, keyname, comment, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_delete_key(fitsfile *fptr, char *keyname, int *status)
{
    ffdkey(fptr, keyname, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_delete_record(fitsfile *fptr, int keypos, int *status)
{
    ffdrec(fptr, keypos, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_hdu_num(fitsfile *fptr, int *hdunum)
{
    /* This routine returns the HDU number instead of status.  Therefore,
     * it should not use any error checking. */
    int status;
    status = ffghdn(fptr, hdunum);
    return status;
}

int FITS_get_hdu_type(fitsfile *fptr, int *exttype, int *status)
{
    ffghdt(fptr, exttype, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_movabs_hdu(fitsfile *fptr, int hdunum, int *exttype, int *status)
{
    ffmahd(fptr, hdunum, exttype, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_movrel_hdu(fitsfile *fptr, int hdumov, int *exttype, int *status)
{
    ffmrhd(fptr, hdumov, exttype, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_movnam_hdu(fitsfile *fptr, int exttype, char *hduname, int hduvers,
                    int *status)
{
    ffmnhd(fptr, exttype, hduname, hduvers,
           status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_num_hdus(fitsfile *fptr, int *nhdu, int *status)
{
    ffthdu(fptr, nhdu, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_create_img(fitsfile *fptr, int bitpix, int naxis, long *naxes,
                    int *status)
{
    ffcrim(fptr, bitpix, naxis, naxes, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_create_tbl(fitsfile *fptr, int tbltype, long naxis2, int tfields,
                    char **ttype, char **tform, char **tunit, char *extname,
                    int *status)
{
    ffcrtb(fptr, tbltype, naxis2, tfields, ttype,
           tform, tunit, extname, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_create_hdu(fitsfile *fptr, int *status)
{
    ffcrhd(fptr, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_img(fitsfile *fptr, int bitpix, int naxis, long *naxes,
                    int *status)
{
    ffiimg(fptr, bitpix, naxis, naxes, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_atbl(fitsfile *fptr, long naxis1, long naxis2, int tfields,
                     char **ttype, long *tbcol, char **tform, char **tunit,
                     char *extname, int *status)
{
    ffitab(fptr, naxis1, naxis2, tfields, ttype,
           tbcol, tform, tunit, extname, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_btbl(fitsfile *fptr,long naxis2, int tfields, char **ttype,
                     char **tform, char **tunit, char *extname, long pcount,
                     int *status)
{
    ffibin(fptr,naxis2, tfields, ttype, tform,
           tunit, extname, pcount, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_delete_hdu(fitsfile *fptr, int *hdutype, int *status)
{
    ffdhdu(fptr, hdutype, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_copy_hdu(fitsfile *infptr, fitsfile *outfptr, int morekeys,
                  int *status)
{
    ffcopy(infptr, outfptr, morekeys, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_copy_header(fitsfile *infptr, fitsfile *outfptr, int *status)
{
    ffcphd(infptr, outfptr, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_copy_data(fitsfile *infptr, fitsfile *outfptr, int *status)
{
    ffcpdt(infptr, outfptr, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_colnum(fitsfile *fptr, int casesen, char *templt, int  *colnum,
                    int *status)
{
    /* force casesen to FALSE, in accordance with the FITS standard */
    casesen = FALSE;
    ffgcno(fptr, casesen, templt, colnum,
           status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_colname(fitsfile *fptr, int casesen, char *templt, char *colname,
                     int *colnum, int *status)
{
    ffgcnn(fptr, casesen, templt, colname,
           colnum, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_coltype(fitsfile *fptr, int colnum, int *typecode, long *repeat,
                     long *width, int *status)
{
    ffgtcl(fptr, colnum, typecode, repeat,
           width, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_num_rows(fitsfile *fptr, long *nrows, int *status)
{
    ffgnrw(fptr, nrows, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_rowsize(fitsfile *fptr, long *nrows, int *status)
{
    ffgrsz(fptr, nrows, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_num_cols(fitsfile *fptr, int  *ncols, int *status)
{
    ffgncl(fptr, ncols, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_acolparms(fitsfile *fptr, int colnum, char *ttype, long *tbcol,
                       char *tunit, char *tform, double *tscal, double *tzero,
                       char *tnull, char *tdisp, int *status)
{
    ffgacl(fptr, colnum, ttype, tbcol,
           tunit, tform, tscal, tzero,
           tnull, tdisp, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_get_bcolparms(fitsfile *fptr, int colnum, char *ttype, char *tunit,
                       char *dtype, long *repeat, double *tscal, double *tzero,
                       long *tnull, char *tdisp, int  *status)
{
    ffgbcl(fptr, colnum, ttype, tunit,
           dtype, repeat, tscal, tzero,
           tnull, tdisp, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_read_img(fitsfile *fptr, int  datatype, long firstelem, long nelem,
                  void *nulval, void *array, int *anynul, int  *status)
{
    ffgpv(fptr,  datatype, firstelem, nelem,
          nulval, array, anynul, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_read_col( fitsfile *fptr, int datatype, int colnum, long firstrow,
                   long firstelem, long nelem, void *nulval, void *array,
                   int *anynul, int  *status)
{
    ffgcv( fptr, datatype, colnum, firstrow,
           firstelem, nelem, nulval, array, anynul,
           status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_img(fitsfile *fptr, int datatype, long  firstelem, long  nelem,
                   void  *array, int  *status)
{
    ffppr(fptr, datatype,  firstelem,  nelem,
          array, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_write_col(fitsfile *fptr, int datatype, int colnum, long firstrow,
                   long firstelem, long nelem, void *array, int *status)
{
    ffpcl(fptr, datatype, colnum, firstrow,
          firstelem, nelem, array, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_rows(fitsfile *fptr, long firstrow, long nrows, int *status)
{
    ffirow(fptr, firstrow, nrows, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_delete_rows(fitsfile *fptr, long firstrow, long nrows, int *status)
{
    ffdrow(fptr, firstrow, nrows, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_delete_rowlist(fitsfile *fptr, long *rownum,  long nrows, int *status)
{
    ffdrws(fptr, rownum,  nrows, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_col(fitsfile *fptr, int numcol, char *ttype, char *tform,
                    int *status)
{
    fficol(fptr, numcol, ttype, tform, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_insert_cols(fitsfile *fptr, int firstcol, int ncols, char **ttype,
                     char **tform, int *status)
{
    fficls(fptr, firstcol, ncols, ttype,
           tform, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_delete_col(fitsfile *fptr, int numcol, int *status)
{
    ffdcol(fptr, numcol, status);
    cf_if_fits_error(*status);
    return *status;
}

int FITS_copy_col(fitsfile *infptr, fitsfile *outfptr, int incol, int outcol, 
                  int create_col, int *status)
{
    ffcpcl(infptr, outfptr, incol, outcol, 
           create_col, status);
    cf_if_fits_error(*status);
    return *status;
}

