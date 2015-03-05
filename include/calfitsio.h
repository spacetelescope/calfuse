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
 *              04/14/99    barrett   Delete calfuse.h from file
 *              04/15/99    barrett   Added FITS_get_rowsize
 *              04/16/99    barrett   Added FITS_get_hdu_num
 *              04/20/99    barrett   Added FITS_insert_key_[str, log, lng,
 *                                      fixflt, flt, fixdbl, dbl, fixcmp,
 *                                      cmp, fixdblcmp, dblcmp]
 *
 ****************************************************************************/

#ifndef FITSIO_H
#define FITSIO_H

#include "fitsio.h"

int FITS_open_file(fitsfile **fptr, const char *filename, int iomode,
                   int *status);
int FITS_reopen_file(fitsfile *openfptr, fitsfile **newfptr, int *status);
int FITS_create_file(fitsfile **fptr, const char *filename, int *status);
int FITS_flush_file(fitsfile *fptr, int *status);
int FITS_close_file(fitsfile *fptr, int *status);
int FITS_delete_file(fitsfile *fptr, int *status);
int FITS_file_name(fitsfile *fptr, char *filename, int *status);
int FITS_file_mode(fitsfile *fptr, int *filemode, int *status);
int FITS_write_record(fitsfile *fptr, const char *card, int *status);
int FITS_write_key(fitsfile *fptr, int datatype, char *keyname, void *value,
                   char *comm, int *status);
int FITS_write_comment(fitsfile *fptr, const char *comm, int *status);
int FITS_write_history(fitsfile *fptr, const char *history, int *status);
int FITS_write_date(fitsfile *fptr, int *status);
int FITS_copy_key(fitsfile *infptr,fitsfile *outfptr,int incol,int outcol,
                  char *rootname, int *status);
int FITS_write_imghdr(fitsfile *fptr, int bitpix, int naxis, long naxes[],
                      int *status);
int FITS_write_grphdr(fitsfile *fptr, int simple, int bitpix, int naxis,
                      long naxes[], long pcount, long gcount, int extend,
                      int *status);
int FITS_write_btblhdr(fitsfile *fptr, long naxis2, int tfields, char **ttype,
                       char **tform, char **tunit, char *extname, long pcount,
                       int *status);
int FITS_get_hdrpos(fitsfile *fptr, int *nexist, int *position, int *status);
int FITS_read_record(fitsfile *fptr, int nrec,      char *card, int *status);
int FITS_read_card(fitsfile *fptr, char *keyname, char *card, int *status);
int FITS_read_key(fitsfile *fptr, int datatype, char *keyname, void *value,
                  char *comm, int *status);
int FITS_read_imghdr(fitsfile *fptr, int maxdim, int *simple, int *bitpix,
                     int *naxis, long naxes[], long *pcount, long *gcount,
                     int *extend, int *status);
int FITS_read_btblhdr(fitsfile *fptr, int maxfield, long *naxis2, int *tfields,
                      char **ttype, char **tform, char **tunit, char *extname,
                      long *pcount, int *status);
int FITS_update_card(fitsfile *fptr, char *keyname, char *card, int *status);
int FITS_update_key(fitsfile *fptr, int datatype, char *keyname, void *value,
                    char *comm, int *status);
int FITS_modify_record(fitsfile *fptr, int nkey, char *card, int *status);
int FITS_modify_card(fitsfile *fptr, char *keyname, char *card, int *status);
int FITS_modify_comment(fitsfile *fptr, char *keyname, char *comm,
                        int *status);
int FITS_insert_record(fitsfile *fptr, int nkey, char *card, int *status);
int FITS_insert_key_str(fitsfile *fptr, char *keyname, char *value,
			char *comment, int *status);
int FITS_insert_key_log(fitsfile *fptr, char *keyname, int value,
			char *comment, int *status);
int FITS_insert_key_lng(fitsfile *fptr, char *keyname, long value,
			char *comment, int *status);
int FITS_insert_key_flt(fitsfile *fptr, char *keyname, float value,
			int decimals, char *comment, int *status);
int FITS_insert_key_fixflt(fitsfile *fptr, char *keyname, float value,
			   int decimals, char *comment, int *status);
int FITS_insert_key_dbl(fitsfile *fptr, char *keyname, double value,
			int decimals, char *comment, int *status);
int FITS_insert_key_fixdbl(fitsfile *fptr, char *keyname, double value,
			   int decimals, char *comment, int *status);
int FITS_insert_key_cmp(fitsfile *fptr, char *keyname, float *value,
			int decimals, char *comment, int *status);
int FITS_insert_key_dblcmp(fitsfile *fptr, char *keyname, double *value,
			   int decimals, char *comment, int *status);
int FITS_insert_key_fixcmp(fitsfile *fptr, char *keyname, float *value,
			   int decimals, char *comment, int *status);
int FITS_insert_key_fixdblcmp(fitsfile *fptr, char *keyname, double *value,
			      int decimals, char *comment, int *status);
int FITS_insert_key_null(fitsfile *fptr, char *keyname, char *comment,
			 int *status);
int FITS_delete_key(fitsfile *fptr, char *keyname, int *status);
int FITS_delete_record(fitsfile *fptr, int keypos, int *status);
int FITS_get_hdu_num(fitsfile *fptr, int *hdunum);
int FITS_get_hdu_type(fitsfile *fptr, int *exttype, int *status);
int FITS_movabs_hdu(fitsfile *fptr, int hdunum, int *exttype, int *status);
int FITS_movrel_hdu(fitsfile *fptr, int hdumov, int *exttype, int *status);
int FITS_movnam_hdu(fitsfile *fptr, int exttype, char *hduname, int hduvers,
                    int *status);
int FITS_get_num_hdus(fitsfile *fptr, int *nhdu, int *status);
int FITS_create_img(fitsfile *fptr, int bitpix, int naxis, long *naxes,
                    int *status);
int FITS_create_tbl(fitsfile *fptr, int tbltype, long naxis2, int tfields,
                    char **ttype, char **tform, char **tunit, char *extname,
                    int *status);
int FITS_create_hdu(fitsfile *fptr, int *status);
int FITS_insert_img(fitsfile *fptr, int bitpix, int naxis, long *naxes,
                    int *status);
int FITS_insert_atbl(fitsfile *fptr, long naxis1, long naxis2, int tfields,
                     char **ttype, long *tbcol, char **tform, char **tunit,
                     char *extname, int *status);
int FITS_insert_btbl(fitsfile *fptr,long naxis2, int tfields, char **ttype,
                     char **tform, char **tunit, char *extname, long pcount,
                     int *status);
int FITS_delete_hdu(fitsfile *fptr, int *hdutype, int *status);
int FITS_copy_hdu(fitsfile *infptr, fitsfile *outfptr, int morekeys,
                  int *status);
int FITS_copy_header(fitsfile *infptr, fitsfile *outfptr, int *status);
int FITS_copy_data(fitsfile *infptr, fitsfile *outfptr, int *status);
int FITS_get_colnum(fitsfile *fptr, int casesen, char *templt, int  *colnum,
                    int *status);
int FITS_get_colname(fitsfile *fptr, int casesen, char *templt, char *colname,
                     int *colnum, int *status);
int FITS_get_coltype(fitsfile *fptr, int colnum, int *typecode, long *repeat,
           long *width, int *status);
int FITS_get_num_rows(fitsfile *fptr, long *nrows, int *status);
int FITS_get_rowsize(fitsfile *fptr, long *nrows, int *status);
int FITS_get_num_cols(fitsfile *fptr, int  *ncols, int *status);
int FITS_get_acolparms(fitsfile *fptr, int colnum, char *ttype, long *tbcol,
                       char *tunit, char *tform, double *tscal, double *tzero,
                       char *tnull, char *tdisp, int *status);
int FITS_get_bcolparms(fitsfile *fptr, int colnum, char *ttype, char *tunit,
                       char *dtype, long *repeat, double *tscal, double *tzero,
                       long *tnull, char *tdisp, int  *status);
int FITS_read_img(fitsfile *fptr, int  datatype, long firstelem, long nelem,
                  void *nulval, void *array, int *anynul, int  *status);
int FITS_read_col(fitsfile *fptr, int datatype, int colnum, long firstrow,
                  long firstelem, long nelem, void *nulval, void *array,
                  int *anynul, int  *status);
int FITS_write_img(fitsfile *fptr, int datatype, long  firstelem, long  nelem,
                   void  *array, int  *status);
int FITS_write_col(fitsfile *fptr, int datatype, int colnum, long firstrow,
                   long firstelem, long nelem, void *array, int *status);
int FITS_insert_rows(fitsfile *fptr, long firstrow, long nrows, int *status);
int FITS_delete_rows(fitsfile *fptr, long firstrow, long nrows, int *status);
int FITS_delete_rowlist(fitsfile *fptr, long *rownum,  long nrows,
                        int *status);
int FITS_insert_col(fitsfile *fptr, int numcol, char *ttype, char *tform,
                    int *status);
int FITS_insert_cols(fitsfile *fptr, int firstcol, int ncols, char **ttype,
                     char **tform, int *status);
int FITS_delete_col(fitsfile *fptr, int numcol, int *status);
int FITS_copy_col(fitsfile *infptr, fitsfile *outfptr, int incol, int outcol, 
                  int create_col, int *status);

#endif
