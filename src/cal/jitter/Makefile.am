SUBDIRS = .
bin_PROGRAMS = cf_jitter_diag \
			   cf_jitter
cf_jitter_diag_SOURCES = cf_jitter_diag.c
cf_jitter_SOURCES = cf_jitter.c
LDADD = ../../libcf/libcf.la ../../slalib/libsla.la $(cfitsio_LIBS)
AM_CFLAGS = -I$(top_srcdir)/include -I./ -DCFORTRAN -Df2cFortran $(cfitsio_CFLAGS)
