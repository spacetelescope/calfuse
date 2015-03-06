SUBDIRS = . src scripts doc
pkgdatadir = $(datadir)/@PACKAGE@-@VERSION@
ACLOCAL_AMFLAGS = -I m4
DATADIRS = calfiles parmfiles idl
INCLUDE_FILES = \
			 src/fes/cf_calfes.h \
			 src/include/calfitsio.h \
			 src/include/calfuse.h \
			 src/include/cfortran.h \
			 src/include/compress.h \
			 src/include/config.h \
			 src/include/drvrsmem.h \
			 src/include/eval_defs.h \
			 src/include/eval_tab.h \
			 src/include/f77_wrap.h \
			 src/include/fusebuf.h \
			 src/include/group.h \
			 src/include/grparser.h \
			 src/include/imcompress.h \
			 src/include/longnam.h \
			 src/include/pctype.h \
			 src/include/region.h \
			 src/include/ricecomp.h \
			 src/include/sgp4.h \
			 src/include/slalib.h \
			 src/include/slamac.h

noinst_HEADERS = $(INCLUDE_FILES)

EXTRA_DIST = $(DATADIRS)

install-data-local: $(DATADIRS)
	mkdir -p $(DESTDIR)$(pkgdatadir)
	cp -a $(DATADIRS) $(DESTDIR)$(pkgdatadir)

