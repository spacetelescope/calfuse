SUBDIRS = .
bin_SCRIPTS = \
			calfuse.csh \
			calfuse.sh \
			calfuse.wavecal \
			cf_make_900_obs.csh \
			cf_make_all.csh \
			cf_make_all_obs.csh \
			cfsetup.csh \
			cfsetup.sh \
			idl_obsplot.pl \
			idlplot_rate.pl \
			idlplot_spex.pl

EXTRA_DIST = \
			cfsetup.csh.in \
			cfsetup.sh.in \
			calfuse.csh \
			calfuse.sh \
			calfuse.wavecal \
			cf_make_900_obs.csh \
			cf_make_all.csh.in \
			cf_make_all_obs.csh \
			idl_obsplot.pl \
			idlplot_rate.pl \
			idlplot_spex.pl

CLEANFILES = cfsetup.csh \
			 cfsetup.sh \
			 cf_make_all.csh

do_subst = sed -e 's,[@]DATADIR[@],$(pkgdatadir)-$(VERSION),g' \
		   -e 's,[@]PREFIX[@],$(prefix),g' \
		   -e 's,[@]PACKAGE[@],$(PACKAGE),g' \
		   -e 's,[@]VERSION[@],$(VERSION),g'

cfsetup.csh: cfsetup.csh.in Makefile
	$(do_subst) < cfsetup.csh.in > cfsetup.csh
	chmod +x cfsetup.csh

cfsetup.sh: cfsetup.sh.in Makefile
	$(do_subst) < cfsetup.sh.in > cfsetup.sh
	chmod +x cfsetup.sh

cf_make_all.csh: cf_make_all.csh.in Makefile
	$(do_subst) < cf_make_all.csh.in > cf_make_all.csh
	chmod +x cf_make_all.csh
