bin_SCRIPTS = \
	add_tle.pl \
	check_tle.pl \
	get_tle.pl \
	get_tle.pl.space-track \
	make_cvzramtool.pl \
	make_orbit.pl \
	test_get_tle.pl \
	update_tle.csh

noinst_SCRIPTS = \
	GetTLE.class \
	GetTLE.java

EXTRA_DIST = \
	GetTLE.class.in \
	GetTLE.java.in \
	add_tle.pl.in \
	check_tle.pl.in \
	get_tle.pl.in \
	get_tle.pl.space-track.in \
	make_cvzramtool.pl.in \
	make_orbit.pl.in \
	test_get_tle.pl.in \
	update_tle.csh.in

CLEANFILES = \
	GetTLE.class \
	GetTLE.java \
	add_tle.pl \
	check_tle.pl \
	get_tle.pl \
	get_tle.pl.space-track \
	make_cvzramtool.pl \
	make_orbit.pl \
	test_get_tle.pl \
	update_tle.csh

do_subst = sed -e 's,[@]DATADIR[@],$(pkgdatadir)-$(VERSION),g' \
		   		   -e 's,[@]PREFIX[@],$(prefix),g' \
				   -e 's,[@]VERSION[@],$(VERSION),g'

GetTLE.class: GetTLE.class.in Makefile
	$(do_subst) < GetTLE.class.in > GetTLE.class
		chmod 644 GetTLE.class

GetTLE.java: GetTLE.java.in Makefile
	$(do_subst) < GetTLE.java.in > GetTLE.java
		chmod 644 GetTLE.java

add_tle.pl: add_tle.pl.in Makefile
	$(do_subst) < add_tle.pl.in > add_tle.pl
		chmod +x add_tle.pl

check_tle.pl: check_tle.pl.in Makefile
	$(do_subst) < check_tle.pl.in > check_tle.pl
		chmod +x check_tle.pl

get_tle.pl: get_tle.pl.in Makefile
	$(do_subst) < get_tle.pl.in > get_tle.pl
		chmod +x get_tle.pl

get_tle.pl.space-track: get_tle.pl.space-track.in Makefile
	$(do_subst) < get_tle.pl.space-track.in > get_tle.pl.space-track
		chmod +x get_tle.pl.space-track

make_cvzramtool.pl: make_cvzramtool.pl.in Makefile
	$(do_subst) < make_cvzramtool.pl.in > make_cvzramtool.pl
		chmod +x make_cvzramtool.pl

make_orbit.pl: make_orbit.pl.in Makefile
	$(do_subst) < make_orbit.pl.in > make_orbit.pl
		chmod +x make_orbit.pl

test_get_tle.pl: test_get_tle.pl.in Makefile
	$(do_subst) < test_get_tle.pl.in > test_get_tle.pl
		chmod +x test_get_tle.pl

update_tle.csh: update_tle.csh.in Makefile
	$(do_subst) < update_tle.csh.in > update_tle.csh
		chmod +x update_tle.csh

install-data-local: GetTLE.class GetTLE.java
	mkdir -p $(DESTDIR)$(pkgdatadir)-$(VERSION)/java
	install -D -m 644 GetTLE.class $(DESTDIR)$(pkgdatadir)-$(VERSION)/java
	install -D -m 644 GetTLE.java $(DESTDIR)$(pkgdatadir)-$(VERSION)/java
