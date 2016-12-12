SUBDIRS = .
bin_PROGRAMS = \
	cf_ttag_to_hist \
	get_shift \
	mjd2hjd \
	cf_reflux \
	idf_screen \
	extract_jitter \
	modhead \
	cf_arith \
	cf_combine \
	cf_nvo \
	cf_pack \
	cf_xcorr \
	idf_combine \
	bpm_combine \
	ttag_combine \
	cf_coadd \
	cf_make_all_exp \
	remove_target_orbital_motion \
	idf_cut \
	ttag_lightcurve \
	ttag_lightcurve_channel_sum \
	ttag_lightcurve_mjd2hmjd \
	ttag_lightcurve_combine \
	ttag_lightcurve_periodogram

bpm_combine_SOURCES = bpm_combine.c
cf_arith_SOURCES = cf_arith.c cf_wrspec7.c cf_wrspec_cf2.c
cf_coadd_SOURCES = cf_coadd.c
cf_combine_SOURCES = cf_combine.c cf_wrspec7.c
cf_make_all_exp_SOURCES = cf_make_all_exp.c
cf_nvo_SOURCES = cf_nvo.c
cf_pack_SOURCES = cf_pack.c
cf_reflux_SOURCES = cf_reflux.c
cf_ttag_to_hist_SOURCES = cf_ttag_to_hist.c
cf_xcorr_SOURCES = cf_xcorr.c
extract_jitter_SOURCES = extract_jitter.c
get_shift_SOURCES = get_shift.c
idf_combine_SOURCES = idf_combine.c
idf_cut_SOURCES = idf_cut.c
idf_screen_SOURCES = idf_screen.c
mjd2hjd_SOURCES = mjd2hjd.c gethmjd.c
modhead_SOURCES = modhead.c
remove_target_orbital_motion_SOURCES = remove_target_orbital_motion.c gethmjd.c
ttag_combine_SOURCES = ttag_combine.c
ttag_lightcurve_SOURCES = ttag_lightcurve.c
ttag_lightcurve_channel_sum_SOURCES = ttag_lightcurve_channel_sum.c
ttag_lightcurve_combine_SOURCES = ttag_lightcurve_combine.c
ttag_lightcurve_mjd2hmjd_SOURCES = ttag_lightcurve_mjd2hmjd.c gethmjd.c
ttag_lightcurve_periodogram_SOURCES = ttag_lightcurve_periodogram.c

LDADD = ../libcf/libcf.la ../slalib/libsla.la $(cfitsio_LIBS)
AM_CFLAGS = $(cfitsio_CFLAGS) -I$(srcdir)/include -I./ -DCFORTRAN -Df2cFortran
AM_LDFLAGS = -lm
