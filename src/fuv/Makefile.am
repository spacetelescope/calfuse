SUBDIRS = .
BUILD_LIBS = ../slalib/libsla.la ../libcf/libcf.la $(cfitsio_LIBS)
BUILD_FLAGS = -I../../include -I./ $(cfitsio_CFLAGS)
bin_PROGRAMS = cf_assign_wavelength \
			   cf_bad_pixels \
			   cf_convert_to_farf \
			   cf_countmap \
			   cf_extract_spectra \
			   cf_flux_calibrate \
			   cf_gainmap \
			   cf_hist_init \
			   cf_remove_motions \
			   cf_screen_photons \
			   cf_ttag_init
cf_assign_wavelength_SOURCES = cf_assign_wavelength.c
cf_bad_pixels_SOURCES = cf_bad_pixels.c
cf_convert_to_farf_SOURCES = cf_convert_to_farf.c
cf_countmap_SOURCES = cf_countmap.c
cf_extract_spectra_SOURCES = cf_extract_spectra.c
cf_flux_calibrate_SOURCES = cf_flux_calibrate.c
cf_gainmap_SOURCES = cf_gainmap.c
cf_hist_init_SOURCES = cf_hist_init.c
cf_remove_motions_SOURCES = cf_remove_motions.c
cf_screen_photons_SOURCES = cf_screen_photons.c
cf_ttag_init_SOURCES = cf_ttag_init.c

cf_assign_wavelength_LDADD = $(BUILD_LIBS)
cf_bad_pixels_LDADD = $(BUILD_LIBS)
cf_convert_to_farf_LDADD = $(BUILD_LIBS)
cf_countmap_LDADD = $(BUILD_LIBS)
cf_extract_spectra_LDADD = $(BUILD_LIBS)
cf_flux_calibrate_LDADD = $(BUILD_LIBS)
cf_gainmap_LDADD = $(BUILD_LIBS)
cf_hist_init_LDADD = $(BUILD_LIBS)
cf_remove_motions_LDADD = $(BUILD_LIBS)
cf_screen_photons_LDADD = $(BUILD_LIBS)
cf_ttag_init_LDADD = $(BUILD_LIBS)

cf_assign_wavelength_CFLAGS = $(BUILD_FLAGS)
cf_bad_pixels_CFLAGS = $(BUILD_FLAGS)
cf_convert_to_farf_CFLAGS = $(BUILD_FLAGS)
cf_countmap_CFLAGS = $(BUILD_FLAGS)
cf_extract_spectra_CFLAGS = $(BUILD_FLAGS)
cf_flux_calibrate_CFLAGS = $(BUILD_FLAGS)
cf_gainmap_CFLAGS = $(BUILD_FLAGS)
cf_hist_init_CFLAGS = $(BUILD_FLAGS)
cf_remove_motions_CFLAGS = $(BUILD_FLAGS)
cf_screen_photons_CFLAGS = $(BUILD_FLAGS)
cf_ttag_init_CFLAGS = $(BUILD_FLAGS)
