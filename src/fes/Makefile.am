SUBDIRS = .
bin_PROGRAMS = cf_calfes
cf_calfes_SOURCES = \
		cf_calfes.c cf_fes_init.c cf_fes_read.c cf_fes_write.c \
		cf_fes_cal.c cf_fes_get_cal_image.c cf_limbang.c \
		cf_fes_apply_bias.c cf_fes_apply_flat.c cf_fes_apply_mask.c
cf_calfes_LDADD = ../slalib/libsla.la ../libcf/libcf.la
cf_calfes_CFLAGS = -I$(top_srcdir)/include -I./
