AM_CPPFLAGS = -I$(top_srcdir)/src

LDADD = \
	$(top_builddir)/src/util/libutil.a

noinst_LIBRARIES = libseq.a
libseq_a_SOURCES = \
	alphabet.cc \
	alignment.cc \
	gff.cc \
	mercator.cc \
	sequence.cc \
	similarity_matrix.cc

noinst_HEADERS = \
	alphabet.h \
	alignment.h \
	gff.h \
	interval.h \
	mercator.h \
	sequence.h \
	similarity_matrix.h
