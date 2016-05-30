AM_CPPFLAGS = -I$(top_srcdir)/src -Wno-deprecated

# AM_CXXFLAGS = -finline-limit=10000 --param inline-unit-growth=70
# Commented out b/c was causing bugs on some older compilers.
#  --RKB & CD 4/14/09

if HAVE_CONDOR
AM_CPPFLAGS += -I$(top_srcdir)/MW/src -I$(top_srcdir)/MW/src/MWControlTasks -I$(top_srcdir)/MW/src/RMComm -I$(top_srcdir)/MW/src/RMComm/MW-Socket
endif

noinst_LIBRARIES = libfsa.a

libfsa_a_SOURCES = \
	algebras.cc \
	aminoacid_indel2dp.cc \
	aminoaciddp.cc \
	anchors.cc \
	constraints.cc \
	fsa.cc \
	model.cc \
	nucleotide_indel2dp.cc \
	nucleotidedp.cc \
	sequence_pair_selector.cc

noinst_HEADERS = \
	algebras.h \
	aminoacid_indel2dp.h \
	aminoaciddp.h \
	anchors.h \
	constraints.h \
	dptables.h \
	fsa.h \
	model.h \
	mybanding.h \
	nucleotide_indel2dp.h \
	nucleotidedp.h \
	sequence_pair_selector.h
