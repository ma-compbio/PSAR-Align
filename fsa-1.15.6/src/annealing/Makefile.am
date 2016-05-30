AM_CPPFLAGS = -I$(top_srcdir)/src

LDADD = \
	libannealing.a \
	$(top_builddir)/src/manager/libmanager.a \
	$(top_builddir)/src/fsa/libfsa.a \
	libannealing.a \
	$(top_builddir)/src/manager/libmanager.a \
	$(top_builddir)/src/seq/libseq.a \
	$(top_builddir)/src/util/libutil.a

if HAVE_CONDOR
AM_CPPFLAGS += -I$(top_srcdir)/MW/src -I$(top_srcdir)/MW/src/MWControlTasks -I$(top_srcdir)/MW/src/RMComm -I$(top_srcdir)/MW/src/RMComm/MW-Socket

AM_LDFLAGS = -static -L$(top_builddir)/MW/lib
LDADD += \
	-lMW \
	-lMWRMComm \
	-lMWutil \
	-lNWS \
	-lMWsocketmaster \
    -lpthread
endif

noinst_LIBRARIES = libannealing.a
libannealing_a_SOURCES = \
	alignment_DAG.cc \
	dotplot.cc \
	tree_weights.cc

noinst_HEADERS = \
	SparseMatrix.h \
	alignment_DAG.h \
	dotplot.h \
	tree_weights.h
