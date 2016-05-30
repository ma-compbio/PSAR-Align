AM_CPPFLAGS = -I$(top_srcdir)/src 

# AM_CXXFLAGS = -finline-limit=10000 --param inline-unit-growth=70
# Commented out b/c was causing bugs on some older compilers.
#  --RKB & CD 4/14/09

LDADD = \
	$(top_builddir)/src/annealing/libannealing.a \
	$(top_builddir)/src/manager/libmanager.a \
	$(top_builddir)/src/fsa/libfsa.a \
	$(top_builddir)/src/annealing/libannealing.a \
	$(top_builddir)/src/manager/libmanager.a \
	$(top_builddir)/src/seq/libseq.a \
	$(top_builddir)/src/util/libutil.a \
	$(top_builddir)/src/math/libmath.a

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

bin_PROGRAMS = \
	fsa \
	gapcleaner \
	isect_mercator_alignment_gff \
	map_coords \
	map_gff_coords \
	percentid \
	prot2codon \
	slice_fasta \
	slice_fasta_gff \
	slice_mercator_alignment \
	translate

fsa_SOURCES = main.cc
gapcleaner_SOURCES = gapcleaner.cc
isect_mercator_alignment_gff_SOURCES = isect_mercator_alignment_gff.cc
map_coords_SOURCES = map_coords.cc
map_gff_coords_SOURCES = map_gff_coords.cc
percentid_SOURCES = percentid.cc
prot2codon_SOURCES = prot2codon.cc
slice_fasta_SOURCES = slice_fasta.cc
slice_fasta_gff_SOURCES = slice_fasta_gff.cc
slice_mercator_alignment_SOURCES = slice_mercator_alignment.cc
translate_SOURCES = translate.cc

noinst_HEADERS = 
