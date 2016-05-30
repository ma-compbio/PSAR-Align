
# Initialize.
m4_include([version.m4])
AC_INIT([FSA], VERSION_NUMBER, [fsa@math.berkeley.edu])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror tar-ustar foreign])

# Set language to c++.
AC_LANG([C++])

 # Make sure CXXFLAGS is defined so that AC_PROG_CXX doesn't set it.
CFLAGS="$CFLAGS"
CXXFLAGS="$CXXFLAGS"

# Check for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL

# Check for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Check for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_UINTMAX_T
AC_CHECK_TYPES([ptrdiff_t])

# Check for unordered_map.
AC_CHECK_HEADERS([tr1/unordered_map tr1/unordered_set])

# Check for library functions.
AC_FUNC_FORK
AC_CHECK_FUNCS([floor memmove pow regcomp sqrt strchr strcspn strspn strstr])

# Check the platform information.
AC_CANONICAL_HOST

# Set CFLAGS and CXXFLAGS.
user_CFLAGS=${CFLAGS}
generic_CFLAGS="-Wall"
ext_CFLAGS=""
debug_CFLAGS=""

# Try to choose appropriate flags for 64-bit compilation.
# (This code is borrowed from Cole Trapnell's TopHat (http://tophat.cbcb.umd.edu/).)
case "${host_cpu}-${host_os}" in
   i*86-*linux*)
    ext_CFLAGS="-march=i686";;
   i*86-darwin*)
    AC_ARG_ENABLE([universal],
		  [AS_HELP_STRING([--enable-universal],
			          [enable universal binary (Intel Apple only)])],
		  [is_uni=1], [is_uni=0])
	case $is_uni in
	  1) CFLAGS="-m64"
	     AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
			       [ext_CFLAGS="-arch x86_64 -arch i386 -arch ppc64 -arch ppc"],
			       [ext_CFLAGS="-arch i386 -arch ppc"]);;
	  0) CFLAGS="-m64"
	     AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
			       [ext_CFLAGS="-arch x86_64"],
			       []);;
	esac;;
  *)
    AC_MSG_CHECKING([if gcc accepts -m64])
    CFLAGS="-m64"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
		      [ext_CFLAGS="-m64"; AC_MSG_RESULT([yes])],
		      [ext_CFLAGS="-D_FILE_OFFSET_BITS=64"; AC_MSG_RESULT([no])]);;
esac

AC_ARG_ENABLE([intel64],
	      [AS_HELP_STRING([--enable-intel64],
		              [optimize for Intel64 CPU such as Xeon and Core2])],
	      [ext_CFLAGS="${ext_CFLAGS} -mtune=nocona"], [])

# Option to manually force 32-bit compilation.
AC_ARG_ENABLE([force-32bit],
	      [AS_HELP_STRING([--enable-force-32bit],
			      [force 32-bit compilation (takes precedence over all other compilation flags)])],
	      [ext_CFLAGS=""], [])

# Add options for compiler flags.
AC_ARG_ENABLE([optim],
              [AS_HELP_STRING([--enable-optim@<:@=0|1|2|3@:>@],
	                      [set optimization level (default is 3)])],
              [if test "x$enable_optim" = xyes; then enable_optim=3; fi],
              [enable_optim=3])
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
	                      [enable debugging info (default is no)])],
              [], [enable_debug=no])
AC_ARG_ENABLE([profile],
              [AS_HELP_STRING([--enable-profile],
	                      [enable profiling (default is no)])],
              [], [enable_profile=no])
AC_ARG_ENABLE([fast-math],
              [AS_HELP_STRING([--disable-fast-math],
	                      [disable fast-math optimization])],
              [], [enable_fast_math=yes])

# Set compiler flags according to specified options.
AS_IF([test "x$enable_optim" != xno], [ext_CFLAGS="$ext_CFLAGS -O$enable_optim"])
AS_IF([test "x$enable_debug" = xyes],
      [debug_CFLAGS="-g -DDEBUG"],
      [debug_CFLAGS="-DNDEBUG"])
AS_IF([test "x$enable_profile" = xyes], [ext_CFLAGS="$ext_CFLAGS -pg"])
AS_IF([test "x$enable_fast_math" = xyes], [ext_CFLAGS="$ext_CFLAGS -ffast-math"])

# Set CFLAGS and CXXFLAGS with all of the specified options.
CFLAGS="${generic_CFLAGS} ${ext_CFLAGS} ${user_CFLAGS} ${debug_CFLAGS}"
CXXFLAGS="$CFLAGS"

# Add option for compiling for condor standard environment.
AC_ARG_ENABLE([condor-compile],
              [AS_HELP_STRING([--enable-condor-compile],
	                      [Build for condor standard environment
                               (default is no)])],
              [], [enable_condor_compile=no])

# Set compiler depending on condor option.
AS_IF([test "x$enable_condor_compile" = xyes],
      [AC_CHECK_PROG([HAVE_CONDOR_COMPILE],[condor_compile],["yes"])
       if test "x$HAVE_CONDOR_COMPILE" = xyes ; then
          CXX="condor_compile $CXX"; CC="condor_compile $CC";
       else
          AC_MSG_FAILURE([--enable-condor given, but condor_compile was not found in path])
       fi])

# Add option for building applet.
AC_ARG_ENABLE([applet],
              [AS_HELP_STRING([--enable-applet],
	                      [Build the MAD GUI as an applet
                               (default is no)])],
              [],
	      [enable_applet=no])

# Define whether to use applet or not.
AS_IF([test "x$enable_applet" = xyes],
      [MAD_MAIN_CLASS=mad.MadApplet],
      [MAD_MAIN_CLASS=mad.MAD])
AC_ARG_VAR([MAD_MAIN_CLASS],[Main class for MAD GUI])

# Add option for using mummer.
AC_ARG_WITH([mummer],
            [AS_HELP_STRING([--with-mummer@<:@=yes|no|check|PATH@:>@],
                            [Use MUMmer for aligning long sequences. Specify
                             PATH to MUMmer executable if not in current path.
                             [check]])],
            [],
	    [with_mummer=check])

# Search for mummer executable if not specified.
AS_IF([test "x$with_mummer" = xcheck || test "x$with_mummer" = xyes],
      [AC_PATH_PROG([MUMMER_EXEC],[mummer],["no"])],
      [MUMMER_EXEC=$with_mummer])

# Check mummer executable and if everything is correct, define macro.
AS_IF([test "x$MUMMER_EXEC" = xno && test "x$with_mummer" = xyes],
      [AC_MSG_FAILURE([--with-mummer given, but mummer was not found in path])])
AS_IF([test "x$MUMMER_EXEC" != xno && test ! -x $MUMMER_EXEC],
      [AC_MSG_FAILURE(["$MUMMER_EXEC" is not a valid mummer executable])])
AS_IF([test "x$MUMMER_EXEC" != xno],
      [AC_DEFINE_UNQUOTED([MUMMER_EXEC],["$MUMMER_EXEC"],[Path to mummer])])

# Add option for using exonerate.
AC_ARG_WITH([exonerate],
            [AS_HELP_STRING([--with-exonerate@<:@=yes|no|check|PATH@:>@],
                            [Use exonerate for aligning long sequences. Specify
                             PATH to exonerate executable if not in current path.
                             [check]])],
            [],
	    [with_exonerate=check])

# Search for exonerate executable if not specified.
AS_IF([test "x$with_exonerate" = xcheck || test "x$with_exonerate" = xyes],
      [AC_PATH_PROG([EXONERATE_EXEC],[exonerate],["no"])],
      [EXONERATE_EXEC=$with_exonerate])

# Check exonerate executable and if everything is correct, define macro.
AS_IF([test "x$EXONERATE_EXEC" = xno && test "x$with_exonerate" = xyes],
      [AC_MSG_FAILURE([--with-exonerate given, but exonerate was not found in path])])
AS_IF([test "x$EXONERATE_EXEC" != xno && test ! -x $EXONERATE_EXEC],
      [AC_MSG_FAILURE(["$EXONERATE_EXEC" is not a valid exonerate executable])])
AS_IF([test "x$EXONERATE_EXEC" != xno],
      [AC_DEFINE_UNQUOTED([EXONERATE_EXEC],["$EXONERATE_EXEC"],[Path to exonerate])])

# Add option for using postgresql.
AC_ARG_WITH([postgresql],
            [AS_HELP_STRING([--with-postgresql@<:@=yes|no|check@:>@],
                            [Enable option to use database mode using
                             PostgreSQL [no]])],
            [],
	    [with_postgresql=no])

# Check for postgresql.
AS_IF([test "x$with_postgresql" = xcheck || test "x$with_postgresql" = xyes],
      [# Check for Postgres library
       AC_CHECK_LIB([pq], [PQconnectdb],
                    # Check for Postgres header
                    [AC_CHECK_HEADER([libpq-fe.h],
                      # Check for a function in the Postgres header
                      [AC_CHECK_DECL([PQgetCopyData],
                                     [HAVE_POSTGRES=yes], [],
                                     [#include <libpq-fe.h>])])])])

AS_IF([test "x$HAVE_POSTGRES" != xyes && test "x$with_postgresql" = xyes],
      [AC_MSG_FAILURE([--with-postgresql given, but valid postgresql installation not found])])

AS_IF([test "x$HAVE_POSTGRES" = xyes],
      [# Check for crypt library, which is required by postgres on Linux
       AC_CHECK_LIB([crypt], [crypt], [LIBS="-lcrypt $LIBS"])
       # Add postgres library
       LIBS="-lpq $LIBS";
       AC_DEFINE_UNQUOTED([HAVE_POSTGRES], [1],
                          [Define to 1 if you have postgres])])
AM_CONDITIONAL([HAVE_POSTGRES], [test "x$HAVE_POSTGRES" = xyes])

# Add option for using condor (parellel mode).
AC_ARG_WITH([condor],
            [AS_HELP_STRING([--with-condor@<:@=yes|no|check@:>@],
                            [Enable option to use parallel mode using
                             condor [no]])],
            [],
	    [with_condor=no])

# Check for condor.
AS_IF([test "x$with_condor" = xcheck || test "x$with_condor" = xyes],
      [AC_CHECK_PROG([HAVE_CONDOR],[condor_submit],["yes"])])

AS_IF([test "x$HAVE_CONDOR" != xyes && test "x$with_condor" = xyes],
      [AC_MSG_FAILURE([--with-condor given, but valid condor installation not found])])

AS_IF([test "x$HAVE_CONDOR" = xyes],
      [AC_DEFINE_UNQUOTED([HAVE_CONDOR],[1],[Define to 1 if you have condor])])

AM_CONDITIONAL([HAVE_CONDOR], [test "x$HAVE_CONDOR" = xyes])

# Build with MW if Condor is found.
AS_IF([test "x$HAVE_CONDOR" = xyes], [AC_CONFIG_SUBDIRS([MW])])

# Check for structures/functions that can be used to determine system memory.
AC_CHECK_MEMBERS([struct sysinfo.totalram, struct sysinfo.mem_unit],
                 [], [], [#include <sys/sysinfo.h>])
AC_CHECK_DECLS([sysctl, CTL_HW, HW_PHYSMEM], [], [], [#include <sys/sysctl.h>])

# Check for javac.  If it exists, set automake conditional.
AC_CHECK_PROG([HAVE_JAVAC], [javac], [yes], [no])
AM_CONDITIONAL([HAVE_JAVAC], [test "x$HAVE_JAVAC" = xyes])

AC_CONFIG_FILES([Makefile
                 src/main/Makefile
                 src/manager/Makefile
                 src/fsa/Makefile
                 src/annealing/Makefile
                 src/math/Makefile
                 src/util/Makefile
                 src/seq/Makefile
                 perl/Makefile
                 examples/Makefile
		 doc/Makefile
		 doc/version.tex
		 html/Makefile
		 tests/Makefile])

AC_OUTPUT


# Dump some configuration confirmations.
echo \
"
-- ${PACKAGE_STRING} Configuration Results --
  C++ compiler:        ${CXX} ${CXXFLAGS}"

if test x"${GCC}" = x"yes" ; then
   gcc_version=`${CC} --version | head -n 1`
   echo "  GCC version:         ${gcc_version}"
else
   gcc_version=''
fi

echo \
"  Host System type:    ${host}
  Install prefix:      ${prefix}

  See config.h for further configuration information.
  Email <${PACKAGE_BUGREPORT}> with questions and bug reports.
"
