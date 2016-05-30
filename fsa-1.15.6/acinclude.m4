
AC_DEFUN([AC_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ac_cv_cxx_namespaces,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([namespace Outer { namespace Inner { int i = 0; }}],
                [using namespace Outer::Inner; return i;],
 ac_cv_cxx_namespaces=yes, ac_cv_cxx_namespaces=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])

AC_DEFUN([AC_CXX_GNUCXX_HASHMAP],[
AC_CACHE_CHECK(whether the compiler supports __gnu_cxx::hash_map,
ac_cv_cxx_gnucxx_hashmap,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <ext/hash_map>
using __gnu_cxx::hash_map;],
 [],
 ac_cv_cxx_gnucxx_hashmap=yes, ac_cv_cxx_gnucxx_hashmap=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_gnucxx_hashmap" = yes; then
  AC_DEFINE(HAVE_GNUCXX_HASHMAP,,[define if the compiler supports __gnu_cxx::hash_map])
fi
])

AC_DEFUN([AC_CXX_HAVE_EXT_HASH_MAP],
[AC_CACHE_CHECK(whether the compiler has ext/hash_map,
ac_cv_cxx_have_ext_hash_map,
[AC_REQUIRE([AC_CXX_NAMESPACES])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([#include <ext/hash_map>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],[hash_map<int, int> t; return 0;],
  ac_cv_cxx_have_ext_hash_map=yes, ac_cv_cxx_have_ext_hash_map=no)
  AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_ext_hash_map" = yes; then
   AC_DEFINE(HAVE_EXT_HASH_MAP,,[define if the compiler has ext/hash_map])
fi
])

AC_DEFUN([AC_CXX_STLPORT_HASHMAP],[
AC_CACHE_CHECK(whether the compiler supports std::hash_map,
ac_cv_cxx_stlport_hashmap,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <hash_map>
using std::hash_map;],
 [],
 ac_cv_cxx_stlport_hashmap=yes, ac_cv_cxx_stlport_hashmap=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_stlport_hashmap" = yes; then
  AC_DEFINE(HAVE_STLPORT_HASHMAP,,[define if the compiler supports std::hash_map])
fi
])


#################################################################################
# The below definitions for AC_TYPE_UINTMAX_T AND AC_TYPE_UNSIGNED_LONG_LONG_INT
# were taken from the Autoconf types.m4 file;
# they're included here since older distributions of autoconf don't
# have these macros defined.
#################################################################################

# AC_TYPE_UINTMAX_T
# -----------------
AC_DEFUN([AC_TYPE_UINTMAX_T],
[
  AC_REQUIRE([AC_TYPE_UNSIGNED_LONG_LONG_INT])
  AC_CHECK_TYPE([uintmax_t],
    [AC_DEFINE([HAVE_UINTMAX_T], 1,
       [Define to 1 if the system has the type `uintmax_t'.])],
    [test $ac_cv_type_unsigned_long_long_int = yes \
       && ac_type='unsigned long long int' \
       || ac_type='unsigned long int'
     AC_DEFINE_UNQUOTED([uintmax_t], [$ac_type],
       [Define to the widest unsigned integer type
	if <stdint.h> and <inttypes.h> do not define.])])
])


# AC_TYPE_UNSIGNED_LONG_LONG_INT
# ------------------------------
AC_DEFUN([AC_TYPE_UNSIGNED_LONG_LONG_INT],
[
  AC_CACHE_CHECK([for unsigned long long int],
    [ac_cv_type_unsigned_long_long_int],
    [AC_LINK_IFELSE(
       [AC_LANG_PROGRAM(
	  [[unsigned long long int ull = 18446744073709551615ULL;
	    typedef int a[(18446744073709551615ULL <= (unsigned long long int) -1
			   ? 1 : -1)];
	   int i = 63;]],
	  [[unsigned long long int ullmax = 18446744073709551615ull;
	    return (ull << 63 | ull >> 63 | ull << i | ull >> i
		    | ullmax / ull | ullmax % ull);]])],
       [ac_cv_type_unsigned_long_long_int=yes],
       [ac_cv_type_unsigned_long_long_int=no])])
  if test $ac_cv_type_unsigned_long_long_int = yes; then
    AC_DEFINE([HAVE_UNSIGNED_LONG_LONG_INT], 1,
      [Define to 1 if the system has the type `unsigned long long int'.])
  fi
])

