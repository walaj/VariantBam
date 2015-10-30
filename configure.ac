# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(variant, 0.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/variant.cpp])
AC_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE([disable])

# Checks for programs.
AC_PROG_CXX ## test for cpp compiler
AC_PROG_CC ## test for C compiler
AC_PROG_RANLIB ## required if libraries are built in package

# Check for headers
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])

# Check for libraries
##AC_SEARCH_LIBS([library],[function], [action-if-found], [action if not found])
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [clock_getttime found])], )

AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
	[specify directory containing the boost library)]))
if test "$with_boost" -a -d "$with_boost"; then
    boost_include="-I$with_boost"
else
    boostr=/xchip/gistic/Jeremiah/boost_1_58_0
    boost_include="-I$boostr"
fi


########
### CHECK FOR SNOW TOOLS
########
#AC_ARG_WITH(snowtools, AS_HELP_STRING([--with-snowtools=PATH],
#	[specify directory containing the SnowTools library]))
#if test "$with_snowtools" -a -d "$with_snowtools"; then
#   a=1
#else
#    with_snowtools="/xchip/gistic/Jeremiah/GIT/SnowTools"
#fi

#with_snowtools="./SnowTools"
#snowtools_ldflags="-L$with_snowtools/src -L$with_snowtools/src/bwa -L$with_snowtools/src/multifast-v1.4.2/ahocorasick/build"
#snowtools_include="-I$with_snowtools/src -I$with_snowtools/src/multifast-v1.4.2 -I$with_snowtools/src/htslib -I$with_snowtools/src/bwa"
snowtools_libs="-lsnowtools $with_snowtools/src/htslib/libhts.a -lbwa -lahocorasick"

# Only fail on warnings when the --enable-development flag is passed into configure
AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-g -std=c++11 -Wall -Wextra $fail_on_warning -Wno-unknown-pragmas")
##AC_SUBST(AM_CXXFLAGS, "-g -Wall -Wextra $fail_on_warning -Wno-unknown-pragmas")
AC_SUBST(CXXFLAGS, "-O3")
AC_SUBST(CFLAGS, "-O3")
AC_SUBST(CPPFLAGS, "$CPPFLAGS")
AC_SUBST(LDFLAGS, "$LDFLAGS -L.")

AC_SUBST(LIBS, "$LIBS -lpthread")

# Make sure the bamtools headers can be found
##AC_CHECK_HEADERS([htslib/hts.h],,[AC_MSG_ERROR([The HTS library must be installed (https://github.com/samtools/htslib). You can specify its path with the --with-hts=PATH option])])

# Make sure the bamtools headers can be found
##AC_CHECK_HEADERS([boost/icl/interval_set.hpp],,[AC_MSG_WARN([The Boost library must be installed for a few of the interval operations. NOT REQUIRED for VariantBam. Specify its path with the --with-boostxs=PATH option])])

## Make sure the snowtools headers can be found
#AC_CHECK_HEADERS([SnowTools/SnowTools.h],,[AC_MSG_ERROR([The SnowTools library must be installed. You can specify its path with the --with-snowtools=PATH option])])

# Make sure the bamtools headers can be found
#AC_CHECK_HEADERS([ahocorasick/ahocorasick.h],,[AC_MSG_WARN([The Aho-Corasick library must be installed (http://sourceforge.net/projects/multifast/), to do string matching rules. Can specify --with-aho=PATH option])])

AC_CONFIG_FILES([Makefile
                 SnowTools/src/Makefile
	         src/Makefile])

AC_OUTPUT
