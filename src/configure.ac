# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(variant, 0.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([variant.cpp])
AC_CONFIG_HEADER([config.h])

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


# Check for the bamtools library pathz
# Bamtools has two different path formats
# If it is built in place the files will be in {path}/lib/ and {path}/include/
# If it is installed (with make install) it will be in {path}/include/bamtools
# and {path}/lib/bamtools. Both variants are handled here
AC_ARG_WITH(bamtools, AS_HELP_STRING([--with-bamtools=PATH],
	[specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]))

if test "$with_bamtools" -a -d "$with_bamtools"; then
    bamtools_ldflags="-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib -Wl,-rpath,$with_bamtools/lib/bamtools -L$with_bamtools/lib/bamtools"
    bamtools_include="-I$with_bamtools/include -I$with_bamtools/include/bamtools"
fi

# Only fail on warnings when the --enable-development flag is passed into configure
AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-g -std=c++11 -Wall -Wextra $fail_on_warning -Wno-unknown-pragmas")
AC_SUBST(CXXFLAGS, "-O3")
AC_SUBST(CFLAGS, "-O3")
AC_SUBST(CPPFLAGS, "$CPPFLAGS $bamtools_include")
AC_SUBST(LDFLAGS, "$bamtools_ldflags $LDFLAGS -L.")

# We always need to specify to link in bamtools and ahocorasick
AC_SUBST(LIBS, "$LIBS -lbamtools -lahocorasick")

# Make sure the bamtools headers can be found
AC_CHECK_HEADERS([api/BamReader.h],,[AC_MSG_ERROR([The bamtools library must be installed (http://github.com/pezmaster31/bamtools). You can specify its path with the --with-bamtools=PATH option])])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
