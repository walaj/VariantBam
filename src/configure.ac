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
##AC_ARG_WITH(bamtools, AS_HELP_STRING([--with-bamtools=PATH],
##	[specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]))

##bamtools_libs=""
##if test "$with_bamtools" -a -d "$with_bamtools"; then
##   bamtools_ldflags="-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib -Wl,-rpath,$with_bamtools/lib/bamtools -L$with_bamtools/lib/bamtools"
##    bamtools_include="-I$with_bamtools/include -I$with_bamtools/include/bamtools"
##    bamtools_libs="-lbamtools"
##fi

########
### CHECK FOR SNOW TOOLS
########
AC_ARG_WITH(snowtools, AS_HELP_STRING([--with-snowtools=PATH],
	[specify directory containing the SnowTools library]))
if test "$with_snowtools" -a -d "$with_snowtools"; then
    snowtools_ldflags="-L$with_snowtools"
    snowtools_include="-Iwith_snowtools"
    snowtools_libs="-lsnowtools"
else
    snowstable="SnowToolsStable/src"
    snowcurrent="/xchip/gistic/Jeremiah/GIT/SnowTools/src"
    snowtools_ldflags="-L$snowcurrent"
    snowtools_include="-I$snowcurrent"
    snowtools_libs="-lsnowtools"
fi  

########
### CHECK FOR AHO-CORASICK LIBRARY
########
# Check for the aho library paths
AC_ARG_WITH(aho, AS_HELP_STRING([--with-aho=PATH],
        [specify directory containing the AhoCorasick string matching algorithm C code (http://sourceforge.net/projects/multifast/)]))

aholib=""
if test "$with_aho" -a -d "$with_aho"; then
    aho_include="-I$with_aho"
    aho_ldflags="-L$with_aho/ahocorasick/build"
    aholib="-lahocorasick"
##else
##    aho_include="-I ThirdParty/multifast-v1.4.2"
##    aho_ldflags="-L ThirdParty/multifast-v1.4.2/ahocorasick/build"
##    aholib="-lahocorasick"
else
    aho=/xchip/gistic/Jeremiah/GIT/isva/ThirdParty/multifast-v1.4.2
    aho_include="-I$aho"
    aho_ldflags="-L$aho/ahocorasick/build"
    aholib="-lahocorasick"
fi   

#AC_ARG_WITH(hts, AS_HELP_STRING([--with-hts=PATH],
#	[specify directory containing the HTS library (https://github.com/samtools/htslib)]))

#if test "$with_hts" -a -d "$with_hts"; then
#    hts_ldflags="-Wl,-rpath,$with_hts -L$with_hts"
#    hts_include="-I$with_hts"
#else
#    htsloc="/xchip/gistic/Jeremiah/GIT/htslib"
#    hts_ldflags="-Wl,-rpath,$htsloc,-L$htsloc"
#    hts_include="-I$htsloc"
#fi


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
AC_SUBST(CPPFLAGS, "$CPPFLAGS $aho_include $hts_include $snowtools_include")
AC_SUBST(LDFLAGS, "$aho_ldflags $hts_ldflags $snowtools_ldflags $LDFLAGS -L.")

# We always need to specify to link in bamtools and ahocorasick
AC_SUBST(LIBS, "$LIBS $aholib $snowtools_libs -lhts -lpthread")

# Make sure the bamtools headers can be found
#AC_CHECK_HEADERS([htslib/hts.h],,[AC_MSG_ERROR([The HTS library must be installed (https://github.com/samtools/htslib). You can specify its path with the --with-hts=PATH option])])

## Make sure the snowtools headers can be found
AC_CHECK_HEADERS([SnowTools/SnowTools.h],,[AC_MSG_ERROR([The SnowTools library must be installed. You can specify its path with the --with-snowtools=PATH option])])

# Make sure the bamtools headers can be found
AC_CHECK_HEADERS([ahocorasick/ahocorasick.h],,[AC_MSG_WARN([The Aho-Corasick library must be installed (http://sourceforge.net/projects/multifast/), to do string matching rules. Can specify --with-aho=PATH option])])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
