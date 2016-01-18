bin_PROGRAMS = variant_test

variant_test_CPPFLAGS = \
     -I$(top_srcdir)/SnowTools/src \
     -I$(top_srcdir)/SnowTools/htslib \
     -I$(top_srcdir)/SnowTools \
     -I$(top_srcdir)/SnowTools/multifast-v1.4.2 \
     -I$(top_srcdir)/src 

variant_test_LDADD = \
	$(top_builddir)/SnowTools/src/libsnowtools.a \
	$(top_builddir)/SnowTools/bwa/libbwa.a \
	$(top_builddir)/SnowTools/htslib/libhts.a \
	@boost_lib@/libboost_regex.a @boost_lib@/libboost_unit_test_framework.a \
	@boost_lib@/libboost_system.a

#variant_test_LDFLAGS = --coverage ##-BOOST_TEST_DYN_LINK

variant_test_SOURCES = variant_test.cpp variant_test_main.cpp 
