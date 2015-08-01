bin_PROGRAMS = variant

variant_CPPFLAGS = \
     -I$(top_srcdir)/SnowTools/src \
     -I$(top_srcdir)/SnowTools/src/htslib \
     -I$(top_srcdir)/SnowTools/src/bwa \
     -I$(top_srcdir)/SnowTools/src/multifast-v1.4.2

variant_LDADD = \
	$(top_builddir)/SnowTools/src/libsnowtools.a \
	$(top_builddir)/SnowTools/src/bwa/libbwa.a \
	$(top_builddir)/SnowTools/src/htslib/libhts.a 

variant_SOURCES = variant.cpp VariantBamWalker.cpp
