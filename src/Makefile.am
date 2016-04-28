bin_PROGRAMS = variant

variant_CPPFLAGS = \
     -I$(top_srcdir)/SnowTools/src \
     -I$(top_srcdir)/SnowTools/htslib \
     -I$(top_srcdir)/SnowTools \
     -I$(top_srcdir)/SnowTools/multifast-v1.4.2 

variant_LDADD = \
	$(top_builddir)/SnowTools/src/libsnowtools.a \
	$(top_builddir)/SnowTools/htslib/libhts.a \
	$(top_builddir)/SnowTools/multifast-v1.4.2/ahocorasick/build/libahocorasick.a

##variant_LDFLAGS = -Wl,-Bstatic -lboost_regex
#variant_LDFLAGS = @boost_lib@/libboost_regex.a

variant_SOURCES = variant.cpp VariantBamWalker.cpp
