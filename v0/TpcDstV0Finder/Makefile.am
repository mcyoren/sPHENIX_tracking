AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  $(ROOTCFLAGS)

lib_LTLIBRARIES = libTpcDstV0Finder.la

libTpcDstV0Finder_la_SOURCES = \
  TpcDstV0Finder.cc \
  TpcDstV0Finder_Dict.cc

libTpcDstV0Finder_la_LDFLAGS = -module -avoid-version
libTpcDstV0Finder_la_LIBADD = \
  -L$(OFFLINE_MAIN)/lib -lfun4all -lphool -lffaobjects -linmoduletracks \
  $(ROOTLIBS)

pkginclude_HEADERS = TpcDstV0Finder.h

BUILT_SOURCES = TpcDstV0Finder_Dict.cc
CLEANFILES = TpcDstV0Finder_Dict.cc TpcDstV0Finder_Dict_rdict.pcm

TpcDstV0Finder_Dict.cc: TpcDstV0Finder.h LinkDef.h
	rootcling -f $@ -c $(AM_CPPFLAGS) TpcDstV0Finder.h LinkDef.h
