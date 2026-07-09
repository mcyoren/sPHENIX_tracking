AUTOMAKE_OPTIONS = foreign

TPC_PATTERN_RECO_INSTALL = /sphenix/user/mitrankova/F4A/TPC_pattern_reco/install
TPC_PATTERN_RECO_SOURCE = /sphenix/user/mitrankova/F4A/TPC_pattern_reco

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(TPC_PATTERN_RECO_SOURCE) \
  -I$(TPC_PATTERN_RECO_INSTALL)/include \
  -I$(TPC_PATTERN_RECO_INSTALL)/include/inmoduletracks \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(TPC_PATTERN_RECO_INSTALL)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  InModuleTrackDisplay.h \
  FullTrackDisplay.h \
  TpcPolyTrackDisplay.h \
  TpcPolyClusterDisplay.h \
  TpcPolyClusterResiduals.h

lib_LTLIBRARIES = \
  libInModuleTrackDisplay.la

libInModuleTrackDisplay_la_SOURCES = \
  InModuleTrackDisplay.cc \
  FullTrackDisplay.cc \
  TpcPolyTrackDisplay.cc \
  TpcPolyClusterDisplay.cc \
  TpcPolyClusterResiduals.cc

libInModuleTrackDisplay_la_LIBADD = \
  -lInModuleTracks \
  -lffamodules \
  -lcdbobjects \
  -lphool \
  -lfun4all \
  -lg4detectors \
  -lg4detectors_io \
  -ltpc_io \
  -ltrack_io \
  -lSubsysReco

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc

testexternals_LDADD = \
  libInModuleTrackDisplay.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f testexternals.cc
