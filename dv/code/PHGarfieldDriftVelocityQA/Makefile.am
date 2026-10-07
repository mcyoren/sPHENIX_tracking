AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkgincludedir = $(includedir)/phgarfielddriftvelocityqa

pkginclude_HEADERS = \
  PHGarfieldDriftVelocityQA.h

lib_LTLIBRARIES = \
  libPHGarfieldDriftVelocityQA.la

libPHGarfieldDriftVelocityQA_la_SOURCES = \
  PHGarfieldDriftVelocityQA.cc

libPHGarfieldDriftVelocityQA_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -ltrack_io \
  -lqautils \
  -lcdbobjects \
  -lffamodules \
  -lTpcConditions \
  -lg4detectors \
  -lPHGarfield

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libPHGarfieldDriftVelocityQA.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
