AM_CPPFLAGS = -I${top_srcdir} -I${top_builddir}

if AUTO_EXAMPLES
bin_PROGRAMS = AdjointGsr
TESTS = AdjointGsr$(EXEEXT)
else
noinst_PROGRAMS = AdjointGsr
endif
AdjointGsr_SOURCES = AdjointGsr.cpp
AdjointGsr_LDADD = ../../ql/libQuantLib.la

EXTRA_DIST =

.PHONY: examples check-examples

examples: AdjointGsr$(EXEEXT)

check-examples: examples
	./AdjointGsr$(EXEEXT) --batch

dist-hook:
	mkdir -p $(distdir)/bin
	mkdir -p $(distdir)/build
