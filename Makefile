LINK = g++
CXX = g++
#CFLAGS = -g -Wall -O3 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF

GSL_LIB = /usr/lib/

SRCDIR := src
OBJDIR := obj
VPATH := $(SRCDIR)

CXXFLAGS := ${CXXFLAGS} ${DEBUGFLAGS} ${OPTFLAGS} ${INCLUDEPATHS}

TARGET := structhdp
CXXSRCS := utils.cpp structhdp.cpp

CXXOBJS := ${addprefix ${OBJDIR}/, ${CXXSRCS:.cpp=.o}}

LIBS := -lgsl -lgslcblas -lm ${LIBS}

optimized:
	$(MAKE) $(TARGET) OPTFLAGS="-O3 -DNDEBUG -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF" DEBUGFLAGS=""
debug:
	$(MAKE) $(TARGET) OPTFLAGS="-O0" DEBUGFLAGS="-g"
fastdebug:
	$(MAKE) $(TARGET) OPTFLAGS="-O3" DEBUGFLAGS="-g -DBOOST_UBLAS_NDEBUG"

${TARGET}: ${OBJDIR} ${CXXOBJS}
	${LINK} ${LINKFLAGS} -o ${TARGET} ${CXXOBJS} ${LIBPATHS} ${LIBS}

${OBJDIR}:
	mkdir -p ${OBJDIR}

${CXXOBJS}: $(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	 ${CXX} ${CXXFLAGS} -o $@ -c $<

clean: 
	-rm -f structhdp
	-rm -f ${OBJDIR}/*.o  ${OBJDIR}/*.mod  
doc:
	cd doc; make

test:
	python tests/cram-0.4/cram.py -q tests/tiny.t
	python tests/cram-0.4/cram.py -q tests/stirling.t


.PHONY: clean optimized debug fastdebug doc test

