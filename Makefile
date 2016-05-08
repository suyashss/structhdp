CC = g++
#CFLAGS = -g -Wall -O3 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
# CFLAGS = -g -Wall
LDFLAGS = -lgsl -lgslcblas -lm

#GSL_INCLUDE = /usr/include/
GSL_LIB = /usr/lib/

LSOURCE =  utils.cpp structhdp.cpp
LHEADER =  utils.h structhdp.h

structhdp: $(LSOURCE) $(HEADER)
	#$(CC) $(LSOURCE) -o $@ $(LDFLAGS)
	#$(CC)  -I$(GSL_INCLUDE) -L$(GSL_LIB) $(LSOURCE) -o $@ $(LDFLAGS)
	#$(CC) -L$(GSL_LIB) $(LSOURCE) -o $@ $(LDFLAGS)
	$(CC) $(CFLAGS) -c utils.cpp
	$(CC) $(CFLAGS) -c structhdp.cpp
	$(CC) -g -o $@ -L$(GSL_LIB) structhdp.o utils.o  $(LDFLAGS) 

clean:
	rm -f *.o structhdp

test:
	python tests/cram-0.4/cram.py -q tests/tiny.t
	python tests/cram-0.4/cram.py -q tests/stirling.t
