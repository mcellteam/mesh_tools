#COMPILE = g++ -O3 -g -o meshrefine
COMPILE = g++ -O3 -g -o meshrefine -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
LINK = g++

meshrefine: classes.cc functions.cc meshrefine.cc  
	$(COMPILE) meshrefine.cc

clean:
	-rm *.o
