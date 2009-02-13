COMPILE = g++ -Wall -Wshadow -W -O3 -g -o measure_ecw_error -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o measure_ecw_error-node56 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc

measure_ecw_error: measure_ecw_error.cc subroutines.cc classes.cc controls.cc
	$(COMPILE) measure_ecw_error.cc

clean:
	-rm *.o
