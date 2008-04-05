#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshmorph -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshmorph-node56 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
#COMPILE = g++ -m64 -Wall -Wshadow -W -O3 -g -o meshmorph64-node56

meshmorph: meshmorph.cc subroutines.cc classes.cc controls.cc
	$(COMPILE) meshmorph.cc

clean:
	-rm *.o
