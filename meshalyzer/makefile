COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshalyzer-lite
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshalyzer -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshalyzer-node56 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
#COMPILE = g++ -m64 -Wall -Wshadow -W -O3 -g -o meshalyzer64 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc

meshalyzer: meshalyzer.cc subroutines.cc classes.cc 
	$(COMPILE) meshalyzer.cc

clean:
	-rm *.o
