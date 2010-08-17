COMPILE = g++ -Wall -Wshadow -W -O3 -g -std=c++0x -o meshrefine
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshrefine -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshrefine-node56 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
LINK = g++

meshrefine: classes.cc functions.cc meshrefine.cc  
	$(COMPILE) meshrefine.cc

clean:
	-rm *.o
