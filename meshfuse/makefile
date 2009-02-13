COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshfuse -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshfuse-node56 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshfuse
LINK = g++

meshfuse: meshfuse.cc functions.cc classes.cc
	$(COMPILE) meshfuse.cc

clean:
	-rm *.o
