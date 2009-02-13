COMPILE = g++ -O3 -o meshsimplify
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshsimplify -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshsimplify-node56 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
LINK = g++

meshsimplify: classes.cc functions.cc meshsimplify.cc  
	$(COMPILE) meshsimplify.cc

clean:
	-rm *.o
