# Author: Justin Kinney
# Date: Feb 2009

#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshalyzer-lite
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshalyzer -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshalyzer-node56 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
#COMPILE = g++ -m64 -Wall -Wshadow -W -O3 -g -o meshalyzer64 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc

CXX = g++
FLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math
#LIBS = -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
LIBS =

SRCS = contour.cc extra_vertex.cc face.cc group.cc meshstitch.cc site.cc vertex.cc
OBJS = $(SRCS:%.cc=%.o)
DEPSINCLUDE = $(addprefix .deps/,$(SRCS:%.cc=%.d))

.deps/%.d: %.cc
	mkdir -p .deps; $(SHELL) -vec '$(CXX) -MM $(FLAGS) $< | sed '\''s@\($*\)\.o[ :]*@\1.o $@ : @g'\'' > $@; [ -s $@ ] || rm -f $@'

meshalyzer: $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o meshalyzer

%.o: %.cc
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm *.o

include $(DEPSINCLUDE)
