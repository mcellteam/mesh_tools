# Author: Justin Kinney
# Date: Sep 2008

CXX = g++
FLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math -std=c++0x
#LIBS = -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
LIBS =

SRCS = vertex.cc face.cc object.cc edge.cc controls.cc meshheal.cc
OBJS = $(SRCS:%.cc=%.o)
DEPSINCLUDE = $(addprefix .deps/,$(SRCS:%.cc=%.d))

.deps/%.d: %.cc
	mkdir -p .deps; $(SHELL) -vec '$(CXX) -MM $(FLAGS) $< | sed '\''s@\($*\)\.o[ :]*@\1.o $@ : @g'\'' > $@; [ -s $@ ] || rm -f $@'

meshheal: $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o meshheal

%.o: %.cc
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm *.o

include $(DEPSINCLUDE)
