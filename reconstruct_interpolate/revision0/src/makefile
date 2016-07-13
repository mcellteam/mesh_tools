# Author: Justin Kinney
# Date: Feb 2009

CXX = g++
#FLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math
FLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math -DNOASSERT=true
LIBS = -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
#LIBS =

SRCS = container.cc contour.cc controls.cc histogram.cc object.cc point.cc \
       reconstruct2contourtiler.cc control_points.cc parameter.cc sim_anneal.cc
OBJS = $(SRCS:%.cc=%.o)
DEPSINCLUDE = $(addprefix .deps/,$(SRCS:%.cc=%.d))

.deps/%.d: %.cc
	mkdir -p .deps; $(SHELL) -vec '$(CXX) -MM $(FLAGS) $< | sed '\''s@\($*\)\.o[ :]*@\1.o $@ : @g'\'' > $@; [ -s $@ ] || rm -f $@'

reconstruct2contourtiler: $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o reconstruct2contourtiler

%.o: %.cc
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm *.o

include $(DEPSINCLUDE)
