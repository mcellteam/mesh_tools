# Author: Justin Kinney
# Date: Feb 2009

CXX = g++
FLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math -std=c++0x
#LIBS = -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
LIBS =

SRCS = boundary.cc box.cc container.cc controls.cc edge.cc face.cc contour.cc \
       face_pair.cc mesh2contours.cc object.cc space.cc stats.cc vertex.cc
OBJS = $(SRCS:%.cc=%.o)
DEPSINCLUDE = $(addprefix .deps/,$(SRCS:%.cc=%.d))

.deps/%.d: %.cc
	mkdir -p .deps; $(SHELL) -vec '$(CXX) -MM $(FLAGS) $< | sed '\''s@\($*\)\.o[ :]*@\1.o $@ : @g'\'' > $@; [ -s $@ ] || rm -f $@'

mesh2contours: $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o mesh2contours

%.o: %.cc
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm *.o

include $(DEPSINCLUDE)
