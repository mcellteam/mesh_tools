# Author: Justin Kinney
# Date: Sep 2008

#COMPILE = g++ -Wall -Wshadow -W -O0 -g -pg -o meshmorph
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -pg -o meshmorph
#COMPILE = g++ -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math -o meshmorph
#COMPILE = g++ -Wall -Wshadow -Wextra -O3 -g -ffast-math -o meshmorph
#COMPILE = g++ -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math -o meshmorph -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
#COMPILE = g++ -Wall -Wshadow -W -O0 -g -ffast-math -o meshmorph -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -ffast-math -o meshmorph-node56 -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
#COMPILE = g++ -m64 -Wall -Wshadow -W -O3 -g -o meshmorph64-node56

CXX = g++
FLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math -std=c++0x
#FLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math -fno-strict-aliasing
#FLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math -fprofile-arcs -ftest-coverage
#LIBS = -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
#LIBS = -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
LIBS =

SRCS = controls.cc log.cc edge.cc face.cc nice.cc state.cc energy.cc grid.cc      \
       Octree.cc OctreeImplementation.cc OctreeAuxiliary.cc Array.cc Vector3r.cc \
       octree_agent_face.cc octree_visitor_face.cc Wm4Quaternion.cc \
       octree_visitor_check_face.cc octree_visitor_update.cc octree_visitor_measure.cc\
       object.cc vertex.cc container.cc meshmorph.cc opttritri.cc     \
       refractory.cc virtual_disp.cc gain_schedule.cc vertex_schedule.cc     \
       intersecting_faces.cc bestfit.cc
OBJS = $(SRCS:%.cc=%.o)
DEPSINCLUDE = $(addprefix .deps/,$(SRCS:%.cc=%.d))

.deps/%.d: %.cc
	mkdir -p .deps; $(SHELL) -vec '$(CXX) -MM $(FLAGS) $< | sed '\''s@\($*\)\.o[ :]*@\1.o $@ : @g'\'' > $@; [ -s $@ ] || rm -f $@'

meshmorph: $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o meshmorph

%.o: %.cc
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm -f *.o
	rm -rf .deps
	rm -f meshmorph

include $(DEPSINCLUDE)
