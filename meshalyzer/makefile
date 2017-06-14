# Author: Justin Kinney
# Date: Feb 2009

CXX = g++
FLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -std=c++11
LIBS = 
INSTALL_DIR = $(HOME)/bin

SRCS = boundary.cc box.cc container.cc controls.cc edge.cc face.cc \
       face_pair.cc meshalyzer.cc object.cc space.cc stats.cc vertex.cc
OBJS = $(SRCS:%.cc=%.o)
DEPSINCLUDE = $(addprefix .deps/,$(SRCS:%.cc=%.d))

.deps/%.d: %.cc
	mkdir -p .deps; $(SHELL) -vec '$(CXX) -MM $(FLAGS) $< | sed '\''s@\($*\)\.o[ :]*@\1.o $@ : @g'\'' > $@; [ -s $@ ] || rm -f $@'

meshalyzer: $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o meshalyzer

%.o: %.cc
	$(CXX) $(FLAGS) -c $< -o $@

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp meshalyzer $(INSTALL_DIR); \
	  cp meshalyzer $(INSTALL_DIR); \
	fi


clean:
	rm -f *.o
	rm -rf .deps
	rm -f meshalyzer

include $(DEPSINCLUDE)
