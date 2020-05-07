# Author: Justin Kinney
# Date: Feb 2009

CXX = g++
CC = gcc
#CFLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math
CFLAGS = -Wall -Wshadow -Wextra -Weffc++ -O3 -g -ffast-math -DNOASSERT=true -fpermissive
#LIBS = -L/home/jkinney/install/googleperf64/lib -lprofiler -ltcmalloc
LIBS =
INSTALL_DIR = $(HOME)/bin

SRCS = container.cc container.h contour.cc contour.h controls.cc controls.h histogram.cc histogram.h object.cc object.h point.cc point.h reconstruct_transform.c reconstruct_transform.h reconstruct_interpolate.cc reconstruct_interpolate.h control_points.cc control_points.h parameter.cc parameter.h sim_anneal.cc sim_anneal.h
OBJS = container.o contour.o controls.o histogram.o object.o point.o reconstruct_transform.o reconstruct_interpolate.o control_points.o parameter.o sim_anneal.o

reconstruct_interpolate: $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) -o reconstruct_interpolate

.cc.o: 
	$(CXX) $(CFLAGS) -c $< -o $@

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

install:
	@if [ "$(INSTALL_DIR)" ]; then \
		if [ ! -s $(INSTALL_DIR) ]; then \
			echo mkdir $(INSTALL_DIR); \
			mkdir $(INSTALL_DIR); \
		fi; \
		echo cp reconstruct_interpolate $(INSTALL_DIR); \
		cp reconstruct_interpolate $(INSTALL_DIR); \
	fi


clean:
	rm -f $(OBJS)
	rm -f reconstruct_interpolate

