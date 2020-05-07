CC = gcc
CXX = g++
DEFS =
FLAGS = -O3 -fpermissive
# FLAGS = -g
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = approximate.o branch.o branch_util.o contour.o contour_read.o correspond.o cross.o CTSliceH2DInit.o CTSliceHartley2D.o CTVolumeHartley3D.o decom_util.o decompose.o filter.o gen_data.o group.o image2cont2D.o legal_region.o main.o math_util.o mymarch.o myutil.o parse_config.o qsort.o read_slice.o scale.o slice.o tile_hash.o tile_util.o tiling_main.o volume.o


contour_tiler: $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o contour_tiler

approximate.o: approximate.cpp approximate.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c approximate.cpp

branch.o: branch.cpp branch.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c branch.cpp

branch_util.o: branch_util.cpp branch_util.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c branch_util.cpp

contour.o: contour.cpp contour.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c contour.cpp

contour_read.o: contour_read.cpp contour_read.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c contour_read.cpp

correspond.o: correspond.cpp correspond.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c correspond.cpp

cross.o: cross.cpp cross.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c cross.cpp

CTSliceH2DInit.o: CTSliceH2DInit.cpp common.h
	$(CXX) $(FLAGS) $(DEFS) -c CTSliceH2DInit.cpp

CTSliceHartley2D.o: CTSliceHartley2D.cpp common.h
	$(CXX) $(FLAGS) $(DEFS) -c CTSliceHartley2D.cpp

CTVolumeHartley3D.o: CTVolumeHartley3D.cpp common.h
	$(CXX) $(FLAGS) $(DEFS) -c CTVolumeHartley3D.cpp

decom_util.o: decom_util.cpp decom_util.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c decom_util.cpp

decompose.o: decompose.cpp decompose.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c decompose.cpp

filter.o: filter.cpp filter.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c filter.cpp

gen_data.o: gen_data.cpp gen_data.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c gen_data.cpp

group.o: group.cpp group.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c group.cpp

image2cont2D.o: image2cont2D.cpp image2cont2D.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c image2cont2D.cpp

legal_region.o: legal_region.cpp legal_region.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c legal_region.cpp

main.o: main.cpp main.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c main.cpp

math_util.o: math_util.cpp math_util.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c math_util.cpp

mymarch.o: mymarch.cpp mymarch.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c mymarch.cpp

myutil.o: myutil.cpp myutil.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c myutil.cpp

parse_config.o: parse_config.cpp parse_config.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c parse_config.cpp

qsort.o: qsort.cpp qsort.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c qsort.cpp

read_slice.o: read_slice.cpp read_slice.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c read_slice.cpp

scale.o: scale.cpp scale.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c scale.cpp

slice.o: slice.cpp common.h
	$(CXX) $(FLAGS) $(DEFS) -c slice.cpp

tile_hash.o: tile_hash.cpp tile_hash.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c tile_hash.cpp

tile_util.o: tile_util.cpp tile_util.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c tile_util.cpp

tiling_main.o: tiling_main.cpp tiling_main.h common.h
	$(CXX) $(FLAGS) $(DEFS) -c tiling_main.cpp

volume.o: volume.cpp common.h
	$(CXX) $(FLAGS) $(DEFS) -c volume.cpp


clean:
	rm -f $(OBJS) contour_tiler

install: contour_tiler
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp contour_tiler $(INSTALL_DIR); \
	  cp contour_tiler $(INSTALL_DIR); \
	fi

