#COMPILE = g++ -W
COMPILE = g++ -O3
#COMPILE = g++-3.4 -m64 -O3
LINK = g++

mesh_manip: mesh_manip.cc subroutines.cc definitions.cc 
	$(COMPILE) mesh_manip.cc

clean:
	-rm *.o
