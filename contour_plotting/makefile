COMPILE = g++ -Wall -Wshadow -W -O3 -g -o mesh2contours

meshalyzer: mesh2contours.cc subroutines.cc classes.cc 
	$(COMPILE) mesh2contours.cc

clean:
	-rm *.o
