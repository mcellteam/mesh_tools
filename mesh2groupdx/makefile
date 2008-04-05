#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o mesh2groupdx
COMPILE = g++ -m64 -Wall -Wshadow -W -O3 -g -o mesh2groupdx64

mesh2groupdx: mesh2groupdx.cc functions.cc 
	$(COMPILE) mesh2groupdx.cc

clean:
	-rm *.o
