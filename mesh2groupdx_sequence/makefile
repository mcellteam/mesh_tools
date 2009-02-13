COMPILE = g++ -Wall -Wshadow -W -O3 -g -o mesh2groupdx_sequence
#COMPILE = g++ -m64 -Wall -Wshadow -W -O3 -g -o mesh2groupdx_sequence64

mesh2groupdx: mesh2groupdx_sequence.cc functions.cc 
	$(COMPILE) mesh2groupdx_sequence.cc

clean:
	-rm *.o
