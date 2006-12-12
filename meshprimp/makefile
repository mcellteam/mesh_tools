COMPILE = g++ -O3 -o meshprimp
#COMPILE = g++-3.4 -m64 -O3 -o meshprimp
LINK = g++

mesprimp: meshprimp.cc subroutines.cc definitions.cc 
	$(COMPILE) meshprimp.cc

clean:
	-rm *.o
