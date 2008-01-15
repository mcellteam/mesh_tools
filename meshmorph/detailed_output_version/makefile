COMPILE = g++ -Wall -Wshadow -W -O3 -g -o meshmorph
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -o porn_viewer
#COMPILE = g++ -Wall -Wshadow -W -O3 -g -L/home/jed/install/lib -ltcmalloc -o meshmorph
#COMPILE = g++-3.4 -m64 -Wall -Wshadow -W -O3 -march=opteron -o meshprimp
#COMPILE = g++-3.4 -m64 -O3 -o meshprimp
#LINK = -I/home/jkinney/packages/boost_1_33_1

mesprimp: meshmorph.cc subroutines.cc classes.cc 
#porn_viewer: meshmorph.cc subroutines.cc classes.cc 
#	$(COMPILE) meshprimp.cc $(LINK)
#	$(COMPILE) meshmorph.cc
	$(COMPILE) meshmorph.cc

clean:
	-rm *.o
