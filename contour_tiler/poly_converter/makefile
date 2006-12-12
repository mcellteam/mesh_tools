CC = gcc
CXX = g++
DEFS =
FLAGS = -O3 -march=pentiumpro 
# FLAGS = -g
LIBS = -L/usr/local/lib -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = header.o iv2raw.o new_adjust.o poly.o vrml.o vrml2raw.o wrl2raw.o


poly_convert: $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o poly_convert

header.o: header.cpp header.h
	$(CXX) $(FLAGS) $(DEFS) -c header.cpp

iv2raw.o: iv2raw.cpp header.h
	$(CXX) $(FLAGS) $(DEFS) -c iv2raw.cpp

new_adjust.o: new_adjust.cpp header.h
	$(CXX) $(FLAGS) $(DEFS) -c new_adjust.cpp

poly.o: poly.cpp header.h
	$(CXX) $(FLAGS) $(DEFS) -c poly.cpp

vrml.o: vrml.cpp header.h
	$(CXX) $(FLAGS) $(DEFS) -c vrml.cpp

vrml2raw.o: vrml2raw.cpp header.h
	$(CXX) $(FLAGS) $(DEFS) -c vrml2raw.cpp

wrl2raw.o: wrl2raw.cpp readwrl.h header.h
	$(CXX) $(FLAGS) $(DEFS) -c wrl2raw.cpp



clean:
	rm -f $(OBJS) poly_convert

install: poly_convert
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp poly_convert $(INSTALL_DIR); \
	  cp poly_convert $(INSTALL_DIR); \
	fi

