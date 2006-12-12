CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vtk2rib.o 

vtk2rib: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o vtk2rib

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y vtk2rib.h
	bison -o parse.c parse.y

parse.o: parse.c vtk2rib.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vtk2rib.o: vtk2rib.c vtk2rib.h
	$(CC) $(FLAGS) $(DEFS) -c vtk2rib.c

clean:
	rm -f $(OBJS) lex.c parse.c vtk2rib

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp vtk2rib $(INSTALL_DIR); \
	  cp vtk2rib $(INSTALL_DIR); \
	fi

tar:
	tar -cvf vtk2rib.tar lex.l parse.y vtk2rib.c vtk2rib.h makefile

gzip:
	tar -cvf - lex.l parse.y vtk2rib.c vtk2rib.h makefile | gzip -9c vtk2rib.tar.gz
