CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o mesh2vtk.o 

mesh2vtk: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2vtk

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2vtk.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2vtk.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

mesh2vtk.o: mesh2vtk.c mesh2vtk.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2vtk.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2vtk

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2vtk $(INSTALL_DIR); \
	  cp mesh2vtk $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2vtk.tar lex.l parse.y mesh2vtk.c mesh2vtk.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2vtk.c mesh2vtk.h makefile | gzip -9c mesh2vtk.tar.gz
