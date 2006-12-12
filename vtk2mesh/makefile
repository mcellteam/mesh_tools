CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vtk2mesh.o 

vtk2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o vtk2mesh

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y vtk2mesh.h
	bison -o parse.c parse.y

parse.o: parse.c vtk2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vtk2mesh.o: vtk2mesh.c vtk2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c vtk2mesh.c

clean:
	rm -f $(OBJS) lex.c parse.c vtk2mesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp vtk2mesh $(INSTALL_DIR); \
	  cp vtk2mesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf vtk2mesh.tar lex.l parse.y vtk2mesh.c vtk2mesh.h makefile

gzip:
	tar -cvf - lex.l parse.y vtk2mesh.c vtk2mesh.h makefile | gzip -9c vtk2mesh.tar.gz
