CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o mesh2smesh.o 

mesh2smesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2smesh

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2smesh.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2smesh.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

mesh2smesh.o: mesh2smesh.c mesh2smesh.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2smesh.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2smesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2smesh $(INSTALL_DIR); \
	  cp mesh2smesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2smesh.tar lex.l parse.y mesh2smesh.c mesh2smesh.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2smesh.c mesh2smesh.h makefile | gzip -9c mesh2smesh.tar.gz
