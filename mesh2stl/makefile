CC = gcc
DEFS =
FLAGS = -O3
#FLAGS = -Wall -ggdb3 -O0 -fno-strict-aliasing -pg
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o mesh2stl.o 

mesh2stl: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2stl

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2stl.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2stl.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.c vector.h 
	$(CC) $(FLAGS) $(DEFS) -c vector.c

mesh2stl.o: mesh2stl.c mesh2stl.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2stl.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2stl

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2stl $(INSTALL_DIR); \
	  cp mesh2stl $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2stl.tar lex.l parse.y mesh2stl.c mesh2stl.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2stl.c mesh2stl.h makefile | gzip -9c mesh2stl.tar.gz
