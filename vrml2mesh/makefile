CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
#FLAGS = -Wall -ggdb3 -O0 -fno-strict-aliasing -pg
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o vrml2mesh.o 

vrml2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o vrml2mesh

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y vrml2mesh.h
	bison -o parse.c parse.y

parse.o: parse.c vrml2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.c vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

vrml2mesh.o: vrml2mesh.c vrml2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c vrml2mesh.c

clean:
	rm -f $(OBJS) lex.c parse.c vrml2mesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp vrml2mesh $(INSTALL_DIR); \
	  cp vrml2mesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf vrml2mesh.tar lex.l parse.y vrml2mesh.c vrml2mesh.h makefile

gzip:
	tar -cvf - lex.l parse.y vrml2mesh.c vrml2mesh.h makefile | gzip -9c vrml2mesh.tar.gz
