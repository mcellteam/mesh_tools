CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o mesh2irit.o 

mesh2irit: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2irit

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2irit.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2irit.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

mesh2irit.o: mesh2irit.c mesh2irit.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2irit.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2irit

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2irit $(INSTALL_DIR); \
	  cp mesh2irit $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2irit.tar lex.l parse.y mesh2irit.c mesh2irit.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2irit.c mesh2irit.h makefile | gzip -9c mesh2irit.tar.gz
