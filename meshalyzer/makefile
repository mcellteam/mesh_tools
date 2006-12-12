CC = gcc
DEFS =
#FLAGS = -O3 -march=pentiumpro
FLAGS = -Wall -ggdb3 -O0 -fno-strict-aliasing -pg
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o meshalyzer.o 

meshalyzer: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o meshalyzer

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y meshalyzer.h
	bison -o parse.c parse.y

parse.o: parse.c meshalyzer.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.c vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

meshalyzer.o: meshalyzer.c meshalyzer.h
	$(CC) $(FLAGS) $(DEFS) -c meshalyzer.c

clean:
	rm -f $(OBJS) lex.c parse.c meshalyzer

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp meshalyzer $(INSTALL_DIR); \
	  cp meshalyzer $(INSTALL_DIR); \
	fi

tar:
	tar -cvf meshalyzer.tar lex.l parse.y meshalyzer.c meshalyzer.h makefile

gzip:
	tar -cvf - lex.l parse.y meshalyzer.c meshalyzer.h makefile | gzip -9c meshalyzer.tar.gz
