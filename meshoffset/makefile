CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
#FLAGS = -Wall -ggdb3 -O0 -fno-strict-aliasing -pg
LIBS = -lm
# LIBS = -lefence -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o meshoffset.o 

meshoffset: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o meshoffset

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y meshoffset.h
	bison -o parse.c parse.y

parse.o: parse.c meshoffset.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

meshoffset.o: meshoffset.c meshoffset.h
	$(CC) $(FLAGS) $(DEFS) -c meshoffset.c

clean:
	rm -f $(OBJS) lex.c parse.c meshoffset

install: meshoffset
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp meshoffset $(INSTALL_DIR); \
	  cp meshoffset $(INSTALL_DIR); \
	fi

tar:
	tar -cvf meshoffset.tar lex.l parse.y meshoffset.c meshoffset.h makefile

gzip:
	tar -cvf - lex.l parse.y meshoffset.c meshoffset.h makefile | gzip -9c meshoffset.tar.gz
