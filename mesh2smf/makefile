CC = gcc
DEFS =
FLAGS = -O2 -m486
LIBS = -lfl -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o mesh2smf.o 

mesh2smf: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2smf

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2smf.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2smf.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

mesh2smf.o: mesh2smf.c mesh2smf.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2smf.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2smf

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2smf $(INSTALL_DIR); \
	  cp mesh2smf $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2smf.tar lex.l parse.y mesh2smf.c mesh2smf.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2smf.c mesh2smf.h makefile | gzip -9c mesh2smf.tar.gz
