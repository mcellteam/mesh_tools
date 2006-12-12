CC = gcc
DEFS =
FLAGS = -O2 -m486
LIBS = -lfl -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o mesh2off.o 

mesh2off: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2off

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2off.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2off.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

mesh2off.o: mesh2off.c mesh2off.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2off.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2off

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2off $(INSTALL_DIR); \
	  cp mesh2off $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2off.tar lex.l parse.y mesh2off.c mesh2off.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2off.c mesh2off.h makefile | gzip -9c mesh2off.tar.gz
