CC = gcc
DEFS =
#FLAGS = -O3 -march=athlon-xp
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o mesh2rib.o 

mesh2rib: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2rib

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2rib.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2rib.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.c vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

mesh2rib.o: mesh2rib.c mesh2rib.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2rib.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2rib

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2rib $(INSTALL_DIR); \
	  cp mesh2rib $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2rib.tar lex.l parse.y vector.c vector.h mesh2rib.c mesh2rib.h makefile

gzip:
	tar -cvf - lex.l parse.y vector.c vector.h mesh2rib.c mesh2rib.h makefile | gzip -9c mesh2rib.tar.gz
