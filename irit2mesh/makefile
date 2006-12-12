CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o irit2mesh.o 

irit2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o irit2mesh

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y irit2mesh.h
	bison -o parse.c parse.y

parse.o: parse.c irit2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.c vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

irit2mesh.o: irit2mesh.c irit2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c irit2mesh.c

clean:
	rm -f $(OBJS) lex.c parse.c irit2mesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp irit2mesh $(INSTALL_DIR); \
	  cp irit2mesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf irit2mesh.tar lex.l parse.y irit2mesh.c irit2mesh.h makefile

gzip:
	tar -cvf - lex.l parse.y irit2mesh.c irit2mesh.h makefile | gzip -9c irit2mesh.tar.gz
