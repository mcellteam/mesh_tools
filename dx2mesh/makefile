CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o dx2mesh.o 

dx2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o dx2mesh

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y dx2mesh.h
	bison -o parse.c parse.y

parse.o: parse.c dx2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

dx2mesh.o: dx2mesh.c dx2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c dx2mesh.c

clean:
	rm -f $(OBJS) lex.c parse.c dx2mesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp dx2mesh $(INSTALL_DIR); \
	  cp dx2mesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf dx2mesh.tar lex.l parse.y dx2mesh.c dx2mesh.h makefile

gzip:
	tar -cvf - lex.l parse.y dx2mesh.c dx2mesh.h makefile | gzip -9c dx2mesh.tar.gz
