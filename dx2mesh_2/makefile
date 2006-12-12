CC = gcc
DEFS =
FLAGS = -O2 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o dx2mesh_2.o 

dx2mesh_2: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o dx2mesh_2

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y dx2mesh_2.h
	bison -o parse.c parse.y

parse.o: parse.c dx2mesh_2.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

dx2mesh_2.o: dx2mesh_2.c dx2mesh_2.h
	$(CC) $(FLAGS) $(DEFS) -c dx2mesh_2.c

clean:
	rm -f $(OBJS) lex.c parse.c dx2mesh_2

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp dx2mesh_2 $(INSTALL_DIR); \
	  cp dx2mesh_2 $(INSTALL_DIR); \
	fi

tar:
	tar -cvf dx2mesh_2.tar lex.l parse.y dx2mesh_2.c dx2mesh_2.h makefile

gzip:
	tar -cvf - lex.l parse.y dx2mesh_2.c dx2mesh_2.h makefile | gzip -9c dx2mesh_2.tar.gz
