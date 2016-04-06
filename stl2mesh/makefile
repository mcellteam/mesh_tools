CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.bison.o vector.o stl2mesh.o 

stl2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o stl2mesh

lex.flex.c: lex.l
	flex -Cfra lex.l

parse.bison.c: lex.flex.c parse.y vector.h stl2mesh.h
	bison parse.y

parse.bison.o: parse.bison.c vector.h stl2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c parse.bison.c

vector.o: vector.c vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

stl2mesh.o: stl2mesh.c vector.h stl2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c stl2mesh.c

clean:
	rm -f $(OBJS) lex.flex.c parse.bison.c stl2mesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp stl2mesh $(INSTALL_DIR); \
	  cp stl2mesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf stl2mesh.tar lex.l parse.y stl2mesh.c stl2mesh.h makefile

gzip:
	tar -cvf - lex.l parse.y stl2mesh.c stl2mesh.h makefile | gzip -9c stl2mesh.tar.gz
