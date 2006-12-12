CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o stl2mesh.o 

stl2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o stl2mesh

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y vector.h stl2mesh.h
	bison -o parse.c parse.y

parse.o: parse.c vector.h stl2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.c vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

stl2mesh.o: stl2mesh.c vector.h stl2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c stl2mesh.c

clean:
	rm -f $(OBJS) lex.c parse.c stl2mesh

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
