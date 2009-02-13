CC = gcc
DEFS =
FLAGS = -O2 -march=pentiumpro -fomit-frame-pointer
LIBS = -lfl -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o mesh2ribwireframe.o 

mesh2ribwireframe: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2ribwireframe

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2ribwireframe.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2ribwireframe.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.c vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

mesh2ribwireframe.o: mesh2ribwireframe.c mesh2ribwireframe.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2ribwireframe.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2ribwireframe

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2ribwireframe $(INSTALL_DIR); \
	  cp mesh2ribwireframe $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2ribwireframe.tar lex.l parse.y mesh2ribwireframe.c mesh2ribwireframe.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2ribwireframe.c mesh2ribwireframe.h makefile | gzip -9c mesh2ribwireframe.tar.gz
