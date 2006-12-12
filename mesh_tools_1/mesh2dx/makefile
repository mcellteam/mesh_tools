CC = gcc
DEFS =
FLAGS = -O2 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o mesh2dx.o 

mesh2dx: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2dx

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2dx.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2dx.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

mesh2dx.o: mesh2dx.c mesh2dx.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2dx.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2dx

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2dx $(INSTALL_DIR); \
	  cp mesh2dx $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2dx.tar lex.l parse.y mesh2dx.c mesh2dx.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2dx.c mesh2dx.h makefile | gzip -9c mesh2dx.tar.gz
