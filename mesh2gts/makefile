CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o hash.o mesh2gts.o 

mesh2gts: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2gts

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2gts.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2gts.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

hash.o: hash.c hash.h
	$(CC) $(FLAGS) $(DEFS) -c hash.c

mesh2gts.o: mesh2gts.c mesh2gts.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2gts.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2gts

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2gts $(INSTALL_DIR); \
	  cp mesh2gts $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2gts.tar lex.l parse.y mesh2gts.c mesh2gts.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2gts.c mesh2gts.h makefile | gzip -9c mesh2gts.tar.gz
