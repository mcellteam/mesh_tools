CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o mesh2mcell.o 

mesh2mcell: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2mcell

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh2mcell.h
	bison -o parse.c parse.y

parse.o: parse.c mesh2mcell.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

mesh2mcell.o: mesh2mcell.c mesh2mcell.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2mcell.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh2mcell

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2mcell $(INSTALL_DIR); \
	  cp mesh2mcell $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2mcell.tar lex.l parse.y mesh2mcell.c mesh2mcell.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2mcell.c mesh2mcell.h makefile | gzip -9c mesh2mcell.tar.gz
