CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.bison.o mesh2mcell.o 

mesh2mcell: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2mcell

lex.flex.c: lex.l
	flex -Cfra lex.l

parse.bison.c: lex.flex.c parse.y mesh2mcell.h
	bison parse.y

parse.bison.o: parse.bison.c mesh2mcell.h
	$(CC) $(FLAGS) $(DEFS) -c parse.bison.c

mesh2mcell.o: mesh2mcell.c mesh2mcell.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2mcell.c

clean:
	rm -f $(OBJS) lex.flex.c parse.bison.c mesh2mcell

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
