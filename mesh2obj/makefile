CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.bison.o mesh2obj.o 

mesh2obj: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2obj

lex.flex.c: lex.l
	flex -Cfra lex.l

parse.bison.c: lex.flex.c parse.y mesh2obj.h
	bison parse.y

parse.bison.o: parse.bison.c mesh2obj.h
	$(CC) $(FLAGS) $(DEFS) -c parse.bison.c

mesh2obj.o: mesh2obj.c mesh2obj.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2obj.c

clean:
	rm -f $(OBJS) lex.flex.c parse.bison.c mesh2obj

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh2obj $(INSTALL_DIR); \
	  cp mesh2obj $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh2obj.tar lex.l parse.y mesh2obj.c mesh2obj.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh2obj.c mesh2obj.h makefile | gzip -9c mesh2obj.tar.gz
