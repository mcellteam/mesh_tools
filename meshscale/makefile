CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o meshscale.o 

meshscale: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o meshscale

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y meshscale.h
	bison -o parse.c parse.y

parse.o: parse.c meshscale.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

meshscale.o: meshscale.c meshscale.h
	$(CC) $(FLAGS) $(DEFS) -c meshscale.c

clean:
	rm -f $(OBJS) lex.c parse.c meshscale

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp meshscale $(INSTALL_DIR); \
	  cp meshscale $(INSTALL_DIR); \
	fi

tar:
	tar -cvf meshscale.tar lex.l parse.y meshscale.c meshscale.h makefile

gzip:
	tar -cvf - lex.l parse.y meshscale.c meshscale.h makefile | gzip -9c meshscale.tar.gz
