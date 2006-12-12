CC = gcc
DEFS =
FLAGS = -O2 -m486
LIBS = -lfl -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o meshfilter.o 

meshfilter: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o meshfilter

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y meshfilter.h
	bison -o parse.c parse.y

parse.o: parse.c meshfilter.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

meshfilter.o: meshfilter.c meshfilter.h
	$(CC) $(FLAGS) $(DEFS) -c meshfilter.c

clean:
	rm -f $(OBJS) lex.c parse.c meshfilter

install: meshfilter
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp meshfilter $(INSTALL_DIR); \
	  cp meshfilter $(INSTALL_DIR); \
	fi

tar:
	tar -cvf meshfilter.tar lex.l parse.y meshfilter.c meshfilter.h makefile

gzip:
	tar -cvf - lex.l parse.y meshfilter.c meshfilter.h makefile | gzip -9c meshfilter.tar.gz
