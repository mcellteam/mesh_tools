CC = gcc
DEFS =
#FLAGS = -O3 -march=pentiumpro
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o meshclip.o 

meshclip: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o meshclip

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y meshclip.h
	bison -o parse.c parse.y

parse.o: parse.c meshclip.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

meshclip.o: meshclip.c meshclip.h
	$(CC) $(FLAGS) $(DEFS) -c meshclip.c

clean:
	rm -f $(OBJS) lex.c parse.c meshclip

install: meshclip
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp meshclip $(INSTALL_DIR); \
	  cp meshclip $(INSTALL_DIR); \
	fi

tar:
	tar -cvf meshclip.tar lex.l parse.y meshclip.c meshclip.h makefile

gzip:
	tar -cvf - lex.l parse.y meshclip.c meshclip.h makefile | gzip -9c meshclip.tar.gz
