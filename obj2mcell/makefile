CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.bison.o obj2mcell.o 

obj2mcell: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o obj2mcell

lex.flex.c: lex.l
	flex -Cfra lex.l

parse.bison.c: lex.flex.c parse.y obj2mcell.h
	bison parse.y

parse.bison.o: parse.bison.c obj2mcell.h
	$(CC) $(FLAGS) $(DEFS) -c parse.bison.c

obj2mcell.o: obj2mcell.c obj2mcell.h
	$(CC) $(FLAGS) $(DEFS) -c obj2mcell.c

clean:
	rm -f $(OBJS) lex.flex.c parse.bison.c obj2mcell

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp obj2mcell $(INSTALL_DIR); \
	  cp obj2mcell $(INSTALL_DIR); \
	fi

tar:
	tar -cvf obj2mcell.tar lex.l parse.y obj2mcell.c obj2mcell.h makefile

gzip:
	tar -cvf - lex.l parse.y obj2mcell.c obj2mcell.h makefile | gzip -9c obj2mcell.tar.gz
