CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.bison.o obj2mesh.o 

obj2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o obj2mesh

lex.flex.c: lex.l
	flex -Cfra lex.l

parse.bison.c: lex.flex.c parse.y obj2mesh.h
	bison parse.y

parse.bison.o: parse.bison.c obj2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c parse.bison.c

obj2mesh.o: obj2mesh.c obj2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c obj2mesh.c

clean:
	rm -f $(OBJS) lex.flex.c parse.bison.c obj2mesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp obj2mesh $(INSTALL_DIR); \
	  cp obj2mesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf obj2mesh.tar lex.l parse.y obj2mesh.c obj2mesh.h makefile

gzip:
	tar -cvf - lex.l parse.y obj2mesh.c obj2mesh.h makefile | gzip -9c obj2mesh.tar.gz
