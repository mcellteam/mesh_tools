CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = netgenparse.o strfunc.o netgen2smesh.o 

netgen2smesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o netgen2smesh

netgenlex.c: netgenlex.l
	flex -Cfra -Pnetgen -onetgenlex.c netgenlex.l

netgenparse.c: netgenlex.c netgenparse.y netgen2smesh.h
	bison -p netgen -o netgenparse.c netgenparse.y

netgenparse.o: netgenparse.c netgen2smesh.h
	$(CC) $(FLAGS) $(DEFS) -c netgenparse.c

strfunc.o: strfunc.c strfunc.h
	$(CC) $(FLAGS) $(DEFS) -c strfunc.c

netgen2smesh.o: netgen2smesh.c netgen2smesh.h
	$(CC) $(FLAGS) $(DEFS) -c netgen2smesh.c

clean:
	rm -f $(OBJS) netgenlex.c netgenparse.c netgen2smesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp netgen2smesh $(INSTALL_DIR); \
	  cp netgen2smesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf netgen2smesh.tar lex.l parse.y netgen2smesh.c netgen2smesh.h makefile

gzip:
	tar -cvf - lex.l parse.y netgen2smesh.c netgen2smesh.h makefile | gzip -9c netgen2smesh.tar.gz
