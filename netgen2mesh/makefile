CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = netgenparse.o strfunc.o netgen2mesh.o 

netgen2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o netgen2mesh

netgenlex.c: netgenlex.l
	flex -Cfra -Pnetgen -onetgenlex.c netgenlex.l

netgenparse.c: netgenlex.c netgenparse.y netgen2mesh.h
	bison -p netgen -o netgenparse.c netgenparse.y

netgenparse.o: netgenparse.c netgen2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c netgenparse.c

strfunc.o: strfunc.c strfunc.h
	$(CC) $(FLAGS) $(DEFS) -c strfunc.c

netgen2mesh.o: netgen2mesh.c netgen2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c netgen2mesh.c

clean:
	rm -f $(OBJS) netgenlex.c netgenparse.c netgen2mesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp netgen2mesh $(INSTALL_DIR); \
	  cp netgen2mesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf netgen2mesh.tar lex.l parse.y netgen2mesh.c netgen2mesh.h makefile

gzip:
	tar -cvf - lex.l parse.y netgen2mesh.c netgen2mesh.h makefile | gzip -9c netgen2mesh.tar.gz
