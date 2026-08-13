CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = meshparse.bison.o mesh2obj.o 

mesh2obj: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh2obj

meshlex.flex.c: meshlex.l
	flex -Cfra meshlex.l

meshparse.bison.c: meshlex.flex.c meshparse.y mesh2obj.h
	bison meshparse.y

meshparse.bison.o: meshparse.bison.c mesh2obj.h
	$(CC) $(FLAGS) $(DEFS) -c meshparse.bison.c

mesh2obj.o: mesh2obj.c mesh2obj.h
	$(CC) $(FLAGS) $(DEFS) -c mesh2obj.c

clean:
	rm -f $(OBJS) meshlex.flex.c meshparse.bison.c mesh2obj

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
	tar -cvf mesh2obj.tar meshlex.l meshparse.y mesh2obj.c mesh2obj.h makefile

gzip:
	tar -cvf - meshlex.l meshparse.y mesh2obj.c mesh2obj.h makefile | gzip -9c mesh2obj.tar.gz
