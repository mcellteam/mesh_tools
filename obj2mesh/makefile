CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = objparse.bison.o obj2mesh.o 

obj2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o obj2mesh

objlex.flex.c: objlex.l
	flex -Cfra objlex.l

objparse.bison.c: objlex.flex.c objparse.y obj2mesh.h
	bison objparse.y

objparse.bison.o: objparse.bison.c obj2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c objparse.bison.c

obj2mesh.o: obj2mesh.c obj2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c obj2mesh.c

clean:
	rm -f $(OBJS) objlex.flex.c objparse.bison.c obj2mesh

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
	tar -cvf obj2mesh.tar objlex.l objparse.y obj2mesh.c obj2mesh.h makefile

gzip:
	tar -cvf - objlex.l objparse.y obj2mesh.c obj2mesh.h makefile | gzip -9c obj2mesh.tar.gz
