CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = meshparse.bison.o meshscale.o 

meshscale: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o meshscale

meshlex.flex.c: meshlex.l
	flex -Cfra meshlex.l

meshparse.bison.c: meshlex.flex.c meshparse.y meshscale.h
	bison meshparse.y

meshparse.bison.o: meshparse.bison.c meshscale.h
	$(CC) $(FLAGS) $(DEFS) -c meshparse.bison.c

meshscale.o: meshscale.c meshscale.h
	$(CC) $(FLAGS) $(DEFS) -c meshscale.c

clean:
	rm -f $(OBJS) meshlex.flex.c meshparse.bison.c meshscale

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
	tar -cvf meshscale.tar meshlex.flex.l meshparse.y meshscale.c meshscale.h makefile

gzip:
	tar -cvf - meshlex.flex.l meshparse.y meshscale.c meshscale.h makefile | gzip -9c meshscale.tar.gz
