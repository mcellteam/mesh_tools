CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = synu2mesh.o 

synu2mesh: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o synu2mesh

synu2mesh.o: synu2mesh.c synu2mesh.h
	$(CC) $(FLAGS) $(DEFS) -c synu2mesh.c

clean:
	rm -f $(OBJS) synu2mesh

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp synu2mesh $(INSTALL_DIR); \
	  cp synu2mesh $(INSTALL_DIR); \
	fi

tar:
	tar -cvf synu2mesh.tar synu2mesh.c synu2mesh.h makefile

gzip:
	tar -cvf - synu2mesh.c synu2mesh.h makefile | gzip -9c synu2mesh.tar.gz
