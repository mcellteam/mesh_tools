CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = stlb2stla.o 

stlb2stla: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o stlb2stla

stlb2stla.o: stlb2stla.c
	$(CC) $(FLAGS) $(DEFS) -c stlb2stla.c

clean:
	rm -f $(OBJS) stlb2stla

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp stl2mesh $(INSTALL_DIR); \
	  cp stlb2stla $(INSTALL_DIR); \
	fi

