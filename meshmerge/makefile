CC = gcc
DEFS =
#FLAGS = -O3 -march=pentiumpro
FLAGS = -O3 -march=pentiumpro -L/home/jkinney/install/googleperf/lib -lprofiler -ltcmalloc
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o meshmerge.o 

meshmerge: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o meshmerge

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y meshmerge.h
	bison -o parse.c parse.y

parse.o: parse.c meshmerge.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

meshmerge.o: meshmerge.c meshmerge.h
	$(CC) $(FLAGS) $(DEFS) -c meshmerge.c

clean:
	rm -f $(OBJS) lex.c parse.c meshmerge

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp meshmerge $(INSTALL_DIR); \
	  cp meshmerge $(INSTALL_DIR); \
	fi

tar:
	tar -cvf meshmerge.tar lex.l parse.y meshmerge.c meshmerge.h makefile

gzip:
	tar -cvf - lex.l parse.y meshmerge.c meshmerge.h makefile | gzip -9c meshmerge.tar.gz
