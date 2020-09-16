CC = gcc
DEFS =
FLAGS = -O3 
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = argparse.bison.o reconparse.bison.o strfunc.o reconstruct_transform.o reconstruct2imod.o 

reconstruct2imod: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o reconstruct2imod

arglex.flex.c: arglex.l
	flex -Cfra arglex.l

argparse.bison.c: arglex.flex.c argparse.y reconstruct2imod.h
	bison argparse.y

reconlex.flex.c: reconlex.l
	flex -Cfra reconlex.l

reconparse.bison.c: reconlex.flex.c reconparse.y reconstruct2imod.h
	bison reconparse.y

reconparse.o: reconparse.c reconstruct2imod.h
	$(CC) $(FLAGS) $(DEFS) -c reconparse.c

strfunc.o: strfunc.c strfunc.h
	$(CC) $(FLAGS) $(DEFS) -c strfunc.c

reconstruct_transform.o: reconstruct_transform.c reconstruct_transform.h
	$(CC) $(FLAGS) $(DEFS) -c reconstruct_transform.c

reconstruct2imod.o: reconstruct2imod.c reconstruct2imod.h
	$(CC) $(FLAGS) $(DEFS) -c reconstruct2imod.c

clean:
	rm -f $(OBJS) arglex.flex.c argparse.bison.c reconlex.flex.c reconparse.bison.c reconstruct2imod

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp reconstruct2imod $(INSTALL_DIR); \
	  cp reconstruct2imod $(INSTALL_DIR); \
	fi

tar:
	tar -cvf reconstruct2imod.tar lex.l parse.y reconstruct2imod.c reconstruct2imod.h makefile

gzip:
	tar -cvf - lex.l parse.y reconstruct2imod.c reconstruct2imod.h makefile | gzip -9c reconstruct2imod.tar.gz
