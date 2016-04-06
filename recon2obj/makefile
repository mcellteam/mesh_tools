CC = gcc
DEFS =
FLAGS = -O3 
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = argparse.bison.o reconparse.bison.o strfunc.o recon2obj.o 

recon2obj: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o recon2obj

arglex.flex.c: arglex.l
	flex -Cfra arglex.l

argparse.bison.c: arglex.flex.c argparse.y recon2obj.h
	bison argparse.y

reconlex.flex.c: reconlex.l
	flex -Cfra reconlex.l

reconparse.bison.c: reconlex.flex.c reconparse.y recon2obj.h
	bison reconparse.y

reconparse.o: reconparse.c recon2obj.h
	$(CC) $(FLAGS) $(DEFS) -c reconparse.c

strfunc.o: strfunc.c strfunc.h
	$(CC) $(FLAGS) $(DEFS) -c strfunc.c

recon2obj.o: recon2obj.c recon2obj.h
	$(CC) $(FLAGS) $(DEFS) -c recon2obj.c

clean:
	rm -f $(OBJS) arglex.flex.c argparse.bison.c reconlex.flex.c reconparse.bison.c recon2obj

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp recon2obj $(INSTALL_DIR); \
	  cp recon2obj $(INSTALL_DIR); \
	fi

tar:
	tar -cvf recon2obj.tar lex.l parse.y recon2obj.c recon2obj.h makefile

gzip:
	tar -cvf - lex.l parse.y recon2obj.c recon2obj.h makefile | gzip -9c recon2obj.tar.gz
