CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.bison.o imod2reconstruct.o 

imod2reconstruct: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o imod2reconstruct

lex.flex.c: lex.l
	flex -Cfra -o lex.flex.c lex.l

parse.bison.c: lex.flex.c parse.y imod2reconstruct.h
	bison -o parse.bison.c parse.y

parse.bison.o: parse.bison.c imod2reconstruct.h
	$(CC) $(FLAGS) $(DEFS) -c parse.bison.c

imod2reconstruct.o: imod2reconstruct.c imod2reconstruct.h
	$(CC) $(FLAGS) $(DEFS) -c imod2reconstruct.c

clean:
	rm -f $(OBJS) lex.flex.c parse.bison.c imod2reconstruct

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp imod2reconstruct $(INSTALL_DIR); \
	  cp imod2reconstruct $(INSTALL_DIR); \
	fi

tar:
	tar -cvf imod2reconstruct.tar lex.l parse.y imod2reconstruct.c imod2reconstruct.h makefile

gzip:
	tar -cvf - lex.l parse.y imod2reconstruct.c imod2reconstruct.h makefile | gzip -9c imod2reconstruct.tar.gz
